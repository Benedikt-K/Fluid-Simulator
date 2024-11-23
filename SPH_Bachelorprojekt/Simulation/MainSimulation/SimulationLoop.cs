using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Numerics;
using SPH_Bachelorprojekt.Simulation.Particles;
using SPH_Bachelorprojekt.Simulation.Kernel_Function;
using SPH_Bachelorprojekt.Simulation.Neighbours;

namespace SPH_Bachelorprojekt.Simulation.MainSimulation
{
    class SimulationLoop
    {
        public List<Particle> Particles;
        public List<Particle> ParticlesUpdated;
        public SpatialHashing SpatialHashing;
        public float Density;
        public float ParticleSizeH;
        public float TimeStep;
        public float ElapsedTime;
        public float Viscosity;
        public float Stiffness;
        public float minDensity;
        public float Gravity;
        public int FluidParticleCount;
        public int DensityErrorIterCollectionCount;

        public float AverageDensity;
        public float MaxCurrentParticlePressure;
        public float MaxVelocity;

        //later do own document
        public float Omega;
        public float Gamma;
        public float Lambda;
        public float LambdaSESPH;
        public float MaxTimestep;
        public float MinTimeStep;
        public float MaxTimestepSESPH;
        public float MinTimeStepSESPH;

        public SimulationLoop(List<Particle> particles, float density, float particleSizeH, float timeStep, float viscosity, float stiffness, float gravity)
        {
            Particles = particles;
            ParticlesUpdated = new List<Particle>();
            SpatialHashing = new SpatialHashing((int) particleSizeH * 2);
            Density = density;
            ParticleSizeH = particleSizeH;
            TimeStep = timeStep;
            ElapsedTime = 0f;
            Viscosity = viscosity;
            Stiffness = stiffness;
            minDensity = density * 0.01f;
            Gravity = gravity;
            Omega = 0.5f;
            Gamma = 0.7f;
            Lambda = 0.7f;
            LambdaSESPH = 0.1f;
            MaxTimestep = 0.005f;
            MinTimeStep = 0.00005f;
            MaxTimestepSESPH = 0.005f;
            MinTimeStepSESPH = 0.0005f;
            //minDensity = 0f;
            int fluidCount = 0;
            foreach (Particle particle in Particles)
            {
                if (!particle.IsBoundaryParticle)
                {
                    fluidCount++;
                }
            }
            FluidParticleCount = fluidCount;
            DensityErrorIterCollectionCount = 0;
    }

    public SimulationLoop(float density, float particleSizeH, float timeStep, float viscosity)
        {
            Particles = new List<Particle>();
            ParticlesUpdated = new List<Particle>();
            Density = density;
            ParticleSizeH = particleSizeH;
            TimeStep = timeStep;
            Viscosity = viscosity;
        }

        public void UpdateAllParticles(float smoothingLength, bool useIISPH, bool useNeighbour, ref List<float> densityErrorData, ref List<float> iterationData, bool CollectDensityIterErr)
        {
            if (useIISPH)
            {
                UpdateAllParticles_IISPH(smoothingLength, useNeighbour, ref densityErrorData, ref iterationData, CollectDensityIterErr);
            }
            else
            {
                UpdateAllParticlesSPH(smoothingLength, useNeighbour);
            }
            ElapsedTime += TimeStep;
        }


        public void UpdateAllParticles_IISPH(float smoothingLength, bool useNeighbour, ref List<float> densityErrorData, ref List<float> iterationData, bool CollectDensityIterErr)
        {
            SPH.NeighbourhoodSearch(ref Particles, ParticleSizeH, useNeighbour);
            // Implementation based on "https://cg.informatik.uni-freiburg.de/publications/2013_TVCG_IISPH.pdf" and notes from Prof. Teschner
            Kernel kernel = new Kernel(smoothingLength);
            PredictAdvection(kernel);
            //IISPH.PressureSolve(ref Particles, min_Iterations, max_Iterations, max_error_Percentage, Density, Gamma, kernel);
            PressureSolve(kernel, ref densityErrorData, ref iterationData, ref iterationData, CollectDensityIterErr);
            UpdateIISPH(kernel);
        }

        
        public void PredictAdvection(Kernel kernel)
        {
            AverageDensity = 0;
            Parallel.ForEach(Particles, particle =>
            {
                if (!particle.IsBoundaryParticle)
                {
                    float density = SPH.CalculateDensityAtParticle(particle, kernel);
                    AverageDensity += density;
                    particle.Density = density;
                    Vector2 nonPressureAcceleration = CalculateViscosityAcceleration(particle, particle.Neighbours, kernel) + GetGravity(); ///add surface tension
                    particle.NonPressureAcceleration = nonPressureAcceleration;
                    particle.PredictedVelocity = particle.Velocity + TimeStep * nonPressureAcceleration;
                }
            });
            AverageDensity /= FluidParticleCount;
            //Console.WriteLine("averageDensity : " + AverageDensity);

            Parallel.ForEach(Particles, particle =>
            {
                if (!particle.IsBoundaryParticle)
                {
                    float densityFac = particle.Density / Density;
                    float densityFac2 = densityFac * densityFac;
                    particle.SourceTerm = IISPH.GetSourceTerm(particle, ParticleSizeH, TimeStep, ElapsedTime, Density, kernel);
                    particle.DiagonalElement = IISPH.GetDiagonalElement(particle, ParticleSizeH, TimeStep, Gamma, Density, kernel);
                    //particle.Pressure = 0.5f * particle.Pressure * (ElapsedTime - TimeStep);
                    //particle.Pressure = Math.Max(0, Omega * (particle.SourceTerm / particle.DiagonalElement));
                    particle.Pressure = 0;
                }
            });
        }

        

        public void PressureSolve(Kernel kernel, ref List<float> densityErrorData, ref List<float> iterationData, ref List<float> iterationCountData, bool collectAverageDensityErrIter)
        {
            ///
            /// calculate Pressures of all particles
            ///
            int min_Iterations = 3;
            int max_Iterations = 30;
            float max_error_Percentage = 0.5f; // given in %
            // dislocate to other file
            int currentIteration = 0;
            float averageDensityError = 0;
            bool continueWhile = true;


            while ((continueWhile || (currentIteration < min_Iterations)) && (currentIteration < max_Iterations))
            {
                averageDensityError = 0;
                DoPressureSolveIteration(kernel, ref averageDensityError);
                float eta = max_error_Percentage * 0.01f * Density;
                float absoluteAverageDensityError = Math.Abs(averageDensityError);
                continueWhile = absoluteAverageDensityError >= eta;
                // add data for graph
                if (collectAverageDensityErrIter && currentIteration > 0)
                {
                    densityErrorData.Add(((absoluteAverageDensityError + Density) - Density) / Density * 100); // get DensityError in %
                    iterationData.Add(currentIteration);
                }
                currentIteration++;
            }
            //Console.WriteLine("iterations needed: " + currentIteration);
            Console.WriteLine("iter: " + currentIteration + "---Err: " + (((Math.Abs(averageDensityError) + Density) - Density) / Density * 100));
        }

        public void DoPressureSolveIteration(Kernel kernel, ref float averageDensityError)
        {
            foreach (Particle particle in Particles)
            {
                if (!particle.IsBoundaryParticle)
                {
                    // compute pressure acceleration
                    float fluidDensity2 = Density * Density;
                    Vector2 pressureAcceleration = Vector2.Zero;
                    foreach (Particle neighbour in particle.Neighbours)
                    {
                        if (neighbour.IsBoundaryParticle)
                        {
                            pressureAcceleration -= Gamma * neighbour.GetMass() * 2 * (particle.Pressure / fluidDensity2) * kernel.GradW(particle.Position, neighbour.Position);
                        }
                        else
                        {
                            float innerTerm = (particle.Pressure / fluidDensity2) + (neighbour.Pressure / (fluidDensity2));
                            pressureAcceleration -= neighbour.GetMass() * innerTerm * kernel.GradW(particle.Position, neighbour.Position);
                        }
                    }
                    particle.PressureAcceleration = pressureAcceleration;
                }
            }

            
            float timeStep2 = TimeStep * TimeStep;
            float omega = Omega;

            foreach (Particle particle in Particles) 
            { 
                if (!particle.IsBoundaryParticle) 
                {
                    float Ap = 0f;
                    // compute divergence of velocity change due to pressureAcc
                    foreach (Particle neighbour in particle.Neighbours)
                    {
                        if (neighbour.IsBoundaryParticle)
                        {
                            float dotProduct = Vector2.Dot(particle.PressureAcceleration, kernel.GradW(particle.Position, neighbour.Position));
                            Ap += neighbour.GetMass() * dotProduct;
                        }
                        else
                        {
                            float dotProduct = Vector2.Dot(particle.PressureAcceleration - neighbour.PressureAcceleration, kernel.GradW(particle.Position, neighbour.Position));
                            Ap += neighbour.GetMass() * dotProduct;
                        }
                    }
                    Ap *= timeStep2;

                    // update pressure
                    if (particle.DiagonalElement != 0)
                    {
                        float innerTerm = particle.Pressure + omega * ((particle.SourceTerm - Ap) / particle.DiagonalElement);
                        particle.Pressure = Math.Max(innerTerm, 0);
                    }
                    if (!Single.IsNaN(Math.Abs(Ap - particle.SourceTerm)))
                    {
                        averageDensityError += Math.Abs(Ap - particle.SourceTerm);
                    }
                }
            }
            averageDensityError /= FluidParticleCount;
        }

        public void UpdateIISPH(Kernel kernel)
        {
            // update timestep
            float maxVelocity = Particles.Max(p => p.Velocity.Length());
            float maxAcceleration = Particles.Max(p => (p.NonPressureAcceleration.Length() + p.PressureAcceleration.Length()));

            float velocityConstrain = ParticleSizeH / maxVelocity;
            float accelerationConstraint = ParticleSizeH / maxAcceleration;

            float newTimeStep = Lambda * Math.Min(maxVelocity, maxAcceleration);
            newTimeStep = Math.Min(MinTimeStep, newTimeStep);
            newTimeStep = Math.Max(MaxTimestep, newTimeStep);

            TimeStep = newTimeStep;
            ///
            /// update particle positions and velocitys
            /// 
            foreach (Particle particle in Particles)
            {
                if (!particle.IsBoundaryParticle)
                {
                    // old update
                    particle.Velocity = TimeStep * particle.PressureAcceleration + particle.PredictedVelocity;
                    particle.Position += TimeStep * particle.Velocity;
                }
            }
        }

            public void UpdateAllParticlesSPH(float smoothingLength, bool useNeighbour)
        { 
            float restDensity = Density;
            float gravity = Gravity;

            ///
            /// Neighbour search
            ///
            if (useNeighbour)
            {
                int CellSize = (int) ParticleSizeH * 2; /// Particle Size should be int here
                SpatialHashing = new SpatialHashing(CellSize);
                foreach (Particle particle in Particles)
                {
                    SpatialHashing.AddParticle(particle);
                }
                foreach (Particle particle in Particles)
                {
                    SpatialHashing.GetNeighbours(particle.Position, ParticleSizeH * 2f, ref particle.Neighbours);
                }
            }
            else
            {
                Quadratic quadraticSolver = new Quadratic(Particles, ParticleSizeH);             
                Parallel.ForEach(Particles, particle =>
                {
                    //Quadratic neighbour
                    List <Particle> neighbours = quadraticSolver.GetNeighboursQuadratic(particle);

                    // null reference check
                    if (neighbours != null)
                    {
                        particle.Neighbours = neighbours;
                        if (!particle.IsBoundaryParticle)
                        {
                            //Console.WriteLine("Neighbours Count: " + particle.Neighbours.Count);
                        }
                    }
                    else
                    {
                        particle.Neighbours = new List<Particle>();
                    }

                });
            }


            // Compute Density and Pressure
            Kernel kernel = new Kernel(smoothingLength);

            // Test kernel
            //kernel.TestKernel();
            //CalculateNewTimeStep();

            AverageDensity = 0f;
            MaxCurrentParticlePressure = 0f;
            int nonBoundaryParticlesCount = 0;
            foreach (Particle particle in Particles)
            {
                if (!particle.IsBoundaryParticle)
                {
                    float newDensity = CalculateDensityAtParticle(particle, particle.Neighbours, kernel);   // calculate density only for fluid particles
                    particle.Density = newDensity;
                }
                //float newDensity = CalculateDensityAtParticle(particle, particle.Neighbours, kernel);
                //particle.Density = newDensity;
                //Console.WriteLine("Density: " + particle.Density);
                float pressure = Math.Max(Stiffness * ((particle.Density / (restDensity)) - 1), 0); // was negative before

                if (pressure > MaxCurrentParticlePressure)
                {
                    MaxCurrentParticlePressure = pressure;
                }

                particle.Pressure = pressure;
                if (!particle.IsBoundaryParticle)
                {
                    AverageDensity += particle.Density;
                    nonBoundaryParticlesCount++;
                    //Console.WriteLine(particle.Pressure);
                }
                //AverageDensity += particle.Density;
            }
            /*Parallel.ForEach (Particles, particle =>
            {
                if (!particle.IsBoundaryParticle)
                {
                    float newDensity = CalculateDensityAtParticle(particle, particle.Neighbours, kernel);   // calculate density only for fluid particles
                    particle.Density = newDensity;
                }
                //float newDensity = CalculateDensityAtParticle(particle, particle.Neighbours, kernel);
                //particle.Density = newDensity;
                //Console.WriteLine("Density: " + particle.Density);
                float pressure = Math.Max(Stiffness * ((particle.Density / (restDensity)) - 1), 0); // was negative before

                if (pressure > MaxCurrentParticlePressure)
                {
                    MaxCurrentParticlePressure = pressure;
                }

                // Other pressure calculation
                //float pressure = Stiffness * (particle.Density - restDensity);

                //float pressure = Stiffness * ((particle.Density / restDensity) - 1);
                //float pressure = -Stiffness * ((particle.Density / restDensity) - 1);
                particle.Pressure = pressure;
                if (!particle.IsBoundaryParticle)
                {
                    AverageDensity += particle.Density;
                    nonBoundaryParticlesCount++;
                    //Console.WriteLine(particle.Pressure);
                }
                //AverageDensity += particle.Density;
            });*/
            AverageDensity /= nonBoundaryParticlesCount;
            //AverageDensity /= Particles.Count;
            Console.WriteLine(Math.Abs(AverageDensity - Density) / Density * 100);

            // Compute accelerations
            Parallel.ForEach(Particles, particle =>
            {
                if (!particle.IsBoundaryParticle)
                {
                    Vector2 acceleration_T = CalculatePressureAcceleration(particle, particle.Neighbours, kernel);
                    /*if (!particle.IsBoundaryParticle)
                    {
                        Console.WriteLine("pressureAcc: " + acceleration_T);
                        Console.WriteLine("viscosity Acc: " + (CalculateViscosityAcceleration(particle, particle.Neighbours, kernel)));
                    }*/

                    acceleration_T += CalculateViscosityAcceleration(particle, particle.Neighbours, kernel);
                    acceleration_T += CalculateSurfaceTensionAcceleration(particle, particle.Neighbours, kernel);
                    if (particle.Mass != 0)
                    {
                        //acceleration_T = acceleration_T / particle.Mass;
                    }
                    else
                    {
                        acceleration_T = Vector2.Zero;
                    }
                    // Add gravity
                    particle.Acceleration = acceleration_T + new Vector2(0, gravity);
                    /*if (!particle.IsBoundaryParticle)
                    {
                        Console.WriteLine("both: " + particle.Acceleration);
                        Console.WriteLine("-------------");
                    }*/
                }
            });

            // Update velocitys and positions
            foreach(Particle particle in Particles)
            {
                // update timestep
                float maxVelocity = Particles.Max(p => p.Velocity.Length());
                float maxAcceleration = Particles.Max(p => (p.Acceleration.Length()));

                float velocityConstrain = ParticleSizeH / maxVelocity;
                float accelerationConstraint = ParticleSizeH / maxAcceleration;

                float newTimeStep = LambdaSESPH * Math.Min(maxVelocity, maxAcceleration);
                newTimeStep = Math.Min(MinTimeStepSESPH, newTimeStep);
                newTimeStep = Math.Max(MaxTimestepSESPH, newTimeStep);

                //deleting and later add hash particles (update basically)
                //SpatialHashing.RemoveObject(particle);
                if (particle.IsBoundaryParticle)
                {
                    //continue;
                }
                else
                {
                    UpdateParticle(particle, particle.Neighbours, kernel);
                    if (particle.Velocity.Length() > MaxVelocity)
                    {
                        MaxVelocity = particle.Velocity.Length();
                    }
                }
                //SpatialHashing.InsertObject(particle);
            }
        }

        public void UpdateAllParticlesSPHNotParallel(float smoothingLength, bool useNeighbour)
        {
            float restDensity = Density;
            float gravity = Gravity;

            ///
            /// Neighbour search
            ///
            if (useNeighbour)
            {
                foreach (Particle particle in Particles) 
                {
                    SpatialHashing.GetNeighbours(particle.Position, ParticleSizeH * 2f, ref particle.Neighbours);
                }
            }
            else
            {
                Quadratic quadraticSolver = new Quadratic(Particles, ParticleSizeH);
                foreach (Particle particle in Particles)
                {
                    // Quadratic neighbour 
                    List<Particle> neighbours = quadraticSolver.GetNeighboursQuadratic(particle);
                    // null reference check
                    if (neighbours != null)
                    {
                        particle.Neighbours = neighbours;
                        if (!particle.IsBoundaryParticle)
                        {
                            //Console.WriteLine("Neighbours Count: " + particle.Neighbours.Count);
                        }
                    }
                    else
                    {
                        particle.Neighbours = new List<Particle>();
                    }
                }
            }

            //Console.WriteLine("Number Particles: " + Particles.Count);
            //Console.WriteLine("Number Particles: " + Particles.Count);


            // Compute Density and Pressure
            Kernel kernel = new Kernel(smoothingLength);

            // Test kernel
            //kernel.TestKernel();
            //CalculateNewTimeStep();

            AverageDensity = 0f;
            MaxCurrentParticlePressure = 0f;
            int nonBoundaryParticlesCount = 0;
            foreach (Particle particle in Particles)
            {
                if (!particle.IsBoundaryParticle)
                {
                    float newDensity = CalculateDensityAtParticle(particle, particle.Neighbours, kernel);   // calculate density only for fluid particles
                    particle.Density = newDensity;
                }
                //float newDensity = CalculateDensityAtParticle(particle, particle.Neighbours, kernel);
                //particle.Density = newDensity;
                //Console.WriteLine("Density: " + particle.Density);
                float pressure = Math.Max(Stiffness * ((particle.Density / (restDensity)) - 1), 0); // was negative before

                if (pressure > MaxCurrentParticlePressure)
                {
                    MaxCurrentParticlePressure = pressure;
                }

                // Other pressure calculation
                //float pressure = Stiffness * (particle.Density - restDensity);

                //float pressure = Stiffness * ((particle.Density / restDensity) - 1);
                //float pressure = -Stiffness * ((particle.Density / restDensity) - 1);
                particle.Pressure = pressure;
                if (!particle.IsBoundaryParticle)
                {
                    //AverageDensity += particle.Density;
                    nonBoundaryParticlesCount++;
                    //Console.WriteLine(particle.Pressure);
                }
                AverageDensity += particle.Density;
            }
            
            //AverageDensity /= nonBoundaryParticlesCount;
            AverageDensity /= Particles.Count;

            // Compute accelerations
            Parallel.ForEach(Particles, particle =>
            {
                if (!particle.IsBoundaryParticle)
                {
                    Vector2 acceleration_T = CalculatePressureAcceleration(particle, particle.Neighbours, kernel);
                    acceleration_T += CalculateViscosityAcceleration(particle, particle.Neighbours, kernel);
                    acceleration_T += CalculateSurfaceTensionAcceleration(particle, particle.Neighbours, kernel);
                    if (particle.Mass != 0)
                    {
                        //acceleration_T = acceleration_T / particle.Mass;
                    }
                    else
                    {
                        acceleration_T = Vector2.Zero;
                    }
                    // Add gravity
                    particle.Acceleration = acceleration_T + GetGravity();
                }
            });

            // Update velocitys and positions
            foreach (Particle particle in Particles)
            {
                //deleting hash particles
                SpatialHashing.RemoveParticle(particle);
                if (particle.IsBoundaryParticle)
                {
                    //continue;
                }
                else
                {
                    UpdateParticle(particle, particle.Neighbours, kernel);
                    if (particle.Velocity.Length() > MaxVelocity)
                    {
                        MaxVelocity = particle.Velocity.Length();
                    }
                }
                SpatialHashing.AddParticle(particle);
            }
        }

        public void UpdateParticleOld(Particle particle, List<Particle> neighbours, Kernel kernel)
        {
            Vector2 newPosition = new Vector2();
            if (particle.IsBoundaryParticle)
            {
                newPosition = particle.Position;
            }
            else 
            {
                Vector2 velocity_TdT = particle.Velocity + TimeStep * particle.Acceleration;
                newPosition = particle.Position + TimeStep * velocity_TdT;         
            }
            ParticlesUpdated.Add(new Particle(newPosition, particle.Velocity, particle.Density, particle.Mass , ParticleSizeH, particle.IsBoundaryParticle));
        }

        public void UpdateParticle(Particle particle, List<Particle> neighbours, Kernel kernel)
        {
            // update timestep
            /*float maxVelocity = Particles.Max(p => p.Velocity.Length());
            float maxAcceleration = Particles.Max(p => (p.Acceleration.Length()));

            float velocityConstrain = ParticleSizeH / maxVelocity;
            float accelerationConstraint = ParticleSizeH / maxAcceleration;

            float newTimeStep = LambdaSESPH * Math.Min(maxVelocity, maxAcceleration);
            newTimeStep = Math.Min(MinTimeStep, newTimeStep);
            newTimeStep = Math.Max(MaxTimestep, newTimeStep);

            TimeStep = newTimeStep;*/
            //SpatialHashing.RemoveObject(particle);
            particle.Velocity += TimeStep * particle.Acceleration;
            particle.Position += TimeStep * particle.Velocity;
            //SpatialHashing.InsertObject(particle);
        }

        public float CalculateDensityAtParticle(Particle particle, List<Particle> neighbours, Kernel kernel)
        {
            float density = 0f;
            if (particle.Neighbours.Count <= 1)
            {
                return Density; /// TODO ? 0f oder RestDensity
            }
            foreach (Particle neighbour in neighbours)
            {
                density += neighbour.Mass * kernel.W(Vector2.Distance(particle.Position, neighbour.Position));
            }
            return density;
        }

        public Vector2 GetGravity()
        {
            return new Vector2(0, Gravity);
        }

        public Vector2 CalculatePressureAcceleration(Particle particle, List<Particle> neighbours, Kernel kernel)
        {
            // boundary handling with mirroring
            Vector2 pressureAcceleration = Vector2.Zero;
            Vector2 pressureAccelerationBoundary = Vector2.Zero;
            float pressureOverDensity2 = particle.Pressure / (particle.Density * particle.Density);
            foreach (Particle neighbour in neighbours)
            {
                Vector2 GradW = kernel.GradW(particle.Position, neighbour.Position);
                if (neighbour.IsBoundaryParticle)
                {
                    pressureAccelerationBoundary += particle.Mass * (pressureOverDensity2 + pressureOverDensity2) * GradW; //changed particle to neighbour
                    continue; 
                }
                float pressureOverDenisty2Neighbour = neighbour.Pressure / (neighbour.Density * neighbour.Density);
                pressureAcceleration += neighbour.Mass * (pressureOverDensity2 + pressureOverDenisty2Neighbour) * GradW;
            }
            Vector2 totalAcceleration = Vector2.Zero - pressureAcceleration - pressureAccelerationBoundary;
            return totalAcceleration;
        }

        public Vector2 CalculateViscosityAcceleration(Particle particle, List<Particle> neighbours, Kernel kernel)
        {
            Vector2 viscostiyAcc = Vector2.Zero;
            foreach (Particle neighbour in particle.Neighbours)
            {
                Vector2 v_ij = particle.Velocity - neighbour.Velocity;
                Vector2 x_ij = particle.Position - neighbour.Position;

                Vector2 gradW = kernel.GradW(particle.Position, neighbour.Position);
                float massDensityRatio = neighbour.Mass / neighbour.Density;
                float scaledParticleSize = 0.01f * ParticleSizeH * ParticleSizeH;
                float v_dot_x = Vector2.Dot(v_ij, x_ij);
                float x_dot_x = Vector2.Dot(x_ij, x_ij);
                viscostiyAcc += massDensityRatio * (v_dot_x / (x_dot_x + scaledParticleSize)) * gradW;
            }
            viscostiyAcc = 2 * Viscosity * viscostiyAcc;

            return viscostiyAcc;
        }

        public Vector2 CalculateSurfaceTensionAcceleration(Particle particle, List<Particle> neighbours, Kernel kernel)
        {
            // Have the feeling that it doesnt work properly
            Vector2 surfaceTensionAcc = Vector2.Zero;

            foreach (Particle neighbor in particle.Neighbours)
            {
                Vector2 gradW = kernel.GradW(particle.Position, neighbor.Position);
                float factor = 0f;         // works with 10, but gets not completely round
                surfaceTensionAcc += gradW * factor;
            }

            return surfaceTensionAcc;
        }

        public float CalculateParticleLambdaCFL(Vector2 velocity)
        {
            float lambda = (TimeStep * velocity.Length()) / ParticleSizeH;
            return lambda;
        }

        public float CalculateParticleLambdaCFL(float velocityLength)
        {
            float lambda = (TimeStep * velocityLength) / ParticleSizeH;
            return lambda;
        }
    }
}
