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
        public float Viscosity;
        public float Stiffness;
        public float minDensity;
        public float Gravity;

        public float AverageDensity;
        public float MaxCurrentParticlePressure;
        public float MaxVelocity;

        //later do own document
        public float RelaxationFactor;
        public float Gamma;

        public SimulationLoop(List<Particle> particles, float density, float particleSizeH, float timeStep, float viscosity, float stiffness, float gravity)
        {
            Particles = particles;
            ParticlesUpdated = new List<Particle>();
            SpatialHashing = new SpatialHashing((int) particleSizeH * 2);
            Density = density;
            ParticleSizeH = particleSizeH;
            TimeStep = timeStep;
            Viscosity = viscosity;
            Stiffness = stiffness;
            minDensity = density * 0.01f;
            Gravity = gravity;
            //minDensity = 0f;
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

        public void UpdateAllParticles(float smoothingLength, bool useIISPH, bool useNeighbour)
        {
            if (useIISPH)
            {
                UpdateAllParticlesIISPH(smoothingLength, useNeighbour);
            }
            else
            {
                UpdateAllParticlesSPH(smoothingLength, useNeighbour);
            }
        }

        public void UpdateAllParticlesIISPH(float smoothingLength, bool useNeighbour)
        {
            /*Kernel kernel = new Kernel(smoothingLength);
            PredictAdvection(kernel);
            PressureSolve(kernel);
            Integrate(kernel);*/
        }

        /*public void PredictAdvection(Kernel kernel)
        {
            foreach (Particle particle in Particles)
            {
                float density = 0.0f;
                foreach (Particle neighbour in particle.Neighbours)
                {
                    density += neighbour.Mass * kernel.W(Vector2.Distance(particle.Position, neighbour.Position));
                }
                particle.PredictedVelocity = particle.Velocity + TimeStep * (particle.PressureForce / particle.Mass); //calculate pressure forces

                Vector2 diagonalTerm = Vector2.Zero;
                foreach (Particle neighbour in particle.Neighbours)
                {
                    diagonalTerm += neighbour.Mass / particle.Density * kernel.GradW(particle.Position, neighbour.Position);
                }
                diagonalTerm *= TimeStep * TimeStep;
            }

            foreach (Particle particle in Particles)
            {
                float predictedDensity = particle.Density;
                foreach (Particle neighbour in particle.Neighbours)
                {
                    predictedDensity += TimeStep * neighbour.Mass * Vector2.Dot()//
                }
            }
        }*/

        /*public void PressureSolve(Kernel kernel)
        {
            int l = 0;
            float avgDensityError = 0;
            while (avgDensityError > RelaxationFactor || l < 2)
            {
                avgDensityError = 0.0f;
                foreach (Particle particle in Particles)
                {
                    Vector2 sum = Vector2.Zero;
                    foreach (Particle neighbour in particle.Neighbours)
                    {
                        sum += TimeStep * TimeStep * neighbour.Mass * neighbour.Pressure / (neighbour.Density * neighbour.Density) * kernel.GradW(particle.Position, neighbour.Position);
                    }
                    particle.NewPressure = particle.Pressure + Gamma * (Density - particle.Density + sum);
                }
            }
        }*/

        public void Integrate(Kernel kernel)
        {

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
                    SpatialHashing.InsertObject(particle);
                }
                foreach (Particle particle in Particles)
                {
                    SpatialHashing.InRadius(particle.Position, ParticleSizeH * 2f, ref particle.Neighbours);
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
            //AverageDensity /= nonBoundaryParticlesCount;
            AverageDensity /= Particles.Count;

            // Compute accelerations
            Parallel.ForEach(Particles, particle =>
            {
                if (!particle.IsBoundaryParticle)
                {
                    Vector2 acceleration_T = CalculatePressureAcceleration(particle, particle.Neighbours, kernel, particle.Density);
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
                        /*Parallel.ForEach(Particles, particle =>
            {
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
            });*/
            //CalculateNewTimeStep();
            // Transfer updated particles and reset 
            //Particles = ParticlesUpdated;
            //ParticlesUpdated = new List<Particle>();

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
                    SpatialHashing.InRadius(particle.Position, ParticleSizeH * 2f, ref particle.Neighbours);
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
                    Vector2 acceleration_T = CalculatePressureAcceleration(particle, particle.Neighbours, kernel, particle.Density);
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
                }
            });

            // Update velocitys and positions
            foreach (Particle particle in Particles)
            {
                //deleting hash particles
                SpatialHashing.RemoveObject(particle);
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
                SpatialHashing.InsertObject(particle);
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
                return 0f;
            }
            foreach (Particle neighbour in neighbours)
            {
                density += neighbour.Mass * kernel.W(Vector2.Distance(particle.Position, neighbour.Position));
            }
            return density;
        }

        public Vector2 CalculatePressureAcceleration(Particle particle, List<Particle> neighbours, Kernel kernel, float density)
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
    }
}
