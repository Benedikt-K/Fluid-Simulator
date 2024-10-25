using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Numerics;
using SPH_Bachelorprojekt.Simulation.Particles;
using SPH_Bachelorprojekt.Simulation.Kernel_Function;

namespace SPH_Bachelorprojekt.Simulation.MainSimulation
{
    /// <summary>
    /// functions to compute the pressure using a global pressure computation
    /// </summary>

    class IISPH
    {
        public static float GetDiagonalElement(Particle particle, float particleSizeH, float TimeStep, float Gamma, Kernel kernel)
        {
            float diagonalElement = 0;
            foreach (Particle neighbour in particle.Neighbours)
            {
                float particleLastDensity2 = particle.LastDensity * particle.LastDensity;
                if (neighbour.IsBoundaryParticle)
                {
                    // BOUNDARY N
                    Vector2 innerTerm = Vector2.Zero;
                    foreach (Particle neighbourInner in particle.Neighbours)
                    {
                        if (neighbourInner.IsBoundaryParticle)
                        {
                            // boundary NN
                            innerTerm -= 2 * Gamma * neighbourInner.Mass / particleLastDensity2 * kernel.GradW(particle.Position, neighbourInner.Position);
                        }
                        else
                        {
                            // fluid NN
                            innerTerm -= neighbourInner.Mass / particleLastDensity2 * kernel.GradW(particle.Position, neighbourInner.Position);
                        }
                    }
                    float dotProduct = Vector2.Dot(innerTerm, kernel.GradW(particle.Position, neighbour.Position));
                    diagonalElement += neighbour.Mass * dotProduct;
                }
                else
                {
                    // FLUID N
                    Vector2 innerTerm = Vector2.Zero;
                    foreach (Particle neighbourInner in particle.Neighbours)
                    {
                        if (neighbourInner.IsBoundaryParticle)
                        {
                            // boundary NN
                            innerTerm -= 2 * Gamma * neighbourInner.Mass / particleLastDensity2 * kernel.GradW(particle.Position, neighbourInner.Position);
                        }
                        else
                        {
                            // fluid NN
                            innerTerm -= neighbourInner.Mass / particleLastDensity2 * kernel.GradW(particle.Position, neighbourInner.Position);
                        }
                    }
                    float dotProduct = Vector2.Dot(innerTerm, kernel.GradW(particle.Position, neighbour.Position));
                    diagonalElement += neighbour.Mass * dotProduct;

                    //second fluid N term
                    Vector2 otherInnerTerm = Vector2.Zero;
                    otherInnerTerm = (particle.Mass / particleLastDensity2) * kernel.GradW(neighbour.Position, particle.Position);
                    float otherdotProduct = Vector2.Dot(otherInnerTerm, kernel.GradW(particle.Position, neighbour.Position));
                    diagonalElement += neighbour.Mass * otherdotProduct;
                }
            }
            diagonalElement *= TimeStep * TimeStep;
            return diagonalElement;
        }


        public static float GetSourceTerm(Particle particle, float particleSizeH, float TimeStep, float ElapsedTime, float fluidDensity, Kernel kernel)
        {
            /////////////// mybe here error ////////// check again
            float predictedDensityError = fluidDensity - particle.Density; // maybe LastDensity
            foreach (Particle neighbour in particle.Neighbours)
            {
                if (neighbour.IsBoundaryParticle)
                {
                    //boundary particles
                    float dotProduct = Vector2.Dot(particle.PredictedVelocity - neighbour.PredictedVelocity * (ElapsedTime + TimeStep), kernel.GradW(particle.Position, neighbour.Position));
                    predictedDensityError -= TimeStep * neighbour.Mass * dotProduct;
                }
                else
                {
                    // FLUID particles
                    float dotProduct = Vector2.Dot(particle.PredictedVelocity - neighbour.PredictedVelocity, kernel.GradW(particle.Position, neighbour.Position));
                    predictedDensityError -= TimeStep * neighbour.Mass * dotProduct;
                }
            }
            return predictedDensityError;
        }

        public static float GetDivergence(Particle particle, float timeStep, Kernel kernel)
        {
            float timeStep2 = timeStep * timeStep;
            float Ap = 0f;
            // compute divergence of velocity change due to pressureAcc
            foreach (Particle neighbour in particle.Neighbours) // OLD = foreach (Particle neighbour in Particles)
            {
                if (neighbour.IsBoundaryParticle)
                {
                    float dotProduct = Vector2.Dot(particle.PressureAcceleration, kernel.GradW(particle.Position, neighbour.Position));
                    Ap += neighbour.Mass * dotProduct;
                }
                else
                {
                    float dotProduct = Vector2.Dot(particle.PressureAcceleration - neighbour.PressureAcceleration, kernel.GradW(particle.Position, neighbour.Position));
                    Ap += neighbour.Mass * dotProduct;
                }
            }
            Ap *= timeStep2;
            return Ap;
        }

        public static void UpdateParticle(ref List<Particle> Particles, float TimeStep, Kernel kernel)
        {
            ///
            /// update particle positions and velocitys
            /// 
            foreach (Particle particle in Particles)
            {
                if (!particle.IsBoundaryParticle)
                {
                    particle.Velocity = TimeStep * particle.PressureAcceleration + particle.PredictedVelocity;
                    particle.Position += TimeStep * particle.Velocity;
                    particle.Pressure = particle.PredictedPressure;
                }
            }
        }

        public static void PressureSolve(ref List<Particle> Particles, int min_Iterations, int max_Iterations, float max_error_Percentage, float Density, float Gamma, Kernel kernel)
        {
            ///
            /// calculate Pressures of all particles
            ///
            int currentIteration = 0;
            float averageDensityError = float.PositiveInfinity;
            bool continueWhile = true;
            float percentageDensityError = float.PositiveInfinity;

            while ((continueWhile || (currentIteration <= min_Iterations)) && (currentIteration < max_Iterations))
            {
                currentIteration++;
                averageDensityError = 0;
                DoPressureSolveIteration(ref Particles, Gamma, kernel, ref averageDensityError);
                percentageDensityError = averageDensityError / Density;
                float eta = max_error_Percentage * 0.01f * Density;
                continueWhile = averageDensityError >= eta;
                Console.WriteLine("iter: " + currentIteration + ", err: " + averageDensityError);

            }
            Console.WriteLine("iterations needed: " + currentIteration);
        }

        public static void DoPressureSolveIteration(ref List<Particle> Particles, float Gamma, Kernel kernel, ref float averageDensityError)
        {
            foreach (Particle particle in Particles)
            {
                if (!particle.IsBoundaryParticle)
                {
                    float particleLastDensity2 = particle.LastDensity * particle.LastDensity;
                    Vector2 pressureAcceleration = Vector2.Zero;
                    foreach (Particle neighbour in particle.Neighbours)
                    {
                        if (particle.IsBoundaryParticle)
                        {
                            pressureAcceleration -= Gamma * neighbour.Mass * 2 * (particle.PredictedPressure / particleLastDensity2) * kernel.GradW(particle.Position, neighbour.Position);
                        }
                        else
                        {
                            float neighbourLastDensity2 = neighbour.LastDensity * neighbour.LastDensity;
                            float innerTerm = (particle.PredictedPressure / particleLastDensity2) + (neighbour.PredictedPressure / neighbourLastDensity2);
                            pressureAcceleration -= neighbour.Mass * innerTerm * kernel.GradW(particle.Position, neighbour.Position);
                        }
                    }
                    particle.PressureAcceleration = pressureAcceleration;
                }
            }
        }
    }
}
