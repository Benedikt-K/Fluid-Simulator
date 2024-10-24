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
    /// functions to compute the pressure using a state equation
    /// </summary>
    class SESPH
    {
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
                    pressureAccelerationBoundary += particle.Mass * (pressureOverDensity2 + pressureOverDensity2) * GradW;
                    continue;
                }
                float pressureOverDenisty2Neighbour = neighbour.Pressure / (neighbour.Density * neighbour.Density);
                pressureAcceleration += neighbour.Mass * (pressureOverDensity2 + pressureOverDenisty2Neighbour) * GradW;
            }
            Vector2 totalAcceleration = Vector2.Zero - pressureAcceleration - pressureAccelerationBoundary;
            return totalAcceleration;
        }
    }
}
