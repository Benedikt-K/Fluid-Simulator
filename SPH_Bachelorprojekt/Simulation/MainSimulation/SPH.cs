using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using SPH_Bachelorprojekt.Simulation.Particles;
using SPH_Bachelorprojekt.Simulation.Kernel_Function;
using System.Numerics;

namespace SPH_Bachelorprojekt.Simulation.MainSimulation
{
    /// <summary>
    /// functions for general SPH (update of particles, density computation...) -> essentially all except pressure computation
    /// </summary>
    class SPH
    {
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

        public Vector2 CalculateViscosityAcceleration(Particle particle, List<Particle> neighbours, float ParticleSizeH, float viscosity, Kernel kernel)
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
            viscostiyAcc = 2 * viscosity * viscostiyAcc;

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

        public float CalculateParticleLambdaCFL(Vector2 velocity, float ParticleSizeH, float TimeStep)
        {
            float lambda = (TimeStep * velocity.Length()) / ParticleSizeH;
            return lambda;
        }
    }
}
