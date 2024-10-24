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
    }
}
