using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using SPH_Bachelorprojekt.Simulation.Particles;
using System.Numerics;

namespace SPH_Bachelorprojekt.Simulation.Kernel_Function
{
    class Quadratic
    {
        public List<Particle> Particles;
        public float ParticleSizeH;

        public Quadratic(List<Particle> particleList, float particleSizeH)
        {
            Particles = particleList;
            ParticleSizeH = particleSizeH;
        }

        public List<Particle> GetNeighboursQuadratic(Particle searchParticle)
        {
            List<Particle> neighbours = new List<Particle>();

            foreach(var potentialNeighbour in Particles)
            {
                if (Vector2.Distance(searchParticle.Position, potentialNeighbour.Position) <= 2.1f * ParticleSizeH)
                {
                    neighbours.Add(potentialNeighbour);
                    if (neighbours.Count >= 13)
                    {
                        //break;
                    }
                }
            }

            return neighbours;
        }
    }
}
