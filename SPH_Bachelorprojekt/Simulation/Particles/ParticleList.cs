using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Numerics;

namespace SPH_Bachelorprojekt.Simulation.Particles
{
    class ParticleList
    {
        public List<Particle> Particles;
        public int Count;

        public ParticleList(List<Particle> particles)
        {
            Particles = particles;
            Count = particles.Count;
        }

        public ParticleList()
        {
            Particles = new List<Particle>();
            Count = 0;
        }

        public void AddParticleToList(Particle particle)
        {
            Particles.Add(particle);
            Count++;
        }

        public void RemoveParticleFromList(Particle particle)
        {
            Particles.Remove(particle);
            Count--;
        }
    }
}
