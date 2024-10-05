using SPH_Bachelorprojekt.Simulation.Particles;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace SPH_Bachelorprojekt.Simulation.Neighbours
{
    internal class NeighbourhoodSearchTest
    {

        public bool TestSpatialHashing()
        {
            int ParticleSize = 2;
            List<Particle> particles = new List<Particle>();
            SpatialHashing spatialHashing = new SpatialHashing(ParticleSize * 2);
            for (int i = 0; i < 5; i++)
            {
                for (int j = 0; j < 5; j++)
                {
                    Particle newParticle = new Particle(new Vector2(i * ParticleSize, j * ParticleSize), 1f, ParticleSize);
                    spatialHashing.InsertObject(newParticle);
                    particles.Add(newParticle);
                }
            }
            foreach (Particle particle in particles)
            {
                spatialHashing.InRadius(particle.Position, ParticleSize * 2f, ref particle.Neighbours);
                //Console.WriteLine("neigbour count: " + particle.Neighbours.Count);
            }
            
            return true;
        }
    }
}
