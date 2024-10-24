using SPH_Bachelorprojekt.Simulation.Kernel_Function;
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
            float ParticleSize = 2;
            List<Particle> particles = new List<Particle>();
            SpatialHashing spatialHashing = new SpatialHashing((int) ParticleSize * 2);
            for (int i = 0; i < 5; i++)
            {
                for (int j = 0; j < 5; j++)
                {
                    Particle newParticle = new Particle(new Vector2(i * ParticleSize, j * ParticleSize), 1f, ParticleSize);
                    spatialHashing.AddParticle(newParticle);
                    particles.Add(newParticle);
                }
            }
            foreach (Particle particle in particles)
            {
                spatialHashing.IsNeighbour(particle.Position, ParticleSize * 2f, ref particle.Neighbours);
                //Console.WriteLine("neighbour count: " + particle.Neighbours.Count);
            }
            
            return true;
        }

        public bool TestQuadratic()
        {
            int ParticleSize = 2;
            List<Particle> particles = new List<Particle>();
            for (int i = 0; i < 5; i++)
            {
                for (int j = 0; j < 5; j++)
                {
                    Particle newParticle = new Particle(new Vector2(i * ParticleSize, j * ParticleSize), 1f, ParticleSize);
                    particles.Add(newParticle);
                }
            }
            Quadratic quadraticNeighbours = new Quadratic(particles, ParticleSize);
            foreach (Particle particle in particles)
            {
                particle.Neighbours = quadraticNeighbours.GetNeighboursQuadratic(particle);
                Console.WriteLine("neighbour count: " + particle.Neighbours.Count);
            }

            return true; 
        }
    }
}
