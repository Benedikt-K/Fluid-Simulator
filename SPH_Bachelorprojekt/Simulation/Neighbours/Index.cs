using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using SPH_Bachelorprojekt.Simulation.Particles;
using System.Numerics;

namespace SPH_Bachelorprojekt.Simulation.Kernel_Function
{
    class Index
    {
        public List<Particle> Particles;
        public float ParticleSizeH;
        public List<Particle>[] Grid;
        public float GridSize;
        public int NumberOfGridCells;
        public float MaxDomainX;
        public float MaxDomainY;
        public float MinDomainX;
        public float MinDomainY;

        public Index(List<Particle> particles, float particleSizeH, float maxDomainX, float maxDomainY, float minDomainX, float minDomainY)
        {
            Particles = particles;
            ParticleSizeH = particleSizeH;
            GridSize = 2 * particleSizeH;
            MaxDomainX = maxDomainX;
            MaxDomainY = maxDomainY;
            MinDomainX = minDomainX;
            MinDomainY = minDomainY;
            //NumberOfGridCells = Math.Ceiling(Convert.ToDouble(maxDomainX - minDomainX * maxDomainY - minDomainY));
        }

        public int ComputeCellIndex(Particle particle)
        // compute cell index, returns int index if successful, returns -1 if not
        {
            int cellIndex = -1;
            float k = (particle.Position.X - MinDomainX) / GridSize;
            float l = (particle.Position.Y - MinDomainY) / GridSize;

            cellIndex = (int) (k + l * GridSize);

            return cellIndex;
        }

        public void SortParticlesInGrid()
        {
            foreach (Particle particle in Particles)
            {
                int cellIndex = ComputeCellIndex(particle);
                if (Grid[cellIndex] == null)
                {
                    Grid[cellIndex] = new List<Particle> { particle };
                }
                else
                {
                    Grid[cellIndex].Add(particle);
                }
            }

        }

        public List<Particle> GetNeighbours(Particle particle)
        {
            List<Particle> neighbours = new List<Particle>();



            return neighbours;
        }

        public void InitializeGrid()
        {
            float TestGridNumber = 50 * ParticleSizeH / 2;

            List<Particle>[] grid = new List<Particle>[(int) (TestGridNumber * TestGridNumber)];
            Grid = grid;
            GridSize = (int) TestGridNumber;
        }
    }
}
