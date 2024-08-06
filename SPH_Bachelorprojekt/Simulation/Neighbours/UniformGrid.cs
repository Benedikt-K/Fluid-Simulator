using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using SPH_Bachelorprojekt.Simulation.Particles;
using System.Numerics;

namespace SPH_Bachelorprojekt.Simulation.Kernel_Function
{
    class UniformGrid
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

        public UniformGrid(List<Particle> particles, float particleSizeH, float maxDomainX, float maxDomainY, float minDomainX, float minDomainY)
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
            float k = 10 * (particle.Position.X - MinDomainX) / GridSize;
            float l = 10 * (particle.Position.Y - MinDomainY) / GridSize;

            cellIndex = (int) (k + l * GridSize);

            return cellIndex;
        }

        public void SortParticlesInGrid()
        {
            Parallel.ForEach(Particles, particle =>
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
            });

        }

        public List<Particle> GetNeighbours(Particle particle)
        {
            List<Particle> neighbours = new List<Particle>();
            float k = (particle.Position.X - MinDomainX) / GridSize;
            float l = (particle.Position.Y - MinDomainY) / GridSize;

            for (int i = -1; i < 2; i++)
            {
                for (int j = -1; j < 2; j++) 
                {

                    int cellIndex = (int) ((k + i) + (l + j) * GridSize);
                    foreach (var potentialNeighbour in neighbours)
                    {
                        if (Vector2.Distance(particle.Position, potentialNeighbour.Position) <= 2f * ParticleSizeH)
                        {
                            neighbours.Add(potentialNeighbour);
                        }
                    }
                }
            }

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
