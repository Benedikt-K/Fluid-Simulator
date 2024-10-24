using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using SPH_Bachelorprojekt.Simulation.Particles;

namespace SPH_Bachelorprojekt.Simulation.Neighbours
{
    public class SpatialHashing
    {
        public int Count;
        public int CellSize;
        private readonly Dictionary<Vector2, HashSet<Particle>> SpatialGrid = new Dictionary<Vector2, HashSet<Particle>>();

        public SpatialHashing(int cellSize)
        {
            CellSize = cellSize;
        }

        public Vector2 Hash(Vector2 vector)
        {
            vector = Vector2.Divide(vector, CellSize);
            Vector2 vectorFloored = new Vector2((float) Math.Floor(vector.X), (float) Math.Floor(vector.Y));
            return vectorFloored;
        }

        public void AddParticle(Particle particle)
        {
            Vector2 hash = Hash(particle.Position);
            if (!SpatialGrid.TryGetValue(hash, out var objectBucket))
            {
                objectBucket = new HashSet<Particle>();
                SpatialGrid[hash] = objectBucket;
            }
            objectBucket.Add(particle);
            Count++;
        }

        public void RemoveParticle(Particle particle)
        {
            var hash = Hash(particle.Position);
            if (!SpatialGrid.TryGetValue(hash, out HashSet<Particle> objectBucket))
            {
                return;
            }
            if (!objectBucket.Remove(particle))
            {
                return;
            }
            Count--;
            SpatialGrid[hash] = objectBucket;
            if (objectBucket.Count == 0)
            {
                SpatialGrid.Remove(hash);
            }
        }

        public void Clear()
        {
            SpatialGrid.Clear();
        }

        public void GetNeighbours(Vector2 position, float radius, ref List<Particle> particleInSearchRadius)
        {
            particleInSearchRadius.Clear();
            int startPos_X = (int) Math.Floor((position.X - radius) / CellSize);
            int endPos_X = (int) Math.Ceiling((position.X + radius) / CellSize);
            int startPos_Y = (int) Math.Floor((position.Y - radius) / CellSize);
            int endPos_Y = (int) Math.Ceiling((position.Y + radius) / CellSize);
            IEnumerable<int> xRange = Enumerable.Range(startPos_X, endPos_X - startPos_X + 1);
            IEnumerable<int> yRange = Enumerable.Range(startPos_Y, endPos_Y - startPos_Y + 1);

            foreach (int x in xRange)
            {
                foreach (int y in yRange)
                {
                    var hash = new Vector2(x, y);
                    if (!SpatialGrid.ContainsKey(hash))
                    {
                        continue;
                    }
                    HashSet<Particle> objectsInBucket = SpatialGrid[hash];
                    foreach (Particle particle in objectsInBucket)
                    {
                        float distance = Vector2.Distance(particle.Position, position);
                        if (distance > radius)
                        {
                            continue;
                        }
                        particleInSearchRadius.Add(particle);
                    }
                }
            }
        }

        
    }
}
