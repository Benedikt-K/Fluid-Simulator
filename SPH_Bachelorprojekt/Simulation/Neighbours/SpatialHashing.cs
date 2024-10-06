using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Numerics;
using SPH_Bachelorprojekt.Simulation.Particles;
using SPH_Bachelorprojekt.Simulation.Kernel_Function;

namespace SPH_Bachelorprojekt.Simulation.Neighbours
{
    public class SpatialHashing
    {
        public int Count { get; private set; }
        public int CellSize;
        private readonly Dictionary<Vector2, HashSet<Particle>> mSpatialGrids = new Dictionary<Vector2, HashSet<Particle>>();

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

        public void InsertObject(Particle particle)
        {
            var hash = Hash(particle.Position);
            if (!mSpatialGrids.TryGetValue(hash, out var objectBucket))
            {
                objectBucket = new HashSet<Particle>();
                mSpatialGrids[hash] = objectBucket;
            }
            objectBucket.Add(particle);
            Count++;
        }

        public void RemoveObject(Particle particle)
        {
            var hash = Hash(particle.Position);
            if (!mSpatialGrids.TryGetValue(hash, out var objectBucket)) return;
            if (!objectBucket.Remove(particle)) return;
            Count--;
            mSpatialGrids[hash] = objectBucket;
            if (objectBucket.Count == 0) mSpatialGrids.Remove(hash);
        }

        public void Clear() => mSpatialGrids.Clear();

        public void InRadius(Vector2 position, float radius, ref List<Particle> particleInRadius)
        {
            particleInRadius.Clear();
            var startX = (int)Math.Floor((position.X - radius) / CellSize);
            var endX = (int)Math.Ceiling((position.X + radius) / CellSize);
            var startY = (int)Math.Floor((position.Y - radius) / CellSize);
            var endY = (int)Math.Ceiling((position.Y + radius) / CellSize);
            var xRange = Enumerable.Range(startX, endX - startX + 1);
            var yRange = Enumerable.Range(startY, endY - startY + 1);

            foreach (var x in xRange)
            {
                foreach (var y in yRange)
                {
                    var hash = new Vector2(x, y);
                    if (!mSpatialGrids.ContainsKey(hash)) continue;
                    var objectsInBucket = mSpatialGrids[hash];
                    foreach (Particle particle in objectsInBucket)
                    {
                        var distance = Vector2.Distance(particle.Position, position);
                        if (distance > radius)
                            continue;
                        particleInRadius.Add(particle);
                    }
                }
            }
        }

        
    }
}
