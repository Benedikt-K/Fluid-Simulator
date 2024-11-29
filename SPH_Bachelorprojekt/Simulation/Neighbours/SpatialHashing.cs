using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using SPH_Bachelorprojekt.Simulation.Particles;

namespace SPH_Bachelorprojekt.Simulation.Neighbours
{
    public class SpatialHashing
    {
        //
        // partially based on https://stackoverflow.com/questions/39580691/how-to-find-neighbor-squares-in-spatial-hashing
        // and https://cg.informatik.uni-freiburg.de/course_notes/sim_03_particleFluids.pdf
        //

        public int Count;
        public int CellSize;
        public Dictionary<Vector2, HashSet<Particle>> DictionarySpatial = new Dictionary<Vector2, HashSet<Particle>>();

        public SpatialHashing(int cellSize)
        {
            CellSize = cellSize;
        }

        public Vector2 HashFunction(Vector2 vector)
        {
            vector = Vector2.Divide(vector, CellSize);
            Vector2 vectorFloored = new Vector2((float) Math.Floor(vector.X), (float) Math.Floor(vector.Y));
            return vectorFloored;
        }

        public void AddParticle(Particle particle)
        {
            Vector2 hash = HashFunction(particle.Position);
            if (!DictionarySpatial.TryGetValue(hash, out HashSet<Particle> newHashSet))
            {
                newHashSet = new HashSet<Particle>();
                DictionarySpatial[hash] = newHashSet;
            }
            newHashSet.Add(particle);
            Count = Count + 1;
        }

        public void RemoveParticle(Particle particle)
        {
            Vector2 hash = HashFunction(particle.Position);
            if (!DictionarySpatial.TryGetValue(hash, out HashSet<Particle> newHashSet))
            {
                return;
            }
            if (!newHashSet.Remove(particle))
            {
                return;
            }
            Count = Count - 1;
            DictionarySpatial[hash] = newHashSet;
            if (newHashSet.Count == 0)
            {
                DictionarySpatial.Remove(hash);
            }
        }

        public void GetNeighbours(Vector2 position, float searchRadius, ref List<Particle> foundNeighbours)
        {
            foundNeighbours.Clear();
            int X_StartPosition = (int) Math.Floor((position.X - searchRadius) / CellSize);
            int X_EndPosition = (int) Math.Ceiling((position.X + searchRadius) / CellSize);
            int Y_StartPosition = (int) Math.Floor((position.Y - searchRadius) / CellSize);
            int Y_EndPosition = (int) Math.Ceiling((position.Y + searchRadius) / CellSize);
            IEnumerable<int> X_PossibleCells = Enumerable.Range(X_StartPosition, X_EndPosition - X_StartPosition + 1);
            IEnumerable<int> Y_PossibleCells = Enumerable.Range(Y_StartPosition, Y_EndPosition - Y_StartPosition + 1);

            foreach (int x in X_PossibleCells)
            {
                foreach (int y in Y_PossibleCells)
                {
                    Vector2 hashValue = new Vector2(x, y);
                    if (DictionarySpatial.ContainsKey(hashValue))
                    {
                        HashSet<Particle> currentCellHasSet = DictionarySpatial[hashValue];
                        foreach (Particle particle in currentCellHasSet)
                        {
                            float distance = Vector2.Distance(particle.Position, position);
                            if (distance <= searchRadius)
                            {
                                foundNeighbours.Add(particle);
                            }
                        }
                    }
                }
            }
        }

        
    }
}
