using System;
using System.Numerics;
using System.Collections.Generic;
using SPH_Bachelorprojekt.Simulation.Particles;
using SPH_Bachelorprojekt.Simulation.MainSimulation;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SPH_Bachelorprojekt.Utilities.ParticleUtils
{
    class ParticleSpawner
    {
        public List<Vector2> Positions;
        public float Density;
        public float ParticleSizeH;
        public Vector2 LeftBottomOfDomain;
        public Vector2 RightTopOfDomain;

        public ParticleSpawner(List<Vector2> positions, float density, float particleSizeH)
        {
            Positions = positions;
            Density = density;
            ParticleSizeH = particleSizeH;
        }

        public ParticleSpawner(float density, float particleSizeH)
        {
            Density = density;
            ParticleSizeH = particleSizeH;
        }

        public List<Particle> SpawnParticles()
        {
            List<Particle> particleList = new List<Particle>();
            for (int i = 0; i < Positions.Count; i++)
            {
                Particle newParticle = new Particle(Positions[i], Density, ParticleSizeH);
                particleList.Add(newParticle);
            }
            return particleList;
        }

        public float GetRandomValue(float maxMagnitude)
        {
            Random rand = new Random();
            float randomNumber = Convert.ToSingle(rand.NextDouble());
            return randomNumber * maxMagnitude;
        }


        /// 
        /// INITIALIZATION METHODS
        /// 
        public List<Particle> TestSetup()
        {
            List<Particle> particleList = new List<Particle>();
            for (float i = 0; i < ParticleSizeH * 50; i += ParticleSizeH)
            {
                for (float j = 0; j < ParticleSizeH * 50; j += ParticleSizeH)
                {
                    if (j == ParticleSizeH * 49)
                    {
                        Particle newParticle = new Particle(new Vector2(i, j), Density, ParticleSizeH, true);
                        particleList.Add(newParticle);
                    }
                    else
                    {
                        Particle newParticle = new Particle(new Vector2(i, j), Density, ParticleSizeH, false);
                        particleList.Add(newParticle);
                    }

                }
            }
            return particleList;
        }

        public List<Particle> TestOneParticleWithTwoLayerBoundarys()
        {
            List<Particle> particleList = new List<Particle>();

            for (float i = 0; i < ParticleSizeH * 20; i += ParticleSizeH)
            {
                Particle Boundarys = new Particle(new Vector2(i, 20f * ParticleSizeH), Density, ParticleSizeH, true);
                particleList.Add(Boundarys);
            }
            for (float i = 0; i < ParticleSizeH * 20; i += ParticleSizeH)
            {
                Particle Boundarys = new Particle(new Vector2(i, 20f * ParticleSizeH + ParticleSizeH), Density, ParticleSizeH, true);
                particleList.Add(Boundarys);
            }
            Particle droppedOne = new Particle(new Vector2(ParticleSizeH * 9, 26f * ParticleSizeH), Density, ParticleSizeH, false);
            particleList.Add(droppedOne);

            return particleList;
        }

        public List<Particle> TestOneParticleWithTwoLayerBoundarysOnOther()
        {
            List<Particle> particleList = new List<Particle>();
            float magitudeOfDeviation = 0.1f;

            for (float i = 0; i < ParticleSizeH * 20; i += ParticleSizeH)
            {
                Particle Boundarys = new Particle(new Vector2(i, 20f * ParticleSizeH), Density, ParticleSizeH, true);
                particleList.Add(Boundarys);
            }
            for (float i = 0; i < ParticleSizeH * 20; i += ParticleSizeH)
            {
                Particle Boundarys = new Particle(new Vector2(i, 20f * ParticleSizeH + ParticleSizeH), Density, ParticleSizeH, true);
                particleList.Add(Boundarys);
            }
            Particle droppedOne = new Particle(new Vector2(ParticleSizeH * 9f + GetRandomValue(magitudeOfDeviation), 22f * ParticleSizeH), Density, ParticleSizeH, false);
            particleList.Add(droppedOne);
            Particle droppedOne2 = new Particle(new Vector2(ParticleSizeH * 9f + GetRandomValue(magitudeOfDeviation), 24f * ParticleSizeH), Density, ParticleSizeH, false);
            particleList.Add(droppedOne2);

            return particleList;
        }

        public List<Particle> TestOneParticleWithOneLayerBoundarys()
        {
            List<Particle> particleList = new List<Particle>();

            for (float i = 0; i < ParticleSizeH * 20; i += ParticleSizeH)
            {
                Particle Boundarys = new Particle(new Vector2(i, 20f * ParticleSizeH), Density, ParticleSizeH, true);
                particleList.Add(Boundarys);
            }
            Particle droppedOne = new Particle(new Vector2(ParticleSizeH * 9, 24f * ParticleSizeH), Density, ParticleSizeH, false);
            particleList.Add(droppedOne);

            return particleList;
        }

        public List<Particle> TestOneParticleOnBoundarys()
        {
            List<Particle> particleList = new List<Particle>();

            for (float i = 0; i < ParticleSizeH * 20; i += ParticleSizeH)
            {
                Particle Boundarys = new Particle(new Vector2(i, 20f * ParticleSizeH), Density, ParticleSizeH, true);
                particleList.Add(Boundarys);
            }
            for (float i = 0; i < ParticleSizeH * 20; i += ParticleSizeH)
            {
                Particle Boundarys = new Particle(new Vector2(i, 20f * ParticleSizeH + ParticleSizeH), Density, ParticleSizeH, true);
                particleList.Add(Boundarys);
            }
            Particle layingOne = new Particle(new Vector2(ParticleSizeH * 9, 20f * ParticleSizeH - ParticleSizeH), Density, ParticleSizeH, false);
            particleList.Add(layingOne);

            return particleList;
        }


        public List<Particle> TestParticlesInJarDroppingTwoLayerBoundary()
        {
            List<Particle> particleList = new List<Particle>();

            // Spawn boundarys
            for (int i = 5; i < 25; i += 1)
            {
                for (int j = 5; j < 45; j += 1)
                {
                    if (i < 7 || i > 22 || j < 7 || j > 42)
                    {
                        Particle particle = new Particle(new Vector2(i * ParticleSizeH, j * ParticleSizeH), Density, ParticleSizeH, true);
                        particleList.Add(particle);
                    }
                }
            }

            // spawn fluid
            for (int i = 10; i < 16; i++)
            {
                for (int j = 30; j < 40; j++)
                {
                    Particle particle1 = new Particle(new Vector2(i * ParticleSizeH, j * ParticleSizeH), Density, ParticleSizeH, false);
                    particleList.Add(particle1);
                }
            }


            return particleList;
        }


        public List<Particle> TestParticlesInJarTwoLayerBoundary()
        {
            List<Particle> particleList = new List<Particle>();
            float magitudeOfDeviation = 0f;

            // Spawn boundarys
            int maxI = 85;
            int maxJ = 45;
            for (int i = 5; i < maxI; i += 1)
            {
                for (int j = 5; j < maxJ; j += 1)
                {
                    if (i < 7 || i > maxI - 3 || j < 7 || j > maxJ - 3)
                    {
                        Particle particle = new Particle(new Vector2(i * ParticleSizeH, j * ParticleSizeH), Density, ParticleSizeH, true);
                        particleList.Add(particle);
                    }
                }
            }

            // spawn fluid
            for (int i = 10; i < 16; i++)
            {
                for (int j = 8; j < 18; j++)
                {
                    /*if (i % 2 == 0 && j % 2 == 0)
                    {
                        Particle particle1 = new Particle(new Vector2(i * ParticleSizeH, j * ParticleSizeH), Density, ParticleSizeH, false);
                        particleList.Add(particle1);
                    }*/
                    Particle particle1 = new Particle(new Vector2(i * ParticleSizeH + GetRandomValue(magitudeOfDeviation), j * ParticleSizeH + GetRandomValue(magitudeOfDeviation)), Density, ParticleSizeH, false);
                    particleList.Add(particle1);
                }
            }
            return particleList;
        }

        public List<Particle> BreakingDam()
        {
            List<Particle> particleList = new List<Particle>();
            float magitudeOfDeviation = 0.1f;

            // Spawn boundarys
            int maxI = 85;
            int maxJ = 45;
            for (int i = 5; i < maxI; i += 1)
            {
                for (int j = 5; j < maxJ; j += 1)
                {
                    if (i < 7 || i > maxI - 3 || j < 7 || j > maxJ - 3)
                    {
                        Particle particle = new Particle(new Vector2(i * ParticleSizeH, j * ParticleSizeH), Density, ParticleSizeH, true);
                        particleList.Add(particle);
                    }
                }
            }

            // spawn fluid
            for (int i = 7; i < 16; i++)
            {
                for (int j = 7; j < 18; j++)
                {
                    Particle particle1 = new Particle(new Vector2(i * ParticleSizeH + GetRandomValue(magitudeOfDeviation), j * ParticleSizeH + GetRandomValue(magitudeOfDeviation)), Density, ParticleSizeH, false);
                    particleList.Add(particle1);
                }
            }
            return particleList;
        }

        public List<Particle> BreakingDamOneLayer()
        {
            List<Particle> particleList = new List<Particle>();
            float magitudeOfDeviation = 0.1f;

            // Spawn boundarys
            int maxI = 85;
            int maxJ = 45;
            for (int i = 5; i < maxI; i += 1)
            {
                for (int j = 5; j < maxJ; j += 1)
                {
                    if (i < 6 || i > maxI - 2 || j < 6 || j > maxJ - 2)
                    {
                        Particle particle = new Particle(new Vector2(i * ParticleSizeH, j * ParticleSizeH), Density, ParticleSizeH, true);
                        particleList.Add(particle);
                    }
                }
            }

            // spawn fluid
            for (int i = 6; i < 16; i++)
            {
                for (int j = 6; j < 18; j++)
                {
                    Particle particle1 = new Particle(new Vector2(i * ParticleSizeH + GetRandomValue(magitudeOfDeviation), j * ParticleSizeH + GetRandomValue(magitudeOfDeviation)), Density, ParticleSizeH, false);
                    particleList.Add(particle1);
                }
            }
            return particleList;
        }

        public List<Particle> BreakingDamBig()
        {
            List<Particle> particleList = new List<Particle>();
            float magitudeOfDeviation = 0.1f;

            // Spawn boundarys
            int maxI = 55;
            int maxJ = 85;
            for (int i = 5; i < maxI; i += 1)
            {
                for (int j = 5; j < maxJ; j += 1)
                {
                    if (i < 7 || i > maxI - 3 || j < 7 || j > maxJ - 3)
                    {
                        Particle particle = new Particle(new Vector2(i * ParticleSizeH, j * ParticleSizeH), Density, ParticleSizeH, true);
                        particleList.Add(particle);
                    }
                }
            }

            // spawn fluid
            for (int i = 7; i < maxI / 3; i++)
            {
                for (int j = 7; j < maxJ / 1.5; j++)
                {
                    Particle particle1 = new Particle(new Vector2(i * ParticleSizeH + GetRandomValue(magitudeOfDeviation), j * ParticleSizeH + GetRandomValue(magitudeOfDeviation)), Density, ParticleSizeH, false);
                    particleList.Add(particle1);
                }
            }
            return particleList;
        }

        public List<Particle> BreakingDamBigAndWideTestLimit()
        {
            List<Particle> particleList = new List<Particle>();
            float magitudeOfDeviation = 0f;

            // Spawn boundarys
            int maxI = 125;
            int maxJ = 85;
            for (int i = 5; i < maxI; i += 1)
            {
                for (int j = 5; j < maxJ; j += 1)
                {
                    if (i < 7 || i > maxI - 3 || j < 7)
                    {
                        Particle particle = new Particle(new Vector2(i * ParticleSizeH, j * ParticleSizeH), Density, ParticleSizeH, true);
                        particleList.Add(particle);
                    }
                }
            }

            // spawn fluid
            for (int i = 7; i < maxI / 3; i++)
            {
                for (int j = 7; j < maxJ / 1.5; j++)
                {
                    Particle particle1 = new Particle(new Vector2(i * ParticleSizeH + GetRandomValue(magitudeOfDeviation), j * ParticleSizeH + GetRandomValue(magitudeOfDeviation)), Density, ParticleSizeH, false);
                    particleList.Add(particle1);
                }
            }
            return particleList;
        }

        public List<Particle> BreakingDamBigAndWideTestLimitOneLayerBoundary()
        {
            List<Particle> particleList = new List<Particle>();
            float magitudeOfDeviation = 0f;

            int maxI = 125;
            int maxJ = 85;

            int fluidStart = maxI / 3;

            // Spawn boundarys
            for (int i = 6; i < maxI; i += 1)
            {
                for (int j = 6; j < maxJ; j += 1)
                {
                    if (i < 7 || i > maxI - 2 || j < 7)
                    {
                        Particle particle = new Particle(new Vector2(i * ParticleSizeH, j * ParticleSizeH), Density, ParticleSizeH, true);
                        particleList.Add(particle);
                    }
                }
            }

            // spawn removeable boundary
            for (int j = 7; j < maxJ; j += 1)
            {
                Particle particle = new Particle(new Vector2(fluidStart * ParticleSizeH, j * ParticleSizeH), Density, ParticleSizeH, true, true);
                particleList.Add(particle);
                
            }
            

            // spawn fluid
            for (int i = 7; i < fluidStart; i++)
            {
                for (int j = 7; j < maxJ / 1.5; j++)
                {
                    Particle particle1 = new Particle(new Vector2(i * ParticleSizeH + GetRandomValue(magitudeOfDeviation), j * ParticleSizeH + GetRandomValue(magitudeOfDeviation)), Density, ParticleSizeH, false);
                    particleList.Add(particle1);
                }
            }
            return particleList;
        }

        public List<Particle> FunnelIntoBox()
        {
            List<Particle> particleList = new List<Particle>();
            float magitudeOfDeviation = 0.1f;

            // Spawn boundarys
            int maxI = 60;
            int maxJ = 80;
            int funnelDownIteration = 0;
            for (int i = 5; i < maxI; i += 1)
            {
                int FunnelIsPlaced = 0;
                for (int j = 5; j < maxJ; j += 1)
                {
                    // outer Box
                    if (i < 7 || i > maxI - 3 || j < 7 || j > maxJ - 3)
                    {
                        Particle particle = new Particle(new Vector2(i * ParticleSizeH, j * ParticleSizeH), Density, ParticleSizeH, true);
                        particleList.Add(particle);
                    }
                    // inner funnel
                    else if (j > maxJ / 2)
                    {
                        if (i <= maxI / 2 && FunnelIsPlaced < 2)
                        {
                            Particle particle = new Particle(new Vector2(i * ParticleSizeH, j * ParticleSizeH - funnelDownIteration), Density, ParticleSizeH, true);
                            particleList.Add(particle);
                            Particle particle2 = new Particle(new Vector2((maxI + 4 - i) * ParticleSizeH, j * ParticleSizeH - funnelDownIteration), Density, ParticleSizeH, true);
                            particleList.Add(particle2);
                            FunnelIsPlaced++;
                            funnelDownIteration++;
                        }
                    }
                }
            }

            // spawn fluid
            for (int i = 8; i < maxI - 3; i++)
            {
                for (int j = maxJ / 2 + 3; j < maxJ - 8; j++)
                {
                    Particle particle1 = new Particle(new Vector2(i * ParticleSizeH + GetRandomValue(magitudeOfDeviation), j * ParticleSizeH + GetRandomValue(magitudeOfDeviation)), Density, ParticleSizeH, false);
                    particleList.Add(particle1);
                }
            }
            return particleList;
        }

        public List<Particle> BreakingDamBigWithHole()
        {
            List<Particle> particleList = new List<Particle>();
            float magitudeOfDeviation = 0.1f;

            // Spawn boundarys
            int maxI = 55;
            int maxJ = 45;
            for (int i = 5; i < maxI; i += 1)
            {
                for (int j = 5; j < maxJ; j += 1)
                {
                    if (i < 7 || i > maxI - 3 || j < 7 || j > maxJ - 3)
                    {
                        if (i == 30 || i == 31)
                        {
                            continue;
                        }
                        else
                        {
                            Particle particle = new Particle(new Vector2(i * ParticleSizeH, j * ParticleSizeH), Density, ParticleSizeH, true);
                            particleList.Add(particle);
                        }
                    }
                }
            }

            // spawn fluid
            for (int i = 7; i < 26; i++)
            {
                for (int j = 7; j < 38; j++)
                {
                    Particle particle1 = new Particle(new Vector2(i * ParticleSizeH + GetRandomValue(magitudeOfDeviation), j * ParticleSizeH + GetRandomValue(magitudeOfDeviation)), Density, ParticleSizeH, false);
                    particleList.Add(particle1);
                }
            }
            return particleList;
        }

        

        public List<Particle> FluidColum()
        {
            List<Particle> particleList = new List<Particle>();
            float magitudeOfDeviation = 0.1f;
            
            // Spawn boundarys
            int maxI = 15;
            int maxJ = 120;
            for (int i = 5; i < maxI; i += 1)
            {
                for (int j = 5; j < maxJ; j += 1)
                {
                    if (i < 7 || i > maxI - 3 || j < 7 || j > maxJ - 3)
                    {
                        Particle particle = new Particle(new Vector2(i * ParticleSizeH, j * ParticleSizeH), Density, ParticleSizeH, true);
                        particleList.Add(particle);
                    }
                }
            }

            //STANDART
            //int maxI = 13;
            //int maxJ = 68;
            // spawn fluid
            for (int i = 7; i < 13; i++)
            {
                for (int j = 7; j < 100; j++)
                {
                    Particle particle1 = new Particle(new Vector2(i * ParticleSizeH + GetRandomValue(magitudeOfDeviation), j * ParticleSizeH), Density, ParticleSizeH, false);
                    particleList.Add(particle1);
                }
            }
            return particleList;
        }

        public List<Particle> FluidColumOneLayerBoundary()
        {
            List<Particle> particleList = new List<Particle>();
            float magitudeOfDeviation = 0.1f;

            // Spawn boundarys
            int maxI = 15;
            int maxJ = 50;
            for (int i = 6; i < maxI; i += 1)
            {
                for (int j = 6; j < maxJ; j += 1)
                {
                    if (i < 7 || i > maxI - 2 || j < 7 || j > maxJ - 2)
                    {
                        Particle particle = new Particle(new Vector2(i * ParticleSizeH, j * ParticleSizeH), Density, ParticleSizeH, true);
                        particleList.Add(particle);
                    }
                }
            }

            //STANDART
            //int maxI = 13;
            //int maxJ = 68;
            // spawn fluid
            for (int i = 7; i < maxI - 1; i++)
            {
                for (int j = 7; j < maxJ - 5; j++)
                {
                    Particle particle1 = new Particle(new Vector2(i * ParticleSizeH + GetRandomValue(magitudeOfDeviation), j * ParticleSizeH), Density, ParticleSizeH, false);
                    particleList.Add(particle1);
                }
            }
            return particleList;
        }

        public List<Particle> FluidColumWithOutRand()
        {
            List<Particle> particleList = new List<Particle>();
            float magitudeOfDeviation = 0.1f;

            // Spawn boundarys
            int maxI = 15;
            int maxJ = 45;
            for (int i = 5; i < maxI; i += 1)
            {
                for (int j = 5; j < maxJ; j += 1)
                {
                    if (i < 7 || i > maxI - 3 || j < 7 || j > maxJ - 3)
                    {
                        Particle particle = new Particle(new Vector2(i * ParticleSizeH, j * ParticleSizeH), Density, ParticleSizeH, true);
                        particleList.Add(particle);
                    }
                }
            }

            // spawn fluid
            for (int i = 7; i < 13; i++)
            {
                for (int j = 7; j < 18; j++)
                {
                    Particle particle1 = new Particle(new Vector2(i * ParticleSizeH, j * ParticleSizeH), Density, ParticleSizeH, false);
                    particleList.Add(particle1);
                }
            }
            return particleList;
        }

        public List<Particle> BreakingDamOnBothSides()
        {
            List<Particle> particleList = new List<Particle>();
            float magitudeOfDeviation = 0f;

            // Spawn boundarys
            int maxI = 105;
            int maxJ = 65;
            for (int i = 5; i < maxI; i += 1)
            {
                for (int j = 5; j < maxJ; j += 1)
                {
                    if (i < 7 || i > maxI - 3 || j < 7)
                    {
                        Particle particle = new Particle(new Vector2(i * ParticleSizeH, j * ParticleSizeH), Density, ParticleSizeH, true);
                        particleList.Add(particle);
                    }
                }
            }

            // spawn fluid LEFT
            for (int i = 7; i < 38; i++)
            {
                for (int j = 7; j < 38; j++)
                {
                    Particle particle1 = new Particle(new Vector2(i * ParticleSizeH + GetRandomValue(magitudeOfDeviation), j * ParticleSizeH + GetRandomValue(magitudeOfDeviation)), Density, ParticleSizeH, false);
                    particleList.Add(particle1);
                }
            }

            // spawn fluid RIGHT
            for (int i = maxI - 34; i < maxI - 2; i++)
            {
                for (int j = 7; j < 38; j++)
                {
                    Particle particle1 = new Particle(new Vector2(i * ParticleSizeH + GetRandomValue(magitudeOfDeviation), j * ParticleSizeH + GetRandomValue(magitudeOfDeviation)), Density, ParticleSizeH, false);
                    particleList.Add(particle1);
                }
            }
            return particleList;
        }

        public List<Particle> DroppingFluidColumn()
        {
            List<Particle> particleList = new List<Particle>();
            float magitudeOfDeviation = 2f;

            // Spawn boundarys
            int maxI = 65;
            int maxJ = 45;
            for (int i = 5; i < maxI; i += 1)
            {
                for (int j = 5; j < maxJ; j += 1)
                {
                    if (i < 7 || i > maxI - 3 || j < 7 || j > maxJ - 3)
                    {
                        Particle particle = new Particle(new Vector2(i * ParticleSizeH, j * ParticleSizeH), Density, ParticleSizeH, true);
                        particleList.Add(particle);
                    }
                }
            }

            // spawn fluid
            for (int i = 15; i < 40; i++)
            {
                for (int j = 15; j < 40; j++)
                {
                    Particle particle1 = new Particle(new Vector2(i * ParticleSizeH + GetRandomValue(magitudeOfDeviation), j * ParticleSizeH + GetRandomValue(magitudeOfDeviation)), Density, ParticleSizeH, false);
                    particleList.Add(particle1);
                }
            }
            return particleList;
        }

        public List<Particle> DroppingFluidColumnBig()
        {
            List<Particle> particleList = new List<Particle>();
            float magitudeOfDeviation = 0.1f;

            // Spawn boundarys
            int maxI = 85;
            int maxJ = 65;
            for (int i = 5; i < maxI; i += 1)
            {
                for (int j = 5; j < maxJ; j += 1)
                {
                    if (i < 7 || i > maxI - 3 || j < 7)
                    {
                        Particle particle = new Particle(new Vector2(i * ParticleSizeH, j * ParticleSizeH), Density, ParticleSizeH, true);
                        particleList.Add(particle);
                    }
                }
            }

            // spawn fluid
            for (int i = 7; i < 50; i++)
            {
                for (int j = 7; j < 50; j++)
                {
                    Particle particle1 = new Particle(new Vector2(i * ParticleSizeH + GetRandomValue(magitudeOfDeviation), j * ParticleSizeH + GetRandomValue(magitudeOfDeviation)), Density, ParticleSizeH, false);
                    particleList.Add(particle1);
                }
            }
            return particleList;
        }

        public List<Particle> TestParticlesInJarTwoLayerBoundaryWithSpacing()
        {
            List<Particle> particleList = new List<Particle>();

            // Spawn boundarys
            int maxI = 85;
            int maxJ = 45;
            for (int i = 5; i < maxI; i += 1)
            {
                for (int j = 5; j < maxJ; j += 1)
                {
                    if (i < 7 || i > maxI - 3 || j < 7 || j > maxJ - 3)
                    {
                        Particle particle = new Particle(new Vector2(i * ParticleSizeH, j * ParticleSizeH), Density, ParticleSizeH, true);
                        particleList.Add(particle);
                    }
                }
            }

            // spawn fluid
            for (int i = 10; i < 16; i++)
            {
                for (int j = 8; j < 18; j++)
                {
                    /*if (i % 2 == 0 && j % 2 == 0)
                    {
                        Particle particle1 = new Particle(new Vector2(i * ParticleSizeH, j * ParticleSizeH), Density, ParticleSizeH, false);
                        particleList.Add(particle1);
                    }*/
                    Particle particle1 = new Particle(new Vector2(i * ParticleSizeH * 1.2f, j * ParticleSizeH * 1.2f), Density, ParticleSizeH, false);
                    particleList.Add(particle1);
                }
            }


            return particleList;
        }

        public List<Particle> TestParticlesBottomOfJar()
        {
            List<Particle> particleList = new List<Particle>();

            for (float i = 5 * ParticleSizeH; i < ParticleSizeH * 105; i += ParticleSizeH)
            {
                for (float j = 5 * ParticleSizeH; j < ParticleSizeH * 105; j += ParticleSizeH)
                {
                    if (i == 5 * ParticleSizeH || i == ParticleSizeH * 104 || j == ParticleSizeH * 104 || i == 6 * ParticleSizeH || i == ParticleSizeH * 103 || j == ParticleSizeH * 103)
                    {
                        Particle particle = new Particle(new Vector2(i, j), Density, ParticleSizeH, true);
                        particleList.Add(particle);
                    }
                    else if (j > 54)
                    {
                        Particle particle = new Particle(new Vector2(i, j), Density, ParticleSizeH, false);
                        particleList.Add(particle);
                    }
                }
            }

            return particleList;
        }



        /// Other spawning functions
        ///
        public Particle SpawnParticleAtPositionWithVelocity(Vector2 position, Vector2 velocity)
        {
            Particle particle = new Particle(position, velocity, Density, ParticleSizeH, false);
            return particle;
        }
    }
}
