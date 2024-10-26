using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Numerics;
using SFML.Graphics;
using SFML.Window;
using SPH_Bachelorprojekt.Simulation.Particles;
using SPH_Bachelorprojekt.Utilities.ParticleUtils;
using SPH_Bachelorprojekt.Simulation.Kernel_Function;
using SPH_Bachelorprojekt.Simulation.MainSimulation;

namespace SPH_Bachelorprojekt
{
    class MainWindow2
    {
        static void Main(string[] args)
        {
            var window = new SimpleWindow();
            window.Run();
        }
        class SimpleWindow
        {
            public ScottPlot.Plot plot = new ScottPlot.Plot();
            List<float> dataX = new List<float>();
            List<float> dataY = new List<float>();
            public bool StartDataCollection = false;
            public int numberOfScreenshots = 0;
            public float StartingDensity;
            public float CurrentTimeStep = 0;
            public bool IsPaused = false;
            public float MaxPressure = 0;
            public float TimeStep;
            public float Viscosity;
            public float Stiffness;
            public SimulationLoop SimulationLoop;
            // what is used for colors
            private readonly bool VelocityColors = true;
            private readonly bool PressureColors = false;
            private readonly bool DensityColors = false;
            // what to use for simulation
            public bool UseIISPH = true;
            public bool UseNeighbour = true;

            public void Run()
            {
                // INITIALIZE IMPORTANT VARIABLES
                float particleSizeH = 5f;                           // works with 8
                float viscosity = 20f;                              // works with 10
                float timeStep = 0.1f;                              // works with 0.2
                float startDensity = 0.3f;                          // works with 0.3
                float gravity = -0.8f;                              // works with -0.4
                float smoothingLength = particleSizeH;

                // ONLY FOR SESPH
                float stiffness = 300f;                             // works with 300  -> größeres k kleinerer TimeStep

                // ONLY FOR VISUALS, scaling
                float scaleFactorDrawing = 2f;

                // for plotting later
                StartingDensity = startDensity;
                TimeStep = timeStep;
                Viscosity = viscosity;
                Stiffness = stiffness;
                
                // initialize window
                uint videoY = 1000;
                var mode = new VideoMode(1800, videoY);
                var window = new RenderWindow(mode, "Bachelorthesis - SPH Simulation - Benedikt Kuss");
                var view = window.GetView();
                window.KeyPressed += Window_KeyPressed;
                ParticleSpawner spawner = new ParticleSpawner(startDensity, particleSizeH);

                // Sim
                //List<Particle> particles = spawner.FluidColum();
                List<Particle> particles = spawner.FluidColumOneLayerBoundary();
                //List<Particle> particles = spawner.FluidColumHighOneLayerBoundary();
                //List<Particle> particles = spawner.FluidColumWithOutRand();
                //List<Particle> particles = spawner.DroppingFluidColumn();
                //List<Particle> particles = spawner.DroppingFluidColumnBig();
                //List<Particle> particles = spawner.BreakingDam();
                //List<Particle> particles = spawner.BreakingDamOneLayer();
                //List<Particle> particles = spawner.BreakingDamBig();
                //List<Particle> particles = spawner.BreakingDamBigAndWide();
                //List<Particle> particles = spawner.BreakingDamBigAndWideTestLimit();
                //List<Particle> particles = spawner.BreakingDamBigAndWideTestLimitOneLayerBoundary();
                //List<Particle> particles = spawner.FunnelIntoBox();
                //List<Particle> particles = spawner.BreakingDamBigWithHole();
                //List<Particle> particles = spawner.BreakingDamOnBothSides();
                // Sim

                // TESTS
                //List<Particle> particles = spawner.TestOneParticleWithTwoLayerBoundarys();
                //List<Particle> particles = spawner.TestOneParticleWithTwoLayerBoundarysOnOther();
                //List<Particle> particles = spawner.TestParticlesInJarTwoLayerBoundaryWithSpacing();
                //List<Particle> particles = spawner.TestParticlesInJarDroppingTwoLayerBoundary();
                //List<Particle> particles = spawner.TestParticlesInJarTwoLayerBoundary();
                //List<Particle> particles = spawner.TestSetup();
                // TESTS

                Console.WriteLine("Particles spawned");
                Quadratic quadraticSolver = new Quadratic(particles, particleSizeH);
                //IndexSort indexSortSolver = new IndexSort(particles, particleSizeH);
                List<Particle> neighbours = new List<Particle>();
                SimulationLoop simulationLoop = new SimulationLoop(particles, startDensity, particleSizeH, timeStep, viscosity, stiffness, gravity);
                SimulationLoop = simulationLoop;
                var circle = new SFML.Graphics.CircleShape((particleSizeH / 2f) * scaleFactorDrawing)
                {
                    FillColor = SFML.Graphics.Color.White
                };


                // RENDER and SIMULATION LOOP
                while (window.IsOpen)
                {
                    // Process events
                    window.DispatchEvents();
                    ////////////////////////
                    //UPDATE TimeStep
                    ////////////////////////
                    if (SimulationLoop.CalculateParticleLambdaCFL(SimulationLoop.MaxVelocity) > 0.5f)
                    {
                        timeStep *= 0.7f;
                    }
                    else if (SimulationLoop.CalculateParticleLambdaCFL(SimulationLoop.MaxVelocity) > 0.2f)
                    {
                        timeStep *= 1.5f;
                    }
                    ////////////////////////
                    //UPDATE PARTICLES
                    ////////////////////////
                    if (!IsPaused) 
                    {
                        SimulationLoop.UpdateAllParticles(smoothingLength, UseIISPH, UseNeighbour);
                    }
                    ///////////////////////
                    //DRAW PARTICLES
                    ///////////////////////
                    window.Clear();
                    MaxPressure = SimulationLoop.MaxCurrentParticlePressure;
                    foreach (var particle in SimulationLoop.Particles)
                    {
                        SFML.System.Vector2f transformedPosition = new SFML.System.Vector2f((particle.Position.X - particleSizeH / 2) * scaleFactorDrawing, videoY - (particle.Position.Y - particleSizeH / 2) * scaleFactorDrawing);
                        circle.Position = new SFML.System.Vector2f(transformedPosition.X, transformedPosition.Y);
                        if (particle.IsBoundaryParticle)
                        {
                            if (!particle.IsRemoveable)
                            {
                                circle.FillColor = SFML.Graphics.Color.White;
                            }
                            else
                            {
                                circle.FillColor = SFML.Graphics.Color.Green;
                            }                            
                        }
                        else
                        {
                            // color coded
                            if (VelocityColors)
                            {
                                // Color for velocity
                                float lambdaParticle = SimulationLoop.CalculateParticleLambdaCFL(particle.Velocity * 2);
                                circle.FillColor = GetColor(lambdaParticle);
                            }
                            else if (PressureColors) 
                            {
                                // Color for pressure
                                if (CurrentTimeStep < timeStep * 2)
                                {
                                    circle.FillColor = SFML.Graphics.Color.Blue;
                                }
                                if (!particle.IsBoundaryParticle)
                                {
                                    //Console.WriteLine(particle.Pressure);
                                }
                                //float pressureFactor = particle.Pressure  / 20f;
                                //float pressureFactor = particle.Pressure / MaxPressure;
                                float maxPressure = 0f;
                                foreach(Particle particle1 in simulationLoop.Particles)
                                {
                                    if (particle1.Pressure < maxPressure)
                                    {
                                        maxPressure = particle1.Pressure;
                                    }
                                }
                                float pressureFactor = particle.Pressure / maxPressure + 0.01f;
                                circle.FillColor = GetColor(pressureFactor);
                            }
                            else if (DensityColors)
                            {
                                // Color for density
                                float densityFactor = (particle.Density - 0.99f * startDensity) / (1.01f * startDensity - 0.99f * startDensity);
                                circle.FillColor = GetColor(densityFactor);
                            }
                            else
                            {
                                circle.FillColor = SFML.Graphics.Color.Blue;
                            }
                        }
                        window.Draw(circle);
                    }
                    // Finally, display the rendered frame on screen
                    window.Display();

                    /////////////////////
                    // DATA COLLECTION
                    /////////////////////
                    // add data to Collection
                    if (StartDataCollection) 
                    {
                        CurrentTimeStep += timeStep;
                        dataX.Add(CurrentTimeStep);
                        dataY.Add((simulationLoop.AverageDensity - StartingDensity) / StartingDensity * 100);
                        if (CurrentTimeStep >= 70)
                        {
                            plot.Add.Scatter(dataX, dataY);
                            // annotation
                            //plot.Add.Annotation("Viscosity: " + Viscosity + Environment.NewLine + "Stiffness: " + Stiffness + Environment.NewLine + "TimeStep: " + TimeStep);
                            // save to png
                            plot.SavePng("AverageDensityOverTime.png", 800, 600);
                            StartDataCollection = false;
                            Console.WriteLine("stopped and saved");
                        }
                    }
                }
            }


            

            private void Window_KeyPressed(object sender, SFML.Window.KeyEventArgs e)
            {
                var window = (SFML.Window.Window)sender;
                if (e.Code == SFML.Window.Keyboard.Key.Escape)
                {
                    window.Close();
                }
                if (e.Code == SFML.Window.Keyboard.Key.P)
                {
                    if (IsPaused)
                    {
                        IsPaused = false;
                    }
                    else
                    {
                        IsPaused = true;
                    }
                }
                if (e.Code == SFML.Window.Keyboard.Key.S)
                {
                    if (StartDataCollection)
                    {
                        // save data to file 
                        plot.Add.Scatter(dataX, dataY);
                        //plot.Add.HorizontalLine(StartingDensity);
                        // annotation
                        plot.Add.Annotation("Viscosity: " + Viscosity + Environment.NewLine + "Stiffness: " + Stiffness + Environment.NewLine + "TimeStep: " + TimeStep);
                        // save to png
                        plot.SavePng("AverageDensityOverTime.png", 800, 600);
                        Console.WriteLine("Saving succsessful");
                    }
                    else
                    {
                        CurrentTimeStep = 0;
                        StartDataCollection = true;
                        Console.WriteLine("Starting data collection");
                    }
                }
                if (e.Code == SFML.Window.Keyboard.Key.I)
                {
                    numberOfScreenshots++;
                    SFML.Graphics.Texture texture = new SFML.Graphics.Texture(window.Size.X, window.Size.Y);
                    texture.Update(window);
                    texture.CopyToImage().SaveToFile("screenshot" + numberOfScreenshots + ".png");
                    Console.WriteLine("Saving succsessful");
                }
                if (e.Code == SFML.Window.Keyboard.Key.X)
                {
                    List<Particle> newP = new List<Particle>();
                    foreach (Particle particle in SimulationLoop.Particles)
                    {
                        if (!particle.IsRemoveable)
                        {
                            newP.Add(particle);
                        }
                    }
                    SimulationLoop.Particles = newP;
                }
            }

            private SFML.Graphics.Color GetColor(float lambdaParticle)
            {
                // make sure value is between 0 and 1
                lambdaParticle = Math.Max(0, Math.Min(1, lambdaParticle));

                if (lambdaParticle < 0.5f)
                {
                    // Interpolate between blue and green
                    float t = lambdaParticle / 0.5f;
                    byte r = (byte)(0 * (1 - t) + 0 * t);
                    byte g = (byte)(0 * (1 - t) + 255 * t);
                    byte b = (byte)(255 * (1 - t) + 0 * t);
                    return new SFML.Graphics.Color(r, g, b);
                }
                else
                {
                    // Interpolate between green and red
                    float t = (lambdaParticle - 0.5f) / 0.5f;
                    byte r = (byte)(0 * (1 - t) + 255 * t);
                    byte g = (byte)(255 * (1 - t) + 0 * t);
                    byte b = (byte)(0 * (1 - t) + 0 * t);
                    return new SFML.Graphics.Color(r, g, b);
                }
            }

        }
    }
}
