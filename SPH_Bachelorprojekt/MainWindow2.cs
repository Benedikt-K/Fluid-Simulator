﻿using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Numerics;
using SFML.Graphics;
using SFML.Window;
using System.Timers;
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
            public float StartingDensity;
            public float CurrentTimeStep = 0;
            public bool IsPaused = false;
            public float MaxPressure = 0;
            public float TimeStep;
            public float Viscosity;
            public float Stiffness;
            public SimulationLoop SimulationLoop;
            // save screenshot on button
            public int numberOfScreenshots = 0;
            // plotting 
            public int ErrorCollectionPlotCount = 0;
            public ScottPlot.Plot plot = new ScottPlot.Plot();
            List<float> dataX = new List<float>();
            List<float> dataY = new List<float>();
            List<float> DensityIterationErrorDataX = new List<float>();
            List<float> DensityIterationErrorDataY = new List<float>();
            public bool StartDensityErrorCollection = false;
            public bool StartDensityErrorAndIterationCollection = false;
            // what is used for colors
            private readonly bool VelocityColors = true;
            private readonly bool PressureColors = false;
            private readonly bool DensityColors = false;
            // runtime measurement
            public float TimertimestepsCount;
            public float TimerTimeStepCurrentValue;
            public bool TimerIsRunning = false;
            public float TimerStoppingTimeStep = 100f;
            public System.Diagnostics.Stopwatch Stopwatch;
            public float AverageDensityErrorForTimer;
            // save Stimestep screen to folder ?
            public int NumerOfTimeStep = 0;
            public string CurrentDate = DateTime.Now.ToString("dd-MM-yyyy-HH-mm");
            public bool SaveSimulationToImages = false;
            public int SaveEvery_X_TimeStep = 2;
            public float SaveThisNextTimeStep = 0f;
            public float SaveEveryX = 1 / 10;
            // what to use for simulation
            public bool UseIISPH = false;
            public bool UseNeighbour = true;

            public void Run()
            {
                // INITIALIZE IMPORTANT VARIABLES
                float particleSizeH = 1f;                           // works with 8 // IISPH 1
                float viscosity = 2f;                              // works with 2
                float timeStep = 0.01f;                              
                float startDensity = 10f;                          // works with 0.5
                float gravity = -1.8f;                              // works with -0.8
                float smoothingLength = particleSizeH;

                // ONLY FOR SESPH
                float stiffness = 10000f;                             // works with 300  -> größeres k kleinerer TimeStep

                // ONLY FOR VISUALS, scaling
                float scaleFactorDrawing = 10f; 

                // for plotting later
                StartingDensity = startDensity;
                TimeStep = timeStep;
                Viscosity = viscosity;
                Stiffness = stiffness;
                
                // initialize window
                uint videoY = 1200;
                uint videoX = 2000;
                var mode = new VideoMode(videoX, videoY);
                var window = new RenderWindow(mode, "Bachelorthesis - SPH Simulation - Benedikt Kuss");
                var view = window.GetView();
                window.KeyPressed += Window_KeyPressed;
                ParticleSpawner spawner = new ParticleSpawner(startDensity, particleSizeH);

                // Sim
                //List<Particle> particles = spawner.FluidColum();
                List<Particle> particles = spawner.FluidColumOneLayerBoundary(15, 50);
                //List<Particle> particles = spawner.FluidColumOneLayerBoundary(15, 200);
                //List<Particle> particles = spawner.BreakingDamBigAndWideTestLimitOneLayerBoundary();
                //List<Particle> particles = spawner.BreakingDamOneLayerBoundary(500, 500);
                //List<Particle> particles = spawner.BreakingDamOneLayerBoundaryBothSides(600, 200); // from size 3 to 1 --> 3x
                //List<Particle> particles = spawner.WaterfallIntoBoxOneLayerBoundary(500, 500); // 300,300
                //List<Particle> particles = spawner.DroppingFluidDropletOnSurface(600, 250);
                //List<Particle> particles = spawner.BreakingDamOneLayerBoundaryBothSidesAndMiddle(700, 300);
                //List<Particle> particles = spawner.FluidColumWithOutRand();
                //List<Particle> particles = spawner.DroppingFluidColumn();
                //List<Particle> particles = spawner.DroppingFluidColumnBig();
                //List<Particle> particles = spawner.BreakingDam();
                //List<Particle> particles = spawner.BreakingDamOneLayer();
                //List<Particle> particles = spawner.BreakingDamBig();
                //List<Particle> particles = spawner.BreakingDamBigAndWide();
                //List<Particle> particles = spawner.BreakingDamBigAndWideTestLimit();
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
                    ///---> is now done in IISPH and SESPH update formulations respectively
                    ////////////////////////
                    //UPDATE PARTICLES
                    ////////////////////////
                    if (!IsPaused) 
                    {
                        SimulationLoop.UpdateAllParticles(smoothingLength, UseIISPH, UseNeighbour, ref DensityIterationErrorDataX, ref DensityIterationErrorDataY, StartDensityErrorAndIterationCollection);
                        CurrentTimeStep += timeStep;
                    }
                    ///////////////////////
                    //DRAW PARTICLES
                    ///////////////////////
                    if (NumerOfTimeStep % SaveEvery_X_TimeStep == SaveEvery_X_TimeStep - 1 || SaveEvery_X_TimeStep == 0 || !SaveSimulationToImages)
                    {
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
                                    float maxPressure = 0f;
                                    foreach (Particle particle1 in simulationLoop.Particles)
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
                                    // standart color is blue
                                    circle.FillColor = SFML.Graphics.Color.Blue;
                                }
                            }
                            window.Draw(circle);
                        }
                        // display the rendered frame on screen
                        window.Display();
                    }
                    /////////////////////
                    // DATA COLLECTION
                    /////////////////////
                    // density error
                    if (StartDensityErrorCollection) 
                    {
                        dataX.Add(CurrentTimeStep);
                        dataY.Add(Math.Abs((simulationLoop.AverageDensity - StartingDensity) / StartingDensity * 100));
                        //dataY.Add(((simulationLoop.AverageDensity - StartingDensity) / StartingDensity * 100));
                        /*if (CurrentTimeStep >= 70)
                        {
                            plot.Add.Scatter(dataX, dataY);
                            // annotation
                            //plot.Add.Annotation("Viscosity: " + Viscosity + Environment.NewLine + "Stiffness: " + Stiffness + Environment.NewLine + "TimeStep: " + TimeStep);
                            // save to png
                            plot.SavePng("AverageDensityOverTime.png", 800, 600);
                            StartDensityErrorCollection = false;
                            Console.WriteLine("stopped and saved");
                        }*/
                    }
                    // display density error
                    // Console.WriteLine((simulationLoop.AverageDensity - StartingDensity) / StartingDensity * 100);
                    // runtime 
                    if (TimerIsRunning)
                    {
                        if (TimerTimeStepCurrentValue >= TimerStoppingTimeStep)
                        {
                            Stopwatch.Stop();
                            AverageDensityErrorForTimer /= TimertimestepsCount;
                            Console.Write("------------------------------------------------\n");
                            Console.WriteLine("Average Density Error: " + AverageDensityErrorForTimer);
                            Console.WriteLine("Time elapsed for 1 simulated second: " + Stopwatch.ElapsedMilliseconds + "ms");
                            Console.Write("------------------------------------------------\n");
                            TimerIsRunning = false;
                        }
                        else
                        {
                            TimerTimeStepCurrentValue += simulationLoop.TimeStep;
                            TimertimestepsCount++;
                            AverageDensityErrorForTimer += Math.Abs((simulationLoop.AverageDensity - StartingDensity) / StartingDensity * 100);
                        }
                    }
                    /////////////////////
                    //SAVE SIMULATION IMAGES
                    /////////////////////
                    if (SaveSimulationToImages && !IsPaused)
                    {
                        NumerOfTimeStep++;
                        SFML.Graphics.Texture screenshot = new SFML.Graphics.Texture(videoX, videoY);
                        screenshot.Update(window);
                        System.IO.Directory.CreateDirectory(CurrentDate);
                        bool saveCurrentTimeStep = NumerOfTimeStep % SaveEvery_X_TimeStep == 0;
                        //bool saveCurrentTimeStep = simulationLoop.ElapsedTime - SaveThisNextTimeStep > SaveEveryX;
                        if (NumerOfTimeStep != 1 && saveCurrentTimeStep) 
                        {
                            string fileName = NumerOfTimeStep.ToString("D5");
                            screenshot.CopyToImage().SaveToFile(CurrentDate + "/" + "screenshot" + fileName + ".png");
                        }
                        //SaveThisNextTimeStep += SaveEveryX;
                    }
                    if (StartDensityErrorAndIterationCollection)
                    {
                        // wrong density graph
                        //List<float> datapoints = new List<float> {0.543215f, 0.34143f, 0.224321f, 0.201f, 0.1822f, 0.17f, 0.1503f, 0.13f, 0.12f, 0.118f, 0.112f, 0.1032f, 0.102f, 0.097f };
                        //DensityIterationErrorDataY.AddRange(datapoints);
                        //
                        var densityIterErrPlot = plot.Add.Scatter(DensityIterationErrorDataY, DensityIterationErrorDataX);
                        plot.XLabel("Iteration");
                        plot.YLabel("Average density error in %");
                        plot.SavePng("DensityErrorPerIteration" + ErrorCollectionPlotCount + ".png", 800, 600);
                        StartDensityErrorAndIterationCollection = false;
                        plot.Remove(densityIterErrPlot);
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
                    if (StartDensityErrorCollection)
                    {
                        // save data to file 
                        plot.Add.Scatter(dataX, dataY);
                        //plot.Add.HorizontalLine(StartingDensity);
                        // annotation
                        //plot.Add.Annotation("Viscosity: " + Viscosity + Environment.NewLine + "Stiffness: " + Stiffness + Environment.NewLine + "TimeStep: " + TimeStep);
                        // save to png
                        plot.SavePng("AverageDensityOverTime.png", 800, 600);
                        Console.WriteLine("Saving succsessful");
                    }
                    else
                    {
                        CurrentTimeStep = 0;
                        StartDensityErrorCollection = true;
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
                    // start timer 
                    if (!TimerIsRunning)
                    {
                        StartDensityErrorCollection = true;
                        Console.WriteLine("Starting data collection");
                        Stopwatch = System.Diagnostics.Stopwatch.StartNew();
                        TimerTimeStepCurrentValue = 0;
                        TimertimestepsCount = 0;
                        AverageDensityErrorForTimer = 0;
                        TimerIsRunning = true;
                    }
                    else
                    {
                        float elapsedTime = Stopwatch.ElapsedMilliseconds;
                        Stopwatch.Stop();
                        float timePerSecond = elapsedTime / TimerTimeStepCurrentValue;
                        Console.WriteLine("needed computation time" + timePerSecond);
                        plot.Add.Scatter(dataX, dataY);
                        plot.XLabel("Simulated time", size: 16);
                        plot.YLabel("Density error in %", size: 16);
                        // tick density
                        ScottPlot.TickGenerators.NumericAutomatic tickGenX = new ScottPlot.TickGenerators.NumericAutomatic();
                        tickGenX.TickDensity = 0.4;
                        plot.Axes.Bottom.TickGenerator = tickGenX;

                        ScottPlot.TickGenerators.NumericAutomatic tickGenY = new ScottPlot.TickGenerators.NumericAutomatic();
                        tickGenY.TickDensity = 0.4;
                        plot.Axes.Left.TickGenerator = tickGenY;

                        plot.SavePng("AverageDensityOverTime.png", 1000, 500);
                    }
                }
                if (e.Code == SFML.Window.Keyboard.Key.C)
                {
                    StartDensityErrorAndIterationCollection = true;
                    ErrorCollectionPlotCount++;
                }
                if (e.Code == SFML.Window.Keyboard.Key.L)
                {
                    // toggle saving to image
                    SaveSimulationToImages = true;
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
