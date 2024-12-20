﻿using System;
using System.Numerics;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;


namespace SPH_Bachelorprojekt.Simulation.Particles
{
    public class Particle
    {
        public Vector2 Position;
        public Vector2 Velocity;
        public Vector2 PredictedVelocity;
        public float DiagonalElement;
        public Vector2 Acceleration;
        public Vector2 PressureAcceleration;
        public Vector2 NonPressureAcceleration;
        public float Density;
        public float Mass;
        public float ParticleSizeH;
        public bool IsBoundaryParticle;
        public bool IsRemoveable;
        public List<Particle> Neighbours;
        public float Pressure;
        public float SourceTerm;

        public Particle(Vector2 position, float density, float particleSizeH, bool boundaryParticle=false, bool isRemoveable=false)
        {
            Position = position;
            Velocity = Vector2.Zero;
            PredictedVelocity = Vector2.Zero;
            Acceleration = Vector2.Zero;
            PressureAcceleration = Vector2.Zero;
            NonPressureAcceleration = Vector2.Zero;
            DiagonalElement = 0;
            Density = density;
            ParticleSizeH = particleSizeH;
            IsBoundaryParticle = boundaryParticle;
            IsRemoveable = isRemoveable;
            // initialize mass
            Mass = density * ParticleSizeH * ParticleSizeH;
            Neighbours = new List<Particle>();
            Pressure = 0f;
            SourceTerm = 0;
        }

        public Particle(Vector2 position, Vector2 velocity, float density, float particleSizeH, bool boundaryParticle = false)
        {
            Position = position;
            Velocity = velocity;
            PredictedVelocity = Vector2.Zero;
            Density = density;
            Pressure = 0f;
            ParticleSizeH = particleSizeH;
            IsBoundaryParticle = boundaryParticle;
            // initialize mass
            Mass = density * ParticleSizeH * ParticleSizeH;
            Neighbours = new List<Particle>();
        }

        public Particle(Vector2 position, Vector2 velocity, float density, float mass, float particleSizeH, bool boundaryParticle = false)
        {
            Position = position;
            Velocity = velocity;
            PredictedVelocity = Vector2.Zero;
            Density = density;
            ParticleSizeH = particleSizeH;
            IsBoundaryParticle = boundaryParticle;
            Mass = mass;
            Neighbours = new List<Particle>();
        }

        public float GetVolume()
        {
            float volume = Mass / Density;
            return volume;
        }

        public float GetMass()
        {
            float mass = Density * ParticleSizeH * ParticleSizeH;
            return Mass;
        }

    }
}
