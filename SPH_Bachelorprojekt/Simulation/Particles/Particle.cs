using System;
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
        public Vector2 PositionNew;
        public Vector2 Velocity;
        public Vector2 PredictedVelocity;
        public Vector2 D_i_i;
        public float A_i_i;
        public Vector2 Dij_Pj;
        public Vector2 Acceleration;
        public Vector2 NonPressureForces;
        public float Density;
        public float Mass;
        public float ParticleSizeH;
        public bool IsBoundaryParticle;
        public List<Particle> Neighbours;
        public float Pressure;
        public float PredictedPressure;

        public Particle(Vector2 position, float density, float particleSizeH, bool boundaryParticle=false)
        {
            Position = position;
            PositionNew = position;
            Velocity = Vector2.Zero;
            PredictedVelocity = Vector2.Zero;
            Acceleration = Vector2.Zero;
            NonPressureForces = Vector2.Zero;
            D_i_i = Vector2.Zero;
            A_i_i = 0.0f;
            Dij_Pj = Vector2.Zero;
            Density = density;
            ParticleSizeH = particleSizeH;
            IsBoundaryParticle = boundaryParticle;
            // initialize mass
            Mass = density * ParticleSizeH * ParticleSizeH;
            Neighbours = new List<Particle>();
            Pressure = 1f;
            PredictedPressure = 0f;
        }

        public Particle(Vector2 position, Vector2 velocity, float density, float particleSizeH, bool boundaryParticle = false)
        {
            Position = position;
            Velocity = velocity;
            PredictedVelocity = Vector2.Zero;
            Density = density;
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

    }
}
