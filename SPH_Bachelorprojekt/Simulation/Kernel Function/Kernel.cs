using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Numerics;

namespace SPH_Bachelorprojekt.Simulation.Kernel_Function
{
    class Kernel
    {
        private float h; // particle Size
        private float h2; 
        private float alpha; // Normalization factor 2D
        private float beta;
        private float factor;

        private const float kernelCorrection = 0.04f / 0.0400344729f;

        public Kernel(float smoothingLength)
        {
            h = smoothingLength; // works with Kernel support of h
            h2 = h * h;
            alpha = 5 / (14 * Convert.ToSingle(Math.PI) * h2); // For 2D
            beta = 10 / (7 * (float)Math.PI * h2);
            double rad9 = Math.Pow((double)smoothingLength, 9.0);
            factor = (float)(315.0 / (64 * Math.PI * rad9));
        }

        /*public float W(float distance)
        {
            float distance2 = distance * distance;
            if (distance2 > h)
            {
                return 0f;
            }
            if (distance2 < float.Epsilon)
            {
                distance2 = float.Epsilon;
            }
            float diffSq = h2 - distance2;
            return factor * diffSq * diffSq * diffSq;
        }

        public Vector2 GradW(Vector2 position_I, Vector2 position_J)
        {
            Vector2 difference = position_I - position_J;
            float distance2 = Vector2.DistanceSquared(position_I, position_J);
            if (distance2 > h2)
            {
                return Vector2.Zero;
            }
            if (distance2 < float.Epsilon)
            {
                distance2 = float.Epsilon;
            }
            float diffSq = h2 - distance2;
            float f = -factor * 6.0f * diffSq * diffSq;
            return new Vector2(difference.X * f, difference.Y * f);
        }*/

        /*public float W(float distance)
        {
            float q = distance / h;
            float q2 = q * q;
            if (0 <= q && q <= 1)
            {
                return beta * (1 - (3 / 2) * q2 * (1 - (q / 2)));
            }
            else if (1 < q && q <= 2)
            {
                return (beta / 4) * (2 - q) * (2 - q) * (2 - q);
            }
            else
            {
                return 0;
            }
        }*/

        public float W(float distance)
        {
            /// implementation from https://cg.informatik.uni-freiburg.de/course_notes/sim_03_particleFluids.pdf slide 60
            float d = distance / h;
            var t1 = Math.Max(1 - d, 0);
            var t2 = Math.Max(2 - d, 0);
            var t3 = (t2 * t2 * t2) - 4 * (t1 * t1 * t1);
            return alpha * t3 * kernelCorrection;
        }

        
        public Vector2 GradW(Vector2 position_I, Vector2 position_J)
        {
            /// most of implementation from https://cg.informatik.uni-freiburg.de/course_notes/sim_03_particleFluids.pdf slide 64
            Vector2 positionDifference = position_I - position_J;
            float distance = Vector2.Distance(position_I, position_J);
            float d = distance / h;

            if (d == 0)
            {
                return Vector2.Zero;
            }
            float t1 = Math.Max(1 - d, 0);
            float t2 = Math.Max(2 - d, 0);
            float t3 = (-3 * t2 * t2) + (12 * t1 * t1);
            return alpha * (positionDifference / (distance * h)) * t3 ; 
        }

        public void TestKernel()
        {
            // Tests from https://cg.informatik.uni-freiburg.de/course_notes/sim_03_particleFluids.pdf
            // Test Kernel - W
            Console.WriteLine("Tests for Kernel:");
            float pi = Convert.ToSingle(Math.PI);

            float test1 = W(0 * h);
            float solution1 = 20 / (14 * pi * h2);

            float test2 = W(1 * h);
            float solution2 = 5 / (14 * pi * h2);

            float test3 = W(Convert.ToSingle(Math.Sqrt(2)) * h);
            float solution3 = 1.005f / (14f * pi * h2);

            float test4 = test1 + 4 * test2 + 4 * test3;
            float solution4 = 1.001f / h2;

            Console.WriteLine(test1 + " " + solution1);
            if (test1 == solution1) { Console.WriteLine("Test 1 Passed"); }
            Console.WriteLine("------------------");
            Console.WriteLine(test2 + " " + solution2);
            if (test2 == solution2) { Console.WriteLine("Test 2 Passed"); }
            Console.WriteLine("------------------");
            Console.WriteLine(test3 + " " + solution3);
            if (Math.Abs(test3 - solution3) < 0.01f) { Console.WriteLine("Test 3 Passed"); }
            Console.WriteLine("------------------");
            Console.WriteLine(test4 + " " + solution4);
            if (Math.Abs(test4 - solution4) < 0.01f) { Console.WriteLine("Test 4 Passed"); }
            Console.WriteLine("------------------");

            // Test Kernel Gradient - GradW
            float beta = -3f * Convert.ToSingle(Math.Pow(2 - Math.Sqrt(2), 2));
            float factor = (1 / (h * Convert.ToSingle(Math.Sqrt(2))));
            Console.WriteLine("Tests for Gradient of Kernel:");
            Console.WriteLine("------------------");

            Vector2 test5 = GradW(Vector2.Zero, Vector2.Zero);
            Vector2 solution5 = new Vector2(0, 0);
            Console.WriteLine("test 5: " + test5 + " I " + solution5);
            Console.WriteLine("------------------");

            Vector2 test6 = GradW(Vector2.Zero, new Vector2(h, 0));
            Vector2 solution6 = new Vector2((3 * alpha) / h, 0);
            Console.WriteLine("test 6: " + test6 + " I " + solution6);
            Console.WriteLine("------------------");

            Vector2 test7 = GradW(Vector2.Zero, new Vector2(0, h));
            Vector2 solution7 = new Vector2(0, (3 * alpha) / h);
            Console.WriteLine("test 7: " + test7 + " I " + solution7);
            Console.WriteLine("------------------");

            Vector2 test8 = GradW(Vector2.Zero, new Vector2(h, h));
            Vector2 solution8 = new Vector2(-factor * alpha * beta, -factor * alpha * beta);
            Console.WriteLine("test 8: " + test8 + " I " + solution8);
            Console.WriteLine("------------------");

            Vector2 test9 = GradW(Vector2.Zero, new Vector2(h, -h));
            Vector2 solution9 = new Vector2(-factor * alpha * beta, factor * alpha * beta);
            Console.WriteLine("test 9: " + test9 + " I " + solution9);
            Console.WriteLine("------------------");
        }
    }
}
