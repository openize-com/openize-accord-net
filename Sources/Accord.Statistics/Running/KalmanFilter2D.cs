// Accord Statistics Library
// The Accord.NET Framework
// http://accord-framework.net
//
// Copyright © Pablo Guzman Sanchez, 2013
// pablogsanchez at gmail.com
//
// Copyright © César Souza, 2009-2017
// cesarsouza at gmail.com
//
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
//
// This code originated as a contribution by Pablo Sanches, originally based on
// Student Dave's tutorial on Object Tracking in Images Using 2D Kalman Filters,
// shared under the LGPL by explicit written permissions from both authors:
//
//   http://studentdavestutorials.weebly.com/object-tracking-2d-kalman-filter.html
//

namespace FileFormat.Accord.Statistics.Running
{
    using System;
    using FileFormat.Accord.Core.AForge.Core;
    using FileFormat.Accord.Core.Exceptions;
    using global::Accord.Math;
    using Math.Core;
    using Math.Matrix;

    /// <summary>
    ///   Kalman filter for 2D coordinate systems.
    /// </summary>
    /// 
    /// <remarks>
    /// <para>
    ///   References:
    ///   <list type="bullet">
    ///     <item><description><a href="http://studentdavestutorials.weebly.com/object-tracking-2d-kalman-filter.html">
    ///       Student Dave's tutorial on Object Tracking in Images Using 2D Kalman Filters.
    ///       Available on: http://studentdavestutorials.weebly.com/object-tracking-2d-kalman-filter.html
    ///       </a></description></item>
    ///   </list></para>
    /// </remarks>
    /// 
    /// <example>
    /// <code source="Unit Tests\Accord.Tests.Statistics\KalmanFilterTest.cs" region="doc_push" />
    /// </example>
    /// 
    [Serializable]
    public class KalmanFilter2D : IRunning<DoublePoint>, IRunning<double[]>
    {

        double samplingRate = 1;

        double acceleration = 0.0005f;
        double accelStdDev = 0.1f;

        double[,] Q_estimate; // (location_0, location_1, vel_0, vel_1)

        double[,] A;
        double[,] B;
        double[,] C;

        double[,] Ez;
        double[,] Ex;
        double[,] P;
        double[,] K;
        double[,] Aux;

        static readonly double[,] diagonal =
        {
            { 1, 0, 0, 0 },
            { 0, 1, 0, 0 },
            { 0, 0, 1, 0 },
            { 0, 0, 0, 1 }
        };

        /// <summary>
        ///   Gets or sets the current X position of the object.
        /// </summary>
        /// 
        public double X
        {
            get { return this.Q_estimate[0, 0]; }
            set { this.Q_estimate[0, 0] = value; }
        }

        /// <summary>
        ///   Gets or sets the current Y position of the object.
        /// </summary>
        /// 
        public double Y
        {
            get { return this.Q_estimate[1, 0]; }
            set { this.Q_estimate[1, 0] = value; }
        }

        /// <summary>
        ///   Gets or sets the current object's velocity in the X axis.
        /// </summary>
        /// 
        public double XAxisVelocity
        {
            get { return this.Q_estimate[2, 0]; }
            set { this.Q_estimate[2, 0] = value; }
        }

        /// <summary>
        ///   Gets or sets the current object's velocity in the Y axis.
        /// </summary>
        /// 
        public double YAxisVelocity
        {
            get { return this.Q_estimate[3, 0]; }
            set { this.Q_estimate[3, 0] = value; }
        }

        /// <summary>
        ///   Gets or sets the observational noise 
        ///   of the current object's in the X axis.
        /// </summary>
        /// 
        public double NoiseX
        {
            get { return this.Ez[0, 0]; }
            set { this.Ez[0, 0] = value; }
        }

        /// <summary>
        ///   Gets or sets the observational noise 
        ///   of the current object's in the Y axis.
        /// </summary>
        /// 
        public double NoiseY
        {
            get { return this.Ez[1, 1]; }
            set { this.Ez[1, 1] = value; }
        }

        /// <summary>
        ///   Initializes a new instance of the <see cref="KalmanFilter2D"/> class.
        /// </summary>
        /// 
        public KalmanFilter2D()
        {
            this.initialize();
        }

        /// <summary>
        ///   Initializes a new instance of the <see cref="KalmanFilter2D"/> class.
        /// </summary>
        /// 
        /// <param name="samplingRate">The sampling rate.</param>
        /// <param name="acceleration">The acceleration.</param>
        /// <param name="accelerationStdDev">The acceleration standard deviation.</param>
        /// 
        public KalmanFilter2D(double samplingRate, double acceleration, double accelerationStdDev)
        {
            this.acceleration = acceleration;
            this.accelStdDev = accelerationStdDev;
            this.samplingRate = samplingRate;

            this.initialize();
        }

        private void initialize()
        {
            double dt = this.samplingRate;

            this.A = new double[,]
            {
                { 1,  0, dt,  0 },
                { 0,  1,  0, dt },
                { 0,  0,  1,  0 },
                { 0,  0,  0,  1 }
            };

            this.B = new double[,]
            {
                { (dt * dt) / 2 },
                { (dt * dt) / 2 },
                {       dt      },
                {       dt      }
            };

            this.C = new double[,]
            {
                { 1, 0, 0, 0 },
                { 0, 1, 0, 0 }
            };

            this.Ez = new double[,] 
            {
                { 1.0, 0.0 }, 
                { 0.0, 1.0 }
            };

            double dt2 = dt * dt;
            double dt3 = dt2 * dt;
            double dt4 = dt2 * dt2;

            double aVar = this.accelStdDev * this.accelStdDev;

            this.Ex = new double[4, 4]
            {
                { dt4 / 4,        0,  dt3 / 2,        0 },
                { 0,        dt4 / 4,        0,  dt3 / 2 },
                { dt3 / 2,        0,      dt2,        0 },
                { 0,        dt3 / 2,        0,      dt2 }
            };

            this.Ex.Multiply(aVar, result: this.Ex);

            this.Q_estimate = new double[4, 1];
            this.P = this.Ex.MemberwiseClone();
        }


        /// <summary>
        ///   Registers the occurrence of a value.
        /// </summary>
        /// 
        /// <param name="value">The value to be registered.</param>
        /// 
        public void Push(double[] value)
        {
            if (value.Length != 2)
                throw new DimensionMismatchException("value");

            this.Push(value[0], value[1]);
        }

        /// <summary>
        ///   Registers the occurrence of a value.
        /// </summary>
        /// 
        /// <param name="value">The value to be registered.</param>
        /// 
        public void Push(DoublePoint value)
        {
            this.Push(value.X, value.Y);
        }

        /// <summary>
        ///   Registers the occurrence of a value.
        /// </summary>
        /// 
        /// <param name="x">The x-coordinate of the value to be registered.</param>
        /// <param name="y">The y-coordinate of the value to be registered.</param>
        /// 
        public void Push(double x, double y)
        {
            double[,] Qloc = { { x }, { y } };

            // Predict next state
            this.Q_estimate = Matrix.Dot(this.A, this.Q_estimate).Add(this.B.Multiply(this.acceleration));

            // Predict Covariances
            this.P = Matrix.Dot(this.A, this.P.DotWithTransposed(this.A)).Add(this.Ex);

            this.Aux = Matrix.Dot(this.C, this.P.DotWithTransposed(this.C)).Add(this.Ez).PseudoInverse();

            // Kalman Gain
            this.K = this.P.Dot(this.C.TransposeAndDot(this.Aux));
            this.Q_estimate = this.Q_estimate.Add(this.K.Dot(Qloc.Subtract(this.C.Dot(this.Q_estimate))));

            // Update P (Covariances)
            this.P = Matrix.Dot(diagonal.Subtract(Matrix.Dot(this.K, this.C)), this.P);
        }


        /// <summary>
        ///   Clears all measures previously computed.
        /// </summary>
        /// 
        public void Clear()
        {
            this.NoiseX = 0;
            this.NoiseY = 0;

            this.XAxisVelocity = 0;
            this.YAxisVelocity = 0;

            this.X = 0;
            this.Y = 0;
        }
    }
}
