// Accord Statistics Library
// The Accord.NET Framework
// http://accord-framework.net
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

namespace Openize.Accord.Statistics.Moving
{
    using System;
    using System.Collections.Generic;

    /// <summary>
    ///   Moving-window circular statistics.
    /// </summary>
    /// 
    [Serializable]
    public class MovingCircularStatistics : IMovingStatistics
    {

        private Queue<double> sines;
        private Queue<double> cosines;

        /// <summary>
        ///   Gets the sum of the sines of the angles within the window.
        /// </summary>
        /// 
        public double SumOfSines { get; private set; }

        /// <summary>
        ///   Gets the sum of the cosines of the angles within the window.
        /// </summary>
        /// 
        public double SumOfCosines { get; private set; }


        /// <summary>
        ///   Gets the size of the window.
        /// </summary>
        /// 
        /// <value>The window's size.</value>
        /// 
        public int Window { get; private set; }

        /// <summary>
        ///   Gets the number of samples within the window.
        /// </summary>
        /// 
        /// <value>The number of samples within the window.</value>
        /// 
        public int Count { get { return this.sines.Count; } }

        /// <summary>
        ///   Gets the mean of the angles within the window.
        /// </summary>
        /// 
        /// <value>The mean.</value>
        /// 
        public double Mean { get; private set; }

        /// <summary>
        ///   Gets the variance of the angles within the window.
        /// </summary>
        /// 
        public double Variance { get; private set; }

        /// <summary>
        ///   Gets the standard deviation of the angles within the window.
        /// </summary>
        /// 
        public double StandardDeviation { get; private set; }

        /// <summary>
        /// Gets the current length of the sample mean resultant vector of the gathered values.
        /// </summary>
        /// 
        public double Rho { get; private set; }

         /// <summary>
        ///   Initializes a new instance of the <see cref="MovingCircularStatistics"/> class.
        /// </summary>
        /// 
        /// <param name="windowSize">The size of the moving window.</param>
        /// 
        public MovingCircularStatistics(int windowSize)
        {
            if (windowSize < 0 || windowSize == int.MaxValue)
                throw new ArgumentOutOfRangeException("windowSize");

            this.Window = windowSize;
            this.sines = new Queue<double>(windowSize + 1);
            this.cosines = new Queue<double>(windowSize + 1);
        }

        /// <summary>
        ///   Registers the occurrence of a value.
        /// </summary>
        /// 
        /// <param name="value">The value to be registered.</param>
        /// 
        public void Push(double value)
        {
            if (this.sines.Count == this.Window)
            {
                this.SumOfSines -= this.sines.Dequeue();
                this.SumOfCosines -= this.cosines.Dequeue();
            }

            double cos = Math.Cos(value);
            double sin = Math.Sin(value);

            this.sines.Enqueue(sin);
            this.cosines.Enqueue(cos);

            this.SumOfSines += sin;
            this.SumOfCosines += cos;

            double N = this.sines.Count;

            this.Rho = Math.Sqrt(this.SumOfSines * this.SumOfSines + this.SumOfCosines * this.SumOfCosines);
            this.Mean = Math.Atan2(this.SumOfSines / N, this.SumOfCosines / N);
            this.Variance = Math.Max(0, 1.0 - this.Rho / N);
            this.StandardDeviation = Math.Sqrt(-2.0 * Math.Log(this.Rho / N));
        }

        /// <summary>
        ///   Clears all measures previously computed.
        /// </summary>
        /// 
        public void Clear()
        {
            this.sines.Clear();
            this.cosines.Clear();

            this.Mean = 0;
            this.Variance = 0;
            this.StandardDeviation = 0;
            this.Rho = 0;

            this.SumOfSines = 0;
            this.SumOfCosines = 0;
        }


    }
}
