﻿// Accord Statistics Library
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
    ///   Moving-window statistics.
    /// </summary>
    /// 
    /// <remarks>
    ///   Provides statistics derived from successive segments of constant, overlapping size ('windowSize') 
    ///   of a series of values. Values are added one at a time to a MovingNormalStatistics instance through 
    ///   <see cref="Push(double)"/> method and are actually kept inside the instance.
    /// </remarks>
    /// 
    /// <example>
    ///   <code source="Unit Tests\Accord.Tests.Statistics\MovingNormalStatisticsTest.cs" region="doc_example" />
    /// </example>
    /// 
    [Serializable]
    public class MovingNormalStatistics : IMovingStatistics, IEnumerable<double>
    {
        private Queue<double> values;
        private Queue<double> squares;

        /// <summary>
        ///   Gets the sum the values within the window.
        /// </summary>
        /// 
        /// <value>The sum of values within the window.</value>
        /// 
        public double Sum { get; private set; }

        /// <summary>
        ///   Gets the sum of squared values within the window.
        /// </summary>
        /// 
        /// <value>The sum of squared values.</value>
        /// 
        public double SumOfSquares { get; private set; }

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
        public int Count { get { return this.values.Count; } }

        /// <summary>
        ///   Gets the mean of the values within the window.
        /// </summary>
        /// 
        /// <value>The mean of the values.</value>
        /// 
        public double Mean { get; private set; }

        /// <summary>
        ///   Gets the variance of the values within the window.
        /// </summary>
        /// 
        /// <value>The variance of the values.</value>
        /// 
        public double Variance { get; private set; }

        /// <summary>
        ///   Gets the standard deviation of the values within the window.
        /// </summary>
        /// 
        /// <value>The standard deviation of the values.</value>
        /// 
        public double StandardDeviation
        {
            get { return System.Math.Sqrt(this.Variance); }
        }



        /// <summary>
        ///   Initializes a new instance of the <see cref="MovingNormalStatistics"/> class.
        /// </summary>
        /// 
        /// <param name="windowSize">The size of the moving window.</param>
        /// 
        public MovingNormalStatistics(int windowSize)
        {
            if (windowSize < 0 || windowSize == int.MaxValue)
                throw new ArgumentOutOfRangeException("windowSize");

            this.Window = windowSize;
            this.values = new Queue<double>(windowSize + 1);
            this.squares = new Queue<double>(windowSize + 1);
        }

        /// <summary>
        ///   Pushes a value into the window.
        /// </summary>
        /// 
        public void Push(double value)
        {
            if (this.values.Count == this.Window)
            {
                this.Sum -= this.values.Dequeue();
                this.SumOfSquares -= this.squares.Dequeue();
            }

            double square = value * value;

            this.values.Enqueue(value);
            this.squares.Enqueue(square);
            
            this.Sum += value;
            this.SumOfSquares += square;

            double N = this.values.Count;

            this.Mean = this.Sum / N;
            this.Variance = (N * this.SumOfSquares - (this.Sum * this.Sum)) / (N * (N - 1));
        }


        /// <summary>
        ///   Returns an enumerator that iterates through the collection.
        /// </summary>
        /// 
        /// <returns>
        ///   A <see cref="T:System.Collections.Generic.IEnumerator`1"/> that can be used to iterate through the collection.
        /// </returns>
        /// 
        public IEnumerator<double> GetEnumerator()
        {
            return this.values.GetEnumerator();
        }


        /// <summary>
        ///   Returns an enumerator that iterates through a collection.
        /// </summary>
        /// 
        /// <returns>
        ///   An <see cref="T:System.Collections.IEnumerator"/> object that can be used to iterate through the collection.
        /// </returns>
        /// 
        System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator()
        {
            return this.values.GetEnumerator();
        }


        /// <summary>
        ///   Removes all elements from the window and resets statistics.
        /// </summary>
        /// 
        public void Clear()
        {
            this.values.Clear();
            this.squares.Clear();

            this.Mean = 0;
            this.Variance = 0;

            this.SumOfSquares = 0;
            this.Sum = 0;
        }
    }
}
