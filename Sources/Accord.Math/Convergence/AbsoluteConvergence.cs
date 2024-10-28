// Accord Math Library
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

namespace FileFormat.Accord.Math.Convergence
{
    using System;
    using Base;

    /// <summary>
    ///   Absolute convergence criteria.
    /// </summary>
    /// 
    /// <remarks>
    ///   This class can be used to track progress and convergence
    ///   of methods which rely on the absolute change of a value.
    /// </remarks>
    /// 
    /// <example>
    /// <code>
    ///   // Create a new convergence criteria for a maximum of 10 iterations
    ///   var criteria = new AbsoluteConvergence(iterations: 10, tolerance: 0.1);
    /// 
    ///   int progress = 1;
    /// 
    ///   do
    ///   {
    ///       // Do some processing...
    /// 
    /// 
    ///       // Update current iteration information:
    ///       criteria.NewValue = 12345.6 / progress++;
    /// 
    ///   } while (!criteria.HasConverged);
    /// 
    ///   
    ///   // The method will converge after reaching the 
    ///   // maximum of 10 iterations with a final value
    ///   // of 1371.73:
    ///   
    ///   int iterations = criteria.CurrentIteration; // 10
    ///   double value = criteria.OldValue; // 1371.7333333
    /// </code>
    /// </example>
    /// 
    public class AbsoluteConvergence : ISingleValueConvergence
    {

        private double tolerance = 0;
        private int maxIterations = 100;
        private double newValue;


        /// <summary>
        ///   Gets or sets the maximum change in the watched value
        ///   after an iteration of the algorithm used to detect 
        ///   convergence. Default is 0.
        /// </summary>
        /// 
        public double Tolerance
        {
            get { return this.tolerance; }
            set
            {
                if (value < 0)
                    throw new ArgumentOutOfRangeException("value", "Tolerance should be positive.");

                this.tolerance = value;
            }
        }

        /// <summary>
        ///   Please use MaxIterations instead.
        /// </summary>
        /// 
        [Obsolete("Please use MaxIterations instead.")]
        public int Iterations
        {
            get { return this.MaxIterations; }
            set { this.MaxIterations = value; }
        }

        /// <summary>
        ///   Gets or sets the maximum number of iterations
        ///   performed by the iterative algorithm. Default 
        ///   is 100.
        /// </summary>
        /// 
        public int MaxIterations
        {
            get { return this.maxIterations; }
            set
            {
                if (value < 0)
                    throw new ArgumentOutOfRangeException("value",
                        "The maximum number of iterations should be positive.");

                this.maxIterations = value;
            }
        }

        /// <summary>
        ///   Initializes a new instance of the <see cref="AbsoluteConvergence"/> class.
        /// </summary>
        /// 
        public AbsoluteConvergence()
        {
        }

        /// <summary>
        ///   Initializes a new instance of the <see cref="AbsoluteConvergence"/> class.
        /// </summary>
        /// 
        /// <param name="iterations">The maximum number of iterations which should be
        ///   performed by the iterative algorithm. Setting to zero indicates there
        ///   is no maximum number of iterations. Default is 0.</param>
        /// <param name="tolerance">The maximum change in the watched value
        ///   after an iteration of the algorithm used to detect convergence.
        ///   Default is 0.</param>
        /// <param name="startValue">The initial value for the <see cref="NewValue"/> and
        ///   <see cref="OldValue"/> properties.</param>
        /// 
        public AbsoluteConvergence(int iterations = 100, double tolerance = 0, double startValue = 0)
        {
            this.MaxIterations = iterations;
            this.tolerance = tolerance;
            this.newValue = startValue;
        }

        /// <summary>
        ///   Gets the watched value before the iteration.
        /// </summary>
        /// 
        public double OldValue { get; private set; }

        /// <summary>
        ///   Gets or sets the watched value after the iteration.
        /// </summary>
        /// 
        public double NewValue
        {
            get { return this.newValue; }
            set
            {
                this.OldValue = this.newValue;
                this.newValue = value;
                this.CurrentIteration++;
            }
        }

        /// <summary>
        ///   Gets or sets the current iteration number.
        /// </summary>
        /// 
        public int CurrentIteration { get; set; }

        /// <summary>
        ///   Gets whether the algorithm has converged.
        /// </summary>
        /// 
        public bool HasConverged
        {
            get
            {
                if (this.maxIterations > 0 && this.CurrentIteration >= this.maxIterations)
                    return true;

                if (this.tolerance > 0)
                {
                    // Stopping criteria is likelihood convergence
                    double delta = Math.Abs(this.OldValue - this.NewValue);

                    if (delta <= this.tolerance)
                        return true;
                }

                // Check if we have reached an invalid or perfectly separable answer
                if (Double.IsNaN(this.NewValue) || (Double.IsInfinity(this.OldValue) && Double.IsInfinity(this.NewValue)))
                    return true;

                return false;
            }
        }


        /// <summary>
        ///   Clears this instance.
        /// </summary>
        /// 
        public void Clear()
        {
            this.CurrentIteration = 0;
            this.NewValue = 0;
            this.OldValue = 0;
        }
    }
}
