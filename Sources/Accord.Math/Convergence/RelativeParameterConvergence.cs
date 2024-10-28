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
    ///   Relative parameter change convergence criteria.
    /// </summary>
    /// 
    /// <remarks>
    ///   This class can be used to track progress and convergence
    ///   of methods which rely on the maximum relative change of
    ///   the values within a parameter vector.
    /// </remarks>
    /// 
    /// <example>
    /// <code>
    ///   // Converge if the maximum change amongst all parameters is less than 0.1:
    ///   var criteria = new RelativeParameterConvergence(iterations: 0, tolerance: 0.1);
    /// 
    ///   int progress = 1;
    ///   double[] parameters = { 12345.6, 952.12, 1925.1 };
    ///   
    ///   do
    ///   {
    ///       // Do some processing...
    /// 
    ///       // Update current iteration information:
    ///       criteria.NewValues = parameters.Divide(progress++);
    /// 
    ///   } while (!criteria.HasConverged);
    /// 
    /// 
    ///   // The method will converge after reaching the 
    ///   // maximum of 11 iterations with a final value
    ///   // of { 1234.56, 95.212, 192.51 }:
    /// 
    ///   int iterations = criteria.CurrentIteration; // 11
    ///   var v = criteria.OldValues; // { 1234.56, 95.212, 192.51 }
    /// 
    /// </code>
    /// </example>
    /// 
    public class RelativeParameterConvergence : IConvergence<double[]>
    {
        private double[] oldValues;
        private double[] newValues;

        private double tolerance = 0;
        private int maxIterations = 100;
        private double maxChange;


        /// <summary>
        ///   Gets or sets the maximum change in the watched value
        ///   after an iteration of the algorithm used to detect convergence.
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
        ///   Gets or sets the maximum number of iterations
        ///   performed by the iterative algorithm.
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
        ///   Initializes a new instance of the <see cref="RelativeParameterConvergence"/> class.
        /// </summary>
        /// 
        public RelativeParameterConvergence()
        {
        }

        /// <summary>
        ///   Initializes a new instance of the <see cref="RelativeParameterConvergence"/> class.
        /// </summary>
        /// 
        /// <param name="iterations">The maximum number of iterations which should be
        ///   performed by the iterative algorithm. Setting to zero indicates there
        ///   is no maximum number of iterations. Default is 0.</param>
        /// <param name="tolerance">The maximum relative change in the watched value
        ///   after an iteration of the algorithm used to detect convergence.
        ///   Default is 0.</param>
        /// 
        public RelativeParameterConvergence(int iterations, double tolerance)
        {
            this.MaxIterations = iterations;
            this.tolerance = tolerance;
        }

        /// <summary>
        ///   Gets the maximum relative parameter
        ///   change after the last iteration.
        /// </summary>
        /// 
        public double Delta
        {
            get { return this.maxChange; }
        }

        /// <summary>
        ///   Gets or sets the watched value before the iteration.
        /// </summary>
        /// 
        public double[] OldValues { get { return this.oldValues; } }


        /// <summary>
        ///   Gets or sets the watched value after the iteration.
        /// </summary>
        /// 
        public double[] NewValues
        {
            get { return this.newValues; }
            set
            {
                this.oldValues = this.newValues;
                this.newValues = (double[])value.Clone();
                this.CurrentIteration++;
            }
        }

        /// <summary>
        ///   Gets or sets the current iteration number.
        /// </summary>
        /// 
        public int CurrentIteration { get; set; }

        /// <summary>
        ///   Gets whether the algorithm has diverged.
        /// </summary>
        /// 
        public bool HasDiverged
        {
            get
            {
                for (int i = 0; i < this.NewValues.Length; i++)
                    if (Double.IsNaN(this.NewValues[i]) || Double.IsInfinity(this.NewValues[i]))
                        return true;
                return false;
            }
        }

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

                if (this.NewValues == null && this.OldValues == null)
                    return true;
                if (this.OldValues == null)
                    return false;
                if (this.NewValues.Length == 0 || this.OldValues.Length == 0)
                    return true;

                // Check if we have reached an invalid or perfectly separable answer
                for (int i = 0; i < this.NewValues.Length; i++)
                    if (Double.IsNaN(this.NewValues[i]) || Double.IsInfinity(this.NewValues[i]))
                        return true;

                // Update and verify stop criteria
                if (this.tolerance > 0)
                {
                    // Stopping criteria is likelihood convergence
                    this.maxChange = Math.Abs(this.OldValues[0] - this.NewValues[0]) / Math.Abs(this.OldValues[0]);

                    if (Double.IsNaN(this.maxChange))
                        this.maxChange = 0;

                    for (int i = 1; i < this.OldValues.Length; i++)
                    {
                        double delta = Math.Abs(this.OldValues[i] - this.NewValues[i]) / Math.Abs(this.OldValues[i]);

                        if (delta > this.maxChange)
                            this.maxChange = delta;
                    }

                    if (Double.IsNaN(this.maxChange))
                        return true;

                    if (this.maxChange <= this.tolerance)
                        return true;
                }

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
            this.newValues = null;
            this.oldValues = null;
        }

        double[] IConvergence<double[]>.NewValue
        {
            // TODO: Remove this explicit implementation.
            get { return this.NewValues; }
            set { this.NewValues = value; }
        }

    }
}
