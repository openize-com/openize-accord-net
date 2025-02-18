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

namespace Openize.Accord.Math.Convergence
{
    using System;
    using Base;

    /// <summary>
    ///   Relative convergence criteria.
    /// </summary>
    /// 
    /// <remarks>
    ///   This class can be used to track progress and convergence
    ///   of methods which rely on the relative change of a value.
    /// </remarks>
    /// 
    /// <example>
    /// <code>
    ///   // Create a new convergence criteria with unlimited iterations
    ///   var criteria = new RelativeConvergence(iterations: 0, tolerance: 0.1);
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
    ///   // maximum of 11 iterations with a final value
    ///   // of 1234.56:
    ///   
    ///   int iterations = criteria.CurrentIteration; // 11
    ///   double value = criteria.OldValue; // 1234.56
    /// </code>
    /// </example>
    /// 
    public class RelativeConvergence : ISingleValueConvergence
    {

        private double tolerance = 0;
        private int maxIterations = 100;
        private double newValue;
        private double startValue = 0;

        private int checks;
        private int maxChecks = 1;


        /// <summary>
        ///   Gets or sets the maximum relative change in the watched value
        ///   after an iteration of the algorithm used to detect convergence.
        ///   Default is zero.
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
        ///   Initializes a new instance of the <see cref="RelativeConvergence"/> class.
        /// </summary>
        /// 
        public RelativeConvergence()
        {
            this.init();
        }

        /// <summary>
        ///   Initializes a new instance of the <see cref="RelativeConvergence"/> class.
        /// </summary>
        /// 
        /// <param name="iterations">The maximum number of iterations which should be
        ///   performed by the iterative algorithm. Setting to zero indicates there
        ///   is no maximum number of iterations. Default is 100.</param>
        /// <param name="tolerance">The maximum relative change in the watched value
        ///   after an iteration of the algorithm used to detect convergence.
        ///   Default is 0.</param>
        /// 
        public RelativeConvergence(int iterations, double tolerance)
        {
            this.init(iterations, tolerance);
        }

        /// <summary>
        ///   Initializes a new instance of the <see cref="RelativeConvergence"/> class.
        /// </summary>
        /// 
        /// <param name="iterations">The maximum number of iterations which should be
        ///   performed by the iterative algorithm. Setting to zero indicates there
        ///   is no maximum number of iterations. Default is 0.</param>
        /// <param name="tolerance">The maximum relative change in the watched value
        ///   after an iteration of the algorithm used to detect convergence.
        ///   Default is 0.</param>
        /// <param name="checks">The minimum number of convergence checks that the
        ///   iterative algorithm should pass before convergence can be declared
        ///   reached.</param>
        /// 
        public RelativeConvergence(int iterations, double tolerance, int checks)
        {
            this.init(iterations, tolerance, checks);
        }

        /// <summary>
        ///   Initializes a new instance of the <see cref="RelativeConvergence"/> class.
        /// </summary>
        /// 
        /// <param name="iterations">The maximum number of iterations which should be
        ///   performed by the iterative algorithm. Setting to zero indicates there
        ///   is no maximum number of iterations. Default is 0.</param>
        /// <param name="tolerance">The maximum relative change in the watched value
        ///   after an iteration of the algorithm used to detect convergence.
        ///   Default is 0.</param>
        /// <param name="checks">The minimum number of convergence checks that the
        ///   iterative algorithm should pass before convergence can be declared
        ///   reached.</param>
        /// <param name="startValue">The initial value for the <see cref="NewValue"/> and
        ///   <see cref="OldValue"/> properties.</param>
        /// 
        public RelativeConvergence(int iterations = 100, double tolerance = 0, int checks = 1, double startValue = 0)
        {
            this.init(iterations, tolerance, checks, startValue);
        }

        private void init(int iterations = 100, double tolerance = 0, int checks = 1, double startValue = 0)
        {
            this.MaxIterations = iterations;
            this.tolerance = tolerance;
            this.maxChecks = checks;
            this.startValue = startValue;

            this.Clear();
        }

        /// <summary>
        ///   Gets or sets the watched value before the iteration.
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
        ///   Gets the current iteration number.
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
                bool converged = this.checkConvergence();

                this.checks = converged ? this.checks + 1 : 0;

                return this.checks >= this.maxChecks;
            }
        }

        private bool checkConvergence()
        {
            if (this.maxIterations > 0 && this.CurrentIteration >= this.maxIterations)
                return true;

            if (this.tolerance > 0)
            {
                // Stopping criteria is likelihood convergence
                if (this.Delta <= this.tolerance * Math.Abs(this.OldValue))
                    return true;
            }

            // Check if we have reached an invalid or perfectly separable answer
            if (Double.IsNaN(this.NewValue) || Double.IsInfinity(this.NewValue))
                return true;

            return false;
        }

        /// <summary>
        ///   Gets the absolute difference between the <see cref="NewValue"/> and <see cref="OldValue"/>
        ///   as as <c>Math.Abs(OldValue - NewValue)</c>.
        /// </summary>
        /// 
        public double Delta
        {
            get { return Math.Abs(this.OldValue - this.NewValue); }
        }

        /// <summary>
        ///   Gets the relative difference between the <see cref="NewValue"/> and <see cref="OldValue"/>
        ///   as <c>Math.Abs(OldValue - NewValue) / Math.Abs(OldValue)</c>.
        /// </summary>
        /// 
        public double RelativeDelta
        {
            get { return this.Delta / Math.Abs(this.OldValue); }
        }

        /// <summary>
        ///   Gets the initial value for the <see cref="NewValue"/> 
        ///   and <see cref="OldValue"/> properties.
        /// </summary>
        /// 
        public double StartValue
        {
            get { return this.startValue; }
        }

        /// <summary>
        ///   Resets this instance, reverting all iteration statistics
        ///   statistics (number of iterations, last error) back to zero.
        /// </summary>
        /// 
        public void Clear()
        {
            this.NewValue = this.startValue;
            this.OldValue = this.startValue;
            this.CurrentIteration = 0;
            this.checks = 0;
        }


        /// <summary>
        /// Returns a <see cref="System.String" /> that represents this instance.
        /// </summary>
        /// <returns>A <see cref="System.String" /> that represents this instance.</returns>
        public override string ToString()
        {
            if (this.HasConverged)
                return String.Format("Converged: {0} < {1}", this.RelativeDelta, this.tolerance);
            return String.Format("Running: {0} > {1}", this.RelativeDelta, this.Tolerance);
        }

    }
}
