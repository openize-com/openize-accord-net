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

namespace FileFormat.Accord.Math.Optimization.Unconstrained
{
    using System;
    using System.ComponentModel;
    using System.Threading.Tasks;
    using Base;
    using FileFormat.Accord.Core;
    using FileFormat.Accord.Math.Convergence;
    using global::Accord.Math;
    using Matrix;

    /// <summary>
    ///   Resilient Backpropagation method for unconstrained optimization.
    /// </summary>
    /// 
    /// <seealso cref="ConjugateGradient"/>
    /// <seealso cref="BoundedBroydenFletcherGoldfarbShanno"/>
    /// <seealso cref="BroydenFletcherGoldfarbShanno"/>
    /// <seealso cref="TrustRegionNewtonMethod"/>
    /// 
    public class ResilientBackpropagation : BaseGradientOptimizationMethod, IGradientOptimizationMethod
    {

        private RelativeConvergence convergence;

        private double initialStep = 0.0125;
        private double deltaMax = 50.0;
        private double deltaMin = 1e-6;

        private double etaMinus = 0.5;
        private double etaPlus = 1.2;

        private double[] gradient;
        private double[] previousGradient;

        // update values, also known as deltas
        private double[] weightsUpdates;



        /// <summary>
        ///   Occurs when the current learning progress has changed.
        /// </summary>
        /// 
        public event EventHandler<ProgressChangedEventArgs> ProgressChanged;


        /// <summary>
        ///   Gets or sets the maximum possible update step,
        ///   also referred as delta min. Default is 50.
        /// </summary>
        /// 
        public double UpdateUpperBound
        {
            get { return this.deltaMax; }
            set { this.deltaMax = value; }
        }

        /// <summary>
        ///   Gets or sets the minimum possible update step,
        ///   also referred as delta max. Default is 1e-6.
        /// </summary>
        /// 
        public double UpdateLowerBound
        {
            get { return this.deltaMin; }
            set { this.deltaMin = value; }
        }

        /// <summary>
        ///   Gets the decrease parameter, also 
        ///   referred as eta minus. Default is 0.5.
        /// </summary>
        /// 
        public double DecreaseFactor
        {
            get { return this.etaMinus; }
            set
            {
                if (value <= 0 || value >= 1)
                    throw new ArgumentOutOfRangeException("value", "Value should be between 0 and 1.");
                this.etaMinus = value;
            }
        }

        /// <summary>
        ///   Gets the increase parameter, also
        ///   referred as eta plus. Default is 1.2.
        /// </summary>
        /// 
        public double IncreaseFactor
        {
            get { return this.etaPlus; }
            set
            {
                if (value <= 1)
                    throw new ArgumentOutOfRangeException("value", "Value should be higher than 1.");
                this.etaPlus = value;
            }
        }

        /// <summary>
        ///   Gets or sets the maximum change in the average log-likelihood
        ///   after an iteration of the algorithm used to detect convergence.
        /// </summary>
        /// 
        public double Tolerance
        {
            get { return this.convergence.Tolerance; }
            set { this.convergence.Tolerance = value; }
        }

        /// <summary>
        ///   Gets or sets the maximum number of iterations
        ///   performed by the learning algorithm.
        /// </summary>
        /// 
        public int Iterations
        {
            get { return this.convergence.MaxIterations; }
            set { this.convergence.MaxIterations = value; }
        }

        /// <summary>
        ///   Creates a new <see cref="ResilientBackpropagation"/> function optimizer.
        /// </summary>
        /// 
        /// <param name="function">The function to be optimized.</param>
        /// 
        public ResilientBackpropagation(NonlinearObjectiveFunction function)
            : this(function.NumberOfVariables, function.Function, function.Gradient)
        {
        }

        /// <summary>
        ///   Creates a new <see cref="ResilientBackpropagation"/> function optimizer.
        /// </summary>
        /// 
        /// <param name="numberOfVariables">The number of free parameters in the function to be optimized.</param>
        /// <param name="function">The function to be optimized.</param>
        /// <param name="gradient">The gradient of the function.</param>
        /// 
        public ResilientBackpropagation(int numberOfVariables,
            Func<double[], double> function, Func<double[], double[]> gradient)
            : base(numberOfVariables, function, gradient)
        {
        }

        /// <summary>
        ///   Creates a new <see cref="ResilientBackpropagation"/> function optimizer.
        /// </summary>
        /// 
        /// <param name="numberOfVariables">The number of parameters in the function to be optimized.</param>
        /// 
        public ResilientBackpropagation(int numberOfVariables)
            : base(numberOfVariables)
        {
        }

        /// <summary>
        /// Called when the <see cref="IOptimizationMethod{TInput, TOutput}.NumberOfVariables" /> property has changed.
        /// </summary>
        /// 
        /// <param name="numberOfVariables">The number of variables.</param>
        /// 
        protected override void OnNumberOfVariablesChanged(int numberOfVariables)
        {
            base.OnNumberOfVariablesChanged(numberOfVariables);

            this.convergence = new RelativeConvergence();

            this.gradient = new double[numberOfVariables];
            this.previousGradient = new double[numberOfVariables];
            this.weightsUpdates = new double[numberOfVariables];

            // Initialize steps
            this.Reset(this.initialStep);
        }

        /// <summary>
        ///   Implements the actual optimization algorithm. This
        ///   method should try to minimize the objective function.
        /// </summary>
        /// 
        protected override bool Optimize()
        {
            this.convergence.Clear();

            do
            {
                this.runEpoch();
                if (this.Token.IsCancellationRequested)
                    break;
            }
            while (!this.convergence.HasConverged);

            return true;
        }




        private double runEpoch()
        {
            // Compute the true gradient
            this.gradient = this.Gradient(this.Solution);

            double[] parameters = this.Solution;

            // Do the Resilient Backpropagation parameter update
            for (int k = 0; k < parameters.Length; k++)
            {
                if (Double.IsInfinity(parameters[k]) || Double.IsNaN(this.gradient[k]))
                    continue;

                double g = this.gradient[k];
                if (g > 1e100) g = 1e100;
                if (g < -1e100) g = -1e100;

                double S = this.previousGradient[k] * g;

                if (S > 0.0)
                {
                    this.weightsUpdates[k] = Math.Min(this.weightsUpdates[k] * this.etaPlus, this.deltaMax);
                    parameters[k] -= Math.Sign(g) * this.weightsUpdates[k];
                    this.previousGradient[k] = g;
                }
                else if (S < 0.0)
                {
                    this.weightsUpdates[k] = Math.Max(this.weightsUpdates[k] * this.etaMinus, this.deltaMin);
                    this.previousGradient[k] = 0.0;
                }
                else
                {
                    parameters[k] -= Math.Sign(g) * this.weightsUpdates[k];
                    this.previousGradient[k] = g;
                }
            }

            Debug.Assert(!parameters.HasNaN());

            double value = this.Function(parameters);

            return this.convergence.NewValue = value;
        }

        /// <summary>
        ///   Raises the <see cref="E:ProgressChanged"/> event.
        /// </summary>
        /// 
        /// <param name="args">The ProgressChangedEventArgs instance containing the event data.</param>
        /// 
        protected void OnProgressChanged(ProgressChangedEventArgs args)
        {
            if (this.ProgressChanged != null)
                this.ProgressChanged(this, args);
        }

        /// <summary>
        ///   Resets the current update steps using the given learning rate.
        /// </summary>
        /// 
        public void Reset(double rate)
        {
            this.convergence.Clear();

            Parallel.For(0, this.weightsUpdates.Length, i =>
            {
                for (int j = 0; j < this.weightsUpdates.Length; j++)
                    this.weightsUpdates[i] = rate;
            });
        }

    }
}
