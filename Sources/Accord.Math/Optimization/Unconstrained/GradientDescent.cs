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
    using Base;
    using FileFormat.Accord.Math.Convergence;

    /// <summary>
    ///   Gradient Descent (GD) for unconstrained optimization.
    /// </summary>
    /// 
    /// <seealso cref="ConjugateGradient"/>
    /// <seealso cref="BoundedBroydenFletcherGoldfarbShanno"/>
    /// <seealso cref="BroydenFletcherGoldfarbShanno"/>
    /// <seealso cref="TrustRegionNewtonMethod"/>
    /// 
    public class GradientDescent : BaseGradientOptimizationMethod, IGradientOptimizationMethod
    {

        private RelativeConvergence convergence = new RelativeConvergence();

        private double eta = 1e-3;
        private int numberOfUpdatesBeforeConvergenceCheck = 1;

        /// <summary>
        ///   Occurs when the current learning progress has changed.
        /// </summary>
        /// 
        public event EventHandler<ProgressChangedEventArgs> ProgressChanged;

        /// <summary>
        ///   Gets or sets the learning rate. Default is 1e-3.
        /// </summary>
        /// 
        public double LearningRate
        {
            get { return this.eta; }
            set
            {
                if (value <= 0)
                    throw new ArgumentOutOfRangeException("value", "Learning rate should be higher than 0.");
                this.eta = value;
            }
        }

        /// <summary>
        ///   Gets or sets the maximum change in the average log-likelihood
        ///   after an iteration of the algorithm used to detect convergence.
        ///   Default is 1e-5.
        /// </summary>
        /// 
        public double Tolerance
        {
            get { return this.convergence.Tolerance; }
            set { this.convergence.Tolerance = value; }
        }

        /// <summary>
        ///   Gets or sets the maximum number of iterations
        ///   performed by the learning algorithm. Default is 0.
        /// </summary>
        /// 
        public int Iterations
        {
            get { return this.convergence.MaxIterations; }
            set { this.convergence.MaxIterations = value; }
        }

        /// <summary>
        ///   Creates a new instance of the GD optimization algorithm.
        /// </summary>
        /// 
        public GradientDescent()
        {
            this.Iterations = 0;
            this.Tolerance = 1e-5;
        }


        /// <summary>
        ///   Implements the actual optimization algorithm. This
        ///   method should try to minimize the objective function.
        /// </summary>
        /// 
        protected override bool Optimize()
        {
            this.convergence.Clear();

            int updates = 0;

            do
            {
                if (this.Token.IsCancellationRequested)
                    break;

                double[] gradient = this.Gradient(this.Solution);
                for (int i = 0; i < this.Solution.Length; i++)
                    this.Solution[i] -= this.eta * gradient[i];

                updates++;

                if (updates >= this.numberOfUpdatesBeforeConvergenceCheck)
                {
                    this.convergence.NewValue = this.Function(this.Solution);
                    updates = 0;
                }
            }
            while (!this.convergence.HasConverged);

            return true;
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

    }
}
