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

#pragma warning disable 612, 618

namespace FileFormat.Accord.Statistics.Models.Fields.Learning.Hidden
{
    using System;
    using System.ComponentModel;
    using System.Threading;
    using System.Threading.Tasks;
    using Accord.MachineLearning.Learning;
    using FileFormat.Accord.Core;
    using FileFormat.Accord.Core.MachineLearning;
    using global::Accord.Math;
    using Gradient;
    using Math.Accord.Statistics;
    using Math.Convergence;
    using Math.Convergence.Base;
    using Math.Matrix;

    /// <summary>
    ///   Resilient Gradient Learning.
    /// </summary>
    /// 
    /// <typeparam name="T">The type of the observations being modeled.</typeparam>
    ///
    /// <example>
    ///   <code source="Unit Tests\Accord.Tests.Statistics\Models\Fields\HiddenConditionalRandomFieldTest.cs" region="doc_learn_1" />
    ///   <code source="Unit Tests\Accord.Tests.Statistics\Models\Fields\HiddenConditionalRandomFieldTest.cs" region="doc_learn_2" />
    ///   <code source="Unit Tests\Accord.Tests.Statistics\Models\Fields\HiddenConditionalRandomFieldTest.cs" region="doc_learn_3" />
    ///   
    ///   <para>
    ///   The next example shows how to use the learning algorithms in a real-world dataset,
    ///   including training and testing in separate sets and evaluating its performance:</para>
    ///   <code source="Unit Tests\Accord.Tests.Statistics\Models\Fields\Learning\ResilientGradientHiddenLearningTest.cs" region="doc_learn_pendigits" />
    /// </example>
    /// 
    /// <seealso cref="HiddenQuasiNewtonLearning{T}"/>
    /// <seealso cref="HiddenGradientDescentLearning{T}"/>
    /// 
    public class HiddenResilientGradientLearning<T> : BaseHiddenConditionalRandomFieldLearning<T>,
        ISupervisedLearning<HiddenConditionalRandomField<T>, T[], int>, IParallel,
        IHiddenConditionalRandomFieldLearning<T>, IConvergenceLearning, IDisposable
    {

        private ForwardBackwardGradient<T> calculator = new ForwardBackwardGradient<T>();
        private ISingleValueConvergence convergence;

        private double initialStep = 0.0125;
        private double deltaMax = 50.0;
        private double deltaMin = 1e-6;

        private double etaMinus = 0.5;
        private double etaPlus = 1.2;
        private bool stochastic = true;

        private double[] gradient;
        private double[] previousGradient;

        private Object lockObj = new Object();

        // update values, also known as deltas
        private double[] weightsUpdates;



        /// <summary>
        ///   Gets or sets a value indicating whether this <see cref="HiddenGradientDescentLearning&lt;T&gt;"/>
        ///   should use stochastic gradient updates. Default is true.
        /// </summary>
        /// 
        /// <value><c>true</c> for stochastic updates; otherwise, <c>false</c>.</value>
        /// 
        public bool Stochastic
        {
            get { return this.stochastic; }
            set { this.stochastic = value; }
        }

        /// <summary>
        ///   Gets or sets the amount of the parameter weights
        ///   which should be included in the objective function.
        ///   Default is 0 (do not include regularization).
        /// </summary>
        /// 
        public double Regularization
        {
            get { return this.calculator.Regularization; }
            set { this.calculator.Regularization = value; }
        }

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
        ///   performed by the learning algorithm.
        /// </summary>
        /// 
        public int MaxIterations
        {
            get { return this.convergence.MaxIterations; }
            set { this.convergence.MaxIterations = value; }
        }

        /// <summary>
        ///   Gets or sets the number of performed iterations.
        /// </summary>
        /// 
        public int CurrentIteration
        {
            get { return this.convergence.CurrentIteration; }
        }

        /// <summary>
        /// Gets or sets whether the algorithm has converged.
        /// </summary>
        /// <value><c>true</c> if this instance has converged; otherwise, <c>false</c>.</value>
        public bool HasConverged
        {
            get { return this.convergence.HasConverged; }
        }

        /// <summary>
        /// Gets or sets the parallelization options for this algorithm.
        /// </summary>
        /// <value>The parallel options.</value>
        public ParallelOptions ParallelOptions
        {
            get { return ((IParallel)this.calculator).ParallelOptions; }
            set { ((IParallel)this.calculator).ParallelOptions = value; }
        }

        /// <summary>
        ///   Initializes a new instance of the <see cref="HiddenResilientGradientLearning{T}"/> class.
        /// </summary>
        /// 
        public HiddenResilientGradientLearning()
        {
            this.convergence = new RelativeConvergence(iterations: 100, tolerance: 0, checks: 3);
        }

        /// <summary>
        ///   Initializes a new instance of the <see cref="HiddenResilientGradientLearning{T}"/> class.
        /// </summary>
        /// 
        /// <param name="model">Model to teach.</param>
        /// 
        public HiddenResilientGradientLearning(HiddenConditionalRandomField<T> model)
            : this()
        {
            this.Model = model;
            this.init();
        }

        private void init()
        {
            this.calculator.Model = this.Model;
            int parameters = this.Model.Function.Weights.Length;
            this.gradient = new double[parameters];
            this.previousGradient = new double[parameters];
            this.weightsUpdates = new double[parameters];
            this.Reset(this.initialStep); // Initialize steps
        }



        /// <summary>
        ///   Runs one iteration of the learning algorithm with the
        ///   specified input training observation and corresponding
        ///   output label.
        /// </summary>
        /// 
        /// <param name="observations">The training observations.</param>
        /// <param name="outputs">The observation's labels.</param>
        /// 
        /// <returns>The error in the last iteration.</returns>
        /// 
        [Obsolete("Please use Learn(x, y) instead.")]
        public double Run(T[][] observations, int[] outputs)
        {
            return this.InnerRun(observations, outputs);
        }

        /// <summary>
        ///   Runs the learning algorithm.
        /// </summary>
        /// 
        protected override double InnerRun(T[][] observations, int[] outputs)
        {
            this.init();
            this.convergence.Clear();

            do
            {
                this.RunEpoch(observations, outputs);
                if (this.Token.IsCancellationRequested)
                    break;
            }
            while (!this.convergence.HasConverged);

            return this.convergence.NewValue;
        }

        /// <summary>
        ///   Runs the learning algorithm with the specified input
        ///   training observations and corresponding output labels.
        /// </summary>
        /// 
        /// <param name="observations">The training observations.</param>
        /// <param name="outputs">The observation's labels.</param>
        /// 
        /// <returns>The error in the last iteration.</returns>
        /// 
        public double RunEpoch(T[][] observations, int[] outputs)
        {
            double error = 0;

            if (this.stochastic)
            {

                // In batch mode, we will use the average of the gradients
                // at each point as a better estimate of the true gradient.
                Array.Clear(this.gradient, 0, this.gradient.Length);

                int progress = 0;


#if SERIAL      // For each training point
                for (int i = 0; i < observations.Length; i++)
#else
                Parallel.For(0, observations.Length, this.ParallelOptions, i =>
#endif
                {
                    this.calculator.Inputs = new[] { observations[i] };
                    this.calculator.Outputs = new[] { outputs[i] };

                    // Compute the estimated gradient
                    double[] estimate = this.calculator.Gradient();

                    lock (this.lockObj)
                    {
                        // Accumulate
                        for (int j = 0; j < estimate.Length; j++)
                            this.gradient[j] += estimate[j];
                        error += this.calculator.LastError;
                    }

                    int current = Interlocked.Increment(ref progress);
                    double percent = current / (double)observations.Length * 100.0;
                    this.OnProgressChanged(new ProgressChangedEventArgs((int)percent, i));

                    Debug.Assert(!this.gradient.HasNaN());
                }
#if !SERIAL
);
#endif
                // Compute the average gradient
                for (int i = 0; i < this.gradient.Length; i++)
                    this.gradient[i] /= observations.Length;
            }
            else
            {
                this.calculator.Inputs = observations;
                this.calculator.Outputs = outputs;

                // Compute the true gradient
                this.gradient = this.calculator.Gradient();

                error = this.calculator.LastError;
            }

            double[] parameters = this.Model.Function.Weights;

            // Do the Resilient Backpropagation parameter update
            for (int k = 0; k < this.calculator.Parameters.Length; k++)
            {
                if (Double.IsInfinity(parameters[k])) continue;

                double S = this.previousGradient[k] * this.gradient[k];

                if (S > 0.0)
                {
                    this.weightsUpdates[k] = Math.Min(this.weightsUpdates[k] * this.etaPlus, this.deltaMax);
                    parameters[k] -= Math.Sign(this.gradient[k]) * this.weightsUpdates[k];
                    this.previousGradient[k] = this.gradient[k];
                }
                else if (S < 0.0)
                {
                    this.weightsUpdates[k] = Math.Max(this.weightsUpdates[k] * this.etaMinus, this.deltaMin);
                    this.previousGradient[k] = 0.0;
                }
                else
                {
                    parameters[k] -= Math.Sign(this.gradient[k]) * this.weightsUpdates[k];
                    this.previousGradient[k] = this.gradient[k];
                }
            }

            Debug.Assert(!this.Model.Function.Weights.HasNaN());

            return this.convergence.NewValue = error;
        }

        /// <summary>
        ///   Runs one iteration of the learning algorithm with the
        ///   specified input training observation and corresponding
        ///   output label.
        /// </summary>
        /// 
        /// <param name="observations">The training observations.</param>
        /// <param name="output">The observation's labels.</param>
        /// 
        /// <returns>The error in the last iteration.</returns>
        /// 
        public double Run(T[] observations, int output)
        {
            return this.RunEpoch(new[] { observations }, new[] { output });
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


        #region IDisposable Members

        /// <summary>
        ///   Performs application-defined tasks associated with freeing,
        ///   releasing, or resetting unmanaged resources.
        /// </summary>
        /// 
        public void Dispose()
        {
            this.Dispose(true);
            GC.SuppressFinalize(this);
        }

        /// <summary>
        ///   Releases unmanaged resources and performs other cleanup operations before
        ///   the <see cref="HiddenResilientGradientLearning{T}"/> is reclaimed by garbage
        ///   collection.
        /// </summary>
        /// 
        ~HiddenResilientGradientLearning()
        {
            this.Dispose(false);
        }

        /// <summary>
        ///   Releases unmanaged and - optionally - managed resources
        /// </summary>
        /// 
        /// <param name="disposing"><c>true</c> to release both managed 
        /// and unmanaged resources; <c>false</c> to release only unmanaged
        /// resources.</param>
        /// 
        protected virtual void Dispose(bool disposing)
        {
            if (disposing)
            {
                // free managed resources
                if (this.calculator != null)
                {
                    this.calculator.Dispose();
                    this.calculator = null;
                }
            }
        }

        #endregion

    }
}