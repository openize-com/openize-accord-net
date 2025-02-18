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

namespace Openize.Accord.Statistics.Models.Fields.Learning.Hidden
{
    using System;
    using System.ComponentModel;
    using System.Threading;
    using System.Threading.Tasks;
    using Accord.MachineLearning.Learning;
    using Openize.Accord.Math.Matrix;
    using Gradient;
    using Math.Convergence;
    using Math.Convergence.Base;
    using Openize.Accord.Core;
    using Openize.Accord.Core.MachineLearning;
    using Openize.Accord.Math.Accord.Statistics;

    /// <summary>
    ///   Stochastic Gradient Descent learning algorithm for <see cref="HiddenConditionalRandomField{T}">
    ///   Hidden Conditional Hidden Fields</see>.
    /// </summary>
    /// 
    /// <example>
    ///   <code source="Unit Tests\Accord.Tests.Statistics\Models\Fields\HiddenConditionalRandomFieldTest.cs" region="doc_learn_1" />
    ///   <code source="Unit Tests\Accord.Tests.Statistics\Models\Fields\HiddenConditionalRandomFieldTest.cs" region="doc_learn_2" />
    ///   <code source="Unit Tests\Accord.Tests.Statistics\Models\Fields\HiddenConditionalRandomFieldTest.cs" region="doc_learn_3" />
    /// </example>
    /// 
    /// <seealso cref="HiddenQuasiNewtonLearning{T}" />
    /// <seealso cref="HiddenResilientGradientLearning{T}"/>
    /// 
    public class HiddenGradientDescentLearning<T> : BaseHiddenConditionalRandomFieldLearning<T>,
        ISupervisedLearning<HiddenConditionalRandomField<T>, T[], int>, IParallel,
        IHiddenConditionalRandomFieldLearning<T>, IConvergenceLearning, IDisposable
    {

        private double learningRate = 100;
        private ISingleValueConvergence convergence;

        //private double decay = 0.9;
        //private double tau = 0.5;
        private double stepSize;

        private bool stochastic = true;
        private double[] gradient;

        private ForwardBackwardGradient<T> calculator = new ForwardBackwardGradient<T>();


        private Object lockObj = new Object();



        /// <summary>
        ///   Gets or sets the learning rate to use as the gradient
        ///   descent step size. Default value is 1e-1.
        /// </summary>
        /// 
        public double LearningRate
        {
            get { return this.learningRate; }
            set { this.learningRate = value; }
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
        public int MaxIterations
        {
            get { return this.convergence.MaxIterations; }
            set { this.convergence.MaxIterations = value; }
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
        ///   Gets or sets a value indicating whether this <see cref="HiddenGradientDescentLearning&lt;T&gt;"/>
        ///   should use stochastic gradient updates.
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
        /// Gets or sets the parallelization options for this algorithm.
        /// </summary>
        /// <value>The parallel options.</value>
        public ParallelOptions ParallelOptions
        {
            get { return ((IParallel)this.calculator).ParallelOptions; }
            set { ((IParallel)this.calculator).ParallelOptions = value; }
        }



        /// <summary>
        ///   Occurs when the current learning progress has changed.
        /// </summary>
        /// 
        public event EventHandler<ProgressChangedEventArgs> ProgressChanged;

        /// <summary>
        ///   Initializes a new instance of the <see cref="HiddenGradientDescentLearning&lt;T&gt;"/> class.
        /// </summary>
        /// 

        public HiddenGradientDescentLearning()
        {
            this.convergence = new RelativeConvergence();
        }

        /// <summary>
        ///   Initializes a new instance of the <see cref="HiddenGradientDescentLearning&lt;T&gt;"/> class.
        /// </summary>
        /// 
        /// <param name="model">The model to be trained.</param>
        /// 
        public HiddenGradientDescentLearning(HiddenConditionalRandomField<T> model)
            : this()
        {
            this.Model = model;
            this.init();
        }

        private void init()
        {
            this.calculator.Model = this.Model;
            this.gradient = new double[this.Model.Function.Weights.Length];
        }

        /// <summary>
        ///   Resets the step size.
        /// </summary>
        /// 
        public void Reset()
        {
            this.convergence.Clear();
            this.stepSize = 0;
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

                // For each training point
                if (this.ParallelOptions.MaxDegreeOfParallelism == 1)
                {
                    for (int i = 0; i < observations.Length; i++)
                        this.iterate(observations, outputs, i, ref error, ref progress);
                }
                else
                {
                    Parallel.For(0, observations.Length, this.ParallelOptions, i =>
                        this.iterate(observations, outputs, i, ref error, ref progress));
                }

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
            this.stepSize = this.learningRate / (this.convergence.CurrentIteration + 1);

            // Update the model using a dynamic step size
            for (int i = 0; i < parameters.Length; i++)
            {
                if (Double.IsInfinity(parameters[i])) continue;

                parameters[i] -= this.stepSize * this.gradient[i];

                Debug.Assert(!Double.IsNaN(parameters[i]));
                Debug.Assert(!Double.IsPositiveInfinity(parameters[i]));
            }


            return this.convergence.NewValue = error;
        }

        private void iterate(T[][] observations, int[] outputs, int i, ref double error, ref int progress)
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
            this.calculator.Inputs = new[] { observations };
            this.calculator.Outputs = new[] { output };

            double[] gradient = this.calculator.Gradient();
            double[] parameters = this.Model.Function.Weights;
            double stepSize = this.learningRate / this.convergence.CurrentIteration;

            // Update the model using a dynamic step size
            for (int i = 0; i < parameters.Length; i++)
            {
                if (Double.IsInfinity(parameters[i])) continue;

                parameters[i] -= stepSize * gradient[i];
            }

            return this.calculator.LastError;
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
        ///   the <see cref="ForwardBackwardGradient{T}"/> is reclaimed by garbage
        ///   collection.
        /// </summary>
        /// 
        ~HiddenGradientDescentLearning()
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
