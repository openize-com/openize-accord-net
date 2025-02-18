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

#define SERIAL

namespace Openize.Accord.Statistics.Models.Markov.Learning.Base
{
    using System;
    using System.Linq;
    using Accord.MachineLearning.Learning;
    using Distributions;
    using Openize.Accord.Math;
    using Openize.Accord.Math.Core;
    using Openize.Accord.Math.Matrix;
    using Math;
    using Math.Convergence;
    using Openize.Accord.Core;
    using Openize.Accord.Core.Exceptions;
    using Openize.Accord.Math.Accord.Statistics;
    using Openize.Accord.Math.Convergence.Base;
    using Openize.Accord.Statistics.Distributions.Fitting.Base;

    /// <summary>
    ///   Base class for implementations of the Baum-Welch learning algorithm.
    ///   This class cannot be instantiated.
    /// </summary>
    /// 
    public abstract class BaseBaumWelchLearningOptions<TModel, TDistribution, TObservation, TOptions> :
        BaseBaumWelchLearning<TModel, TDistribution, TObservation, TOptions>,
        IUnsupervisedLearning<TModel, TObservation[], int[]>,
        IConvergenceLearning
        where TModel : HiddenMarkovModel<TDistribution, TObservation>
        where TDistribution : IFittableDistribution<TObservation, TOptions>
        where TOptions : class, IFittingOptions
    {

        /// <summary>
        ///   Creates a new instance of the Baum-Welch learning algorithm.
        /// </summary>
        /// 
        public BaseBaumWelchLearningOptions()
        {
        }

        /// <summary>
        ///   Creates a new instance of the Baum-Welch learning algorithm.
        /// </summary>
        /// 
        public BaseBaumWelchLearningOptions(TModel model)
            : base(model)
        {
        }

        /// <summary>
        /// Fits one emission distribution. This method can be override in a
        /// base class in order to implement special fitting options.
        /// </summary>
        protected override void Fit(int index, TObservation[] values, double[] weights)
        {
            this.Model.Emissions[index].Fit(values, weights, this.FittingOptions);
        }
    }

    /// <summary>
    ///   Base class for implementations of the Baum-Welch learning algorithm.
    ///   This class cannot be instantiated.
    /// </summary>
    /// 
    public abstract class BaseBaumWelchLearning<TModel, TDistribution, TObservation, TOptions> :
        BaseHiddenMarkovModelLearning<TModel, TObservation>,
        IUnsupervisedLearning<TModel, TObservation[], int[]>,
        IConvergenceLearning
        where TOptions : class, IFittingOptions
        where TModel : HiddenMarkovModel<TDistribution, TObservation>
        where TDistribution : IFittableDistribution<TObservation>
    {

        private RelativeConvergence convergence = new RelativeConvergence(
            iterations: 0, tolerance: 1e-5,
            startValue: double.NegativeInfinity);

        private TObservation[][] vectorObservations;
        private TObservation[] samples;

#if SERIAL
        double[,] lnFwd;
        double[,] lnBwd;
        double[] sampleWeights;
#endif

        /// <summary>
        ///   Gets all observations as a single vector.
        /// </summary>
        /// 
        protected TObservation[][] Observations { get { return this.vectorObservations; } }

        /// <summary>
        /// Gets or sets convergence parameters.
        /// </summary>
        /// <value>The convergence parameters.</value>
        public IConvergence Convergence { get { return this.convergence; } }

        /// <summary>
        ///   Gets or sets the distribution fitting options
        ///   to use when estimating distribution densities
        ///   during learning.
        /// </summary>
        /// <value>The distribution fitting options.</value>
        /// 
        public TOptions FittingOptions { get; set; }

        /// <summary>
        ///   Gets the log-likelihood of the model at the last iteration.
        /// </summary>
        /// 
        public double LogLikelihood { get; set; }

        /// <summary>
        ///   Gets or sets the function that initializes the emission
        ///   distributions in the hidden Markov Models.
        /// </summary>
        /// 
        public Func<int, TDistribution> Emissions { get; set; }

        /// <summary>
        ///   Creates a new instance of the Baum-Welch learning algorithm.
        /// </summary>
        /// 
        public BaseBaumWelchLearning(TModel model)
            : base(model)
        {
        }

        /// <summary>
        ///   Creates a new instance of the Baum-Welch learning algorithm.
        /// </summary>
        /// 
        public BaseBaumWelchLearning()
        {
            this.Emissions = (stateIndex) =>
            {
                try
                {
                    return Activator.CreateInstance<TDistribution>();
                }
                catch
                {
                    throw new InvalidOperationException("Please set the Emissions property to specify how the initial distributions should be created.");
                }
            };
        }

        /// <summary>
        ///   Gets or sets the maximum change in the average log-likelihood
        ///   after an iteration of the algorithm used to detect convergence.
        /// </summary>
        /// 
        /// <remarks>
        ///   This is the likelihood convergence limit L between two iterations of the algorithm. The
        ///   algorithm will stop when the change in the likelihood for two consecutive iterations
        ///   has not changed by more than L percent of the likelihood. If left as zero, the
        ///   algorithm will ignore this parameter and iterate over a number of fixed iterations
        ///   specified by the previous parameter.
        /// </remarks>
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
        /// <remarks>
        ///   This is the maximum number of iterations to be performed by the learning algorithm. If
        ///   specified as zero, the algorithm will learn until convergence of the model average
        ///   likelihood respecting the desired limit.
        /// </remarks>
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
            set { this.convergence.CurrentIteration = value; }
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
        ///   Gets the Ksi matrix of log probabilities created during
        ///   the last iteration of the Baum-Welch learning algorithm.
        /// </summary>
        /// 
        public double[][][,] LogKsi { get; protected set; }

        /// <summary>
        ///   Gets the Gamma matrix of log probabilities created during
        ///   the last iteration of the Baum-Welch learning algorithm.
        /// </summary>
        /// 
        public double[][,] LogGamma { get; protected set; }

        /// <summary>
        ///   Gets the sample weights in the last iteration of the
        ///   Baum-Welch learning algorithm.
        /// </summary>
        /// 
        public double[] LogWeights { get; protected set; }


        /// <summary>
        ///   Learns a model that can map the given inputs to the desired outputs.
        /// </summary>
        /// 
        /// <param name="x">The model inputs.</param>
        /// <param name="weights">The weight of importance for each input sample.</param>
        /// 
        /// <returns>A model that has learned how to produce suitable outputs
        ///   given the input data <paramref name="x"/>.</returns>
        /// 
        public TModel Learn(TObservation[][] x, double[] weights = null)
        {
            // Initial argument checks
            CheckArgs(x, weights);

            // Baum-Welch algorithm.

            // The Baum–Welch algorithm is a particular case of a generalized expectation-maximization
            // (GEM) algorithm. It can compute maximum likelihood estimates and posterior mode estimates
            // for the parameters (transition and emission probabilities) of an HMM, when given only
            // emissions as training data.

            // The algorithm has two steps:
            //  - Calculating the forward probability and the backward probability for each HMM state;
            //  - On the basis of this, determining the frequency of the transition-emission pair values
            //    and dividing it by the probability of the entire string. This amounts to calculating
            //    the expected count of the particular transition-emission pair. Each time a particular
            //    transition is found, the value of the quotient of the transition divided by the probability
            //    of the entire string goes up, and this value can then be made the new value of the transition.


            this.samples = x.Concatenate();
            this.vectorObservations = x;

            if (this.Model == null)
                this.Model = this.Create(x);

            MarkovHelperMethods.CheckObservationDimensions(x, this.Model);

            if (this.MaxIterations > 0 && this.CurrentIteration >= this.MaxIterations)
                return this.Model;

            // Grab model information
            int states = this.Model.NumberOfStates;
            var logA = this.Model.LogTransitions;
            var logP = this.Model.LogInitial;

            // Initialize the algorithm
            int N = x.Length;
            double logN = Math.Log(N);
            this.LogKsi = new double[N][][,];
            this.LogGamma = new double[N][,];
            this.LogWeights = new double[N];
            if (weights != null)
                weights.Log(result: this.LogWeights);

            for (int i = 0; i < x.Length; i++)
            {
                int T = x[i].Length;

                this.LogKsi[i] = new double[T][,];
                this.LogGamma[i] = new double[T, states];

                for (int t = 0; t < this.LogKsi[i].Length; t++)
                    this.LogKsi[i][t] = new double[states, states];
            }

            int TMax = x.Max(x_i => x_i.Length);
#if SERIAL
            this.lnFwd = new double[TMax, states];
            this.lnBwd = new double[TMax, states];
            this.sampleWeights = new double[this.samples.Length];
#endif

            this.convergence.CurrentIteration--;
            bool hasUpdated = false;

            do
            {
                if (this.Token.IsCancellationRequested)
                    break;

                // Initialize the model log-likelihood
                this.LogLikelihood = this.Expectation(x, TMax);
                this.convergence.NewValue = this.LogLikelihood;

                // Check for convergence
                if (hasUpdated && this.convergence.HasConverged)
                    break;

                if (this.Token.IsCancellationRequested)
                    break;

                // 3. Continue with parameter re-estimation
                // 3.1 Re-estimation of initial state probabilities 
#if SERIAL
                for (int i = 0; i < logP.Length; i++)
#else
                Parallel.For(0, logP.Length, ParallelOptions, i =>
#endif
                {
                    double lnsum = Double.NegativeInfinity;
                    for (int k = 0; k < this.LogGamma.Length; k++)
                        lnsum = Special.LogSum(lnsum, this.LogGamma[k][0, i]);
                    logP[i] = lnsum - logN;
                }
#if !SERIAL
                );
#endif

                // 3.2 Re-estimation of transition probabilities 
#if SERIAL
                for (int i = 0; i < states; i++)
#else
                Parallel.For(0, states, ParallelOptions, i =>
#endif
                {
                    for (int j = 0; j < states; j++)
                    {
                        double lnnum = Double.NegativeInfinity;
                        double lnden = Double.NegativeInfinity;

                        for (int k = 0; k < this.LogGamma.Length; k++)
                        {
                            int T = x[k].Length;

                            for (int t = 0; t < T - 1; t++)
                            {
                                lnnum = Special.LogSum(lnnum, this.LogKsi[k][t][i, j]);
                                lnden = Special.LogSum(lnden, this.LogGamma[k][t, i]);
                            }
                        }

                        logA[i][j] = (lnnum == lnden) ? 0 : lnnum - lnden;
                        Debug.Assert(!Double.IsNaN(logA[i][j]));
                    }
                }
#if !SERIAL
                );
#endif

                // 3.3 Re-estimation of emission probabilities
                this.UpdateEmissions(); // discrete and continuous
                hasUpdated = true;

            } while (true);

            return this.Model;
        }

        private double Expectation(TObservation[][] x, int maxLength)
        {
            int states = this.Model.NumberOfStates;
            double logLikelihood = 0;

            // For each sequence in the observations input
#if SERIAL
            for (int i = 0; i < x.Length; i++)
            {
#else
            object lockObj = new object();

            Parallel.For(0, x.Length, ParallelOptions,
                () => new
                {
                    lnFwd = new double[maxLength, states],
                    lnBwd = new double[maxLength, states],
                    partialSum = new[] { 0.0 }
                },

                (i, loopState, local) =>
                {
                    double[,] lnFwd = local.lnFwd;
                    double[,] lnBwd = local.lnBwd;
                    double[] partialSum = local.partialSum;
#endif
                int T = x[i].Length;
                double[,] logGamma = this.LogGamma[i];
                double w = this.LogWeights[i];


                // 1st step - Calculating the forward probability and the
                //            backward probability for each HMM state.
                this.ComputeForwardBackward(i, this.lnFwd, this.lnBwd);


                // 2nd step - Determining the frequency of the transition-emission pair values
                //            and dividing it by the probability of the entire string.

                // Calculate gamma values for next computations
#if SERIAL
                for (int t = 0; t < T; t++)
#else
                Parallel.For(0, T, ParallelOptions, t =>
#endif
                {
                    double lnsum = Double.NegativeInfinity;
                    for (int k = 0; k < states; k++)
                    {
                        logGamma[t, k] = this.lnFwd[t, k] + this.lnBwd[t, k] + w;
                        lnsum = Special.LogSum(lnsum, logGamma[t, k]);
                    }

                    Debug.Assert(!Double.IsNaN(lnsum));

                    // Normalize if different from zero
                    if (lnsum != Double.NegativeInfinity)
                        for (int k = 0; k < states; k++)
                            logGamma[t, k] = logGamma[t, k] - lnsum;
                }
#if !SERIAL
);
#endif
                // Calculate ksi values for next computations
                this.ComputeKsi(i, this.lnFwd, this.lnBwd);

                double ll = Double.NegativeInfinity;
                for (int j = 0; j < states; j++)
                    ll = Special.LogSum(ll, this.lnFwd[T - 1, j]);

#if !SERIAL
                    partialSum[0] += ll;
                    return local;
                },

                (local) =>
                {
                    lock (lockObj)
                        logLikelihood += local.partialSum[0];
                });
#else
                logLikelihood += ll;
            }
#endif

            return logLikelihood / (double)x.Length;
        }

        private static void CheckArgs(TObservation[][] observations, double[] weights)
        {
            for (int i = 0; i < observations.Length; i++)
                if (observations[i].Length == 0)
                    throw new ArgumentException("The observation at position {0} has zero length.".Format(i), "observations");

            if (weights != null)
                if (weights.Length != observations.Length)
                    throw new DimensionMismatchException("weights");
        }


        /// <summary>
        ///   Computes the ksi matrix of probabilities for a given observation
        ///   referenced by its index in the input training data.
        /// </summary>
        /// <param name="index">The index of the observation in the input training data.</param>
        /// <param name="lnFwd">The matrix of forward probabilities for the observation.</param>
        /// <param name="lnBwd">The matrix of backward probabilities for the observation.</param>
        /// 
        protected void ComputeKsi(int index, double[,] lnFwd, double[,] lnBwd)
        {
            int states = this.Model.NumberOfStates;
            double[][] logA = this.Model.LogTransitions;
            TDistribution[] B = this.Model.Emissions;

            var sequence = this.vectorObservations[index];

            int T = sequence.Length;
            var logKsi = this.LogKsi[index];
            var w = this.LogWeights[index];

            // Note: cannot be parallelized due to 
            // dependency on previous values of t
            for (int t = 0; t < T - 1; t++)
            {
                double lnsum = Double.NegativeInfinity;
                TObservation x = sequence[t + 1];

                for (int i = 0; i < states; i++)
                {
                    for (int j = 0; j < states; j++)
                    {
                        double b = B[j].LogProbabilityFunction(x);
                        logKsi[t][i, j] = lnFwd[t, i] + lnBwd[t + 1, j] + logA[i][j] + b + w;
                        lnsum = Special.LogSum(lnsum, logKsi[t][i, j]);
                    }
                }

                Debug.Assert(!Double.IsNaN(lnsum));

                // Normalize if different from zero
                if (lnsum != Double.NegativeInfinity)
                    for (int i = 0; i < states; i++)
                        for (int j = 0; j < states; j++)
                            logKsi[t][i, j] = logKsi[t][i, j] - lnsum;
            }
        }

        /// <summary>
        ///   Updates the emission probability matrix.
        /// </summary>
        /// <remarks>
        ///   Implementations of this method should use the observations
        ///   in the training data and the Gamma probability matrix to
        ///   update the probability distributions of symbol emissions.
        /// </remarks>
        /// 
        protected void UpdateEmissions()
        {
            var B = this.Model.Emissions;

            // For each state i in the model
#if SERIAL
            for (int i = 0; i < B.Length; i++)
#else
            Parallel.For(0, B.Length, ParallelOptions,
                () => new double[this.samples.Length],
                (i, loopState, sampleWeights) =>
#endif
            {
                double lnsum = Double.NegativeInfinity;

                // For each observation sequence k
                for (int k = 0, w = 0; k < this.vectorObservations.Length; k++)
                {
                    int T = this.vectorObservations[k].Length;

                    // For each observation t in k
                    for (int t = 0; t < T; t++, w++)
                    {
                        this.sampleWeights[w] = this.LogGamma[k][t, i];
                        lnsum = Special.LogSum(lnsum, this.sampleWeights[w]);
                    }
                }

                Debug.Assert(!Double.IsNaN(lnsum));

                if (lnsum != Double.NegativeInfinity)
                {
                    for (int w = 0; w < this.sampleWeights.Length; w++)
                        this.sampleWeights[w] = this.sampleWeights[w] - lnsum;

                    // Convert to probabilities
                    for (int w = 0; w < this.sampleWeights.Length; w++)
                    {
                        double p = Math.Exp(this.sampleWeights[w]);
                        this.sampleWeights[w] = (Double.IsNaN(p) || Double.IsInfinity(p)) ? 0.0 : p;
                    }

                    // Estimate the distribution for state i
                    this.Fit(i, this.samples, this.sampleWeights);
                }
#if SERIAL
            }
#else
                return weights;
        },

                (weights) => { }
            );
#endif
        }

        /// <summary>
        ///   Fits one emission distribution. This method can be override in a
        ///   base class in order to implement special fitting options.
        /// </summary>
        /// 
        protected virtual void Fit(int index, TObservation[] values, double[] weights)
        {
#if DEBUG
            if (this.Model.Emissions[index] is IFittableDistribution<TObservation, TOptions>)
                throw new Exception("A more specialized method should have been called.");
#endif
            IFittableDistribution<TObservation> dist = this.Model.Emissions[index];
            if (this.FittingOptions == null)
                dist.Fit(values, weights);
            else
                dist.Fit(values, weights, this.FittingOptions);
        }


        /// <summary>
        ///   Computes the forward and backward probabilities matrices
        ///   for a given observation referenced by its index in the
        ///   input training data.
        /// </summary>
        /// <param name="index">The index of the observation in the input training data.</param>
        /// <param name="lnFwd">Returns the computed forward probabilities matrix.</param>
        /// <param name="lnBwd">Returns the computed backward probabilities matrix.</param>
        /// 
        protected void ComputeForwardBackward(int index, double[,] lnFwd, double[,] lnBwd)
        {
            int states = this.Model.NumberOfStates;
            int T = this.vectorObservations[index].Length;

            Debug.Assert(lnBwd.GetLength(0) >= T);
            Debug.Assert(lnBwd.GetLength(1) == states);
            Debug.Assert(lnFwd.GetLength(0) >= T);
            Debug.Assert(lnFwd.GetLength(1) == states);

            ForwardBackwardAlgorithm.LogForward(this.Model, this.vectorObservations[index], lnFwd);
            ForwardBackwardAlgorithm.LogBackward(this.Model, this.vectorObservations[index], lnBwd);
        }

    }
}
