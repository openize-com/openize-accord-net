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

namespace Openize.Accord.Statistics.Models.Markov.Learning
{
#pragma warning disable 612, 618
    using System;
    using Base;
    using Distributions;
    using Openize.Accord.Math;
    using Openize.Accord.Math.Matrix;
    using Math;
    using Openize.Accord.Core;
    using Openize.Accord.Math.Accord.Statistics;
    using Openize.Accord.Statistics.Distributions.Fitting.Base;
    using Openize.Accord.Statistics.Distributions.Univariate.Base;

    /// <summary>
    ///   Obsolete. Please use <see cref="BaumWelchLearning{TDistribution, TObservation}"/> instead.
    /// </summary>
    /// 
    [Obsolete("Please use BaumWelchLearning<TDistribution, TObservation> instead.")]
    public class BaumWelchLearning<TDistribution> : BaseBaumWelchLearning, IUnsupervisedLearning, IConvergenceLearning
        where TDistribution : IDistribution
    {

        private HiddenMarkovModel<TDistribution> model;

        private IFittingOptions fittingOptions;

        private double[][][] vectorObservations;
        private Array samples;
        private double[] weights;

        /// <summary>
        ///   Gets the model being trained.
        /// </summary>
        /// 
        public HiddenMarkovModel<TDistribution> Model
        {
            get { return this.model; }
        }

        /// <summary>
        ///   Gets or sets the distribution fitting options
        ///   to use when estimating distribution densities
        ///   during learning.
        /// </summary>
        /// <value>The distribution fitting options.</value>
        /// 
        public IFittingOptions FittingOptions
        {
            get { return this.fittingOptions; }
            set { this.fittingOptions = value; }
        }

        /// <summary>
        ///   Creates a new instance of the Baum-Welch learning algorithm.
        /// </summary>
        /// 
        public BaumWelchLearning(HiddenMarkovModel<TDistribution> model)
            : base(model)
        {
            this.model = model;
        }


        /// <summary>
        ///   Runs the Baum-Welch learning algorithm for hidden Markov models.
        /// </summary>
        /// <remarks>
        ///   Learning problem. Given some training observation sequences O = {o1, o2, ..., oK}
        ///   and general structure of HMM (numbers of hidden and visible states), determine
        ///   HMM parameters M = (A, B, pi) that best fit training data. 
        /// </remarks>
        /// <param name="observations">
        ///   The sequences of univariate or multivariate observations used to train the model.
        ///   Can be either of type double[] (for the univariate case) or double[][] for the
        ///   multivariate case.
        /// </param>
        /// <returns>
        ///   The average log-likelihood for the observations after the model has been trained.
        /// </returns>
        /// 
        public new double Run(params Array[] observations)
        {

            // Convert the generic representation to a vector of multivariate sequences
            this.vectorObservations = new double[observations.Length][][];
            for (int i = 0; i < this.vectorObservations.Length; i++)
                this.vectorObservations[i] = MarkovHelperMethods.checkAndConvert(observations[i], this.model.Dimension);


            // Sample array, used to store all observations as a sample vector
            //   will be useful when fitting the distribution models.
            if (this.model.Emissions[0] is IUnivariateDistribution)
            {
                this.samples = Matrix.Concatenate(
                    Matrix.Stack(this.vectorObservations));
            }
            else
            {
                this.samples = Matrix.Stack(this.vectorObservations);
            }

            this.weights = new double[this.samples.Length];


            return base.Run(observations);
        }

        /// <summary>
        ///   Computes the ksi matrix of probabilities for a given observation
        ///   referenced by its index in the input training data.
        /// </summary>
        /// <param name="index">The index of the observation in the input training data.</param>
        /// <param name="lnFwd">The matrix of forward probabilities for the observation.</param>
        /// <param name="lnBwd">The matrix of backward probabilities for the observation.</param>
        /// 
        protected override void ComputeKsi(int index, double[,] lnFwd, double[,] lnBwd)
        {
            int states = this.model.States;
            double[,] logA = this.model.Transitions;
            TDistribution[] B = this.model.Emissions;

            var sequence = this.vectorObservations[index];

            int T = sequence.Length;
            var logKsi = this.LogKsi[index];
            var w = this.LogWeights[index];


            for (int t = 0; t < T - 1; t++)
            {
                double lnsum = Double.NegativeInfinity;
                double[] x = sequence[t + 1];

                for (int i = 0; i < states; i++)
                {
                    for (int j = 0; j < states; j++)
                    {
                        double b = B[j].LogProbabilityFunction(x);
                        logKsi[t][i, j] = lnFwd[t, i] + lnBwd[t + 1, j] + logA[i, j] + b + w;
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
        protected override void UpdateEmissions()
        {
            var B = this.model.Emissions;

            // For each state i in the model
            for (int i = 0; i < B.Length; i++)
            {
                double lnsum = Double.NegativeInfinity;

                // For each observation sequence k
                for (int k = 0, w = 0; k < this.vectorObservations.Length; k++)
                {
                    int T = this.vectorObservations[k].Length;

                    // For each observation t in k
                    for (int t = 0; t < T; t++, w++)
                    {
                        this.weights[w] = this.LogGamma[k][t, i];
                        lnsum = Special.LogSum(lnsum, this.weights[w]);
                    }
                }

                Debug.Assert(!Double.IsNaN(lnsum));

                if (lnsum != Double.NegativeInfinity)
                {
                    for (int w = 0; w < this.weights.Length; w++)
                        this.weights[w] = this.weights[w] - lnsum;

                    // Convert to probabilities
                    for (int w = 0; w < this.weights.Length; w++)
                    {
                        double p = Math.Exp(this.weights[w]);
                        this.weights[w] = (Double.IsNaN(p) || Double.IsInfinity(p)) ? 0.0 : p;
                    }

                    // Estimate the distribution for state i
                    B[i].Fit(this.samples, this.weights, this.fittingOptions);
                }
            }
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
        protected override void ComputeForwardBackward(int index, double[,] lnFwd, double[,] lnBwd)
        {
            int states = this.model.States;
            int T = this.vectorObservations[index].Length;

            Debug.Assert(lnBwd.GetLength(0) >= T);
            Debug.Assert(lnBwd.GetLength(1) == states);
            Debug.Assert(lnFwd.GetLength(0) >= T);
            Debug.Assert(lnFwd.GetLength(1) == states);

            ForwardBackwardAlgorithm.LogForward(this.model, this.vectorObservations[index], lnFwd);
            ForwardBackwardAlgorithm.LogBackward(this.model, this.vectorObservations[index], lnBwd);
        }


    }
}
