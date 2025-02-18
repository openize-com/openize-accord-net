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

// #define SERIAL

namespace Openize.Accord.Statistics.Models.Fields.Gradient
{
    using System;
    using System.Threading;
    using System.Threading.Tasks;
    using Features.Base;
    using Openize.Accord.Math;
    using Openize.Accord.Math.Matrix;
    using Math;
    using Openize.Accord.Core;
    using Openize.Accord.Math.Accord.MachineLearning.Learning;
    using Potential_Functions.Base;
    using Potential_Functions.Clique_Factor_Potentials;

    /// <summary>
    ///   Linear Gradient calculator class for <see cref="HiddenConditionalRandomField{T}">
    ///   Hidden Conditional Random Fields</see>.
    /// </summary>
    /// 
    /// <typeparam name="T">The type of the observations being modeled.</typeparam>
    /// 
    public class ForwardBackwardGradient<T> : ParallelLearningBase,
        IHiddenRandomFieldGradient, IDisposable
    {
        private HiddenConditionalRandomField<T> model;
        private IPotentialFunction<T> function;

        private double sigma = 0;

        private ThreadLocal<T[][]> inputs;
        private ThreadLocal<int[]> outputs;
        private ThreadLocal<double[][]> logLikelihoods;

        private ThreadLocal<double[]> gradient;
        private ThreadLocal<double[]> lnZx, lnZxy;
        private ThreadLocal<double> error;

        /// <summary>
        ///   Gets or sets the inputs to be used in the next
        ///   call to the Objective or Gradient functions.
        /// </summary>
        /// 
        public T[][] Inputs
        {
            get { return this.inputs.Value; }
            set
            {
                if (this.inputs.Value != value)
                {
                    this.inputs.Value = value;

                    this.gradient.Value = new double[this.model.Function.Weights.Length];
                    this.lnZx.Value = new double[value.Length];
                    this.lnZxy.Value = new double[value.Length];
                }
            }
        }

        /// <summary>
        ///   Gets or sets the outputs to be used in the next
        ///   call to the Objective or Gradient functions.
        /// </summary>
        /// 
        public int[] Outputs
        {
            get { return this.outputs.Value; }
            set { this.outputs.Value = value; }
        }

        /// <summary>
        ///   Gets or sets the current parameter 
        ///   vector for the model being learned.
        /// </summary>
        /// 
        public double[] Parameters
        {
            get { return this.model.Function.Weights; }
            set { this.model.Function.Weights = value; }
        }

        /// <summary>
        ///   Gets the error computed in the last call
        ///   to the gradient or objective functions.
        /// </summary>
        /// 
        public double LastError
        {
            get { return this.error.Value; }
        }

        /// <summary>
        ///   Gets or sets the amount of the parameter weights
        ///   which should be included in the objective function.
        ///   Default is 0 (do not include regularization).
        /// </summary>
        /// 
        public double Regularization
        {
            get { return this.sigma; }
            set { this.sigma = value; }
        }

        /// <summary>
        ///   Gets the model being trained.
        /// </summary>
        /// 
        public HiddenConditionalRandomField<T> Model
        {
            get { return this.model; }
            set
            {
                this.model = value;
                this.function = this.model.Function;
            }
        }

        /// <summary>
        ///   Initializes a new instance of the <see cref="ForwardBackwardGradient{T}"/> class.
        /// </summary>
        /// 
        public ForwardBackwardGradient()
        {
            this.gradient = new ThreadLocal<double[]>();
            this.lnZx = new ThreadLocal<double[]>();
            this.lnZxy = new ThreadLocal<double[]>();
            this.inputs = new ThreadLocal<T[][]>();
            this.outputs = new ThreadLocal<int[]>();
            this.logLikelihoods = new ThreadLocal<double[][]>();
            this.error = new ThreadLocal<double>();
        }

        /// <summary>
        ///   Initializes a new instance of the <see cref="ForwardBackwardGradient{T}"/> class.
        /// </summary>
        /// 
        /// <param name="model">The model to be trained.</param>
        /// 
        public ForwardBackwardGradient(HiddenConditionalRandomField<T> model)
            : this()
        {
            this.model = model;
            this.function = model.Function;
        }


        /// <summary>
        ///   Computes the gradient (vector of derivatives) vector for
        ///   the cost function, which may be used to guide optimization.
        /// </summary>
        /// 
        /// <param name="parameters">The parameter vector lambda to use in the model.</param>
        /// <param name="input">The inputs to compute the cost function.</param>
        /// <param name="output">The respective outputs to compute the cost function.</param>
        /// <returns>The value of the gradient vector for the given parameters.</returns>
        /// 
        public double[] Gradient(double[] parameters, T[] input, int output)
        {
            this.Inputs = new[] { input };
            this.Outputs = new[] { output };
            this.Parameters = parameters;
            return this.Gradient();
        }

        /// <summary>
        ///   Computes the gradient (vector of derivatives) vector for
        ///   the cost function, which may be used to guide optimization.
        /// </summary>
        /// 
        /// <param name="parameters">The parameter vector lambda to use in the model.</param>
        /// <param name="inputs">The inputs to compute the cost function.</param>
        /// <param name="outputs">The respective outputs to compute the cost function.</param>
        /// <returns>The value of the gradient vector for the given parameters.</returns>
        /// 
        public double[] Gradient(double[] parameters, T[][] inputs, int[] outputs)
        {
            this.Inputs = inputs;
            this.Outputs = outputs;
            this.Parameters = parameters;
            return this.Gradient();
        }

        /// <summary>
        ///   Computes the gradient using the input/outputs stored in this object.
        ///   This method is <b>not</b> thread safe.
        /// </summary>
        /// 
        /// <param name="parameters">The parameter vector lambda to use in the model.</param>
        /// <returns>The value of the gradient vector for the given parameters.</returns>
        /// 
        public double[] Gradient(double[] parameters)
        {
            this.Parameters = parameters;
            return this.Gradient();
        }

        /// <summary>
        ///   Computes the gradient using the input/outputs stored in this object.
        ///   This method is thread-safe.
        /// </summary>
        /// 
        /// <returns>The value of the gradient vector for the given parameters.</returns>
        /// 
        public double[] Gradient()
        {
            // Localize thread locals
            double[][] logLikelihoods = this.logLikelihoods.Value;
            T[][] inputs = this.inputs.Value;
            int[] outputs = this.outputs.Value;
            double[] lnZx = this.lnZx.Value;
            double[] lnZxy = this.lnZxy.Value;
            double[] gradient = this.gradient.Value;

            double error = 0;

            // The previous call to Objective could have computed
            // the log-likelihoods for all input values. However, if
            // this hasn't been the case, compute them now:

            if (logLikelihoods == null)
                this.model.LogLikelihood(inputs, outputs, out logLikelihoods);

            // Compute the partition function using the previously
            // computed likelihoods. Also compute the total error

            // For each x, compute lnZ(x) and lnZ(x,y)
            for (int i = 0; i < inputs.Length; i++)
            {
                double[] lli = logLikelihoods[i];

                // Compute the marginal likelihood
                double sum = Double.NegativeInfinity;
                for (int j = 0; j < lli.Length; j++)
                    sum = Special.LogSum(sum, lli[j]);

                lnZx[i] = sum;
                lnZxy[i] = lli[outputs[i]];

                // compute and return the negative
                // log-likelihood as error function
                error -= lnZxy[i] - lnZx[i];

                Debug.Assert(!Double.IsNaN(error));
            }

            // Now start computing the gradient w.r.t to the
            // feature functions. Each feature function belongs
            // to a factor potential function, so:

            // For each clique potential (factor potential function)
            if (this.ParallelOptions.MaxDegreeOfParallelism == 1)
            {
                Parallel.For(0, this.function.Factors.Length, this.ParallelOptions, c =>
                {
                    this.InnerGradient(this.function.Factors[c], inputs, outputs, lnZx, lnZxy, gradient);
                });
            }
            else
            {
                for (int c = 0; c < this.function.Factors.Length; c++)
                {
                    this.InnerGradient(this.function.Factors[c], inputs, outputs, lnZx, lnZxy, gradient);
                }
            }

            // Reset log-likelihoods so they are recomputed in the next run,
            // either by the Objective function or by the Gradient calculation.

            this.logLikelihoods.Value = null;
            this.error.Value = error;

            Debug.Assert(!Double.IsNaN(error));

            return gradient; // return the gradient.
        }

        private void InnerGradient(FactorPotential<T> factor, T[][] inputs, int[] outputs, double[] lnZx, double[] lnZxy, double[] gradient)
        {
            int factorIndex = factor.Index;

            // Compute all forward and backward matrices to be
            //  used in the feature functions marginal computations.

            double[][,] lnFwds = new double[inputs.Length][,];
            double[][,] lnBwds = new double[inputs.Length][,];
            for (int i = 0; i < inputs.Length; i++)
            {
                lnFwds[i] = ForwardBackwardAlgorithm.LogForward(factor, inputs[i], factorIndex);
                lnBwds[i] = ForwardBackwardAlgorithm.LogBackward(factor, inputs[i], factorIndex);
            }

            double[] marginals = new double[this.function.Outputs];

            // For each feature in the factor potential function
            int end = factor.FactorParameters.Offset + factor.FactorParameters.Count;
            for (int k = factor.FactorParameters.Offset; k < end; k++)
            {
                IFeature<T> feature = this.function.Features[k];
                double parameter = this.function.Weights[k];

                if (Double.IsInfinity(parameter))
                {
                    gradient[k] = 0; continue;
                }


                // Compute the two marginal sums for the gradient calculation
                // as given in eq. 1.52 of Sutton, McCallum; "An introduction to
                // Conditional Random Fields for Relational Learning". The sums
                // will be computed in the log domain for numerical stability.

                double lnsum1 = Double.NegativeInfinity;
                double lnsum2 = Double.NegativeInfinity;

                // For each training sample (sequences)
                for (int i = 0; i < inputs.Length; i++)
                {
                    T[] x = inputs[i]; // training input
                    int y = outputs[i];  // training output

                    // Compute marginals for all possible outputs
                    for (int j = 0; j < marginals.Length; j++)
                        marginals[j] = Double.NegativeInfinity;

                    // However, making the assumption that each factor is responsible for only 
                    // one output label, we can compute the marginal only for the current factor
                    marginals[factorIndex] = feature.LogMarginal(lnFwds[i], lnBwds[i], x, factorIndex);

                    // The first term contains a marginal probability p(w|x,y), which is
                    // exactly a marginal distribution of the clamped CRF (eq. 1.46).
                    lnsum1 = Special.LogSum(lnsum1, (marginals[y] == lnZxy[i]) ? 0 : marginals[y] - lnZxy[i]);

                    // The second term contains a different marginal p(w,y|x) which is the
                    // same marginal probability required in as fully-observed CRF.
                    for (int j = 0; j < marginals.Length; j++)
                        lnsum2 = Special.LogSum(lnsum2, marginals[j] - lnZx[i]);

                    Debug.Assert(!marginals.HasNaN());
                    Debug.Assert(!Double.IsNaN(lnsum1));
                    Debug.Assert(!Double.IsNaN(lnsum2));
                }

                // Compute the current derivative
                double sum1 = Math.Exp(lnsum1);
                double sum2 = Math.Exp(lnsum2);
                double derivative = sum1 - sum2;

                if (sum1 == sum2) derivative = 0;

                Debug.Assert(!Double.IsNaN(derivative));

                // Include regularization derivative if required
                if (this.sigma != 0) derivative -= parameter / this.sigma;

                gradient[k] = -derivative;
            }
        }

        /// <summary>
        ///   Computes the objective (cost) function for the Hidden
        ///   Conditional Random Field (negative log-likelihood).
        /// </summary>
        /// 
        /// <param name="parameters">The parameter vector lambda to use in the model.</param>
        /// <param name="inputs">The inputs to compute the cost function.</param>
        /// <param name="outputs">The respective outputs to compute the cost function.</param>
        /// <returns>The value of the objective function for the given parameters.</returns>
        /// 
        public double Objective(double[] parameters, T[][] inputs, int[] outputs)
        {
            this.Inputs = inputs;
            this.Outputs = outputs;
            return this.Objective(parameters);
        }

        /// <summary>
        ///   Computes the objective (cost) function for the Hidden
        ///   Conditional Random Field (negative log-likelihood) using
        ///   the input/outputs stored in this object.
        /// </summary>
        /// 
        /// <param name="parameters">The parameter vector lambda to use in the model.</param>
        /// 
        public double Objective(double[] parameters)
        {
            this.Parameters = parameters;
            return this.Objective();
        }

        /// <summary>
        ///   Computes the objective (cost) function for the Hidden
        ///   Conditional Random Field (negative log-likelihood) using
        ///   the input/outputs stored in this object.
        /// </summary>
        /// 
        public double Objective()
        {
            var inputs = this.inputs.Value;
            var outputs = this.outputs.Value;
            var parameters = this.Parameters;
            double sumSquaredWeights = 0;

            // Regularization
            if (this.sigma != 0)
            {
                for (int i = 0; i < parameters.Length; i++)
                    if (!(Double.IsInfinity(parameters[i]) || Double.IsNaN(parameters[i])))
                        sumSquaredWeights += parameters[i] * parameters[i];
                sumSquaredWeights = sumSquaredWeights / (2.0 * this.sigma);
            }

            double[][] logLikelihoods;
            double logLikelihood = this.model.LogLikelihood(inputs, outputs, out logLikelihoods);

            this.logLikelihoods.Value = logLikelihoods;

            // Maximize the log-likelihood and minimize
            // a portion of the sum of squared weights
            double objective = logLikelihood - sumSquaredWeights;

            return -objective; // convert to a minimization problem
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
        ~ForwardBackwardGradient()
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
                if (this.inputs != null)
                    this.inputs.Dispose();
                if (this.outputs != null)
                    this.outputs.Dispose();
                if (this.logLikelihoods != null)
                    this.logLikelihoods.Dispose();
                if (this.gradient != null)
                    this.gradient.Dispose();
                if (this.lnZx != null)
                    this.lnZx.Dispose();
                if (this.lnZxy != null)
                    this.lnZxy.Dispose();
                if (this.error != null)
                    this.error.Dispose();

                this.inputs = null;
                this.outputs = null;
                this.gradient = null;
                this.logLikelihoods = null;
                this.lnZxy = null;
                this.lnZx = null;
                this.error = null;
            }
        }

        #endregion

    }
}
