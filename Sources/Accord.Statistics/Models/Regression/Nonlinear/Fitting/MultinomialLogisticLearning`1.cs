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

namespace FileFormat.Accord.Statistics.Models.Regression.Nonlinear.Fitting
{
    using System;
    using System.Threading;
    using Accord.MachineLearning.Learning;
    using FileFormat.Accord.Core;
    using FileFormat.Accord.Core.MachineLearning;
    using FileFormat.Accord.Core.Ranges;
    using global::Accord.Math;
    using Math;
    using Math.Core;
    using Math.Matrix;
    using Math.Optimization.Base;

    /// <summary>
    ///   Gradient optimization for Multinomial logistic regression fitting.
    /// </summary>
    /// 
    /// <examples>
    ///   <para>
    ///     The gradient optimization class allows multinomial logistic regression models to be learnt
    ///     using any mathematical optimization algorithm that implements the <see cref="IFunctionOptimizationMethod{TInput, TOutput}"/>
    ///     interface. </para>
    ///   <code source = "Unit Tests\Accord.Tests.Statistics\Models\Regression\MultinomialLogisticGradientDescentTest.cs" region="doc_learn_0" />
    ///   
    /// <para>Using Conjugate Gradient:</para>
    ///   <code source = "Unit Tests\Accord.Tests.Statistics\Models\Regression\MultinomialLogisticGradientDescentTest.cs" region="doc_learn_cg" />
    ///   
    /// <para>Using Gradient Descent:</para>
    ///   <code source = "Unit Tests\Accord.Tests.Statistics\Models\Regression\MultinomialLogisticGradientDescentTest.cs" region="doc_learn_gd" />
    ///   
    /// <para>Using BFGS:</para>
    ///   <code source = "Unit Tests\Accord.Tests.Statistics\Models\Regression\MultinomialLogisticGradientDescentTest.cs" region="doc_learn_bfgs" />
    /// </examples>
    /// 
    public class MultinomialLogisticLearning<TMethod> :
        ISupervisedLearning<MultinomialLogisticRegression, double[], int>,
        ISupervisedLearning<MultinomialLogisticRegression, double[], int[]>,
        ISupervisedLearning<MultinomialLogisticRegression, double[], double[]>,
        ISupervisedLearning<MultinomialLogisticRegression, double[], bool[]>
        where TMethod : IFunctionOptimizationMethod<double[], double>, new()
    {
        [NonSerialized]
        CancellationToken token = new CancellationToken();

        private MultinomialLogisticRegression regression;
        private TMethod method = new TMethod();

        private double[][] inputs;
        private int[] outputs;

        private double[] gradient;
        private double[] log_y_hat;

        int miniBatchSize;
        IntRange[] miniBatches;
        int current = 0;

        /// <summary>
        /// Gets or sets a cancellation token that can be used to
        /// stop the learning algorithm while it is running.
        /// </summary>
        public CancellationToken Token
        {
            get { return this.token; }
            set { this.token = value; }
        }

        /// <summary>
        ///   Gets or sets the optimization method used to optimize
        ///   the parameters (learn) the <see cref="MultinomialLogisticRegression"/>.
        /// </summary>
        /// 
        public TMethod Method
        {
            get { return this.method; }
            set { this.method = value; }
        }

        /// <summary>
        ///   Gets or sets the number of samples to be used as the mini-batch.
        ///   If set to 0 (or a negative number) the total number of training
        ///   samples will be used as the mini-batch.
        /// </summary>
        /// 
        /// <value>The size of the mini batch.</value>
        /// 
        public int MiniBatchSize
        {
            get { return this.miniBatchSize; }
            set { this.miniBatchSize = value; }
        }

        /// <summary>
        ///   Creates a new <see cref="MultinomialLogisticLearning{TMethod}"/>.
        /// </summary>
        /// 
        public MultinomialLogisticLearning()
        {
        }

        /// <summary>
        ///   Creates a new <see cref="MultinomialLogisticLearning{TMethod}"/>.
        /// </summary>
        /// <param name="regression">The regression to estimate.</param>
        /// 
        public MultinomialLogisticLearning(MultinomialLogisticRegression regression)
            : this()
        {
            this.init(regression);
        }

        private void init(MultinomialLogisticRegression regression)
        {
            this.regression = regression;
        }

        private void compute(double[] w, double[] x, double[] log_responses)
        {
            log_responses[0] = 0;
            for (int j = 1, c = 0; j < log_responses.Length; j++)
            {
                double logit = w[c++]; // intercept
                for (int k = 0; k < x.Length; k++)
                    logit += w[c++] * x[k];

                log_responses[j] = logit;
            }

            double sum = Special.LogSumExp(log_responses);

            // Normalize the probabilities
            for (int j = 0; j < log_responses.Length; j++)
                log_responses[j] -= sum;

#if DEBUG
            double[] exp = log_responses.Exp();
            double one = exp.Sum();
            Debug.Assert(one.IsEqual(1, atol: 1e-5));
#endif
        }

        internal double crossEntropy(double[] w)
        {
            double sum = 0;

            // Negative error log-likelihood / cross-entropy error function
            for (int j = 0; j < this.inputs.Length; j++)
            {
                this.compute(w, this.inputs[j], this.log_y_hat);
                sum -= this.log_y_hat[this.outputs[j]];
            }

            return sum / (double)this.inputs.Length;
        }

        internal double[] crossEntropyGradient(double[] w)
        {
            this.gradient.Clear();

            IntRange miniBatch = this.miniBatches[this.current++];
            if (this.current >= this.miniBatches.Length)
                this.current = 0;

            for (int i = miniBatch.Min; i < miniBatch.Max; i++)
            {
                double[] x = this.inputs[i];
                int y = this.outputs[i];

                this.compute(w, x, this.log_y_hat);

                for (int s = 1, c = 0; s < this.log_y_hat.Length; s++)
                {
                    double h = Math.Exp(this.log_y_hat[s]);

                    if (s == y)
                    {
                        this.gradient[c++] += 1 * h - 1;
                        for (int p = 0; p < x.Length; p++)
                            this.gradient[c++] += x[p] * h - x[p];
                    }
                    else
                    {
                        this.gradient[c++] += h;
                        for (int p = 0; p < x.Length; p++)
                            this.gradient[c++] += x[p] * h;
                    }
                }
            }

            for (int i = 0; i < this.gradient.Length; i++)
                this.gradient[i] /= (double)miniBatch.Length;

            return this.gradient;
        }





        /// <summary>
        /// Learns a model that can map the given inputs to the given outputs.
        /// </summary>
        /// <param name="x">The model inputs.</param>
        /// <param name="y">The desired outputs associated with each <paramref name="x">inputs</paramref>.</param>
        /// <param name="weights">The weight of importance for each input-output pair (if supported by the learning algorithm).</param>
        /// <returns>
        /// A model that has learned how to produce <paramref name="y" /> given <paramref name="x" />.
        /// </returns>
        public MultinomialLogisticRegression Learn(double[][] x, int[][] y, double[] weights = null)
        {
            return this.Learn(x, y.ArgMax(dimension: 0), weights);
        }

        /// <summary>
        /// Learns a model that can map the given inputs to the given outputs.
        /// </summary>
        /// <param name="x">The model inputs.</param>
        /// <param name="y">The desired outputs associated with each <paramref name="x">inputs</paramref>.</param>
        /// <param name="weights">The weight of importance for each input-output pair (if supported by the learning algorithm).</param>
        /// <returns>
        /// A model that has learned how to produce <paramref name="y" /> given <paramref name="x" />.
        /// </returns>
        public MultinomialLogisticRegression Learn(double[][] x, double[][] y, double[] weights = null)
        {
            return this.Learn(x, y.ArgMax(dimension: 0), weights);
        }

        /// <summary>
        /// Learns a model that can map the given inputs to the given outputs.
        /// </summary>
        /// <param name="x">The model inputs.</param>
        /// <param name="y">The desired outputs associated with each <paramref name="x">inputs</paramref>.</param>
        /// <param name="weights">The weight of importance for each input-output pair (if supported by the learning algorithm).</param>
        /// <returns>
        /// A model that has learned how to produce <paramref name="y" /> given <paramref name="x" />.
        /// </returns>
        public MultinomialLogisticRegression Learn(double[][] x, bool[][] y, double[] weights = null)
        {
            return this.Learn(x, y.ArgMax(dimension: 0), weights);
        }

        /// <summary>
        /// Learns a model that can map the given inputs to the given outputs.
        /// </summary>
        /// <param name="x">The model inputs.</param>
        /// <param name="y">The desired outputs associated with each <paramref name="x">inputs</paramref>.</param>
        /// <param name="weights">The weight of importance for each input-output pair (if supported by the learning algorithm).</param>
        /// <returns>
        /// A model that has learned how to produce <paramref name="y" /> given <paramref name="x" />.
        /// </returns>
        public MultinomialLogisticRegression Learn(double[][] x, int[] y, double[] weights = null)
        {
            if (weights != null)
                throw new ArgumentException(global::Accord.Properties.Resources.NotSupportedWeights, "weights");

            this.inputs = x;
            this.outputs = y;

            if (this.regression == null)
                this.regression = new MultinomialLogisticRegression(x.Columns(), y.Max() + 1);

            if (this.method.NumberOfVariables != this.regression.NumberOfParameters)
                this.method.NumberOfVariables = this.regression.NumberOfParameters;

            this.method.Function = this.crossEntropy;

            var gom = this.method as IGradientOptimizationMethod;
            if (gom != null)
                gom.Gradient = this.crossEntropyGradient;

            var sc = this.method as ISupportsCancellation;
            if (sc != null)
                sc.Token = this.Token;

            if (this.miniBatchSize <= 0)
                this.miniBatchSize = x.Length;

            this.gradient = new double[this.regression.NumberOfParameters];
            this.log_y_hat = new double[this.regression.NumberOfOutputs];

            this.current = 0;
            this.miniBatches = new IntRange[(int)Math.Floor(x.Length / (double)this.miniBatchSize)];
            for (int i = 0; i < this.miniBatches.Length; i++)
                this.miniBatches[i] = new IntRange(i, Math.Min(i + this.miniBatchSize, x.Length));

            //bool success = method.Minimize();

            for (int i = 0, k = 0; i < this.regression.Coefficients.Length; i++)
                for (int j = 0; j < this.regression.Coefficients[i].Length; j++, k++)
                    this.regression.Coefficients[i][j] = this.method.Solution[k];

            return this.regression;
        }
    }
}
