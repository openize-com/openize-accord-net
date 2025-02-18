// Accord Statistics Library
// The Accord.NET Framework
// http://accord-framework.net
//
// Copyright © Clement Schiano di Colella, 2015-2016
// clement.schiano@gmail.com
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

namespace Openize.Accord.Statistics.Models.Regression.Nonlinear.Fitting
{
    using System;
    using System.Collections.Generic;
    using System.Threading;
    using Accord.MachineLearning.Learning;
    using Base;
    using Openize.Accord.Math.Core;
    using Openize.Accord.Math.Matrix;
    using Linear;
    using Linear.Fitting;
    using Openize.Accord.Math.Optimization.Losses;

    /// <summary>
    ///   Non-negative Least Squares for <see cref="MultipleLinearRegression"/> optimization.
    /// </summary>
    /// 
    /// <remarks>
    /// <para>
    ///   References:
    ///   <list type="bullet">
    ///     <item><description>
    ///      Donghui Chen and Robert J.Plemmons, Nonnegativity Constraints in Numerical Analysis.
    ///      Available on: http://users.wfu.edu/plemmons/papers/nonneg.pdf </description></item>
    ///   </list></para>
    /// </remarks>
    /// 
    /// <example>
    ///  <para>
    ///   The following example shows how to fit a multiple linear regression model with the 
    ///   additional constraint that none of its coefficients should be negative. For this
    ///   we can use the <see cref="NonNegativeLeastSquares"/> learning algorithm instead of
    ///   the <see cref="OrdinaryLeastSquares"/> used above.</para>
    ///   
    /// <code source="Unit Tests\Accord.Tests.Statistics\Models\Regression\NonNegativeLeastSquaresTest.cs" region="doc_learn" />
    /// </example>
    /// 
    /// <seealso cref="OrdinaryLeastSquares"/>
    /// <seealso cref="MultipleLinearRegression"/>
    /// 
#pragma warning disable 612, 618
    public class NonNegativeLeastSquares : IRegressionFitting,
        ISupervisedLearning<MultipleLinearRegression, double[], double>
#pragma warning restore 612, 618
    {
        [NonSerialized]
        CancellationToken token = new CancellationToken();

        MultipleLinearRegression regression;

        List<int> p = new List<int>();
        List<int> r = new List<int>();
        double[] s;
        double tolerance = 0.001;
        double[][] scatter;
        double[] vector;
        double[] W;
        double[][] X;
        int cols;
        int maxIter;

        /// <summary>
        ///   Gets the coefficient vector being fitted.
        /// </summary>
        /// 
        public double[] Coefficients { get { return this.regression.Weights; } }

        /// <summary>
        ///   Gets or sets the maximum number of iterations to be performed.
        /// </summary>
        /// 
        public int MaxIterations
        {
            get { return this.maxIter; }
            set { this.maxIter = value; }
        }

        /// <summary>
        ///   Gets or sets the tolerance for detecting
        ///   convergence. Default is 0.001.
        /// </summary>
        /// 
        public double Tolerance
        {
            get { return this.tolerance; }
            set { this.tolerance = value; }
        }

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
        ///   Initializes a new instance of the <see cref="NonNegativeLeastSquares"/> class.
        /// </summary>
        /// 
        public NonNegativeLeastSquares()
        {
        }

        /// <summary>
        ///   Initializes a new instance of the <see cref="NonNegativeLeastSquares"/> class.
        /// </summary>
        /// 
        /// <param name="regression">The regression to be fitted.</param>
        /// 
        public NonNegativeLeastSquares(MultipleLinearRegression regression)
        {
            this.init(regression);
        }

        private void init(MultipleLinearRegression regression)
        {
            this.regression = regression;
            this.cols = regression.Weights.Length;
            this.s = new double[this.cols];
            this.W = new double[this.cols];
        }


        /// <summary>
        ///   Runs the fitting algorithm.
        /// </summary>
        /// 
        /// <param name="inputs">The input training data.</param>
        /// <param name="outputs">The output associated with each of the outputs.</param>
        /// 
        /// <returns>
        ///   The sum of squared errors after the learning.
        /// </returns>
        /// 
        [Obsolete("Please use the Learn() method instead.")]
        public double Run(double[][] inputs, double[] outputs)
        {
            this.regression = this.Learn(inputs, outputs);
            return new SquareLoss(inputs).Loss(this.regression.Transform(inputs));
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
        public MultipleLinearRegression Learn(double[][] x, double[] y, double[] weights = null)
        {
            if (weights != null)
                throw new ArgumentException(global::Accord.Properties.Resources.NotSupportedWeights, "weights");

            if (this.regression == null)
                this.init(new MultipleLinearRegression { NumberOfInputs = x.Columns() });

            this.X = x;
            this.scatter = this.X.TransposeAndDot(this.X);
            this.vector = this.X.TransposeAndDot(y);

            // Initialization
            this.p.Clear();
            this.r.Clear();
            for (var i = 0; i < this.cols; i++)
                this.r.Add(i);

            var w = this.Coefficients;

            this.ComputeWeights(w);
            var iter = 0;
            int maxWeightIndex;
            this.W.Max(out maxWeightIndex);

            while (this.r.Count > 0 && this.W[maxWeightIndex] > this.tolerance && iter < this.maxIter)
            {
                if (this.Token.IsCancellationRequested)
                    break;

                // Include the index j in P and remove it from R
                if (!this.p.Contains(maxWeightIndex))
                    this.p.Add(maxWeightIndex);

                if (this.r.Contains(maxWeightIndex))
                    this.r.Remove(maxWeightIndex);

                this.GetSP();
                int iter2 = 0;

                while (GetElements(this.s, this.p).Min() <= 0 && iter2 < this.maxIter)
                {
                    this.InnerLoop(w);
                    iter2++;
                }

                // 5
                Array.Copy(this.s, w, this.s.Length);

                // 6
                this.ComputeWeights(w);

                this.W.Max(out maxWeightIndex);
                iter++;
            }

            //Coefficients = x;
            return this.regression;
        }

        private void InnerLoop(double[] x)
        {
            var alpha = double.PositiveInfinity;
            foreach (int i in this.p)
            {
                if (this.s[i] <= 0)
                    alpha = System.Math.Min(alpha, x[i] / (x[i] - this.s[i]));
            }

            if (System.Math.Abs(alpha) < 0.001 || double.IsNaN(alpha))
                return;

            x = (this.s.Subtract(x)).Multiply(alpha).Add(x);

            // 4.4 Update R and P
            for (var i = 0; i < this.p.Count;)
            {
                int pItem = this.p[i];
                if (System.Math.Abs(x[pItem]) < double.Epsilon)
                {
                    this.r.Add(pItem);
                    this.p.RemoveAt(i);
                }
                else
                {
                    i++;
                }
            }

            // 4.5 
            this.GetSP();

            // 4.6
            foreach (var i in this.r)
                this.s[i] = 0;
        }

        private void ComputeWeights(double[] x)
        {
            this.W = this.vector.Subtract(this.scatter.Dot(x));
        }

        private void GetSP()
        {
            int[] array = this.p.ToArray();
            double[][] left = this.scatter
                .GetColumns(array)
                .GetRows(array)
                .PseudoInverse();

            double[] columnVector = GetElements(this.vector, this.p);
            double[] result = left.Dot(columnVector);
            for (int i = 0; i < this.p.Count; i++)
                this.s[this.p[i]] = result[i];
        }

        private static double[] GetElements(double[] vector, List<int> elementsIndex)
        {
            var z = new double[elementsIndex.Count];
            for (var i = 0; i < elementsIndex.Count; i++)
                z[i] = vector[elementsIndex[i]];
            return z;
        }

    }
}
