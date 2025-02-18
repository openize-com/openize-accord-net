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

namespace Openize.Accord.Statistics.Models.Regression.Nonlinear.Fitting
{
    using System;
    using System.Threading;
    using Accord.MachineLearning.Learning;
    using Base;
    using Openize.Accord.Math.Core;
    using Openize.Accord.Math.Matrix;
    using Math.Decompositions;
    using Math.Decompositions.Base;
    using Openize.Accord.Math.Accord.Statistics;
    using Openize.Accord.Math.Convergence;
    using Openize.Accord.Math.Convergence.Base;

    /// <summary>
    ///   Lower-Bound Newton-Raphson for Multinomial logistic regression fitting.
    /// </summary>
    /// 
    /// <remarks>
    /// <para>
    ///   The Lower Bound principle consists of replacing the second derivative
    ///   matrix by a global lower bound in the Leowner ordering [Böhning, 92].
    ///   In the case of multinomial logistic regression estimation, the Hessian
    ///   of the negative log-likelihood function can be replaced by one of those
    ///   lower bounds, leading to a monotonically converging sequence of iterates.
    ///   Furthermore, [Krishnapuram, Carin, Figueiredo and Hartemink, 2005] also
    ///   have shown that a lower bound can be achieved which does not depend on 
    ///   the coefficients for the current iteration.</para>
    /// 
    /// <para>
    ///   References:
    ///   <list type="bullet">
    ///     <item><description>
    ///      B. Krishnapuram, L. Carin, M.A.T. Figueiredo, A. Hartemink. Sparse Multinomial 
    ///      Logistic Regression: Fast Algorithms and Generalization Bounds. 2005. Available on: 
    ///      http://www.lx.it.pt/~mtf/Krishnapuram_Carin_Figueiredo_Hartemink_2005.pdf </description></item>
    ///     <item><description>
    ///       D. Böhning. Multinomial logistic regression algorithm. Annals of the Institute
    ///       of Statistical Mathematics, 44(9):197 ˝U200, 1992. 2. M. Corney.</description></item>
    ///     <item><description>
    ///       Bishop, Christopher M.; Pattern Recognition and Machine Learning. 
    ///       Springer; 1st ed. 2006.</description></item>
    ///   </list></para>
    /// </remarks>
    ///
    /// <example>
    ///   <code source="Unit Tests\Accord.Tests.Statistics\Models\Regression\MultinomialLogisticRegressionTest.cs" region="doc_learn" />
    /// </example>
    ///
#pragma warning disable 612, 618
    public class LowerBoundNewtonRaphson :
        ISupervisedLearning<MultinomialLogisticRegression, double[], int>,
        ISupervisedLearning<MultinomialLogisticRegression, double[], int[]>,
        ISupervisedLearning<MultinomialLogisticRegression, double[], double[]>,
        IMultipleRegressionFitting, IConvergenceLearning
#pragma warning restore 612, 618
    {
        [NonSerialized]
        CancellationToken token = new CancellationToken();

        private MultinomialLogisticRegression regression;

        private double[] previous;
        private double[] solution;
        private double[] deltas;

        private double[] errors;
        private double[] output;

        private double[,] weights;

        private double[,] lowerBound;
        private double[] gradient;
        private double[,] xxt;

        private int K;
        private int M;
        private int parameterCount;

        private bool computeStandardErrors = true;
        private bool updateLowerBound = true;
        private ISolverMatrixDecomposition<double> decomposition = null;

        private IConvergence<double> convergence;

        /// <summary>
        ///   Gets the previous values for the coefficients which were
        ///   in place before the last learning iteration was performed.
        /// </summary>
        /// 
        public double[] Previous
        {
            get { return this.previous; }
        }

        /// <summary>
        ///   Gets the current values for the coefficients.
        /// </summary>
        /// 
        public double[] Solution
        {
            get { return this.solution; }
        }

        /// <summary>
        ///   Gets or sets a value indicating whether the lower bound 
        ///   should be updated using new data. Default is <c>true</c>.
        /// </summary>
        /// 
        /// <value>
        ///   <c>true</c> if the lower bound should be
        ///   updated; otherwise, <c>false</c>.</value>
        ///   
        public bool UpdateLowerBound
        {
            get { return this.updateLowerBound; }
            set { this.updateLowerBound = value; }
        }

        /// <summary>
        ///   Gets the Lower-Bound matrix being used in place of
        ///   the Hessian matrix in the Newton-Raphson iterations.
        /// </summary>
        /// 
        public double[,] HessianLowerBound
        {
            get { return this.lowerBound; }
        }

        /// <summary>
        ///   Gets the Gradient vector computed in
        ///   the last Newton-Raphson iteration.
        /// </summary>
        /// 
        public double[] Gradient
        {
            get { return this.gradient; }
        }

        /// <summary>
        ///   Gets the total number of parameters in the model.
        /// </summary>
        /// 
        public int Parameters { get { return this.parameterCount; } }


        /// <summary>
        ///   Gets or sets a value indicating whether standard
        ///   errors should be computed in the next iteration.
        /// </summary>
        /// <value>
        /// 	<c>true</c> to compute standard errors; otherwise, <c>false</c>.
        /// </value>
        /// 
        public bool ComputeStandardErrors
        {
            get { return this.computeStandardErrors; }
            set { this.computeStandardErrors = value; }
        }

        /// <summary>
        ///   Please use MaxIterations instead.
        /// </summary>
        [Obsolete("Please use MaxIterations instead.")]
        public int Iterations
        {
            get { return this.convergence.MaxIterations; }
            set { this.convergence.MaxIterations = value; }
        }

        /// <summary>
        /// Gets or sets the maximum number of iterations
        /// performed by the learning algorithm.
        /// </summary>
        public int MaxIterations
        {
            get { return this.convergence.MaxIterations; }
            set { this.convergence.MaxIterations = value; }
        }

        /// <summary>
        /// Gets or sets the tolerance value used to determine
        /// whether the algorithm has converged.
        /// </summary>
        public double Tolerance
        {
            get { return this.convergence.Tolerance; }
            set { this.convergence.Tolerance = value; }
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
        ///   Gets the vector of parameter updates in the last iteration.
        /// </summary>
        /// 
        /// <value>How much each parameter changed after the last update.</value>
        /// 
        public double[] ParameterChange {  get { return this.deltas; } }

        /// <summary>
        ///   Gets the maximum <see cref="ParameterChange">parameter change</see> in the last 
        ///   iteration. If this value is less than <see cref="Tolerance"/>, the algorithm
        ///   has converged.
        /// </summary>
        /// 
        public double MaximumChange {  get { return this.convergence.NewValue; } }

        /// <summary>
        ///   Creates a new <see cref="LowerBoundNewtonRaphson"/>.
        /// </summary>
        /// 
        public LowerBoundNewtonRaphson()
        {
            this.convergence = new AbsoluteConvergence();
        }

        /// <summary>
        ///   Creates a new <see cref="LowerBoundNewtonRaphson"/>.
        /// </summary>
        /// <param name="regression">The regression to estimate.</param>
        /// 
        public LowerBoundNewtonRaphson(MultinomialLogisticRegression regression)
            : this()
        {
            this.init(regression);
        }

        private void init(MultinomialLogisticRegression regression)
        {
            this.regression = regression;

            this.K = regression.NumberOfOutputs - 1;
            this.M = regression.NumberOfInputs + 1;
            this.parameterCount = this.K * this.M;

            this.solution = regression.Coefficients.Reshape();

            this.xxt = new double[this.M, this.M];
            this.errors = new double[this.K];
            this.output = new double[this.K];

            this.lowerBound = new double[this.parameterCount, this.parameterCount];
            this.gradient = new double[this.parameterCount];

            // Differently from the IRLS iteration, the weight matrix can be fixed
            // as it does not depend on the current coefficients anymore [I - 11/m]

            // TODO: Avoid the multiple allocations in the line below
            this.weights = (-0.5).Multiply(Matrix.Identity(this.K).Subtract(Matrix.Create(this.K, this.K, 1.0 / this.M)));
        }


        /// <summary>
        ///   Runs one iteration of the Lower-Bound Newton-Raphson iteration.
        /// </summary>
        /// <param name="inputs">The input data.</param>
        /// <param name="classes">The outputs associated with each input vector.</param>
        /// <returns>The maximum relative change in the parameters after the iteration.</returns>
        /// 
        [Obsolete("Please use the Learn() method instead.")]
        public double Run(double[][] inputs, int[] classes)
        {
            int old = this.Iterations;
            this.Iterations = 1;
            this.Learn(inputs, classes);
            this.Iterations = old;
            return this.deltas.Max();
        }

        /// <summary>
        ///   Runs one iteration of the Lower-Bound Newton-Raphson iteration.
        /// </summary>
        /// <param name="inputs">The input data.</param>
        /// <param name="outputs">The outputs associated with each input vector.</param>
        /// <returns>The maximum relative change in the parameters after the iteration.</returns>
        /// 
        [Obsolete("Please use the Learn() method instead.")]
        public double Run(double[][] inputs, double[][] outputs)
        {
            int old = this.Iterations;
            this.Iterations = 1;
            this.Learn(inputs, outputs);
            this.Iterations = old;
            return this.deltas.Max();
        }



        private void compute(double[] x, double[] responses)
        {
            double[][] coefficients = this.regression.Coefficients;

            double sum = 1;

            // For each category (except the first)
            for (int j = 0; j < coefficients.Length; j++)
            {
                // Get category coefficients
                double[] c = coefficients[j];

                double logit = c[0]; // intercept

                for (int i = 0; i < x.Length; i++)
                    logit += c[i + 1] * x[i];

                sum += responses[j] = Math.Exp(logit);
            }

            // Normalize the probabilities
            for (int j = 0; j < responses.Length; j++)
                responses[j] /= sum;
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
            return this.Learn(x, Jagged.OneHot(y), weights);
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
            return this.Learn(x, y.ToDouble(), weights);
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
            if (weights != null)
                throw new ArgumentException(global::Accord.Properties.Resources.NotSupportedWeights, "weights");

            // Regress using Lower-Bound Newton-Raphson estimation
            //
            // The main idea is to replace the Hessian matrix with a 
            //   suitable lower bound. Indeed, the Hessian is lower
            //   bounded by a negative definite matrix that does not
            //   even depend on w [Krishnapuram et al].
            //
            //   - http://www.lx.it.pt/~mtf/Krishnapuram_Carin_Figueiredo_Hartemink_2005.pdf
            //

            if (this.regression == null)
                this.init(new MultinomialLogisticRegression(x.Columns(), y.Columns()));

            // Initial definitions and memory allocations
            int N = x.Length;

            double[][] design = x.InsertColumn(value: 1, index: 0);
            double[][] coefficients = this.regression.Coefficients;

            // In the multinomial logistic regression, the objective
            // function is the log-likelihood function l(w). As given
            // by Krishnapuram et al and Böhning, this is a concave 
            // function with Hessian given by:
            //
            //       H(w) = -sum(P(w) - p(w)p(w)')  (x)  xx'
            //      (see referenced paper for proper indices)
            //       
            // In which (x) denotes the Kronecker product. By using
            // the lower bound principle, Krishnapuram has shown that
            // we can replace H(w) with a lower bound approximation B
            // which does not depend on w (eq. 8 on aforementioned paper):
            // 
            //      B = -(1/2) [I - 11/M]  (x)  sum(xx')
            //
            // Thus we can compute and invert this matrix only once.
            //

            do
            {
                if (this.Token.IsCancellationRequested)
                    break;

                // Reset Hessian matrix and gradient
                Array.Clear(this.gradient, 0, this.gradient.Length);

                if (this.UpdateLowerBound)
                    Array.Clear(this.lowerBound, 0, this.lowerBound.Length);

                // For each input sample in the dataset
                for (int i = 0; i < x.Length; i++)
                {
                    // Grab variables related to the sample
                    double[] rx = design[i];
                    double[] ry = y[i];

                    // Compute and estimate outputs
                    this.compute(x[i], this.output);

                    // Compute errors for the sample
                    for (int j = 0; j < this.errors.Length; j++)
                        this.errors[j] = ry[j + 1] - this.output[j];


                    // Compute current gradient and Hessian
                    //   We can take advantage of the block structure of the 
                    //   Hessian matrix and gradient vector by employing the
                    //   Kronecker product. See [Böhning, 1992]

                    // (Re-) Compute error gradient
                    double[] g = Matrix.Kronecker(this.errors, rx);
                    for (int j = 0; j < g.Length; j++)
                        this.gradient[j] += g[j];

                    if (this.UpdateLowerBound)
                    {
                        // Compute xxt matrix
                        for (int k = 0; k < rx.Length; k++)
                            for (int j = 0; j < rx.Length; j++)
                                this.xxt[k, j] = rx[k] * rx[j];

                        // (Re-) Compute weighted "Hessian" matrix 
                        double[,] h = Matrix.Kronecker(this.weights, this.xxt);
                        for (int j = 0; j < this.parameterCount; j++)
                            for (int k = 0; k < this.parameterCount; k++)
                                this.lowerBound[j, k] += h[j, k];
                    }
                }


                if (this.UpdateLowerBound)
                {
                    this.UpdateLowerBound = false;

                    // Decompose to solve the linear system. Usually the Hessian will
                    // be invertible and LU will succeed. However, sometimes the Hessian
                    // may be singular and a Singular Value Decomposition may be needed.

                    // The SVD is very stable, but is quite expensive, being on average
                    // about 10-15 times more expensive than LU decomposition. There are
                    // other ways to avoid a singular Hessian. For a very interesting 
                    // reading on the subject, please see:
                    //
                    //  - Jeff Gill & Gary King, "What to Do When Your Hessian Is Not Invertible",
                    //    Sociological Methods & Research, Vol 33, No. 1, August 2004, 54-87.
                    //    Available in: http://gking.harvard.edu/files/help.pdf
                    //

                    // Moreover, the computation of the inverse is optional, as it will
                    // be used only to compute the standard errors of the regression.


                    this.decomposition = new SingularValueDecomposition(this.lowerBound);
                    this.deltas = this.decomposition.Solve(this.gradient);
                }
                else
                {
                    this.deltas = this.decomposition.Solve(this.gradient);
                }


                this.previous = coefficients.Reshape();

                // Update coefficients using the calculated deltas
                for (int i = 0, k = 0; i < coefficients.Length; i++)
                    for (int j = 0; j < coefficients[i].Length; j++)
                        coefficients[i][j] -= this.deltas[k++];

                this.solution = coefficients.Reshape();

                // Return the relative maximum parameter change
                for (int i = 0; i < this.deltas.Length; i++)
                    this.deltas[i] = Math.Abs(this.deltas[i]) / Math.Abs(this.previous[i]);

                this.convergence.NewValue = this.deltas.Max();

                if (this.Token.IsCancellationRequested)
                    break;

            } while (!this.convergence.HasConverged);

            if (this.computeStandardErrors)
            {
                // Grab the regression information matrix
                double[,] inverse = this.decomposition.Inverse();

                // Calculate coefficients' standard errors
                double[][] standardErrors = this.regression.StandardErrors;
                for (int i = 0, k = 0; i < standardErrors.Length; i++)
                    for (int j = 0; j < standardErrors[i].Length; j++, k++)
                        standardErrors[i][j] = Math.Sqrt(Math.Abs(inverse[k, k]));
            }

            return this.regression;
        }
    }
}
