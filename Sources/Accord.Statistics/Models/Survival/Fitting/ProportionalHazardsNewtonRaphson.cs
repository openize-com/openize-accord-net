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

namespace Openize.Accord.Statistics.Models.Survival.Fitting
{
    using System;
    using System.Threading;
    using Accord.MachineLearning.Learning;
    using Distributions;
    using Openize.Accord.Math.Core;
    using Openize.Accord.Math.Matrix;
    using Math.Convergence;
    using Math.Decompositions;
    using Math.Decompositions.Base;
    using Openize.Accord.Core;
    using Openize.Accord.Core.Exceptions;
    using Openize.Accord.Math.Accord.Statistics;
    using Openize.Accord.Statistics.Analysis;
    using Openize.Accord.Statistics.Distributions.Fitting;
    using Openize.Accord.Statistics.Distributions.Univariate.Continuous;

    /// <summary>
    ///   Newton-Raphson learning updates for Cox's Proportional Hazards models.
    /// </summary>
    /// 
    /// <example>
    ///   <code source="Unit Tests\Accord.Tests.Statistics\Models\Regression\CoxProportionalHazardsTest.cs" region="doc_learn" />
    /// </example>
    /// 
    /// <seealso cref="ProportionalHazards"/>
    /// <seealso cref="ProportionalHazardsAnalysis"/>
    /// 
#pragma warning disable 612, 618
    public class ProportionalHazardsNewtonRaphson :
        ISupervisedLearning<ProportionalHazards, Tuple<double[], double>, int>,
        ISupervisedLearning<ProportionalHazards, Tuple<double[], double>, bool>,
        ISurvivalFitting, IConvergenceLearning
#pragma warning disable 612, 618
    {
        [NonSerialized]
        CancellationToken token = new CancellationToken();

        private ProportionalHazards regression;
        private int parameterCount;

        private double[,] hessian;
        private double[] gradient;

        private double[] partialGradient;
        private double[,] partialHessian;

        private bool computeStandardErrors = true;
        private bool computeBaselineFunction = true;
        private bool normalize = true;

        private ISolverMatrixDecomposition<double> decomposition;
        private RelativeParameterConvergence convergence;


        /// <summary>
        ///   Gets or sets the maximum absolute parameter change detectable
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
        ///   Please use MaxIterations instead.
        /// </summary>
        [Obsolete("Please use MaxIterations instead.")]
        public int Iterations
        {
            get { return this.convergence.MaxIterations; }
            set { this.convergence.MaxIterations = value; }
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
        ///   Gets or sets the hazard estimator that should be used by the
        ///   proportional hazards learning algorithm. Default is to use
        ///   <see cref="HazardEstimator.BreslowNelsonAalen"/>.
        /// </summary>
        /// 
        public HazardEstimator Estimator { get; set; }

        /// <summary>
        ///   Gets or sets the ties handling method to be used by the
        ///   proportional hazards learning algorithm. Default is to use
        ///   <see cref="HazardTiesMethod.Efron"/>'s method.
        /// </summary>
        /// 
        public HazardTiesMethod Ties { get; set; }


        /// <summary>
        ///   Gets the previous values for the coefficients which were
        ///   in place before the last learning iteration was performed.
        /// </summary>
        /// 
        public double[] Previous { get { return this.convergence.OldValues; } }

        /// <summary>
        ///   Gets the current values for the coefficients.
        /// </summary>
        /// 
        public double[] Solution { get { return this.regression.Coefficients; } }

        /// <summary>
        ///   Gets the Hessian matrix computed in 
        ///   the last Newton-Raphson iteration.
        /// </summary>
        /// 
        public double[,] Hessian { get { return this.hessian; } }

        /// <summary>
        ///   Gets the Gradient vector computed in
        ///   the last Newton-Raphson iteration.
        /// </summary>
        /// 
        public double[] Gradient { get { return this.gradient; } }

        /// <summary>
        ///   Gets the total number of parameters in the model.
        /// </summary>
        /// 
        public int Parameters { get { return this.parameterCount; } }

        /// <summary>
        ///   Gets or sets the regression model being learned.
        /// </summary>
        /// 
        public ProportionalHazards Model
        {
            get { return this.regression; }
            set
            {
                this.regression = value;
                if (this.regression != null)
                    this.init(this.regression);
            }
        }

        /// <summary>
        ///   Gets or sets a value indicating whether standard
        ///   errors should be computed at the end of the next 
        ///   iterations.
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
        ///   Gets or sets a value indicating whether an estimate
        ///   of the baseline hazard function should be computed
        ///   at the end of the next iterations.
        /// </summary>
        /// <value>
        /// 	<c>true</c> to compute the baseline function; otherwise, <c>false</c>.
        /// </value>
        /// 
        public bool ComputeBaselineFunction
        {
            get { return this.computeBaselineFunction; }
            set { this.computeBaselineFunction = value; }
        }

        /// <summary>
        ///   Gets or sets a value indicating whether the Cox model should
        ///   be computed using the mean-centered version of the covariates.
        ///   Default is true.
        /// </summary>
        /// 
        public bool Normalize
        {
            get { return this.normalize; }
            set { this.normalize = value; }
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
        ///   Gets or sets the smoothing factor used to avoid numerical
        ///   problems in the beginning of the training. Default is 0.1.
        /// </summary>
        /// 
        public double Lambda { get; set; }

        /// <summary>
        ///   Constructs a new Newton-Raphson learning algorithm
        ///   for Cox's Proportional Hazards models.
        /// </summary>
        /// 
        public ProportionalHazardsNewtonRaphson()
        {
            this.convergence = new RelativeParameterConvergence()
            {
                Iterations = 0,
                Tolerance = 1e-5
            };

            this.Estimator = HazardEstimator.BreslowNelsonAalen;
            this.Ties = HazardTiesMethod.Efron;
            this.Lambda = 0.1;
            this.Token = new CancellationToken();
        }

        /// <summary>
        ///   Constructs a new Newton-Raphson learning algorithm
        ///   for Cox's Proportional Hazards models.
        /// </summary>
        /// 
        /// <param name="hazards">The model to estimate.</param>
        /// 
        public ProportionalHazardsNewtonRaphson(ProportionalHazards hazards)
            : this()
        {
            this.init(hazards);
        }

        private void init(ProportionalHazards hazards)
        {
            this.regression = hazards;
            this.parameterCount = hazards.Coefficients.Length;

            this.hessian = new double[this.parameterCount, this.parameterCount];
            this.gradient = new double[this.parameterCount];

            this.partialHessian = new double[this.parameterCount, this.parameterCount];
            this.partialGradient = new double[this.parameterCount];
        }



        /// <summary>
        ///   Runs the Newton-Raphson update for Cox's hazards learning until convergence.
        /// </summary>
        /// 
        /// <param name="inputs">The input data.</param>
        /// <param name="time">The time-to-event for the training samples.</param>
        /// 
        /// <returns>The maximum relative change in the parameters after the iteration.</returns>
        /// 
        [Obsolete("Please use Learn(x, y) instead.")]
        public double Run(double[][] inputs, double[] time)
        {
            var censor = new SurvivalOutcome[time.Length];

            Debug.Assert(censor[0] == SurvivalOutcome.Failed);

            return this.Run(inputs, time, censor);
        }

        /// <summary>
        ///   Runs the Newton-Raphson update for Cox's hazards learning until convergence.
        /// </summary>
        /// 
        /// <param name="inputs">The input data.</param>
        /// <param name="censor">The output (event) associated with each input vector.</param>
        /// <param name="time">The time-to-event for the non-censored training samples.</param>
        /// 
        /// <returns>The maximum relative change in the parameters after the iteration.</returns>
        /// 
        [Obsolete("Please use Learn(x, y) instead.")]
        public double Run(double[][] inputs, double[] time, int[] censor)
        {
            return this.Run(inputs, time, censor.To<SurvivalOutcome[]>());
        }

        /// <summary>
        ///   Runs the Newton-Raphson update for Cox's hazards learning until convergence.
        /// </summary>
        /// 
        /// <param name="inputs">The input data.</param>
        /// <param name="censor">The output (event) associated with each input vector.</param>
        /// <param name="time">The time-to-event for the non-censored training samples.</param>
        /// 
        /// <returns>The maximum relative change in the parameters after the iteration.</returns>
        /// 
        [Obsolete("Please use Learn(x, y) instead.")]
        public double Run(double[][] inputs, double[] time, SurvivalOutcome[] censor)
        {
            this.innerLearn(inputs, time, censor, null);
            return this.regression.GetPartialLogLikelihood(inputs, time, censor);
        }

        /// <summary>
        ///   Runs the Newton-Raphson update for Cox's hazards learning until convergence.
        /// </summary>
        /// 
        /// <param name="censor">The output (event) associated with each input vector.</param>
        /// <param name="time">The time-to-event for the non-censored training samples.</param>
        /// 
        /// <returns>The maximum relative change in the parameters after the iteration.</returns>
        /// 
        [Obsolete("Please use Learn(x, y) instead.")]
        public double Run(double[] time, int[] censor)
        {
            return this.Run(time, censor.To<SurvivalOutcome[]>());
        }

        /// <summary>
        ///   Runs the Newton-Raphson update for Cox's hazards learning until convergence.
        /// </summary>
        /// 
        /// <param name="censor">The output (event) associated with each input vector.</param>
        /// <param name="time">The time-to-event for the non-censored training samples.</param>
        /// 
        /// <returns>The maximum relative change in the parameters after the iteration.</returns>
        /// 
        [Obsolete("Please use Learn(x, y) instead.")]
        public double Run(double[] time, SurvivalOutcome[] censor)
        {
            if (time.Length != censor.Length)
            {
                throw new DimensionMismatchException("time",
                    "The time and output vector must have the same length.");
            }

            // Sort data by time to accelerate performance
            EmpiricalHazardDistribution.Sort(ref time, ref censor);

            this.createBaseline(time, censor);

            return this.regression.GetPartialLogLikelihood(time, censor);
        }

        private void createBaseline(double[] time, SurvivalOutcome[] censor, double[] output = null)
        {
            if (this.regression.BaselineHazard == null)
                return;

            var hazard = this.regression.BaselineHazard as IFittableDistribution<double, EmpiricalHazardOptions>;
            if (hazard != null)
            {
                // Compute an estimate of the cumulative Hazard
                //   function using the Nelson-Aalen estimator
                hazard.Fit(time, output, new EmpiricalHazardOptions()
                {
                    Outcome = censor,
                    Estimator = this.Estimator,
                    Ties = this.Ties
                });
                return;
            }

            var survival = this.regression.BaselineHazard as IFittableDistribution<double, SurvivalOptions>;
            if (survival != null)
            {
                // Compute an estimate of the cumulative Hazard
                //   function using the Kaplan-Meier estimator
                survival.Fit(time, new SurvivalOptions()
                {
                    Outcome = censor,
                });
            }
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
        public ProportionalHazards Learn(Tuple<double[], double>[] x, bool[] y, double[] weights = null)
        {
            double[][] inputs = x.Apply(x_i => x_i.Item1);
            double[] time = x.Apply(x_i => x_i.Item2);
            return this.Learn(inputs, time, y.To<SurvivalOutcome[]>(), weights);
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
        public ProportionalHazards Learn(Tuple<double[], double>[] x, int[] y, double[] weights = null)
        {
            double[][] inputs = x.Apply(x_i => x_i.Item1);
            double[] time = x.Apply(x_i => x_i.Item2);
            return this.Learn(inputs, time, y.To<SurvivalOutcome[]>(), weights);
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
        public ProportionalHazards Learn(Tuple<double[], double>[] x, SurvivalOutcome[] y, double[] weights = null)
        {
            double[][] inputs = x.Apply(x_i => x_i.Item1);
            double[] time = x.Apply(x_i => x_i.Item2);
            SurvivalOutcome[] censor = y;
            return this.Learn(inputs, time, censor, weights);
        }

        /// <summary>
        /// Learns a model that can map the given inputs to the given outputs.
        /// </summary>
        /// <param name="inputs">The model inputs.</param>
        /// <param name="censor">The output (event) associated with each input vector.</param>
        /// <param name="time">The time-to-event for the non-censored training samples.</param>
        /// <param name="weights">The weight of importance for each input-output pair (if supported by the learning algorithm).</param>
        /// <returns>
        /// A model that has learned how to produce <paramref name="censor" /> given <paramref name="inputs" /> and <paramref name="time" />.
        /// </returns>
        public ProportionalHazards Learn(double[][] inputs, double[] time, int[] censor, double[] weights = null)
        {
            return this.Learn(inputs, time, censor.To<SurvivalOutcome[]>(), weights);
        }

        /// <summary>
        /// Learns a model that can map the given inputs to the given outputs.
        /// </summary>
        /// <param name="inputs">The model inputs.</param>
        /// <param name="censor">The output (event) associated with each input vector.</param>
        /// <param name="time">The time-to-event for the non-censored training samples.</param>
        /// <param name="weights">The weight of importance for each input-output pair (if supported by the learning algorithm).</param>
        /// <returns>
        /// A model that has learned how to produce <paramref name="censor" /> given <paramref name="inputs" /> and <paramref name="time" />.
        /// </returns>
        public ProportionalHazards Learn(double[][] inputs, double[] time, SurvivalOutcome[] censor, double[] weights = null)
        {
            var regression = this.innerLearn(inputs, time, censor, weights);

            // Disable deprecated mechanism
            regression.Intercept = -regression.Coefficients.Dot(regression.Offsets);
            regression.Offsets.Clear();

            return regression;
        }

        private ProportionalHazards innerLearn(double[][] inputs, double[] time, SurvivalOutcome[] censor, double[] weights)
        {
            if (weights != null)
                throw new ArgumentException(global::Accord.Properties.Resources.NotSupportedWeights, "weights");

            if (inputs.Length != time.Length || time.Length != censor.Length)
            {
                throw new DimensionMismatchException("time",
                    "The inputs, time and output vector must have the same length.");
            }

            if (this.regression == null)
                this.init(new ProportionalHazards(inputs.Columns()));

            // Sort data by time to accelerate performance
            EmpiricalHazardDistribution.Sort(ref time, ref censor, ref inputs);


            var means = new double[this.parameterCount];
            var sdev = new double[this.parameterCount];
            for (int i = 0; i < sdev.Length; i++)
                sdev[i] = 1;

            if (this.normalize)
            {
                // Store means as regression centers
                means = inputs.Mean(dimension: 0);
                for (int i = 0; i < means.Length; i++)
                    this.regression.Offsets[i] = means[i];

                // Convert to unit scores for increased accuracy
                sdev = Measures.StandardDeviation(inputs);
                inputs = Elementwise.Divide(inputs.Subtract(means, 0), sdev, 0);

                for (int i = 0; i < this.regression.Coefficients.Length; i++)
                    this.regression.Coefficients[i] *= sdev[i];
            }

            // Compute actual outputs
            var output = new double[inputs.Length];
            for (int i = 0; i < output.Length; i++)
            {
                double sum = 0;
                for (int j = 0; j < this.regression.Coefficients.Length; j++)
                    sum += this.regression.Coefficients[j] * inputs[i][j];
                output[i] = Math.Exp(sum);
            }

            // Compute ties
            int[] ties = new int[inputs.Length];
            for (int i = 0; i < inputs.Length; i++)
                for (int j = 0; j < time.Length; j++)
                    if (time[j] == time[i]) ties[i]++;

            if (this.parameterCount == 0)
            {
                this.createBaseline(time, censor, output);
                return this.regression;
            }

            this.CurrentIteration = 0;
            double smooth = this.Lambda;

            do
            {
                if (this.Token.IsCancellationRequested)
                    break;

                // learning iterations until convergence
                // or maximum number of iterations reached

                this.CurrentIteration++;

                // Reset Hessian matrix and gradient
                Array.Clear(this.gradient, 0, this.gradient.Length);
                Array.Clear(this.hessian, 0, this.hessian.Length);

                // For each observation instance
                for (int i = 0; i < inputs.Length; i++)
                {
                    // Check if we should censor
                    if (censor[i] == SurvivalOutcome.Censored)
                        continue;

                    // Compute partials 
                    double den = 0;
                    Array.Clear(this.partialGradient, 0, this.partialGradient.Length);
                    Array.Clear(this.partialHessian, 0, this.partialHessian.Length);

                    for (int j = 0; j < inputs.Length; j++)
                    {
                        if (time[j] >= time[i])
                            den += output[j];
                    }

                    for (int j = 0; j < inputs.Length; j++)
                    {
                        if (time[j] >= time[i])
                        {
                            // Compute partial gradient
                            for (int k = 0; k < this.partialGradient.Length; k++)
                                this.partialGradient[k] += inputs[j][k] * output[j] / den;

                            // Compute partial Hessian
                            for (int ii = 0; ii < inputs[j].Length; ii++)
                                for (int jj = 0; jj < inputs[j].Length; jj++)
                                    this.partialHessian[ii, jj] += inputs[j][ii] * inputs[j][jj] * output[j] / den;
                        }
                    }

                    // Compute gradient vector
                    for (int j = 0; j < this.gradient.Length; j++)
                        this.gradient[j] += inputs[i][j] - this.partialGradient[j];

                    // Compute Hessian matrix
                    for (int j = 0; j < this.partialGradient.Length; j++)
                        for (int k = 0; k < this.partialGradient.Length; k++)
                            this.hessian[j, k] -= this.partialHessian[j, k] - this.partialGradient[j] * this.partialGradient[k];
                }


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

                this.decomposition = new SingularValueDecomposition(this.hessian);
                double[] deltas = this.decomposition.Solve(this.gradient);

                if (this.convergence.Iterations > 0 || this.convergence.Tolerance > 0)
                {
                    // Update coefficients using the calculated deltas
                    for (int i = 0; i < this.regression.Coefficients.Length; i++)
                        this.regression.Coefficients[i] -= smooth * deltas[i];
                }

                smooth += this.Lambda;
                if (smooth > 1)
                    smooth = 1;

                // Check relative maximum parameter change
                this.convergence.NewValues = this.regression.Coefficients;


                if (this.convergence.HasDiverged)
                {
                    // Restore previous coefficients
                    for (int i = 0; i < this.regression.Coefficients.Length; i++)
                        this.regression.Coefficients[i] = this.convergence.OldValues[i];
                }

                // Recompute current outputs
                for (int i = 0; i < output.Length; i++)
                {
                    double sum = 0;
                    for (int j = 0; j < this.regression.Coefficients.Length; j++)
                        sum += this.regression.Coefficients[j] * inputs[i][j];
                    output[i] = Math.Exp(sum);
                }

                if (this.Token.IsCancellationRequested)
                    return this.regression;

            } while (!this.convergence.HasConverged);


            for (int i = 0; i < this.regression.Coefficients.Length; i++)
                this.regression.Coefficients[i] /= sdev[i];

            if (this.computeStandardErrors)
            {
                // Grab the regression information matrix
                double[,] inverse = this.decomposition.Inverse();

                // Calculate coefficients' standard errors
                double[] standardErrors = this.regression.StandardErrors;
                for (int i = 0; i < standardErrors.Length; i++)
                    standardErrors[i] = Math.Sqrt(Math.Abs(inverse[i, i])) / sdev[i];
            }

            if (this.computeBaselineFunction)
                this.createBaseline(time, censor, output);

            return this.regression;
        }
    }
}
