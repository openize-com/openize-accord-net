﻿// Accord Math Library
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

namespace Openize.Accord.Math.Optimization.Unconstrained.Least_Squares
{
    using System.Threading.Tasks;
    using Decompositions;
    using Openize.Accord.Math.Core;
    using Openize.Accord.Math.Matrix;
    using Openize.Accord.Math.Accord.Statistics;
    using Openize.Accord.Math.Differentiation;
    using Openize.Accord.Math.Optimization.Base;

    /// <summary>
    ///   Gauss-Newton algorithm for solving Least-Squares problems.
    /// </summary>
    /// 
    /// <remarks>
    ///   This class isn't suitable for most real-world problems. Instead, this class
    ///   is intended to be use as a baseline comparison to help debug and check other
    ///   optimization methods, such as <see cref="LevenbergMarquardt"/>.
    /// </remarks>
    /// 
    /// <example>
    /// <para>
    ///   While it is possible to use the <see cref="GaussNewton"/> class as a standalone
    ///   method for solving least squares problems, this class is intended to be used as
    ///   a strategy for NonlinearLeastSquares, as shown in the example below:</para>
    ///   <code source="Unit Tests\Accord.Tests.Statistics\Models\Regression\NonlinearLeastSquaresTest.cs" region="doc_learn_gn" lang="cs"/>
    ///   <code source="Unit Tests\Accord.Tests.Statistics.VB\Models\Regression\NonlinearLeastSquaresTest.vb" region="doc_learn_gn" lang="vb"/>
    ///   
    /// <para>
    ///   However, as mentioned above it is also possible to use <see cref="GaussNewton"/> as
    ///   a standalone class, as shown in the example below:</para>
    ///   <code source="Unit Tests\Accord.Tests.Statistics\Models\Regression\GaussNewtonTest.cs" region="doc_minimize"/>
    /// </example>
    /// 
    /// <seealso cref="LevenbergMarquardt"/>
    /// <seealso cref="FiniteDifferences"/>
    /// 
    public class GaussNewton : BaseLeastSquaresMethod, ILeastSquaresMethod, IConvergenceLearning
    {

        private double[] gradient; // this is just a cached vector to avoid memory allocations
        private double[] errors;
        private double[] deltas;

        private double[][] hessian;
        private double[][] jacobian;

        private JaggedSingularValueDecomposition decomposition;

        /// <summary>
        ///   Gets the approximate Hessian matrix of second derivatives
        ///   created during the last algorithm iteration.
        /// </summary>
        /// 
        /// <remarks>
        /// <para>
        ///   Please note that this value is actually just an approximation to the
        ///   actual Hessian matrix using the outer Jacobian approximation (H ~ J'J).
        /// </para>
        /// </remarks>
        /// 
        public double[][] Hessian
        {
            get { return this.hessian; }
        }

        /// <summary>
        ///   Gets the vector of residuals computed in the last iteration.
        ///   The residuals are computed as <c>(y - f(w, x))</c>, in which 
        ///   <c>y</c> are the expected output values, and <c>f</c> is the
        ///   parameterized model function.
        /// </summary>
        /// 
        public double[] Residuals
        {
            get { return this.errors; }
        }

        /// <summary>
        ///   Gets the Jacobian matrix of first derivatives computed in the
        ///   last iteration.
        /// </summary>
        /// 
        public double[][] Jacobian
        {
            get { return this.jacobian; }
        }

        /// <summary>
        ///   Gets the vector of coefficient updates computed in the last iteration.
        /// </summary>
        /// 
        public double[] Deltas
        {
            get { return this.deltas; }
        }

        /// <summary>
        ///   Gets standard error for each parameter in the solution.
        /// </summary>
        /// 
        public double[] StandardErrors
        {
            get { return this.decomposition.Inverse().Diagonal().Sqrt(); }
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="GaussNewton" /> class.
        /// </summary>
        /// 
        public GaussNewton()
        {

        }

        /// <summary>
        ///   Initializes a new instance of the <see cref="GaussNewton"/> class.
        /// </summary>
        /// 
        /// <param name="parameters">The number of variables (free parameters)
        ///   in the objective function.</param>
        /// 
        public GaussNewton(int parameters)
            : this()
        {
            this.NumberOfParameters = parameters;
        }

        /// <summary>
        /// This method should be implemented by child classes to initialize
        /// their fields once the <see cref="BaseLeastSquaresMethod.NumberOfParameters" /> is known.
        /// </summary>
        /// 
        protected override void Initialize()
        {
            this.hessian = Jagged.Zeros(this.NumberOfParameters, this.NumberOfParameters);
            this.gradient = new double[this.NumberOfParameters];
            this.jacobian = new double[this.NumberOfParameters][];
        }


        /// <summary>
        ///   Attempts to find the best values for the parameter vector
        ///   minimizing the discrepancy between the generated outputs
        ///   and the expected outputs for a given set of input data.
        /// </summary>
        /// 
        /// <param name="inputs">A set of input data.</param>
        /// <param name="outputs">The values associated with each 
        ///   vector in the <paramref name="inputs"/> data.</param>
        /// 
        public double Minimize(double[][] inputs, double[] outputs)
        {
            this.Convergence.CurrentIteration = 0;

            do
            {
                this.Convergence.NewValue = this.iterate(inputs, outputs);
            } while (!this.Convergence.HasConverged);

            return this.Value = this.Convergence.NewValue;
        }

        private double iterate(double[][] inputs, double[] outputs)
        {
            if (this.errors == null || inputs.Length != this.errors.Length)
            {
                this.errors = new double[inputs.Length];
                for (int i = 0; i < this.jacobian.Length; i++)
                    this.jacobian[i] = new double[inputs.Length];
            }

            for (int i = 0; i < inputs.Length; i++)
                this.errors[i] = outputs[i] - this.Function(this.Solution, inputs[i]);

            if (this.ParallelOptions.MaxDegreeOfParallelism == 1)
            {
                for (int i = 0; i < inputs.Length; i++)
                {
                    this.Gradient(this.Solution, inputs[i], result: this.gradient);

                    for (int j = 0; j < this.gradient.Length; j++)
                        this.jacobian[j][i] = -this.gradient[j];

                    if (this.Token.IsCancellationRequested)
                        break;
                }
            }
            else
            {
                Parallel.For(0, inputs.Length, this.ParallelOptions,

                    () => new double[this.NumberOfParameters],

                    (i, state, grad) =>
                    {
                        this.Gradient(this.Solution, inputs[i], result: grad);

                        for (int j = 0; j < grad.Length; j++)
                            this.jacobian[j][i] = -grad[j];

                        return grad;
                    },

                    (grad) => { }
                );
            }


            // Compute error gradient using Jacobian
            this.jacobian.Dot(this.errors, result: this.gradient);

            // Compute Quasi-Hessian Matrix approximation
            this.jacobian.DotWithTransposed(this.jacobian, result: this.hessian);

            this.decomposition = new JaggedSingularValueDecomposition(this.hessian,
                computeLeftSingularVectors: true, computeRightSingularVectors: true, autoTranspose: true);

            this.deltas = this.decomposition.Solve(this.gradient);

            for (int i = 0; i < this.deltas.Length; i++)
                this.Solution[i] -= this.deltas[i];

            return this.ComputeError(inputs, outputs);
        }

    }
}

