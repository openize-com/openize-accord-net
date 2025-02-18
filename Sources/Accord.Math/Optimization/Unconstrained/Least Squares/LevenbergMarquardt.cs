// Accord Math Library
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
    using System;
    using System.Threading.Tasks;
    using Decompositions;
    using Openize.Accord.Math.Core;
    using Openize.Accord.Math.Matrix;
    using Openize.Accord.Math.Accord.Statistics;
    using Openize.Accord.Math.Differentiation;
    using Openize.Accord.Math.Optimization.Base;
    using Vector = Openize.Accord.Math.Vector.Vector;

    /// <summary>
    ///   Levenberg-Marquardt algorithm for solving Least-Squares problems.
    /// </summary>
    /// 
    /// <example>
    /// <para>
    ///   While it is possible to use the <see cref="LevenbergMarquardt"/> class as a standalone
    ///   method for solving least squares problems, this class is intended to be used as
    ///   a strategy for NonlinearLestSquares, as shown in the example below:</para>
    ///   <code source="Unit Tests\Accord.Tests.Statistics\Models\Regression\NonlinearLeastSquaresTest.cs" region="doc_learn_lm" lang="cs"/>
    ///   <code source="Unit Tests\Accord.Tests.Statistics.VB\Models\Regression\NonlinearLeastSquaresTest.vb" region="doc_learn_lm" lang="vb"/>
    ///   
    /// <para>
    ///   However, as mentioned above it is also possible to use <see cref="LevenbergMarquardt"/> 
    ///   as a standalone class, as shown in the example below:</para>
    ///   <code source="Unit Tests\Accord.Tests.Statistics\Models\Regression\LevenbergMarquardtTest.cs" region="doc_minimize"/>
    /// </example>
    /// 
    /// <seealso cref="GaussNewton"/>
    /// <seealso cref="FiniteDifferences"/>
    /// 
    public class LevenbergMarquardt : BaseLeastSquaresMethod, ILeastSquaresMethod, IConvergenceLearning
    {
        private const double lambdaMax = 1e25;
        private double eps = 1e-12;

        // Levenberg-Marquardt variables
        private double[][] jacobian;
        private double[][] hessian;

        private double[] diagonal;
        private double[] gradient;
        private double[] weights;
        private double[] deltas;
        private double[] errors;


        // Levenberg damping factor
        private double lambda = 0.1;

        // The amount the damping factor is adjusted
        // when searching the minimum error surface
        private double v = 10.0;


        private int blocks = 1;
        private int outputCount = 1;

        JaggedCholeskyDecomposition decomposition;



        /// <summary>
        ///   Levenberg's damping factor, also known as lambda.
        /// </summary>
        /// 
        /// <remarks><para>The value determines speed of learning.</para>
        /// 
        /// <para>Default value is <b>0.1</b>.</para>
        /// </remarks>
        ///
        public double LearningRate
        {
            get { return this.lambda; }
            set { this.lambda = value; }
        }

        /// <summary>
        ///   Learning rate adjustment. 
        /// </summary>
        /// 
        /// <remarks><para>The value by which the learning rate
        /// is adjusted when searching for the minimum cost surface.</para>
        /// 
        /// <para>Default value is <b>10</b>.</para>
        /// </remarks>
        ///
        public double Adjustment
        {
            get { return this.v; }
            set { this.v = value; }
        }


        /// <summary>
        ///   Gets or sets the number of blocks to divide the 
        ///   Jacobian matrix in the Hessian calculation to
        ///   preserve memory. Default is 1.
        /// </summary>
        /// 
        public int Blocks
        {
            get { return this.blocks; }
            set { this.blocks = value; }
        }

        /// <summary>
        ///   Gets or sets a small epsilon value to be added to the
        ///   diagonal of the Hessian matrix. Default is 1e-12.
        /// </summary>
        /// 
        public double Epsilon
        {
            get { return this.eps; }
            set { this.eps = value; }
        }

        /// <summary>
        ///   Gets the approximate Hessian matrix of second derivatives 
        ///   generated in the last algorithm iteration. The Hessian is 
        ///   stored in the upper triangular part of this matrix. See 
        ///   remarks for details.
        ///   </summary>
        ///   
        /// <remarks>
        /// <para>
        ///   The Hessian needs only be upper-triangular, since
        ///   it is symmetric. The Cholesky decomposition will
        ///   make use of this fact and use the lower-triangular
        ///   portion to hold the decomposition, conserving memory</para>
        ///   
        /// <para>
        ///   Thus said, this property will hold the Hessian matrix
        ///   in the upper-triangular part of this matrix, and store
        ///   its Cholesky decomposition on its lower triangular part.</para>
        ///   
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
        ///   Gets standard error for each parameter in the solution.
        /// </summary>
        /// 
        public double[] StandardErrors
        {
            get { return this.decomposition.InverseDiagonal().Sqrt(); }
        }



        /// <summary>
        /// Initializes a new instance of the <see cref="LevenbergMarquardt" /> class.
        /// </summary>
        public LevenbergMarquardt()
        {
        }

        /// <summary>
        ///   Initializes a new instance of the <see cref="LevenbergMarquardt"/> class.
        /// </summary>
        /// 
        /// <param name="parameters">The number of free parameters in the optimization problem.</param>
        /// 
        public LevenbergMarquardt(int parameters)
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
            this.weights = new double[this.NumberOfParameters];
            this.diagonal = new double[this.NumberOfParameters];
            this.gradient = new double[this.NumberOfParameters];

            this.jacobian = new double[this.NumberOfParameters][];
            this.hessian = Jagged.Zeros(this.NumberOfParameters, this.NumberOfParameters);
            for (int i = 0; i < this.hessian.Length; i++)
                this.hessian[i] = new double[this.NumberOfParameters];
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
            if (this.NumberOfParameters == 0)
                throw new InvalidOperationException("Please set the NumberOfVariables property first.");

            // Divide the problem into blocks. Instead of computing
            // a single Jacobian and a single error vector, we will
            // be computing multiple Jacobians for smaller problems
            // and then sum all blocks into the final Hessian matrix
            // and gradient vector.

            int blockSize = inputs.Length / this.Blocks;
            int finalBlock = inputs.Length % this.Blocks;
            int jacobianSize = blockSize * this.outputCount;

            // Re-allocate the partial Jacobian matrix only if needed
            if (this.jacobian[0] == null || this.jacobian[0].Length < jacobianSize)
            {
                for (int i = 0; i < this.jacobian.Length; i++)
                    this.jacobian[i] = new double[jacobianSize];
            }

            // Re-allocate error vector only if needed
            if (this.errors == null || this.errors.Length < jacobianSize)
                this.errors = new double[jacobianSize];

            this.Convergence.CurrentIteration = 0;

            do
            {
                this.Convergence.NewValue = this.iterate(inputs, outputs, blockSize, finalBlock, jacobianSize);
            } while (!this.Convergence.HasConverged);


            return this.Value = this.Convergence.NewValue;
        }

        private double iterate(double[][] inputs, double[] outputs, int blockSize, int finalBlock, int jacobianSize)
        {
            double sumOfSquaredErrors = 0;

            // Set upper triangular Hessian to zero
            for (int i = 0; i < this.hessian.Length; i++)
                Array.Clear(this.hessian[i], i, this.hessian.Length - i);

            // Set Gradient vector to zero
            Array.Clear(this.gradient, 0, this.gradient.Length);

            // For each block
            for (int s = 0; s <= this.Blocks; s++)
            {
                if (s == this.Blocks && finalBlock == 0)
                    continue;

                int B = (s == this.Blocks) ? finalBlock : blockSize;
                int[] block = Vector.Range(s * blockSize, s * blockSize + B);

                // Compute the partial residuals vector
                sumOfSquaredErrors += this.computeErrors(inputs, outputs, block);

                // Compute the partial Jacobian
                this.computeJacobian(inputs, block);

                if (Double.IsNaN(sumOfSquaredErrors))
                {
                    throw new ArithmeticException("Error calculation has produced a non-finite number."
                        + " Please make sure that there are no constant columns in the input data.");
                }


                // Compute error gradient using Jacobian
                for (int i = 0; i < this.jacobian.Length; i++)
                {
                    double sum = 0;
                    for (int j = 0; j < jacobianSize; j++)
                        sum += this.jacobian[i][j] * this.errors[j];
                    this.gradient[i] += sum;
                }


                // Compute Quasi-Hessian Matrix approximation
                //  using the outer product Jacobian (H ~ J'J)
                //
                Parallel.For(0, this.jacobian.Length, this.ParallelOptions, i =>
                {
                    double[] ji = this.jacobian[i];
                    double[] hi = this.hessian[i];

                    for (int j = i; j < hi.Length; j++)
                    {
                        double[] jj = this.jacobian[j];

                        double sum = 0;
                        for (int k = 0; k < jj.Length; k++)
                            sum += ji[k] * jj[k];

                        // The Hessian need only be upper-triangular, since
                        // it is symmetric. The Cholesky decomposition will
                        // make use of this fact and use the lower-triangular
                        // portion to hold the decomposition, conserving memory.

                        hi[j] += 2 * sum;
                    }
                });
            }


            // Store the Hessian's diagonal for future computations. The
            // diagonal will be destroyed in the decomposition, so it can
            // still be updated on every iteration by restoring this copy.
            //
            for (int i = 0; i < this.hessian.Length; i++)
                this.diagonal[i] = this.hessian[i][i] + this.eps;

            // Create the initial weights vector
            for (int i = 0; i < this.Solution.Length; i++)
                this.weights[i] = this.Solution[i];


            // Define the objective function:
            double objective = sumOfSquaredErrors;
            double current = objective + 1.0;


            // Begin of the main Levenberg-Marquardt method
            this.lambda /= this.v;

            // We'll try to find a direction with less error
            //  (or where the objective function is smaller)
            while (current >= objective && this.lambda < lambdaMax)
            {
                if (this.Token.IsCancellationRequested)
                    break;

                this.lambda *= this.v;

                // Update diagonal (Levenberg-Marquardt)
                for (int i = 0; i < this.diagonal.Length; i++)
                    this.hessian[i][i] = this.diagonal[i] * (1 + this.lambda);


                // Decompose to solve the linear system. The Cholesky decomposition
                // is done in place, occupying the Hessian's lower-triangular part.
                this.decomposition = new JaggedCholeskyDecomposition(this.hessian, robust: true, inPlace: true);


                // Check if the decomposition exists
                if (this.decomposition.IsUndefined)
                {
                    // The Hessian is singular. Continue to the next
                    // iteration until the diagonal update transforms
                    // it back to non-singular.
                    continue;
                }


                // Solve using Cholesky decomposition
                this.deltas = this.decomposition.Solve(this.gradient);


                // Update weights using the calculated deltas
                for (int i = 0; i < this.Solution.Length; i++)
                    this.Solution[i] = this.weights[i] + this.deltas[i];


                // Calculate the new error
                sumOfSquaredErrors = this.ComputeError(inputs, outputs);

                // Update the objective function
                current = sumOfSquaredErrors;

                // If the object function is bigger than before, the method
                // is tried again using a greater damping factor.
            }

            // If this iteration caused a error drop, then next iteration
            //  will use a smaller damping factor.
            this.lambda /= this.v;

            return sumOfSquaredErrors;
        }



        private double computeErrors(double[][] input, double[] output, int[] block)
        {
            double sumOfSquaredErrors = 0.0;

            // for each input sample
            foreach (int i in block)
            {
                double actual = this.Function(this.Solution, input[i]);
                double expected = output[i];

                double e = expected - actual;
                sumOfSquaredErrors += e * e;

                this.errors[i] = e;
            }

            return sumOfSquaredErrors / 2.0;
        }

        private void computeJacobian(double[][] input, int[] block)
        {
            double[] derivatives = new double[this.NumberOfParameters];

            // for each input sample
            foreach (int i in block)
            {
                this.Gradient(this.Solution, input[i], derivatives);

                // copy the gradient vector into the Jacobian
                for (int j = 0; j < derivatives.Length; j++)
                    this.jacobian[j][i] = derivatives[j];
            }
        }

    }
}
