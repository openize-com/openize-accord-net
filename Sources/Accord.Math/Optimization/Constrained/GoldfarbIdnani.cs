// Accord Math Library
// The Accord.NET Framework
// http://accord-framework.net
//
// Copyright © 1995, 1996, 1997, 1998
// Berwin A. Turlach <bturlach@stats.adelaide.edu.au>
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
// This work is based on the original Fortran implementation by Berwin Turlach,
// also shared under the LGPL license.
//

namespace FileFormat.Accord.Math.Optimization.Constrained
{
    using System;
    using System.Collections.Generic;
    using Base;
    using Constraints;
    using FileFormat.Accord.Core.Exceptions;
    using FileFormat.Accord.Math;
    using FileFormat.Accord.Math.Vector;
    using Matrix;

    /// <summary>
    ///   Status codes for the <see cref="GoldfarbIdnani"/>
    ///   constrained quadratic programming solver.
    /// </summary>
    /// 
    public enum GoldfarbIdnaniStatus
    {
        /// <summary>
        ///   Convergence was attained.
        /// </summary>
        /// 
        Success,

        /// <summary>
        ///   The quadratic problem matrix is not positive definite.
        /// </summary>
        /// 
        NonPositiveDefinite,

        /// <summary>
        ///   The posed constraints cannot be fulfilled.
        /// </summary>
        /// 
        NoPossibleSolution
    }

    /// <summary>
    ///   Goldfarb-Idnani Quadratic Programming Solver.
    /// </summary>
    /// 
    /// <remarks>
    /// <para>
    ///   References:
    ///   <list type="bullet">
    ///     <item><description><a href="http://www.javaquant.net/papers/GoldfarbIdnani.pdf">
    ///       Goldfarb D., Idnani A. (1982) Dual and Primal-Dual Methods for Solving Strictly Convex Quadratic Programs.
    ///       Available on: http://www.javaquant.net/papers/GoldfarbIdnani.pdf .</a></description></item>
    ///     <item><description><a href="http://www.javaquant.net/papers/GoldfarbIdnani.pdf">
    ///       Berwin A Turlach. QuadProg, Quadratic Programming Solver (implementation in Fortran).
    ///       Available on:  http://school.maths.uwa.edu.au/~berwin/software/quadprog.html .</a></description></item>
    ///   </list>
    /// </para>   
    /// </remarks>
    /// 
    /// <example>
    /// <para>
    ///   There are three ways to state a quadratic programming problem in this framework.</para>
    ///   
    /// <list type="bullet">
    ///   <item><description>
    ///   The first is to state the problem in its canonical form, explicitly stating the
    ///   matrix Q and vector d specifying the quadratic function and the matrices A and
    ///   vector b specifying the problem constraints.</description></item>
    ///   <item><description>
    ///   The second is to state the problem with lambda expressions using symbolic variables.</description></item>
    ///   <item><description>
    ///   The third is to state the problem using text strings.</description></item>
    /// </list>
    ///   
    /// <para>  
    ///   In the following section we will provide examples for those ways. 
    /// </para>
    /// 
    /// <para>
    ///   This is an example stating the problem using lambdas:</para>
    ///   <code source="Unit Tests\Accord.Tests.Math\Optimization\GoldfarbIdnaniTest.cs" region="doc_lambdas" />
    /// 
    /// <para>
    ///   This is an example stating the problem using strings:</para>
    ///   <code source="Unit Tests\Accord.Tests.Math\Optimization\GoldfarbIdnaniTest.cs" region="doc_string" />
    ///   
    /// <para>
    ///   And finally, an example stating the problem using matrices:</para>
    ///   <code source="Unit Tests\Accord.Tests.Math\Optimization\GoldfarbIdnaniTest.cs" region="doc_matrix" />
    /// </example>
    /// 
    /// <seealso cref="AugmentedLagrangian"/>
    /// 
    public class GoldfarbIdnani : BaseGradientOptimizationMethod,
        IOptimizationMethod, IOptimizationMethod<GoldfarbIdnaniStatus>
    {
        private double[,] hessian;
        private double[] linearTerms;

        private double[,] constraintMatrix;
        private double[] constraintValues;
        private double[] constraintTolerances;

        //private double[] work;
        private int r;
        private double[] work;
        private double[] iwrv;
        private double[] iwzv;
        private double[] iwuv;
        private double[] iwrm;
        private double[] iwsv;
        private double[] iwnbv;


        /// <summary>
        ///   Gets the total number of constraints in the problem.
        /// </summary>
        /// 
        public int NumberOfConstraints { get; private set; }

        /// <summary>
        ///   Gets how many constraints are inequality constraints.
        /// </summary>
        /// 
        public int NumberOfEqualities { get; private set; }

        /// <summary>
        ///   Gets the total number of iterations performed on the
        ///   last call to the <see cref="Minimize"/> or <see cref="Maximize"/> methods.
        /// </summary>
        /// 
        public int Iterations { get; set; }

        /// <summary>
        ///   Gets or sets the maximum number of iterations that should be 
        ///   performed before the method terminates. If set to zero, the 
        ///   method will run to completion. Default is 0.
        /// </summary>
        /// 
        public int MaxIterations { get; set; }

        /// <summary>
        ///   Gets the total number of constraint removals performed
        ///   on the last call to the <see cref="Minimize"/> or <see cref="Maximize"/> methods.
        /// </summary>
        /// 
        public int Deletions { get; set; }

        /// <summary>
        ///   Gets the Lagrangian multipliers for the
        ///   last solution found.
        /// </summary>
        /// 
        public double[] Lagrangian { get; private set; }


        /// <summary>
        ///   Gets the indices of the active constraints
        ///   found during the last call of the 
        ///   <see cref="Minimize"/> or <see cref="Maximize"/>
        ///   methods.
        /// </summary>
        /// 
        public int[] ActiveConstraints { get; private set; }

        /// <summary>
        ///   Gets the constraint matrix <c>A</c> for the problem.
        /// </summary>
        /// 
        public double[,] ConstraintMatrix
        {
            get { return this.constraintMatrix; }
        }

        /// <summary>
        ///   Gets the constraint values <c>b</c> for the problem.
        /// </summary>
        /// 
        public double[] ConstraintValues
        {
            get { return this.constraintValues; }
        }

        /// <summary>
        ///   Gets the constraint tolerances <c>b</c> for the problem.
        /// </summary>
        /// 
        public double[] ConstraintTolerances
        {
            get { return this.constraintTolerances; }
        }

        /// <summary>
        ///   Gets the matrix of quadratic terms of
        ///   the quadratic optimization problem.
        /// </summary>
        /// 
        public double[,] QuadraticTerms { get { return this.hessian; } }

        /// <summary>
        ///   Gets the vector of linear terms of the
        ///   quadratic optimization problem.
        /// </summary>
        /// 
        public double[] LinearTerms { get { return this.linearTerms; } }

        /// <summary>
        ///   Get the exit code returned in the last call to the
        ///   <see cref="IOptimizationMethod{TInput, TOutput}.Maximize()"/> or 
        ///   <see cref="IOptimizationMethod{TInput, TOutput}.Minimize()"/> methods.
        /// </summary>
        /// 
        public GoldfarbIdnaniStatus Status { get; private set; }

        /// <summary>
        ///   Constructs a new <see cref="GoldfarbIdnani"/> class.
        /// </summary>
        /// 
        /// <param name="function">The objective function to be optimized.</param>
        /// <param name="constraints">The problem's constraints.</param>
        /// 
        public GoldfarbIdnani(QuadraticObjectiveFunction function, IEnumerable<LinearConstraint> constraints)
            : this(function, new LinearConstraintCollection(constraints))
        {

        }

        /// <summary>
        ///   Constructs a new <see cref="GoldfarbIdnani"/> class.
        /// </summary>
        /// 
        /// <param name="function">The objective function to be optimized.</param>
        /// <param name="constraints">The problem's constraints.</param>
        /// 
        public GoldfarbIdnani(QuadraticObjectiveFunction function, LinearConstraintCollection constraints)
            : base(function.NumberOfVariables, function.Function, function.Gradient)
        {
            int equalities;

            // Create the constraint matrix A from the specified constraint list
            double[,] A = constraints.CreateMatrix(function.NumberOfVariables,
                out this.constraintValues, out this.constraintTolerances, out equalities);

            System.Diagnostics.Debug.Assert(A.GetLength(1) == function.NumberOfVariables);

            this.initialize(function.NumberOfVariables,
                function.QuadraticTerms, function.LinearTerms,
                A, this.constraintValues, equalities);
        }


        /// <summary>
        ///   Constructs a new instance of the <see cref="GoldfarbIdnani"/> class.
        /// </summary>
        /// 
        /// <param name="function">The objective function to be optimized.</param>
        /// <param name="constraintMatrix">The constraints matrix <c>A</c>.</param>
        /// <param name="constraintValues">The constraints values <c>b</c>.</param>
        /// <param name="numberOfEqualities">The number of equalities in the constraints.</param>
        /// 
        public GoldfarbIdnani(QuadraticObjectiveFunction function, double[,] constraintMatrix,
            double[] constraintValues, int numberOfEqualities = 0)
            : base(function.NumberOfVariables, function.Function, function.Gradient)
        {
            if (function.NumberOfVariables != constraintMatrix.GetLength(1))
            {
                throw new ArgumentException("The number of columns in the constraint matrix A "
                    + "should equal the number of variables in the problem.", "constraintMatrix");
            }

            if (constraintValues.Length != constraintMatrix.GetLength(0))
                throw new DimensionMismatchException("constraintValues");

            if (numberOfEqualities < 0 || numberOfEqualities > constraintValues.Length)
                throw new ArgumentOutOfRangeException("numberOfEqualities");

            this.constraintTolerances = Vector.Create(constraintValues.Length, LinearConstraint.DefaultTolerance);

            this.initialize(function.NumberOfVariables, function.QuadraticTerms,
                function.LinearTerms, constraintMatrix, constraintValues, numberOfEqualities);
        }

        /// <summary>
        ///   Constructs a new instance of the <see cref="GoldfarbIdnani"/> class.
        /// </summary>
        /// 
        /// <param name="quadratic">The symmetric matrix of quadratic terms defining the objective function.</param>
        /// <param name="linear">The vector of linear terms defining the objective function.</param>
        /// <param name="constraintMatrix">The constraints matrix <c>A</c>.</param>
        /// <param name="constraintValues">The constraints values <c>b</c>.</param>
        /// <param name="numberOfEqualities">The number of equalities in the constraints.</param>
        /// 
        public GoldfarbIdnani(double[,] quadratic, double[] linear,
            double[,] constraintMatrix, double[] constraintValues, int numberOfEqualities = 0)
            : this(new QuadraticObjectiveFunction(quadratic, linear), constraintMatrix, constraintValues, numberOfEqualities)
        {
        }

        private void initialize(int numberOfVariables, double[,] hessian, double[] linearTerms,
            double[,] constraintMatrix, double[] b, int numberOfEqualities)
        {
            this.NumberOfVariables = numberOfVariables;
            this.linearTerms = linearTerms;
            this.hessian = hessian;

            this.constraintMatrix = constraintMatrix;
            this.constraintValues = b;

            this.NumberOfEqualities = numberOfEqualities;
            this.NumberOfConstraints = constraintMatrix.GetLength(0);
            this.r = Math.Min(this.NumberOfVariables, this.NumberOfConstraints);

            this.ActiveConstraints = new int[this.NumberOfConstraints];
            this.Lagrangian = new double[this.NumberOfConstraints];
            this.Solution = new double[this.NumberOfVariables];

            // initialize work vector
            this.work = new double[this.NumberOfVariables];
            this.iwzv = new double[this.NumberOfVariables];
            this.iwrv = new double[this.r];
            this.iwuv = new double[this.r + 1];
            this.iwrm = new double[this.r * (this.r + 1) / 2];
            this.iwsv = new double[this.NumberOfConstraints];
            this.iwnbv = new double[this.NumberOfConstraints];
        }


        /// <summary>
        ///   Finds the minimum value of a function. The solution vector
        ///   will be made available at the <see cref="IOptimizationMethod{TInput, TOutput}.Solution"/> property.
        /// </summary>
        /// 
        /// <returns>
        ///   Returns <c>true</c> if the method converged to a <see cref="IOptimizationMethod{TInput, TOutput}.Solution"/>.
        ///   In this case, the found value will also be available at the <see cref="IOptimizationMethod{TInput, TOutput}.Value"/>
        ///   property.
        /// </returns>
        /// 
        public override bool Minimize()
        {
            double[,] h = new double[this.NumberOfVariables, this.NumberOfVariables];
            double[] d = new double[this.NumberOfVariables];

            // Prepare a minimization problem
            for (int i = 0; i < this.NumberOfVariables; i++)
                for (int j = 0; j < this.NumberOfVariables; j++)
                    h[i, j] = this.hessian[i, j];

            for (int i = 0; i < this.linearTerms.Length; i++)
                d[i] = -this.linearTerms[i];

            this.Status = this.minimize(h, d);

            this.Value = this.Function(this.Solution);

            return (this.Status == GoldfarbIdnaniStatus.Success);
        }


        /// <summary>
        ///   Finds the maximum value of a function. The solution vector
        ///   will be made available at the <see cref="IOptimizationMethod{TInput, TOutput}.Solution"/> property.
        /// </summary>
        /// <returns>
        ///   Returns <c>true</c> if the method converged to a <see cref="IOptimizationMethod{TInput, TOutput}.Solution"/>.
        ///   In this case, the found value will also be available at the <see cref="IOptimizationMethod{TInput, TOutput}.Value"/>
        ///   property.
        /// </returns>
        /// 
        public override bool Maximize()
        {
            double[,] h = new double[this.NumberOfVariables, this.NumberOfVariables];
            double[] d = new double[this.NumberOfVariables];

            // Prepare a maximization problem
            for (int i = 0; i < this.NumberOfVariables; i++)
                for (int j = 0; j < this.NumberOfVariables; j++)
                    h[i, j] = -this.hessian[i, j];

            for (int i = 0; i < d.Length; i++)
                d[i] = this.linearTerms[i];

            this.Status = this.minimize(h, d);

            this.Value = this.Function(this.Solution);

            return (this.Status == GoldfarbIdnaniStatus.Success);
        }

        /// <summary>
        ///   Not available.
        /// </summary>
        /// 
        protected override bool Optimize()
        {
            throw new NotImplementedException();
        }

        private GoldfarbIdnaniStatus minimize(double[,] D, double[] d)
        {
            int numberOfActiveConstraints;
            int[] activeConstraints = new int[this.NumberOfConstraints];

            int ierr = 0;

            // Call qpgen2 routine from Turlach's Fortran code
            this.qpgen2(D, d, activeConstraints, out numberOfActiveConstraints, ref ierr);

            if (ierr != 0)
            {
                if (ierr == 1)
                    return GoldfarbIdnaniStatus.NoPossibleSolution;

                if (ierr == 2)
                    return GoldfarbIdnaniStatus.NonPositiveDefinite;

                throw new InvalidOperationException("Unexpected error.");
            }

            // Extract Lagrange multipliers from the work vector
            this.ActiveConstraints = Matrix.First(activeConstraints, numberOfActiveConstraints);

            for (int i = 0; i < this.ActiveConstraints.Length; i++)
                this.Lagrangian[this.ActiveConstraints[i]] = this.iwuv[i];

            return GoldfarbIdnaniStatus.Success;
        }



        //
        // This routine uses the Goldfarb/Idnani algorithm to solve the
        // following minimization problem:
        //
        //       minimize 1/2 * x^T D x + d^T x
        //       where   A1 x  = b1
        //               A2 x >= b2
        //
        // the matrix D is assumed to be positive definite.  Especially,
        // w.l.o.g. D is assumed to be symmetric. This is slightly different
        // from the original implementation by Berwin A. Turlach.
        // 
        // Input parameter:
        // dmat   nxn matrix, the matrix D from above (dp)
        //        *** WILL BE DESTROYED ON EXIT ***
        //        The user has two possibilities:
        //        a) Give D (ierr=0), in this case we use routines from LINPACK
        //           to decompose D.
        //        b) To get the algorithm started we need R^-1, where D=R^TR.
        //           So if it is cheaper to calculate R^-1 in another way (D may
        //           be a band matrix) then with the general routine, the user
        //           may pass R^{-1}.  Indicated by ierr not equal to zero.
        // dvec   nx1 vector, the vector d from above (dp)
        //        *** WILL BE DESTROYED ON EXIT ***
        //        contains on exit the solution to the initial, i.e.,
        //        unconstrained problem
        // fddmat scalar, the leading dimension of the matrix dmat
        // n      the dimension of dmat and dvec (int)
        // amat   nxq matrix, the matrix A from above (dp) [ A=(A1 A2) ]
        //        *** ENTRIES CORRESPONDING TO EQUALITY CONSTRAINTS MAY HAVE
        //            CHANGED SIGNES ON EXIT ***
        // bvec   qx1 vector, the vector of constants b in the constraints (dp)
        //        [ b = (b1^T b2^T)^T ]
        //        *** ENTRIES CORRESPONDING TO EQUALITY CONSTRAINTS MAY HAVE
        //            CHANGED SIGNES ON EXIT ***
        // fdamat the first dimension of amat as declared in the calling program. 
        //        fdamat >= n !!
        // q      int, the number of constraints.
        // meq    int, the number of equality constraints, 0 <= meq <= q.
        // ierr   int, code for the status of the matrix D:
        //           ierr =  0, we have to decompose D
        //           ierr != 0, D is already decomposed into D=R^TR and we were
        //                      given R^{-1}.
        //
        // Output parameter:
        // sol   nx1 the final solution (x in the notation above)
        // crval scalar, the value of the criterion at the minimum      
        // iact  qx1 vector, the constraints which are active in the final
        //       fit (int)
        // nact  scalar, the number of constraints active in the final fit (int)
        // iter  2x1 vector, first component gives the number of "main" 
        //       iterations, the second one says how many constraints were
        //        deleted after they became active
        //  ierr  int, error code on exit, if
        //          ierr = 0, no problems
        //          ierr = 1, the minimization problem has no solution
        //          ierr = 2, problems with decomposing D, in this case sol
        //                     contains garbage!!
        //
        //  Working space:
        //  work  vector with length at least 2n+r*(r+5)/2 + 2q +1
        //        where r=min(n,q)
        //
        private void qpgen2(double[,] dmat, double[] dvec, int[] iact, out int nact, ref int ierr)
        {

            int n = this.NumberOfVariables;
            int q = this.NumberOfConstraints;
            int meq = this.NumberOfEqualities;
            double[,] amat = (double[,])this.constraintMatrix.Clone();
            double[] bvec = (double[])this.constraintValues.Clone();
            double[] sol = this.Solution;

            int l1;
            double gc, gs, tt, sum;

            double f = 0;
            nact = 0;


            // Store the initial dvec to calculate below the
            //  unconstrained minima of the critical value.

            Array.Clear(this.iwzv, 0, this.iwzv.Length);
            Array.Clear(this.iwrv, 0, this.iwrv.Length);
            Array.Clear(this.iwuv, 0, this.iwuv.Length);
            Array.Clear(this.iwrm, 0, this.iwrm.Length);
            Array.Clear(this.iwsv, 0, this.iwsv.Length);
            Array.Clear(this.iwnbv, 0, this.iwnbv.Length);

            for (int i = 0; i < dvec.Length; i++)
                this.work[i] = dvec[i];

            for (int i = 0; i < iact.Length; i++)
                iact[i] = -1;


            // Get the initial solution
            if (ierr == 0)
            {
                // L'L = chol(D)
                bool success = this.dpofa(dmat);

                if (!success)
                {
                    ierr = 2;
                    return;
                }

                // L*x = d
                this.dposl(dmat, dvec);

                // D = inv(L)
                this.dpori(dmat);
            }
            else
            {
                // Matrix D is already factorized, so we have to multiply d first with 
                // R^-T and then with R^-1.  R^-1 is stored in the upper half of the 
                // array dmat. 

                for (int j = 0; j < sol.Length; j++)
                {
                    sol[j] = 0.0;

                    for (int i = 0; i < j; i++)
                        sol[j] += dmat[j, i] * dvec[i];
                }

                for (int j = 0; j < dvec.Length; j++)
                {
                    dvec[j] = 0.0;

                    for (int i = j; i < sol.Length; i++)
                        dvec[j] += dmat[i, j] * sol[i];
                }
            }

            // Set upper triangular of dmat to zero, store dvec in sol and 
            //   calculate value of the criterion at unconstrained minima

            f = 0.0;

            // calculate some constants, i.e., from which index on 
            // the different quantities are stored in the work matrix 

            for (int j = 0; j < sol.Length; j++)
            {
                sol[j] = dvec[j];
                f += this.work[j] * sol[j];
                this.work[j] = 0.0;

                for (int i = j + 1; i < n; i++)
                    dmat[j, i] = 0.0;
            }

            f = -f / 2.0;
            ierr = 0;


            // calculate the norm of each column of the A matrix 

            for (int i = 0; i < this.iwnbv.Length; i++)
            {
                sum = 0.0;
                for (int j = 0; j < n; j++)
                    sum += amat[i, j] * amat[i, j];

                this.iwnbv[i] = Math.Sqrt(sum);
            }

            nact = 0;
            this.Iterations = 0;
            this.Deletions = 0;


        L50: // start a new iteration 

            if (this.Token.IsCancellationRequested)
                return;

            this.Iterations++;

            if (this.MaxIterations > 0 && this.Iterations > this.MaxIterations)
                return;

            // calculate all constraints and check which are still violated 
            // for the equality constraints we have to check whether the normal 
            // vector has to be negated (as well as bvec in that case) 

            int l = 0;

            for (int i = 0; i < bvec.Length; i++)
            {
                sum = -bvec[i];

                for (int j = 0; j < sol.Length; j++)
                    sum += amat[i, j] * sol[j];

                if (Math.Abs(sum) < Constants.DoubleEpsilon)
                    sum = 0.0;

                if (Double.IsNaN(sum))
                    sum = 0.0;

                if (i >= meq)
                {
                    // this is an inequality constraint
                    this.iwsv[l] = sum;
                }
                else
                {
                    // this is an equality constraint
                    this.iwsv[l] = -Math.Abs(sum);

                    if (sum > 0.0)
                    {
                        for (int j = 0; j < n; j++)
                            amat[i, j] = -amat[i, j];
                        bvec[i] = -bvec[i];
                    }
                }

                l++;
            }

            // as safeguard against rounding errors set 
            // already active constraints explicitly to zero 

            for (int i = 0; i < nact; i++)
                this.iwsv[iact[i]] = 0.0;

            // We weight each violation by the number of non-zero elements in the 
            // corresponding row of A. then we choose the violated constraint which 
            // has maximal absolute value, i.e., the minimum. By obvious commenting
            // and uncommenting we can choose the strategy to take always the first
            // constraint which is violated. ;-) 

            int nvl = -1;
            double temp = 0.0;
            for (int i = 0; i < this.iwnbv.Length; i++)
            {
                double w = temp * this.iwnbv[i];
                double tol = this.constraintTolerances[i];

                if (this.iwsv[i] < w - tol)
                {
                    nvl = i;
                    temp = this.iwsv[i] / this.iwnbv[i];
                }

                // if (work(iwsv+i) .LT. 0.d0) then 
                //     nvl = i 
                //     goto 72 
                // endif 
            }

            if (nvl == -1)
                return;


        L55:

            // calculate d=J^Tn^+ where n^+ is the normal vector of the violated 
            // constraint. J is stored in dmat in this implementation!! 
            // if we drop a constraint, we have to jump back here. 

            for (int i = 0; i < this.work.Length; i++)
            {
                sum = 0.0;
                for (int j = 0; j < n; j++)
                    sum += dmat[i, j] * amat[nvl, j];

                this.work[i] = sum;
            }

            // Now calculate z = J_2 d_2 ...

            for (int i = 0; i < this.iwzv.Length; i++)
                this.iwzv[i] = 0.0;

            for (int j = nact; j < this.work.Length; j++)
                for (int i = 0; i < this.iwzv.Length; i++)
                    this.iwzv[i] += dmat[j, i] * this.work[j];

            // ... and r = R^{-1} d_1, check also if r has positive elements
            // (among the entries corresponding to inequalities constraints). 

            l1 = 0;
            int it1 = 0;
            double t1 = 0;
            bool t1inf = true;

            for (int i = nact - 1; i >= 0; i--)
            {
                sum = this.work[i];
                l = ((i + 1) * (i + 4)) / 2 - 1;
                l1 = l - i - 1;

                for (int j = i + 1; j < nact; j++)
                {
                    sum -= this.iwrm[l] * this.iwrv[j];
                    l += j + 1;
                }

                sum /= this.iwrm[l1];

                this.iwrv[i] = sum;

                if (iact[i] < meq)
                    continue;

                if (sum <= 0.0)
                    continue;

                if (Double.IsNaN(sum))
                    continue;

                t1inf = false;
                it1 = i + 1;
            }


            // if r has positive elements, find the partial step length t1, which is 
            // the maximum step in dual space without violating dual feasibility. 
            // it1 stores in which component t1, the min of u/r, occurs. 

            if (!t1inf)
            {
                t1 = this.iwuv[it1 - 1] / this.iwrv[it1 - 1];

                for (int i = 0; i < nact; i++)
                {
                    if (iact[i] < meq)
                        continue;

                    if (this.iwrv[i] <= 0.0)
                        continue;

                    temp = this.iwuv[i] / this.iwrv[i];

                    if (temp < t1)
                    {
                        t1 = temp;
                        it1 = i + 1;
                    }
                }
            }


            // test if the z vector is equal to zero 

            sum = 0.0;
            for (int i = 0; i < this.iwzv.Length; i++)
                sum += this.iwzv[i] * this.iwzv[i];

            if (Math.Abs(sum) <= Constants.DoubleEpsilon)
            {
                // No step in primal space such that the new constraint becomes 
                // feasible. Take step in dual space and drop a constant. 

                if (t1inf)
                {
                    // No step in dual space possible 
                    // either, problem is not solvable
                    ierr = 1;
                    return;
                }
                else
                {
                    // we take a partial step in dual space and drop constraint it1, 
                    // that is, we drop the it1-th active constraint. 
                    // then we continue at step 2(a) (marked by label 55) 

                    for (int i = 0; i < nact; i++)
                        this.iwuv[i] -= t1 * this.iwrv[i];

                    this.iwuv[nact] += t1;
                    goto L700;
                }
            }
            else
            {
                // compute full step length t2, minimum step in primal space such that 
                // the constraint becomes feasible. 
                // keep sum (which is z^Tn^+) to update crval below! 

                sum = 0.0;
                for (int i = 0; i < this.iwzv.Length; i++)
                    sum += this.iwzv[i] * amat[nvl, i];

                tt = -this.iwsv[nvl] / sum;
                bool t2min = true;

                if (!t1inf)
                {
                    if (t1 < tt)
                    {
                        tt = t1;
                        t2min = false;
                    }
                }

                // take step in primal and dual space 
                for (int i = 0; i < sol.Length; i++)
                    sol[i] += tt * this.iwzv[i];

                f += tt * sum * (tt / 2.0 + this.iwuv[nact]);

                for (int i = 0; i < nact; i++)
                    this.iwuv[i] -= tt * this.iwrv[i];

                this.iwuv[nact] += tt;

                // if it was a full step, then we check whether further constraints are 
                // violated otherwise we can drop the current constraint and iterate once 
                // more 

                if (t2min)
                {
                    // we took a full step. Thus add constraint nvl to the list of active 
                    // constraints and update J and R 

                    iact[nact++] = nvl;


                    // to update R we have to put the first nact-1 components of the d vector 
                    // into column (nact) of R 

                    l = ((nact - 1) * (nact)) / 2;
                    for (int i = 0; i < nact - 1; i++, l++)
                        this.iwrm[l] = this.work[i];

                    // if now nact=n, then we just have to add the last element to the new 
                    // row of R. 

                    // Otherwise we use Givens transformations to turn the vector d(nact:n) 
                    // into a multiple of the first unit vector. That multiple goes into the 
                    // last element of the new row of R and J is accordingly updated by the 
                    // Givens transformations. 

                    if (nact == n)
                    {
                        this.iwrm[l] = this.work[n - 1];
                    }
                    else
                    {
                        for (int i = n - 1; i >= nact; i--)
                        {
                            // We have to find the Givens rotation which will reduce the element 
                            // (l1) of d to zero. If it is already zero we don't have to do anything,
                            // except of decreasing l1 

                            if (this.work[i] == 0.0)
                                continue;

                            gc = Math.Max(Math.Abs(this.work[i - 1]), Math.Abs(this.work[i]));
                            gs = Math.Min(Math.Abs(this.work[i - 1]), Math.Abs(this.work[i]));
                            temp = Special.Sign(gc * Math.Sqrt(1.0 + (gs * gs) / (gc * gc)), this.work[i - 1]);
                            gc = this.work[i - 1] / temp;
                            gs = this.work[i] / temp;

                            // The Givens rotation is done with the matrix (gc gs, gs -gc). If
                            // gc is one, then element (i) of d is zero compared with element 
                            // (l1-1). Hence we don't have to do anything. If gc is zero, then
                            // we just have to switch column (i) and column (i-1) of J. Since 
                            // we only switch columns in J, we have to be careful how we update
                            // d depending on the sign of gs. Otherwise we have to apply the
                            // Givens rotation to these columns. The i-1 element of d has to be
                            // updated to temp. 

                            if (gc == 1.0)
                                continue;

                            if (gc == 0.0)
                            {
                                this.work[i - 1] = gs * temp;

                                for (int j = 0; j < n; j++)
                                {
                                    temp = dmat[i - 1, j];
                                    dmat[i - 1, j] = dmat[i, j];
                                    dmat[i, j] = temp;
                                }
                            }
                            else
                            {
                                this.work[i - 1] = temp;
                                double nu = gs / (gc + 1.0);

                                for (int j = 0; j < n; j++)
                                {
                                    temp = gc * dmat[i - 1, j] + gs * dmat[i, j];
                                    dmat[i, j] = nu * (dmat[i - 1, j] + temp) - dmat[i, j];
                                    dmat[i - 1, j] = temp;
                                }
                            }
                        }

                        // l is still pointing to element (nact,nact) of
                        // the matrix R. So store d(nact) in R(nact,nact) 
                        this.iwrm[l] = this.work[nact - 1];
                    }
                }
                else
                {

                    // We took a partial step in dual space. Thus drop constraint it1, 
                    // that is, we drop the it1-th active constraint. Then we continue
                    // at step 2(a) (marked by label 55) but since the fit changed, we
                    // have to recalculate now "how much" the fit violates the chosen
                    // constraint now. 

                    sum = -bvec[nvl];

                    for (int j = 0; j < sol.Length; j++)
                        sum += sol[j] * amat[nvl, j];

                    if (nvl + 1 > meq)
                    {
                        this.iwsv[nvl] = sum;
                    }
                    else
                    {
                        this.iwsv[nvl] = -Math.Abs(sum);

                        if (sum > 0.0)
                        {
                            for (int j = 0; j < n; j++)
                                amat[nvl, j] = -amat[nvl, j];

                            bvec[nvl] = -bvec[nvl];
                        }
                    }

                    goto L700;
                }
            }

            goto L50;


        L700: // Drop constraint it1 

            // if it1 = nact it is only necessary
            // to update the vector u and nact

            if (it1 == nact)
                goto L799;


        L797: // After updating one row of R (column of J) we will also come back here

            // We have to find the Givens rotation which will reduce the element
            // (it1+1,it1+1) of R to zero. If it is already zero we don't have to
            // do anything except of updating u, iact, and shifting column (it1+1)
            // of R to column (it1). Then l  will point to element (1,it1+1) of R 
            // and l1 will point to element (it1+1,it1+1) of R.

            l = it1 * (it1 + 1) / 2;
            l1 = l + it1;

            if (this.iwrm[l1] == 0.0)
                goto L798;

            gc = Math.Max(Math.Abs(this.iwrm[l1 - 1]), Math.Abs(this.iwrm[l1]));
            gs = Math.Min(Math.Abs(this.iwrm[l1 - 1]), Math.Abs(this.iwrm[l1]));
            temp = Special.Sign(gc * Math.Sqrt(1.0 + (gs * gs) / (gc * gc)), this.iwrm[l1 - 1]);
            gc = this.iwrm[l1 - 1] / temp;
            gs = this.iwrm[l1] / temp;


            // The Givens rotation is done with the matrix (gc gs, gs -gc). If gc is
            // one, then element (it1+1,it1+1) of R is zero compared with element 
            // (it1,it1+1). Hence we don't have to do anything. if gc is zero, then
            // we just have to switch row (it1) and row (it1+1) of R and column (it1)
            // and column (it1+1) of J. Since we switch rows in R and columns in J,
            // we can ignore the sign of gs. Otherwise we have to apply the Givens 
            // rotation to these rows/columns. 

            if (gc == 1.0)
                goto L798;

            if (gc == 0.0)
            {
                for (int i = it1; i < nact; i++)
                {
                    temp = this.iwrm[l1 - 1];
                    this.iwrm[l1 - 1] = this.iwrm[l1];
                    this.iwrm[l1] = temp;
                    l1 += i + 1;
                }

                for (int i = 0; i < n; i++)
                {
                    temp = dmat[it1 - 1, i];
                    dmat[it1 - 1, i] = dmat[it1, i];
                    dmat[it1, i] = temp;
                }
            }
            else
            {
                double nu = gs / (gc + 1.0);

                for (int i = it1; i < nact; i++)
                {
                    temp = gc * this.iwrm[l1 - 1] + gs * this.iwrm[l1];
                    this.iwrm[l1] = nu * (this.iwrm[l1 - 1] + temp) - this.iwrm[l1];
                    this.iwrm[l1 - 1] = temp;
                    l1 += i + 1;
                }

                for (int i = 0; i < n; i++)
                {
                    temp = gc * dmat[it1 - 1, i] + gs * dmat[it1, i];
                    dmat[it1, i] = nu * (dmat[it1 - 1, i] + temp) - dmat[it1, i];
                    dmat[it1 - 1, i] = temp;
                }
            }

        L798:

            // shift column (it1+1) of R to column (it1) (that is, the first it1 
            // elements). The position of element (1,it1+1) of R was calculated 
            // above and stored in l. 

            l1 = l - it1;
            for (int i = 0; i < it1; i++, l++, l1++)
            {
                this.iwrm[l1] = this.iwrm[l];
            }

            // update vector u and iact as necessary 
            // Continue with updating the matrices J and R 

            this.iwuv[it1 - 1] = this.iwuv[it1];
            iact[it1 - 1] = iact[it1];
            it1++;

            if (it1 < nact)
                goto L797;

        L799:

            this.iwuv[nact - 1] = this.iwuv[nact];
            this.iwuv[nact] = 0.0;
            iact[nact - 1] = -1;

            nact--;
            this.Deletions++;

            goto L55;
        }

        private void dpori(double[,] a)
        {
            int n = this.NumberOfVariables;

            for (int k = 0; k < n; k++)
            {
                a[k, k] = 1.0 / a[k, k];
                double t = -a[k, k];

                for (int j = 0; j < k; j++)
                    a[k, j] = a[k, j] * t;

                for (int j = k + 1; j < n; j++)
                {
                    t = a[j, k];
                    a[j, k] = 0.0;

                    for (int i = 0; i <= k; i++)
                        a[j, i] = t * a[k, i] + a[j, i];
                }
            }
        }

        private void dposl(double[,] a, double[] b)
        {
            int n = this.NumberOfVariables;

            // solve trans(r)*y = b
            for (int k = 0; k < b.Length; k++)
            {
                for (int i = 0; i < k; i++)
                    b[k] -= b[i] * a[k, i];
                b[k] = b[k] / a[k, k];
            }

            // solve r*x = y
            for (int k = n - 1; k >= 0; k--)
            {
                for (int i = k + 1; i < n; i++)
                    b[k] -= b[i] * a[i, k];
                b[k] /= a[k, k];
            }
        }

        private bool dpofa(double[,] a)
        {
            int n = this.NumberOfVariables;

            for (int j = 0; j < n; j++)
            {
                double s = 0;
                for (int k = 0; k < j; k++)
                {
                    double t = a[j, k];
                    for (int i = 0; i < k; i++)
                        t -= a[j, i] * a[k, i];
                    t = t / a[k, k];

                    a[j, k] = t;
                    s += t * t;
                }

                s = a[j, j] - s;

                // Use a tolerance for positive-definiteness
                if (s <= 1e-14 * Math.Abs(a[j, j]))
                    return false;

                a[j, j] = Math.Sqrt(s);
            }

            return true;
        }


    }
}
