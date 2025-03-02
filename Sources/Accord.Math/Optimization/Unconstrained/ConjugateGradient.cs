﻿// Accord Math Library
// The Accord.NET Framework
// http://accord-framework.net
//
// Copyright © César Souza, 2009-2017
// cesarsouza at gmail.com
//
// Copyright © Jorge Nocedal, 1990
// http://users.eecs.northwestern.edu/~nocedal/
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

namespace Openize.Accord.Math.Optimization.Unconstrained
{
    using System;
    using Openize.Accord.Math;
    using Openize.Accord.Math.Matrix;
    using Openize.Accord.Core.MachineLearning;
    using Openize.Accord.Math.Optimization.Base;

    /// <summary>
    ///   Conjugate gradient direction update formula.
    /// </summary>
    /// 
    public enum ConjugateGradientMethod
    {
        /// <summary>
        ///   Fletcher-Reeves formula.
        /// </summary>
        /// 
        FletcherReeves = 1,

        /// <summary>
        ///   Polak-Ribière formula.
        /// </summary>
        /// 
        /// <remarks>
        ///   The Polak-Ribière is known to perform better for non-quadratic functions.
        /// </remarks>
        /// 
        PolakRibiere = 2,

        /// <summary>
        ///   Polak-Ribière formula.
        /// </summary>
        /// 
        /// <remarks>
        ///   The Polak-Ribière is known to perform better for non-quadratic functions.
        ///   The positive version B=max(0,Bpr) provides a direction reset automatically.
        /// </remarks>
        /// 
        PositivePolakRibiere = 3,
    }

    /// <summary>
    ///   Conjugate Gradient exit codes.
    /// </summary>
    /// 
    public enum ConjugateGradientCode
    {
        /// <summary>
        ///   Success.
        /// </summary>
        /// 
        Success,

        /// <summary>
        ///   Invalid step size.
        /// </summary>
        /// 
        StepSize = 1,

        /// <summary>
        ///   Descent direction was not obtained.
        /// </summary>
        /// 
        DescentNotObtained = -2,

        /// <summary>
        ///   Rounding errors prevent further progress. There may not be a step 
        ///   which satisfies the sufficient decrease and curvature conditions. 
        ///   Tolerances may be too small.
        /// </summary>
        /// 
        RoundingErrors = 6,

        /// <summary>
        ///   The step size has reached the upper bound.
        /// </summary>
        /// 
        StepHigh = 5,

        /// <summary>
        ///   The step size has reached the lower bound.
        /// </summary>
        /// 
        StepLow = 4,

        /// <summary>
        ///   Maximum number of function evaluations has been reached.
        /// </summary>
        /// 
        MaximumEvaluations = 3,

        /// <summary>
        ///   Relative width of the interval of uncertainty is at machine precision.
        /// </summary>
        /// 
        Precision = 2,
    }

    /// <summary>
    ///   Conjugate Gradient (CG) optimization method.
    /// </summary>
    /// 
    /// <remarks>
    /// <para>
    ///   In mathematics, the conjugate gradient method is an algorithm for the numerical solution of
    ///   particular systems of linear equations, namely those whose matrix is symmetric and positive-
    ///   definite. The conjugate gradient method is an iterative method, so it can be applied to sparse
    ///   systems that are too large to be handled by direct methods. Such systems often arise when
    ///   numerically solving partial differential equations. The nonlinear conjugate gradient method 
    ///   generalizes the conjugate gradient method to nonlinear optimization (Wikipedia, 2011).</para>
    /// <para>
    ///   T</para>
    /// <para>
    ///   The framework implementation of this method is based on the original FORTRAN source code
    ///   by Jorge Nocedal (see references below). The original FORTRAN source code of CG+ (for large
    ///   scale unconstrained problems) is available at http://users.eecs.northwestern.edu/~nocedal/CG+.html
    ///   and had been made freely available for educational or commercial use. The original authors
    ///   expect that all publications describing work using this software quote the (Gilbert and Nocedal, 1992)
    ///   reference given below.</para>
    /// 
    /// <para>
    ///   References:
    ///   <list type="bullet">
    ///     <item><description><a href="http://users.eecs.northwestern.edu/~nocedal/CG+.html">
    ///        J. C. Gilbert and J. Nocedal. Global Convergence Properties of Conjugate Gradient
    ///        Methods for Optimization, (1992) SIAM J. on Optimization, 2, 1.</a></description></item>
    ///     <item><description>
    ///        Wikipedia contributors, "Nonlinear conjugate gradient method," Wikipedia, The Free 
    ///        Encyclopedia, http://en.wikipedia.org/w/index.php?title=Nonlinear_conjugate_gradient_method
    ///        (accessed December 22, 2011).</description></item>
    ///     <item><description>
    ///        Wikipedia contributors, "Conjugate gradient method," Wikipedia, The Free Encyclopedia,
    ///        http://en.wikipedia.org/w/index.php?title=Conjugate_gradient_method 
    ///        (accessed December 22, 2011).</description></item>
    ///    </list></para>
    /// </remarks>
    /// 
    /// <seealso cref="BroydenFletcherGoldfarbShanno"/>
    /// <seealso cref="ResilientBackpropagation"/>
    /// <seealso cref="BoundedBroydenFletcherGoldfarbShanno"/>
    /// <seealso cref="TrustRegionNewtonMethod"/>
    /// 
    public class ConjugateGradient : BaseGradientOptimizationMethod,
        IGradientOptimizationMethod, IOptimizationMethod<ConjugateGradientCode>, 
        ISupportsCancellation
    {

        private double[] g; // gradient at current solution

        private double[] d;
        private double[] gold;
        private double[] w;

        private double ftol = 1e-4;
        private double gtol = 0.1;
        private double xtol = 1e-17;
        private double stpmin = 1e-20;
        private double stpmax = 1e20;
        private double epsilon = 1e-5;
        private int maxfev = 40;
        private int iterations;
        private int evaluations;
        private int searches;
        private int maxIterations;
        private double tolerance = 1e-5;


        /// <summary>
        ///   Gets or sets the relative difference threshold
        ///   to be used as stopping criteria between two
        ///   iterations. Default is 0 (iterate until convergence). 
        /// </summary>
        /// 
        public double Tolerance
        {
            get { return this.tolerance; }
            set { this.tolerance = value; }
        }

        /// <summary>
        ///   Gets or sets the maximum number of iterations
        ///   to be performed during optimization. Default
        ///   is 0 (iterate until convergence).
        /// </summary>
        /// 
        public int MaxIterations
        {
            get { return this.maxIterations; }
            set { this.maxIterations = value; }
        }

        /// <summary>
        ///   Gets or sets the conjugate gradient update 
        ///   method to be used during optimization.
        /// </summary>
        /// 
        public ConjugateGradientMethod Method { get; set; }

        /// <summary>
        ///   Gets the number of iterations performed 
        ///   in the last call to <see cref="IOptimizationMethod{TInput, TOutput}.Minimize()"/>.
        /// </summary>
        /// 
        /// <value>
        ///   The number of iterations performed
        ///   in the previous optimization.</value>
        ///   
        public int Iterations
        {
            get { return this.iterations; }
        }

        /// <summary>
        ///   Gets the number of function evaluations performed
        ///   in the last call to <see cref="IOptimizationMethod{TInput, TOutput}.Minimize()"/>.
        /// </summary>
        /// 
        /// <value>
        ///   The number of evaluations performed
        ///   in the previous optimization.</value>
        ///   
        public int Evaluations
        {
            get { return this.evaluations; }
        }

        /// <summary>
        ///   Gets the number of linear searches performed
        ///   in the last call to <see cref="IOptimizationMethod{TInput, TOutput}.Minimize()"/>.
        /// </summary>
        /// 
        public int Searches
        {
            get { return this.searches; }
        }

        /// <summary>
        ///   Get the exit code returned in the last call to the
        ///   <see cref="IOptimizationMethod{TInput, TOutput}.Maximize()"/> or 
        ///   <see cref="IOptimizationMethod{TInput, TOutput}.Minimize()"/> methods.
        /// </summary>
        /// 
        public ConjugateGradientCode Status { get; private set; }

        /// <summary>
        ///   Occurs when progress is made during the optimization.
        /// </summary>
        /// 
        public event EventHandler<OptimizationProgressEventArgs> Progress;

        /// <summary>
        ///   Creates a new instance of the CG optimization algorithm.
        /// </summary>
        /// 
        public ConjugateGradient()
            : base()
        {
        }

        /// <summary>
        ///   Creates a new instance of the CG optimization algorithm.
        /// </summary>
        /// 
        /// <param name="numberOfVariables">The number of free parameters in the optimization problem.</param>
        /// 
        public ConjugateGradient(int numberOfVariables)
            : base(numberOfVariables)
        {
        }

        /// <summary>
        ///   Creates a new instance of the CG optimization algorithm.
        /// </summary>
        /// 
        /// <param name="numberOfVariables">The number of free parameters in the function to be optimized.</param>
        /// <param name="function">The function to be optimized.</param>
        /// <param name="gradient">The gradient of the function.</param>
        /// 
        public ConjugateGradient(int numberOfVariables,
            Func<double[], double> function, Func<double[], double[]> gradient)
            : base(numberOfVariables, function, gradient)
        {
        }

        /// <summary>
        /// Called when the <see cref="IOptimizationMethod{TInput, TOutput}.NumberOfVariables" /> property has changed.
        /// </summary>
        /// 
        /// <param name="numberOfVariables">The number of variables.</param>
        /// 
        protected override void OnNumberOfVariablesChanged(int numberOfVariables)
        {
            base.OnNumberOfVariablesChanged(numberOfVariables);

            this.d = new double[numberOfVariables];
            this.gold = new double[numberOfVariables];
            this.w = new double[numberOfVariables];
        }

        /// <summary>
        ///   Implements the actual optimization algorithm. This
        ///   method should try to minimize the objective function.
        /// </summary>
        /// 
        protected override bool Optimize()
        {
            // This code has been adapted from the original
            // FORTRAN function CGFAM by Jorge Nocedal, 1992.

            int irest = 1;
            int n = this.NumberOfVariables;
            double[] x = this.Solution;

            double f = this.Function(x);
            this.g = this.Gradient(x);

            int method = (int)this.Method;

            this.iterations = 0;
            this.evaluations = 1;
            bool bnew = true;
            int nrst = 0;
            int im = 0;   // Number of times betapr was negative for method 2 or 3
            this.searches = 0; // Number of line search iterations after Wolfe conditions were satisfied.
            double dg0 = 0;

            for (int i = 0; i < this.g.Length; ++i)
                this.d[i] = -this.g[i];

            double gnorm = Norm.Euclidean(this.g);
            double xnorm = Math.Max(1.0, Norm.Euclidean(x));
            double stp1 = 1.0 / gnorm;
            double f_old = f;


            bool finish = false;

            // Make initial progress report with initialization parameters
            if (this.Progress != null) this.Progress(this, new OptimizationProgressEventArgs
                (this.iterations, this.evaluations, this.g, gnorm, this.Solution, xnorm, f, stp1, finish));


            // Main iteration
            while (!finish)
            {
                if (this.Token.IsCancellationRequested)
                    break;

                this.iterations++;
                nrst++;

                // Call the line search routine of Mor'e and Thuente
                // (modified for Nocedal's CG method) 
                // ------------------------------------------------- 
                //
                //  J.J. Mor'e and D. Thuente, "Linesearch Algorithms with Guaranteed 
                //  Sufficient Decrease". ACM Transactions on Mathematical 
                //  Software 20 (1994), pp 286-307. 
                //

                int nfev = 0;
                int info = 0;
                double dgout = 0;

                // Save original gradient
                for (int i = 0; i < this.g.Length; i++)
                    this.gold[i] = this.g[i];

                double dg = this.d.Dot(this.g);
                double dgold = dg;
                double stp = 1.0;

                // Shanno-Phua's formula for trial step
                if (!bnew) stp = dg0 / dg;

                if (this.iterations == 1)
                    stp = stp1;

                int ides = 0;
                bnew = false;

            L72:

                // Call to the line search subroutine
                this.Status = this.cvsmod(ref f, this.d, ref stp, ref info, ref nfev, this.w, ref dg, ref dgout);

                if (this.Status != ConjugateGradientCode.Success)
                    return false;

                // Test if descent direction is obtained for methods 2 and 3
                double gg = Matrix.Dot(this.g, this.g);
                double gg0 = Matrix.Dot(this.g, this.gold);
                double betapr = (gg - gg0) / (gnorm * gnorm);

                // When nrst > n and irest == 1 then restart.
                if (irest == 1 && nrst > n)
                {
                    nrst = 0;
                    bnew = true;
                }
                else
                {
                    if (method != 1)
                    {
                        double dg1 = -gg + betapr * dgout;

                        if (dg1 >= 0.0)
                        {
                            ides++;

                            if (ides > 5)
                            {
                                this.Status = ConjugateGradientCode.DescentNotObtained;
                                return false;
                            }

                            goto L72; // retry
                        }
                    }
                }

                this.evaluations += nfev;
                this.searches += ides;

                // Determine correct beta value for method chosen
                double betafr = gg / (gnorm * gnorm);
                double beta = 0;

                if (nrst == 0)
                {
                    beta = 0.0;
                }
                else
                {
                    if (method == 1)
                    {
                        beta = betafr;
                    }
                    else if (method == 2)
                    {
                        beta = betapr;
                    }
                    else if (method == 3)
                    {
                        beta = Math.Max(0.0, betapr);
                    }

                    if ((method == 2 || method == 3) && betapr < 0.0)
                    {
                        im++;
                    }
                }

                // Compute the new direction
                for (int i = 0; i < this.g.Length; i++)
                    this.d[i] = -this.g[i] + beta * this.d[i];

                dg0 = dgold * stp;

                // Check for termination
                gnorm = Norm.Euclidean(this.g);
                xnorm = Math.Max(1.0, Norm.Euclidean(x));

                // Convergence test
                if (gnorm / xnorm <= this.epsilon)
                    finish = true;

                // Stopping criteria by function delta
                if (this.tolerance > 0 && this.iterations > 1)
                {
                    double delta = (f_old - f) / f;
                    f_old = f;

                    if (delta < this.tolerance)
                        finish = true;
                }

                // Stopping criteria by max iterations
                if (this.maxIterations > 0)
                {
                    if (this.iterations > this.maxIterations)
                        finish = true;
                }

                if (this.Progress != null)
                {
                    this.Progress(this, new OptimizationProgressEventArgs(this.iterations,
                        this.evaluations, this.g, gnorm, this.Solution, xnorm, f, stp, finish));
                }
            }

            return this.Status == ConjugateGradientCode.Success;
        }



        // TODO: Move to separate classes

        bool brackt;
        bool stage1;
        double finit;
        double dgtest;
        double width;
        double width1;

        double stmin;
        double stmax;
        double dg2;
        double dgx;
        double dgy;
        int infoc;
        double stx;
        double fx;
        double sty;
        double fy;
        double ftest1;

        private unsafe ConjugateGradientCode cvsmod(ref double f, double[] s, ref double stp, ref int info,
             ref int nfev, double[] wa, ref double dginit, ref double dgout)
        {
            int n = this.NumberOfVariables;

            double[] x = this.Solution;

            if (info == 1)
                goto L321;

            this.infoc = 1;

            if (stp <= 0) // Check the input parameters for errors
                return ConjugateGradientCode.StepSize;

            // Compute the initial gradient in the search direction
            // and check that S is a descent direction.

            if (dginit >= 0)
                throw new LineSearchFailedException(0, "The search direction is not a descent direction.");

            // Initialize local variables
            this.brackt = false;
            this.stage1 = true;
            nfev = 0;
            this.finit = f;
            this.dgtest = this.ftol * dginit;
            this.width = this.stpmax - this.stpmin;
            this.width1 = this.width / 0.5;

            for (int j = 0; j < x.Length; ++j)
                wa[j] = x[j];

            // The variables STX, FX, DGX contain the values of the step,
            //   function, and directional derivative at the best step.
            // The variables STY, FY, DGY contain the value of the step,
            //   function, and derivative at the other endpoint of the interval
            //   of uncertainty.
            // The variables STP, F, DG contain the values of the step,
            //   function, and derivative at the current step.

            this.stx = 0;
            this.fx = this.finit;
            this.dgx = dginit;
            this.sty = 0;
            this.fy = this.finit;
            this.dgy = dginit;


        L30: // Start of iteration.

            // Set the minimum and maximum steps to correspond
            // to the present interval of uncertainty.

            if (this.brackt)
            {
                this.stmin = Math.Min(this.stx, this.sty);
                this.stmax = Math.Max(this.stx, this.sty);
            }
            else
            {
                this.stmin = this.stx;
                this.stmax = stp + 4 * (stp - this.stx);
            }

            // Force the step to be within
            // the bounds STPMAX and STPMIN.

            stp = Math.Max(stp, this.stpmin);
            stp = Math.Min(stp, this.stpmax);

            // If an unusual termination is to occur then 
            // let STP be the lowest point obtained so far.

            if (this.brackt && (stp <= this.stmin || stp >= this.stmax) || nfev >= this.maxfev - 1 ||
                this.infoc == 0 || this.brackt && this.stmax - this.stmin <= this.xtol * this.stmax)
            {
                stp = this.stx;
            }

            // Evaluate the function and gradient at STP
            // and compute the directional derivative.

            for (int j = 0; j < s.Length; ++j)
                x[j] = wa[j] + stp * s[j];

            // Fetch function and gradient
            f = this.Function(x);
            this.g = this.Gradient(x);

            info = 0;
            nfev++;
            this.dg2 = 0;

            for (int j = 0; j < this.g.Length; ++j)
                this.dg2 += this.g[j] * s[j];

            this.ftest1 = this.finit + stp * this.dgtest;

            if ((this.brackt && (stp <= this.stmin || stp >= this.stmax)) || this.infoc == 0)
                return ConjugateGradientCode.RoundingErrors;

            if (stp == this.stpmax && f <= this.ftest1 && this.dg2 <= this.dgtest)
                return ConjugateGradientCode.StepHigh;

            if (stp == this.stpmin && (f > this.ftest1 || this.dg2 >= this.dgtest))
                return ConjugateGradientCode.StepLow;

            if (nfev >= this.maxfev)
                return ConjugateGradientCode.MaximumEvaluations;

            if (this.brackt && this.stmax - this.stmin <= this.xtol * this.stmax)
                return ConjugateGradientCode.Precision;


            // More's code has been modified so that at least one new 
            //  function value is computed during the line search (enforcing 
            //  at least one interpolation is not easy, since the code may 
            //  override an interpolation) 

            if (f <= this.ftest1 && Math.Abs(this.dg2) <= this.gtol * (-dginit) && nfev > 1)
            {
                info = 1;
                dgout = this.dg2;
                return ConjugateGradientCode.Success;
            }


        L321:

            // In the first stage we seek a step for which the modified
            // function has a nonpositive value and nonnegative derivative.
            if (this.stage1 && f <= this.ftest1 && this.dg2 >= Math.Min(this.ftol, this.gtol) * dginit)
            {
                this.stage1 = false;
            }

            // A modified function is used to predict the step only if
            // we have not obtained a step for which the modified function
            // has a nonpositive function value and nonnegative derivative,
            // and if a lower function value has been obtained but the 
            // decrease is not sufficient.

            if (this.stage1 && f <= this.fx && f > this.ftest1)
            {

                // Define the modified function and derivative values
                double fm = f - stp * this.dgtest;
                double fxm = this.fx - this.stx * this.dgtest;
                double fym = this.fy - this.sty * this.dgtest;
                double dgm = this.dg2 - this.dgtest;
                double dgxm = this.dgx - this.dgtest;
                double dgym = this.dgy - this.dgtest;

                // Call CSTEPM to update the interval of
                // uncertainty and to compute the new step.

                BoundedBroydenFletcherGoldfarbShanno.dcstep(ref this.stx, ref fxm, ref dgxm,
                    ref this.sty, ref fym, ref dgym, ref stp, fm, dgm, ref this.brackt, this.stpmin, this.stpmax);

                // Reset the function and gradient values for f.
                this.fx = fxm + this.stx * this.dgtest;
                this.fy = fym + this.sty * this.dgtest;
                this.dgx = dgxm + this.dgtest;
                this.dgy = dgym + this.dgtest;
            }
            else
            {
                // Call CSTEPM to update the interval of
                // uncertainty and to compute the new step.
                BoundedBroydenFletcherGoldfarbShanno.dcstep(ref this.stx, ref this.fx, ref this.dgx,
                    ref this.sty, ref this.fy, ref this.dgy, ref stp, f, this.dg2, ref this.brackt, this.stpmin, this.stpmax);
            }

            // Force a sufficient decrease in the 
            // size of the interval of uncertainty.

            if (this.brackt)
            {
                if ((Math.Abs(this.sty - this.stx)) >= 0.66 * this.width1)
                    stp = this.stx + 0.5 * (this.sty - this.stx);

                this.width1 = this.width;
                this.width = Math.Abs(this.sty - this.stx);
            }

            goto L30;
        }

    }
}
