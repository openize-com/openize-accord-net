// Accord Math Library
// The Accord.NET Framework
// http://accord-framework.net
//
// Copyright © Steven G. Johnson, 2008
// stevenj@alum.mit.edu
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
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
//
//
// The source code presented in this file has been adapted from the
// Sbplx method (based on Nelder-Mead's Simplex) given in the NLopt
// Numerical Optimization Library. Original license is given below.
//
//    Copyright (c) 2007-2011 Massachusetts Institute of Technology
//
//    Permission is hereby granted, free of charge, to any person obtaining
//    a copy of this software and associated documentation files (the
//    "Software"), to deal in the Software without restriction, including
//    without limitation the rights to use, copy, modify, merge, publish,
//    distribute, sublicense, and/or sell copies of the Software, and to
//    permit persons to whom the Software is furnished to do so, subject to
//    the following conditions:
// 
//    The above copyright notice and this permission notice shall be
//    included in all copies or substantial portions of the Software.
// 
//    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
//    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
//    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
//    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
//    LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
//    OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
//    WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
//

namespace Openize.Accord.Math.Optimization.Unconstrained
{
    using System;
    using System.Diagnostics;
    using Convergence;
    using Openize.Accord.Math;
    using Openize.Accord.Math.Optimization.Base;

    /// <summary>
    ///   Subplex
    /// </summary>
    /// 
    /// <remarks>
    /// <para>
    ///   The source code presented in this file has been adapted from the
    ///   Sbplx method (based on Nelder-Mead's Simplex) given in the NLopt
    ///   Numerical Optimization Library, created by Steven G. Johnson.</para>
    ///   
    /// <para>
    ///   References:
    ///   <list type="bullet">
    ///     <item><description><a href="http://ab-initio.mit.edu/nlopt">
    ///       Steven G. Johnson, The NLopt nonlinear-optimization package, 
    ///       http://ab-initio.mit.edu/nlopt </a></description></item>
    ///     <item><description><a href="http://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method">
    ///       Wikipedia, The Free Encyclopedia. Nelder Mead method. Available on:
    ///       http://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method </a></description></item>
    ///    </list></para>
    /// </remarks>
    /// 
    public class Subplex : BaseOptimizationMethod, IOptimizationMethod<NelderMeadStatus>
    {

        private int n;

        private double minf_max;
        private GeneralConvergence stop;


        // subplex strategy constants:
        const double psi = 0.25;
        const double omega = 0.1;
        const int nsmin = 2;
        const int nsmax = 5;

        // bounds
        private double[] lb;
        private double[] ub;

        private double[] xprev;
        private double[] dx;
        private double[] absdx;

        private double[] xstep; // initial step sizes
        private double[] xstep0;


        private int[] p; // subspace index permutation
        private int sindex; // starting index for this subspace


        NelderMead nelderMead;
        NelderMeadStatus status;


        /// <summary>
        ///   Creates a new <see cref="Subplex"/> optimization algorithm.
        /// </summary>
        /// 
        /// <param name="numberOfVariables">The number of free parameters in the optimization problem.</param>
        /// 
        public Subplex(int numberOfVariables)
            : base(numberOfVariables)
        {
            this.init(numberOfVariables);
        }

        /// <summary>
        ///   Creates a new <see cref="Subplex"/> optimization algorithm.
        /// </summary>
        /// 
        /// <param name="numberOfVariables">The number of free parameters in the optimization problem.</param>
        /// <param name="function">The objective function whose optimum values should be found.</param>
        /// 
        public Subplex(int numberOfVariables, Func<double[], double> function)
            : base(numberOfVariables, function)
        {
            this.init(numberOfVariables);
        }

        /// <summary>
        ///   Creates a new <see cref="Subplex"/> optimization algorithm.
        /// </summary>
        /// 
        /// <param name="function">The objective function whose optimum values should be found.</param>
        /// 
        public Subplex(NonlinearObjectiveFunction function)
            : base(function)
        {
            this.init(function.NumberOfVariables);
        }

        private void init(int n)
        {
            this.n = n;
            this.stop = new GeneralConvergence(n);

            this.nelderMead = new NelderMead(nsmax, this.subspace_func);
            this.nelderMead.Convergence = this.stop;

            this.xstep = new double[n];
            this.xstep0 = new double[n];
            for (int i = 0; i < this.xstep0.Length; i++)
                this.xstep0[i] = 1e-5;

            this.p = new int[n];
            this.dx = new double[n];
            this.xprev = new double[n];
            this.absdx = new double[n];

            this.lb = new double[n];
            for (int i = 0; i < this.lb.Length; i++)
                this.lb[i] = Double.NegativeInfinity;

            this.ub = new double[n];
            for (int i = 0; i < this.ub.Length; i++)
                this.ub[i] = Double.PositiveInfinity;
        }

        /// <summary>
        ///   Get the exit code returned in the last call to the
        ///   <see cref="IOptimizationMethod{TInput, TOutput}.Maximize()"/> or 
        ///   <see cref="IOptimizationMethod{TInput, TOutput}.Minimize()"/> methods.
        /// </summary>
        /// 
        public NelderMeadStatus Status
        {
            get { return this.status; }
        }

        /// <summary>
        ///   Gets or sets the maximum value that the objective 
        ///   function could produce before the algorithm could
        ///   be terminated as if the solution was good enough.
        /// </summary>
        /// 
        public double MaximumValue
        {
            get { return this.minf_max; }
            set { this.minf_max = value; }
        }

        /// <summary>
        ///   Gets the step sizes to be used by the optimization
        ///   algorithm. Default is to initialize each with 1e-5.
        /// </summary>
        /// 
        public double[] StepSize
        {
            get { return this.xstep0; }
        }

        /// <summary>
        ///   Gets or sets multiple convergence options to 
        ///   determine when the optimization can terminate.
        /// </summary>
        /// 
        public GeneralConvergence Convergence
        {
            get { return this.stop; }
            set { this.stop = value; }
        }

        /// <summary>
        ///   Gets the lower bounds that should be respected in this 
        ///   optimization problem. Default is to initialize this vector
        ///   with <see cref="Double.NegativeInfinity"/>.
        /// </summary>
        /// 
        public double[] LowerBounds
        {
            get { return this.lb; }
        }

        /// <summary>
        ///   Gets the upper bounds that should be respected in this 
        ///   optimization problem. Default is to initialize this vector
        ///   with <see cref="Double.PositiveInfinity"/>.
        /// </summary>
        /// 
        public double[] UpperBounds
        {
            get { return this.ub; }
        }



        /// <summary>
        ///   Implements the actual optimization algorithm. This
        ///   method should try to minimize the objective function.
        /// </summary>
        /// 
        protected override bool Optimize()
        {
            this.status = this.sbplx_minimize();

            return (this.status == NelderMeadStatus.Success ||
                this.status == NelderMeadStatus.FunctionToleranceReached ||
                this.status == NelderMeadStatus.SolutionToleranceReached);
        }

        NelderMeadStatus sbplx_minimize()
        {
            var ret = NelderMeadStatus.Success;

            double[] x = this.Solution;
            this.Value = this.Function(x);

            this.stop.Evaluations++;
            if (NelderMead.nlopt_stop_forced(this.stop))
                return NelderMeadStatus.ForcedStop;
            if (this.Value < this.minf_max)
                return NelderMeadStatus.MinimumAllowedValueReached;
            if (NelderMead.nlopt_stop_evals(this.stop))
                return NelderMeadStatus.MaximumEvaluationsReached;
            if (NelderMead.nlopt_stop_time(this.stop))
                return NelderMeadStatus.MaximumTimeReached;


            Array.Copy(this.xstep0, this.xstep, this.xstep.Length);


            while (true)
            {
                double normi = 0;
                double normdx = 0;
                int ns, nsubs = 0;
                int nevals = this.stop.Evaluations;
                double fdiff, fdiff_max = 0;

                Array.Copy(x, this.xprev, x.Length);

                double fprev = this.Value;

                // sort indices into the progress vector dx
                // by decreasing order of magnitude abs(dx)
                //
                for (int i = 0; i < this.p.Length; i++)
                    this.p[i] = i;

                for (int j = 0; j < this.absdx.Length; j++)
                    this.absdx[j] = Math.Abs(this.dx[j]);

                Array.Sort(this.p, this.absdx);


                // find the subspaces, and perform nelder-mead on each one
                for (int i = 0; i < this.absdx.Length; i++)
                    normdx += this.absdx[i]; // L1 norm

                int last = 0;
                for (int i = 0; i + nsmin < this.n; i += ns)
                {
                    last = i;

                    // find subspace starting at index i
                    double ns_goodness = -Double.MaxValue;
                    double norm = normi;
                    int nk = i + nsmax > this.n ? this.n : i + nsmax; // max k for this subspace

                    for (int k = i; k < i + nsmin - 1; k++)
                        norm += this.absdx[this.p[k]];

                    ns = nsmin;
                    for (int k = i + nsmin - 1; k < nk; k++)
                    {
                        double goodness;
                        norm += this.absdx[this.p[k]];

                        // remaining subspaces must be big enough to partition
                        if (this.n - (k + 1) < nsmin)
                            continue;

                        // maximize figure of merit defined by Rowan thesis:
                        // look for sudden drops in average |dx|

                        if (k + 1 < this.n)
                        {
                            goodness = norm / (k + 1) - (normdx - norm) / (this.n - (k + 1));
                        }
                        else
                        {
                            goodness = normdx / this.n;
                        }

                        if (goodness > ns_goodness)
                        {
                            ns_goodness = goodness;
                            ns = (k + 1) - i;
                        }
                    }

                    for (int k = i; k < i + ns; ++k)
                        normi += this.absdx[this.p[k]];

                    // do nelder-mead on subspace of dimension ns starting w/i 
                    this.sindex = i;
                    for (int k = i; k < i + ns; ++k)
                    {
                        this.nelderMead.Solution[k - i] = x[this.p[k]];
                        this.nelderMead.StepSize[k - i] = this.xstep[this.p[k]];
                        this.nelderMead.LowerBounds[k - i] = this.lb[this.p[k]];
                        this.nelderMead.UpperBounds[k - i] = this.ub[this.p[k]];
                    }

                    nsubs++;
                    nevals = this.stop.Evaluations;

                    this.nelderMead.NumberOfVariables = ns;
                    this.nelderMead.DiameterTolerance = psi;
                    ret = this.nelderMead.Minimize(this.Value);

                    fdiff = this.nelderMead.Difference;
                    this.Value = this.nelderMead.Value;

                    if (fdiff > fdiff_max)
                        fdiff_max = fdiff;

                    Trace.WriteLine(String.Format("{0} NM iterations for ({1},{2}) subspace",
                       this.stop.Evaluations - nevals, this.sindex, ns));

                    for (int k = i; k < i + ns; k++)
                        x[this.p[k]] = this.nelderMead.Solution[k - i];

                    if (ret == NelderMeadStatus.Failure)
                        return NelderMeadStatus.SolutionToleranceReached;

                    if (ret != NelderMeadStatus.SolutionToleranceReached)
                        return ret;
                }

                // nelder-mead on last subspace 
                ns = this.n - last;
                this.sindex = last;
                for (int i = last; i < this.n; i++)
                {
                    this.nelderMead.Solution[i - this.sindex] = x[this.p[i]];
                    this.nelderMead.StepSize[i - this.sindex] = this.xstep[this.p[i]];
                    this.nelderMead.LowerBounds[i - this.sindex] = this.lb[this.p[i]];
                    this.nelderMead.UpperBounds[i - this.sindex] = this.ub[this.p[i]];
                }

                nsubs++;
                nevals = this.stop.Evaluations;

                this.nelderMead.NumberOfVariables = ns;
                this.nelderMead.DiameterTolerance = psi;
                ret = this.nelderMead.Minimize(this.Value);

                fdiff = this.nelderMead.Difference;
                this.Value = this.nelderMead.Value;

                if (fdiff > fdiff_max)
                    fdiff_max = fdiff;

                Trace.WriteLine(String.Format("sbplx: {0} NM iterations for ({1},{2}) subspace",
                   this.stop.Evaluations - nevals, this.sindex, ns));


                for (int i = this.sindex; i < this.p.Length; i++)
                    x[this.p[i]] = this.nelderMead.Solution[i - this.sindex];

                if (ret == NelderMeadStatus.Failure)
                    return NelderMeadStatus.SolutionToleranceReached;

                if (ret != NelderMeadStatus.SolutionToleranceReached)
                    return ret;

                // termination tests:
                if (NelderMead.nlopt_stop_ftol(this.stop, this.Value, this.Value + fdiff_max))
                    return NelderMeadStatus.FunctionToleranceReached;

                if (NelderMead.nlopt_stop_xtol(this.stop, x, this.xprev, this.n))
                {
                    int j;

                    // as explained in Rowan's thesis, it is important
                    // to check |xstep| as well as |x-xprev|, since if
                    // the step size is too large (in early iterations),
                    // the inner Nelder-Mead may not make much progress 
                    //
                    for (j = 0; j < this.xstep.Length; j++)
                    {
                        if (Math.Abs(this.xstep[j]) * psi > this.stop.AbsoluteParameterTolerance[j]
                         && Math.Abs(this.xstep[j]) * psi > this.stop.RelativeParameterTolerance * Math.Abs(x[j]))
                            break;
                    }

                    if (j == this.n)
                    {
                        return NelderMeadStatus.SolutionToleranceReached;
                    }
                }

                // compute change in optimal point
                for (int i = 0; i < x.Length; i++)
                    this.dx[i] = x[i] - this.xprev[i];

                // setting step sizes
                {
                    double scale;
                    if (nsubs == 1)
                    {
                        scale = psi;
                    }
                    else
                    {
                        double stepnorm = 0, dxnorm = 0;
                        for (int i = 0; i < this.dx.Length; i++)
                        {
                            stepnorm += Math.Abs(this.xstep[i]);
                            dxnorm += Math.Abs(this.dx[i]);
                        }

                        scale = dxnorm / stepnorm;

                        if (scale < omega)
                            scale = omega;

                        if (scale > 1 / omega)
                            scale = 1 / omega;
                    }


                    Trace.WriteLine("sbplx: stepsize scale factor = " + scale);


                    for (int i = 0; i < this.xstep.Length; i++)
                    {
                        this.xstep[i] = (this.dx[i] == 0) ?
                            -(this.xstep[i] * scale) : Special.Sign(this.xstep[i] * scale, this.dx[i]);
                    }
                }
            }
        }


        /// <summary>
        ///   Wrapper around objective function for subspace optimization.
        /// </summary>
        /// 
        double subspace_func(double[] xs)
        {
            double[] x = this.Solution;

            for (int i = this.sindex; i < this.sindex + this.n; i++)
                x[this.p[i]] = xs[i - this.sindex];

            return this.Function(x);
        }

    }
}
