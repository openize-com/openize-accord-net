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

// A note on compatibility: Up to version 3.5, users were supposed to implement their own probability
// distributions by inheriting from this class and overriding the public members ProbabilityDensityFunction,
// DistributionFunction, etc. However, since those were public methods this meant that users (and I) had to
// write validation checks on every method override, resulting in lots of duplicated code. Starting from version
// 3.6, users should override methods that start with "Inner" in their name, such as InnerProbabilityDensityFunction 
// and InnerDistributionFunction. The framework will have already validated the inputs of those functions, and
// will also take care to check whether the implementation of those functions is correct.

// For now, compatibility mode is enabled for release builds, meaning that old code that has been written
// using the old way will keep working. However, debug (development) builds will have this feature turned
// off to force new classes to be implemented using this new way.

# if !DEBUG
#define COMPATIBILITY
#endif

namespace Openize.Accord.Statistics.Distributions.Univariate.Base
{
    using System;
    using System.ComponentModel.DataAnnotations;
    using Discrete;
    using Openize.Accord.Math.Matrix;
    using Openize.Accord.Core.Ranges;
    using Openize.Accord.Math.Optimization;
    using Openize.Accord.Math.Random;
    using Openize.Accord.Statistics.Distributions.Fitting.Base;

    /// <summary>
    ///   Abstract class for univariate discrete probability distributions.
    /// </summary>
    /// <remarks>
    /// <para>
    ///   A probability distribution identifies either the probability of each value of an
    ///   unidentified random variable (when the variable is discrete), or the probability
    ///   of the value falling within a particular interval (when the variable is continuous).</para>
    /// <para>
    ///   The probability distribution describes the range of possible values that a random
    ///   variable can attain and the probability that the value of the random variable is
    ///   within any (measurable) subset of that range.</para>  
    /// <para>
    ///   The function describing the probability that a given discrete value will
    ///   occur is called the probability function (or probability mass function,
    ///   abbreviated PMF), and the function describing the cumulative probability
    ///   that a given value or any value smaller than it will occur is called the
    ///   distribution function (or cumulative distribution function, abbreviated CDF).</para>  
    ///   
    /// <para>    
    ///   References:
    ///   <list type="bullet">
    ///     <item><description><a href="http://en.wikipedia.org/wiki/Probability_distribution">
    ///       Wikipedia, The Free Encyclopedia. Probability distribution. Available on:
    ///       http://en.wikipedia.org/wiki/Probability_distribution </a></description></item>
    ///     <item><description><a href="http://mathworld.wolfram.com/StatisticalDistribution.html">
    ///       Weisstein, Eric W. "Statistical Distribution." From MathWorld--A Wolfram Web Resource.
    ///       http://mathworld.wolfram.com/StatisticalDistribution.html </a></description></item>
    ///   </list></para>
    /// </remarks>
    /// 
    /// <seealso cref="BernoulliDistribution"/>
    /// <seealso cref="GeometricDistribution"/>
    /// <seealso cref="PoissonDistribution"/>
    /// 
    [Serializable]
    public abstract class UnivariateDiscreteDistribution : DistributionBase,
        IUnivariateDistribution<int>, IUnivariateDistribution,
        IUnivariateDistribution<double>, IDistribution<double[]>,
        IDistribution<double>,
        ISampleableDistribution<double>, ISampleableDistribution<int>,
        IFormattable
    {

        double? median;
        double? stdDev;
        DoubleRange? quartiles;


        /// <summary>
        ///   Constructs a new UnivariateDistribution class.
        /// </summary>
        /// 
        protected UnivariateDiscreteDistribution()
        {
        }

        /// <summary>
        ///   Gets the mean for this distribution.
        /// </summary>
        /// 
        /// <value>The distribution's mean value.</value>
        /// 
        public abstract double Mean { get; }

        /// <summary>
        ///   Gets the variance for this distribution.
        /// </summary>
        /// 
        /// <value>The distribution's variance.</value>
        /// 
        public abstract double Variance { get; }

        /// <summary>
        ///   Gets the entropy for this distribution.
        /// </summary>
        /// 
        /// <value>The distribution's entropy.</value>
        /// 
        public abstract double Entropy { get; }

        /// <summary>
        ///   Gets the support interval for this distribution.
        /// </summary>
        /// 
        /// <value>A <see cref="IntRange"/> containing
        ///  the support interval for this distribution.</value>
        ///  
        public abstract IntRange Support { get; }

        /// <summary>
        ///   Gets the mode for this distribution.
        /// </summary>
        /// 
        /// <value>The distribution's mode value.</value>
        /// 
        public virtual double Mode
        {
            get { return this.Mean; }
        }

        /// <summary>
        ///   Gets the median for this distribution.
        /// </summary>
        /// 
        /// <value>The distribution's median value.</value>
        /// 
        public virtual double Median
        {
            get
            {
                if (this.median == null)
                    this.median = this.InverseDistributionFunction(0.5);

#if DEBUG
                double expected = this.BaseInverseDistributionFunction(0.5);
                if (this.median != expected)
                    throw new Exception();
#endif

                return this.median.Value;
            }
        }

        /// <summary>
        ///   Gets the Standard Deviation (the square root of
        ///   the variance) for the current distribution.
        /// </summary>
        /// 
        /// <value>The distribution's standard deviation.</value>
        /// 
        public virtual double StandardDeviation
        {
            get
            {
                if (!this.stdDev.HasValue)
                    this.stdDev = Math.Sqrt(this.Variance);
                return this.stdDev.Value;
            }
        }


        /// <summary>
        ///   Gets the Quartiles for this distribution.
        /// </summary>
        /// 
        /// <value>A <see cref="DoubleRange"/> object containing the first quartile
        /// (Q1) as its minimum value, and the third quartile (Q2) as the maximum.</value>
        /// 
        public virtual DoubleRange Quartiles
        {
            get
            {
                if (this.quartiles == null)
                {
                    double min = this.InverseDistributionFunction(0.25);
                    double max = this.InverseDistributionFunction(0.75);
                    this.quartiles = new DoubleRange(min, max);
                }

                return this.quartiles.Value;
            }
        }

        /// <summary>
        ///   Gets the distribution range within a given percentile.
        /// </summary>
        /// 
        /// <remarks>
        ///   If <c>0.25</c> is passed as the <paramref name="percentile"/> argument, 
        ///   this function returns the same as the <see cref="Quartiles"/> function.
        /// </remarks>
        /// 
        /// <param name="percentile">
        ///   The percentile at which the distribution ranges will be returned.</param>
        /// 
        /// <value>A <see cref="DoubleRange"/> object containing the minimum value
        /// for the distribution value, and the third quartile (Q2) as the maximum.</value>
        /// 
        public virtual IntRange GetRange(double percentile)
        {
            if (percentile <= 0 || percentile > 1)
                throw new ArgumentOutOfRangeException("percentile", "The percentile must be between 0 and 1.");

            int a = this.InverseDistributionFunction(1.0 - percentile);
            int b = this.InverseDistributionFunction(percentile);

            if (b > a)
                return new IntRange(a, b);
            return new IntRange(b, a);
        }

        #region IDistribution explicit members

        /// <summary>
        ///   Gets the support interval for this distribution.
        /// </summary>
        /// 
        /// <value>
        ///   A <see cref="DoubleRange" /> containing 
        ///   the support interval for this distribution.
        /// </value>
        /// 
        DoubleRange IUnivariateDistribution.Support
        {
            get { return new DoubleRange(this.Support.Min, this.Support.Max); }
        }

        /// <summary>
        ///   Gets the support interval for this distribution.
        /// </summary>
        /// 
        /// <value>
        ///   A <see cref="DoubleRange" /> containing 
        ///   the support interval for this distribution.
        /// </value>
        /// 
        DoubleRange IUnivariateDistribution<int>.Support
        {
            get { return new DoubleRange(this.Support.Min, this.Support.Max); }
        }

        /// <summary>
        ///   Gets the support interval for this distribution.
        /// </summary>
        /// 
        /// <value>
        ///   A <see cref="DoubleRange" /> containing 
        ///   the support interval for this distribution.
        /// </value>
        /// 
        DoubleRange IUnivariateDistribution<double>.Support
        {
            get { return new DoubleRange(this.Support.Min, this.Support.Max); }
        }

        /// <summary>
        ///   Gets the distribution range within a given percentile.
        /// </summary>
        /// 
        /// <remarks>
        ///   If <c>0.25</c> is passed as the <paramref name="percentile"/> argument, 
        ///   this function returns the same as the <see cref="Quartiles"/> function.
        /// </remarks>
        /// 
        /// <param name="percentile">
        ///   The percentile at which the distribution ranges will be returned.</param>
        /// 
        /// <value>A <see cref="DoubleRange"/> object containing the minimum value
        /// for the distribution value, and the third quartile (Q2) as the maximum.</value>
        /// 
        DoubleRange IUnivariateDistribution.GetRange(double percentile)
        {
            IntRange range = this.GetRange(percentile);
            return new DoubleRange(range.Min, range.Max);
        }

        /// <summary>
        ///   Gets the probability density function (pdf) for
        ///   this distribution evaluated at point <c>x</c>.
        /// </summary>
        /// 
        /// <param name="x">
        ///   A single point in the distribution range. For a 
        ///   univariate distribution, this should be a single
        ///   double value. For a multivariate distribution,
        ///   this should be a double array.</param>
        /// <remarks>
        ///   The Probability Density Function (PDF) describes the
        ///   probability that a given value <c>x</c> will occur.
        /// </remarks>
        /// 
        /// <returns>
        ///   The probability of <c>x</c> occurring
        ///   in the current distribution.</returns>
        ///   
        double IDistribution.DistributionFunction(double[] x)
        {
            if (double.IsNegativeInfinity(x[0]))
                return this.DistributionFunction(int.MinValue);
            else if (double.IsPositiveInfinity(x[0]))
                return this.DistributionFunction(int.MaxValue);

            return this.DistributionFunction((int)x[0]);
        }

        double IUnivariateDistribution.DistributionFunction(double x)
        {
            if (double.IsNegativeInfinity(x))
                return this.DistributionFunction(int.MinValue);
            else if (double.IsPositiveInfinity(x))
                return this.DistributionFunction(int.MaxValue);

            return this.DistributionFunction((int)x);
        }

        double IUnivariateDistribution.DistributionFunction(double a, double b)
        {
            int ia;
            if (double.IsNegativeInfinity(a))
                ia = int.MinValue;
            else if (double.IsPositiveInfinity(a))
                ia = int.MaxValue;
            else ia = (int)a;

            int ib;
            if (double.IsNegativeInfinity(b))
                ib = int.MinValue;
            else if (double.IsPositiveInfinity(b))
                ib = int.MaxValue;
            else ib = (int)b;

            return this.DistributionFunction(ia, ib);
        }

        double IUnivariateDistribution.InverseDistributionFunction(double p)
        {
            return this.InverseDistributionFunction(p);
        }

        double IDistribution.ComplementaryDistributionFunction(double[] x)
        {
            if (double.IsNegativeInfinity(x[0]))
                return this.ComplementaryDistributionFunction(int.MinValue);
            else if (double.IsPositiveInfinity(x[0]))
                return this.ComplementaryDistributionFunction(int.MaxValue);

            return this.ComplementaryDistributionFunction((int)x[0]);
        }

        double IDistribution.ProbabilityFunction(double[] x)
        {
            if (double.IsNegativeInfinity(x[0]))
                return this.ProbabilityMassFunction(int.MinValue);
            else if (double.IsPositiveInfinity(x[0]))
                return this.ProbabilityMassFunction(int.MaxValue);

            return this.ProbabilityMassFunction((int)x[0]);
        }

        double IDistribution.LogProbabilityFunction(double[] x)
        {
            return this.LogProbabilityMassFunction((int)x[0]);
        }

        double IUnivariateDistribution.ProbabilityFunction(double x)
        {
            return this.ProbabilityMassFunction((int)x);
        }

        double IUnivariateDistribution.LogProbabilityFunction(double x)
        {
            return this.LogProbabilityMassFunction((int)x);
        }

        double IUnivariateDistribution.ComplementaryDistributionFunction(double x)
        {
            return this.ComplementaryDistributionFunction((int)x);
        }

        double IUnivariateDistribution.HazardFunction(double x)
        {
            return this.HazardFunction((int)x);
        }

        double IUnivariateDistribution.CumulativeHazardFunction(double x)
        {
            return this.CumulativeHazardFunction((int)x);
        }

        double IUnivariateDistribution.LogCumulativeHazardFunction(double x)
        {
            return this.LogCumulativeHazardFunction((int)x);
        }

        void IDistribution.Fit(Array observations)
        {
            (this as IDistribution).Fit(observations, (IFittingOptions)null);
        }

        void IDistribution.Fit(Array observations, double[] weights)
        {
            (this as IDistribution).Fit(observations, weights, (IFittingOptions)null);
        }

        void IDistribution.Fit(Array observations, int[] weights)
        {
            (this as IDistribution).Fit(observations, weights, (IFittingOptions)null);
        }

        void IDistribution.Fit(Array observations, IFittingOptions options)
        {
            (this as IDistribution).Fit(observations, (double[])null, options);
        }

        void IDistribution.Fit(Array observations, double[] weights, IFittingOptions options)
        {
            double[] univariate = observations as double[];
            if (univariate != null)
            {
                this.Fit(univariate, weights, options);
                return;
            }

            double[][] multivariate = observations as double[][];
            if (multivariate != null)
            {
                var concat = Matrix.Concatenate(multivariate);
                this.Fit(concat, weights, options);
                return;
            }

            int[] iunivariate = observations as int[];
            if (iunivariate != null)
            {
                this.Fit(iunivariate, weights, options);
                return;
            }

            int[][] imultivariate = observations as int[][];
            if (imultivariate != null)
            {
                var concat = Matrix.Concatenate(imultivariate);
                this.Fit(concat, weights, options);
                return;
            }

            throw new ArgumentException("Invalid input type.", "observations");
        }

        void IDistribution.Fit(Array observations, int[] weights, IFittingOptions options)
        {
            double[] univariate = observations as double[];
            if (univariate != null)
            {
                this.Fit(univariate, weights, options);
                return;
            }

            double[][] multivariate = observations as double[][];
            if (multivariate != null)
            {
                this.Fit(Matrix.Concatenate(multivariate), weights, options);
                return;
            }

            int[] iunivariate = observations as int[];
            if (iunivariate != null)
            {
                this.Fit(iunivariate, weights, options);
                return;
            }

            int[][] imultivariate = observations as int[][];
            if (imultivariate != null)
            {
                this.Fit(Matrix.Concatenate(imultivariate), weights, options);
                return;
            }

            throw new ArgumentException("Invalid input type.", "observations");
        }
        #endregion


        /// <summary>
        ///   Gets P(X ≤ k), the cumulative distribution function
        ///   (cdf) for this distribution evaluated at point <c>k</c>.
        /// </summary>
        /// 
        /// <param name="k">
        ///   A single point in the distribution range.</param>
        ///   
        /// <remarks>
        ///   The Cumulative Distribution Function (CDF) describes the cumulative
        ///   probability that a given value or any value smaller than it will occur.
        /// </remarks>
        /// 
        public
#if COMPATIBILITY
        virtual
#endif
        double DistributionFunction(int k)
        {
            if (k < this.Support.Min)
                return 0;

            if (k >= this.Support.Max)
                return 1;

            double result = this.InnerDistributionFunction(k);

            if (Double.IsNaN(result))
                throw new InvalidOperationException("CDF computation generated NaN values.");
            if (result < 0 || result > 1)
                new InvalidOperationException("CDF computation generated values out of the [0,1] range.");

            return result;
        }

        /// <summary>
        ///   Gets P(X ≤ k), the cumulative distribution function
        ///   (cdf) for this distribution evaluated at point <c>k</c>.
        /// </summary>
        /// 
        /// <param name="k">
        ///   A single point in the distribution range.</param>
        ///   
        /// <remarks>
        ///   The Cumulative Distribution Function (CDF) describes the cumulative
        ///   probability that a given value or any value smaller than it will occur.
        /// </remarks>
        /// 
#if COMPATIBILITY
        protected internal virtual double InnerDistributionFunction(int k)
        {
            throw new NotImplementedException();
        }
#else
        protected internal abstract double InnerDistributionFunction(int k);
#endif


        /// <summary>
        ///   Gets P(X ≤ k) or P(X &lt; k), the cumulative distribution function
        ///   (cdf) for this distribution evaluated at point <c>k</c>, depending
        ///   on the value of the <paramref name="inclusive"/> parameter.
        /// </summary>
        /// 
        /// <param name="k">
        ///   A single point in the distribution range.</param>
        /// <param name="inclusive">
        ///   True to return P(X ≤ x), false to return P(X &lt; x). Default is true.</param>
        ///   
        /// <remarks>
        ///   The Cumulative Distribution Function (CDF) describes the cumulative
        ///   probability that a given value or any value smaller than it will occur.
        /// </remarks>
        /// 
        /// <example>
        /// <code>
        /// // Compute P(X = k) 
        /// double equal = dist.ProbabilityMassFunction(k: 1);
        /// 
        /// // Compute P(X &lt; k) 
        /// double less = dist.DistributionFunction(k: 1, inclusive: false);
        /// 
        /// // Compute P(X ≤ k) 
        /// double lessThanOrEqual = dist.DistributionFunction(k: 1, inclusive: true);
        /// 
        /// // Compute P(X > k) 
        /// double greater = dist.ComplementaryDistributionFunction(k: 1);
        /// 
        /// // Compute P(X ≥ k) 
        /// double greaterThanOrEqual = dist.ComplementaryDistributionFunction(k: 1, inclusive: true);
        /// </code>
        /// </example>
        /// 
        /// 
        public
#if COMPATIBILITY
        virtual
#endif
        double DistributionFunction(int k, bool inclusive)
        {
            if (inclusive)
                return this.DistributionFunction(k);
            return this.DistributionFunction(k) - this.ProbabilityMassFunction(k);
        }

        /// <summary>
        ///   Gets the cumulative distribution function (cdf) for this
        ///   distribution in the semi-closed interval (a; b] given as
        ///   <c>P(a &lt; X ≤ b)</c>.
        /// </summary>
        /// 
        /// <param name="a">The start of the semi-closed interval (a; b].</param>
        /// <param name="b">The end of the semi-closed interval (a; b].</param>
        /// 
        /// <remarks>
        ///   The Cumulative Distribution Function (CDF) describes the cumulative
        ///   probability that a given value or any value smaller than it will occur.
        /// </remarks>
        /// 
        public
#if COMPATIBILITY
        virtual
#endif
        double DistributionFunction(int a, int b)
        {
            if (a >= b)
            {
                throw new ArgumentOutOfRangeException("b",
                    "The start of the interval a must be smaller than b.");
            }

            return this.DistributionFunction(b, inclusive: true) - this.DistributionFunction(a);
        }

        /// <summary>
        ///   Gets the inverse of the cumulative distribution function (icdf) for
        ///   this distribution evaluated at probability <c>p</c>. This function 
        ///   is also known as the Quantile function.
        /// </summary>
        /// 
        /// <remarks>
        ///   The Inverse Cumulative Distribution Function (ICDF) specifies, for
        ///   a given probability, the value which the random variable will be at,
        ///   or below, with that probability.
        /// </remarks>
        /// 
        /// <param name="p">A probability value between 0 and 1.</param>
        /// 
        /// <returns>A sample which could original the given probability 
        ///   value when applied in the <see cref="DistributionFunction(int)"/>.</returns>
        /// 
        public int InverseDistributionFunction(
#if !NET35
[Range(0, 1)]
#endif 
            double p)
        {
            if (p < 0.0 || p > 1.0)
                throw new ArgumentOutOfRangeException("p", "Value must be between 0 and 1.");

            if (Double.IsNaN(p))
                throw new ArgumentOutOfRangeException("p", "Value is Not-a-Number (NaN).");

            if (p == 0)
            {
                if (this.Support.Min == this.Support.Max) // Needed by Degenerate
                    return this.Support.Min - 1;
                return this.Support.Min;
            }


            if (p == 1)
                return this.Support.Max;

            int result = this.InnerInverseDistributionFunction(p);

            if (result < this.Support.Min || result > this.Support.Max)
                new InvalidOperationException("invCDF computation generated values outside the distribution supported range.");

            return result;
        }

        /// <summary>
        ///   Gets the inverse of the cumulative distribution function (icdf) for
        ///   this distribution evaluated at probability <c>p</c>. This function 
        ///   is also known as the Quantile function.
        /// </summary>
        /// 
        /// <remarks>
        ///   The Inverse Cumulative Distribution Function (ICDF) specifies, for
        ///   a given probability, the value which the random variable will be at,
        ///   or below, with that probability.
        /// </remarks>
        /// 
        /// <param name="p">A probability value between 0 and 1.</param>
        /// 
        /// <returns>A sample which could original the given probability 
        ///   value when applied in the <see cref="DistributionFunction(int)"/>.</returns>
        /// 
        protected virtual int InnerInverseDistributionFunction(double p)
        {
            return this.BaseInverseDistributionFunction(p);
        }

        /// <summary>
        ///   Gets the inverse of the cumulative distribution function (icdf) for
        ///   this distribution evaluated at probability <c>p</c> using a numerical
        ///   approximation based on binary search.
        /// </summary>
        /// 
        /// <remarks>
        ///   The Inverse Cumulative Distribution Function (ICDF) specifies, for
        ///   a given probability, the value which the random variable will be at,
        ///   or below, with that probability.
        /// </remarks>
        /// 
        /// <param name="p">A probability value between 0 and 1.</param>
        /// 
        /// <returns>A sample which could original the given probability 
        ///   value when applied in the <see cref="DistributionFunction(int)"/>.</returns>
        /// 
        protected int BaseInverseDistributionFunction(double p)
        {
            bool lowerBounded = !Double.IsInfinity(this.Support.Min);
            bool upperBounded = !Double.IsInfinity(this.Support.Max);

            if (lowerBounded && upperBounded)
            {
                return new BinarySearch(this.DistributionFunction, (int)this.Support.Min, (int)this.Support.Max).Find(p);
            }

            try
            {
                checked
                {
                    if (lowerBounded && !upperBounded)
                    {
                        int lower = (int)this.Support.Min;
                        int upper = lower + 1;

                        double f = this.DistributionFunction(lower);

                        if (f > p)
                        {
                            while (f > p)
                            {
                                upper = 2 * upper;
                                f = this.DistributionFunction(upper);
                            }
                        }
                        else
                        {
                            while (f < p)
                            {
                                upper = 2 * upper;
                                f = this.DistributionFunction(upper);
                            }
                        }

                        return new BinarySearch(this.DistributionFunction, lower, upper).Find(p);
                    }

                    if (!lowerBounded && upperBounded)
                    {
                        int upper = (int)this.Support.Max;
                        int lower = upper - 1;

                        double f = this.DistributionFunction(upper);

                        if (f > p)
                        {
                            while (f > p)
                            {
                                lower = lower - 2 * lower;
                                f = this.DistributionFunction(lower);
                            }
                        }
                        else
                        {
                            while (f < p)
                            {
                                lower = lower - 2 * lower;
                                f = this.DistributionFunction(lower);
                            }
                        }

                        return new BinarySearch(this.DistributionFunction, lower, upper).Find(p);
                    }

                    // completely unbounded
                    return this.UnboundedBaseInverseDistributionFunction(p, lower: -1, upper: +1, start: 0);
                }
            }
            catch (OverflowException)
            {
                return 0;
            }
        }

        private int UnboundedBaseInverseDistributionFunction(double p, int lower, int upper, int start)
        {
            checked
            {
                double f = this.DistributionFunction(start);

                if (f > p)
                {
                    while (f > p)
                    {
                        upper = lower;
                        lower = 2 * lower - 1;
                        f = this.DistributionFunction(lower);
                    }
                }
                else
                {
                    while (f < p)
                    {
                        lower = upper;
                        upper = 2 * upper + 1;
                        f = this.DistributionFunction(upper);
                    }
                }

                return new BinarySearch(this.DistributionFunction, lower, upper).Find(p);
            }
        }

        /// <summary>
        ///   Computes the cumulative distribution function by summing the outputs of the <see cref="ProbabilityMassFunction"/>
        ///   for all elements in the distribution domain. Note that this method should not be used in case there is a more
        ///   efficient formula for computing the CDF of a distribution.
        /// </summary>
        /// 
        /// <param name="k">A single point in the distribution range.</param>
        /// 
        protected double BaseDistributionFunction(int k)
        {
            double sum = 0;
            for (int i = this.Support.Min; i <= k; i++)
                sum += this.ProbabilityMassFunction(i);
            return sum;
        }

        /// <summary>
        ///   Gets the first derivative of the <see cref="InverseDistributionFunction">
        ///   inverse distribution function</see> (icdf) for this distribution evaluated
        ///   at probability <c>p</c>. 
        /// </summary>
        /// 
        /// <param name="p">A probability value between 0 and 1.</param>
        /// 
        public virtual double QuantileDensityFunction(double p)
        {
            return 1.0 / this.ProbabilityMassFunction(this.InverseDistributionFunction(p));
        }

        /// <summary>
        ///   Gets the complementary cumulative distribution function
        ///   (ccdf) for this distribution evaluated at point <c>k</c>.
        ///   This function is also known as the Survival function.
        /// </summary>
        /// 
        /// <param name="k">
        ///   A single point in the distribution range.</param>
        /// <param name="inclusive">
        ///   True to return P(X &gt;= x), false to return P(X &gt; x). Default is false.</param>
        ///
        /// <remarks>
        ///   The Complementary Cumulative Distribution Function (CCDF) is
        ///   the complement of the Cumulative Distribution Function, or 1
        ///   minus the CDF.
        /// </remarks>
        /// 
        /// <example>
        /// <code>
        /// // Compute P(X = k) 
        /// double equal = dist.ProbabilityMassFunction(k: 1);
        /// 
        /// // Compute P(X &lt; k) 
        /// double less = dist.DistributionFunction(k: 1, inclusive: false);
        /// 
        /// // Compute P(X ≤ k) 
        /// double lessThanOrEqual = dist.DistributionFunction(k: 1, inclusive: true);
        /// 
        /// // Compute P(X > k) 
        /// double greater = dist.ComplementaryDistributionFunction(k: 1);
        /// 
        /// // Compute P(X ≥ k) 
        /// double greaterThanOrEqual = dist.ComplementaryDistributionFunction(k: 1, inclusive: true);
        /// </code>
        /// </example>
        /// 
        public
#if COMPATIBILITY
        virtual
#endif
        double ComplementaryDistributionFunction(int k, bool inclusive)
        {
            if (inclusive)
                return this.ComplementaryDistributionFunction(k) + this.ProbabilityMassFunction(k);
            return this.ComplementaryDistributionFunction(k);
        }

        /// <summary>
        ///   Gets P(X &gt; k) the complementary cumulative distribution function
        ///   (ccdf) for this distribution evaluated at point <c>k</c>.
        ///   This function is also known as the Survival function.
        /// </summary>
        /// 
        /// <param name="k">
        ///   A single point in the distribution range.</param>
        ///   
        /// <remarks>
        ///   The Complementary Cumulative Distribution Function (CCDF) is
        ///   the complement of the Cumulative Distribution Function, or 1
        ///   minus the CDF.
        /// </remarks>
        /// 
        public
#if COMPATIBILITY
        virtual
#endif
        double ComplementaryDistributionFunction(int k)
        {
            if (k < this.Support.Min)
                return 1;
            if (k >= this.Support.Max)
                return 0;

            double result = this.InnerComplementaryDistributionFunction(k);

            if (Double.IsNaN(result))
                throw new InvalidOperationException("CCDF computation generated NaN values.");
            if (result < 0 || result > 1)
                new InvalidOperationException("CCDF computation generated values out of the [0,1] range.");

            return result;
        }

        /// <summary>
        ///   Gets P(X &gt; k) the complementary cumulative distribution function
        ///   (ccdf) for this distribution evaluated at point <c>k</c>.
        ///   This function is also known as the Survival function.
        /// </summary>
        /// 
        /// <param name="k">
        ///   A single point in the distribution range.</param>
        ///   
        /// <remarks>
        ///   The Complementary Cumulative Distribution Function (CCDF) is
        ///   the complement of the Cumulative Distribution Function, or 1
        ///   minus the CDF.
        /// </remarks>
        /// 
        protected internal virtual double InnerComplementaryDistributionFunction(int k)
        {
            return 1.0 - this.DistributionFunction(k);
        }

        /// <summary>
        ///   Gets the probability mass function (pmf) for
        ///   this distribution evaluated at point <c>x</c>.
        /// </summary>
        /// 
        /// <param name="k">
        ///   A single point in the distribution range.</param>
        ///   
        /// <remarks>
        ///   The Probability Mass Function (PMF) describes the
        ///   probability that a given value <c>x</c> will occur.
        /// </remarks>
        /// 
        /// <returns>
        ///   The probability of <c>k</c> occurring
        ///   in the current distribution.</returns>
        ///   
        public
#if COMPATIBILITY
        virtual
#endif
        double ProbabilityMassFunction(int k)
        {
            if (k < this.Support.Min)
                return 0;
            if (k > this.Support.Max)
                return 0;

            double result = this.InnerProbabilityMassFunction(k);

            if (Double.IsNaN(result))
                throw new InvalidOperationException("PMF computation generated NaN values.");

            return result;
        }

        /// <summary>
        ///   Gets the probability mass function (pmf) for
        ///   this distribution evaluated at point <c>x</c>.
        /// </summary>
        /// 
        /// <param name="k">
        ///   A single point in the distribution range.</param>
        ///   
        /// <remarks>
        ///   The Probability Mass Function (PMF) describes the
        ///   probability that a given value <c>x</c> will occur.
        /// </remarks>
        /// 
        /// <returns>
        ///   The probability of <c>k</c> occurring
        ///   in the current distribution.</returns>
        ///   
#if COMPATIBILITY
        protected internal virtual double InnerProbabilityMassFunction(int k)
        {
            throw new NotImplementedException();
        }
#else
        protected internal abstract double InnerProbabilityMassFunction(int k);
#endif


        /// <summary>
        ///   Gets the log-probability mass function (pmf) for
        ///   this distribution evaluated at point <c>x</c>.
        /// </summary>
        /// 
        /// <param name="k">
        ///   A single point in the distribution range.</param>
        ///   
        /// <remarks>
        ///   The Probability Mass Function (PMF) describes the
        ///   probability that a given value <c>k</c> will occur.
        /// </remarks>
        /// 
        /// <returns>
        ///   The logarithm of the probability of <c>x</c>
        ///   occurring in the current distribution.</returns>
        ///   
        public
#if COMPATIBILITY
        virtual
#endif
        double LogProbabilityMassFunction(int k)
        {
            if (k < this.Support.Min || k > this.Support.Max)
                return Double.NegativeInfinity;

            double result = this.InnerLogProbabilityMassFunction(k);

            if (Double.IsNaN(result))
                throw new InvalidOperationException("LogPDF computation generated NaN values.");

            return result;
        }

        /// <summary>
        ///   Gets the log-probability mass function (pmf) for
        ///   this distribution evaluated at point <c>x</c>.
        /// </summary>
        /// 
        /// <param name="k">
        ///   A single point in the distribution range.</param>
        ///   
        /// <remarks>
        ///   The Probability Mass Function (PMF) describes the
        ///   probability that a given value <c>k</c> will occur.
        /// </remarks>
        /// 
        /// <returns>
        ///   The logarithm of the probability of <c>x</c>
        ///   occurring in the current distribution.</returns>
        ///   
        protected internal virtual double InnerLogProbabilityMassFunction(int k)
        {
            return Math.Log(this.ProbabilityMassFunction(k));
        }

        /// <summary>
        ///   Gets the hazard function, also known as the failure rate or
        ///   the conditional failure density function for this distribution
        ///   evaluated at point <c>x</c>.
        /// </summary>
        /// 
        /// <remarks>
        ///   The hazard function is the ratio of the probability
        ///   density function f(x) to the survival function, S(x).
        /// </remarks>
        /// 
        /// <param name="x">
        ///   A single point in the distribution range.</param>
        ///   
        /// <returns>
        ///   The conditional failure density function <c>h(x)</c>
        ///   evaluated at <c>x</c> in the current distribution.</returns>
        /// 
        public virtual double HazardFunction(int x)
        {
            return this.ProbabilityMassFunction(x) / this.ComplementaryDistributionFunction(x);
        }

        /// <summary>
        ///   Gets the cumulative hazard function for this
        ///   distribution evaluated at point <c>x</c>.
        /// </summary>
        /// 
        /// <param name="x">
        ///   A single point in the distribution range.</param>
        /// 
        /// <returns>
        ///   The cumulative hazard function <c>H(x)</c>  
        ///   evaluated at <c>x</c> in the current distribution.</returns>
        /// 
        public virtual double CumulativeHazardFunction(int x)
        {
            return -Math.Log(this.ComplementaryDistributionFunction(x));
        }

        /// <summary>
        ///   Gets the log-cumulative hazard function for this
        ///   distribution evaluated at point <c>x</c>.
        /// </summary>
        /// 
        /// <param name="x">
        ///   A single point in the distribution range.</param>
        /// 
        /// <returns>
        ///   The logarithm of the cumulative hazard function <c>H(x)</c>  
        ///   evaluated at <c>x</c> in the current distribution.</returns>
        /// 
        public virtual double LogCumulativeHazardFunction(int x)
        {
            return Math.Log(-Math.Log(this.ComplementaryDistributionFunction(x)));
        }

        /// <summary>
        ///   Fits the underlying distribution to a given set of observations.
        /// </summary>
        /// 
        /// <param name="observations">
        ///   The array of observations to fit the model against. The array
        ///   elements can be either of type double (for univariate data) or
        ///   type double[] (for multivariate data).</param>
        ///   
        public virtual void Fit(double[] observations)
        {
            this.Fit(observations, (IFittingOptions)null);
        }

        /// <summary>
        ///   Fits the underlying distribution to a given set of observations.
        /// </summary>
        /// 
        /// <param name="observations">
        ///   The array of observations to fit the model against. The array
        ///   elements can be either of type double (for univariate data) or
        ///   type double[] (for multivariate data).</param>
        ///   
        /// <param name="weights">
        ///   The weight vector containing the weight for each of the samples.</param>
        ///   
        public virtual void Fit(double[] observations, double[] weights)
        {
            this.Fit(observations, weights, (IFittingOptions)null);
        }

        /// <summary>
        ///   Fits the underlying distribution to a given set of observations.
        /// </summary>
        /// 
        /// <param name="observations">
        ///   The array of observations to fit the model against. The array
        ///   elements can be either of type double (for univariate data) or
        ///   type double[] (for multivariate data).</param>
        ///   
        /// <param name="weights">
        ///   The weight vector containing the weight for each of the samples.</param>
        ///   
        public virtual void Fit(double[] observations, int[] weights)
        {
            this.Fit(observations, weights, (IFittingOptions)null);
        }

        /// <summary>
        ///   Fits the underlying distribution to a given set of observations.
        /// </summary>
        /// 
        /// <param name="observations">
        ///   The array of observations to fit the model against. The array
        ///   elements can be either of type double (for univariate data) or
        ///   type double[] (for multivariate data).</param>
        ///   
        /// <param name="options">
        ///   Optional arguments which may be used during fitting, such
        ///   as regularization constants and additional parameters.</param>
        ///   
        public virtual void Fit(double[] observations, IFittingOptions options)
        {
            this.Fit(observations, (double[])null, options);
        }

        /// <summary>
        ///   Fits the underlying distribution to a given set of observations.
        /// </summary>
        /// 
        /// <param name="observations">
        ///   The array of observations to fit the model against. The array
        ///   elements can be either of type double (for univariate data) or
        ///   type double[] (for multivariate data).
        /// </param>
        /// 
        /// <param name="weights">
        ///   The weight vector containing the weight for each of the samples.</param>
        ///   
        /// <param name="options">
        ///   Optional arguments which may be used during fitting, such
        ///   as regularization constants and additional parameters.</param>
        ///   
        public virtual void Fit(double[] observations, double[] weights, IFittingOptions options)
        {
            throw new NotSupportedException();
        }

        /// <summary>
        ///   Fits the underlying distribution to a given set of observations.
        /// </summary>
        /// 
        /// <param name="observations">
        ///   The array of observations to fit the model against. The array
        ///   elements can be either of type double (for univariate data) or
        ///   type double[] (for multivariate data).
        /// </param>
        /// 
        /// <param name="weights">
        ///   The weight vector containing the weight for each of the samples.</param>
        ///   
        /// <param name="options">
        ///   Optional arguments which may be used during fitting, such
        ///   as regularization constants and additional parameters.</param>
        ///   
        public virtual void Fit(double[] observations, int[] weights, IFittingOptions options)
        {
            if (weights == null)
            {
                this.Fit(observations, (double[])null, options);
            }
            else
            {
                throw new NotSupportedException();
            }
        }

        /// <summary>
        ///   Fits the underlying distribution to a given set of observations.
        /// </summary>
        /// 
        /// <param name="observations">
        ///   The array of observations to fit the model against. The array
        ///   elements can be either of type double (for univariate data) or
        ///   type double[] (for multivariate data).</param>
        ///   
        public virtual void Fit(int[] observations)
        {
            this.Fit(observations, (IFittingOptions)null);
        }

        /// <summary>
        ///   Fits the underlying distribution to a given set of observations.
        /// </summary>
        /// 
        /// <param name="observations">
        ///   The array of observations to fit the model against. The array
        ///   elements can be either of type double (for univariate data) or
        ///   type double[] (for multivariate data).</param>
        ///   
        /// <param name="options">
        ///   Optional arguments which may be used during fitting, such
        ///   as regularization constants and additional parameters.</param>
        ///   
        public virtual void Fit(int[] observations, IFittingOptions options)
        {
            this.Fit(observations, (double[])null, options);
        }

        /// <summary>
        ///   Fits the underlying distribution to a given set of observations.
        /// </summary>
        /// 
        /// <param name="observations">
        ///   The array of observations to fit the model against. The array
        ///   elements can be either of type double (for univariate data) or
        ///   type double[] (for multivariate data).
        /// </param>
        /// 
        /// <param name="weights">
        ///   The weight vector containing the weight for each of the samples.</param>
        ///   
        /// <param name="options">
        ///   Optional arguments which may be used during fitting, such
        ///   as regularization constants and additional parameters.</param>
        ///   
        public virtual void Fit(int[] observations, double[] weights, IFittingOptions options)
        {
            throw new NotSupportedException();
        }

        /// <summary>
        ///   Fits the underlying distribution to a given set of observations.
        /// </summary>
        /// 
        /// <param name="observations">
        ///   The array of observations to fit the model against. The array
        ///   elements can be either of type double (for univariate data) or
        ///   type double[] (for multivariate data).
        /// </param>
        /// 
        /// <param name="weights">
        ///   The weight vector containing the weight for each of the samples.</param>
        ///   
        /// <param name="options">
        ///   Optional arguments which may be used during fitting, such
        ///   as regularization constants and additional parameters.</param>
        ///   
        public virtual void Fit(int[] observations, int[] weights, IFittingOptions options)
        {
            if (weights != null)
                throw new NotSupportedException();

            this.Fit(observations, (double[])null, options);
        }



        /// <summary>
        ///   Generates a random vector of observations from the current distribution.
        /// </summary>
        /// 
        /// <param name="samples">The number of samples to generate.</param>
        /// <returns>A random vector of observations drawn from this distribution.</returns>
        /// 
        public int[] Generate(int samples)
        {
            return this.Generate(samples, global::Openize.Accord.Math.Random.Generator.Random);
        }

        /// <summary>
        ///   Generates a random vector of observations from the current distribution.
        /// </summary>
        /// 
        /// <param name="samples">The number of samples to generate.</param>
        /// <param name="result">The location where to store the samples.</param>
        /// 
        /// <returns>A random vector of observations drawn from this distribution.</returns>
        /// 
        public int[] Generate(int samples, int[] result)
        {
            return this.Generate(samples, result, global::Openize.Accord.Math.Random.Generator.Random);
        }

        /// <summary>
        ///   Generates a random vector of observations from the current distribution.
        /// </summary>
        /// 
        /// <param name="samples">The number of samples to generate.</param>
        /// <param name="result">The location where to store the samples.</param>
        ///
        /// <returns>A random vector of observations drawn from this distribution.</returns>
        /// 
        public double[] Generate(int samples, double[] result)
        {
            return this.Generate(samples, result, global::Openize.Accord.Math.Random.Generator.Random);
        }

        /// <summary>
        ///   Generates a random observation from the current distribution.
        /// </summary>
        /// 
        /// <returns>A random observations drawn from this distribution.</returns>
        /// 
        public int Generate()
        {
            return this.Generate(global::Openize.Accord.Math.Random.Generator.Random);
        }



        /// <summary>
        ///   Generates a random vector of observations from the current distribution.
        /// </summary>
        /// 
        /// <param name="samples">The number of samples to generate.</param>
        /// <param name="source">The random number generator to use as a source of randomness. 
        ///   Default is to use <see cref="Generator.Random"/>.</param>
        ///   
        /// <returns>A random vector of observations drawn from this distribution.</returns>
        /// 
        public int[] Generate(int samples, Random source)
        {
            return this.Generate(samples, new int[samples]);
        }

        /// <summary>
        ///   Generates a random vector of observations from the current distribution.
        /// </summary>
        /// 
        /// <param name="samples">The number of samples to generate.</param>
        /// <param name="result">The location where to store the samples.</param>
        /// <param name="source">The random number generator to use as a source of randomness. 
        ///   Default is to use <see cref="Generator.Random"/>.</param>
        /// 
        /// <returns>A random vector of observations drawn from this distribution.</returns>
        /// 
        public virtual int[] Generate(int samples, int[] result, Random source)
        {
            for (int i = 0; i < samples; i++)
                result[i] = this.InverseDistributionFunction(source.NextDouble());
            return result;
        }

        /// <summary>
        ///   Generates a random vector of observations from the current distribution.
        /// </summary>
        /// 
        /// <param name="samples">The number of samples to generate.</param>
        /// <param name="result">The location where to store the samples.</param>
        /// <param name="source">The random number generator to use as a source of randomness. 
        ///   Default is to use <see cref="Generator.Random"/>.</param>
        ///   
        /// <returns>A random vector of observations drawn from this distribution.</returns>
        /// 
        public virtual double[] Generate(int samples, double[] result, Random source)
        {
            for (int i = 0; i < samples; i++)
                result[i] = this.InverseDistributionFunction(source.NextDouble());
            return result;
        }

        /// <summary>
        ///   Generates a random observation from the current distribution.
        /// </summary>
        /// 
        /// <param name="source">The random number generator to use as a source of randomness. 
        ///   Default is to use <see cref="Generator.Random"/>.</param>
        ///   
        /// <returns>A random observations drawn from this distribution.</returns>
        /// 
        public virtual int Generate(Random source)
        {
            return this.InverseDistributionFunction(source.NextDouble());
        }




        double[] IRandomNumberGenerator<double>.Generate(int samples)
        {
            return this.Generate(samples, new double[samples]);
        }

        double[] ISampleableDistribution<double>.Generate(int samples, Random source)
        {
            return this.Generate(samples, new double[samples], source);
        }

        double ISampleableDistribution<double>.Generate(double result)
        {
            return this.Generate();
        }

        int ISampleableDistribution<int>.Generate(int result)
        {
            return this.Generate();
        }

        double ISampleableDistribution<double>.Generate(double result, Random source)
        {
            return this.Generate();
        }

        int ISampleableDistribution<int>.Generate(int result, Random source)
        {
            return this.Generate();
        }

        double IRandomNumberGenerator<double>.Generate()
        {
            return (double)this.Generate();
        }

        double ISampleableDistribution<double>.Generate(Random source)
        {
            return (double)this.Generate(source);
        }


        double IDistribution<int>.ProbabilityFunction(int x)
        {
            return this.ProbabilityMassFunction(x);
        }

        double IDistribution<int>.LogProbabilityFunction(int x)
        {
            return this.LogProbabilityMassFunction(x);
        }

        double IDistribution<double[]>.DistributionFunction(double[] x)
        {
            return (this as IDistribution).DistributionFunction(x);
        }

        double IDistribution<double[]>.ProbabilityFunction(double[] x)
        {
            return (this as IDistribution).ProbabilityFunction(x);
        }

        double IDistribution<double[]>.LogProbabilityFunction(double[] x)
        {
            return (this as IDistribution).LogProbabilityFunction(x);
        }

        double IDistribution<double[]>.ComplementaryDistributionFunction(double[] x)
        {
            return (this as IDistribution).ComplementaryDistributionFunction(x);
        }

        double IDistribution<double>.DistributionFunction(double x)
        {
            return this.DistributionFunction((int)x);
        }

        double IDistribution<double>.ProbabilityFunction(double x)
        {
            return this.ProbabilityMassFunction((int)x);
        }

        double IDistribution<double>.LogProbabilityFunction(double x)
        {
            return this.LogProbabilityMassFunction((int)x);
        }

        double IDistribution<double>.ComplementaryDistributionFunction(double x)
        {
            return this.ComplementaryDistributionFunction((int)x);
        }

        double IUnivariateDistribution<double>.InverseDistributionFunction(double p)
        {
            return this.InverseDistributionFunction(p);
        }

        double IUnivariateDistribution<double>.HazardFunction(double x)
        {
            return this.HazardFunction((int)x);
        }

        double IUnivariateDistribution<double>.CumulativeHazardFunction(double x)
        {
            return this.CumulativeHazardFunction((int)x);
        }

    }
}