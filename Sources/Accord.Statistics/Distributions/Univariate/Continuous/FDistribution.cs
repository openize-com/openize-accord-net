﻿// Accord Statistics Library
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

namespace Openize.Accord.Statistics.Distributions.Univariate.Continuous
{
    using System;
    using Openize.Accord.Core.Attributes;
    using Openize.Accord.Core.Ranges;
    using Openize.Accord.Math.Functions;
    using Openize.Accord.Math.Random;
    using Openize.Accord.Statistics.Distributions.Fitting.Base;
    using Openize.Accord.Statistics.Distributions.Univariate.Base;
    using Testing;

    /// <summary>
    ///   F (Fisher-Snedecor) distribution.
    /// </summary>
    /// 
    /// <remarks>
    /// <para>
    ///   In probability theory and statistics, the F-distribution is a continuous
    ///   probability distribution. It is also known as Snedecor's F distribution
    ///   or the Fisher-Snedecor distribution (after R.A. Fisher and George W. Snedecor). 
    ///   The F-distribution arises frequently as the null distribution of a test statistic,
    ///   most notably in the analysis of variance; see <see cref="FTest"/>.</para>
    ///   
    /// <para>
    ///   References:
    ///   <list type="bullet">
    ///     <item><description><a href="http://en.wikipedia.org/wiki/F-distribution">
    ///       Wikipedia, The Free Encyclopedia. F-distribution. Available on: 
    ///       http://en.wikipedia.org/wiki/F-distribution </a></description></item>
    ///   </list></para>     
    /// </remarks>
    /// 
    /// <example>
    /// <para>
    ///   The following example shows how to construct a Fisher-Snedecor's F-distribution
    ///   with 8 and 5 degrees of freedom, respectively.</para>
    ///
    /// <code>
    ///   // Create a Fisher-Snedecor's F distribution with 8 and 5 d.f.
    ///   FDistribution F = new FDistribution(degrees1: 8, degrees2: 5);
    /// 
    ///   // Common measures
    ///   double mean = F.Mean;     // 1.6666666666666667
    ///   double median = F.Median; // 1.0545096252132447
    ///   double var = F.Variance;  // 7.6388888888888893
    ///   
    ///   // Cumulative distribution functions
    ///   double cdf = F.DistributionFunction(x: 0.27); // 0.049463408057268315
    ///   double ccdf = F.ComplementaryDistributionFunction(x: 0.27); // 0.95053659194273166
    ///   double icdf = F.InverseDistributionFunction(p: cdf); // 0.27
    ///   
    ///   // Probability density functions
    ///   double pdf = F.ProbabilityDensityFunction(x: 0.27); // 0.45120469723580559
    ///   double lpdf = F.LogProbabilityDensityFunction(x: 0.27); // -0.79583416831212883
    ///   
    ///   // Hazard (failure rate) functions
    ///   double hf = F.HazardFunction(x: 0.27); // 0.47468419528555084
    ///   double chf = F.CumulativeHazardFunction(x: 0.27); // 0.050728620222091653
    ///   
    ///   // String representation
    ///   string str = F.ToString(CultureInfo.InvariantCulture); // F(x; df1 = 8, df2 = 5)
    /// </code>
    /// </example>
    /// 
    [Serializable]
    public class FDistribution : UnivariateContinuousDistribution,
        ISampleableDistribution<double>
    {

        // Distribution parameters
        private int d1;
        private int d2;

        // Derived values
        private double b;

        private double? mean;
        private double? variance;

        /// <summary>
        ///   Constructs a F-distribution with
        ///   the given degrees of freedom.
        /// </summary>
        /// 
        public FDistribution()
            : this(1, 1)
        {
        }

        /// <summary>
        ///   Constructs a F-distribution with
        ///   the given degrees of freedom.
        /// </summary>
        /// 
        /// <param name="degrees1">The first degree of freedom. Default is 1.</param>
        /// <param name="degrees2">The second degree of freedom. Default is 1.</param>
        /// 
        public FDistribution([PositiveInteger] int degrees1, [PositiveInteger] int degrees2)
        {
            if (degrees1 <= 0)
                throw new ArgumentOutOfRangeException("degrees1", "Degrees of freedom must be positive.");

            if (degrees2 <= 0)
                throw new ArgumentOutOfRangeException("degrees1", "Degrees of freedom must be positive.");

            this.d1 = degrees1;
            this.d2 = degrees2;

            this.b = Beta.Function(degrees1 * 0.5, degrees2 * 0.5);
        }

        /// <summary>
        ///   Gets the first degree of freedom.
        /// </summary>
        /// 
        public int DegreesOfFreedom1
        {
            get { return this.d1; }
        }

        /// <summary>
        ///   Gets the second degree of freedom.
        /// </summary>
        /// 
        public int DegreesOfFreedom2
        {
            get { return this.d2; }
        }

        /// <summary>
        ///   Gets the mean for this distribution.
        /// </summary>
        /// 
        public override double Mean
        {
            get
            {
                if (!this.mean.HasValue)
                {
                    if (this.d2 <= 2)
                    {
                        this.mean = Double.NaN;
                    }
                    else
                    {
                        this.mean = this.d2 / (this.d2 - 2.0);
                    }
                }

                return this.mean.Value;
            }
        }

        /// <summary>
        ///   Gets the variance for this distribution.
        /// </summary>
        /// 
        public override double Variance
        {
            get
            {
                if (!this.variance.HasValue)
                {
                    if (this.d2 <= 4)
                    {
                        this.variance = Double.NaN;
                    }
                    else
                    {
                        this.variance = (2.0 * this.d2 * this.d2 * (this.d1 + this.d2 - 2)) /
                                        (this.d1 * (this.d2 - 2) * (this.d2 - 2) * (this.d2 - 4));
                    }
                }

                return this.variance.Value;
            }
        }

        /// <summary>
        ///   Gets the mode for this distribution.
        /// </summary>
        /// 
        /// <value>
        ///   The distribution's mode value.
        /// </value>
        /// 
        public override double Mode
        {
            get
            {
                if (this.d1 > 2)
                {
                    double a = (this.d1 - 2.0) / this.d1;
                    double b = this.d2 / (this.d2 + 2.0);
                    return a * b;
                }

                return Double.NaN;
            }
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
        public override DoubleRange Support
        {
            get { return new DoubleRange(0, Double.PositiveInfinity); }
        }

        /// <summary>
        ///   Gets the entropy for this distribution.
        /// </summary>
        /// 
        public override double Entropy
        {
            get { throw new NotSupportedException(); }
        }

        /// <summary>
        ///   Gets the cumulative distribution function (cdf) for
        ///   the F-distribution evaluated at point <c>x</c>.
        /// </summary>
        /// 
        /// <param name="x">A single point in the distribution range.</param>
        /// 
        /// <remarks>
        /// <para>
        ///   The Cumulative Distribution Function (CDF) describes the cumulative
        ///   probability that a given value or any value smaller than it will occur.</para>
        /// <para>
        ///   The F-distribution CDF is computed in terms of the <see cref="Beta.Incomplete">
        ///   Incomplete Beta function Ix(a,b)</see> as CDF(x) = Iu(d1/2, d2/2) in which
        ///   u is given as u = (d1 * x) / (d1 * x + d2)</para>.
        /// </remarks>
        /// 
        protected internal override double InnerDistributionFunction(double x)
        {
            double u = (this.d1 * x) / (this.d1 * x + this.d2);
            return Beta.Incomplete(this.d1 * 0.5, this.d2 * 0.5, u);
        }

        /// <summary>
        ///   Gets the complementary cumulative distribution
        ///   function evaluated at point <c>x</c>.
        /// </summary>
        /// 
        /// <remarks>
        /// <para>
        ///   The F-distribution complementary CDF is computed in terms of the <see cref="Beta.Incomplete">
        ///   Incomplete Beta function Ix(a,b)</see> as CDFc(x) = Iu(d2/2, d1/2) in which
        ///   u is given as u = (d2 * x) / (d2 * x + d1)</para>.
        /// </remarks>
        /// 
        protected internal override double InnerComplementaryDistributionFunction(double x)
        {
            double u = this.d2 / (this.d2 + this.d1 * x);
            return Beta.Incomplete(this.d2 * 0.5, this.d1 * 0.5, u);
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
        protected internal override double InnerInverseDistributionFunction(double p)
        {
            // Cephes Math Library Release 2.8:  June, 2000
            // Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
            // Adapted under the LGPL with permission of original author.

            double d1 = this.d1;
            double d2 = this.d2;

            double x;

            double w = Beta.Incomplete(0.5 * d2, 0.5 * d1, 0.5);

            if (w > p || p < 0.001)
            {
                w = Beta.IncompleteInverse(0.5 * d1, 0.5 * d2, p);
                x = d2 * w / (d1 * (1.0 - w));
            }
            else
            {
                w = Beta.IncompleteInverse(0.5 * d2, 0.5 * d1, 1.0 - p);
                x = (d2 - d2 * w) / (d1 * w);
            }

            return x;
        }

        /// <summary>
        ///   Gets the probability density function (pdf) for
        ///   the F-distribution evaluated at point <c>x</c>.
        /// </summary>
        /// 
        /// <param name="x">A single point in the distribution range.</param>
        /// 
        /// <returns>
        ///   The probability of <c>x</c> occurring
        ///   in the current distribution.
        /// </returns>
        /// 
        /// <remarks>
        ///   The Probability Density Function (PDF) describes the
        ///   probability that a given value <c>x</c> will occur.
        /// </remarks>
        /// 
        protected internal override double InnerProbabilityDensityFunction(double x)
        {
            double u = Math.Pow(this.d1 * x, this.d1) * Math.Pow(this.d2, this.d2) /
                Math.Pow(this.d1 * x + this.d2, this.d1 + this.d2);
            return Math.Sqrt(u) / (x * this.b);
        }

        /// <summary>
        ///   Gets the log-probability density function (pdf) for
        ///   this distribution evaluated at point <c>x</c>.
        /// </summary>
        /// 
        /// <param name="x">A single point in the distribution range.</param>
        /// 
        /// <returns>
        ///   The logarithm of the probability of <c>x</c>
        ///   occurring in the current distribution.
        /// </returns>
        /// 
        /// <remarks>
        ///   The Probability Density Function (PDF) describes the
        ///   probability that a given value <c>x</c> will occur.
        /// </remarks>
        /// 
        protected internal override double InnerLogProbabilityDensityFunction(double x)
        {
            double lnu = this.d1 * Math.Log(this.d1 * x) + this.d2 * Math.Log(this.d2) -
                (this.d1 + this.d2) * Math.Log(this.d1 * x + this.d2);
            return 0.5 * lnu - Math.Log(x * this.b);
        }


        /// <summary>
        ///   Not available.
        /// </summary>
        /// 
        public override void Fit(double[] observations, double[] weights, IFittingOptions options)
        {
            throw new NotSupportedException();
        }

        /// <summary>
        ///   Creates a new object that is a copy of the current instance.
        /// </summary>
        /// <returns>
        ///   A new object that is a copy of this instance.
        /// </returns>
        /// 
        public override object Clone()
        {
            return new FDistribution(this.d1, this.d2);
        }

        #region ISamplableDistribution<double> Members

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
        public override double[] Generate(int samples, double[] result, Random source)
        {
            return Random(this.d1, this.d2, samples, result, source);
        }

        /// <summary>
        ///   Generates a random observation from the current distribution.
        /// </summary>
        /// 
        /// <returns>A random observations drawn from this distribution.</returns>
        /// 
        public override double Generate(Random source)
        {
            return Random(this.d1, this.d2, source);
        }

        /// <summary>
        ///   Generates a random vector of observations from the 
        ///   F-distribution with the given parameters.
        /// </summary>
        /// 
        /// <param name="d1">The first degree of freedom.</param>
        /// <param name="d2">The second degree of freedom.</param>
        /// <param name="samples">The number of samples to generate.</param>
        ///
        /// <returns>An array of double values sampled from the specified F-distribution.</returns>
        /// 
        public static double[] Random(int d1, int d2, int samples)
        {
            return Random(d1, d2, samples, global::Openize.Accord.Math.Random.Generator.Random);
        }

        /// <summary>
        ///   Generates a random vector of observations from the 
        ///   F-distribution with the given parameters.
        /// </summary>
        /// 
        /// <param name="d1">The first degree of freedom.</param>
        /// <param name="d2">The second degree of freedom.</param>
        /// <param name="samples">The number of samples to generate.</param>
        /// <param name="result">The location where to store the samples.</param>
        ///
        /// <returns>An array of double values sampled from the specified F-distribution.</returns>
        /// 
        public static double[] Random(int d1, int d2, int samples, double[] result)
        {
            return Random(d1, d2, samples, result, global::Openize.Accord.Math.Random.Generator.Random);
        }

        /// <summary>
        ///   Generates a random observation from the 
        ///   F-distribution with the given parameters.
        /// </summary>
        /// 
        /// <param name="d1">The first degree of freedom.</param>
        /// <param name="d2">The second degree of freedom.</param>
        /// 
        /// <returns>A random double value sampled from the specified F-distribution.</returns>
        /// 
        public static double Random(int d1, int d2)
        {
            return Random(d1, d2, global::Openize.Accord.Math.Random.Generator.Random);
        }





        /// <summary>
        ///   Generates a random vector of observations from the 
        ///   F-distribution with the given parameters.
        /// </summary>
        /// 
        /// <param name="d1">The first degree of freedom.</param>
        /// <param name="d2">The second degree of freedom.</param>
        /// <param name="samples">The number of samples to generate.</param>
        /// <param name="source">The random number generator to use as a source of randomness. 
        ///   Default is to use <see cref="Generator.Random"/>.</param>
        ///
        /// <returns>An array of double values sampled from the specified F-distribution.</returns>
        /// 
        public static double[] Random(int d1, int d2, int samples, Random source)
        {
            return Random(d1, d2, samples, new double[samples], source);
        }

        /// <summary>
        ///   Generates a random vector of observations from the 
        ///   F-distribution with the given parameters.
        /// </summary>
        /// 
        /// <param name="d1">The first degree of freedom.</param>
        /// <param name="d2">The second degree of freedom.</param>
        /// <param name="samples">The number of samples to generate.</param>
        /// <param name="result">The location where to store the samples.</param>
        /// <param name="source">The random number generator to use as a source of randomness. 
        ///   Default is to use <see cref="Generator.Random"/>.</param>
        ///
        /// <returns>An array of double values sampled from the specified F-distribution.</returns>
        /// 
        public static double[] Random(int d1, int d2, int samples, double[] result, Random source)
        {
            double[] x = GammaDistribution.Random(shape: d1 / 2.0, scale: 2, samples: samples, result: result, source: source);
            double[] y = GammaDistribution.Random(shape: d2 / 2.0, scale: 2, samples: samples, source: source);

            for (int i = 0; i < x.Length; i++)
                x[i] /= y[i];
            return x;
        }

        /// <summary>
        ///   Generates a random observation from the 
        ///   F-distribution with the given parameters.
        /// </summary>
        /// 
        /// <param name="d1">The first degree of freedom.</param>
        /// <param name="d2">The second degree of freedom.</param>
        /// <param name="source">The random number generator to use as a source of randomness. 
        ///   Default is to use <see cref="Generator.Random"/>.</param>
        /// 
        /// <returns>A random double value sampled from the specified F-distribution.</returns>
        /// 
        public static double Random(int d1, int d2, Random source)
        {
            double x = GammaDistribution.Random(shape: d1 / 2.0, scale: 2, source: source);
            double y = GammaDistribution.Random(shape: d2 / 2.0, scale: 2, source: source);
            return x / y;
        }


        #endregion


        /// <summary>
        ///   Returns a <see cref="System.String"/> that represents this instance.
        /// </summary>
        /// 
        /// <returns>
        ///   A <see cref="System.String"/> that represents this instance.
        /// </returns>
        /// 
        public override string ToString(string format, IFormatProvider formatProvider)
        {
            return String.Format("F(x; df1 = {0}, df2 = {1})",
                this.d1.ToString(format, formatProvider),
                this.d2.ToString(format, formatProvider));
        }
    }
}
