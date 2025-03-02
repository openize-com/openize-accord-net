﻿// Accord Statistics Library
// The Accord.NET Framework
// http://accord-framework.net
//
// Copyright © Fredrik Enqvist, 2016
// fredrikenqvist at hotmail.com
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

namespace Openize.Accord.Statistics.Distributions.Univariate
{
    using System;
    using System.ComponentModel;
    using Openize.Accord.Core.Attributes;
    using Openize.Accord.Core.Ranges;
    using Openize.Accord.Math.Random;
    using Openize.Accord.Statistics.Distributions.Univariate.Base;
    using Openize.Accord.Statistics.Distributions.Univariate.Continuous;

    /// <summary>
    ///   Generalized Pareto distribution (three parameters).
    /// </summary>
    /// 
    /// <remarks>
    /// <para>
    ///   In statistics, the generalized Pareto distribution (GPD) is a family of 
    ///   continuous probability distributions. It is often used to model the tails
    ///   of another distribution. It is specified by three parameters: location μ,
    ///   scale σ, and shape ξ. Sometimes it is specified by only scale and shape
    ///   and sometimes only by its shape parameter. Some references give the shape 
    ///   parameter as κ = − ξ.
    /// </para>
    /// 
    /// <para>
    ///   References:
    ///   <list type="bullet">
    ///     <item><description><a href="https://en.wikipedia.org/wiki/Generalized_Pareto_distribution">
    ///       Wikipedia, The Free Encyclopedia. Generalized Pareto distribution. 
    ///       Available from: https://en.wikipedia.org/wiki/Generalized_Pareto_distribution </a></description></item>
    ///   </list></para>
    /// </remarks>
    /// 
    /// <example>
    /// <code source="Unit Tests\Accord.Tests.Statistics\Distributions\Univariate\Continuous\GeneralizedParetoDistributionTest.cs" region="doc_example1" />
    /// </example>
    /// 
    /// <seealso cref="ParetoDistribution" />
    /// <seealso cref="UnivariateContinuousDistribution" />
    /// 
    [Serializable]
    public class GeneralizedParetoDistribution : UnivariateContinuousDistribution, IFormattable
    {

        double mu;
        double sigma;
        double ksi;


        /// <summary>
        /// Initializes a new instance of the <see cref="GeneralizedParetoDistribution"/> class.
        /// </summary>
        /// 
        public GeneralizedParetoDistribution()
            : this(1, 1, 2)
        {
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="GeneralizedParetoDistribution"/> class.
        /// </summary>
        /// 
        /// <param name="location">The location parameter μ (mu). Default is 0.</param>
        /// <param name="scale">The scale parameter σ (sigma). Must be > 0. Default is 1.</param>
        /// <param name="shape">The shape parameter ξ (Xi). Default is 2.</param>
        /// 
        public GeneralizedParetoDistribution([Real, DefaultValue(0)] double location,
            [Positive] double scale,
            [Real, DefaultValue(0)] double shape)
        {
            if (scale <= 0)
                throw new ArgumentOutOfRangeException("scale", "Scale must be positive.");

            this.init(location, scale, shape);
        }


        private void init(double location, double scale, double shape)
        {
            this.mu = location;
            this.sigma = scale;
            this.ksi = shape;
        }

        /// <summary>
        ///   Gets the scale parameter σ (sigma).
        /// </summary>
        /// 
        public double Scale
        {
            get { return this.sigma; }
        }

        /// <summary>
        ///   Gets shape parameter ξ (Xi).
        /// </summary>
        /// 
        public double Shape
        {
            get { return this.ksi; }
        }

        /// <summary>
        ///   Gets the location parameter μ (mu).
        /// </summary>
        /// 
        public double Location
        {
            get { return this.mu; }
        }

        /// <summary>
        /// Gets the variance for this distribution.
        /// </summary>
        /// <value>The distribution's variance.</value>
        public override double Variance
        {
            get
            {
                double den = ((1.0 - this.ksi) * (1.0 - this.ksi) * 2.0 * (0.5 - this.ksi));
                double num = (this.sigma * this.sigma);
                return num / den;
            }
        }


        /// <summary>
        /// Gets the entropy for this distribution.
        /// </summary>
        /// <value>The distribution's entropy.</value>
        public override double Entropy
        {
            get { return Double.NaN; }
        }


        /// <summary>
        /// Gets the support interval for this distribution.
        /// </summary>
        /// <value>A <see cref="DoubleRange" /> containing
        /// the support interval for this distribution.</value>
        public override DoubleRange Support
        {
            get
            {
                if (this.ksi >= 0)
                    return new DoubleRange(this.mu, Double.PositiveInfinity);
                return new DoubleRange(this.mu, this.mu - this.sigma / this.ksi);
            }
        }


        /// <summary>
        /// Gets the mean for this distribution.
        /// </summary>
        /// <value>The distribution's mean value.</value>
        public override double Mean
        {
            get { return this.mu + (this.sigma / (1 - this.ksi)); }
        }


        /// <summary>
        /// Gets the median for this distribution.
        /// </summary>
        /// <value>The distribution's median value.</value>
        public override double Median
        {
            get
            {
                if (this.ksi == 0)
                    return -Math.Log(0.5) * this.sigma + this.mu;

                double m = Math.Pow(2, this.ksi) - 1;
                return this.mu + (this.sigma * m / this.ksi);
            }
        }

        /// <summary>
        /// Gets the inverse of the cumulative distribution function (icdf) for
        /// this distribution evaluated at probability <c>p</c>. This function
        /// is also known as the Quantile function.
        /// </summary>
        /// <param name="p">A probability value between 0 and 1.</param>
        /// <returns>A sample which could original the given probability
        /// value when applied in the <see cref="UnivariateContinuousDistribution.DistributionFunction(double)" />.</returns>
        /// <remarks>The Inverse Cumulative Distribution Function (ICDF) specifies, for
        /// a given probability, the value which the random variable will be at,
        /// or below, with that probability.</remarks>
        protected internal override double InnerInverseDistributionFunction(double p)
        {
            if (this.ksi == 0)
                return -Math.Log(1 - p) * this.sigma + this.mu;

            double m = Math.Exp(-this.ksi * Math.Log(1 - p)) - 1;
            return m * this.sigma / this.ksi + this.mu;
        }

        /// <summary>
        /// Gets the cumulative distribution function (cdf) for
        /// this distribution evaluated at point <c>x</c>.
        /// </summary>
        /// <param name="x">A single point in the distribution range.</param>
        /// <returns>System.Double.</returns>
        /// <remarks>The Cumulative Distribution Function (CDF) describes the cumulative
        /// probability that a given value or any value smaller than it will occur.</remarks>
        protected internal override double InnerDistributionFunction(double x)
        {
            double m = (x - this.mu) / this.sigma;

            if (this.ksi == 0)
                return 1 - Math.Exp(-m);

            return 1 - Math.Pow(1 + (this.ksi * m), -1 / this.ksi);
        }


        /// <summary>
        /// Gets the probability density function (pdf) for
        /// this distribution evaluated at point <c>x</c>.
        /// </summary>
        /// <param name="x">A single point in the distribution range.</param>
        /// <returns>The probability of <c>x</c> occurring
        /// in the current distribution.</returns>
        /// <remarks>The Probability Density Function (PDF) describes the
        /// probability that a given value <c>x</c> will occur.</remarks>
        protected internal override double InnerProbabilityDensityFunction(double x)
        {
            double m = (x - this.mu) / this.sigma;

            if (this.ksi == 0)
                return Math.Exp(-m) / this.sigma;

            return Math.Pow(1 + (this.ksi * m), -((1 / this.ksi) + 1)) / this.sigma;
        }





        /// <summary>
        /// Creates a new object that is a copy of the current instance.
        /// </summary>
        /// <returns>A new object that is a copy of this instance.</returns>
        public override object Clone()
        {
            return new GeneralizedParetoDistribution(this.mu, this.sigma, this.ksi);
        }


        /// <summary>
        /// Returns a <see cref="System.String" /> that represents this instance.
        /// </summary>
        /// <param name="format">The format.</param>
        /// <param name="formatProvider">The format provider.</param>
        /// <returns>A <see cref="System.String" /> that represents this instance.</returns>
        public override string ToString(string format, IFormatProvider formatProvider)
        {
            return String.Format(formatProvider, "Pareto(x; μ = {0}, σ = {1}, ξ = {2})",
                this.mu.ToString(format, formatProvider),
                this.sigma.ToString(format, formatProvider),
                this.ksi.ToString(format, formatProvider));
        }

        /// <summary>
        /// Generates a random vector of observations from the current distribution.
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
            UniformContinuousDistribution.Random(samples, result, source);

            if (this.ksi == 0)
            {
                for (int i = 0; i < result.Length; i++)
                    result[i] = this.mu - this.sigma * Math.Log(result[i]);
            }
            else
            {
                for (int i = 0; i < result.Length; i++)
                    result[i] = this.mu + this.sigma * (Math.Pow(result[i], -this.ksi) - 1) / this.ksi;
            }

            return result;
        }

        /// <summary>
        /// Generates a random observation from the current distribution.
        /// </summary>
        /// 
        /// <param name="source">The random number generator to use as a source of randomness. 
        ///   Default is to use <see cref="Generator.Random"/>.</param>
        ///   
        /// <returns>A random observations drawn from this distribution.</returns>
        /// 
        public override double Generate(Random source)
        {
            double U = UniformContinuousDistribution.Random(source);

            if (this.ksi == 0)
                return this.mu - this.sigma * Math.Log(U);
            return this.mu + this.sigma * (Math.Pow(U, -this.ksi) - 1) / this.ksi;
        }
    }
}