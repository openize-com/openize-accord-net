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

namespace Openize.Accord.Statistics.Distributions.Multivariate.Base
{
    using System;
    using Openize.Accord.Math.Matrix;
    using Openize.Accord.Math.Random;
    using Openize.Accord.Statistics.Distributions.Fitting.Base;

    /// <summary>
    ///   Abstract class for Matrix Probability Distributions.
    /// </summary>
    /// 
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
    ///   The function describing the probability that a given value will occur is called
    ///   the probability function (or probability density function, abbreviated PDF), and
    ///   the function describing the cumulative probability that a given value or any value
    ///   smaller than it will occur is called the distribution function (or cumulative
    ///   distribution function, abbreviated CDF).</para>  
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
    [Serializable]
    public abstract class MatrixContinuousDistribution : DistributionBase,
        IMultivariateDistribution,
        IMultivariateDistribution<double[,]>,
        IFittableDistribution<double[,]>,
        ISampleableDistribution<double[,]>,
        IMultivariateDistribution<double[]>,
        IFittableDistribution<double[]>,
        ISampleableDistribution<double[]>,
        IFormattable
    {

        private int rows;
        private int cols;


        /// <summary>
        ///   Constructs a new MultivariateDistribution class.
        /// </summary>
        /// 
        /// <param name="rows">The number of rows for matrices modeled by the distribution.</param>
        /// <param name="cols">The number of rows for matrices modeled by the distribution.</param>
        /// 
        protected MatrixContinuousDistribution(int rows, int cols)
        {
            this.rows = rows;
            this.cols = cols;
        }

        /// <summary>
        ///   Gets the number of variables for this distribution.
        /// </summary>
        /// 
        public int Dimension { get { return this.rows * this.cols; } }

        /// <summary>
        ///   Gets the number of rows that matrices from this distribution should have.
        /// </summary>
        /// 
        public int NumberOfRows { get { return this.rows; } }

        /// <summary>
        ///   Gets the number of columns that matrices from this distribution should have.
        /// </summary>
        /// 
        public int NumberOfColumns { get { return this.cols; } }

        /// <summary>
        ///   Gets the mean for this distribution.
        /// </summary>
        /// 
        /// <value>A vector containing the mean values for the distribution.</value>
        /// 
        public abstract double[,] Mean { get; }

        /// <summary>
        ///   Gets the variance for this distribution.
        /// </summary>
        /// 
        /// <value>A vector containing the variance values for the distribution.</value>
        /// 
        public abstract double[] Variance { get; }

        /// <summary>
        ///   Gets the variance-covariance matrix for this distribution.
        /// </summary>
        /// 
        /// <value>A matrix containing the covariance values for the distribution.</value>
        /// 
        public abstract double[,] Covariance { get; }

        /// <summary>
        ///   Gets the mode for this distribution.
        /// </summary>
        /// 
        /// <value>A vector containing the mode values for the distribution.</value>
        /// 
        public virtual double[,] Mode
        {
            get { return this.Mean; }
        }

        /// <summary>
        ///   Gets the median for this distribution.
        /// </summary>
        /// 
        /// <value>A vector containing the median values for the distribution.</value>
        /// 
        public virtual double[,] Median
        {
            get { return this.Mean; }
        }

        double[] IMultivariateDistribution.Mean
        {
            get { return this.Mean.Reshape(); }
        }

        double[] IMultivariateDistribution.Median
        {
            get { return this.Median.Reshape(); }
        }

        double[] IMultivariateDistribution.Mode
        {
            get { return this.Mode.Reshape(); }
        }

        double[] IMultivariateDistribution<double[,]>.Mean
        {
            get { return this.Mean.Reshape(); }
        }

        double[] IMultivariateDistribution<double[,]>.Median
        {
            get { return this.Median.Reshape(); }
        }

        double[] IMultivariateDistribution<double[,]>.Mode
        {
            get { return this.Mode.Reshape(); }
        }

        double[] IMultivariateDistribution<double[]>.Mean
        {
            get { return this.Mean.Reshape(); }
        }

        double[] IMultivariateDistribution<double[]>.Median
        {
            get { return this.Median.Reshape(); }
        }

        double[] IMultivariateDistribution<double[]>.Mode
        {
            get { return this.Mode.Reshape(); }
        }



        #region IDistribution explicit members

        double IDistribution.ProbabilityFunction(double[] x)
        {
            return this.ProbabilityDensityFunction(this.reshape(x));
        }

        double IDistribution.LogProbabilityFunction(double[] x)
        {
            return this.LogProbabilityDensityFunction(this.reshape(x));
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
            double[][,] multivariate = observations as double[][,];
            if (multivariate != null)
            {
                this.Fit(multivariate, weights, options);
                return;
            }

            throw new ArgumentException("Unsupported parameter type.", "observations");
        }

        void IDistribution.Fit(Array observations, int[] weights, IFittingOptions options)
        {
            double[][,] multivariate = observations as double[][,];
            if (multivariate != null)
            {
                this.Fit(multivariate, weights, options);
                return;
            }

            throw new ArgumentException("Unsupported parameter type.", "observations");
        }
        #endregion


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
        ///   
        /// <remarks>
        ///   The Probability Density Function (PDF) describes the
        ///   probability that a given value <c>x</c> will occur.
        /// </remarks>
        /// 
        /// <returns>
        ///   The probability of <c>x</c> occurring
        ///   in the current distribution.</returns>
        ///   
        public double DistributionFunction(double[,] x)
        {
            if (x == null)
                throw new ArgumentNullException("x");

            return this.InnerDistributionFunction(x);
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
        ///   
        /// <remarks>
        ///   The Probability Density Function (PDF) describes the
        ///   probability that a given value <c>x</c> will occur.
        /// </remarks>
        /// 
        /// <returns>
        ///   The probability of <c>x</c> occurring
        ///   in the current distribution.</returns>
        ///   
        protected internal abstract double InnerDistributionFunction(double[,] x);

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
        ///   
        /// <remarks>
        ///   The Probability Density Function (PDF) describes the
        ///   probability that a given value <c>x</c> will occur.
        /// </remarks>
        /// 
        /// <returns>
        ///   The probability of <c>x</c> occurring
        ///   in the current distribution.</returns>
        ///   
        public double ProbabilityDensityFunction(double[,] x)
        {
            if (x == null)
                throw new ArgumentNullException("x");

            return this.InnerProbabilityDensityFunction(x);
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
        ///   
        /// <remarks>
        ///   The Probability Density Function (PDF) describes the
        ///   probability that a given value <c>x</c> will occur.
        /// </remarks>
        /// 
        /// <returns>
        ///   The probability of <c>x</c> occurring
        ///   in the current distribution.</returns>
        ///   
        protected internal abstract double InnerProbabilityDensityFunction(double[,] x);

        /// <summary>
        ///   Gets the log-probability density function (pdf)
        ///   for this distribution evaluated at point <c>x</c>.
        /// </summary>
        /// 
        /// <param name="x">
        ///   A single point in the distribution range. For a 
        ///   univariate distribution, this should be a single
        ///   double value. For a multivariate distribution,
        ///   this should be a double array.</param>
        ///   
        /// <returns>
        ///   The logarithm of the probability of <c>x</c>
        ///   occurring in the current distribution.</returns>
        ///   
        public double LogProbabilityDensityFunction(double[,] x)
        {
            if (x == null)
                throw new ArgumentNullException("x");

            return this.InnerLogProbabilityDensityFunction(x);
        }

        /// <summary>
        ///   Gets the log-probability density function (pdf)
        ///   for this distribution evaluated at point <c>x</c>.
        /// </summary>
        /// 
        /// <param name="x">
        ///   A single point in the distribution range. For a 
        ///   univariate distribution, this should be a single
        ///   double value. For a multivariate distribution,
        ///   this should be a double array.</param>
        ///   
        /// <returns>
        ///   The logarithm of the probability of <c>x</c>
        ///   occurring in the current distribution.</returns>
        ///   
        protected internal virtual double InnerLogProbabilityDensityFunction(double[,] x)
        {
            return Math.Log(this.ProbabilityDensityFunction(x));
        }

        /// <summary>
        ///   Gets the complementary cumulative distribution function
        ///   (ccdf) for this distribution evaluated at point <c>x</c>.
        ///   This function is also known as the Survival function.
        /// </summary>
        /// 
        /// <remarks>
        ///   The Complementary Cumulative Distribution Function (CCDF) is
        ///   the complement of the Cumulative Distribution Function, or 1
        ///   minus the CDF.
        /// </remarks>
        /// 
        public double ComplementaryDistributionFunction(double[,] x)
        {
            if (x == null)
                throw new ArgumentNullException("x");

            return this.InnerComplementaryDistributionFunction(x);
        }

        /// <summary>
        ///   Gets the complementary cumulative distribution function
        ///   (ccdf) for this distribution evaluated at point <c>x</c>.
        ///   This function is also known as the Survival function.
        /// </summary>
        /// 
        /// <remarks>
        ///   The Complementary Cumulative Distribution Function (CCDF) is
        ///   the complement of the Cumulative Distribution Function, or 1
        ///   minus the CDF.
        /// </remarks>
        /// 
        protected internal virtual double InnerComplementaryDistributionFunction(double[,] x)
        {
            return 1.0 - this.DistributionFunction(x);
        }

        /// 
        /// <summary>
        ///   Fits the underlying distribution to a given set of observations.
        /// </summary>
        /// 
        /// <param name="observations">
        ///   The array of observations to fit the model against. The array
        ///   elements can be either of type double (for univariate data) or
        ///   type double[] (for multivariate data).</param>
        ///   
        public virtual void Fit(double[][,] observations)
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
        /// <param name="weights">
        ///   The weight vector containing the weight for each of the samples.</param>
        ///   
        public virtual void Fit(double[][,] observations, double[] weights)
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
        /// <param name="weights">
        ///   The weight vector containing the weight for each of the samples.</param>
        ///   
        public virtual void Fit(double[][,] observations, int[] weights)
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
        /// <param name="options">
        ///   Optional arguments which may be used during fitting, such
        ///   as regularization constants and additional parameters.</param>
        ///   
        public virtual void Fit(double[][,] observations, IFittingOptions options)
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
        ///   type double[] (for multivariate data). </param>
        /// <param name="weights">
        ///   The weight vector containing the weight for each of the samples.</param>
        /// <param name="options">
        ///   Optional arguments which may be used during fitting, such
        ///   as regularization constants and additional parameters.</param>
        ///   
        public virtual void Fit(double[][,] observations, double[] weights, IFittingOptions options)
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
        ///   type double[] (for multivariate data). </param>
        /// <param name="weights">
        ///   The weight vector containing the weight for each of the samples.</param>
        /// <param name="options">
        ///   Optional arguments which may be used during fitting, such
        ///   as regularization constants and additional parameters.</param>
        ///   
        public virtual void Fit(double[][,] observations, int[] weights, IFittingOptions options)
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
        public double[][,] Generate(int samples)
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
        public double[][,] Generate(int samples, double[][,] result)
        {
            return this.Generate(samples, result, global::Openize.Accord.Math.Random.Generator.Random);
        }

        /// <summary>
        ///   Generates a random observation from the current distribution.
        /// </summary>
        /// 
        /// <returns>A random observations drawn from this distribution.</returns>
        /// 
        public double[,] Generate()
        {
            return this.Generate(global::Openize.Accord.Math.Random.Generator.Random);
        }

        /// <summary>
        ///   Generates a random observation from the current distribution.
        /// </summary>
        /// 
        /// <returns>A random observations drawn from this distribution.</returns>
        /// 
        public double[,] Generate(double[,] result)
        {
            return this.Generate(result, global::Openize.Accord.Math.Random.Generator.Random);
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
        public double[][,] Generate(int samples, Random source)
        {
            return this.Generate(samples, new double[samples].Apply(x => new double[this.rows, this.cols]), source);
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
        public virtual double[][,] Generate(int samples, double[][,] result, Random source)
        {
            throw new NotSupportedException();
        }

        /// <summary>
        ///   Generates a random observation from the current distribution.
        /// </summary>
        /// 
        /// <returns>A random observations drawn from this distribution.</returns>
        /// 
        public double[,] Generate(Random source)
        {
            return this.Generate(1, source)[0];
        }


        /// <summary>
        /// Generates a random observation from the current distribution.
        /// </summary>
        /// <param name="result">The location where to store the sample.</param>
        /// <param name="source">The random number generator to use as a source of randomness.
        /// Default is to use <see cref="P:Accord.Math.Random.Generator.Random" />.</param>
        /// <returns>A random observation drawn from this distribution.</returns>
        public double[,] Generate(double[,] result, Random source)
        {
            return this.Generate(1, result: new[] { result }, source: source)[0];
        }



        /// <summary>
        /// Generates a random observation from the current distribution.
        /// </summary>
        /// <param name="result">The location where to store the sample.</param>
        /// <returns>A random observation drawn from this distribution.</returns>
        public virtual double[] Generate(double[] result)
        {
            return this.Generate(result, global::Openize.Accord.Math.Random.Generator.Random);
        }


        /// <summary>
        /// Generates a random observation from the current distribution.
        /// </summary>
        /// <param name="result">The location where to store the sample.</param>
        /// <param name="source">The random number generator to use as a source of randomness.
        /// Default is to use <see cref="P:Accord.Math.Random.Generator.Random" />.</param>
        /// <returns>A random observation drawn from this distribution.</returns>
        public double[] Generate(double[] result, Random source)
        {
            return this.Generate(1, result: new[] { result }, source: source)[0];
        }

        double[][] ISampleableDistribution<double[]>.Generate(int samples, Random source)
        {
            return this.Generate(samples, Jagged.Create<double>(samples, this.Dimension));
        }

        /// <summary>
        /// Generates a random vector of observations from the current distribution.
        /// </summary>
        /// <param name="samples">The number of samples to generate.</param>
        /// <param name="result">The location where to store the samples.</param>
        /// <param name="source">The random number generator to use as a source of randomness.
        /// Default is to use <see cref="P:Accord.Math.Random.Generator.Random" />.</param>
        /// <returns>A random vector of observations drawn from this distribution.</returns>
        public double[][] Generate(int samples, double[][] result, Random source)
        {
            var buffer = new double[this.rows, this.cols];
            for (int i = 0; i < samples; i++)
            {
                this.Generate(buffer, source);
                Buffer.BlockCopy(buffer, 0, result, 0, this.Dimension);
            }

            return result;
        }

        double[] ISampleableDistribution<double[]>.Generate(Random source)
        {
            return this.Generate(source).Reshape();
        }

        double[][] IRandomNumberGenerator<double[]>.Generate(int samples)
        {
            return this.Generate(samples, Jagged.Create<double>(samples, this.Dimension));
        }

        /// <summary>
        /// Generates a random vector of observations from the current distribution.
        /// </summary>
        /// <param name="samples">The number of samples to generate.</param>
        /// <param name="result">The location where to store the samples.</param>
        /// <returns>A random vector of observations drawn from this distribution.</returns>
        public double[][] Generate(int samples, double[][] result)
        {
            return this.Generate(samples, result, global::Openize.Accord.Math.Random.Generator.Random);
        }

        double[] IRandomNumberGenerator<double[]>.Generate()
        {
            return this.Generate().Reshape();
        }




        double IDistribution<double[,]>.ProbabilityFunction(double[,] x)
        {
            return this.ProbabilityDensityFunction(x);
        }

        double IDistribution<double[,]>.LogProbabilityFunction(double[,] x)
        {
            return this.LogProbabilityDensityFunction(x);
        }

        double IDistribution<double[]>.ProbabilityFunction(double[] x)
        {
            return this.ProbabilityDensityFunction(this.reshape(x));
        }

        double IDistribution<double[]>.LogProbabilityFunction(double[] x)
        {
            return this.LogProbabilityDensityFunction(this.reshape(x));
        }


        private double[,] reshape(double[] x)
        {
            return x.Reshape(this.rows, this.cols);
        }

        /// <summary>
        /// Gets the cumulative distribution function (cdf) for
        /// this distribution evaluated at point <c>x</c>.
        /// </summary>
        /// <param name="x">The x.</param>
        /// <returns>System.Double.</returns>
        /// <remarks>The Cumulative Distribution Function (CDF) describes the cumulative
        /// probability that a given value or any value smaller than it will occur.</remarks>
        public double DistributionFunction(params double[] x)
        {
            return this.DistributionFunction(this.reshape(x));
        }

        /// <summary>
        /// Gets the complementary cumulative distribution function
        /// (ccdf) for this distribution evaluated at point <c>x</c>.
        /// This function is also known as the Survival function.
        /// </summary>
        /// <param name="x">The x.</param>
        /// <returns>System.Double.</returns>
        /// <remarks>The Complementary Cumulative Distribution Function (CCDF) is
        /// the complement of the Cumulative Distribution Function, or 1
        /// minus the CDF.</remarks>
        public double ComplementaryDistributionFunction(params double[] x)
        {
            return this.ComplementaryDistributionFunction(this.reshape(x));
        }

        /// <summary>
        ///   Fits the underlying distribution to a given set of observations.
        /// </summary>
        /// 
        /// <param name="observations">
        ///   The array of observations to fit the model against. The array
        ///   elements can be either of type double (for univariate data) or
        ///   type double[] (for multivariate data). </param>
        ///   
        public void Fit(double[][] observations)
        {
            this.Fit(observations, null);
        }

        /// <summary>
        ///   Fits the underlying distribution to a given set of observations.
        /// </summary>
        /// 
        /// <param name="observations">
        ///   The array of observations to fit the model against. The array
        ///   elements can be either of type double (for univariate data) or
        ///   type double[] (for multivariate data). </param>
        /// <param name="weights">
        ///   The weight vector containing the weight for each of the samples.</param>
        ///   
        public virtual void Fit(double[][] observations, double[] weights)
        {
            this.Fit(observations.Apply(x => x.Reshape(this.rows, this.cols)), weights);
        }

    }
}