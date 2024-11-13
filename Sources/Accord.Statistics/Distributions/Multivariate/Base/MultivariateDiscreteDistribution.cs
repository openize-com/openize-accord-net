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

namespace FileFormat.Accord.Statistics.Distributions.Multivariate.Base
{
    using System;
    using FileFormat.Accord.Core.Ranges;
    using Fitting.Base;
    using global::Accord.Math;
    using Math;
    using Math.Matrix;
    using Math.Random;
    using Sampling;

    /// <summary>
    ///   Abstract class for multivariate discrete probability distributions.
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
    [Serializable]
    public abstract class MultivariateDiscreteDistribution : DistributionBase,
        IMultivariateDistribution, IMultivariateDistribution<int[]>,
        IMultivariateDistribution<double[]>,
        ISampleableDistribution<int[]>, ISampleableDistribution<double[]>,
        IDistribution<double[]>, IFormattable
    {

        private int dimension;

        [NonSerialized]
        private MetropolisHasting<int> generator;


        /// <summary>
        ///   Constructs a new MultivariateDiscreteDistribution class.
        /// </summary>
        /// 
        protected MultivariateDiscreteDistribution(int dimension)
        {
            this.dimension = dimension;
        }

        /// <summary>
        ///   Gets the number of variables for this distribution.
        /// </summary>
        /// 
        public int Dimension { get { return this.dimension; } }

        /// <summary>
        ///   Gets the mean for this distribution.
        /// </summary>
        /// 
        /// <value>
        ///   An array of double-precision values containing
        ///   the mean values for this distribution.
        /// </value>
        /// 
        public abstract double[] Mean { get; }

        /// <summary>
        ///   Gets the mean for this distribution.
        /// </summary>
        /// 
        /// <value>
        ///   An array of double-precision values containing
        ///   the variance values for this distribution.
        /// </value>
        /// 
        public abstract double[] Variance { get; }

        /// <summary>
        ///   Gets the variance for this distribution.
        /// </summary>
        /// 
        /// <value>
        ///   An multidimensional array of double-precision values
        ///   containing the covariance values for this distribution.
        /// </value>
        /// 
        public abstract double[,] Covariance { get; }

        /// <summary>
        ///   Gets the mode for this distribution.
        /// </summary>
        /// 
        /// <value>
        ///   An array of double-precision values containing
        ///   the mode values for this distribution.
        /// </value>
        /// 
        public virtual double[] Mode
        {
            get { return this.Mean; }
        }

        /// <summary>
        ///   Gets the median for this distribution.
        /// </summary>
        /// 
        /// <value>
        ///   An array of double-precision values containing
        ///   the median values for this distribution.
        /// </value>
        /// 
        public virtual double[] Median
        {
            get { return this.Mean; }
        }

        /// <summary>
        ///   Gets the support interval for this distribution.
        /// </summary>
        /// 
        /// <value>A <see cref="IntRange"/> containing
        ///  the support interval for this distribution.</value>
        ///  
        public abstract IntRange[] Support { get; }

        #region IDistribution explicit members

        double IDistribution.DistributionFunction(double[] x)
        {
            return this.DistributionFunction(x.Apply(Convert.ToInt32));
        }

        double IDistribution.ProbabilityFunction(double[] x)
        {
            return this.ProbabilityMassFunction(x.Apply(Convert.ToInt32));
        }

        double IDistribution.LogProbabilityFunction(double[] x)
        {
            return this.LogProbabilityMassFunction(x.Apply(Convert.ToInt32));
        }

        double IDistribution.ComplementaryDistributionFunction(double[] x)
        {
            return this.ComplementaryDistributionFunction(x.Apply(Convert.ToInt32));
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
            double[][] multivariate = observations as double[][];
            if (multivariate != null)
            {
                this.Fit(multivariate, weights, options);
                return;
            }

            double[] univariate = observations as double[];
            if (univariate != null)
            {
                this.Fit(univariate.Split(this.dimension), weights, options);
                return;
            }

            throw new ArgumentException("Unsupported parameter type.", "observations");
        }

        void IDistribution.Fit(Array observations, int[] weights, IFittingOptions options)
        {
            double[][] multivariate = observations as double[][];
            if (multivariate != null)
            {
                this.Fit(multivariate, weights, options);
                return;
            }

            double[] univariate = observations as double[];
            if (univariate != null)
            {
                this.Fit(univariate.Split(this.dimension), weights, options);
                return;
            }

            throw new ArgumentException("Unsupported parameter type.", "observations");
        }
        #endregion


        /// <summary>
        ///   Gets the cumulative distribution function (cdf) for
        ///   this distribution evaluated at point <c>x</c>.
        /// </summary>
        /// 
        /// <param name="x">
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
         double DistributionFunction(int[] x)
        {
            if (x == null)
                throw new ArgumentNullException("x");

            return this.InnerDistributionFunction(x);
        }

        /// <summary>
        ///   Gets the cumulative distribution function (cdf) for
        ///   this distribution evaluated at point <c>x</c>.
        /// </summary>
        /// 
        /// <param name="x">
        ///   A single point in the distribution range.</param>
        ///   
        /// <remarks>
        ///   The Cumulative Distribution Function (CDF) describes the cumulative
        ///   probability that a given value or any value smaller than it will occur.
        /// </remarks>
        /// 
#if COMPATIBILITY
        protected internal virtual double InnerDistributionFunction(int[] x)
        {
            throw new NotImplementedException();
        }
#else
        protected internal abstract double InnerDistributionFunction(int[] x);
#endif

        /// <summary>
        ///   Gets the probability mass function (pmf) for
        ///   this distribution evaluated at point <c>x</c>.
        /// </summary>
        /// 
        /// <param name="x">
        ///   A single point in the distribution range.</param>
        ///   
        /// <remarks>
        ///   The Probability Mass Function (PMF) describes the
        ///   probability that a given value <c>x</c> will occur.
        /// </remarks>
        /// 
        /// <returns>
        ///   The probability of <c>x</c> occurring
        ///   in the current distribution.</returns>
        ///   
        public
#if COMPATIBILITY
        virtual
#endif
         double ProbabilityMassFunction(int[] x)
        {
            if (x == null)
                throw new ArgumentNullException("x");

            return this.InnerProbabilityMassFunction(x);
        }

        /// <summary>
        ///   Gets the probability mass function (pmf) for
        ///   this distribution evaluated at point <c>x</c>.
        /// </summary>
        /// 
        /// <param name="x">
        ///   A single point in the distribution range.</param>
        ///   
        /// <remarks>
        ///   The Probability Mass Function (PMF) describes the
        ///   probability that a given value <c>x</c> will occur.
        /// </remarks>
        /// 
        /// <returns>
        ///   The probability of <c>x</c> occurring
        ///   in the current distribution.</returns>
        ///   
#if COMPATIBILITY
        protected internal virtual double InnerProbabilityMassFunction(int[] x)
        {
            throw new NotImplementedException();
        }
#else
        protected internal abstract double InnerProbabilityMassFunction(int[] x);
#endif

        /// <summary>
        ///   Gets the log-probability mass function (pmf) for
        ///   this distribution evaluated at point <c>x</c>.
        /// </summary>
        /// 
        /// <param name="x">
        ///   A single point in the distribution range.</param>
        ///   
        /// <remarks>
        ///   The Probability Mass Function (PMF) describes the
        ///   probability that a given value <c>x</c> will occur.
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
        double LogProbabilityMassFunction(int[] x)
        {
            if (x == null)
                throw new ArgumentNullException("x");

            return this.InnerLogProbabilityMassFunction(x);
        }

        /// <summary>
        ///   Gets the log-probability mass function (pmf) for
        ///   this distribution evaluated at point <c>x</c>.
        /// </summary>
        /// 
        /// <param name="x">
        ///   A single point in the distribution range.</param>
        ///   
        /// <remarks>
        ///   The Probability Mass Function (PMF) describes the
        ///   probability that a given value <c>x</c> will occur.
        /// </remarks>
        /// 
        /// <returns>
        ///   The logarithm of the probability of <c>x</c> 
        ///   occurring in the current distribution.</returns>
        ///   
#if COMPATIBILITY
        protected internal virtual double InnerLogProbabilityMassFunction(int[] x)
        {
            throw new NotImplementedException();
        }
#else
        protected internal virtual double InnerLogProbabilityMassFunction(int[] x)
        {
            return Math.Log(this.ProbabilityMassFunction(x));
        }
#endif

        /// <summary>
        ///   Not supported.
        /// </summary>
        /// 
        public
#if COMPATIBILITY
        virtual
#endif
         int[] InverseDistributionFunction(double p)
        {
            return this.InnerInverseDistributionFunction(p);
        }

        /// <summary>
        ///   Not supported.
        /// </summary>
        protected internal int[] InnerInverseDistributionFunction(double p)
        {
            throw new NotImplementedException();
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
        public
#if COMPATIBILITY
        virtual
#endif
         double ComplementaryDistributionFunction(int[] x)
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
        protected internal double InnerComplementaryDistributionFunction(int[] x)
        {
            return 1.0 - this.DistributionFunction(x);
        }

        /// <summary>
        ///   Gets the marginal distribution of a given variable.
        /// </summary>
        /// 
        /// <param name="index">The variable index.</param>
        /// 
        public double[] MarginalDistributionFunction(int index)
        {
            double[] num = new double[this.Support[index].Length];
            double den = 0;

            int[] symbols = this.Support.Apply(x => x.Length);
            int[] probe = new int[this.Dimension];

            foreach (int[] input in Combinatorics.Sequences(symbols, inPlace: true))
            {
                for (int i = 0; i < probe.Length; i++)
                    probe[i] = input[i] + this.Support[i].Min;

                double p = this.ProbabilityMassFunction(probe);

                for (int k = 0; k < probe.Length; k++)
                {
                    if (k == index)
                        continue;

                    for (int j = 0; j < this.Support[k].Length; j++)
                    {
                        if (probe[k] == this.Support[k].Min + j)
                            num[j] += p;
                    }
                }

                den += p;
            }

            return num.Apply(v => v == 0 ? 0 : v / den);
        }

        /// <summary>
        ///   Gets the marginal distribution of a given variable evaluated at a given value.
        /// </summary>
        /// 
        /// <param name="index">The variable index.</param>
        /// <param name="value">The variable value.</param>
        /// 
        public double MarginalDistributionFunction(int index, int value)
        {
            double num = 0;
            double den = 0;

            int[] symbols = this.Support.Apply(x => x.Length);
            int[] probe = new int[this.Dimension];

            foreach (int[] input in Combinatorics.Sequences(symbols, inPlace: true))
            {
                for (int i = 0; i < probe.Length; i++)
                    probe[i] = input[i] + this.Support[i].Min;

                double p = this.ProbabilityMassFunction(probe);

                for (int k = 0; k < probe.Length; k++)
                {
                    if (k == index)
                        continue;

                    if (probe[k] == this.Support[k].Min + value)
                        num += p;
                }

                den += p;
            }

            return num == 0 ? 0 : num / den;
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
        /// <remarks>
        ///   Although both double[] and double[][] arrays are supported,
        ///   providing a double[] for a multivariate distribution or a
        ///   double[][] for a univariate distribution may have a negative
        ///   impact in performance.
        /// </remarks>
        ///   
        public virtual void Fit(double[][] observations)
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
        /// <remarks>
        ///   Although both double[] and double[][] arrays are supported,
        ///   providing a double[] for a multivariate distribution or a
        ///   double[][] for a univariate distribution may have a negative
        ///   impact in performance.
        /// </remarks>
        /// 
        public virtual void Fit(double[][] observations, double[] weights)
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
        /// <remarks>
        ///   Although both double[] and double[][] arrays are supported,
        ///   providing a double[] for a multivariate distribution or a
        ///   double[][] for a univariate distribution may have a negative
        ///   impact in performance.
        /// </remarks>
        /// 
        public virtual void Fit(double[][] observations, int[] weights)
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
        /// <remarks>
        ///   Although both double[] and double[][] arrays are supported,
        ///   providing a double[] for a multivariate distribution or a
        ///   double[][] for a univariate distribution may have a negative
        ///   impact in performance.
        /// </remarks>
        /// 
        public virtual void Fit(double[][] observations, IFittingOptions options)
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
        /// <param name="weights">
        ///   The weight vector containing the weight for each of the samples.</param>
        /// <param name="options">
        ///   Optional arguments which may be used during fitting, such
        ///   as regularization constants and additional parameters.</param>
        ///   
        /// <remarks>
        ///   Although both double[] and double[][] arrays are supported,
        ///   providing a double[] for a multivariate distribution or a
        ///   double[][] for a univariate distribution may have a negative
        ///   impact in performance.
        /// </remarks>
        /// 
        public virtual void Fit(double[][] observations, double[] weights, IFittingOptions options)
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
        /// <param name="weights">
        ///   The weight vector containing the weight for each of the samples.</param>
        /// <param name="options">
        ///   Optional arguments which may be used during fitting, such
        ///   as regularization constants and additional parameters.</param>
        ///   
        /// <remarks>
        ///   Although both double[] and double[][] arrays are supported,
        ///   providing a double[] for a multivariate distribution or a
        ///   double[][] for a univariate distribution may have a negative
        ///   impact in performance.
        /// </remarks>
        /// 
        public virtual void Fit(double[][] observations, int[] weights, IFittingOptions options)
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
        public int[][] Generate(int samples)
        {
            return this.Generate(samples, global::FileFormat.Accord.Math.Random.Generator.Random);
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
        public double[][] Generate(int samples, double[][] result)
        {
            return this.Generate(samples, result, global::FileFormat.Accord.Math.Random.Generator.Random);
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
        public int[][] Generate(int samples, int[][] result)
        {
            return this.Generate(samples, result, global::FileFormat.Accord.Math.Random.Generator.Random);
        }

        /// <summary>
        ///   Generates a random observation from the current distribution.
        /// </summary>
        /// 
        /// <returns>A random observations drawn from this distribution.</returns>
        /// 
        public int[] Generate()
        {
            return this.Generate(global::FileFormat.Accord.Math.Random.Generator.Random);
        }

        /// <summary>
        ///   Generates a random observation from the current distribution.
        /// </summary>
        /// 
        /// <returns>A random observations drawn from this distribution.</returns>
        /// 
        public double[] Generate(double[] result)
        {
            return this.Generate(result, global::FileFormat.Accord.Math.Random.Generator.Random);
        }

        /// <summary>
        ///   Generates a random observation from the current distribution.
        /// </summary>
        /// 
        /// <returns>A random observations drawn from this distribution.</returns>
        /// 
        public int[] Generate(int[] result)
        {
            return this.Generate(result, global::FileFormat.Accord.Math.Random.Generator.Random);
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
        public int[][] Generate(int samples, Random source)
        {
            return this.Generate(samples, Jagged.Create<int>(samples, this.dimension), source);
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
        public virtual double[][] Generate(int samples, double[][] result, Random source)
        {
            return this.Generate(samples, source).ToDouble(result: result);
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
        public virtual int[][] Generate(int samples, int[][] result, Random source)
        {
            if (this.generator == null)
                this.generator = MetropolisHasting.Discrete(this.Dimension, this);

            if (this.generator.RandomSource != source)
                this.generator.RandomSource = source;

            return this.generator.Generate(samples, result);
        }

        /// <summary>
        ///   Generates a random observation from the current distribution.
        /// </summary>
        /// 
        /// <returns>A random observations drawn from this distribution.</returns>
        /// 
        public int[] Generate(Random source)
        {
            return this.Generate(1, source: source)[0];
        }

        /// <summary>
        ///   Generates a random observation from the current distribution.
        /// </summary>
        /// 
        /// <returns>A random observations drawn from this distribution.</returns>
        /// 
        public double[] Generate(double[] result, Random source)
        {
            return this.Generate(1, result: new[] { result }, source: source)[0];
        }

        /// <summary>
        ///   Generates a random observation from the current distribution.
        /// </summary>
        /// 
        /// <returns>A random observations drawn from this distribution.</returns>
        /// 
        public int[] Generate(int[] result, Random source)
        {
            return this.Generate(1, result: new[] { result }, source: source)[0];
        }



        double IDistribution<int[]>.ProbabilityFunction(int[] x)
        {
            return this.ProbabilityMassFunction(x);
        }

        double IDistribution<int[]>.LogProbabilityFunction(int[] x)
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

        double[] IRandomNumberGenerator<double[]>.Generate()
        {
            return this.Generate(new double[this.dimension]);
        }

        double[][] IRandomNumberGenerator<double[]>.Generate(int samples)
        {
            return this.Generate(samples, Jagged.Create<double>(samples, this.dimension));
        }

        double[][] ISampleableDistribution<double[]>.Generate(int samples, Random source)
        {
            return this.Generate(samples, Jagged.Create<double>(samples, this.dimension), source);
        }

        double[] ISampleableDistribution<double[]>.Generate(Random source)
        {
            return this.Generate(new double[this.dimension], source);
        }

    }
}