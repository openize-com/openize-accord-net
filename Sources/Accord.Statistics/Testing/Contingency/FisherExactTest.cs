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

namespace Openize.Accord.Statistics.Testing.Contingency
{
    using System;
    using System.Diagnostics;
    using Base;
    using Openize.Accord.Statistics.Analysis.Performance;
    using Openize.Accord.Statistics.Distributions.Univariate.Discrete;

    /// <summary>
    ///   Fisher's exact test for contingency tables.
    /// </summary>
    /// 
    /// <remarks>
    /// <para>
    ///   This test statistic distribution is the <see cref="HypergeometricDistribution">
    ///   Hypergeometric</see>.</para>
    /// </remarks>
    /// 
    /// <seealso cref="HypergeometricDistribution"/>
    /// 
    /// <example>
    /// <code source="Unit Tests\Accord.Tests.Statistics\Testing\Contingency\FisherExactTestTest.cs" region="doc_test" />
    /// </example>
    /// 
    [Serializable]
    public class FisherExactTest : HypothesisTest<HypergeometricDistribution>
    {

        /// <summary>
        ///   Gets the alternative hypothesis under test. If the test is
        ///   <see cref="IHypothesisTest.Significant"/>, the null hypothesis
        ///   can be rejected in favor of this alternative hypothesis.
        /// </summary>
        /// 
        public OneSampleHypothesis Hypothesis { get; protected set; }

        /// <summary>
        ///   Constructs a new Fisher's exact test.
        /// </summary>
        /// 
        /// <param name="matrix">The matrix to be tested.</param>
        /// <param name="alternate">The alternative hypothesis (research hypothesis) to test.</param>
        /// 
        public FisherExactTest(ConfusionMatrix matrix,
            OneSampleHypothesis alternate = OneSampleHypothesis.ValueIsDifferentFromHypothesis)
        {
            if (matrix == null)
                throw new ArgumentNullException("matrix");

            int[,] mat = matrix.Matrix;

            int a = mat[0, 0];
            int b = mat[0, 1];
            int c = mat[1, 0];
            int d = mat[1, 1];

            int N = matrix.NumberOfSamples; // total number in population
            int k = a;                      // number of success in sample
            int m = a + b;                  // number of success in population
            int n = a + c;                  // sample size

            Debug.Assert(a + b + c + d == N);

            this.Compute(statistic: k, populationSize: N, n: n, m: m, alternate: alternate);
        }

        /// <summary>
        ///   Computes the Fisher's exact test.
        /// </summary>
        /// 
        protected void Compute(double statistic, int populationSize, int n, int m, OneSampleHypothesis alternate)
        {
            this.Statistic = statistic;
            this.StatisticDistribution = new HypergeometricDistribution(populationSize, m, n);

            this.Hypothesis = alternate;
            this.Tail = (DistributionTail)alternate;
            this.PValue = this.StatisticToPValue(this.Statistic);
        }

        /// <summary>
        ///   Converts a given test statistic to a p-value.
        /// </summary>
        /// 
        /// <param name="x">The value of the test statistic.</param>
        /// 
        /// <returns>The p-value for the given statistic.</returns>
        /// 
        public override double StatisticToPValue(double x)
        {
            double p;
            switch (this.Tail)
            {
                case DistributionTail.TwoTail:
                    p = 2.0 * this.StatisticDistribution.DistributionFunction((int)x);
                    break;

                case DistributionTail.OneUpper:
                    p = this.StatisticDistribution.ComplementaryDistributionFunction((int)x, inclusive: true);
                    break;

                case DistributionTail.OneLower:
                    p = this.StatisticDistribution.DistributionFunction((int)x);
                    break;

                default:
                    throw new InvalidOperationException();
            }
            return p;
        }

        /// <summary>
        ///   Converts a given p-value to a test statistic.
        /// </summary>
        /// 
        /// <param name="p">The p-value.</param>
        /// 
        /// <returns>The test statistic which would generate the given p-value.</returns>
        /// 
        public override double PValueToStatistic(double p)
        {
            throw new System.NotImplementedException();
        }
    }
}
