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

namespace Openize.Accord.Statistics.Testing.Contingency
{
    using System;
    using Openize.Accord.Math.Matrix;
    using Openize.Accord.Statistics.Analysis.Performance;

    /// <summary>
    ///   Bhapkar test of homogeneity for contingency tables.
    /// </summary>
    /// 
    /// <remarks>
    ///   The Bhapkar test is a more powerful alternative to the
    ///   <see cref="StuartMaxwellTest">Stuart-Maxwell test</see>.
    ///   
    /// <para>
    ///   This is a <see cref="ChiSquareTest">Chi-square kind of test</see>.</para>
    ///   
    /// <para>
    ///   References:
    ///   <list type="bullet">
    ///     <item><description>
    ///      Bhapkar, V.P. (1966). A note on the equivalence of two test criteria
    ///      for hypotheses in categorical data. Journal of the American Statistical
    ///      Association, 61, 228-235.</description></item>
    ///    </list></para>
    /// </remarks>
    /// 
    /// <seealso cref="ChiSquareTest"/>
    ///
    [Serializable]
    public class BhapkarTest : ChiSquareTest
    {

        double[] d;
        double[,] S;
        double[,] invS;

        /// <summary>
        ///   Gets the delta vector <c>d</c> used
        ///   in the test calculations.
        /// </summary>
        /// 
        public double[] Delta
        {
            get { return this.d; }
        }

        /// <summary>
        ///   Gets the covariance matrix <c>S</c>
        ///   used in the test calculations.
        /// </summary>
        /// 
        public double[,] Covariance
        {
            get { return this.S; }
        }

        /// <summary>
        ///   Gets the inverse covariance matrix
        ///   <c>S^-1</c> used in the calculations.
        /// </summary>
        /// 
        public double[,] Precision
        {
            get { return this.invS; }
        }


        /// <summary>
        ///   Creates a new Bhapkar test.
        /// </summary>
        /// 
        /// <param name="matrix">The contingency table to test.</param>
        /// 
        public BhapkarTest(GeneralConfusionMatrix matrix)
        {
            int classes = matrix.NumberOfClasses;
            int samples = matrix.NumberOfSamples;

            int df = classes - 1;

            int[] rowMarginals = matrix.RowTotals;
            int[] colMarginals = matrix.ColumnTotals;

            this.d = new double[df];
            for (int i = 0; i < this.d.Length; i++)
                this.d[i] = rowMarginals[i] - colMarginals[i];

            this.S = new double[df, df];

            for (int i = 0; i < df; i++)
            {
                for (int j = 0; j < df; j++)
                {
                    if (i == j)
                    {
                        double u = (rowMarginals[i] - colMarginals[i]);
                        double pii = matrix.Matrix[i, i];

                        this.S[i, i] = rowMarginals[i] + colMarginals[i] - 2.0 * pii - u * u / (double)samples;
                    }
                    else
                    {
                        double pij = matrix.Matrix[i, j];
                        double pji = matrix.Matrix[j, i];

                        this.S[i, j] = -(pij + pji) - (rowMarginals[i] - colMarginals[i]) * (rowMarginals[j] - colMarginals[j]) / (double)samples;
                    }
                }
            }

            this.invS = this.S.PseudoInverse();

            double chiSquare = this.d.DotAndDot(this.invS, this.d);
            
            this.Compute(chiSquare, df);
        }

      

    }
}
