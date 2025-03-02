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

namespace Openize.Accord.Statistics.Distributions.Multivariate.Continuous
{
    using System;
    using Base;
    using Openize.Accord.Math.Core;
    using Openize.Accord.Math.Matrix;
    using Openize.Accord.Core.Exceptions;
    using Openize.Accord.Math.Decompositions;
    using Openize.Accord.Math.Functions;
    using Openize.Accord.Statistics.Distributions.Fitting.Base;

    /// <summary>
    ///   Inverse Wishart Distribution.
    /// </summary>
    ///
    /// <remarks>
    /// <para>
    ///   The inverse Wishart distribution, also called the inverted Wishart distribution,
    ///   is a probability distribution defined on real-valued positive-definite matrices.
    ///   In Bayesian statistics it is used as the conjugate prior for the covariance matrix
    ///   of a multivariate normal distribution.</para>
    ///
    /// <para>
    ///   References:
    ///   <list type="bullet">
    ///     <item><description><a href="http://en.wikipedia.org/wiki/Inverse-Wishart_distribution">
    ///       Wikipedia, The Free Encyclopedia. Inverse Wishart distribution. 
    ///       Available from: http://en.wikipedia.org/wiki/Inverse-Wishart_distribution </a></description></item>
    ///   </list></para>
    /// </remarks>
    /// 
    /// <example>
    /// <code>
    ///   // Create a Inverse Wishart with the parameters
    ///   var invWishart = new InverseWishartDistribution(
    ///   
    ///       // Degrees of freedom
    ///       degreesOfFreedom: 4,
    ///   
    ///       // Scale parameter
    ///       inverseScale: new double[,] 
    ///       {
    ///            {  1.7, -0.2 },
    ///            { -0.2,  5.3 },
    ///       }
    ///   );
    ///   
    ///   // Common measures
    ///   double[] var = invWishart.Variance;  // { -3.4, -10.6 }
    ///   double[,] cov = invWishart.Covariance;  // see below
    ///   double[,] mmean = invWishart.MeanMatrix; // see below
    ///   
    ///   //        cov                mean
    ///   //   -5.78   -4.56        1.7  -0.2 
    ///   //   -4.56  -56.18       -0.2   5.3 
    ///   
    ///   // (the above matrix representations have been transcribed to text using)
    ///   string scov = cov.ToString(DefaultMatrixFormatProvider.InvariantCulture);
    ///   string smean = mmean.ToString(DefaultMatrixFormatProvider.InvariantCulture);
    ///   
    ///   // For compatibility reasons, .Mean stores a flattened mean matrix
    ///   double[] mean = invWishart.Mean; // { 1.7, -0.2, -0.2, 5.3 }
    ///   
    ///   
    ///   // Probability density functions
    ///   double pdf = invWishart.ProbabilityDensityFunction(new double[,] 
    ///   {
    ///       {  5.2,  0.2 }, // 0.000029806281690351203
    ///       {  0.2,  4.2 },
    ///   });
    ///   
    ///   double lpdf = invWishart.LogProbabilityDensityFunction(new double[,] 
    ///   {
    ///       {  5.2,  0.2 }, // -10.420791391688828
    ///       {  0.2,  4.2 },
    ///   });
    /// </code>
    /// </example>
    /// 
    /// <seealso cref="WishartDistribution"/>
    /// 
    [Serializable]
    public class InverseWishartDistribution : MatrixContinuousDistribution
    {
        int size;
        double v; // degrees of freedom
        double[,] inverseScaleMatrix; // Ψ (psi)

        double constant;
        double power;

        double[,] mean;
        double[] variance;
        double[,] covariance;

        /// <summary>
        ///   Creates a new Inverse Wishart distribution.
        /// </summary>
        /// 
        /// <param name="degreesOfFreedom">The degrees of freedom v.</param>
        /// <param name="inverseScale">The inverse scale matrix Ψ (psi).</param>
        /// 
        public InverseWishartDistribution(double degreesOfFreedom, double[,] inverseScale)
            : base(inverseScale.Rows(), inverseScale.Columns())
        {

            if (inverseScale.GetLength(0) != inverseScale.GetLength(1))
                throw new DimensionMismatchException("inverseScale", "Matrix must be square.");

            this.inverseScaleMatrix = inverseScale;
            this.size = inverseScale.GetLength(0);
            this.v = degreesOfFreedom;

            if (degreesOfFreedom <= this.size - 1)
                throw new ArgumentOutOfRangeException("degreesOfFreedom", "Degrees of freedom must be greater "
                + "than or equal to the number of rows in the inverse scale matrix.");

            var chol = new CholeskyDecomposition(inverseScale);

            if (!chol.IsPositiveDefinite)
                throw new NonPositiveDefiniteMatrixException("scale");
            //if (!chol.Symmetric)
            //    throw new NonSymmetricMatrixException("scale");

            double a = Math.Pow(chol.Determinant, this.v / 2.0);
            double b = Math.Pow(2, (this.v * this.size) / 2.0);
            double c = Gamma.Multivariate(this.v / 2.0, this.size);

            this.constant = a / (b * c);
            this.power = -(this.v + this.size + 1) / 2.0;
        }


        /// <summary>
        ///   Gets the mean for this distribution.
        /// </summary>
        /// 
        /// <value>A vector containing the mean values for the distribution.</value>
        /// 
        public override double[,] Mean
        {
            get
            {
                if (this.mean == null)
                    this.mean = this.inverseScaleMatrix.Divide(this.v - this.size - 1);
                return this.mean;
            }
        }

        /// <summary>
        ///   Gets the variance for this distribution.
        /// </summary>
        /// 
        /// <value>A vector containing the variance values for the distribution.</value>
        /// 
        public override double[] Variance
        {
            get
            {
                if (this.variance == null)
                {
                    this.variance = new double[this.size];
                    for (int i = 0; i < this.variance.Length; i++)
                    {
                        double num = 2 * this.inverseScaleMatrix[i, i];
                        double den = (this.v - this.size - 1) * (this.v - this.size - 1) * (this.v - this.size - 3);
                        this.variance[i] = num / den;
                    }
                }
                return this.variance;
            }
        }


        /// <summary>
        ///   Gets the variance-covariance matrix for this distribution.
        /// </summary>
        /// 
        /// <value>A matrix containing the covariance values for the distribution.</value>
        /// 
        public override double[,] Covariance
        {
            get
            {
                if (this.covariance == null)
                {
                    this.covariance = new double[this.size, this.size];
                    for (int i = 0; i < this.size; i++)
                    {
                        double pii = this.inverseScaleMatrix[i, i];

                        for (int j = 0; j < this.size; j++)
                        {
                            double pij = this.inverseScaleMatrix[i, j];
                            double pjj = this.inverseScaleMatrix[j, j];

                            double num = (this.v - this.size + 1) * (pij * pij) + (this.v - this.size - 1) * (pii * pjj);
                            double den = (this.v - this.size) * (this.v - this.size - 1) * (this.v - this.size - 1) * (this.v - this.size - 3);

                            this.covariance[i, j] = num / den;
                        }
                    }
                }

                return this.covariance;
            }
        }

        /// <summary>
        ///   Not supported.
        /// </summary>
        /// 
        protected internal override double InnerDistributionFunction(double[,] x)
        {
            throw new NotSupportedException();
        }

        /// <summary>
        ///   Gets the probability density function (pdf) for
        ///   this distribution evaluated at point <c>x</c>.
        /// </summary>
        /// 
        /// <param name="x">A single point in the distribution range.
        ///   For a matrix distribution, such as the Wishart's, this
        ///   should be a positive-definite matrix or a matrix written
        ///   in flat vector form.
        /// </param>
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
        protected internal override double InnerProbabilityDensityFunction(double[,] x)
        {
            var chol = new CholeskyDecomposition(x);

            double det = chol.Determinant;
            double[,] Vx = chol.Solve(this.inverseScaleMatrix);

            double z = -0.5 * Vx.Trace();
            double a = Math.Pow(det, this.power);
            double b = Math.Exp(z);

            return this.constant * a * b;
        }

        /// <summary>
        ///   Gets the probability density function (pdf) for
        ///   this distribution evaluated at point <c>x</c>.
        /// </summary>
        /// 
        /// <param name="x">A single point in the distribution range.
        ///   For a matrix distribution, such as the Wishart's, this
        ///   should be a positive-definite matrix or a matrix written
        ///   in flat vector form.
        /// </param>
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
        protected internal override double InnerLogProbabilityDensityFunction(double[,] x)
        {
            var chol = new CholeskyDecomposition(x);

            double det = chol.Determinant;
            double[,] Vx = chol.Solve(this.inverseScaleMatrix);

            double z = -0.5 * Vx.Trace();

            return Math.Log(this.constant) + this.power * Math.Log(det) + z;
        }

        /// <summary>
        ///   Not supported.
        /// </summary>
        /// 
        public override void Fit(double[][,] observations, double[] weights, IFittingOptions options)
        {
            throw new NotSupportedException();
        }

        /// <summary>
        ///   Creates a new object that is a copy of the current instance.
        /// </summary>
        /// 
        /// <returns>
        ///   A new object that is a copy of this instance.
        /// </returns>
        /// 
        public override object Clone()
        {
            return new InverseWishartDistribution(this.v, this.inverseScaleMatrix);
        }


        /// <summary>
        ///   Returns a <see cref="System.String" /> that represents this instance.
        /// </summary>
        /// 
        /// <param name="format">The format.</param>
        /// <param name="formatProvider">The format provider.</param>
        /// 
        /// <returns>
        ///   A <see cref="System.String" /> that represents this instance.
        /// </returns>
        /// 
        public override string ToString(string format, IFormatProvider formatProvider)
        {
            return String.Format(formatProvider, "Wishart^-1(X)");
        }
    }
}
