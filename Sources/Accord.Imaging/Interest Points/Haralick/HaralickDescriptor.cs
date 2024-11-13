// Accord Imaging Library
// The Accord.NET Framework
// http://accord-framework.net
//
// Copyright © Diego Catalano, 2013
// diego.catalano at live.com
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

namespace FileFormat.Accord.Imaging.Interest_Points.Haralick
{
    using System;
    using Core.Exceptions;
    using global::Accord;
    using global::Accord.Math;
    using Math;
    using Math.Accord.Statistics;
    using Math.Decompositions;
    using Math.Matrix;

    /// <summary>
    ///   Haralick's Texture Features.
    /// </summary>
    /// 
    /// <remarks>
    /// <para>
    ///   Haralick's texture features are based on measures derived from
    ///   <see cref="GrayLevelCooccurrenceMatrix">Gray-level Co-occurrence 
    ///   matrices (GLCM)</see>.</para>
    /// <para>
    ///   Whether considering the intensity or grayscale values of the image 
    ///   or various dimensions of color, the co-occurrence matrix can measure
    ///   the texture of the image. Because co-occurrence matrices are typically
    ///   large and sparse, various metrics of the matrix are often taken to get
    ///   a more useful set of features. Features generated using this technique
    ///   are usually called Haralick features, after R. M. Haralick, attributed to
    ///   his paper Textural features for image classification (1973).</para>
    ///   
    /// <para>
    ///   This class encompasses most of the features derived on Haralick's original
    ///   paper. All features are lazy-evaluated until needed; but may also be
    ///   combined in a single feature vector by calling <see cref="GetVector(int)"/>.</para>
    ///   
    /// <para>
    ///   References:
    ///   <list type="bullet">
    ///     <item><description>
    ///       Wikipedia Contributors, "Co-occurrence matrix". Available at
    ///       http://en.wikipedia.org/wiki/Co-occurrence_matrix </description></item>
    ///     <item><description>
    ///       Robert M Haralick, K Shanmugam, Its'hak Dinstein; "Textural 
    ///       Features for Image Classification". IEEE Transactions on Systems, Man,
    ///       and Cybernetics. SMC-3 (6): 610–621, 1973. Available at:
    ///       <a href="http://www.makseq.com/materials/lib/Articles-Books/Filters/Texture/Co-occurrence/haralick73.pdf">
    ///       http://www.makseq.com/materials/lib/Articles-Books/Filters/Texture/Co-occurrence/haralick73.pdf </a>
    ///       </description></item>
    ///   </list></para>
    /// </remarks>
    /// 
    /// <example>
    /// <para>
    ///   For a complete example on how to use <see cref="Haralick"/>, please refer to
    ///   the <see cref="Haralick">documentation of the main Haralick class</see>.</para>
    /// </example>
    /// 
    /// <seealso cref="GrayLevelCooccurrenceMatrix"/>
    /// <seealso cref="Haralick"/>
    /// 
    [Serializable]
    public class HaralickDescriptor
    {
        double[,] matrix;


        int size; // Ng, number of gray levels

        double[] px; // row marginal
        double[] py; // col marginal

        double? sum;
        double? mean;  // μ (mu), mean of matrix
        double? xmean; // μ_x (mu x), mean of px
        double? ymean; // μ_y (mu y), mean of py
        double? xdev;  // σ_x (sigma y), standard deviation of py
        double? ydev;  // σ_y (sigma x), standard deviation of px
        double? xentropy;    // H_x, entropy of px
        double? yentropy;    // H_y, entropy of py

        double[] xysum;  // p_x+y
        double[] xydiff; // p_x-y

        double? angular;       // f1: energy / angular second moment
        double? contrast;      // f2: contrast
        double? correlation;   // f3: correlation
        double? variance;      // f4: sum of squares: variance
        double? inverse;       // f5: inverse difference moment
        double? sumAverage;    // f6: sum average
        double? sumVariance;   // f7: sum variance
        double? sumEntropy;    // f8: sum entropy
        double? entropy;       // f9: matrix entropy
        double? diffVariance;  // f10: difference variance
        double? diffEntropy;   // f11: difference entropy
        double? information1;  // f12: first information measure of correlation
        double? information2;  // f13: second information measure of correlation
        double? maximal;       // f14: maximal correlation coefficient

        // bonus features
        double? laplace;  // contrast using absolute value instead of square
        double? textureHomogeneity;
        double? clusterShade;
        double? clusterProminence;

        const double epsilon = Constants.DoubleSmall;


        /// <summary>
        ///   Initializes a new instance of the <see cref="HaralickDescriptor"/> class.
        /// </summary>
        /// 
        /// <param name="cooccurrenceMatrix">The co-occurrence matrix to compute features from.</param>
        /// 
        public HaralickDescriptor(double[,] cooccurrenceMatrix)
        {
            if (cooccurrenceMatrix == null)
                throw new ArgumentNullException("cooccurrenceMatrix");
            if (cooccurrenceMatrix.Rows() != cooccurrenceMatrix.Columns())
                throw new DimensionMismatchException("cooccurrenceMatrix", "Matrix must be square");

            this.matrix = cooccurrenceMatrix;
            this.size = cooccurrenceMatrix.GetLength(0);
        }


        #region Common calculations

        /// <summary>
        ///   Gets the number of gray levels in the 
        ///   original image. This is the number of
        ///   dimensions of the co-occurrence matrix.
        /// </summary>
        /// 
        public int GrayLevels { get { return this.size; } }

        /// <summary>
        ///   Gets the matrix sum.
        /// </summary>
        /// 
        public double Sum
        {
            get
            {
                if (this.sum == null)
                {
                    double s = 0;
                    foreach (double v in this.matrix)
                        s += v;
                    this.sum = s;
                }
                return this.sum.Value;
            }
        }

        /// <summary>
        ///   Gets the matrix mean μ.
        /// </summary>
        /// 
        public double Mean
        {
            get
            {
                if (this.mean == null)
                    this.mean = this.Sum / this.matrix.Length;
                return this.mean.Value;
            }
        }

        /// <summary>
        ///   Gets the marginal probability vector
        ///   obtained by summing the rows of p(i,j),
        ///   given as p<sub>x</sub>(i) = Σ<sub>j</sub> p(i,j).
        /// </summary>
        /// 
        public double[] RowMarginal
        {
            get
            {
                if (this.px == null)
                {
                    this.px = new double[this.size];
                    for (int i = 0; i < this.px.Length; i++)
                        for (int j = 0; j < this.size; j++)
                            this.px[i] += this.matrix[i, j];
                }
                return this.px;
            }
        }

        /// <summary>
        ///   Gets the marginal probability vector
        ///   obtained by summing the columns of p(i,j),
        ///   given as p<sub>y</sub>(j) = Σ<sub>i</sub> p(i,j).
        /// </summary>
        /// 
        public double[] ColumnMarginal
        {
            get
            {
                if (this.py == null)
                {
                    this.py = new double[this.size];
                    for (int j = 0; j < this.py.Length; j++)
                        for (int i = 0; i < this.size; i++)
                            this.py[j] += this.matrix[i, j];
                }
                return this.py;
            }
        }

        /// <summary>
        ///   Gets μ<sub>x</sub>, the mean value of the 
        ///   <see cref="RowMarginal"/> vector.
        /// </summary>
        /// 
        public double RowMean
        {
            get
            {
                if (this.xmean == null)
                    this.xmean = this.RowMarginal.Mean();
                return this.xmean.Value;
            }
        }

        /// <summary>
        ///   Gets μ_y, the mean value of the 
        ///   <see cref="ColumnMarginal"/> vector.
        /// </summary>
        /// 
        public double ColumnMean
        {
            get
            {
                if (this.ymean == null)
                    this.ymean = this.ColumnMarginal.Mean();
                return this.ymean.Value;
            }
        }

        /// <summary>
        ///   Gets σ<sub>x</sub>, the variance of the 
        ///   <see cref="RowMarginal"/> vector.
        /// </summary>
        /// 
        public double RowStandardDeviation
        {
            get
            {
                if (this.xdev == null)
                    this.xdev = this.RowMarginal.StandardDeviation(this.RowMean);
                return this.xdev.Value;
            }
        }

        /// <summary>
        ///   Gets σ<sub>y</sub>, the variance of the 
        ///   <see cref="ColumnMarginal"/> vector.
        /// </summary>
        /// 
        public double ColumnStandardDeviation
        {
            get
            {
                if (this.ydev == null)
                    this.ydev = this.ColumnMarginal.StandardDeviation(this.ColumnMean);
                return this.ydev.Value;
            }
        }

        /// <summary>
        ///   Gets H<sub>x</sub>, the entropy of the 
        ///   <see cref="RowMarginal"/> vector.
        /// </summary>
        /// 
        public double RowEntropy
        {
            get
            {
                if (this.xentropy == null)
                    this.xentropy = this.RowMarginal.Entropy(epsilon);
                return this.xentropy.Value;
            }
        }

        /// <summary>
        ///   Gets H<sub>y</sub>, the entropy of the 
        ///   <see cref="ColumnMarginal"/> vector.
        /// </summary>
        /// 
        public double ColumnEntropy
        {
            get
            {
                if (this.yentropy == null)
                    this.yentropy = this.ColumnMarginal.Entropy(epsilon);
                return this.yentropy.Value;
            }
        }

        /// <summary>
        ///   Gets p<sub>(x+y)</sub>(k), the sum 
        ///   of elements whose indices sum to k.
        /// </summary>
        /// 
        public double[] Sums
        {
            get
            {
                if (this.xysum == null)
                {
                    this.xysum = new double[2 * this.size];
                    for (int i = 0; i < this.size; i++)
                        for (int j = 0; j < this.size; j++)
                            this.xysum[i + j] += this.matrix[i, j];
                }
                return this.xysum;
            }
        }

        /// <summary>
        ///   Gets p<sub>(x-y)</sub> (k), the sum of elements 
        ///   whose absolute indices diferences equals to k.
        /// </summary>
        /// 
        public double[] Differences
        {
            get
            {
                if (this.xydiff == null)
                {
                    this.xydiff = new double[this.size];
                    for (int i = 0; i < this.size; i++)
                        for (int j = 0; j < this.size; j++)
                            this.xydiff[Math.Abs(i - j)] += this.matrix[i, j];
                }
                return this.xydiff;
            }
        }
        #endregion


        #region Feature Numbers

        /// <summary>
        ///   Gets Haralick's first textural feature,
        ///   the Angular Second Momentum.
        /// </summary>
        /// 
        public double F01 { get { return this.AngularSecondMomentum; } }

        /// <summary>
        ///   Gets Haralick's second textural feature,
        ///   the Contrast.
        /// </summary>
        /// 
        public double F02 { get { return this.Contrast; } }

        /// <summary>
        ///   Gets Haralick's third textural feature,
        ///   the Correlation.
        /// </summary>
        /// 
        public double F03 { get { return this.Correlation; } }

        /// <summary>
        ///   Gets Haralick's fourth textural feature,
        ///   the Sum of Squares: Variance.
        /// </summary>
        /// 
        public double F04 { get { return this.SumOfSquares; } }

        /// <summary>
        ///   Gets Haralick's fifth textural feature,
        ///   the Inverse Difference Moment.
        /// </summary>
        ///
        public double F05 { get { return this.InverseDifferenceMoment; } }

        /// <summary>
        ///   Gets Haralick's sixth textural feature,
        ///   the Sum Average.
        /// </summary>
        /// 
        public double F06 { get { return this.SumAverage; } }

        /// <summary>
        ///   Gets Haralick's seventh textural feature,
        ///   the Sum Variance.
        /// </summary>
        /// 
        public double F07 { get { return this.SumVariance; } }

        /// <summary>
        ///   Gets Haralick's eighth textural feature,
        ///   the Sum Entropy.
        /// </summary>
        /// 
        public double F08 { get { return this.SumEntropy; } }

        /// <summary>
        ///   Gets Haralick's ninth textural feature,
        ///   the Entropy.
        /// </summary>
        /// 
        public double F09 { get { return this.Entropy; } }

        /// <summary>
        ///   Gets Haralick's tenth textural feature,
        ///   the Difference Variance.
        /// </summary>
        /// 
        public double F10 { get { return this.DifferenceVariance; } }

        /// <summary>
        ///   Gets Haralick's eleventh textural feature,
        ///   the Difference Entropy.
        /// </summary>
        /// 
        public double F11 { get { return this.DifferenceEntropy; } }

        /// <summary>
        ///   Gets Haralick's twelfth textural feature,
        ///   the First Information Measure.
        /// </summary>
        /// 
        public double F12 { get { return this.FirstInformationMeasure; } }

        /// <summary>
        ///   Gets Haralick's thirteenth textural feature,
        ///   the Second Information Measure.
        /// </summary>
        /// 
        public double F13 { get { return this.SecondInformationMeasure; } }

        /// <summary>
        ///   Gets Haralick's fourteenth textural feature,
        ///   the Maximal Correlation Coefficient.
        /// </summary>
        /// 
        public double F14 { get { return this.MaximalCorrelationCoefficient; } }

        #endregion


        #region Features

        /// <summary>
        ///   Gets Haralick's first textural feature, the
        ///   Angular Second Momentum, also known as Energy
        ///   or Homogeneity.
        /// </summary>
        /// 
        public double AngularSecondMomentum
        {
            get
            {
                if (this.angular == null)
                {
                    double s = 0;
                    foreach (double v in this.matrix)
                        s += v * v;
                    this.angular = s;
                }
                return this.angular.Value;
            }
        }

        /// <summary>
        ///   Gets a variation of Haralick's second textural feature,
        ///   the Contrast with Absolute values (instead of squares).
        /// </summary>
        /// 
        public double LaplaceContrast
        {
            get
            {
                if (this.laplace == null)
                {
                    double[] p_xmy = this.Differences;

                    double s = 0;
                    for (int n = 0; n < p_xmy.Length; n++)
                        s += Math.Abs(n) * p_xmy[n];

                    this.laplace = s;
                }
                return this.laplace.Value;
            }
        }

        /// <summary>
        ///   Gets Haralick's second textural feature,
        ///   the Contrast.
        /// </summary>
        /// 
        public double Contrast
        {
            get
            {
                if (this.contrast == null)
                {
                    double[] p_xmy = this.Differences;

                    double s = 0;
                    for (int n = 0; n < p_xmy.Length; n++)
                        s += (n * n) * p_xmy[n];

                    this.contrast = s;
                    global::FileFormat.Accord.Core.Debug.Assert(!Double.IsNaN(this.contrast.Value));
                }
                return this.contrast.Value;
            }
        }

        /// <summary>
        ///   Gets Haralick's third textural feature,
        ///   the Correlation.
        /// </summary>
        /// 
        public double Correlation
        {
            get
            {
                if (this.correlation == null)
                {
                    double mx = this.RowMean;
                    double my = this.ColumnMean;
                    double sx = this.RowStandardDeviation;
                    double sy = this.ColumnStandardDeviation;
                    this.correlation = 0;

                    if (sx * sy > 0)
                    {
                        double s = 0;
                        for (int i = 0; i < this.size; i++)
                            for (int j = 0; j < this.size; j++)
                                s += (i * j) * this.matrix[i, j];

                        this.correlation = (s - mx * my) / (sx * sy);
                    }

                    global::FileFormat.Accord.Core.Debug.Assert(!Double.IsNaN(this.correlation.Value));
                }

                return this.correlation.Value;
            }
        }

        /// <summary>
        ///   Gets Haralick's fourth textural feature,
        ///   the Sum of Squares: Variance.
        /// </summary>
        /// 
        public double SumOfSquares
        {
            get
            {
                if (this.variance == null)
                {
                    double s = 0, mu = this.Mean;
                    for (int i = 0; i < this.size; i++)
                        for (int j = 0; j < this.size; j++)
                            s += (i - mu) * this.matrix[i, j];

                    this.variance = s;
                }
                return this.variance.Value;
            }
        }

        /// <summary>
        ///   Gets Haralick's fifth textural feature, the Inverse
        ///   Difference Moment, also known as Local Homogeneity.
        ///   Can be regarded as a complement to <see cref="Contrast"/>.
        /// </summary>
        ///
        public double InverseDifferenceMoment
        {
            get
            {
                if (this.inverse == null)
                {
                    double s = 0;
                    for (int i = 0; i < this.size; i++)
                        for (int j = 0; j < this.size; j++)
                            s += this.matrix[i, j] / (double)(1 + (i - j) * (i - j));
                    this.inverse = s;
                }
                return this.inverse.Value;
            }
        }

        /// <summary>
        ///   Gets a variation of Haralick's fifth textural feature,
        ///   the Texture Homogeneity. Can be regarded as a complement
        ///   to <see cref="LaplaceContrast"/>.
        /// </summary>
        /// 
        public double TextureHomogeneity
        {
            get
            {
                if (this.textureHomogeneity == null)
                {
                    double s = 0;
                    for (int i = 0; i < this.size; i++)
                        for (int j = 0; j < this.size; j++)
                            s += this.matrix[i, j] / (double)(1 + Math.Abs(i - j));
                    this.textureHomogeneity = s;
                }
                return this.textureHomogeneity.Value;
            }
        }

        /// <summary>
        ///   Gets Haralick's sixth textural feature,
        ///   the Sum Average.
        /// </summary>
        /// 
        public double SumAverage
        {
            get
            {
                if (this.sumAverage == null)
                {
                    double[] sums = this.Sums;
                    double s = 0;
                    for (int i = 0; i < sums.Length; i++)
                        s += i * sums[i];
                    this.sumAverage = s;

                    global::FileFormat.Accord.Core.Debug.Assert(!Double.IsNaN(this.sumAverage.Value));
                }
                return this.sumAverage.Value;
            }
        }

        /// <summary>
        ///   Gets Haralick's seventh textural feature,
        ///   the Sum Variance.
        /// </summary>
        /// 
        public double SumVariance
        {
            get
            {
                if (this.sumVariance == null)
                {
                    double[] sums = this.Sums;
                    double f8 = this.F08;
                    double s = 0;
                    for (int i = 0; i < sums.Length; i++)
                        s += (i - f8) * (i - f8) * sums[i];
                    this.sumVariance = s;

                    global::FileFormat.Accord.Core.Debug.Assert(!Double.IsNaN(this.sumVariance.Value));
                }
                return this.sumVariance.Value;
            }
        }

        /// <summary>
        ///   Gets Haralick's eighth textural feature,
        ///   the Sum Entropy.
        /// </summary>
        /// 
        public double SumEntropy
        {
            get
            {
                if (this.sumEntropy == null)
                    this.sumEntropy = this.Sums.Entropy(epsilon);
                global::FileFormat.Accord.Core.Debug.Assert(!Double.IsNaN(this.sumEntropy.Value));
                return this.sumEntropy.Value;
            }
        }

        /// <summary>
        ///   Gets Haralick's ninth textural feature,
        ///   the Entropy.
        /// </summary>
        /// 
        public double Entropy
        {
            get
            {
                if (this.entropy == null)
                    this.entropy = this.matrix.Entropy(epsilon);
                global::FileFormat.Accord.Core.Debug.Assert(!Double.IsNaN(this.entropy.Value));
                return this.entropy.Value;
            }
        }

        /// <summary>
        ///   Gets Haralick's tenth textural feature,
        ///   the Difference Variance.
        /// </summary>
        /// 
        public double DifferenceVariance
        {
            get
            {
                if (this.diffVariance == null)
                    this.diffVariance = this.Differences.Variance();
                global::FileFormat.Accord.Core.Debug.Assert(!Double.IsNaN(this.diffVariance.Value));
                return this.diffVariance.Value;
            }
        }

        /// <summary>
        ///   Gets Haralick's eleventh textural feature,
        ///   the Difference Entropy.
        /// </summary>
        /// 
        public double DifferenceEntropy
        {
            get
            {
                if (this.diffEntropy == null)
                    this.diffEntropy = this.Differences.Entropy(epsilon);
                global::FileFormat.Accord.Core.Debug.Assert(!Double.IsNaN(this.diffEntropy.Value));
                return this.diffEntropy.Value;
            }
        }

        /// <summary>
        ///   Gets Haralick's twelfth textural feature,
        ///   the First Information Measure.
        /// </summary>
        /// 
        public double FirstInformationMeasure
        {
            get
            {
                if (this.information1 == null)
                {
                    double hxy = this.Entropy;
                    double hxy1 = 0;

                    double[] px = this.RowMarginal;
                    double[] py = this.ColumnMarginal;

                    double hx = this.RowEntropy;
                    double hy = this.ColumnEntropy;

                    this.information1 = 0;

                    if (hx + hy > 0)
                    {
                        for (int i = 0; i < this.size; i++)
                            for (int j = 0; j < this.size; j++)
                                hxy1 -= this.matrix[i, j] * Math.Log(px[i] * py[j] + epsilon);

                        this.information1 = (hxy - hxy1) / Math.Max(hx, hy);
                    }

                    global::FileFormat.Accord.Core.Debug.Assert(!Double.IsNaN(this.information1.Value));
                }

                return this.information1.Value;
            }
        }

        /// <summary>
        ///   Gets Haralick's thirteenth textural feature,
        ///   the Second Information Measure.
        /// </summary>
        /// 
        public double SecondInformationMeasure
        {
            get
            {
                if (this.information2 == null)
                {
                    double hxy = this.Entropy;

                    double hxy2 = 0;

                    double[] px = this.RowMarginal;
                    double[] py = this.ColumnMarginal;

                    double hx = this.RowEntropy;
                    double hy = this.ColumnEntropy;

                    for (int i = 0; i < this.size; i++)
                        for (int j = 0; j < this.size; j++)
                            hxy2 -= px[i] * py[j] * Math.Log(px[i] * py[j] + epsilon);

                    this.information2 = Math.Sqrt(1.0 - Math.Exp(-2 * (hxy2 - hxy)));
                    global::FileFormat.Accord.Core.Debug.Assert(!Double.IsNaN(this.information2.Value));
                }

                return this.information2.Value;
            }
        }

        /// <summary>
        ///   Gets Haralick's fourteenth textural feature,
        ///   the Maximal Correlation Coefficient.
        /// </summary>
        /// 
        public double MaximalCorrelationCoefficient
        {
            get
            {
                if (this.maximal == null)
                {
                    double[] px = this.RowMarginal;
                    double[] py = this.ColumnMarginal;

                    double[,] Q = new double[this.size, this.size];
                    for (int i = 0; i < this.size; i++)
                    {
                        for (int j = 0; j < this.size; j++)
                        {
                            for (int k = 0; k < this.size; k++)
                            {
                                double num = this.matrix[i, k] * this.matrix[j, k];
                                double den = px[i] * py[i];
                                Q[i, j] += num / den;
                            }
                        }
                    }

                    var evd = new EigenvalueDecomposition(Q,
                        assumeSymmetric: false, inPlace: true);

                    this.maximal = Math.Sqrt(evd.RealEigenvalues[1]);
                }

                return this.maximal.Value;
            }
        }

        /// <summary>
        ///   Gets the Cluster Shade textural feature.
        /// </summary>
        /// 
        public double ClusterShade
        {
            get
            {
                if (this.clusterShade == null)
                {
                    double mx = this.RowMean;
                    double my = this.ColumnMean;

                    double s = 0;
                    for (int i = 0; i < this.size; i++)
                    {
                        for (int j = 0; j < this.size; j++)
                        {
                            double v = (i + j - mx - my);
                            s += (v * v * v) * this.matrix[i, j];
                        }
                    }

                    this.clusterShade = s;
                }

                return this.clusterShade.Value;
            }
        }

        /// <summary>
        ///   Gets the Cluster Prominence textural feature.
        /// </summary>
        /// 
        public double ClusterProminence
        {
            get
            {
                if (this.clusterProminence == null)
                {
                    double mx = this.RowMean;
                    double my = this.ColumnMean;

                    double s = 0;
                    for (int i = 0; i < this.size; i++)
                    {
                        for (int j = 0; j < this.size; j++)
                        {
                            double v = (i + j - mx - my);
                            s += (v * v * v * v) * this.matrix[i, j];
                        }
                    }

                    this.clusterProminence = s;
                }

                return this.clusterProminence.Value;
            }
        }

        #endregion


        /// <summary>
        ///   Creates a feature vector with 
        ///   the chosen feature functions.
        /// </summary>
        /// 
        /// <param name="features">How many features to include in the vector. Default is 13.</param>
        /// 
        /// <returns>A vector with Haralick's features up 
        /// to the given number passed as input.</returns>
        /// 
        public double[] GetVector(int features = 13)
        {
            double[] vector = new double[features];

            int i = features;

            switch (features)
            {
                case 01: vector[--i] = this.F01; break;
                case 02: vector[--i] = this.F02; goto case 1;
                case 03: vector[--i] = this.F03; goto case 2;
                case 04: vector[--i] = this.F04; goto case 3;
                case 05: vector[--i] = this.F05; goto case 4;
                case 06: vector[--i] = this.F06; goto case 5;
                case 07: vector[--i] = this.F07; goto case 6;
                case 08: vector[--i] = this.F08; goto case 7;
                case 09: vector[--i] = this.F09; goto case 8;
                case 10: vector[--i] = this.F10; goto case 9;
                case 11: vector[--i] = this.F11; goto case 10;
                case 12: vector[--i] = this.F12; goto case 11;
                case 13: vector[--i] = this.F13; goto case 12;
                case 14: vector[--i] = this.F14; goto case 13;

                default:
                    throw new ArgumentOutOfRangeException("features");
            }

            return vector;
        }
    }
}