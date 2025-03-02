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

namespace Openize.Accord.Statistics.Analysis.Performance
{
    using System;
    using System.ComponentModel;
    using Testing.Contingency;

    /// <summary>
    ///   Weighted confusion matrix for multi-class decision problems.
    /// </summary>
    /// 
    /// <remarks>
    /// <para>
    ///   References:
    ///   <list type="bullet">
    ///     <item><description>
    ///       <a href="http://uwf.edu/zhu/evr6930/2.pdf">
    ///       R.  G.  Congalton. A Review  of Assessing  the Accuracy  of Classifications 
    ///       of Remotely  Sensed  Data. Available on: http://uwf.edu/zhu/evr6930/2.pdf </a></description></item>
    ///     <item><description>
    ///       <a href="http://www.iiasa.ac.at/Admin/PUB/Documents/IR-98-081.pdf">
    ///       G. Banko. A Review of Assessing the Accuracy of Classiﬁcations of Remotely Sensed Data and
    ///       of Methods Including Remote Sensing Data in Forest Inventory. Interim report. Available on:
    ///       http://www.iiasa.ac.at/Admin/PUB/Documents/IR-98-081.pdf </a></description></item>
    ///     </list></para>  
    /// </remarks>
    /// 
    [Serializable]
    public class WeightedConfusionMatrix : GeneralConfusionMatrix
    {

        private double[,] weights;

        private double[] weightedRowProportion;
        private double[] weightedColProportion;

        private double? kappa;
        private double? kappaStdError;
        private double? kappaVariance;

        private double? kappaStdError0;
        private double? kappaVariance0;

        /// <summary>
        ///   Gets the Weights matrix.
        /// </summary>
        /// 
        public double[,] Weights
        {
            get { return this.weights; }
        }



        /// <summary>
        ///   Creates a new Confusion Matrix.
        /// </summary>
        /// 
        public WeightedConfusionMatrix(double[,] matrix, double[,] weights, int samples)
            : base(matrix, samples)
        {
            this.weights = weights;
        }

        /// <summary>
        ///   Creates a new Confusion Matrix.
        /// </summary>
        /// 
        public WeightedConfusionMatrix(int[,] matrix, double[,] weights)
            : base(matrix)
        {
            this.weights = weights;
        }

        /// <summary>
        ///   Creates a new Confusion Matrix.
        /// </summary>
        /// 
        public WeightedConfusionMatrix(double[,] weights, int[] expected, int[] predicted)
            : base(weights.GetLength(0), expected, predicted)
        {
            this.weights = weights;
        }

        /// <summary>
        ///   Gets the row marginals (proportions).
        /// </summary>
        /// 
        [DisplayName("Weighted Row Proportions")]
        public double[] WeightedRowProportions
        {
            get
            {
                if (this.weightedRowProportion == null)
                {
                    double[] rowProportion = this.RowProportions;
                    this.weightedRowProportion = new double[this.NumberOfClasses];
                    for (int i = 0; i < rowProportion.Length; i++)
                        for (int j = 0; j < rowProportion.Length; j++)
                            this.weightedRowProportion[i] += rowProportion[j] * this.weights[i, j];
                }

                return this.weightedRowProportion;
            }
        }

        /// <summary>
        ///   Gets the column marginals (proportions).
        /// </summary>
        /// 
        [DisplayName("Weighted Column Proportions")]
        public double[] WeightedColumnProportions
        {
            get
            {
                if (this.weightedColProportion == null)
                {
                    double[] colProportion = this.ColumnProportions;
                    this.weightedColProportion = new double[this.NumberOfClasses];
                    for (int i = 0; i < colProportion.Length; i++)
                        for (int j = 0; j < colProportion.Length; j++)
                            this.weightedColProportion[i] += colProportion[j] * this.weights[i, j];
                }

                return this.weightedColProportion;
            }
        }

        /// <summary>
        ///   Gets the Kappa coefficient of performance.
        /// </summary>
        /// 
        [DisplayName("Kappa Coefficient (κ)")]
        public double WeightedKappa
        {
            get
            {
                if (this.kappa == null)
                {
                    double Po = this.WeightedOverallAgreement;
                    double Pc = this.WeightedChanceAgreement;

                    this.kappa = (Po - Pc) / (1.0 - Pc);
                }

                return this.kappa.Value;
            }
        }

        /// <summary>
        ///   Gets the standard error of the <see cref="WeightedKappa"/>
        ///   coefficient of performance. 
        /// </summary>
        /// 
        [DisplayName("Kappa (κ) Std. Error")]
        public double WeightedStandardError
        {
            get
            {
                if (this.kappaStdError == null)
                {
                    double se;
                    this.kappaVariance = KappaTest.AsymptoticKappaVariance(this, out se);
                    this.kappaStdError = se;
                }

                return this.kappaStdError.Value;
            }
        }

        /// <summary>
        ///   Gets the variance of the <see cref="WeightedKappa"/>
        ///   coefficient of performance. 
        /// </summary>
        /// 
        [DisplayName("Kappa (κ) Variance")]
        public double WeightedVariance
        {
            get
            {
                if (this.kappaVariance == null)
                {
                    double se;
                    this.kappaVariance = KappaTest.AsymptoticKappaVariance(this, out se);
                    this.kappaStdError = se;
                }

                return this.kappaVariance.Value;
            }
        }

        /// <summary>
        ///   Gets the variance of the <see cref="WeightedKappa"/>
        ///   under the null hypothesis that the underlying
        ///   Kappa value is 0. 
        /// </summary>
        /// 
        [DisplayName("Kappa (κ) H₀ Variance")]
        public double WeightedVarianceUnderNull
        {
            get
            {
                if (this.kappaVariance0 == null)
                {
                    double se;
                    this.kappaVariance0 = KappaTest.AsymptoticKappaVariance(this,
                        out se, nullHypothesis: true);
                    this.kappaStdError0 = se;
                }

                return this.kappaVariance0.Value;
            }
        }

        /// <summary>
        ///   Gets the standard error of the <see cref="WeightedKappa"/>
        ///   under the null hypothesis that the underlying Kappa
        ///   value is 0. 
        /// </summary>
        /// 
        [DisplayName("Kappa (κ) H₀ Std. Error")]
        public double WeightedStandardErrorUnderNull
        {
            get
            {
                if (this.kappaStdError0 == null)
                {
                    double se;
                    this.kappaVariance0 = KappaTest.AsymptoticKappaVariance(this,
                        out se, nullHypothesis: true);
                    this.kappaStdError0 = se;
                }

                return this.kappaStdError0.Value;
            }
        }



        /// <summary>
        ///   Overall agreement.
        /// </summary>
        /// 
        [DisplayName("Overall Agreement")]
        public double WeightedOverallAgreement
        {
            get
            {
                double sum = 0;
                for (int i = 0; i < this.NumberOfClasses; i++)
                    for (int j = 0; j < this.NumberOfClasses; j++)
                        sum += this.weights[i, j] * this.Matrix[i, j];

                return sum / (double)this.NumberOfSamples;
            }
        }

        /// <summary>
        ///   Chance agreement.
        /// </summary>
        /// 
        /// <remarks>
        ///   The chance agreement tells how many samples
        ///   were correctly classified by chance alone.
        /// </remarks>
        /// 
        [DisplayName("Chance Agreement")]
        public double WeightedChanceAgreement
        {
            get
            {
                double chance = 0;
                for (int i = 0; i < this.NumberOfClasses; i++)
                    for (int j = 0; j < this.NumberOfClasses; j++)
                        chance += this.weights[i, j] 
                            * this.RowTotals[i] 
                            * this.ColumnTotals[j];
                return chance / (this.NumberOfSamples * this.NumberOfSamples);
            }
        }

        /// <summary>
        ///   Creates a new Weighted Confusion Matrix with linear weighting.
        /// </summary>
        /// 
        public static WeightedConfusionMatrix LinearWeighting(int[,] matrix)
        {
            // Create matrix of weights with linear weighting.
            int classes = matrix.GetLength(0);

            double[,] weights = new double[classes, classes];
            for (int i = 0; i < classes; i++)
            {
                for (int j = 0; j < classes; j++)
                {
                    double num = Math.Abs(i - j);
                    double den = classes-1;
                    weights[i, j] = 1.0 - num / den;
                }
            }

            return new WeightedConfusionMatrix(matrix, weights);
        }


        /// <summary>
        ///   Creates a new Weighted Confusion Matrix with linear weighting.
        /// </summary>
        /// 
        public static WeightedConfusionMatrix QuadraticWeighting(int[,] matrix)
        {
            // Create matrix of weights with quadratic weighting.
            int classes = matrix.GetLength(0);

            double[,] weights = new double[classes, classes];
            for (int i = 0; i < classes; i++)
            {
                for (int j = 0; j < classes; j++)
                {
                    double num = Math.Abs(i - j);
                    double den = classes - 1;
                    weights[i, j] = 1.0 - (num / den) * (num / den);
                }
            }

            return new WeightedConfusionMatrix(matrix, weights);
        }
    }
}
