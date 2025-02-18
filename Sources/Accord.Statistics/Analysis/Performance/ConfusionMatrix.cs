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

namespace Openize.Accord.Statistics.Analysis.Performance
{
    using System;
    using System.ComponentModel;
    using Openize.Accord.Math.Matrix;
    using Openize.Accord.Core.Exceptions;
    using Openize.Accord.Core.MachineLearning.Classifiers;
    using Openize.Accord.Math.Accord.Statistics;
    using Testing.Contingency;

    /// <summary>
    ///   Binary confusion matrix for binary decision problems. For multi-class 
    ///   decision problems, please see <see cref="GeneralConfusionMatrix"/>.
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
    /// <example>
    /// <para>
    ///   The first example shows how a confusion matrix can be constructed from a vector of expected (ground-truth)
    ///   values and their associated predictions (as done by a test, procedure or machine learning classifier):</para>
    ///   <code source="Unit Tests\Accord.Tests.Statistics\Analysis\GeneralConfusionMatrixTest.cs" region="doc_ctor_values" />
    /// <para>
    ///   The second example shows how to construct a binary confusion matrix directly from a classifier and a dataset:</para>
    ///   <code source="Unit Tests\Accord.Tests.MachineLearning\Boosting\AdaBoostTest.cs" region="doc_learn" />
    /// </example>
    /// 
    /// <seealso cref="GeneralConfusionMatrix"/>
    /// 
    [Serializable]
    public class ConfusionMatrix // TODO: Rename ConfusionMatrix to BinaryConfusionMatrix to make explicit it is only valid for binary problems
    {

        //  2x2 confusion matrix
        //
        //         | a(TP)    b(FN) |
        //   A =   |                |
        //         | c(FP)    d(TN) |
        //

        private int truePositives;  // a
        private int falseNegatives; // b
        private int falsePositives; // c
        private int trueNegatives;  // d


        // Association measures
        private double? kappaVariance;
        private double? kappaStdError;
        private double? kappaVariance0;
        private double? kappaStdError0;

        private double? chiSquare;



        /// <summary>
        ///   Constructs a new Confusion Matrix.
        /// </summary>
        /// 
        public ConfusionMatrix(
            int truePositives, int falseNegatives,
            int falsePositives, int trueNegatives)
        {
            this.truePositives = truePositives;
            this.trueNegatives = trueNegatives;
            this.falsePositives = falsePositives;
            this.falseNegatives = falseNegatives;
        }

        /// <summary>
        ///   Constructs a new Confusion Matrix.
        /// </summary>
        /// 
        public ConfusionMatrix(int[,] matrix)
        {
            if (matrix.Rows() != 2 || matrix.Columns() != 2)
                throw new DimensionMismatchException("matrix");

            this.truePositives = matrix[0, 0];
            this.falseNegatives = matrix[0, 1];

            this.falsePositives = matrix[1, 0];
            this.trueNegatives = matrix[1, 1];
        }

        /// <summary>
        ///   Constructs a new Confusion Matrix.
        /// </summary>
        /// 
        /// <param name="predicted">The values predicted by the model.</param>
        /// <param name="expected">The actual, truth values from the data.</param>
        /// 
        public ConfusionMatrix(bool[] predicted, bool[] expected)
        {
            // Initial argument checking
            if (predicted == null)
                throw new ArgumentNullException("predicted");
            if (expected == null)
                throw new ArgumentNullException("expected");
            if (predicted.Length != expected.Length)
                throw new DimensionMismatchException("expected", "The size of the predicted and expected arrays must match.");


            // For each of the predicted values,
            for (int i = 0; i < predicted.Length; i++)
            {
                bool prediction = predicted[i];
                bool expectation = expected[i];


                // If the prediction equals the true measured value
                if (expectation == prediction)
                {
                    // We have a hit. Now we have to see
                    //  if it was a positive or negative hit
                    if (prediction == true)
                    {
                        this.truePositives++; // Positive hit
                    }
                    else
                    {
                        this.trueNegatives++; // Negative hit
                    }
                }
                else
                {
                    // We have a miss. Now we have to see
                    //  if it was a positive or negative miss
                    if (prediction == true)
                    {
                        this.falsePositives++; // Positive hit
                    }
                    else
                    {
                        this.falseNegatives++; // Negative hit
                    }
                }
            }
        }

        /// <summary>
        ///   Constructs a new Confusion Matrix.
        /// </summary>
        /// 
        /// <param name="predicted">The values predicted by the model.</param>
        /// <param name="expected">The actual, truth values from the data.</param>
        /// 
        public ConfusionMatrix(int[] predicted, int[] expected)
        {
            int[] expectedSymbols = expected.Distinct();
            if (expectedSymbols.Length != 2)
            {
                throw new Exception("The vector of expected values contains more than two possible symbols. Please make sure the " +
                    "expected vector contains only values 0 or 1, or -1 and +1. If you have a multi-class decision problem, please " +
                    "use GeneralConfusionMatrix instead.");
            }

            if (predicted.DistinctCount() > 2)
                throw new Exception("The vector of predicted values contains more than two possible symbols. Please make sure the predicted vector contains only values 0 or 1, or -1 and +1.");

            predicted = Classes.ToZeroOne(predicted);
            expected = Classes.ToZeroOne(expected);

            this.compute(predicted, expected, positiveValue: 1);
        }

        /// <summary>
        ///   Constructs a new Confusion Matrix.
        /// </summary>
        /// 
        /// <param name="predicted">The values predicted by the model.</param>
        /// <param name="expected">The actual, truth values from the data.</param>
        /// <param name="positiveValue">The integer label which identifies a value as positive.</param>
        /// 
        public ConfusionMatrix(int[] predicted, int[] expected, int positiveValue = 1)
        {
            this.compute(predicted, expected, positiveValue);
        }

        private void compute(int[] predicted, int[] expected, int positiveValue)
        {
            // Initial argument checking
            if (predicted == null)
                throw new ArgumentNullException("predicted");
            if (expected == null)
                throw new ArgumentNullException("expected");
            if (predicted.Length != expected.Length)
                throw new DimensionMismatchException("expected", "The size of the predicted and expected arrays must match.");


            for (int i = 0; i < predicted.Length; i++)
            {
                bool prediction = predicted[i] == positiveValue;
                bool expectation = expected[i] == positiveValue;


                // If the prediction equals the true measured value
                if (expectation == prediction)
                {
                    // We have a hit. Now we have to see
                    //  if it was a positive or negative hit
                    if (prediction == true)
                    {
                        this.truePositives++; // Positive hit
                    }
                    else
                    {
                        this.trueNegatives++; // Negative hit
                    }
                }
                else
                {
                    // We have a miss. Now we have to see
                    //  if it was a positive or negative miss
                    if (prediction == true)
                    {
                        this.falsePositives++; // Positive hit
                    }
                    else
                    {
                        this.falseNegatives++; // Negative hit
                    }
                }
            }
        }


        /// <summary>
        ///   Constructs a new Confusion Matrix.
        /// </summary>
        /// 
        /// <param name="predicted">The values predicted by the model.</param>
        /// <param name="expected">The actual, truth values from the data.</param>
        /// <param name="positiveValue">The integer label which identifies a value as positive.</param>
        /// <param name="negativeValue">The integer label which identifies a value as negative.</param>
        /// 
        public ConfusionMatrix(int[] predicted, int[] expected, int positiveValue, int negativeValue)
        {
            // Initial argument checking
            if (predicted == null) throw new ArgumentNullException("predicted");
            if (expected == null) throw new ArgumentNullException("expected");
            if (predicted.Length != expected.Length)
                throw new DimensionMismatchException("expected", "The size of the predicted and expected arrays must match.");


            for (int i = 0; i < predicted.Length; i++)
            {

                // If the prediction equals the true measured value
                if (predicted[i] == expected[i])
                {
                    // We have a hit. Now we have to see
                    //  if it was a positive or negative hit
                    if (predicted[i] == positiveValue)
                    {
                        this.truePositives++; // Positive hit
                    }
                    else if (predicted[i] == negativeValue)
                    {
                        this.trueNegatives++; // Negative hit
                    }
                }
                else
                {
                    // We have a miss. Now we have to see
                    //  if it was a positive or negative miss
                    if (predicted[i] == positiveValue)
                    {
                        this.falsePositives++; // Positive hit
                    }
                    else if (predicted[i] == negativeValue)
                    {
                        this.falseNegatives++; // Negative hit
                    }
                }
            }
        }


        /// <summary>
        ///   Gets the confusion matrix in count matrix form.
        /// </summary>
        /// 
        /// <remarks>
        ///   The table is listed as true positives, false negatives
        ///   on its first row, false positives and true negatives in
        ///   its second row.
        /// </remarks>
        /// 
        public int[,] Matrix
        {
            get
            {
                return new int[,]
                {
                    { this.truePositives, this.falseNegatives },
                    { this.falsePositives, this.trueNegatives },
                };
            }
        }

        /// <summary>
        ///   Gets the marginal sums for table rows.
        /// </summary>
        /// 
        /// <value>
        ///   Returns a vector with the sum of true positives and 
        ///   false negatives on its first position, and the sum
        ///   of false positives and true negatives on the second.
        /// </value>
        /// 
        [DisplayName("Row Totals")]
        public int[] RowTotals
        {
            get
            {
                return new int[]
                {
                    this.truePositives + this.falseNegatives, // ActualPositives
                    this.falsePositives + this.trueNegatives, // ActualNegatives
                };
            }
        }

        /// <summary>
        ///   Gets the marginal sums for table columns.
        /// </summary>
        /// 
        /// <value>
        ///   Returns a vector with the sum of true positives and
        ///   false positives on its first position, and the sum
        ///   of false negatives and true negatives on the second.
        /// </value>
        /// 
        [DisplayName("Column Totals")]
        public int[] ColumnTotals
        {
            get
            {
                return new int[]
                {
                    this.truePositives + this.falsePositives, // PredictedPositives
                    this.falseNegatives + this.trueNegatives, // PredictedNegatives
                };
            }
        }


        /// <summary>
        ///   Gets the number of observations for this matrix.
        /// </summary>
        /// 
        [DisplayName("Number of samples")]
        public int NumberOfSamples
        {
            get
            {
                return this.trueNegatives + this.truePositives +
                    this.falseNegatives + this.falsePositives;
            }
        }

        /// <summary>
        ///   Gets the number of classes in this decision problem.
        /// </summary>
        /// 
        [DisplayName("Number of classes")]
        public int NumberOfClasses
        {
            get { return 2; }
        }

        /// <summary>
        ///   Gets the number of observations for this matrix
        /// </summary>
        /// 
        [Obsolete("Please use NumberOfSamples instead.")]
        [DisplayName("Number of samples")]
        public int Samples
        {
            get { return this.NumberOfSamples; }
        }

        /// <summary>
        ///   Gets the number of actual positives.
        /// </summary>
        /// 
        /// <remarks>
        ///   The number of positives cases can be computed by
        ///   taking the sum of true positives and false negatives.
        /// </remarks>
        /// 
        [DisplayName("Actual Positives")]
        public int ActualPositives
        {
            get { return this.truePositives + this.falseNegatives; }
        }

        /// <summary>
        ///   Gets the number of actual negatives
        /// </summary>
        /// 
        /// <remarks>
        ///   The number of negatives cases can be computed by
        ///   taking the sum of true negatives and false positives.
        /// </remarks>
        /// 
        [DisplayName("Actual Negatives")]
        public int ActualNegatives
        {
            get { return this.trueNegatives + this.falsePositives; }
        }

        /// <summary>
        ///   Gets the number of predicted positives.
        /// </summary>
        /// 
        /// <remarks>
        ///   The number of cases predicted as positive by the
        ///   test. This value can be computed by adding the
        ///   true positives and false positives.
        /// </remarks>
        /// 
        [DisplayName("Predicted Positives")]
        public int PredictedPositives
        {
            get { return this.truePositives + this.falsePositives; }
        }

        /// <summary>
        ///   Gets the number of predicted negatives.
        /// </summary>
        /// 
        /// <remarks>
        ///   The number of cases predicted as negative by the
        ///   test. This value can be computed by adding the
        ///   true negatives and false negatives.
        /// </remarks>
        /// 
        [DisplayName("Predicted Negatives")]
        public int PredictedNegatives
        {
            get { return this.trueNegatives + this.falseNegatives; }
        }

        /// <summary>
        ///   Cases correctly identified by the system as positives.
        /// </summary>
        /// 
        [DisplayName("True Positives (TP)")]
        public int TruePositives
        {
            get { return this.truePositives; }
        }

        /// <summary>
        ///   Cases correctly identified by the system as negatives.
        /// </summary>
        /// 
        [DisplayName("True Negatives (TN)")]
        public int TrueNegatives
        {
            get { return this.trueNegatives; }
        }

        /// <summary>
        ///   Cases incorrectly identified by the system as positives.
        /// </summary>
        /// 
        [DisplayName("False Positives (FP)")]
        public int FalsePositives
        {
            get { return this.falsePositives; }
        }

        /// <summary>
        ///   Cases incorrectly identified by the system as negatives.
        /// </summary>
        /// 
        [DisplayName("False Negatives (FN)")]
        public int FalseNegatives
        {
            get { return this.falseNegatives; }
        }

        /// <summary>
        ///   Sensitivity, also known as True Positive Rate
        /// </summary>
        /// 
        /// <remarks>
        ///   The Sensitivity is calculated as <c>TPR = TP / (TP + FN)</c>.
        /// </remarks>
        /// 
        public double Sensitivity
        {
            get
            {
                return (this.truePositives == 0) ?
                    0 : (double)this.truePositives / (this.truePositives + this.falseNegatives);
            }
        }

        /// <summary>
        ///   Specificity, also known as True Negative Rate
        /// </summary>
        /// 
        /// <remarks>
        ///   The Specificity is calculated as <c>TNR = TN / (FP + TN)</c>.
        ///   It can also be calculated as: <c>TNR = (1-False Positive Rate)</c>.
        /// </remarks>
        /// 
        public double Specificity
        {
            get
            {
                return (this.trueNegatives == 0) ?
                    0 : (double)this.trueNegatives / (this.trueNegatives + this.falsePositives);
            }
        }

        /// <summary>
        ///  Efficiency, the arithmetic mean of sensitivity and specificity
        /// </summary>
        /// 
        public double Efficiency
        {
            get { return (this.Sensitivity + this.Specificity) / 2.0; }
        }

        /// <summary>
        ///   Gets the number of errors between the expected and predicted values.
        /// </summary>
        /// 
        public int Errors
        {
            get { return this.falseNegatives + this.falsePositives; }
        }

        /// <summary>
        ///   Gets the number of hits between the expected and predicted values.
        /// </summary>
        /// 
        public int Hits
        {
            get { return this.truePositives + this.trueNegatives; }
        }

        /// <summary>
        ///   Accuracy, or raw performance of the system
        /// </summary>
        /// 
        /// <remarks>
        ///   The Accuracy is calculated as 
        ///   <c>ACC = (TP + TN) / (P + N).</c>
        /// </remarks>
        /// 
        public double Accuracy
        {
            get { return (this.truePositives + this.trueNegatives) / (double)this.NumberOfSamples; }
        }

        /// <summary>
        ///   Error rate, or 1 - accuracy.
        /// </summary>
        /// 
        /// <remarks>
        ///   The Error is calculated as 
        ///   <c>ACC = (FP + FN) / (P + N).</c>
        /// </remarks>
        /// 
        public double Error
        {
            get { return (this.falsePositives + this.falseNegatives) / (double)this.NumberOfSamples; }
        }

        /// <summary>
        ///  Prevalence of outcome occurrence.
        /// </summary>
        /// 
        public double Prevalence
        {
            get { return this.ActualPositives / (double)this.NumberOfSamples; }
        }

        /// <summary>
        ///   Gets the diagonal of the confusion matrix.
        /// </summary>
        /// 
        public int[] Diagonal
        {
            get { return new int[] { this.truePositives, this.trueNegatives }; }
        }

        /// <summary>
        ///   Positive Predictive Value, also known as Positive Precision
        /// </summary>
        /// 
        /// <remarks>
        /// <para>
        ///   The Positive Predictive Value tells us how likely is 
        ///   that a patient has a disease, given that the test for
        ///   this disease is positive.</para>
        /// <para>
        ///   The Positive Predictive Rate is calculated as
        ///   <c>PPV = TP / (TP + FP)</c>.</para>
        /// </remarks>
        /// 
        [DisplayName("Positive Predictive Value (PPV)")]
        public double PositivePredictiveValue
        {
            get
            {
                double f = this.truePositives + this.FalsePositives;

                if (f != 0)
                    return this.truePositives / f;

                return 1.0;
            }
        }

        /// <summary>
        ///   Negative Predictive Value, also known as Negative Precision
        /// </summary>
        /// 
        /// <remarks>
        /// <para>
        ///   The Negative Predictive Value tells us how likely it is
        ///   that the disease is NOT present for a patient, given that
        ///   the patient's test for the disease is negative.</para>
        /// <para>
        ///   The Negative Predictive Value is calculated as 
        ///   <c>NPV = TN / (TN + FN)</c>.</para> 
        /// </remarks>
        /// 
        [DisplayName("Negative Predictive Value (NPV)")]
        public double NegativePredictiveValue
        {
            get
            {
                double f = (this.trueNegatives + this.falseNegatives);

                if (f != 0)
                    return this.trueNegatives / f;

                return 1.0;
            }
        }


        /// <summary>
        ///   False Positive Rate, also known as false alarm rate.
        /// </summary>
        /// 
        /// <remarks>
        /// <para>
        ///   The rate of false alarms in a test.</para>
        /// <para>
        ///   The False Positive Rate can be calculated as
        ///   <c>FPR = FP / (FP + TN)</c> or <c>FPR = (1-specificity)</c>.
        /// </para>
        /// </remarks>
        /// 
        [DisplayName("False Positive Rate")]
        public double FalsePositiveRate
        {
            get
            {
                return (double)this.falsePositives / (this.falsePositives + this.trueNegatives);
            }
        }

        /// <summary>
        ///   False Discovery Rate, or the expected false positive rate.
        /// </summary>
        /// 
        /// <remarks>
        /// <para>
        ///   The False Discovery Rate is actually the expected false positive rate.</para>
        /// <para>
        ///   For example, if 1000 observations were experimentally predicted to
        ///   be different, and a maximum FDR for these observations was 0.10, then
        ///   100 of these observations would be expected to be false positives.</para>
        /// <para>
        ///   The False Discovery Rate is calculated as
        ///   <c>FDR = FP / (FP + TP)</c>.</para>
        /// </remarks>
        /// 
        [DisplayName("False Discovery Rate")]
        public double FalseDiscoveryRate
        {
            get
            {
                double d = this.falsePositives + this.truePositives;

                if (d != 0.0)
                    return this.falsePositives / d;

                return 1.0;
            }
        }

        /// <summary>
        ///   Matthews Correlation Coefficient, also known as Phi coefficient 
        /// </summary>
        /// 
        /// <remarks>
        ///   A coefficient of +1 represents a perfect prediction, 0 an
        ///   average random prediction and −1 an inverse prediction.
        /// </remarks>
        /// 
        [DisplayName("Mattews Correlation (φ)")]
        public double MatthewsCorrelationCoefficient
        {
            get
            {
                double tp = this.truePositives;
                double tn = this.trueNegatives;
                double fp = this.falsePositives;
                double fn = this.falseNegatives;

                double den = System.Math.Sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn));

                if (den == 0.0)
                    return 0;

                double num = (tp * tn) - (fp * fn);

                return num / den;
            }
        }

        /// <summary>
        ///   Pearson's contingency coefficient C.
        /// </summary>
        /// 
        /// <remarks>
        /// <para>
        ///   Pearson's C measures the degree of association between the two variables. However,
        ///   C suffers from the disadvantage that it does not reach a maximum of 1 or the minimum 
        ///   of -1; the highest it can reach in a 2 x 2 table is .707; the maximum it can reach in
        ///   a 4 × 4 table is 0.870. It can reach values closer to 1 in contingency tables with more
        ///   categories. It should, therefore, not be used to compare associations among tables with
        ///   different numbers of categories.</para>
        ///   
        /// <para>
        ///   References:
        ///     http://en.wikipedia.org/wiki/Contingency_table </para>
        /// </remarks>
        /// 
        [DisplayName("Pearson's C")]
        public double Pearson
        {
            get { return Math.Sqrt(this.ChiSquare / (this.ChiSquare + this.NumberOfSamples)); }
        }

        /// <summary>
        ///   Geometric agreement.
        /// </summary>
        /// 
        /// <remarks>
        ///   The geometric agreement is the geometric mean of the
        ///   diagonal elements of the confusion matrix.
        /// </remarks>
        /// 
        [DisplayName("Geometric Agreement")]
        public double GeometricAgreement
        {
            get { return Math.Exp(Measures.LogGeometricMean(this.Diagonal)); }
        }

        /// <summary>
        ///   Odds-ratio.
        /// </summary>
        /// 
        /// <remarks>
        ///   References: http://www.iph.ufrgs.br/corpodocente/marques/cd/rd/presabs.htm
        /// </remarks>
        /// 
        [DisplayName("Odds Ratio")]
        public double OddsRatio
        {
            get
            {
                return (double)(this.truePositives * this.trueNegatives) / (this.falsePositives * this.falseNegatives);
            }
        }

        /// <summary>
        ///   Overall agreement.
        /// </summary>
        /// 
        /// <remarks>
        ///   The overall agreement is the sum of the diagonal elements
        ///   of the contingency table divided by the number of samples.
        /// </remarks>
        /// 
        [DisplayName("Overall Agreement")]
        public double OverallAgreement
        {
            get
            {
                double a = this.truePositives;
                double d = this.trueNegatives;
                double N = this.NumberOfSamples;

                return (a + d) / N;
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
        public double ChanceAgreement
        {
            get
            {
                double a = this.truePositives;
                double b = this.falsePositives;
                double c = this.falseNegatives;
                double d = this.trueNegatives;
                double N = this.NumberOfSamples;

                return (((a + c) * (a + b) + (b + d) * (c + d)) / (N * N));
            }
        }

        /// <summary>
        ///   Kappa coefficient.
        /// </summary>
        ///
        /// <remarks>
        ///   References: http://www.iph.ufrgs.br/corpodocente/marques/cd/rd/presabs.htm
        /// </remarks>
        ///
        [DisplayName("Kappa Coefficient (κ)")]
        public double Kappa
        {
            get
            {
                double Po = this.OverallAgreement;
                double Pc = this.ChanceAgreement;
                return (Po - Pc) / (1.0 - Pc);
            }
        }

        /// <summary>
        ///   Gets the standard error of the <see cref="Kappa"/>
        ///   coefficient of performance. 
        /// </summary>
        /// 
        [DisplayName("Kappa (κ) Std. Error")]
        public double StandardError
        {
            get
            {
                if (this.kappaStdError == null)
                {
                    double se;
                    this.kappaVariance = KappaTest.AsymptoticKappaVariance(this.ToGeneralMatrix(), out se);
                    this.kappaStdError = se;
                }

                return this.kappaStdError.Value;
            }
        }

        /// <summary>
        ///   Gets the variance of the <see cref="Kappa"/>
        ///   coefficient of performance. 
        /// </summary>
        /// 
        [DisplayName("Kappa (κ) Variance")]
        public double Variance
        {
            get
            {
                if (this.kappaVariance == null)
                {
                    double se;
                    this.kappaVariance = KappaTest.AsymptoticKappaVariance(this.ToGeneralMatrix(), out se);
                    this.kappaStdError = se;
                }

                return this.kappaVariance.Value;
            }
        }

        /// <summary>
        ///   Gets the variance of the <see cref="Kappa"/>
        ///   under the null hypothesis that the underlying
        ///   Kappa value is 0. 
        /// </summary>
        /// 
        [DisplayName("Kappa (κ) H₀ Variance")]
        public double VarianceUnderNull
        {
            get
            {
                if (this.kappaVariance0 == null)
                {
                    double se;
                    this.kappaVariance0 = KappaTest.AsymptoticKappaVariance(this.ToGeneralMatrix(),
                        out se, nullHypothesis: true);
                    this.kappaStdError0 = se;
                }

                return this.kappaVariance0.Value;
            }
        }

        /// <summary>
        ///   Gets the standard error of the <see cref="Kappa"/>
        ///   under the null hypothesis that the underlying Kappa
        ///   value is 0. 
        /// </summary>
        /// 
        [DisplayName("Kappa (κ) H₀ Std. Error")]
        public double StandardErrorUnderNull
        {
            get
            {
                if (this.kappaStdError0 == null)
                {
                    double se;
                    this.kappaVariance0 = KappaTest.AsymptoticKappaVariance(this.ToGeneralMatrix(),
                        out se, nullHypothesis: true);
                    this.kappaStdError0 = se;
                }

                return this.kappaStdError0.Value;
            }
        }

        /// <summary>
        ///   Diagnostic power.
        /// </summary>
        /// 
        [DisplayName("Overall Diagnostic Power")]
        public double OverallDiagnosticPower
        {
            get { return (double)(this.falsePositives + this.trueNegatives) / this.NumberOfSamples; }
        }

        /// <summary>
        ///   Normalized Mutual Information.
        /// </summary>
        /// 
        [DisplayName("Normalized Mutual Information")]
        public double NormalizedMutualInformation
        {
            get
            {
                double a = this.truePositives;
                double b = this.falsePositives;
                double c = this.falseNegatives;
                double d = this.trueNegatives;
                double N = this.NumberOfSamples;

                double num = a * Math.Log(a) + b * Math.Log(b) + c * Math.Log(c) + d * Math.Log(d)
                           - (a + b) * Math.Log(a + b) - (c + d) * Math.Log(c + d);

                double den = N * Math.Log(N) - ((a + c) * Math.Log(a + c) + (b + d) * Math.Log(b + d));

                return 1 + num / den;
            }
        }

        /// <summary>
        ///   Precision, same as the <see cref="PositivePredictiveValue"/>.
        /// </summary>
        /// 
        public double Precision
        {
            get { return this.PositivePredictiveValue; }
        }

        /// <summary>
        ///   Recall, same as the <see cref="Sensitivity"/>.
        /// </summary>
        /// 
        public double Recall
        {
            get { return this.Sensitivity; }
        }

        /// <summary>
        ///   F-Score, computed as the harmonic mean of
        ///   <see cref="Precision"/> and <see cref="Recall"/>.
        /// </summary>
        /// 
        [DisplayName("F-Score")]
        public double FScore
        {
            get { return 2.0 * (this.Precision * this.Recall) / (this.Precision + this.Recall); }
        }

        /// <summary>
        ///   Expected values, or values that could
        ///   have been generated just by chance.
        /// </summary>
        /// 
        [DisplayName("Expected Values")]
        public double[,] ExpectedValues
        {
            get
            {
                var row = this.RowTotals;
                var col = this.ColumnTotals;

                var expected = new double[2, 2];

                for (int i = 0; i < row.Length; i++)
                    for (int j = 0; j < col.Length; j++)
                        expected[i, j] = (col[j] * row[i]) / (double)this.NumberOfSamples;

                return expected;
            }
        }

        /// <summary>
        ///   Gets the Chi-Square statistic for the contingency table.
        /// </summary>
        /// 
        [DisplayName("Chi-Square (χ²)")]
        public double ChiSquare
        {
            get
            {
                if (this.chiSquare == null)
                {
                    var row = this.RowTotals;
                    var col = this.ColumnTotals;

                    double x = 0;
                    for (int i = 0; i < row.Length; i++)
                    {
                        for (int j = 0; j < col.Length; j++)
                        {
                            double e = (row[i] * col[j]) / (double)this.NumberOfSamples;
                            double o = this.Matrix[i, j];

                            double num = (o - e) * (o - e);
                            if (num != 0)
                                x += num / e;
                        }
                    }

                    this.chiSquare = x;
                }

                return this.chiSquare.Value;
            }
        }

        /// <summary>
        ///   Returns a <see cref="System.String"/> representing this confusion matrix.
        /// </summary>
        /// 
        /// <returns>
        ///   A <see cref="System.String"/> representing this confusion matrix.
        /// </returns>
        /// 
        public override string ToString()
        {
            return String.Format(System.Globalization.CultureInfo.CurrentCulture,
                "TP:{0} FP:{2}, FN:{3} TN:{1}",
                this.truePositives, this.trueNegatives, this.falsePositives, this.falseNegatives);
        }

        /// <summary>
        ///   Converts this matrix into a <see cref="GeneralConfusionMatrix"/>.
        /// </summary>
        /// 
        /// <returns>
        ///   A <see cref="GeneralConfusionMatrix"/> with the same contents as this matrix.
        /// </returns>
        /// 
        public GeneralConfusionMatrix ToGeneralMatrix()
        {
            // Create a new matrix assuming negative instances 
            // are class 0, and positive instances are class 1.

            int[,] matrix =
            {
                //   class 0          class 1
                { this.trueNegatives,  this.falsePositives }, // class 0
                { this.falseNegatives, this.truePositives  }, // class 1
            };

            return new GeneralConfusionMatrix(matrix);
        }

        /// <summary>
        ///   Combines several confusion matrices into one single matrix.
        /// </summary>
        /// 
        /// <param name="matrices">The matrices to combine.</param>
        /// 
        public static ConfusionMatrix Combine(params ConfusionMatrix[] matrices)
        {
            if (matrices == null)
                throw new ArgumentNullException("matrices");
            if (matrices.Length == 0)
                throw new ArgumentException("At least one confusion matrix is required.");

            int[,] total = new int[2, 2];

            foreach (var matrix in matrices)
            {
                for (int j = 0; j < 2; j++)
                    for (int k = 0; k < 2; k++)
                        total[j, k] += matrix.Matrix[j, k];
            }

            return new ConfusionMatrix(total);
        }

        /// <summary>
        ///   Estimates a <see cref="ConfusionMatrix"/> directly from a classifier, a set of inputs and its expected outputs.
        /// </summary>
        /// 
        /// <typeparam name="TInput">The type of the inputs accepted by the classifier.</typeparam>
        /// 
        /// <param name="classifier">The classifier.</param>
        /// <param name="inputs">The input vectors.</param>
        /// <param name="expected">The expected outputs associated with each input vector.</param>
        /// 
        /// <returns>A <see cref="ConfusionMatrix"/> capturing the performance of the classifier when
        ///   trying to predict the outputs <paramref name="expected"/> from the <paramref name="inputs"/>.</returns>
        ///   
        public static ConfusionMatrix Estimate<TInput>(IClassifier<TInput, int> classifier, TInput[] inputs, int[] expected)
        {
            return new ConfusionMatrix(expected: expected, predicted: classifier.Decide(inputs));
        }

        /// <summary>
        ///   Estimates a <see cref="ConfusionMatrix"/> directly from a classifier, a set of inputs and its expected outputs.
        /// </summary>
        /// 
        /// <typeparam name="TInput">The type of the inputs accepted by the classifier.</typeparam>
        /// 
        /// <param name="classifier">The classifier.</param>
        /// <param name="inputs">The input vectors.</param>
        /// <param name="expected">The expected outputs associated with each input vector.</param>
        /// 
        /// <returns>A <see cref="ConfusionMatrix"/> capturing the performance of the classifier when
        ///   trying to predict the outputs <paramref name="expected"/> from the <paramref name="inputs"/>.</returns>
        ///   
        public static ConfusionMatrix Estimate<TInput>(IClassifier<TInput, bool> classifier, TInput[] inputs, bool[] expected)
        {
            return new ConfusionMatrix(expected: expected, predicted: classifier.Decide(inputs));
        }
    }
}
