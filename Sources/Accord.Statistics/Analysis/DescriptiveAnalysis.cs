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

namespace FileFormat.Accord.Statistics.Analysis
{
    using System;
    using System.ComponentModel;
    using Accord.MachineLearning.Learning;
    using Base;
    using Distributions.Univariate.Continuous;
    using FileFormat.Accord.Core.Collections;
    using FileFormat.Accord.Core.Ranges;
    using Math.Accord.Statistics;
    using Measures;
    using FileFormat.Accord.Math.Matrix;

    /// <summary>
    ///   Descriptive statistics analysis.
    /// </summary>
    /// 
    /// <remarks>
    /// <para>
    ///   Descriptive statistics are used to describe the basic features of the data
    ///   in a study. They provide simple summaries about the sample and the measures.
    ///   Together with simple graphics analysis, they form the basis of virtually
    ///   every quantitative analysis of data.</para>
    ///   
    /// <para>
    ///   This class can also be bound to standard controls such as the 
    ///   <a href="http://msdn.microsoft.com/en-us/library/system.windows.forms.datagridview.aspx">DataGridView</a>
    ///   by setting their DataSource property to the analysis' <see cref="Measures"/> property.</para>
    ///   
    /// <para>
    ///   References:
    ///   <list type="bullet">
    ///     <item><description>
    ///        Wikipedia, The Free Encyclopedia. Descriptive Statistics. Available on:
    ///        http://en.wikipedia.org/wiki/Descriptive_statistics </description></item>
    ///   </list></para>
    /// </remarks>
    ///
    /// <example>
    /// <code source="Unit Tests\Accord.Tests.Statistics\Analysis\DescriptiveAnalysisTest.cs" region="doc_learn" />
    /// </example>
    /// 
    /// <seealso cref="Tools"/>
    /// <seealso cref="DescriptiveMeasures"/>
    ///
    [Serializable]
#pragma warning disable 612, 618
    public class DescriptiveAnalysis : IMultivariateAnalysis,
        IDescriptiveLearning<DescriptiveAnalysis, double[]>
#pragma warning restore 612, 618
    {

        private int samples;
        private int variables;

        private double[] sums;
        private double[] means;
        private double[] standardDeviations;
        private double[] variances;
        private double[] medians;
        private double[] modes;
        private int[] distinct;

        private string[] columnNames;

        private QuantileMethod quantileMethod = QuantileMethod.Default;
        private DoubleRange[] ranges;
        private DoubleRange[] quartiles;
        private DoubleRange[] innerFences;
        private DoubleRange[] outerFences;
        private DoubleRange[] confidence;
        private DoubleRange[] deviance;

        private double[] kurtosis;
        private Double[] skewness;
        private double[] standardErrors;


        private double[,] covarianceMatrix;
        private double[,] correlationMatrix;

        private double[,] zScores;
        private double[,] dScores;

        private double[,] sourceMatrix;
        private double[][] sourceArray;

        private DescriptiveMeasureCollection measuresCollection;

        private bool lazy = true;

        /// <summary>
        ///   Constructs the Descriptive Analysis.
        /// </summary>
        /// 
        public DescriptiveAnalysis()
        {
        }

        /// <summary>
        ///   Constructs the Descriptive Analysis.
        /// </summary>
        /// 
        /// <param name="columnNames">Names for the analyzed variables.</param>
        /// 
        public DescriptiveAnalysis(string[] columnNames)
        {
            if (columnNames == null)
                throw new ArgumentNullException("columnNames");

            this.init(null, null, columnNames);
        }

        /// <summary>
        ///   Constructs the Descriptive Analysis.
        /// </summary>
        /// 
        /// <param name="data">The source data to perform analysis.</param>
        /// 
        [Obsolete("Please call the Learn() method passing the data to be analyzed.")]
        public DescriptiveAnalysis(double[] data)
        {
            if (data == null)
                throw new ArgumentNullException("data");

            double[,] matrix = new double[data.Length, 1];

            System.Buffer.BlockCopy(data, 0, matrix, 0, data.Length * sizeof(double));

            this.init(matrix, null, null);
        }

        /// <summary>
        ///   Constructs the Descriptive Analysis.
        /// </summary>
        /// 
        /// <param name="data">The source data to perform analysis.</param>
        /// 
        [Obsolete("Please call the Learn() method passing the data to be analyzed.")]
        public DescriptiveAnalysis(double[,] data)
        {
            if (data == null)
                throw new ArgumentNullException("data");

            this.init(data, null, null);
        }

        /// <summary>
        ///   Constructs the Descriptive Analysis.
        /// </summary>
        /// 
        /// <param name="data">The source data to perform analysis.</param>
        /// <param name="columnNames">Names for the analyzed variables.</param>
        /// 
        [Obsolete("Please pass only columnNames and call the Learn() method passing the data to be analyzed.")]
        public DescriptiveAnalysis(double[,] data, string[] columnNames)
        {
            if (data == null)
                throw new ArgumentNullException("data");

            if (columnNames == null)
                throw new ArgumentNullException("columnNames");

            this.init(data, null, columnNames);
        }

        /// <summary>
        ///   Constructs the Descriptive Analysis.
        /// </summary>
        /// 
        /// <param name="data">The source data to perform analysis.</param>
        /// 
        [Obsolete("Please call the Learn() method passing the data to be analyzed.")]
        public DescriptiveAnalysis(double[][] data)
        {
            // Initial argument checking
            if (data == null)
                throw new ArgumentNullException("data");

            this.init(null, data, null);
        }

        /// <summary>
        ///   Constructs the Descriptive Analysis.
        /// </summary>
        /// 
        /// <param name="data">The source data to perform analysis.</param>
        /// <param name="columnNames">Names for the analyzed variables.</param>
        /// 
        [Obsolete("Please pass only columnNames and call the Learn() method passing the data to be analyzed.")]
        public DescriptiveAnalysis(double[][] data, string[] columnNames)
        {
            // Initial argument checking
            if (data == null)
                throw new ArgumentNullException("data");

            if (columnNames == null)
                throw new ArgumentNullException("columnNames");

            this.init(null, data, columnNames);
        }

        private void init(double[,] matrix, double[][] array, string[] columnNames)
        {
            this.sourceMatrix = matrix;
            this.sourceArray = array;
            this.columnNames = columnNames;

            if (matrix != null)
            {
                this.samples = matrix.GetLength(0);
                this.variables = matrix.GetLength(1);
            }
            else if (array != null)
            {
                this.samples = array.Length;
                this.variables = array[0].Length;
            }
            else
            {
                return;
            }

            // Create object-oriented structure to access data
            DescriptiveMeasures[] measures = new DescriptiveMeasures[this.variables];
            for (int i = 0; i < measures.Length; i++)
                measures[i] = new DescriptiveMeasures(this, i);
            this.measuresCollection = new DescriptiveMeasureCollection(measures);
        }


        /// <summary>
        ///   Computes the analysis using given source data and parameters.
        /// </summary>
        /// 
        [Obsolete("Please use Learn() instead.")]
        public void Compute()
        {
            // Clear analysis
            this.reset();

            this.sums = this.Sums;
            this.means = this.Means;
            this.standardDeviations = this.StandardDeviations;
            this.ranges = this.Ranges;
            this.kurtosis = this.Kurtosis;
            this.skewness = this.Skewness;
            this.medians = this.Medians;
            this.modes = this.Modes;
            this.variances = this.Variances;
            this.standardErrors = this.StandardErrors;
            this.distinct = this.Distinct;
            this.quartiles = this.Quartiles;
            this.innerFences = this.InnerFences;
            this.outerFences = this.OuterFences;

            // Mean centered and standardized data
            this.dScores = this.DeviationScores;
            this.zScores = this.StandardScores;

            // Covariance and correlation
            this.covarianceMatrix = this.CovarianceMatrix;
            this.correlationMatrix = this.CorrelationMatrix;

            this.confidence = this.Confidence;
            this.deviance = this.Deviance;
        }

        private void reset()
        {
            this.sums = null;
            this.means = null;
            this.standardDeviations = null;
            this.ranges = null;
            this.kurtosis = null;
            this.skewness = null;
            this.medians = null;
            this.modes = null;
            this.variances = null;
            this.standardErrors = null;
            this.distinct = null;
            this.dScores = null;
            this.zScores = null;
            this.covarianceMatrix = null;
            this.correlationMatrix = null;
            this.deviance = null;
            this.quartiles = null;
            this.innerFences = null;
            this.outerFences = null;
        }

        /// <summary>
        /// Learns a model that can map the given inputs to the desired outputs.
        /// </summary>
        /// <param name="x">The model inputs.</param>
        /// <returns>
        /// A model that has learned how to produce suitable outputs
        /// given the input data <paramref name="x" />.
        /// </returns>
        public DescriptiveAnalysis Learn(double[][] x)
        {
            this.reset();

            this.init(null, x, this.columnNames);

            if (!this.lazy)
            {
                this.sums = this.Sums;
                this.means = this.Means;
                this.standardDeviations = this.StandardDeviations;
                this.ranges = this.Ranges;
                this.kurtosis = this.Kurtosis;
                this.skewness = this.Skewness;
                this.medians = this.Medians;
                this.modes = this.Modes;
                this.variances = this.Variances;
                this.standardErrors = this.StandardErrors;
                this.distinct = this.Distinct;
                this.quartiles = this.Quartiles;
                this.innerFences = this.InnerFences;
                this.outerFences = this.OuterFences;

                // Mean centered and standardized data
                this.dScores = this.DeviationScores;
                this.zScores = this.StandardScores;

                // Covariance and correlation
                this.covarianceMatrix = this.CovarianceMatrix;
                this.correlationMatrix = this.CorrelationMatrix;

                this.confidence = this.Confidence;
                this.deviance = this.Deviance;

                this.sourceArray = null;
                this.sourceMatrix = null;
            }

            return this;
        }

        /// <summary>
        ///   Gets or sets whether the properties of this class should
        ///   be computed only when necessary. If set to true, a copy
        ///   of the input data will be maintained inside an instance
        ///   of this class, using more memory.
        /// </summary>
        /// 
        public bool Lazy
        {
            get { return this.lazy; }
            set { this.lazy = value; }
        }

        /// <summary>
        ///   Gets the source matrix from which the analysis was run.
        /// </summary>
        /// 
        [Obsolete("This property will be removed.")]
        public double[,] Source
        {
            get
            {
                if (this.sourceMatrix == null)
                    this.sourceMatrix = Matrix.ToMatrix(this.sourceArray);
                return this.sourceMatrix;
            }
        }

        /// <summary>
        ///   Gets the source matrix from which the analysis was run.
        /// </summary>
        /// 
        [Obsolete("This property will be removed.")]
        public double[][] Array
        {
            get
            {
                if (this.sourceArray == null)
                    this.sourceArray = Matrix.ToJagged(this.sourceMatrix);
                return this.sourceArray;
            }
        }

        /// <summary>
        ///   Gets or sets the method to be used when computing quantiles (median and quartiles).
        /// </summary>
        /// 
        /// <value>The quantile method.</value>
        /// 
        public QuantileMethod QuantileMethod
        {
            get { return this.quantileMethod; }
            set { this.quantileMethod = value; }
        }

        /// <summary>
        ///   Gets the column names from the variables in the data.
        /// </summary>
        /// 
        public String[] ColumnNames
        {
            get
            {
                if (this.columnNames == null)
                {
                    this.columnNames = new string[this.variables];
                    for (int i = 0; i < this.columnNames.Length; i++)
                        this.columnNames[i] = "Column " + i;
                }

                return this.columnNames;
            }
            set
            {
                this.columnNames = value;
            }
        }

        /// <summary>
        ///   Gets the mean subtracted data.
        /// </summary>
        /// 
        public double[,] DeviationScores
        {
            get
            {
                if (this.dScores == null)
                {
                    if (this.sourceMatrix != null)
                        this.dScores = this.sourceMatrix.Center(this.Means, inPlace: false);
                    else this.dScores = Matrix.ToMatrix(this.sourceArray.Center(this.Means, inPlace: false));
                }

                return this.dScores;
            }
        }

        /// <summary>
        /// Gets the mean subtracted and deviation divided data. Also known as Z-Scores.
        /// </summary>
        /// 
        public double[,] StandardScores
        {
            get
            {
                if (this.zScores == null)
                    this.zScores = global::FileFormat.Accord.Statistics.Measures.Tools.Standardize(this.DeviationScores, this.StandardDeviations, inPlace: false);
                return this.zScores;
            }
        }

        /// <summary>
        ///   Gets the Covariance Matrix
        /// </summary>
        /// 
        public double[,] CovarianceMatrix
        {
            get
            {
                if (this.covarianceMatrix == null)
                {
                    if (this.sourceMatrix != null)
                        this.covarianceMatrix = this.sourceMatrix.Covariance(this.Means);
                    else this.covarianceMatrix = Matrix.ToMatrix(this.sourceArray.Covariance(this.Means));
                }

                return this.covarianceMatrix;
            }
        }

        /// <summary>
        ///   Gets the Correlation Matrix
        /// </summary>
        /// 
        public double[,] CorrelationMatrix
        {
            get
            {
                if (this.correlationMatrix == null)
                {
                    if (this.sourceMatrix != null)
                        this.correlationMatrix = this.sourceMatrix.Correlation(this.Means, this.StandardDeviations);
                    else this.correlationMatrix = Matrix.ToMatrix(this.sourceArray.Correlation(this.Means, this.StandardDeviations));
                }

                return this.correlationMatrix;
            }
        }

        /// <summary>
        ///   Gets a vector containing the Mean of each data column.
        /// </summary>
        /// 
        public double[] Means
        {
            get
            {
                if (this.means == null)
                {
                    if (this.sourceMatrix != null)
                        this.means = this.sourceMatrix.Mean(this.Sums);
                    else this.means = this.sourceArray.Mean(this.Sums);
                }
                return this.means;
            }
        }

        /// <summary>
        ///   Gets a vector containing the Standard Deviation of each data column.
        /// </summary>
        /// 
        public double[] StandardDeviations
        {
            get
            {
                if (this.standardDeviations == null)
                {
                    if (this.sourceMatrix != null)
                        this.standardDeviations = this.sourceMatrix.StandardDeviation(this.Means);
                    else this.standardDeviations = this.sourceArray.StandardDeviation(this.Means);
                }

                return this.standardDeviations;
            }
        }

        /// <summary>
        ///   Gets a vector containing the Standard Error of the Mean of each data column.
        /// </summary>
        /// 
        public double[] StandardErrors
        {
            get
            {
                if (this.standardErrors == null)
                    this.standardErrors = global::FileFormat.Accord.Math.Accord.Statistics.Measures.StandardError(this.samples, this.StandardDeviations);

                return this.standardErrors;
            }
        }

        /// <summary>
        ///   Gets the 95% confidence intervals for the <see cref="Means"/>.
        /// </summary>
        /// 
        public DoubleRange[] Confidence
        {
            get
            {
                if (this.confidence == null)
                {
                    this.confidence = new DoubleRange[this.variables];
                    for (int i = 0; i < this.confidence.Length; i++)
                        this.confidence[i] = this.GetConfidenceInterval(i);
                }

                return this.confidence;
            }
        }

        /// <summary>
        ///   Gets the 95% deviance intervals for the <see cref="Means"/>.
        /// </summary>
        /// 
        /// <remarks>
        ///   A deviance interval uses the standard deviation rather 
        ///   than the standard error to compute the range interval 
        ///   for a variable.
        /// </remarks>
        /// 
        public DoubleRange[] Deviance
        {
            get
            {
                if (this.deviance == null)
                {
                    this.deviance = new DoubleRange[this.variables];
                    for (int i = 0; i < this.deviance.Length; i++)
                        this.deviance[i] = this.GetDevianceInterval(i);
                }

                return this.deviance;
            }
        }

        /// <summary>
        ///   Gets a vector containing the Mode of each data column.
        /// </summary>
        /// 
        public double[] Modes
        {
            get
            {
                if (this.modes == null)
                {
                    if (this.sourceMatrix != null)
                        this.modes = this.sourceMatrix.Mode();
                    else this.modes = this.sourceArray.Mode();
                }

                return this.modes;
            }
        }

        /// <summary>
        ///   Gets a vector containing the Median of each data column.
        /// </summary>
        /// 
        public double[] Medians
        {
            get
            {
                if (this.medians == null)
                {
                    if (this.sourceMatrix != null)
                        this.medians = this.sourceMatrix.Median(type: this.quantileMethod);
                    else this.medians = this.sourceArray.Median(type: this.quantileMethod);
                }

                return this.medians;
            }
        }

        /// <summary>
        ///   Gets a vector containing the Variance of each data column.
        /// </summary>
        /// 
        public double[] Variances
        {
            get
            {
                if (this.variances == null)
                {
                    if (this.sourceMatrix != null)
                        this.variances = this.sourceMatrix.Variance(this.Means);
                    else this.variances = this.sourceArray.Variance(this.Means);
                }

                return this.variances;
            }
        }

        /// <summary>
        ///   Gets a vector containing the number of distinct elements for each data column.
        /// </summary>
        /// 
        public int[] Distinct
        {
            get
            {
                if (this.distinct == null)
                {
                    if (this.sourceMatrix != null)
                        this.distinct = Matrix.DistinctCount(this.sourceMatrix);
                    else this.distinct = Matrix.DistinctCount(this.sourceArray);
                }

                return this.distinct;
            }
        }

        /// <summary>
        ///   Gets an array containing the Ranges of each data column.
        /// </summary>
        /// 
        public DoubleRange[] Ranges
        {
            get
            {
                if (this.ranges == null)
                {
                    if (this.sourceMatrix != null)
                        this.ranges = Matrix.GetRange(this.sourceMatrix, 0);
                    else this.ranges = Matrix.GetRange(this.sourceArray, 0);
                }

                return this.ranges;
            }
        }

        /// <summary>
        ///   Gets an array containing the interquartile range of each data column.
        /// </summary>
        /// 
        public DoubleRange[] Quartiles
        {
            get
            {
                if (this.quartiles == null)
                {
                    if (this.sourceMatrix != null)
                        this.medians = this.sourceMatrix.Quartiles(out this.quartiles, type: this.quantileMethod);
                    else this.medians = this.sourceArray.Quartiles(out this.quartiles, type: this.quantileMethod);
                }

                return this.quartiles;
            }
        }

        /// <summary>
        ///   Gets an array containing the inner fences of each data column.
        /// </summary>
        /// 
        public DoubleRange[] InnerFences
        {
            get
            {
                if (this.innerFences == null)
                {
                    DoubleRange[] Q = this.Quartiles;

                    this.innerFences = new DoubleRange[this.variables];
                    for (int i = 0; i < this.innerFences.Length; i++)
                        this.innerFences[i] = global::FileFormat.Accord.Statistics.Measures.Tools.InnerFence(Q[i]);
                }

                return this.innerFences;
            }
        }

        /// <summary>
        ///   Gets an array containing the outer fences of each data column.
        /// </summary>
        /// 
        public DoubleRange[] OuterFences
        {
            get
            {
                if (this.outerFences == null)
                {
                    DoubleRange[] Q = this.Quartiles;

                    this.outerFences = new DoubleRange[this.variables];
                    for (int i = 0; i < this.outerFences.Length; i++)
                        this.outerFences[i] = global::FileFormat.Accord.Statistics.Measures.Tools.OuterFence(Q[i]);
                }

                return this.outerFences;
            }
        }

        /// <summary>
        ///   Gets an array containing the sum of each data column.
        /// </summary>
        /// 
        public double[] Sums
        {
            get
            {
                if (this.sums == null)
                {
                    if (this.sourceMatrix != null)
                        this.sums = Matrix.Sum(this.sourceMatrix, 0);
                    else this.sums = Matrix.Sum(this.sourceArray, 0);
                }

                return this.sums;
            }
        }

        /// <summary>
        /// Gets an array containing the skewness for of each data column.
        /// </summary>
        /// 
        public double[] Skewness
        {
            get
            {
                if (this.skewness == null)
                {
                    if (this.sourceMatrix != null)
                        this.skewness = this.sourceMatrix.Skewness();
                    else this.skewness = this.sourceArray.Skewness();
                }

                return this.skewness;
            }
        }

        /// <summary>
        /// Gets an array containing the kurtosis for of each data column.
        /// </summary>
        /// 
        public double[] Kurtosis
        {
            get
            {
                if (this.kurtosis == null)
                {
                    if (this.sourceMatrix != null)
                        this.kurtosis = this.sourceMatrix.Kurtosis();
                    else this.kurtosis = this.sourceArray.Kurtosis();
                }

                return this.kurtosis;
            }
        }

        /// <summary>
        ///   Gets the number of samples (or observations) in the data.
        /// </summary>
        /// 
        public int Samples
        {
            get { return this.samples; }
        }

        /// <summary>
        ///   Gets the number of variables (or features) in the data.
        /// </summary>
        /// 
        public int Variables
        {
            get { return this.variables; }
        }

        /// <summary>
        /// Gets a collection of DescriptiveMeasures objects that can be bound to a DataGridView.
        /// </summary>
        /// 
        public DescriptiveMeasureCollection Measures
        {
            get { return this.measuresCollection; }
        }


        /// <summary>
        ///   Gets a confidence interval for the <see cref="Means"/>
        ///   within the given confidence level percentage.
        /// </summary>
        /// 
        /// <param name="percent">The confidence level. Default is 0.95.</param>
        /// <param name="index">The index of the data column whose confidence
        ///   interval should be calculated.</param>
        /// 
        /// <returns>A confidence interval for the estimated value.</returns>
        /// 
        public DoubleRange GetConfidenceInterval(int index, double percent = 0.95)
        {
            double z = NormalDistribution.Standard
                .InverseDistributionFunction(0.5 + percent / 2.0);

            return new DoubleRange(
                this.Means[index] - z * this.StandardErrors[index],
                this.Means[index] + z * this.StandardErrors[index]);
        }

        /// <summary>
        ///   Gets a deviance interval for the <see cref="Means"/>
        ///   within the given confidence level percentage (i.e. uses
        ///   the standard deviation rather than the standard error to
        ///   compute the range interval for the variable).
        /// </summary>
        /// 
        /// <param name="percent">The confidence level. Default is 0.95.</param>
        /// <param name="index">The index of the data column whose confidence
        ///   interval should be calculated.</param>
        /// 
        /// <returns>A confidence interval for the estimated value.</returns>
        /// 
        public DoubleRange GetDevianceInterval(int index, double percent = 0.95)
        {
            double z = NormalDistribution.Standard
                .InverseDistributionFunction(0.5 + percent / 2.0);

            return new DoubleRange(
                this.Means[index] - z * this.StandardDeviations[index],
                this.Means[index] + z * this.StandardDeviations[index]);
        }

    }

    /// <summary>
    ///   Descriptive measures for a variable.
    /// </summary>
    /// 
    /// <seealso cref="DescriptiveAnalysis"/>
    /// 
    [Serializable]
    public class DescriptiveMeasures : IDescriptiveMeasures
    {

        private DescriptiveAnalysis analysis;
        private int index;

        internal DescriptiveMeasures(DescriptiveAnalysis analysis, int index)
        {
            this.analysis = analysis;
            this.index = index;
        }

        /// <summary>
        ///   Gets the descriptive analysis 
        ///   that originated this measure.
        /// </summary>
        /// 
        [Browsable(false)]
        public DescriptiveAnalysis Analysis
        {
            get { return this.analysis; }
        }

        /// <summary>
        ///   Gets the variable's index.
        /// </summary>
        /// 
        public int Index
        {
            get { return this.index; }
        }

        /// <summary>
        ///   Gets the variable's name
        /// </summary>
        /// 
        public string Name
        {
            get { return this.analysis.ColumnNames[this.index]; }
        }

        /// <summary>
        ///   Gets the variable's total sum.
        /// </summary>
        /// 
        public double Sum
        {
            get { return this.analysis.Sums[this.index]; }
        }

        /// <summary>
        ///   Gets the variable's mean.
        /// </summary>
        /// 
        public double Mean
        {
            get { return this.analysis.Means[this.index]; }
        }

        /// <summary>
        ///   Gets the variable's standard deviation.
        /// </summary>
        /// 
        public double StandardDeviation
        {
            get { return this.analysis.StandardDeviations[this.index]; }
        }

        /// <summary>
        ///   Gets the variable's median.
        /// </summary>
        /// 
        public double Median
        {
            get { return this.analysis.Medians[this.index]; }
        }

        /// <summary>
        ///   Gets the variable's outer fences range.
        /// </summary>
        /// 
        public DoubleRange OuterFence
        {
            get { return this.analysis.OuterFences[this.index]; }
        }

        /// <summary>
        ///   Gets the variable's inner fence range.
        /// </summary>
        /// 
        public DoubleRange InnerFence
        {
            get { return this.analysis.InnerFences[this.index]; }
        }

        /// <summary>
        ///   Gets the variable's interquartile range.
        /// </summary>
        /// 
        public DoubleRange Quartiles
        {
            get { return this.analysis.Quartiles[this.index]; }
        }

        /// <summary>
        ///   Gets the variable's mode.
        /// </summary>
        /// 
        public double Mode
        {
            get { return this.analysis.Modes[this.index]; }
        }

        /// <summary>
        ///   Gets the variable's variance.
        /// </summary>
        /// 
        public double Variance
        {
            get { return this.analysis.Variances[this.index]; }
        }

        /// <summary>
        ///   Gets the variable's skewness.
        /// </summary>
        /// 
        public double Skewness
        {
            get { return this.analysis.Skewness[this.index]; }
        }

        /// <summary>
        ///   Gets the variable's kurtosis.
        /// </summary>
        /// 
        public double Kurtosis
        {
            get { return this.analysis.Kurtosis[this.index]; }
        }

        /// <summary>
        ///   Gets the variable's standard error of the mean.
        /// </summary>
        /// 
        public double StandardError
        {
            get { return this.analysis.StandardErrors[this.index]; }
        }

        /// <summary>
        ///   Gets the variable's maximum value.
        /// </summary>
        /// 
        public double Max
        {
            get { return this.analysis.Ranges[this.index].Max; }
        }

        /// <summary>
        ///   Gets the variable's minimum value.
        /// </summary>
        /// 
        public double Min
        {
            get { return this.analysis.Ranges[this.index].Min; }
        }

        /// <summary>
        ///   Gets the variable's length.
        /// </summary>
        /// 
        public double Length
        {
            get { return this.analysis.Ranges[this.index].Length; }
        }

        /// <summary>
        ///   Gets the number of distinct values for the variable.
        /// </summary>
        /// 
        public int Distinct
        {
            get { return this.analysis.Distinct[this.index]; }
        }

        /// <summary>
        ///   Gets the number of samples for the variable.
        /// </summary>
        /// 
        public int Count
        {
            get { return this.analysis.Samples; }
        }

        /// <summary>
        ///   Gets the 95% confidence interval around the <see cref="Mean"/>.
        /// </summary>
        /// 
        public DoubleRange Confidence
        {
            get { return this.analysis.Confidence[this.index]; }
        }

        /// <summary>
        ///   Gets the 95% deviance interval around the <see cref="Mean"/>.
        /// </summary>
        /// 
        public DoubleRange Deviance
        {
            get { return this.analysis.Deviance[this.index]; }
        }

        /// <summary>
        ///   Gets the variable's observations.
        /// </summary>
        /// 
        [Obsolete("This property will be removed.")]
        public double[] Samples
        {
            get { return Matrix.GetColumn(this.analysis.Source, this.index); }
        }

        /// <summary>
        ///   Gets a confidence interval for the <see cref="Mean"/>
        ///   within the given confidence level percentage.
        /// </summary>
        /// 
        /// <param name="percent">The confidence level. Default is 0.95.</param>
        /// 
        /// <returns>A confidence interval for the estimated value.</returns>
        /// 
        public DoubleRange GetConfidenceInterval(double percent = 0.95)
        {
            return this.analysis.GetConfidenceInterval(this.index, percent);
        }

        /// <summary>
        ///   Gets a deviance interval for the <see cref="Mean"/>
        ///   within the given confidence level percentage (i.e. uses
        ///   the standard deviation rather than the standard error to
        ///   compute the range interval for the variable).
        /// </summary>
        /// 
        /// <param name="percent">The confidence level. Default is 0.95.</param>
        /// 
        /// <returns>A confidence interval for the estimated value.</returns>
        /// 
        public DoubleRange GetDevianceInterval(double percent = 0.95)
        {
            return this.analysis.GetDevianceInterval(this.index, percent);
        }

    }

    /// <summary>
    ///   Collection of descriptive measures.
    /// </summary>
    /// 
    /// <seealso cref="DescriptiveMeasures"/>
    /// <seealso cref="DescriptiveAnalysis"/>
    /// 
    [Serializable]
    public class DescriptiveMeasureCollection : ReadOnlyKeyedCollection<string, DescriptiveMeasures>
    {
        internal DescriptiveMeasureCollection(DescriptiveMeasures[] components)
            : base(components)
        {
        }

        /// <summary>
        ///   Gets the key for item.
        /// </summary>
        /// 
        protected override string GetKeyForItem(DescriptiveMeasures item)
        {
            return item.Name;
        }
    }

}
