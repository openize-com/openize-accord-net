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

namespace Openize.Accord.Statistics.Analysis
{
    using System;
    using System.ComponentModel;
    using Accord.MachineLearning.Learning;
    using Base;
    using Distributions.Univariate.Continuous;
    using Openize.Accord.Math.Matrix;
    using Openize.Accord.Statistics;
    using Measures;
    using Openize.Accord.Core.Collections;
    using Openize.Accord.Core.Exceptions;
    using Openize.Accord.Core.Ranges;
    using Openize.Accord.Math.Accord.Statistics;

    /// <summary>
    ///   Descriptive statistics analysis for circular data.
    /// </summary>
    /// 
    /// 
    /// <seealso cref="Tools"/>
    /// <seealso cref="DescriptiveAnalysis"/>
    /// <seealso cref="DescriptiveMeasures"/>
    ///
    [Serializable]
#pragma warning disable 612, 618
    public class CircularDescriptiveAnalysis : IMultivariateAnalysis,
        IDescriptiveLearning<CircularDescriptiveAnalysis, double[]>
#pragma warning restore 612, 618
    {

        private int samples;
        private int variables;

        private double[][] angles;
        private double[] lengths;

        private double[] sums;
        private double[] sin;
        private double[] cos;

        private double[] means;
        private double[] standardDeviations;
        private double[] standardErrors;
        private double[] variances;
        private double[] medians;
        private double[] modes;
        private double[] kurtosis;
        private double[] skewness;
        private int[] distinct;

        private double[] angularMeans;
        private double[] angularMedians;

        private string[] columnNames;

        private QuantileMethod quantileMethod = QuantileMethod.Default;
        private DoubleRange[] ranges;
        private DoubleRange[] quartiles;
        private DoubleRange[] innerFences;
        private DoubleRange[] outerFences;

        DoubleRange[] confidences;
        DoubleRange[] deviances;

        private double[] concentration;


        private double[,] sourceMatrix;
        private double[][] sourceArray;
        private double[] sourceRow;

        private CircularDescriptiveMeasureCollection measuresCollection;

        private bool useStrictRanges = true;
        private bool lazy = true;

        /// <summary>
        ///   Constructs the Circular Descriptive Analysis.
        /// </summary>
        /// 
        /// <param name="data">The source data to perform analysis.</param>
        /// <param name="length">The length of each circular variable (i.e. 24 for hours).</param>
        /// <param name="columnName">The names for the analyzed variable.</param>
        /// <param name="inPlace">
        ///   Whether the analysis should conserve memory by doing 
        ///   operations over the original <paramref name="data"/> array.
        /// </param>
        /// 
        [Obsolete("Please pass only the lengths and columnNames parameters and call the Learn() method passing the data to be analyzed.")]
        public CircularDescriptiveAnalysis(double[] data, double length, String columnName, bool inPlace = false)
        {
            if (data == null)
                throw new ArgumentNullException("data");

            this.compute(data, null, null, new[] { length }, new[] { columnName }, inPlace);
        }

        /// <summary>
        ///   Constructs the Circular Descriptive Analysis.
        /// </summary>
        /// 
        /// <param name="data">The source data to perform analysis.</param>
        /// <param name="length">The length of each circular variable (i.e. 24 for hours).</param>
        /// <param name="inPlace">
        ///   Whether the analysis should conserve memory by doing 
        ///   operations over the original <paramref name="data"/> array.
        /// </param>
        /// 
        [Obsolete("Please pass only the lengths and columnNames parameters and call the Learn() method passing the data to be analyzed.")]
        public CircularDescriptiveAnalysis(double[] data, double length, bool inPlace = false)
        {
            if (data == null)
                throw new ArgumentNullException("data");

            this.compute(data, null, null, new[] { length }, null, inPlace);
        }

        /// <summary>
        ///   Constructs the Circular Descriptive Analysis.
        /// </summary>
        /// 
        /// <param name="data">The source data to perform analysis.</param>
        /// <param name="length">The length of each circular variable (i.e. 24 for hours).</param>
        /// 
        [Obsolete("Please pass only the lengths and columnNames parameters and call the Learn() method passing the data to be analyzed.")]
        public CircularDescriptiveAnalysis(double[,] data, double[] length)
        {
            if (data == null)
                throw new ArgumentNullException("data");

            if (length == null)
                throw new ArgumentNullException("length");

            this.compute(null, data, null, length, null);
        }

        /// <summary>
        ///   Constructs the Circular Descriptive Analysis.
        /// </summary>
        /// 
        /// <param name="data">The source data to perform analysis.</param>
        /// <param name="length">The length of each circular variable (i.e. 24 for hours).</param>
        /// <param name="columnNames">Names for the analyzed variables.</param>
        /// 
        [Obsolete("Please pass only the lengths and columnNames parameters and call the Learn() method passing the data to be analyzed.")]
        public CircularDescriptiveAnalysis(double[,] data, double[] length, string[] columnNames)
        {
            if (data == null)
                throw new ArgumentNullException("data");

            if (length == null)
                throw new ArgumentNullException("length");

            if (columnNames == null)
                throw new ArgumentNullException("columnNames");

            this.compute(null, data, null, length, columnNames);
        }

        /// <summary>
        ///   Constructs the Circular Descriptive Analysis.
        /// </summary>
        /// 
        /// <param name="data">The source data to perform analysis.</param>
        /// <param name="length">The length of each circular variable (i.e. 24 for hours).</param>
        /// 
        [Obsolete("Please pass only the lengths and columnNames parameters and call the Learn() method passing the data to be analyzed.")]
        public CircularDescriptiveAnalysis(double[][] data, double[] length)
        {
            // Initial argument checking
            if (data == null)
                throw new ArgumentNullException("data");

            if (length == null)
                throw new ArgumentNullException("length");

            this.compute(null, null, data, length, null);
        }

        /// <summary>
        ///   Constructs the Circular Descriptive Analysis.
        /// </summary>
        /// 
        /// <param name="data">The source data to perform analysis.</param>
        /// <param name="length">The length of each circular variable (i.e. 24 for hours).</param>
        /// <param name="columnNames">Names for the analyzed variables.</param>
        /// 
        [Obsolete("Please pass only the lengths and columnNames parameters and call the Learn() method passing the data to be analyzed.")]
        public CircularDescriptiveAnalysis(double[][] data, double[] length, string[] columnNames)
        {
            // Initial argument checking
            if (data == null)
                throw new ArgumentNullException("data");

            if (length == null)
                throw new ArgumentNullException("length");

            if (columnNames == null)
                throw new ArgumentNullException("columnNames");

            this.compute(null, null, data, length, columnNames);
        }

        /// <summary>
        ///   Constructs the Circular Descriptive Analysis.
        /// </summary>
        /// 
        /// <param name="length">The length of each circular variable (i.e. 24 for hours).</param>
        /// <param name="columnNames">Names for the analyzed variables.</param>
        /// 
        public CircularDescriptiveAnalysis(double[] length, string[] columnNames)
        {
            this.compute(null, null, null, length, columnNames);
        }

        /// <summary>
        ///   Constructs the Circular Descriptive Analysis.
        /// </summary>
        /// 
        /// <param name="length">The length of each circular variable (i.e. 24 for hours).</param>
        /// 
        public CircularDescriptiveAnalysis(double[] length)
        {
            this.compute(null, null, null, length, null);
        }


        private void compute(double[] row, double[,] matrix, double[][] array,
            double[] length, string[] columnNames, bool inPlace = false)
        {
            this.columnNames = columnNames;
            this.sourceArray = array;
            this.sourceMatrix = matrix;
            this.sourceRow = row;

            if (matrix != null)
            {
                this.samples = Matrix.Rows(matrix);
                this.variables = Matrix.Columns(matrix);
                if (this.lengths == null)
                    this.lengths = Matrix.Max(matrix, (int)0);
                if (this.lengths.Length != this.variables)
                    throw new DimensionMismatchException("length");
                this.lengths = length;
                this.angles = new double[this.variables][];
                for (int i = 0; i < this.angles.Length; i++)
                {
                    this.angles[i] = new double[this.samples];
                    for (int j = 0; j < this.angles[i].Length; j++)
                        this.angles[i][j] = Circular.ToRadians(matrix[j, i], this.lengths[i]);
                }
            }
            else if (array != null)
            {
                this.samples = array.Length;
                this.variables = array[0].Length;
                this.angles = new double[this.variables][];
                if (this.lengths == null)
                    this.lengths = Matrix.Max(array, (int)1);
                if (this.lengths.Length != this.variables)
                    throw new DimensionMismatchException("length");
                this.lengths = length;
                for (int i = 0; i < this.angles.Length; i++)
                {
                    this.angles[i] = new double[this.samples];
                    for (int j = 0; j < this.angles[i].Length; j++)
                        this.angles[i][j] = Circular.ToRadians(array[j][i], length[i]);
                }
            }
            else
            {
                this.samples = this.sourceRow.Length;
                this.variables = 1;
                this.angles = new double[this.variables][];
                if (this.lengths == null)
                    this.lengths = new[] { Matrix.Max(this.sourceRow) };
                if (this.lengths.Length != this.variables)
                    throw new DimensionMismatchException("length");
                this.lengths = length;
                this.angles[0] = inPlace ? this.sourceRow : new double[this.samples];
                for (int j = 0; j < this.angles[0].Length; j++)
                    this.angles[0][j] = Circular.ToRadians(this.sourceRow[j], length[0]);
            }


            // Create object-oriented structure to access data
            var measures = new CircularDescriptiveMeasures[this.variables];
            for (int i = 0; i < measures.Length; i++)
                measures[i] = new CircularDescriptiveMeasures(this, i);
            this.measuresCollection = new CircularDescriptiveMeasureCollection(measures);
        }


        /// <summary>
        ///   Computes the analysis using given source data and parameters.
        /// </summary>
        /// 
        [Obsolete("Please use Learn() instead.")]
        public void Compute()
        {
            this.reset();

            this.sums = this.Sums;
            this.means = this.Means;
            this.standardDeviations = this.StandardDeviations;
            this.ranges = this.Ranges;
            this.medians = this.Medians;
            this.variances = this.Variances;
            this.distinct = this.Distinct;
            this.quartiles = this.Quartiles;
            this.innerFences = this.InnerFences;
            this.outerFences = this.OuterFences;
            this.modes = this.Modes;
            this.cos = this.CosineSum;
            this.sin = this.SineSum;
            this.skewness = this.Skewness;
            this.kurtosis = this.Kurtosis;
            this.concentration = this.Concentration;
            this.deviances = this.Deviance;
            this.confidences = this.Confidence;
        }

        private void reset()
        {
            this.sums = null;
            this.means = null;
            this.standardDeviations = null;
            this.ranges = null;
            this.medians = null;
            this.variances = null;
            this.standardErrors = null;
            this.distinct = null;
            this.modes = null;
            this.deviances = null;
            this.confidences = null;
            this.quartiles = null;
            this.innerFences = null;
            this.outerFences = null;
            this.sin = null;
            this.cos = null;
            this.skewness = null;
            this.kurtosis = null;
            this.concentration = null;
            this.standardErrors = null;
        }

        /// <summary>
        /// Learns a model that can map the given inputs to the desired outputs.
        /// </summary>
        /// <param name="x">The model inputs.</param>
        /// <returns>
        /// A model that has learned how to produce suitable outputs
        /// given the input data <paramref name="x" />.
        /// </returns>
        public CircularDescriptiveAnalysis Learn(double[][] x)
        {
            this.reset();

            this.compute(null, null, x, this.lengths, this.columnNames);

            if (!this.lazy)
            {
                this.sums = this.Sums;
                this.means = this.Means;
                this.standardDeviations = this.StandardDeviations;
                this.ranges = this.Ranges;
                this.medians = this.Medians;
                this.variances = this.Variances;
                this.distinct = this.Distinct;
                this.quartiles = this.Quartiles;
                this.innerFences = this.InnerFences;
                this.outerFences = this.OuterFences;
                this.modes = this.Modes;
                this.cos = this.CosineSum;
                this.sin = this.SineSum;
                this.skewness = this.Skewness;
                this.kurtosis = this.Kurtosis;
                this.concentration = this.Concentration;
                this.deviances = this.Deviance;
                this.confidences = this.Confidence;
                this.sourceArray = null;
                this.sourceMatrix = null;
                this.sourceRow = null;
            }

            return this;
        }

        /// <summary>
        ///   Gets or sets whether all reported statistics should respect the circular 
        ///   interval. For example, setting this property to <c>false</c> would allow
        ///   the <see cref="Confidence"/>, <see cref="Deviance"/>, <see cref="InnerFences"/>
        ///   and <see cref="OuterFences"/> properties report minimum and maximum values 
        ///   outside the variable's allowed circular range. Default is <c>true</c>.
        /// </summary>
        /// 
        public bool UseStrictRanges
        {
            get { return this.useStrictRanges; }
            set
            {
                if (this.useStrictRanges != value)
                {
                    this.useStrictRanges = value;
                    this.innerFences = null;
                    this.outerFences = null;
                    this.deviances = null;
                    this.confidences = null;
                }
            }
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
                {
                    if (this.sourceArray == null)
                    {
                        this.sourceMatrix = Matrix.ToMatrix(this.sourceRow, (bool)true);
                    }
                    else
                    {
                        this.sourceMatrix = Matrix.ToMatrix(this.sourceArray);
                    }
                }

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
                {
                    if (this.sourceRow == null)
                    {
                        this.sourceArray = Matrix.ToJagged(this.sourceMatrix);
                    }
                    else
                    {
                        this.sourceArray = Matrix.ToJagged(this.sourceRow, (bool)true);
                    }
                }

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
        public double[][] Angles
        {
            get { return this.angles; }
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
        }

        /// <summary>
        ///   Gets a vector containing the length of
        ///   the circular domain for each data column.
        /// </summary>
        /// 
        public double[] Lengths
        {
            get { return this.lengths; }
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
                    double[] c = this.CosineSum, s = this.SineSum;

                    this.means = new double[this.variables];
                    this.angularMeans = new double[this.variables];
                    for (int i = 0; i < this.means.Length; i++)
                    {
                        this.angularMeans[i] = Circular.Mean(this.samples, c[i], s[i]);
                        this.means[i] = Circular.ToCircular(this.angularMeans[i], this.lengths[i]);
                    }
                }

                return this.means;
            }
            private set { this.means = value; }
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
                    else if (this.sourceArray != null)
                        this.modes = this.sourceArray.Mode();
                    else this.modes = new[] { this.sourceRow.Mode() };
                }

                return this.modes;
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
                    double[] c = this.CosineSum, s = this.SineSum;

                    this.standardDeviations = new double[this.variables];
                    for (int i = 0; i < this.standardDeviations.Length; i++)
                    {
                        this.standardDeviations[i] = Circular.StandardDeviation(this.samples, c[i], s[i])
                            * this.lengths[i] / (2 * Math.PI);
                    }
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
                {
                    double[] c = this.CosineSum, s = this.SineSum;

                    this.standardErrors = new double[this.variables];
                    for (int i = 0; i < this.standardErrors.Length; i++)
                    {
                        this.standardErrors[i] = Circular.StandardError(this.samples, c[i], s[i], 0.05)
                            * this.lengths[i] / (2 * Math.PI);
                    }
                }

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
                if (this.confidences == null)
                {
                    this.confidences = new DoubleRange[this.variables];
                    for (int i = 0; i < this.confidences.Length; i++)
                        this.confidences[i] = this.GetConfidenceInterval(i);
                }

                return this.confidences;
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
                if (this.deviances == null)
                {
                    this.deviances = new DoubleRange[this.variables];
                    for (int i = 0; i < this.deviances.Length; i++)
                        this.deviances[i] = this.GetDevianceInterval(i);
                }

                return this.deviances;
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
                    this.computeMedians();

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
                    double[] c = this.CosineSum, s = this.SineSum;

                    this.variances = new double[this.variables];
                    for (int i = 0; i < this.variances.Length; i++)
                    {
                        double scale = this.lengths[i] / (2 * Math.PI);
                        this.variances[i] = Circular.Variance(this.samples, c[i], s[i]) * scale * scale;
                    }
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
                    else if (this.sourceArray != null)
                        this.distinct = Matrix.DistinctCount(this.sourceArray);
                    else
                        this.distinct = new[] { Matrix.DistinctCount(this.sourceRow) };
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
                    this.ranges = new DoubleRange[this.variables];
                    for (int i = 0; i < this.ranges.Length; i++)
                        this.ranges[i] = new DoubleRange(0, this.lengths[i]);
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
                    if (this.medians == null)
                        this.computeMedians();

                    this.quartiles = new DoubleRange[this.variables];
                    for (int i = 0; i < this.variances.Length; i++)
                    {
                        Circular.Quartiles(this.angles[i], out this.quartiles[i], this.angularMedians[i], wrap: this.useStrictRanges, type: this.quantileMethod);

                        this.quartiles[i].Min = Circular.ToCircular(this.quartiles[i].Min, this.lengths[i], this.useStrictRanges);
                        this.quartiles[i].Max = Circular.ToCircular(this.quartiles[i].Max, this.lengths[i], this.useStrictRanges);
                    }
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
                    {
                        this.innerFences[i] = global::Openize.Accord.Statistics.Measures.Tools.InnerFence(Q[i]);

                        if (this.useStrictRanges)
                        {
                            this.innerFences[i].Min = global::Openize.Accord.Math.Tools.Mod(this.innerFences[i].Min, this.lengths[i]);
                            this.innerFences[i].Max = global::Openize.Accord.Math.Tools.Mod(this.innerFences[i].Max, this.lengths[i]);
                        }
                    }
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
                    {
                        this.outerFences[i] = global::Openize.Accord.Statistics.Measures.Tools.OuterFence(Q[i]);

                        if (this.useStrictRanges)
                        {
                            this.outerFences[i].Min = global::Openize.Accord.Math.Tools.Mod(this.outerFences[i].Min, this.lengths[i]);
                            this.outerFences[i].Max = global::Openize.Accord.Math.Tools.Mod(this.outerFences[i].Max, this.lengths[i]);
                        }
                    }
                }

                return this.outerFences;
            }
        }

        /// <summary>
        ///   Gets an array containing the sum of each data column. If 
        ///   the analysis has been computed in place, this will contain 
        ///   the sum of the transformed angle values instead.
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
                    else if (this.sourceArray != null)
                        this.sums = Matrix.Sum(this.sourceArray, 0);
                    else this.sums = new[] { Matrix.Sum(this.sourceRow) };
                }

                return this.sums;
            }
        }

        /// <summary>
        ///   Gets an array containing the sum of cosines for each data column.
        /// </summary>
        /// 
        public double[] CosineSum
        {
            get
            {
                if (this.cos == null)
                    this.computeSums();

                return this.cos;
            }
        }

        /// <summary>
        ///   Gets an array containing the sum of sines for each data column.
        /// </summary>
        /// 
        public double[] SineSum
        {
            get
            {
                if (this.sin == null)
                    this.computeSums();

                return this.sin;
            }
        }

        /// <summary>
        ///   Gets an array containing the circular concentration for each data column.
        /// </summary>
        /// 
        public double[] Concentration
        {
            get
            {
                if (this.concentration == null)
                {
                    var m = this.Means;

                    this.concentration = new double[this.variables];
                    for (int i = 0; i < this.variances.Length; i++)
                        this.concentration[i] = Circular.Concentration(this.angles[i], m[i]);
                }

                return this.variances;
            }
        }

        /// <summary>
        ///   Gets an array containing the skewness for of each data column.
        /// </summary>
        /// 
        public double[] Skewness
        {
            get
            {
                if (this.skewness == null)
                {
                    this.skewness = new double[this.variables];
                    for (int i = 0; i < this.skewness.Length; i++)
                    {
                        this.skewness[i] = Circular.Skewness(this.angles[i])
                            * this.lengths[i] / (2 * Math.PI);
                    }
                }

                return this.skewness;
            }
        }

        /// <summary>
        ///   Gets an array containing the kurtosis for of each data column.
        /// </summary>
        /// 
        public double[] Kurtosis
        {
            get
            {
                if (this.kurtosis == null)
                {
                    this.kurtosis = new double[this.variables];
                    for (int i = 0; i < this.kurtosis.Length; i++)
                    {
                        this.kurtosis[i] = Circular.Kurtosis(this.angles[i])
                            * this.lengths[i] / (2 * Math.PI);
                    }
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
        public CircularDescriptiveMeasureCollection Measures
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
            double[] c = this.CosineSum, s = this.SineSum;

            double t = Circular.StandardError(this.samples, c[index], s[index], 1 - percent);

            t *= this.lengths[index] / (2 * Math.PI);

            double min = this.Means[index] - t;
            double max = this.Means[index] + t;

            if (this.useStrictRanges)
            {
                min = global::Openize.Accord.Math.Tools.Mod(min, this.lengths[index]);
                max = global::Openize.Accord.Math.Tools.Mod(max, this.lengths[index]);
            }

            return new DoubleRange(min, max);
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

            double min = this.Means[index] - z * this.StandardDeviations[index];
            double max = this.Means[index] + z * this.StandardDeviations[index];

            if (this.useStrictRanges)
            {
                min = global::Openize.Accord.Math.Tools.Mod(min, this.lengths[index]);
                max = global::Openize.Accord.Math.Tools.Mod(max, this.lengths[index]);
            }

            return new DoubleRange(min, max);
        }



        private void computeSums()
        {
            this.cos = new double[this.variables];
            this.sin = new double[this.variables];
            for (int i = 0; i < this.angles.Length; i++)
                Circular.Sum(this.angles[i], out this.cos[i], out this.sin[i]);
        }

        private void computeMedians()
        {
            this.medians = new double[this.variables];
            this.angularMedians = new double[this.variables];
            for (int i = 0; i < this.medians.Length; i++)
            {
                this.angularMedians[i] = Circular.Median(this.angles[i]);
                this.medians[i] = Circular.ToCircular(this.angularMedians[i], this.lengths[i]);
            }
        }

    }

    /// <summary>
    ///   Circular descriptive measures for a variable.
    /// </summary>
    /// 
    /// <seealso cref="CircularDescriptiveAnalysis"/>
    /// 
    [Serializable]
    public class CircularDescriptiveMeasures : IDescriptiveMeasures
    {

        private CircularDescriptiveAnalysis analysis;
        private int index;

        internal CircularDescriptiveMeasures(CircularDescriptiveAnalysis analysis, int index)
        {
            this.analysis = analysis;
            this.index = index;
        }

        /// <summary>
        ///   Gets the circular analysis 
        ///   that originated this measure.
        /// </summary>
        /// 
        [Browsable(false)]
        public CircularDescriptiveAnalysis Analysis
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
        ///   Gets the variable's mode.
        /// </summary>
        /// 
        public double Mode
        {
            get { return this.analysis.Modes[this.index]; }
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
        ///   Gets the variable's variance.
        /// </summary>
        /// 
        public double Variance
        {
            get { return this.analysis.Variances[this.index]; }
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
        ///   Gets the sum of cosines for the variable.
        /// </summary>
        /// 
        public double CosineSum
        {
            get { return this.analysis.CosineSum[this.index]; }
        }

        /// <summary>
        ///   Gets the sum of sines for the variable.
        /// </summary>
        /// 
        public double SineSum
        {
            get { return this.analysis.SineSum[this.index]; }
        }

        /// <summary>
        ///   Gets the transformed variable's observations.
        /// </summary>
        /// 
        [Obsolete("This property will be removed.")]
        public double[] Angles
        {
            get { return this.analysis.Angles[this.index]; }
        }

        /// <summary>
        ///   Gets the variable's standard error of the mean.
        /// </summary>
        /// 
        public double StandardError
        {
            get { return this.analysis.StandardErrors[this.index]; ; }
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
        ///   Gets the variable <see cref="Circular.Skewness">skewness</see>.
        /// </summary>
        /// 
        public double Skewness
        {
            get { return this.analysis.Skewness[this.index]; }
        }

        /// <summary>
        ///   Gets the variable <see cref="Circular.Skewness">kurtosis</see>.
        /// </summary>
        /// 
        public double Kurtosis
        {
            get { return this.analysis.Kurtosis[this.index]; }
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
    /// <seealso cref="CircularDescriptiveMeasures"/>
    /// <seealso cref="CircularDescriptiveAnalysis"/>
    /// 
    [Serializable]
    public class CircularDescriptiveMeasureCollection : ReadOnlyKeyedCollection<string, CircularDescriptiveMeasures>
    {
        internal CircularDescriptiveMeasureCollection(CircularDescriptiveMeasures[] components)
            : base(components)
        {
        }

        /// <summary>
        ///   Gets the key for item.
        /// </summary>
        /// 
        protected override string GetKeyForItem(CircularDescriptiveMeasures item)
        {
            return item.Name;
        }
    }
}
