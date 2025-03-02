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

namespace Openize.Accord.Statistics.Visualizations
{
    using System;
    using Openize.Accord.Math.Matrix;
    using Openize.Accord.Core.Ranges;
    using Openize.Accord.Math.Accord.Statistics;
    using Vector = Openize.Accord.Math.Vector.Vector;

    /// <summary>
    ///   Optimum histogram bin size adjustment rule.
    /// </summary>
    /// 
    public enum BinAdjustmentRule
    {
        /// <summary>
        ///   Does not attempts to automatically calculate 
        ///   an optimum bin width and preserves the current
        ///   histogram organization.
        /// </summary>
        /// 
        None,

        /// <summary>
        ///   Calculates the optimum bin width as 3.49σN, where σ 
        ///   is the sample standard deviation and N is the number
        ///   of samples.
        /// </summary>
        /// <remarks>
        ///   Scott, D. 1979. On optimal and data-based histograms. Biometrika, 66:605-610.
        /// </remarks>
        /// 
        Scott,

        /// <summary>
        ///   Calculates the optimum bin width as <c>ceiling(log2(N) + 1)</c>
        ///   where N is the number of samples. The rule implicitly bases
        ///   the bin sizes on the range of the data, and can perform poorly
        ///   if n &lt; 30.
        /// </summary>
        /// 
        Sturges,

        /// <summary>
        ///   Calculates the optimum bin width as the square root of the
        ///   number of samples. This is the same rule used by Microsoft (c)
        ///   Excel and many others.
        /// </summary>
        /// 
        SquareRoot,
    }

    /// <summary>
    ///   Histogram.
    /// </summary>
    /// 
    /// <remarks>
    ///  <para>
    ///   In a more general mathematical sense, a histogram is a mapping Mi
    ///   that counts the number of observations that fall into various 
    ///   disjoint categories (known as bins).</para>
    ///  <para>
    ///   This class represents a Histogram mapping of Discrete or Continuous
    ///   data. To use it as a discrete mapping, pass a bin size (length) of 1.
    ///   To use it as a continuous mapping, pass any real number instead.</para>
    ///  <para>
    ///   Currently, only a constant bin width is supported.</para>
    /// </remarks>
    /// 
    /// <example>
    /// <code source="Unit Tests\Accord.Tests.Statistics\Visualizations\HistogramTest.cs" region="doc_example1" />
    /// </example>
    /// 
    [Serializable]
    public class Histogram : ICloneable
    {

        string title = "Histogram";

        private int[] values;
        private double[] ranges;
        private bool cumulative;
        private bool inclusiveUpperBound = true;

        private HistogramBinCollection binCollection;
        private BinAdjustmentRule rule = BinAdjustmentRule.SquareRoot;

        private bool uniform;
        private double? mean;
        private double? stdDev;
        private int? median;
        private int? min;
        private int? max;
        private long? total;


        /// <summary>
        ///   Constructs an empty histogram
        /// </summary>
        /// 
        public Histogram()
            : this(new int[0])
        {
        }



        /// <summary>
        ///   Constructs an empty histogram
        /// </summary>
        /// 
        /// <param name="values">The values to be binned in the histogram.</param>
        /// 
        [Obsolete("Please use Histogram.FromData instead.")]
        public Histogram(double[] values)
        {
            this.Compute(values);
        }

        /// <summary>
        ///   Constructs an empty histogram
        /// </summary>
        /// 
        /// <param name="values">The values to be binned in the histogram.</param>
        /// 
        public Histogram(int[] values)
        {
            this.values = values;
            this.ranges = Vector.Ones(values.Length + 1);
            this.uniform = true;
        }


        /// <summary>
        ///   Initializes the histogram's bins.
        /// </summary>
        /// 
        private void initialize(int numberOfBins)
        {
            this.values = new int[numberOfBins];
            this.ranges = new double[numberOfBins + 1];
            this.Update();
        }

        /// <summary>
        ///   Sets the histogram's bin ranges (edges).
        /// </summary>
        /// 
        private void initialize(double startValue, double width)
        {
            this.ranges[0] = startValue;
            for (int i = 1; i < this.ranges.Length; i++)
                this.ranges[i] = this.ranges[i - 1] + width;
        }

        /// <summary>
        /// Update statistical value of the histogram.
        /// </summary>
        /// 
        /// <remarks>The method recalculates statistical values of the histogram, like mean,
        /// standard deviation, etc., in the case if histogram's values were changed directly.
        /// The method should be called only in the case if histogram's values were retrieved
        /// through <see cref="Values"/> property and updated after that.
        /// </remarks>
        /// 
        public void Update()
        {
            this.mean = null;
            this.stdDev = null;
            this.median = null;
            this.min = null;
            this.max = null;
            this.total = null;
        }



        /// <summary>
        ///   Gets the Bin values of this Histogram.
        /// </summary>
        /// 
        /// <param name="index">Bin index.</param>
        /// 
        /// <returns>The number of hits of the selected bin.</returns>
        /// 
        public int this[int index]
        {
            get { return this.values[index]; }
            set { this.values[index] = value; }
        }

        /// <summary>
        /// Gets or sets the title of this histogram. Default value is "Histogram".
        /// </summary>
        /// 
        public string Title
        {
            get { return this.title; }
            set { this.title = value; }
        }

        /// <summary>
        ///   Gets the Bin values for this Histogram.
        /// </summary>
        /// 
        public int[] Values
        {
            get { return this.values; }
        }

        /// <summary>
        ///   Gets the Range of the values in this Histogram.
        /// </summary>
        /// 
        public DoubleRange Range
        {
            get { return new DoubleRange(this.ranges[0], this.ranges[this.ranges.Length - 1]); }
        }

        /// <summary>
        ///   Gets the edges of each bin in this Histogram.
        /// </summary>
        /// 
        public double[] Edges
        {
            get { return this.ranges; }
        }

        /// <summary>
        ///   Gets the collection of bins of this Histogram.
        /// </summary>
        /// 
        public HistogramBinCollection Bins
        {
            get
            {
                if (this.binCollection == null)
                {
                    var bins = new HistogramBin[this.values.Length];
                    for (int i = 0; i < bins.Length; i++)
                        bins[i] = new HistogramBin(this, i);
                    this.binCollection = new HistogramBinCollection(bins);
                }
                return this.binCollection;
            }
        }

        /// <summary>
        ///   Gets or sets whether this histogram represents a cumulative distribution.
        /// </summary>
        /// 
        public bool Cumulative
        {
            get { return this.cumulative; }
            set { this.cumulative = value; }
        }


        /// <summary>
        ///   Gets or sets the bin size auto adjustment rule
        ///   to be used when computing this histogram from
        ///   new data. Default is <see cref="BinAdjustmentRule.SquareRoot"/>.
        /// </summary>
        /// 
        /// <value>The bin size auto adjustment rule.</value>
        /// 
        public BinAdjustmentRule AutoAdjustmentRule
        {
            get { return this.rule; }
            set { this.rule = value; }
        }

        /// <summary>
        /// Mean value.
        /// </summary>
        /// 
        /// <remarks><para>The property allows to retrieve mean value of the histogram.</para>
        /// </remarks>
        /// 
        public double Mean
        {
            get
            {
                if (this.mean == null && this.uniform)
                    this.mean = this.values.HistogramMean();
                return this.mean.Value;
            }
        }

        /// <summary>
        /// Standard deviation.
        /// </summary>
        /// 
        /// <remarks><para>The property allows to retrieve standard deviation value of the histogram.</para>
        /// </remarks>
        /// 
        public double StdDev
        {
            get
            {
                if (this.stdDev == null && this.uniform)
                    this.stdDev = this.values.HistogramStandardDeviation(this.Mean);
                return this.stdDev.Value;
            }
        }

        /// <summary>
        /// Median value.
        /// </summary>
        /// 
        /// <remarks><para>The property allows to retrieve median value of the histogram.</para>
        /// </remarks>
        /// 
        public int Median
        {
            get
            {
                if (this.median == null && this.uniform)
                    this.median = this.values.HistogramMedian();
                return this.median.Value;
            }
        }

        /// <summary>
        /// Minimum value.
        /// </summary>
        /// 
        /// <remarks><para>The property allows to retrieve minimum value of the histogram with non zero
        /// hits count.</para>
        /// </remarks>
        /// 
        public int Min
        {
            get
            {
                if (this.min == null && this.uniform)
                    this.min = this.values.HistogramMin();
                return this.min.Value;
            }
        }

        /// <summary>
        /// Maximum value.
        /// </summary>
        /// 
        /// <remarks><para>The property allows to retrieve maximum value of the histogram with non zero
        /// hits count.</para>
        /// </remarks>
        /// 
        public int Max
        {
            get
            {
                if (this.max == null && this.uniform)
                    this.max = this.values.HistogramMax();
                return this.max.Value;
            }
        }

        /// <summary>
        /// Total count of values.
        /// </summary>
        /// 
        /// <remarks><para>The property represents total count of values contributed to the histogram, which is
        /// essentially sum of the <see cref="Values"/> array.</para>
        /// </remarks>
        /// 
        public long TotalCount
        {
            get
            {
                if (this.total == null && this.uniform)
                    this.total = this.values.HistogramSum();
                return this.total.Value;
            }
        }

        /// <summary>
        ///   Gets or sets a value indicating whether the last bin
        ///   should have an inclusive upper bound. Default is <c>true</c>.
        /// </summary>
        /// 
        /// <remarks>
        ///   If set to <c>false</c>, the last bin's range will be defined
        ///   as Edge[i] &lt;= x &lt; Edge[i+1]. If set to <c>true</c>, the
        ///   last bin will have an inclusive upper bound and be defined as
        ///   Edge[i] &lt;= x &lt;= Edge[i+1] instead.
        /// </remarks>
        /// 
        /// <value>
        ///   <c>true</c> if the last bin should have an inclusive upper bound;
        ///   <c>false</c> otherwise.
        /// </value>
        /// 
        public bool InclusiveUpperBound
        {
            get { return this.inclusiveUpperBound; }
            set { this.inclusiveUpperBound = value; }
        }

        /// <summary>
        /// Get range around median containing specified percentage of values.
        /// </summary>
        /// 
        /// <param name="percent">Values percentage around median.</param>
        /// 
        /// <returns>Returns the range which containes specifies percentage of values.</returns>
        /// 
        /// <remarks><para>The method calculates range of stochastic variable, which summary probability
        /// comprises the specified percentage of histogram's hits.</para>
        /// 
        /// <para>Sample usage:</para>
        /// <code>
        /// // create histogram
        /// Histogram histogram = new Histogram( new int[10] { 0, 0, 1, 3, 6, 8, 11, 0, 0, 0 } );
        /// // get 50% range
        /// IntRange range = histogram.GetRange( 0.5 );
        /// // show the range ([4, 6])
        /// Console.WriteLine( "50% range = [" + range.Min + ", " + range.Max + "]" );
        /// </code>
        /// </remarks>
        /// 
        public IntRange GetRange(double percent)
        {
            return this.values.GetHistogramRange(percent);
        }



        /// <summary>
        ///   Computes (populates) an Histogram mapping with values from a sample. 
        /// </summary>
        /// 
        /// <param name="values">The values to be binned in the histogram.</param>
        /// <param name="binWidth">The desired width for the histogram's bins.</param>
        /// 
        public void Compute(double[] values, double binWidth)
        {
            if (values == null)
                throw new ArgumentNullException("values");

            if (binWidth <= 0.0)
                throw new ArgumentOutOfRangeException("binWidth");

            // Compute values' range
            DoubleRange range = values.GetRange();

            // Determine number of bins based on the given width
            int numberOfBins = (int)Math.Ceiling(range.Length / binWidth);

            // Create bin structure
            this.initialize(numberOfBins);

            // Create ranges w/ fixed width
            this.initialize(range.Min, binWidth);

            // Create histogram
            this.compute(values);
        }

        /// <summary>
        ///   Computes (populates) an Histogram mapping with values from a sample. 
        /// </summary>
        /// 
        /// <param name="values">The values to be binned in the histogram.</param>
        /// <param name="numberOfBins">The desired number of histogram's bins.</param>
        /// 
        public void Compute(double[] values, int numberOfBins)
        {
            this.Compute(values, numberOfBins, false);
        }

        /// <summary>
        ///   Computes (populates) an Histogram mapping with values from a sample. 
        /// </summary>
        /// 
        /// <param name="values">The values to be binned in the histogram.</param>
        /// <param name="numberOfBins">The desired number of histogram's bins.</param>
        /// <param name="extraUpperBin">Whether to include an extra upper bin going to infinity.</param>
        /// 
        public void Compute(double[] values, int numberOfBins, bool extraUpperBin)
        {
            if (values == null)
                throw new ArgumentNullException("values");

            if (numberOfBins <= 0)
                throw new ArgumentOutOfRangeException("numberOfBins");

            // Compute values' range
            DoubleRange range = values.GetRange();

            // Determine bin width based on the given number of bins
            double binWidth = range.Length / (double)numberOfBins;

            if (extraUpperBin) numberOfBins++;

            // Create bin structure
            this.initialize(numberOfBins);

            // Create ranges w/ fixed width
            this.initialize(range.Min, binWidth);

            // Create histogram
            this.compute(values);
        }

        /// <summary>
        ///   Computes (populates) an Histogram mapping with values from a sample. 
        /// </summary>
        /// 
        /// <param name="values">The values to be binned in the histogram.</param>
        /// <param name="numberOfBins">The desired number of histogram's bins.</param>
        /// <param name="binWidth">The desired width for the histogram's bins.</param>
        /// 
        public void Compute(double[] values, int numberOfBins, double binWidth)
        {
            if (values == null)
                throw new ArgumentNullException("values");

            if (numberOfBins <= 0)
                throw new ArgumentOutOfRangeException("numberOfBins");

            if (binWidth <= 0.0)
                throw new ArgumentOutOfRangeException("binWidth");

            // Compute values' range
            DoubleRange range = values.GetRange();

            // Create bin structure
            this.initialize(numberOfBins);

            // Create ranges w/ fixed width
            this.initialize(range.Min, binWidth);

            // Create histogram
            this.compute(values);
        }

        /// <summary>
        ///   Computes (populates) an Histogram mapping with values from a sample. 
        /// </summary>
        /// 
        /// <param name="values">The values to be binned in the histogram.</param>
        /// 
        public void Compute(double[] values)
        {
            // Compute values' range
            DoubleRange range = values.GetRange();

            // Check if there are no values
            if (values.Length == 0)
            {
                this.initialize(0);
            }

            // Check if we have a constant value
            else if (range.Length == 0)
            {
                // Yes, we will create a special histogram
                // bin to accommodate those constant values.
                this.initialize(1);

                this.ranges[0] = range.Min;
                this.ranges[1] = range.Max;
            }

            // Check if we have to auto-adjust the histogram
            //  bins according to some selection choice.
            else if (this.rule != BinAdjustmentRule.None)
            {
                // Yes, we will be recomputing the optimal number of bins
                int numberOfBins = NumberOfBins(values, range, this.rule);

                // Determine bin width based on the given number of bins
                double binWidth = range.Length / (double)numberOfBins;

                // Create bin structure
                this.initialize(numberOfBins);

                // Create ranges w/ fixed width
                this.initialize(range.Min, binWidth);
            }

            // Create histogram
            this.compute(values);
        }






        /// <summary>
        ///   Actually computes the histogram.
        /// </summary>
        /// 
        private void compute(double[] values)
        {
            int numberOfBins = this.values.Length;
            DoubleRange scale = new DoubleRange(0, numberOfBins);
            DoubleRange range = this.Range;

            if (range.Length == 0)
            {
                if (values.Length > 0)
                    this.values[0] = values.Length;
                return;
            }

            // Populate Bins
            for (int i = 0; i < values.Length; i++)
            {
                double v = values[i];

                // Convert the value to the range of histogram
                //  bins to check which bin the value belongs.
                int index = (int)Vector.Scale(v, range, scale);

                if (index < numberOfBins)
                    this.values[index]++;
                else if (this.inclusiveUpperBound)
                    this.values[numberOfBins - 1]++;
            }


            // If this is a cumulative histogram,
            //  accumulate values in the bins.
            if (this.cumulative)
            {
                for (int i = 1; i < this.values.Length; i++)
                    this.values[i] += this.values[i - 1];
            }
        }





        /// <summary>
        ///   Computes the optimum number of bins based on a <see cref="BinAdjustmentRule"/>.
        /// </summary>
        /// 
        public static int NumberOfBins(double[] values, DoubleRange range, BinAdjustmentRule rule)
        {
            switch (rule)
            {
                case BinAdjustmentRule.None:
                    return 0;

                case BinAdjustmentRule.Scott:
                    double h = (3.49 * Measures.StandardDeviation(values))
                        / System.Math.Pow(values.Length, 1.0 / 3.0);
                    return (int)Math.Ceiling(range.Length / h);

                case BinAdjustmentRule.Sturges:
                    return (int)Math.Ceiling(Math.Log(values.Length, 2));

                case BinAdjustmentRule.SquareRoot:
                    return (int)Math.Floor(Math.Sqrt(values.Length));

                default:
                    goto case BinAdjustmentRule.SquareRoot;
            }
        }


        /// <summary>
        ///   Integer array implicit conversion.
        /// </summary>
        /// 
        public static implicit operator int[](Histogram value)
        {
            return value.values;
        }

        /// <summary>
        ///   Converts this histogram into an integer array representation.
        /// </summary>
        /// 
        public int[] ToArray()
        {
            return this.values;
        }

        /// <summary>
        ///   Creates a histogram of values from a sample.
        /// </summary>
        /// 
        /// <param name="values">The values to be binned in the histogram.</param>
        /// 
        /// <returns>A histogram reflecting the distribution of values in the sample.</returns>
        /// 
        public Histogram FromData(double[] values)
        {
            var hist = new Histogram();
            hist.Compute(values);
            return hist;
        }

        /// <summary>
        ///   Subtracts one histogram from the other, storing
        ///   results in a new histogram, without changing the
        ///   current instance.
        /// </summary>
        /// 
        /// <param name="histogram">The histogram whose bin values will be subtracted.</param>
        /// 
        /// <returns>A new <see cref="Histogram"/> containing the result of this operation.</returns>
        /// 
        public Histogram Subtract(Histogram histogram)
        {
            return this.Subtract(histogram.Values);
        }

        /// <summary>
        ///   Subtracts one histogram from the other, storing
        ///   results in a new histogram, without changing the
        ///   current instance.
        /// </summary>
        /// 
        /// <param name="histogram">The histogram whose bin values will be subtracted.</param>
        /// 
        /// <returns>A new <see cref="Histogram"/> containing the result of this operation.</returns>
        /// 
        public Histogram Subtract(int[] histogram)
        {
            Histogram clone = this.Clone() as Histogram;
            for (int i = 0; i < this.values.Length; i++)
                clone.values[i] -= histogram[i];
            return clone;
        }

        /// <summary>
        ///   Adds one histogram from the other, storing
        ///   results in a new histogram, without changing the
        ///   current instance.
        /// </summary>
        /// 
        /// <param name="histogram">The histogram whose bin values will be added.</param>
        /// 
        /// <returns>A new <see cref="Histogram"/> containing the result of this operation.</returns>
        /// 
        public Histogram Add(Histogram histogram)
        {
            return this.Add(histogram.Values);
        }

        /// <summary>
        ///   Adds one histogram from the other, storing
        ///   results in a new histogram, without changing the
        ///   current instance.
        /// </summary>
        /// 
        /// <param name="histogram">The histogram whose bin values will be added.</param>
        /// 
        /// <returns>A new <see cref="Histogram"/> containing the result of this operation.</returns>
        /// 
        public Histogram Add(int[] histogram)
        {
            Histogram clone = this.Clone() as Histogram;
            for (int i = 0; i < this.values.Length; i++)
                clone.values[i] += histogram[i];
            return clone;
        }

        /// <summary>
        ///   Multiplies one histogram from the other, storing
        ///   results in a new histogram, without changing the
        ///   current instance.
        /// </summary>
        /// 
        /// <param name="histogram">The histogram whose bin values will be multiplied.</param>
        /// 
        /// <returns>A new <see cref="Histogram"/> containing the result of this operation.</returns>
        /// 
        public Histogram Multiply(int[] histogram)
        {
            Histogram clone = this.Clone() as Histogram;
            for (int i = 0; i < this.values.Length; i++)
                clone.values[i] *= histogram[i];
            return clone;
        }

        /// <summary>
        ///   Multiplies one histogram from the other, storing
        ///   results in a new histogram, without changing the
        ///   current instance.
        /// </summary>
        /// 
        /// <param name="value">The value to be multiplied.</param>
        /// 
        /// <returns>A new <see cref="Histogram"/> containing the result of this operation.</returns>
        /// 
        public Histogram Multiply(int value)
        {
            Histogram clone = this.Clone() as Histogram;
            for (int i = 0; i < this.values.Length; i++)
                clone.values[i] *= value;
            return clone;
        }

        /// <summary>
        ///   Adds a value to each histogram bin.
        /// </summary>
        /// 
        /// <param name="value">The value to be added.</param>
        /// 
        /// <returns>A new <see cref="Histogram"/> containing the result of this operation.</returns>
        /// 
        public Histogram Add(int value)
        {
            Histogram clone = this.Clone() as Histogram;
            for (int i = 0; i < this.values.Length; i++)
                clone.values[i] += value;
            return clone;
        }

        /// <summary>
        ///   Subtracts a value to each histogram bin.
        /// </summary>
        /// 
        /// <param name="value">The value to be subtracted.</param>
        /// 
        /// <returns>A new <see cref="Histogram"/> containing the result of this operation.</returns>
        /// 
        public Histogram Subtract(int value)
        {
            Histogram clone = this.Clone() as Histogram;
            for (int i = 0; i < this.values.Length; i++)
                clone.values[i] -= value;
            return clone;
        }

        /// <summary>
        ///   Creates a new object that is a copy of the current instance.
        /// </summary>
        /// 
        /// <returns>
        ///   A new object that is a copy of this instance.
        /// </returns>
        /// 
        public object Clone()
        {
            Histogram clone = new Histogram();
            clone.initialize(this.values.Length);

            for (int i = 0; i < clone.ranges.Length; i++)
                clone.ranges[i] = this.ranges[i];

            for (int i = 0; i < clone.values.Length; i++)
                clone.values[i] = this.values[i];

            clone.cumulative = this.cumulative;
            clone.inclusiveUpperBound = this.inclusiveUpperBound;
            clone.rule = this.rule;

            return clone;
        }
    }
}