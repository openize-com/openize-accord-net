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

#if !NETSTANDARD1_4
namespace FileFormat.Accord.Statistics.Filters
{
    using System;
    using System.Data;
    using Base;
    using global::Accord.Math;
    using Math.Matrix;

    /// <summary>
    ///   Time-series windowing filter.
    /// </summary>
    /// 
    /// <remarks>
    ///   This filter splits a time-series into overlapping time
    ///   windows, with optional associated output values. This
    ///   filter can be used to create time-window databases for
    ///   time-series regression and latent-state identification.
    /// </remarks>
    /// 
    [Serializable]
    public class Windowing : BaseFilter<Windowing.Options, Windowing>
    {

        /// <summary>
        ///   Gets or sets the length of the time-windows
        ///   that should be extracted from the sequences.
        /// </summary>
        /// 
        public int WindowSize { get; set; }

        /// <summary>
        ///   Gets or sets the step size that should be used
        ///   when extracting windows. If set to the same number
        ///   as the <see cref="WindowSize"/>, windows will not
        ///   overlap. Default is 1.
        /// </summary>
        /// 
        public int StepSize { get; set; }

        /// <summary>
        ///   Creates a new time segmentation filter.
        /// </summary>
        /// 
        public Windowing()
        {
            this.init(1, 1);
        }

        /// <summary>
        ///   Creates a new time segmentation filter.
        /// </summary>
        /// 
        /// <param name="windowSize">The size of the time windows to be extracted.</param>
        /// 
        public Windowing(int windowSize)
        {
            this.init(windowSize, 1);
        }

        /// <summary>
        ///   Creates a new time segmentation filter.
        /// </summary>
        /// 
        /// <param name="windowSize">The size of the time windows to be extracted.</param>
        /// <param name="steps">The number of elements between two taken windows. If set to
        ///   the same number of <paramref name="windowSize"/>, the windows will not overlap. 
        ///   Default is 1.</param>
        /// 
        public Windowing(int windowSize, int steps)
        {
            this.init(windowSize, steps);
        }

        private void init(int inputs, int stepSize)
        {
            this.WindowSize = inputs;
            this.StepSize = stepSize;
        }



        /// <summary>
        ///   Processes the current filter.
        /// </summary>
        /// 
        protected override DataTable ProcessFilter(DataTable data)
        {
            if (this.Columns.Count > 1)
                throw new InvalidOperationException();

            Options column = this.Columns[0];

            string name = column.ColumnName;

            double[] outputs;
            double[][] inputs;

            // Get entire time series
            double[] sequence = data.Columns[name].ToArray();

            inputs = this.Apply(sequence, out outputs);

            return Matrix.ToTable(inputs.InsertColumn(outputs));
        }

        /// <summary>
        ///   Applies the filter to a time series.
        /// </summary>
        /// 
        /// <param name="samples">The source time series.</param>
        /// 
        /// <returns>The time-windows extracted from the time-series.</returns>
        /// 
        public double[][] Apply(double[] samples)
        {
            int m = (samples.Length - this.WindowSize);
            int n = (int)System.Math.Floor(m / (double)this.StepSize);

            double[][] inputs = new double[n][];
            for (int i = 0; i < inputs.Length; i++)
                inputs[i] = this.rectangular(samples, i * this.StepSize);

            return inputs;
        }

        /// <summary>
        ///   Applies the filter to a time series.
        /// </summary>
        /// 
        /// <param name="samples">The source time series.</param>
        /// <param name="outputs">The output associated with each time-window.</param>
        /// 
        /// <returns>The time-windows extracted from the time-series.</returns>
        /// 
        public double[][] Apply(double[] samples, out double[] outputs)
        {
            int OutputCount = 1;
            int m = (samples.Length - this.WindowSize - OutputCount);
            int n = (int)System.Math.Floor(m / (double)this.StepSize);

            double[][] inputs = new double[n][];
            outputs = new double[n];

            for (int i = 0; i < inputs.Length; i++)
            {
                int j = i * this.StepSize;
                inputs[i] = this.rectangular(samples, j);
                outputs[i] = samples[j + this.WindowSize];
            }

            return inputs;
        }


        private double[] rectangular(double[] signal, int sampleIndex)
        {
            double[] window = new double[this.WindowSize];
            Array.Copy(signal, sampleIndex, window, 0, this.WindowSize);
            return window;
        }






        /// <summary>
        ///   Options for segmenting a time-series contained inside a column.
        /// </summary>
        ///
        [Serializable]
        public class Options : ColumnOptionsBase<Windowing>
        {

            /// <summary>
            ///   Constructs a new Options object.
            /// </summary>
            /// 
            public Options()
                : base("New column")
            {
            }

        }

    }
}
#endif