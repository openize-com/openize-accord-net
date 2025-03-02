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
    using Openize.Accord.Core.Ranges;

    /// <summary>
    ///   Histogram Bin
    /// </summary>
    /// 
    /// <remarks>
    /// <para>
    ///   A "bin" is a container, where each element stores the total number of observations of a sample
    ///   whose values lie within a given range. A histogram of a sample consists of a list of such bins
    ///   whose range does not overlap with each other; or in other words, bins that are mutually exclusive.</para>
    /// <para>
    ///   Unless <see cref="Histogram.InclusiveUpperBound"/> is true, the ranges of all bins <c>i</c> are
    ///   defined as Edge[i] &lt;= x &lt; Edge[i+1]. Otherwise, the last bin will have an inclusive upper
    ///   bound (i.e. will be defined as Edge[i] &lt;= x &lt;= Edge[i+1].</para>  
    /// </remarks>
    /// 
    [Serializable]
    public class HistogramBin
    {

        private int index;
        private Histogram histogram;


        internal HistogramBin(Histogram histogram, int index)
        {
            this.index = index;
            this.histogram = histogram;
        }

        /// <summary>
        ///   Gets the actual range of data this bin represents.
        /// </summary>
        /// 
        public DoubleRange Range
        {
            get
            {
                double min = this.histogram.Edges[this.index];
                double max = this.histogram.Edges[this.index + 1];
                return new DoubleRange(min, max);
            }
        }

        /// <summary>
        ///   Gets the Width (range) for this histogram bin.
        /// </summary>
        /// 
        public double Width
        {
            get { return this.histogram.Edges[this.index + 1] - this.histogram.Edges[this.index]; }
        }

        /// <summary>
        ///   Gets the Value (number of occurrences of a variable in a range)
        ///   for this histogram bin.
        /// </summary>
        /// 
        public int Value
        {
            get { return this.histogram.Values[this.index]; }
            set { this.histogram.Values[this.index] = value; }
        }

        /// <summary>
        ///   Gets whether the Histogram Bin contains the given value.
        /// </summary>
        /// 
        public bool Contains(double value)
        {
            if (this.histogram.InclusiveUpperBound &&
                this.index == this.histogram.Bins.Count - 1)
            {
                return value >= this.histogram.Edges[this.index] &&
                       value <= this.histogram.Edges[this.index + 1];
            }
            else
            {
                return value >= this.histogram.Edges[this.index] &&
                       value < this.histogram.Edges[this.index + 1];
            }
        }
    }
}
