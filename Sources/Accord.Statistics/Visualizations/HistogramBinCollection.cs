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

    /// <summary>
    ///   Collection of Histogram bins. This class cannot be instantiated.
    /// </summary>
    /// 
    [Serializable]
    public class HistogramBinCollection : System.Collections.ObjectModel.ReadOnlyCollection<HistogramBin>
    {
        internal HistogramBinCollection(HistogramBin[] objects)
            : base(objects)
        {

        }

        /// <summary>
        ///   Searches for a bin containing the specified value.
        /// </summary>
        /// 
        /// <param name="value">The value to search for.</param>
        /// 
        /// <returns>The histogram bin containing the searched value.</returns>
        /// 
        public HistogramBin Search(double value)
        {
            foreach (HistogramBin bin in this)
            {
                if (bin.Contains(value))
                    return bin;
            }

            return null;
        }

        /// <summary>
        ///   Attempts to find the index of the bin containing the specified value.
        /// </summary>
        /// 
        /// <param name="value">The value to search for.</param>
        /// 
        /// <returns>The index of the bin containing the specified value.</returns>
        /// 
        public int SearchIndex(double value)
        {
            for (int i = 0; i < this.Count; i++)
            {
                if (this[i].Contains(value))
                    return i;
            }
            return -1;
        }

    }
  
}
