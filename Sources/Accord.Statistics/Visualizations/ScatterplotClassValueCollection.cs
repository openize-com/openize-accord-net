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

namespace Openize.Accord.Statistics.Visualizations
{
    using System;
    using System.Collections.Generic;
    using Openize.Accord.Math.Matrix;
    using Openize.Accord.Core.AForge.Core;

    /// <summary>
    ///   Scatter Plot class.
    /// </summary>
    /// 
    [Serializable]
    public class ScatterplotClassValueCollection : IEnumerable<DoublePoint>
    {
        private Scatterplot parent;

        private int index;

        /// <summary>
        ///   Gets the integer label associated with this class.
        /// </summary>
        /// 
        public int Label { get { return this.parent.LabelValues[this.index]; } }

        /// <summary>
        ///   Gets the indices of all points of this class.
        /// </summary>
        /// 
        public int[] Indices
        {
            get { return this.parent.LabelAxis.Find(x => x == this.Label); }
        }

        /// <summary>
        ///   Gets all X values of this class.
        /// </summary>
        /// 
        public double[] XAxis
        {
            get { return this.parent.XAxis.Get(this.Indices); }
        }

        /// <summary>
        ///   Gets all Y values of this class.
        /// </summary>
        /// 
        public double[] YAxis
        {
            get { return this.parent.YAxis.Get(this.Indices); }
        }

        /// <summary>
        ///   Gets or sets the class' text label.
        /// </summary>
        /// 
        public string Text
        {
            get { return this.parent.LabelNames[this.index]; }
            set { this.parent.LabelNames[this.index] = value; }
        }

        internal ScatterplotClassValueCollection(Scatterplot parent, int index)
        {
            this.parent = parent;
            this.index = index;
        }

        /// <summary>
        ///   Returns an enumerator that iterates through a collection.
        /// </summary>
        /// <returns>
        ///   An <see cref="T:System.Collections.IEnumerator"/> object that can be used to iterate through the collection.
        /// </returns>
        /// 
        public IEnumerator<DoublePoint> GetEnumerator()
        {
            for (int i = 0; i < this.parent.LabelAxis.Length; i++)
            {
                if (this.parent.LabelAxis[i] == this.Label)
                    yield return new DoublePoint(this.parent.XAxis[i], this.parent.YAxis[i]);
            }

            yield break;
        }

        /// <summary>
        ///   Returns an enumerator that iterates through a collection.
        /// </summary>
        /// <returns>
        ///   An <see cref="T:System.Collections.IEnumerator"/> object that can be used to iterate through the collection.
        /// </returns>
        /// 
        System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator()
        {
            return this.GetEnumerator();
        }

    }
}
