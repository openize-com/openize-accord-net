// Accord Imaging Library
// The Accord.NET Framework
// http://accord-framework.net
//
// Copyright © Diego Catalano, 2013
// diego.catalano at live.com
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

namespace FileFormat.Accord.Imaging.Filters
{
    using System;
    using AForge.Imaging.Filters.Convolution;

    /// <summary>
    ///   High boost filter.
    /// </summary>
    /// 
    /// <remarks>
    /// <para>
    ///   The High-boost filter can be used to emphasize high frequency
    ///   components (i.e. points of contrast) without removing the low
    ///   frequency ones.</para>
    ///   
    /// <para>
    ///  This filter implementation has been contributed by Diego Catalano.</para>
    /// </remarks>
    /// 
    public class HighBoost : Convolution
    {
        private int boost = 8;
        private int size = 3;


        /// <summary>
        /// Kernel size, [3, 21].
        /// </summary>
        /// 
        /// <remarks><para>Size of Gaussian kernel.</para>
        /// 
        /// <para>Default value is set to <b>5</b>.</para>
        /// </remarks>
        /// 
        public int Size
        {
            get { return this.size; }
            set
            {
                this.size = Math.Max(3, Math.Min(21, value | 1));
                this.createFilter();
            }
        }

        /// <summary>
        ///   Gets or sets the boost value. Default is 9.
        /// </summary>
        /// 
        public int Boost
        {
            get { return this.boost; }
            set
            {
                this.boost = value;
                this.Kernel[this.size / 2, this.size / 2] = this.boost;
            }
        }

        /// <summary>
        ///   Initializes a new instance of the <see cref="HighBoost"/> class.
        /// </summary>
        /// 
        public HighBoost()
        {
            this.createFilter();
            base.ProcessAlpha = true;
        }

        /// <summary>
        ///   Initializes a new instance of the <see cref="HighBoost"/> class.
        /// </summary>
        /// 
        /// <param name="boost">The boost value. Default is 8.</param>
        /// 
        public HighBoost(int boost)
        {
            this.Boost = boost;
            base.ProcessAlpha = true;
        }

        /// <summary>
        ///   Initializes a new instance of the <see cref="HighBoost"/> class.
        /// </summary>
        /// 
        /// <param name="boost">The boost value. Default is 8.</param>
        /// <param name="size">The kernel size. Default is 3.</param>
        /// 
        public HighBoost(int boost, int size)
        {
            this.Size = size;
            this.Boost = boost;
            base.ProcessAlpha = true;
        }


        private void createFilter()
        {
            // integer kernel
            int[,] kernel = new int[this.size, this.size];

            for (int i = 0; i < this.size; i++)
                for (int j = 0; j < this.size; j++)
                    kernel[i, j] = -1;

            kernel[this.size / 2, this.size / 2] = this.boost;

            this.Kernel = kernel;
        }
    }
}
