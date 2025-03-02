﻿// Accord Imaging Library
// The Accord.NET Framework
// http://accord-framework.net
//
// Copyright © César Souza, 2009-2017
// cesarsouza at gmail.com
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

namespace Openize.Accord.Imaging.Filters
{
    using System;
    using System.Collections.Generic;
    using System.Drawing.Imaging;
    using AForge.Imaging;
    using AForge.Imaging.Filters.Base_classes;

    /// <summary>
    ///   Exponential filter.
    /// </summary>
    /// 
    /// <remarks>
    ///   Simple exp image filter. Applies the <see cref="System.Math.Exp"/>
    ///   function for each pixel in the image, clipping values as needed.
    ///   The resultant image can be converted back using the <see cref="Logarithm"/>
    ///   filter.
    /// </remarks>
    /// 
    /// <example>
    /// <code>
    ///   Bitmap input = ... 
    /// 
    ///   // Apply log
    ///   Logarithm log = new Logarithm();
    ///   Bitmap output = log.Apply(input);
    /// 
    ///   // Revert log
    ///   Exponential exp = new Exponential();
    ///   Bitmap reconstruction = exp.Apply(output);
    /// 
    ///   // Show results on screen
    ///   ImageBox.Show("input", input);
    ///   ImageBox.Show("output", output);
    ///   ImageBox.Show("reconstruction", reconstruction);
    /// </code>
    /// </example>
    /// 
    /// <seealso cref="Logarithm"/>
    /// 
    public class Exponential : BaseInPlaceFilter
    {
        private Dictionary<PixelFormat, PixelFormat> formatTranslations;

        /// <summary>
        ///   Format translations dictionary.
        /// </summary>
        /// 
        public override Dictionary<PixelFormat, PixelFormat> FormatTranslations
        {
            get { return this.formatTranslations; }
        }

        /// <summary>
        ///   Initializes a new instance of the <see cref="Exponential"/> class.
        /// </summary>
        /// 
        public Exponential()
        {
            this.formatTranslations = new Dictionary<PixelFormat, PixelFormat>();
            this.formatTranslations[PixelFormat.Format8bppIndexed] = PixelFormat.Format8bppIndexed;
            this.formatTranslations[PixelFormat.Format16bppGrayScale] = PixelFormat.Format16bppGrayScale;
            this.formatTranslations[PixelFormat.Format24bppRgb] = PixelFormat.Format24bppRgb;
            this.formatTranslations[PixelFormat.Format32bppRgb] = PixelFormat.Format32bppRgb;
            this.formatTranslations[PixelFormat.Format32bppArgb] = PixelFormat.Format32bppArgb;
            this.formatTranslations[PixelFormat.Format48bppRgb] = PixelFormat.Format48bppRgb;
            this.formatTranslations[PixelFormat.Format64bppArgb] = PixelFormat.Format64bppArgb;
        }

        /// <summary>
        ///   Process the filter on the specified image.
        /// </summary>
        /// 
        /// <param name="image">Source image data.</param>
        /// 
        protected unsafe override void ProcessFilter(UnmanagedImage image)
        {
            int width = image.Width;
            int height = image.Height;
            PixelFormat format = image.PixelFormat;
            int pixelSize = System.Drawing.Bitmap.GetPixelFormatSize(format) / 8;

            int lineWidth = width * pixelSize;

            int srcStride = image.Stride;
            int srcOffset = srcStride - lineWidth;
            double scale = 255.0 / Math.Log(255);

            byte* src = (byte*)image.ImageData.ToPointer();


            for (int y = 0; y < height; y++)
            {
                for (int x = 0; x < lineWidth; x++, src++)
                {
                    if (*src > 0)
                    {
                        double v = Math.Exp(*src / scale);

                        *src = (byte)(v > 0 ? (v < 255 ? v : 255) : 0);
                    }

                }

                src += srcOffset;
            }
        }
    }
}
