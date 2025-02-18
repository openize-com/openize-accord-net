﻿// Accord Imaging Library
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

namespace Openize.Accord.Imaging.Filters
{
    using System;
    using System.Collections.Generic;
    using System.Drawing.Imaging;
    using AForge.Imaging;
    using AForge.Imaging.Filters.Base_classes;

    /// <summary>
    ///   Combine channel filter.
    /// </summary>
    /// 
    public class CombineChannel : BaseInPlaceFilter
    {

        private Dictionary<PixelFormat, PixelFormat> formatTranslations = new Dictionary<PixelFormat, PixelFormat>();

        private int baseWidth;
        private int baseHeight;
        private UnmanagedImage[] channels;

        /// <summary>
        ///   Format translations dictionary.
        /// </summary>
        /// 
        /// <remarks>
        ///   <para>The dictionary defines, which pixel formats are supported for
        ///   source images and which pixel format will be used for resulting image.</para>
        /// 
        ///   <para>See <see cref="P:AForge.Imaging.Filters.IFilterInformation.FormatTranslations"/>
        ///   for more information.</para>
        /// </remarks>
        /// 
        public override Dictionary<PixelFormat, PixelFormat> FormatTranslations
        {
            get { return this.formatTranslations; }
        }

        /// <summary>
        ///   Constructs a new CombineChannel filter.
        /// </summary>
        /// 
        public CombineChannel(params UnmanagedImage[] channels)
        {
            if (channels == null)
                throw new ArgumentNullException("channels");
            if (channels.Length < 2)
                throw new ArgumentException("There must be at least two channels to be combined.", "channels");

            this.channels = channels;
            this.baseWidth = channels[0].Width;
            this.baseHeight = channels[0].Height;

            foreach (var c in channels)
                if (c.Width != this.baseWidth || c.Height != this.baseHeight)
                    throw new ArgumentException("All images must have the same dimensions.", "channels");

            this.formatTranslations[PixelFormat.Format64bppArgb] = PixelFormat.Format64bppArgb;
            this.formatTranslations[PixelFormat.Format48bppRgb] = PixelFormat.Format48bppRgb;
            this.formatTranslations[PixelFormat.Format32bppRgb] = PixelFormat.Format32bppRgb;
            this.formatTranslations[PixelFormat.Format32bppArgb] = PixelFormat.Format32bppArgb;
            this.formatTranslations[PixelFormat.Format24bppRgb] = PixelFormat.Format24bppRgb;
        }


        /// <summary>
        ///   Process the filter on the specified image.
        /// </summary>
        /// 
        /// <param name="image">Source image data.</param>
        /// 
        protected unsafe override void ProcessFilter(UnmanagedImage image)
        {
            int pixelSize = System.Drawing.Image.GetPixelFormatSize(image.PixelFormat) / 8;

            // get source image size
            int width = image.Width;
            int height = image.Height;

            // check is the same size
            if (image.Height != this.baseHeight || image.Width != this.baseWidth)
                throw new InvalidImagePropertiesException("Image does not have expected dimensions.", "image");

            if (pixelSize == 8)
            {
                // for each channel
                for (int c = 0; c < this.channels.Length; c++)
                {
                    byte* dst = (byte*)((int)image.ImageData + c);
                    byte* src = (byte*)this.channels[c].ImageData;

                    // copy channel to image
                    for (int y = 0; y < height; y++)
                    {
                        for (int x = 0; x < width; x++)
                        {
                            *(dst += pixelSize) = *(src++);
                        }
                    }
                }
            }
            else if (pixelSize == 16)
            {
                // for each channel
                for (int c = 0; c < this.channels.Length; c++)
                {
                    short* dst = (short*)((int)image.ImageData + c);
                    short* src = (short*)this.channels[c].ImageData;

                    // copy channel to image
                    for (int y = 0; y < height; y++)
                    {
                        for (int x = 0; x < width; x++)
                        {
                            *(dst += pixelSize) = *(src++);
                        }
                    }
                }
            }
            else
            {
                throw new UnsupportedImageFormatException("Unsupported pixel size.");
            }
        }

    }
}
