// AForge Image Processing Library
// AForge.NET framework
// http://www.aforgenet.com/framework/
//
// Copyright © Andrew Kirillov, 2005-2009
// andrew.kirillov@aforgenet.com
//
// Accord Imaging Library
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

namespace Openize.Accord.Imaging.Blob_Processing
{
    using System.Drawing;
    using System.Drawing.Imaging;
    using AForge.Imaging;
    using Colors;
    using Math.AForge.Math.Geometry;

    /// <summary>
    /// Blob counter - counts objects in image, which are separated by black background.
    /// </summary>
    /// 
    /// <remarks><para>The class counts and extracts stand alone objects in
    /// images using connected components labeling algorithm.</para>
    /// 
    /// <para><note>The algorithm treats all pixels with values less or equal to <see cref="BackgroundThreshold"/>
    /// as background, but pixels with higher values are treated as objects' pixels.</note></para>
    /// 
    /// <para>For blobs' searching the class supports 8 bpp indexed grayscale images and
    /// 24/32 bpp color images that are at least two pixels wide. Images that are one
    /// pixel wide can be processed if they are rotated first, or they can be processed
    /// with <see cref="RecursiveBlobCounter"/>.
    /// See documentation about <see cref="BlobCounterBase"/> for information about which
    /// pixel formats are supported for extraction of blobs.</para>
    /// </remarks>
    /// 
    /// <example>
    /// <para>
    ///   For a more complete example, please take a look at the 
    ///   <a href="https://github.com/accord-net/framework/wiki/Sample-applications#blobs-detection">
    ///   Blob detection sample application</a>.</para>
    /// 
    /// <para>
    ///   A simplified code example is shown below:</para>
    ///   
    /// <code source="Unit Tests\Accord.Tests.Controls\BlobCounterTest.cs" region="doc_process" />
    /// <para>
    ///   The input image is shown below:</para>
    ///   <img src="..\images\imaging\blobs\blobs-input.png" />
    /// <para>
    ///   And the output image should be:</para>
    ///   <img src="..\images\imaging\blobs\blobs-output.png" />
    /// </example>
    /// 
    /// <seealso cref="GrahamConvexHull"/>
    /// <seealso cref="PointsCloud"/>
    /// 
    public class BlobCounter : BlobCounterBase
    {
        private byte backgroundThresholdR = 0;
        private byte backgroundThresholdG = 0;
        private byte backgroundThresholdB = 0;

        /// <summary>
        /// Background threshold's value.
        /// </summary>
        /// 
        /// <remarks><para>The property sets threshold value for distinguishing between background
        /// pixel and objects' pixels. All pixel with values less or equal to this property are
        /// treated as background, but pixels with higher values are treated as objects' pixels.</para>
        /// 
        /// <para><note>In the case of colour images a pixel is treated as objects' pixel if <b>any</b> of its
        /// RGB values are higher than corresponding values of this threshold.</note></para>
        /// 
        /// <para><note>For processing grayscale image, set the property with all RGB components eqaul.</note></para>
        ///
        /// <para>Default value is set to <b>(0, 0, 0)</b> - black colour.</para></remarks>
        /// 
        public Color BackgroundThreshold
        {
            get { return Color.FromArgb(this.backgroundThresholdR, this.backgroundThresholdG, this.backgroundThresholdB); }
            set
            {
                this.backgroundThresholdR = value.R;
                this.backgroundThresholdG = value.G;
                this.backgroundThresholdB = value.B;
            }
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="BlobCounter"/> class.
        /// </summary>
        /// 
        /// <remarks>Creates new instance of the <see cref="BlobCounter"/> class with
        /// an empty objects map. Before using methods, which provide information about blobs
        /// or extract them, the <see cref="BlobCounterBase.ProcessImage(Bitmap)"/>,
        /// <see cref="BlobCounterBase.ProcessImage(BitmapData)"/> or <see cref="BlobCounterBase.ProcessImage(UnmanagedImage)"/>
        /// method should be called to collect objects map.</remarks>
        /// 
        public BlobCounter()
        {
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="BlobCounter"/> class.
        /// </summary>
        /// 
        /// <param name="image">Image to look for objects in.</param>
        /// 
        public BlobCounter(Bitmap image)
            : base(image)
        {
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="BlobCounter"/> class.
        /// </summary>
        /// 
        /// <param name="imageData">Image data to look for objects in.</param>
        /// 
        public BlobCounter(BitmapData imageData)
            : base(imageData)
        {
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="BlobCounter"/> class.
        /// </summary>
        /// 
        /// <param name="image">Unmanaged image to look for objects in.</param>
        /// 
        public BlobCounter(UnmanagedImage image)
            : base(image)
        {
        }

        /// <summary>
        /// Actual objects map building.
        /// </summary>
        /// 
        /// <param name="image">Unmanaged image to process.</param>
        /// 
        /// <remarks>The method supports 8 bpp indexed grayscale images and 24/32 bpp color images.</remarks>
        /// 
        /// <exception cref="UnsupportedImageFormatException">Unsupported pixel format of the source image.</exception>
        /// <exception cref="InvalidImagePropertiesException">Cannot process images that are one pixel wide. Rotate the image
        /// or use <see cref="RecursiveBlobCounter"/>.</exception>
        /// 
        protected override void BuildObjectsMap(UnmanagedImage image)
        {
            int stride = image.Stride;

            // check pixel format
            if ((image.PixelFormat != PixelFormat.Format8bppIndexed) &&
                 (image.PixelFormat != PixelFormat.Format24bppRgb) &&
                 (image.PixelFormat != PixelFormat.Format32bppRgb) &&
                 (image.PixelFormat != PixelFormat.Format32bppArgb) &&
                 (image.PixelFormat != PixelFormat.Format32bppPArgb))
            {
                throw new UnsupportedImageFormatException("Unsupported pixel format of the source image.");
            }

            // we don't want one pixel width images
            if (this.ImageWidth == 1)
                throw new InvalidImagePropertiesException("BlobCounter cannot process images that are one pixel wide. Rotate the image or use RecursiveBlobCounter.");

            int imageWidthM1 = this.ImageWidth - 1;

            // allocate labels array
            this.ObjectLabels = new int[this.ImageWidth * this.ImageHeight];

            // initial labels count
            int labelsCount = 0;

            // create map
            int maxObjects = ((this.ImageWidth / 2) + 1) * ((this.ImageHeight / 2) + 1) + 1;
            int[] map = new int[maxObjects];

            // initially map all labels to themself
            for (int i = 0; i < maxObjects; i++)
                map[i] = i;

            // do the job
            unsafe
            {
                byte* src = (byte*)image.ImageData.ToPointer();
                int p = 0;

                if (image.PixelFormat == PixelFormat.Format8bppIndexed)
                {
                    int offset = stride - this.ImageWidth;

                    image.CheckBounds(src);

                    // 1 - for pixels of the first row
                    if (*src > this.backgroundThresholdG)
                        this.ObjectLabels[p] = ++labelsCount;
                    ++src;
                    ++p;

                    // process the rest of the first row
                    for (int x = 1; x < this.ImageWidth; x++, src++, p++)
                    {
                    image.CheckBounds(src);

                        // check if we need to label current pixel
                        if (*src > this.backgroundThresholdG)
                        {
                            image.CheckBounds(src - 1);

                            // check if the previous pixel already was labeled
                            if (src[-1] > this.backgroundThresholdG)
                            {
                                // label current pixel, as the previous
                                this.ObjectLabels[p] = this.ObjectLabels[p - 1];
                            }
                            else
                            {
                                // create new label
                                this.ObjectLabels[p] = ++labelsCount;
                            }
                        }
                    }
                    src += offset;

                    // 2 - for other rows
                    // for each row
                    for (int y = 1; y < this.ImageHeight; y++)
                    {
                        // for the first pixel of the row, we need to check
                        // only upper and upper-right pixels

                        image.CheckBounds(src );

                        if (*src > this.backgroundThresholdG)
                        {
                            image.CheckBounds(src - stride);
                            image.CheckBounds(src + 1 - stride);

                            // check surrounding pixels
                            if (src[-stride] > this.backgroundThresholdG)
                            {
                                // label current pixel, as the above
                                this.ObjectLabels[p] = this.ObjectLabels[p - this.ImageWidth];
                            }
                            else if (src[1 - stride] > this.backgroundThresholdG)
                            {
                                // label current pixel, as the above right
                                this.ObjectLabels[p] = this.ObjectLabels[p + 1 - this.ImageWidth];
                            }
                            else
                            {
                                // create new label
                                this.ObjectLabels[p] = ++labelsCount;
                            }
                        }
                        ++src;
                        ++p;

                        // check left pixel and three upper pixels for the rest of pixels
                        for (int x = 1; x < imageWidthM1; x++, src++, p++)
                        {
                            image.CheckBounds(src );

                            if (*src > this.backgroundThresholdG)
                            {
                                image.CheckBounds(src - 1);
                                image.CheckBounds(src - 1 - stride);
                                image.CheckBounds(src - stride);
                                image.CheckBounds(src +1 - stride);

                                // check surrounding pixels
                                if (src[-1] > this.backgroundThresholdG)
                                {
                                    // label current pixel, as the left
                                    this.ObjectLabels[p] = this.ObjectLabels[p - 1];
                                }
                                else if (src[-1 - stride] > this.backgroundThresholdG)
                                {
                                    // label current pixel, as the above left
                                    this.ObjectLabels[p] = this.ObjectLabels[p - 1 - this.ImageWidth];
                                }
                                else if (src[-stride] > this.backgroundThresholdG)
                                {
                                    // label current pixel, as the above
                                    this.ObjectLabels[p] = this.ObjectLabels[p - this.ImageWidth];
                                }

                                if (src[1 - stride] > this.backgroundThresholdG)
                                {
                                    if (this.ObjectLabels[p] == 0)
                                    {
                                        // label current pixel, as the above right
                                        this.ObjectLabels[p] = this.ObjectLabels[p + 1 - this.ImageWidth];
                                    }
                                    else
                                    {
                                        int l1 = this.ObjectLabels[p];
                                        int l2 = this.ObjectLabels[p + 1 - this.ImageWidth];

                                        if ((l1 != l2) && (map[l1] != map[l2]))
                                        {
                                            // merge
                                            if (map[l1] == l1)
                                            {
                                                // map left value to the right
                                                map[l1] = map[l2];
                                            }
                                            else if (map[l2] == l2)
                                            {
                                                // map right value to the left
                                                map[l2] = map[l1];
                                            }
                                            else
                                            {
                                                // both values already mapped
                                                map[map[l1]] = map[l2];
                                                map[l1] = map[l2];
                                            }

                                            // reindex
                                            for (int i = 1; i <= labelsCount; i++)
                                            {
                                                if (map[i] != i)
                                                {
                                                    // reindex
                                                    int j = map[i];
                                                    while (j != map[j])
                                                    {
                                                        j = map[j];
                                                    }
                                                    map[i] = j;
                                                }
                                            }
                                        }
                                    }
                                }

                                // label the object if it is not yet
                                if (this.ObjectLabels[p] == 0)
                                {
                                    // create new label
                                    this.ObjectLabels[p] = ++labelsCount;
                                }
                            }
                        }

                        // for the last pixel of the row, we need to check
                        // only upper and upper-left pixels
                        image.CheckBounds(src);

                        if (*src > this.backgroundThresholdG)
                        {
                            image.CheckBounds(src - 1);
                            image.CheckBounds(src - 1 - stride);
                            image.CheckBounds(src - stride);

                            // check surrounding pixels
                            if (src[-1] > this.backgroundThresholdG)
                            {
                                // label current pixel, as the left
                                this.ObjectLabels[p] = this.ObjectLabels[p - 1];
                            }
                            else if (src[-1 - stride] > this.backgroundThresholdG)
                            {
                                // label current pixel, as the above left
                                this.ObjectLabels[p] = this.ObjectLabels[p - 1 - this.ImageWidth];
                            }
                            else if (src[-stride] > this.backgroundThresholdG)
                            {
                                // label current pixel, as the above
                                this.ObjectLabels[p] = this.ObjectLabels[p - this.ImageWidth];
                            }
                            else
                            {
                                // create new label
                                this.ObjectLabels[p] = ++labelsCount;
                            }
                        }
                        ++src;
                        ++p;

                        src += offset;
                    }
                }
                else
                {
                    // color images
                    int pixelSize = Bitmap.GetPixelFormatSize(image.PixelFormat) / 8;
                    int offset = stride - this.ImageWidth * pixelSize;

                    int strideM1 = stride - pixelSize;
                    int strideP1 = stride + pixelSize;

                    // 1 - for pixels of the first row
                    if ((src[RGB.R] | src[RGB.G] | src[RGB.B]) != 0)
                        this.ObjectLabels[p] = ++labelsCount;
                    src += pixelSize;
                    ++p;

                    // process the rest of the first row
                    for (int x = 1; x < this.ImageWidth; x++, src += pixelSize, p++)
                    {
                        // check if we need to label current pixel
                        if ((src[RGB.R] > this.backgroundThresholdR) ||
                             (src[RGB.G] > this.backgroundThresholdG) ||
                             (src[RGB.B] > this.backgroundThresholdB))
                        {
                            // check if the previous pixel already was labeled
                            if ((src[RGB.R - pixelSize] > this.backgroundThresholdR) ||
                                 (src[RGB.G - pixelSize] > this.backgroundThresholdG) ||
                                 (src[RGB.B - pixelSize] > this.backgroundThresholdB))
                            {
                                // label current pixel, as the previous
                                this.ObjectLabels[p] = this.ObjectLabels[p - 1];
                            }
                            else
                            {
                                // create new label
                                this.ObjectLabels[p] = ++labelsCount;
                            }
                        }
                    }
                    src += offset;

                    // 2 - for other rows
                    // for each row
                    for (int y = 1; y < this.ImageHeight; y++)
                    {
                        // for the first pixel of the row, we need to check
                        // only upper and upper-right pixels
                        if ((src[RGB.R] > this.backgroundThresholdR) ||
                             (src[RGB.G] > this.backgroundThresholdG) ||
                             (src[RGB.B] > this.backgroundThresholdB))
                        {
                            // check surrounding pixels
                            if ((src[RGB.R - stride] > this.backgroundThresholdR) ||
                                 (src[RGB.G - stride] > this.backgroundThresholdG) ||
                                 (src[RGB.B - stride] > this.backgroundThresholdB))
                            {
                                // label current pixel, as the above
                                this.ObjectLabels[p] = this.ObjectLabels[p - this.ImageWidth];
                            }
                            else if ((src[RGB.R - strideM1] > this.backgroundThresholdR) ||
                                      (src[RGB.G - strideM1] > this.backgroundThresholdG) ||
                                      (src[RGB.B - strideM1] > this.backgroundThresholdB))
                            {
                                // label current pixel, as the above right
                                this.ObjectLabels[p] = this.ObjectLabels[p + 1 - this.ImageWidth];
                            }
                            else
                            {
                                // create new label
                                this.ObjectLabels[p] = ++labelsCount;
                            }
                        }
                        src += pixelSize;
                        ++p;

                        // check left pixel and three upper pixels for the rest of pixels
                        for (int x = 1; x < this.ImageWidth - 1; x++, src += pixelSize, p++)
                        {
                            if ((src[RGB.R] > this.backgroundThresholdR) ||
                                 (src[RGB.G] > this.backgroundThresholdG) ||
                                 (src[RGB.B] > this.backgroundThresholdB))
                            {
                                // check surrounding pixels
                                if ((src[RGB.R - pixelSize] > this.backgroundThresholdR) ||
                                     (src[RGB.G - pixelSize] > this.backgroundThresholdG) ||
                                     (src[RGB.B - pixelSize] > this.backgroundThresholdB))
                                {
                                    // label current pixel, as the left
                                    this.ObjectLabels[p] = this.ObjectLabels[p - 1];
                                }
                                else if ((src[RGB.R - strideP1] > this.backgroundThresholdR) ||
                                          (src[RGB.G - strideP1] > this.backgroundThresholdG) ||
                                          (src[RGB.B - strideP1] > this.backgroundThresholdB))
                                {
                                    // label current pixel, as the above left
                                    this.ObjectLabels[p] = this.ObjectLabels[p - 1 - this.ImageWidth];
                                }
                                else if ((src[RGB.R - stride] > this.backgroundThresholdR) ||
                                          (src[RGB.G - stride] > this.backgroundThresholdG) ||
                                          (src[RGB.B - stride] > this.backgroundThresholdB))
                                {
                                    // label current pixel, as the above
                                    this.ObjectLabels[p] = this.ObjectLabels[p - this.ImageWidth];
                                }

                                if ((src[RGB.R - strideM1] > this.backgroundThresholdR) ||
                                     (src[RGB.G - strideM1] > this.backgroundThresholdG) ||
                                     (src[RGB.B - strideM1] > this.backgroundThresholdB))
                                {
                                    if (this.ObjectLabels[p] == 0)
                                    {
                                        // label current pixel, as the above right
                                        this.ObjectLabels[p] = this.ObjectLabels[p + 1 - this.ImageWidth];
                                    }
                                    else
                                    {
                                        int l1 = this.ObjectLabels[p];
                                        int l2 = this.ObjectLabels[p + 1 - this.ImageWidth];

                                        if ((l1 != l2) && (map[l1] != map[l2]))
                                        {
                                            // merge
                                            if (map[l1] == l1)
                                            {
                                                // map left value to the right
                                                map[l1] = map[l2];
                                            }
                                            else if (map[l2] == l2)
                                            {
                                                // map right value to the left
                                                map[l2] = map[l1];
                                            }
                                            else
                                            {
                                                // both values already mapped
                                                map[map[l1]] = map[l2];
                                                map[l1] = map[l2];
                                            }

                                            // reindex
                                            for (int i = 1; i <= labelsCount; i++)
                                            {
                                                if (map[i] != i)
                                                {
                                                    // reindex
                                                    int j = map[i];
                                                    while (j != map[j])
                                                    {
                                                        j = map[j];
                                                    }
                                                    map[i] = j;
                                                }
                                            }
                                        }
                                    }
                                }

                                // label the object if it is not yet
                                if (this.ObjectLabels[p] == 0)
                                {
                                    // create new label
                                    this.ObjectLabels[p] = ++labelsCount;
                                }
                            }
                        }

                        // for the last pixel of the row, we need to check
                        // only upper and upper-left pixels
                        if ((src[RGB.R] > this.backgroundThresholdR) ||
                             (src[RGB.G] > this.backgroundThresholdG) ||
                             (src[RGB.B] > this.backgroundThresholdB))
                        {
                            // check surrounding pixels
                            if ((src[RGB.R - pixelSize] > this.backgroundThresholdR) ||
                                 (src[RGB.G - pixelSize] > this.backgroundThresholdG) ||
                                 (src[RGB.B - pixelSize] > this.backgroundThresholdB))
                            {
                                // label current pixel, as the left
                                this.ObjectLabels[p] = this.ObjectLabels[p - 1];
                            }
                            else if ((src[RGB.R - strideP1] > this.backgroundThresholdR) ||
                                      (src[RGB.G - strideP1] > this.backgroundThresholdG) ||
                                      (src[RGB.B - strideP1] > this.backgroundThresholdB))
                            {
                                // label current pixel, as the above left
                                this.ObjectLabels[p] = this.ObjectLabels[p - 1 - this.ImageWidth];
                            }
                            else if ((src[RGB.R - stride] > this.backgroundThresholdR) ||
                                      (src[RGB.G - stride] > this.backgroundThresholdG) ||
                                      (src[RGB.B - stride] > this.backgroundThresholdB))
                            {
                                // label current pixel, as the above
                                this.ObjectLabels[p] = this.ObjectLabels[p - this.ImageWidth];
                            }
                            else
                            {
                                // create new label
                                this.ObjectLabels[p] = ++labelsCount;
                            }
                        }
                        src += pixelSize;
                        ++p;

                        src += offset;
                    }
                }
            }

            // allocate remapping array
            int[] reMap = new int[map.Length];

            // count objects and prepare remapping array
            this.ObjectsCount = 0;
            for (int i = 1; i <= labelsCount; i++)
            {
                if (map[i] == i)
                {
                    // increase objects count
                    reMap[i] = ++this.ObjectsCount;
                }
            }
            // second pass to complete remapping
            for (int i = 1; i <= labelsCount; i++)
            {
                if (map[i] != i)
                    reMap[i] = reMap[map[i]];
            }

            // repair object labels
            for (int i = 0, n = this.ObjectLabels.Length; i < n; i++)
            {
                this.ObjectLabels[i] = reMap[this.ObjectLabels[i]];
            }
        }
    }
}
