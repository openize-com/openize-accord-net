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

namespace Openize.Accord.Imaging.Moments
{
    using System;
    using System.Drawing;
    using System.Drawing.Imaging;
    using AForge.Imaging;
    using Base;
    using Colors;

    /// <summary>
    ///   Raw image moments.
    /// </summary>
    /// 
    /// <remarks>
    /// <para>
    ///   In image processing, computer vision and related fields, an image moment is
    ///   a certain particular weighted average (moment) of the image pixels' intensities,
    ///   or a function of such moments, usually chosen to have some attractive property 
    ///   or interpretation.</para>
    ///
    /// <para>
    ///   Image moments are useful to describe objects after segmentation. Simple properties 
    ///   of the image which are found via image moments include area (or total intensity), 
    ///   its centroid, and information about its orientation.</para>
    ///   
    /// <para>
    ///   The raw moments are the most basic moments which can be computed from an image,
    ///   and can then be further processed to achieve <see cref="CentralMoments"/> or even
    ///   <see cref="HuMoments"/>.</para>
    ///   
    /// <para>
    ///   References:
    ///   <list type="bullet">
    ///     <item><description>
    ///       Wikipedia contributors. "Image moment." Wikipedia, The Free Encyclopedia. Wikipedia,
    ///       The Free Encyclopedia. Available at http://en.wikipedia.org/wiki/Image_moment </description></item>
    ///   </list>
    /// </para>
    /// </remarks>
    /// 
    /// <example>
    /// <code>
    /// Bitmap image = ...;
    ///
    /// // Compute the raw moments of up to third order
    /// RawMoments m = new RawMoments(image, order: 3);
    /// </code>
    /// </example>
    /// 
    /// <seealso cref="HuMoments"/>
    /// <seealso cref="CentralMoments"/>
    /// 
    [Serializable]
    public class RawMoments : MomentsBase, IMoments
    {

        /// <summary>
        ///   Gets the default maximum moment order.
        /// </summary>
        /// 
        public const int DefaultOrder = 3;


        /// <summary>
        ///   Raw moment of order (0,0).
        /// </summary>
        /// 
        public float M00 { get; private set; }

        /// <summary>
        ///   Raw moment of order (1,0).
        /// </summary>
        /// 
        public float M10 { get; private set; }

        /// <summary>
        ///   Raw moment of order (0,1).
        /// </summary>
        /// 
        public float M01 { get; private set; }

        /// <summary>
        ///   Raw moment of order (1,1).
        /// </summary>
        /// 
        public float M11 { get; private set; }

        /// <summary>
        ///   Raw moment of order (2,0).
        /// </summary>
        /// 
        public float M20 { get; private set; }

        /// <summary>
        ///   Raw moment of order (0,2).
        /// </summary>
        /// 
        public float M02 { get; private set; }

        /// <summary>
        ///   Raw moment of order (2,1).
        /// </summary>
        /// 
        public float M21 { get; private set; }

        /// <summary>
        ///   Raw moment of order (1,2).
        /// </summary>
        /// 
        public float M12 { get; private set; }

        /// <summary>
        ///   Raw moment of order (3,0).
        /// </summary>
        /// 
        public float M30 { get; private set; }

        /// <summary>
        ///   Raw moment of order (0,3).
        /// </summary>
        /// 
        public float M03 { get; private set; }


        /// <summary>
        ///   Inverse raw moment of order (0,0).
        /// </summary>
        /// 
        public float InvM00 { get; private set; }



        /// <summary>
        ///   Gets the X centroid of the image.
        /// </summary>
        /// 
        public float CenterX { get; private set; }

        /// <summary>
        ///   Gets the Y centroid of the image.
        /// </summary>
        /// 
        public float CenterY { get; private set; }

        /// <summary>
        ///   Gets the area (for binary images) or sum of
        ///   gray level (for grayscale images).
        /// </summary>
        /// 
        public float Area { get { return this.M00; } }


        /// <summary>
        ///   Initializes a new instance of the <see cref="Moments"/> class.
        /// </summary>
        /// 
        public RawMoments(int order = DefaultOrder)
            : base(order) { }

        /// <summary>
        ///   Initializes a new instance of the <see cref="Moments"/> class.
        /// </summary>
        /// 
        public RawMoments(Bitmap image, Rectangle area, int order = DefaultOrder)
            : base(order)
        {
            this.Compute(image, area);
        }

        /// <summary>
        ///   Initializes a new instance of the <see cref="Moments"/> class.
        /// </summary>
        /// 
        public RawMoments(UnmanagedImage image, Rectangle area, int order = DefaultOrder)
            : base(image, area, order) { }

        /// <summary>
        ///   Initializes a new instance of the <see cref="Moments"/> class.
        /// </summary>
        /// 
        public RawMoments(float[,] image, Rectangle area, int order = DefaultOrder)
            : base(image, area, order) { }

        /// <summary>
        ///   Initializes a new instance of the <see cref="Moments"/> class.
        /// </summary>
        /// 
        public RawMoments(UnmanagedImage image, int order = DefaultOrder)
            : base(image, order) { }

        /// <summary>
        ///   Initializes a new instance of the <see cref="Moments"/> class.
        /// </summary>
        /// 
        public RawMoments(Bitmap image, int order = DefaultOrder)
            : base(image, order) { }

        /// <summary>
        ///   Initializes a new instance of the <see cref="Moments"/> class.
        /// </summary>
        /// 
        public RawMoments(float[,] image, int order = DefaultOrder)
            : base(image, order) { }



        /// <summary>
        ///   Computes the raw moments for the specified image.
        /// </summary>
        /// 
        /// <param name="image">The image whose moments should be computed.</param>
        /// <param name="area">The region of interest in the image to compute moments for.</param>
        /// <param name="secondOrder"><c>True</c> to compute second order moments, <c>false</c> otherwise.</param>
        /// 
        [Obsolete("Use the Order property to determine the maximum order to be computed.")]
        public void Compute(float[,] image, Rectangle area, bool secondOrder)
        {
            this.Order = secondOrder ? 2 : 1;
            this.Compute(image, area);
        }

        /// <summary>
        ///   Computes the raw moments for the specified image.
        /// </summary>
        /// 
        /// <param name="image">The image.</param>
        /// <param name="area">The region of interest in the image to compute moments for.</param>
        /// 
        public override void Compute(float[,] image, Rectangle area)
        {
            int height = image.GetLength(0);
            int width = image.GetLength(1);

            if (area == Rectangle.Empty)
                area = new Rectangle(0, 0, width, height);

            this.Reset();

            // stay inside the image
            int windowX = Math.Max(area.X, 0);
            int windowY = Math.Max(area.Y, 0);
            int windowWidth = Math.Min(windowX + area.Width, width);
            int windowHeight = Math.Min(windowY + area.Height, height);

            int offset = width - (windowWidth - windowX);

            unsafe
            {
                fixed (float* ptrImage = image)
                {
                    float* src = ptrImage + windowY * width + windowX;

                    // TODO: Walk using x and y directly instead of i and j.

                    for (int j = windowY; j < windowHeight; j++)
                    {
                        float y = j - windowY;

                        for (int i = windowX; i < windowWidth; i++, src++)
                        {
                            float x = i - windowX;

                            float v = *src;

                            this.M00 += v;
                            this.M01 += y * v;
                            this.M10 += x * v;

                            if (this.Order >= 2)
                            {
                                this.M11 += x * y * v;
                                this.M02 += y * y * v;
                                this.M20 += x * x * v;
                            }

                            if (this.Order >= 3)
                            {
                                this.M12 += x * y * y * v;
                                this.M21 += x * x * y * v;

                                this.M30 += x * x * x * v;
                                this.M03 += y * y * y * v;
                            }
                        }

                        src += offset;
                    }
                }
            }

            this.InvM00 = 1f / this.M00;
            this.CenterX = this.M10 * this.InvM00;
            this.CenterY = this.M01 * this.InvM00;
        }

        /// <summary>
        ///   Computes the raw moments for the specified image.
        /// </summary>
        /// 
        /// <param name="image">The image.</param>
        /// <param name="area">The region of interest in the image to compute moments for.</param>
        /// 
        public override void Compute(UnmanagedImage image, Rectangle area)
        {
            int height = image.Height;
            int width = image.Width;
            int stride = image.Stride;

            if (area == Rectangle.Empty)
                area = new Rectangle(0, 0, width, height);

            this.Reset();


            // stay inside the image
            int windowX = Math.Max(area.X, 0);
            int windowY = Math.Max(area.Y, 0);
            int windowWidth = Math.Min(windowX + area.Width, width);
            int windowHeight = Math.Min(windowY + area.Height, height);

            unsafe
            {
                if (image.PixelFormat == PixelFormat.Format8bppIndexed)
                {
                    int offset = stride - (windowWidth - windowX);

                    byte* src = (byte*)image.ImageData.ToPointer() + windowY * stride + windowX;

                    for (int j = windowY; j < windowHeight; j++)
                    {
                        float y = j - windowY;

                        for (int i = windowX; i < windowWidth; i++, src++)
                        {
                            float x = i - windowX;

                            float v = *src;

                            this.M00 += v;
                            this.M01 += y * v;
                            this.M10 += x * v;

                            if (this.Order >= 2)
                            {
                                this.M11 += x * y * v;
                                this.M02 += y * y * v;
                                this.M20 += x * x * v;
                            }

                            if (this.Order >= 3)
                            {
                                this.M12 += x * y * y * v;
                                this.M21 += x * x * y * v;

                                this.M30 += x * x * x * v;
                                this.M03 += y * y * y * v;
                            }
                        }

                        src += offset;
                    }
                }
                else
                {
                    // color images
                    int pixelSize = Bitmap.GetPixelFormatSize(image.PixelFormat) / 8;
                    int offset = stride - (windowWidth - windowX) * pixelSize;

                    byte* src = (byte*)image.ImageData.ToPointer() + windowY * stride + windowX * pixelSize;

                    for (int j = windowY; j < windowHeight; j++)
                    {
                        float y = j - windowY;

                        for (int i = windowX; i < windowWidth; i++, src += pixelSize)
                        {
                            float x = i - windowX;

                            // BT709 - 0.2125, 0.7154, 0.0721 
                            float v = (float)(0.2125 * src[RGB.R] + 0.7154 * src[RGB.G] + 0.0721 * src[RGB.B]);

                            this.M00 += v;
                            this.M01 += y * v;
                            this.M10 += x * v;

                            if (this.Order >= 2)
                            {
                                this.M11 += x * y * v;
                                this.M02 += y * y * v;
                                this.M20 += x * x * v;
                            }

                            if (this.Order >= 3)
                            {
                                this.M12 += x * y * y * v;
                                this.M21 += x * x * y * v;

                                this.M30 += x * x * x * v;
                                this.M03 += y * y * y * v;
                            }
                        }

                        src += offset;
                    }
                }
            }


            this.InvM00 = 1f / this.M00;
            this.CenterX = this.M10 * this.InvM00;
            this.CenterY = this.M01 * this.InvM00;
        }


        /// <summary>
        ///   Resets all moments to zero.
        /// </summary>
        /// 
        protected void Reset()
        {
            this.M00 = this.M10 = this.M01 = 0;
            this.M11 = this.M20 = this.M02 = 0;

            this.M21 = this.M12 = 0;
            this.M30 = this.M03 = 0;

            this.InvM00 = 0;

            this.CenterX = this.CenterY = 0;
        }

    }
}
