// AForge Image Processing Library
// AForge.NET framework
// http://www.aforgenet.com/framework/
//
// Copyright � AForge.NET, 2005-2011
// contacts@aforgenet.com
//

namespace Openize.Accord.Imaging.AForge.Imaging.Filters.HSL_Filters
{
    using System;
    using System.Collections.Generic;
    using System.Drawing;
    using System.Drawing.Imaging;
    using Base_classes;
    using Colors;

    /// <summary>
    /// Hue modifier.
    /// </summary>
    /// 
    /// <remarks><para>The filter operates in <b>HSL</b> color space and updates
    /// pixels' hue values setting it to the specified value (luminance and
    /// saturation are kept unchanged). The result of the filter looks like the image
    /// is observed through a glass of the given color.</para>
    ///
    /// <para>The filter accepts 24 and 32 bpp color images for processing.</para>
    /// <para>Sample usage:</para>
    /// <code>
    /// // create filter
    /// HueModifier filter = new HueModifier( 180 );
    /// // apply the filter
    /// filter.ApplyInPlace( image );
    /// </code>
    /// 
    /// <para><b>Initial image:</b></para>
    /// <img src="..\images\imaging\sample1.jpg" width="480" height="361" />
    /// <para><b>Result image:</b></para>
    /// <img src="..\images\imaging\hue_modifier.jpg" width="480" height="361" />
    /// </remarks>
    /// 
    public class HueModifier : BaseInPlacePartialFilter
    {
        private int hue = 0;

        // private format translation dictionary
        private Dictionary<PixelFormat, PixelFormat> formatTranslations = new Dictionary<PixelFormat, PixelFormat>();

        /// <summary>
        /// Format translations dictionary.
        /// </summary>
        public override Dictionary<PixelFormat, PixelFormat> FormatTranslations
        {
            get { return this.formatTranslations; }
        }

        /// <summary>
        /// Hue value to set, [0, 359].
        /// </summary>
        /// 
        /// <remarks><para>Default value is set to <b>0</b>.</para></remarks>
        /// 
        public int Hue
        {
            get { return this.hue; }
            set { this.hue = Math.Max(0, Math.Min(359, value)); }
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="HueModifier"/> class.
        /// </summary>
        /// 
        public HueModifier()
        {
            this.formatTranslations[PixelFormat.Format24bppRgb] = PixelFormat.Format24bppRgb;
            this.formatTranslations[PixelFormat.Format32bppRgb] = PixelFormat.Format32bppRgb;
            this.formatTranslations[PixelFormat.Format32bppArgb] = PixelFormat.Format32bppArgb;
            this.formatTranslations[PixelFormat.Format32bppPArgb] = PixelFormat.Format32bppPArgb;
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="HueModifier"/> class.
        /// </summary>
        /// 
        /// <param name="hue">Hue value to set.</param>
        /// 
        public HueModifier(int hue) : this()
        {
            this.hue = hue;
        }

        /// <summary>
        /// Process the filter on the specified image.
        /// </summary>
        /// 
        /// <param name="image">Source image data.</param>
        /// <param name="rect">Image rectangle for processing by the filter.</param>
        ///
        protected override unsafe void ProcessFilter(UnmanagedImage image, Rectangle rect)
        {
            int pixelSize = Bitmap.GetPixelFormatSize(image.PixelFormat) / 8;

            int startX = rect.Left;
            int startY = rect.Top;
            int stopX = startX + rect.Width;
            int stopY = startY + rect.Height;
            int offset = image.Stride - rect.Width * pixelSize;

            RGB rgb = new RGB();
            HSL hsl = new HSL();

            // do the job
            byte* ptr = (byte*)image.ImageData.ToPointer();

            // allign pointer to the first pixel to process
            ptr += (startY * image.Stride + startX * pixelSize);

            // for each row
            for (int y = startY; y < stopY; y++)
            {
                // for each pixel
                for (int x = startX; x < stopX; x++, ptr += pixelSize)
                {
                    rgb.Red = ptr[RGB.R];
                    rgb.Green = ptr[RGB.G];
                    rgb.Blue = ptr[RGB.B];

                    // convert to HSL
                    global::Openize.Accord.Imaging.Colors.HSL.FromRGB(rgb, ref hsl);

                    // modify hue value
                    hsl.Hue = this.hue;

                    // convert back to RGB
                    global::Openize.Accord.Imaging.Colors.HSL.ToRGB(hsl, ref rgb);

                    ptr[RGB.R] = rgb.Red;
                    ptr[RGB.G] = rgb.Green;
                    ptr[RGB.B] = rgb.Blue;
                }
                ptr += offset;
            }
        }
    }
}
