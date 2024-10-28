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

namespace FileFormat.Accord.Imaging.Filters
{
    using System.Collections.Generic;
    using System.Drawing;
    using System.Drawing.Imaging;
    using System.Linq;
    using AForge.Imaging;
    using AForge.Imaging.Filters.Base_classes;
    using Core.AForge.Core;
    using global::Accord;
    using Interest_Points.Base;

    /// <summary>
    ///   Filter to mark (highlight) points in a image.
    /// </summary>
    /// 
    /// <remarks>
    /// <para>The filter highlights points on the image using a given set of points.</para>
    /// 
    /// <para>The filter accepts 8 bpp grayscale, 24 and 32 bpp color images for processing.</para>
    /// </remarks>
    /// 
    /// <example>
    /// <para>Sample usage:
    /// <code>
    /// // Create a blob contour's instance
    /// BlobCounter bc = new BlobCounter(image);
    /// 
    /// // Extract blobs
    /// Blob[] blobs = bc.GetObjectsInformation();
    /// bc.ExtractBlobsImage(bmp, blobs[0], true);
    /// 
    /// // Extract blob's edge points
    /// List&lt;IntPoint&gt; contour = bc.GetBlobsEdgePoints(blobs[0]);
    /// 
    /// // Create a green, 2 pixel width points marker's instance
    /// PointsMarker marker = new PointsMarker(contour, Color.Green, 2);
    /// 
    /// // Apply the filter in a given color image
    /// marker.ApplyInPlace(colorBlob);
    /// </code>
    /// </para>
    /// </example>
    //
    public class PointsMarker : BaseInPlaceFilter
    {
        private int width = 3;
        private Color markerColor = Color.White;
        private IEnumerable<IntPoint> points;
        private Dictionary<PixelFormat, PixelFormat> formatTranslations = new Dictionary<PixelFormat, PixelFormat>();
        private bool connect;

        /// <summary>
        ///   Format translations dictionary.
        /// </summary>
        /// 
        public override Dictionary<PixelFormat, PixelFormat> FormatTranslations
        {
            get { return this.formatTranslations; }
        }

        /// <summary>
        ///   Color used to mark corners.
        /// </summary>
        /// 
        public Color MarkerColor
        {
            get { return this.markerColor; }
            set { this.markerColor = value; }
        }

        /// <summary>
        ///   Gets or sets the set of points to mark.
        /// </summary>
        /// 
        public IEnumerable<IntPoint> Points
        {
            get { return this.points; }
            set { this.points = value; }
        }

        /// <summary>
        /// Gets or sets a value indicating whether the points 
        /// should be connected by line segments. Default is false.
        /// </summary>
        /// 
        public bool Connect
        {
            get { return this.connect; }
            set { this.connect = value; }
        }

        /// <summary>
        ///   Gets or sets the width of the points to be drawn.
        /// </summary>
        /// 
        public int Width
        {
            get { return this.width; }
            set { this.width = value; }
        }

        /// <summary>
        ///   Initializes a new instance of the <see cref="PointsMarker"/> class.
        /// </summary>
        /// 
        public PointsMarker(IEnumerable<IFeaturePoint> points)
            : this(points, Color.White, 3)
        {
        }


        /// <summary>
        ///   Initializes a new instance of the <see cref="PointsMarker"/> class.
        /// </summary>
        /// 
        public PointsMarker(IEnumerable<IFeaturePoint> points, Color markerColor)
            : this(points, markerColor, 3)
        {
        }

        /// <summary>
        ///   Initializes a new instance of the <see cref="PointsMarker"/> class.
        /// </summary>
        /// 
        public PointsMarker(IEnumerable<IFeaturePoint> points, Color markerColor, int width)
        {
            markerColor = this.init(points.Select(x => new IntPoint((int)x.X, (int)x.Y)), markerColor, width);
        }

        /// <summary>
        ///   Initializes a new instance of the <see cref="PointsMarker"/> class.
        /// </summary>
        /// 
        public PointsMarker(IEnumerable<IntPoint> points)
            : this(points, Color.White, 3)
        {
        }


        /// <summary>
        ///   Initializes a new instance of the <see cref="PointsMarker"/> class.
        /// </summary>
        /// 
        public PointsMarker(IEnumerable<IntPoint> points, Color markerColor)
            : this(points, markerColor, 3)
        {
        }

        /// <summary>
        ///   Initializes a new instance of the <see cref="PointsMarker"/> class.
        /// </summary>
        /// 
        public PointsMarker(IEnumerable<IntPoint> points, Color markerColor, int width)
        {
            markerColor = this.init(points, markerColor, width);
        }

        private Color init(IEnumerable<IntPoint> points, Color markerColor, int width)
        {
            this.points = points;
            this.markerColor = markerColor;
            this.width = width;

            this.formatTranslations[PixelFormat.Format8bppIndexed] = PixelFormat.Format8bppIndexed;
            this.formatTranslations[PixelFormat.Format24bppRgb] = PixelFormat.Format24bppRgb;
            this.formatTranslations[PixelFormat.Format32bppArgb] = PixelFormat.Format32bppArgb;
            return markerColor;
        }

        /// <summary>
        ///   Initializes a new instance of the <see cref="PointsMarker"/> class.
        /// </summary>
        /// 
        public PointsMarker(Color markerColor)
            : this((IEnumerable<IntPoint>)null, markerColor, 3)
        {
        }

        /// <summary>
        ///   Initializes a new instance of the <see cref="PointsMarker"/> class.
        /// </summary>
        /// 
        public PointsMarker(Color markerColor, int width)
        {
            this.markerColor = markerColor;
            this.width = width;

            this.formatTranslations[PixelFormat.Format8bppIndexed] = PixelFormat.Format8bppIndexed;
            this.formatTranslations[PixelFormat.Format24bppRgb] = PixelFormat.Format24bppRgb;
            this.formatTranslations[PixelFormat.Format32bppArgb] = PixelFormat.Format32bppArgb;
        }

        /// <summary>
        ///   Process the filter on the specified image.
        /// </summary>
        /// 
        /// <param name="image">Source image data.</param>
        ///
        protected override unsafe void ProcessFilter(UnmanagedImage image)
        {
            if (this.Connect)
            {
                Drawing.Polygon(image, this.points.ToList(), this.markerColor);
            }
            else
            {
                // mark all points
                foreach (IntPoint p in this.points)
                {
                    Drawing.FillRectangle(image, new Rectangle(p.X - this.width / 2, p.Y - this.width / 2, this.width, this.width), this.markerColor);
                }
            }
        }
    }
}