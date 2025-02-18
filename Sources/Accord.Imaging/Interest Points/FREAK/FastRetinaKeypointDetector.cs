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

namespace Openize.Accord.Imaging.Interest_Points.FREAK
{
    using System;
    using System.Collections.Generic;
    using System.Drawing.Imaging;
    using AForge.Imaging;
    using AForge.Imaging.Filters.Color_Filters;
    using Base;
    using Openize.Accord.Core.AForge.Core;

    /// <summary>
    ///   SURF Feature descriptor types.
    /// </summary>
    /// 
    public enum FastRetinaKeypointDescriptorType
    {
        /// <summary>
        ///   Do not compute descriptors.
        /// </summary>
        /// 
        None,

        /// <summary>
        ///   Compute standard 512-bit descriptors.
        /// </summary>
        /// 
        Standard,

        /// <summary>
        ///   Compute extended 1024-bit descriptors.
        /// </summary>
        /// 
        Extended,
    }

    /// <summary>
    ///   Fast Retina Keypoint (FREAK) detector.
    /// </summary>
    /// 
    /// <remarks>
    ///   The FREAK algorithm is a binary based interest point descriptor algorithm
    ///   that relies in another corner 
    /// </remarks>
    /// 
    /// <example>
    /// <para>
    ///   In the following example, we will see how can we extract binary descriptor
    ///   vectors from a given image using the Fast Retina Keypoint Detector together
    ///   a FAST corners detection algorithm.</para>
    /// 
    /// <code>
    /// Bitmap lena = Resources.lena512;
    /// 
    /// // The freak detector can be used with any other corners detection
    /// // algorithm. The default corners detection method used is the FAST
    /// // corners detection. So, let's start creating this detector first:
    /// //
    /// var detector = new FastCornersDetector(60);
    /// 
    /// // Now that we have a corners detector, we can pass it to the FREAK
    /// // feature extraction algorithm. Please note that if we leave this
    /// // parameter empty, FAST will be used by default.
    /// //
    /// var freak = new FastRetinaKeypointDetector(detector);
    /// 
    /// // Now, all we have to do is to process our image:
    /// List&lt;FastRetinaKeypoint> points = freak.ProcessImage(lena);
    /// 
    /// // Afterwards, we should obtain 83 feature points. We can inspect
    /// // the feature points visually using the FeaturesMarker class as
    /// //
    /// FeaturesMarker marker = new FeaturesMarker(points, scale: 20);
    ///             
    /// // And showing it on screen with
    /// ImageBox.Show(marker.Apply(lena));
    /// 
    /// // We can also inspect the feature vectors (descriptors) associated
    /// // with each feature point. In order to get a descriptor vector for
    /// // any given point, we can use
    /// //
    /// byte[] feature = points[42].Descriptor;
    ///             
    /// // By default, feature vectors will have 64 bytes in length. We can also
    /// // display those vectors in more readable formats such as HEX or base64
    /// //
    /// string hex = points[42].ToHex();
    /// string b64 = points[42].ToBase64();
    /// 
    /// // The above base64 result should be:
    /// //
    /// //  "3W8M/ev///ffbr/+v3f34vz//7X+f0609v//+++/1+jfq/e83/X5/+6ft3//b4uaPZf7ePb3n/P93/rIbZlf+g=="
    /// //
    /// </code>
    /// 
    /// <para>
    ///   The resulting image is shown below:</para>
    /// 
    ///   <img src="..\images\freak.png" />
    ///   
    /// </example>
    /// 
    /// <seealso cref="FastCornersDetector"/>
    /// <seealso cref="HistogramsOfOrientedGradients"/>
    /// <seealso cref="LocalBinaryPattern"/>
    /// 
    [Serializable]
    public class FastRetinaKeypointDetector : BaseSparseFeatureExtractor<FastRetinaKeypoint>
    {

        private FastRetinaKeypointDescriptorType featureType = FastRetinaKeypointDescriptorType.Standard;
        private float scale = 22.0f;
        private int octaves = 4;


        [NonSerialized]
        private IntegralImage integral;

        [NonSerialized]
        private UnmanagedImage grayImage;

        [NonSerialized]
        private FastRetinaKeypointPattern pattern;

        [NonSerialized]
        private FastRetinaKeypointDescriptor descriptor;



        /// <summary>
        ///   Gets the corners detector used to generate features.
        /// </summary>
        /// 
        public ICornersDetector Detector { get; private set; }


        /// <summary>
        ///   Gets or sets a value indicating whether all feature points
        ///   should have their descriptors computed after being detected.
        ///   Default is to compute standard descriptors.
        /// </summary>
        /// 
        /// <value><c>true</c> if to compute orientation; otherwise, <c>false</c>.</value>
        /// 
        public FastRetinaKeypointDescriptorType ComputeDescriptors
        {
            get { return this.featureType; }
            set { this.featureType = value; }
        }

        /// <summary>
        ///   Gets or sets the number of octaves to use when 
        ///   building the feature descriptor. Default is 4.
        /// </summary>
        /// 
        public int Octaves
        {
            get { return this.pattern.Octaves; }
            set
            {
                if (value != this.octaves)
                {
                    this.octaves = value;
                    this.pattern = null;
                }
            }
        }

        /// <summary>
        ///   Gets or sets the scale used when building 
        ///   the feature descriptor. Default is 22.
        /// </summary>
        /// 
        public float Scale
        {
            get { return this.pattern.Scale; }
            set
            {
                if (value != this.scale)
                {
                    this.scale = value;
                    this.pattern = null;
                }
            }
        }

        /// <summary>
        ///   Initializes a new instance of the <see cref="FastRetinaKeypointDetector"/> class.
        /// </summary>
        /// 
        /// <param name="threshold">The detection threshold for the 
        /// <see cref="FastCornersDetector">FAST detector</see>.</param>
        /// 
        public FastRetinaKeypointDetector(int threshold)
        {
            this.init(new FastCornersDetector(threshold));
        }

        /// <summary>
        ///   Initializes a new instance of the <see cref="FastRetinaKeypointDetector"/> class.
        /// </summary>
        /// 
        public FastRetinaKeypointDetector()
        {
            this.init(new FastCornersDetector());
        }

        /// <summary>
        ///   Initializes a new instance of the <see cref="FastRetinaKeypointDetector"/> class.
        /// </summary>
        /// 
        /// <param name="detector">A corners detector.</param>
        /// 
        public FastRetinaKeypointDetector(ICornersDetector detector)
        {
            this.init(detector);
        }

        private void init(ICornersDetector detector)
        {
            this.Detector = detector;

            base.SupportedFormats.UnionWith(new[]
            {
                PixelFormat.Format8bppIndexed,
                PixelFormat.Format24bppRgb,
                PixelFormat.Format32bppRgb,
                PixelFormat.Format32bppArgb,
            });
        }

        /// <summary>
        ///   This method should be implemented by inheriting classes to implement the 
        ///   actual feature extraction, transforming the input image into a list of features.
        /// </summary>
        /// 
        protected override IEnumerable<FastRetinaKeypoint> InnerTransform(UnmanagedImage image)
        {
            // make sure we have grayscale image
            if (image.PixelFormat == PixelFormat.Format8bppIndexed)
            {
                this.grayImage = image;
            }
            else
            {
                // create temporary grayscale image
                this.grayImage = Grayscale.CommonAlgorithms.BT709.Apply(image);
            }

            // 1. Extract corners points from the image.
            List<IntPoint> corners = this.Detector.ProcessImage(this.grayImage);

            var features = new List<FastRetinaKeypoint>();
            for (int i = 0; i < corners.Count; i++)
                features.Add(new FastRetinaKeypoint(corners[i].X, corners[i].Y));

            // 2. Compute the integral for the given image
            this.integral = IntegralImage.FromBitmap(this.grayImage);

            // 3. Compute feature descriptors if required
            this.descriptor = null;
            if (this.featureType != FastRetinaKeypointDescriptorType.None)
            {
                this.descriptor = this.GetDescriptor();
                this.descriptor.Compute(features);
            }

            return features;
        }

        /// <summary>
        ///   Gets the <see cref="FastRetinaKeypointDescriptor">
        ///   feature descriptor</see> for the last processed image.
        /// </summary>
        /// 
        public FastRetinaKeypointDescriptor GetDescriptor()
        {
            if (this.descriptor == null || this.pattern == null)
            {
                if (this.pattern == null)
                    this.pattern = new FastRetinaKeypointPattern(this.octaves, this.scale);

                this.descriptor = new FastRetinaKeypointDescriptor(this.grayImage, this.integral, this.pattern);
                this.descriptor.Extended = this.featureType == FastRetinaKeypointDescriptorType.Extended;
            }

            return this.descriptor;
        }

        /// <summary>
        /// Creates a new object that is a copy of the current instance.
        /// </summary>
        /// 
        protected override object Clone(ISet<PixelFormat> supportedFormats)
        {
            var clone = new FastRetinaKeypointDetector();
            clone.featureType = this.featureType;
            clone.octaves = this.octaves;
            clone.scale = this.scale;
            clone.SupportedFormats = supportedFormats;

            if (this.descriptor != null)
            {
                // clone.descriptor = (FastRetinaKeypointDescriptor)descriptor.Clone();
                // clone.grayImage = grayImage.Clone();
                // clone.integral = (IntegralImage)integral.Clone();
                // clone.pattern = (FastRetinaKeypointPattern)pattern.Clone();
            }

            return clone;
        }

    }
}
