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

namespace Openize.Accord.Imaging.Interest_Points.Base
{
    using System;
    using System.Collections.Generic;
    using System.Drawing;
    using System.Drawing.Imaging;
    using System.Linq;
    using AForge.Imaging;
    using Openize.Accord.Core.AForge.Core;

    /// <summary>
    ///   Base class for image feature extractors that implement the <see cref="IImageFeatureExtractor{TPoint}"/> interface.
    /// </summary>
    /// 
    /// <typeparam name="TPoint">The type of the descriptor vector for the feature (e.g. double[]).</typeparam>
    /// 
    /// <seealso cref="IImageFeatureExtractor{TFeature}" />
    /// <seealso cref="BaseSparseFeatureExtractor{T}"/>
    /// 
    [Serializable]
    public abstract class BaseSparseFeatureExtractor<TPoint> : BaseFeatureExtractor<TPoint>, ICornersDetector
        where TPoint : IFeaturePoint<double[]>
    {

        /// <summary>
        ///   Obsolete. Please use the <see cref="BaseFeatureExtractor{T}.Transform(Bitmap)"/> method instead.
        /// </summary>
        /// 
        [Obsolete("Please use the Transform(input) method instead.")]
        public new List<TPoint> ProcessImage(Bitmap image)
        {
            return new List<TPoint>(this.Transform(image));
        }

        /// <summary>
        ///   Obsolete. Please use the <see cref="BaseFeatureExtractor{T}.Transform(Bitmap)"/> method instead.
        /// </summary>
        /// 
        [Obsolete("Please use the Transform(input) method instead.")]
        public new List<TPoint> ProcessImage(BitmapData imageData)
        {
            return new List<TPoint>(this.Transform(new UnmanagedImage(imageData)));
        }

        /// <summary>
        ///   Obsolete. Please use the <see cref="BaseFeatureExtractor{T}.Transform(Bitmap)"/> method instead.
        /// </summary>
        /// 
        [Obsolete("Please use the Transform(input) method instead.")]
        public new List<TPoint> ProcessImage(UnmanagedImage image)
        {
            return new List<TPoint>(this.Transform(image));
        }


        List<IntPoint> ICornersDetector.ProcessImage(Bitmap image)
        {
            return this.Transform(image).Select(x => new IntPoint((int)x.X, (int)x.Y)).ToList();
        }

        List<IntPoint> ICornersDetector.ProcessImage(BitmapData imageData)
        {
            return this.Transform(new UnmanagedImage(imageData)).Select(x => new IntPoint((int)x.X, (int)x.Y)).ToList();
        }

        List<IntPoint> ICornersDetector.ProcessImage(UnmanagedImage image)
        {
            return this.Transform(image).Select(x => new IntPoint((int)x.X, (int)x.Y)).ToList();
        }
     
    }
}