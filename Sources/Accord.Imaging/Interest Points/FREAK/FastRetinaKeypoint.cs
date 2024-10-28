// Accord Imaging Library
// The Accord.NET Framework
// http://accord-framework.net
//
// Copyright © Christopher Evans, 2009-2011
// http://www.chrisevansdev.com/
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

namespace FileFormat.Accord.Imaging.Interest_Points.FREAK
{
    using System;
    using System.Drawing;
    using System.Text;
    using Base;
    using Core.AForge.Core;
    using Core.Features;
    using global::Accord;
    using Point = System.Drawing.Point;

    /// <summary>
    ///   Fast Retina Keypoint (FREAK) point.
    /// </summary>
    /// 
    /// <remarks>
    ///   In order to extract feature points from an image using FREAK,
    ///   please take a look on the <see cref="FastRetinaKeypointDetector"/>
    ///   documentation page.
    /// </remarks>
    /// 
    /// <seealso cref="FastRetinaKeypointDescriptor"/>
    /// <seealso cref="FastRetinaKeypointDetector"/>
    /// 
    [Serializable]
    public class FastRetinaKeypoint : IFeaturePoint<byte[]>, IFeaturePoint<double[]>
    {

        /// <summary>
        ///   Initializes a new instance of the <see cref="FastRetinaKeypoint"/> class.
        /// </summary>
        /// 
        /// <param name="x">The x-coordinate of the point in the image.</param>
        /// <param name="y">The y-coordinate of the point in the image.</param>
        /// 
        public FastRetinaKeypoint(double x, double y)
        {
            this.X = x;
            this.Y = y;
            this.Scale = 1;
        }

        /// <summary>
        ///   Gets or sets the x-coordinate of this point.
        /// </summary>
        /// 
        public double X { get; set; }

        /// <summary>
        ///   Gets or sets the y-coordinate of this point.
        /// </summary>
        /// 
        public double Y { get; set; }

        /// <summary>
        ///   Gets or sets the scale of the point.
        /// </summary>
        /// 
        public double Scale { get; set; }

        /// <summary>
        ///   Gets or sets the orientation of this point in angles.
        /// </summary>
        /// 
        public double Orientation { get; set; }


        /// <summary>
        ///   Gets or sets the descriptor vector
        ///   associated with this point.
        /// </summary>
        /// 
        public byte[] Descriptor { get; set; }


        /// <summary>
        ///   Converts the binary descriptor to 
        ///   string of hexadecimal values.
        /// </summary>
        /// 
        /// <returns>A string containing an hexadecimal
        /// value representing this point's descriptor.</returns>
        /// 
        public string ToHex()
        {
            if (this.Descriptor == null)
                return String.Empty;

            StringBuilder hex = new StringBuilder(this.Descriptor.Length * 2);
            foreach (byte b in this.Descriptor)
                hex.AppendFormat("{0:x2}", b);
            return hex.ToString();
        }

        /// <summary>
        ///   Converts the binary descriptor
        ///   to a string of binary values.
        /// </summary>
        /// 
        /// <returns>A string containing a binary value
        /// representing this point's descriptor.</returns>
        /// 
        public string ToBinary()
        {
            StringBuilder sb = new StringBuilder();

            for (int i = 0; i < this.Descriptor.Length; i++)
            {
                for (int j = 0; j < 8; j++)
                {
                    bool set = (this.Descriptor[i] & (1 << j)) != 0;
                    sb.Append(set ? "1" : "0");
                }
            }

            return sb.ToString();
        }

        /// <summary>
        ///   Converts the binary descriptor to base64.
        /// </summary>
        /// 
        /// <returns>A string containing the base64 
        /// representation of the descriptor.</returns>
        /// 
        public string ToBase64()
        {
            if (this.Descriptor == null)
                return String.Empty;

            return Convert.ToBase64String(this.Descriptor);
        }




        /// <summary>
        ///   Converts the feature point to a <see cref="IntPoint"/>.
        /// </summary>
        /// 
        public IntPoint ToIntPoint()
        {
            return new IntPoint((int)this.X, (int)this.Y);
        }

        /// <summary>
        ///   Converts this object into a <see cref="IntPoint"/>.
        /// </summary>
        /// 
        /// <returns>
        ///   The result of the conversion.
        /// </returns>
        /// 
        public Point ToPoint()
        {
            return new Point((int)this.X, (int)this.Y);
        }

        /// <summary>
        ///   Converts this object into a <see cref="System.Drawing.PointF"/>.
        /// </summary>
        /// 
        /// <returns>
        ///   The result of the conversion.
        /// </returns>
        /// 
        public PointF ToPointF()
        {
            return new PointF((float)this.X, (float)this.Y);
        }

        /// <summary>
        ///   Performs an implicit conversion from <see cref="FastRetinaKeypoint"/>
        ///   to <see cref="System.Drawing.Point"/>.
        /// </summary>
        /// 
        /// <param name="point">The point to be converted.</param>
        /// 
        /// <returns>
        ///   The result of the conversion.
        /// </returns>
        /// 
        public static implicit operator Point(FastRetinaKeypoint point)
        {
            return point.ToPoint();
        }

        /// <summary>
        ///   Performs an implicit conversion from <see cref="FastRetinaKeypoint"/>
        ///   to <see cref="System.Drawing.PointF"/>.
        /// </summary>
        /// 
        /// <param name="point">The point to be converted.</param>
        /// 
        /// <returns>
        ///   The result of the conversion.
        /// </returns>
        /// 
        public static implicit operator PointF(FastRetinaKeypoint point)
        {
            return point.ToPointF();
        }

        /// <summary>
        ///   Performs an implicit conversion from <see cref="FastRetinaKeypoint"/>
        ///   to <see cref="IntPoint"/>.
        /// </summary>
        /// 
        /// <param name="point">The point to be converted.</param>
        /// 
        /// <returns>
        ///   The result of the conversion.
        /// </returns>
        /// 
        public static implicit operator global::FileFormat.Accord.Core.AForge.Core.IntPoint(FastRetinaKeypoint point)
        {
            return point.ToIntPoint();
        }

        double[] IFeatureDescriptor<double[]>.Descriptor
        {
            get
            {
                double[] r = new double[this.Descriptor.Length];
                for (int i = 0; i < this.Descriptor.Length; i++)
                    r[i] = this.Descriptor[i];
                return r;
            }
        }
    }
}
