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

namespace Openize.Accord.Imaging.Interest_Points
{
    using System;
    using Base;
    using Openize.Accord.Core.AForge.Core;

    /// <summary>
    ///   Corner feature point.
    /// </summary>
    /// 
    [Serializable]
    public class CornerFeaturePoint : IFeaturePoint<double[]>
    {

        /// <summary>
        ///   Initializes a new instance of the <see cref="CornerFeaturePoint"/> class.
        /// </summary>
        /// 
        public CornerFeaturePoint(IntPoint point)
        {
            this.Descriptor = new double[] { point.X, point.Y };
        }

        /// <summary>
        ///   Gets the X position of the point.
        /// </summary>
        /// 
        public double X { get { return this.Descriptor[0]; } }

        /// <summary>
        ///   Gets the Y position of the point.
        /// </summary>
        /// 
        public double Y { get { return this.Descriptor[1]; } }

        /// <summary>
        ///   Gets the descriptor vector
        ///   associated with this point.
        /// </summary>
        /// 
        public double[] Descriptor { get; private set; }
    }
}
