﻿// Accord Unit Tests
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

namespace Openize.Accord.Tests.Imaging
{
    using System.Collections.Generic;
    using NUnit.Framework;
    using Image = Openize.Accord.Imaging.AForge.Imaging.Image;
    using Openize.Accord.Tests.Imaging.Properties;
    using Openize.Accord.Core.AForge.Core;
    using Openize.Accord.Imaging.AForge.Imaging;
    using Openize.Accord.Imaging.Interest_Points;

    [TestFixture]
    public class HarrisCornersDetectorTest
    {

        [Test]
        public void ProcessImageTest()
        {
            UnmanagedImage image = UnmanagedImage.FromManagedImage(Image.Clone(Resources.image1));

            HarrisCornersDetector target = new HarrisCornersDetector(0.04f, 1000f, 1.4);
            target.Suppression = 1;

            List<IntPoint> actual = target.ProcessImage(image);

         /*   
                        PointsMarker marker = new PointsMarker(actual.ToArray());
                        marker.Width = 1;
                        marker.MarkerColor = Color.FromArgb(128, 128, 128);
                        var markers = marker.Apply(image);
                        ImageBox.Show(markers.ToManagedImage(), PictureBoxSizeMode.Zoom);
           */ 

            /*
                        Assert.AreEqual(4, actual.Count);
                        Assert.IsTrue(actual.Contains(new IntPoint(3, 3)));
                        Assert.IsTrue(actual.Contains(new IntPoint(14, 3)));
                        Assert.IsTrue(actual.Contains(new IntPoint(3, 14)));
                        Assert.IsTrue(actual.Contains(new IntPoint(14, 14)));
            */

            Assert.AreEqual(4, actual.Count);
            Assert.IsTrue(actual.Contains(new IntPoint(3, 3)));
            Assert.IsTrue(actual.Contains(new IntPoint(12, 3)));
            Assert.IsTrue(actual.Contains(new IntPoint(3, 12)));
            Assert.IsTrue(actual.Contains(new IntPoint(12, 12)));
        }


        [Test]
        public void ProcessImageTest2()
        {
            UnmanagedImage image = UnmanagedImage.FromManagedImage(Image.Clone(Resources.sample_black));

            HarrisCornersDetector target = new HarrisCornersDetector(HarrisCornerMeasure.Noble, 700f, 1.4, 1);

            List<IntPoint> actual = target.ProcessImage(image);
            
            /*
            PointsMarker marker = new PointsMarker(actual.ToArray());
            marker.Width = 3;
            marker.MarkerColor = Color.FromArgb(255, 0, 0);
            var markers = marker.Apply(image);
            ImageBox.Show(markers.ToManagedImage(), PictureBoxSizeMode.Zoom);
             */

            Assert.AreEqual(actual.Count, 42);
        }
 
    }
}
