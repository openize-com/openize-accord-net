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
    using System.Drawing;
    using NUnit.Framework;
    using Image = Openize.Accord.Imaging.AForge.Imaging.Image;
    using Openize.Accord.Tests.Imaging.Properties;
    using Openize.Accord.Imaging.Moments;

    [TestFixture]
    public class CentralMomentsTest
    {

        [Test]
        public void ComputeTest()
        {
            Bitmap image = Image.Clone(Resources.hu);

            CentralMoments target = new CentralMoments(image, order: 3);

            Assert.AreEqual(86424.0 / target.Mu00, 1, 1e-4);
            Assert.AreEqual(0, target.Mu01);
            Assert.AreEqual(0, target.Mu10);
            Assert.AreEqual(5.868206472635379E8 / target.Mu02, 1, 1e-2);

            Assert.AreEqual(6348920.945848465 / target.Mu11, 1, 1e-2);
            Assert.AreEqual(9.084235762166061E8 / target.Mu20, 1, 1e-3);

            Assert.AreEqual(-2.155191E9 / target.Mu12, 1, 1e-4);
            Assert.AreEqual(7.125893E8 / target.Mu21, 1, 1e-3);

            Assert.AreEqual(-1.26244547E10 / target.Mu30, 1, 1e-4);
            Assert.AreEqual(1.71818829E9 / target.Mu03, 1, 1e-4);

            SizeF size = target.GetSize();
            float angle = target.GetOrientation();

            Assert.AreEqual(410.207916f, size.Height);
            Assert.AreEqual(329.534637f, size.Width);
            Assert.AreEqual(0.0196384024f, angle, 1e-4);
        }

        [Test]
        public void ComputeTest2()
        {
            // 0 and 1 are only translated
            var cm0 = new CentralMoments(Image.Clone(Resources.hu0), 3);
            var cm1 = new CentralMoments(Image.Clone(Resources.hu1), 3);

            Assert.AreEqual(cm0.Mu00, cm1.Mu00);
            Assert.AreEqual(cm0.Mu01, cm1.Mu01);
            Assert.AreEqual(cm0.Mu10, cm1.Mu10);
            Assert.AreNotEqual(cm0.Mu11, cm1.Mu11);
        }
    }
}
