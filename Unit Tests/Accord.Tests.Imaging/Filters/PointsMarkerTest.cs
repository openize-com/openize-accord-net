// Accord Unit Tests
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

namespace FileFormat.Accord.Tests.Imaging.Filters
{
    using System.Collections.Generic;
    using System.Drawing;
    using Accord.Imaging;
    using FileFormat.Accord.Imaging.Filters;
    using FileFormat.Accord.Imaging.Interest_Points.FREAK;
    using FileFormat.Accord.Math.Matrix;
    using NUnit.Framework;

    [TestFixture]
    public class PointsMarkerTest
    {

        [Test]
        public void PointsMarkerTest1()
        {
            IEnumerable<FastRetinaKeypoint> points = new FastRetinaKeypoint[] 
            {
                new FastRetinaKeypoint(1, 2),
                new FastRetinaKeypoint(3, 4),
            };

            var marker = new PointsMarker(points);

            double[,] m = Matrix.Zeros(5, 5);
            Bitmap bmp = m.ToBitmap();

            marker.ApplyInPlace(bmp);

            double[,] actual = bmp.ToMatrix(0);

            double[,] expected =
            {
                { 0, 0, 0, 0, 0 },
                { 1, 1, 1, 0, 0 },
                { 1, 1, 1, 0, 0 },
                { 1, 1, 1, 1, 1 },
                { 0, 0, 1, 1, 1 },
            };

            Assert.AreEqual(expected, actual);
        }

    }
}
