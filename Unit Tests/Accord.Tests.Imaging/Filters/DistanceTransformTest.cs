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
    using System.Drawing;
    using System.IO;
    using FileFormat.Accord.Core.AForge.Core;
    using FileFormat.Accord.Imaging.AForge.Imaging.Filters.Binarization;
    using FileFormat.Accord.Imaging.AForge.Imaging.Filters.Color_Filters;
    using FileFormat.Accord.Imaging.Converters;
    using FileFormat.Accord.Imaging.Filters;
    using FileFormat.Accord.Math.Matrix;
    using NUnit.Framework;
    using Image = Accord.Imaging.AForge.Imaging.Image;
    using FileFormat.Accord.Tests.Imaging.Properties;

    [TestFixture]
    public class DistanceTransformTest
    {

        [Test]
        public void ApplyTest1()
        {
            string localPath = Path.Combine(NUnit.Framework.TestContext.CurrentContext.TestDirectory, "distance");
            Directory.CreateDirectory(localPath);

            // Bitmap cards = new TestImages().GetImage("cards.jpg");
            Bitmap cards = Image.Clone(Resources.cards);
            cards.Save(Path.Combine(localPath, "input.jpg"));

            Bitmap grayCards = Grayscale.CommonAlgorithms.BT709.Apply(cards);
            grayCards.Save(Path.Combine(localPath, "grayscale.jpg"));

            Bitmap binaryCards = new Threshold(200).Apply(grayCards);
            binaryCards.Save(Path.Combine(localPath, "binary.jpg"));

            var dt = new DistanceTransform();
            Bitmap result = dt.Apply(binaryCards);
            Assert.AreEqual(35.805027, dt.MaximumDistance, 1e-5);
            Assert.AreEqual(new IntPoint(254, 129), dt.UltimateErodedPoint);
            result.Save(Path.Combine(localPath, "distance_transform.jpg"));

            Assert.IsNotNull(result);
        }

        [Test]
        public void water_test()
        {
            Bitmap water = Image.Clone(Resources.water);

            var dt = new DistanceTransform();
            Bitmap result = dt.Apply(water);
            Assert.AreEqual(111.364265, dt.MaximumDistance, 1e-5);
            Assert.AreEqual(new IntPoint(436, 310), dt.UltimateErodedPoint);
            Assert.IsNotNull(result);

            float[] distance1D;
            new ImageToArray().Convert(result, out distance1D);

            float[][] distance;
            new ImageToMatrix().Convert(result, out distance);

            Assert.AreEqual(distance1D.Sum(), distance.Sum());

            float[][] distances2 = Jagged.Reshape(dt.Pixels, result.Height, result.Width);
            float actual = distances2[84][533];
            float expected = 1.4142135f;
            Assert.AreEqual(expected, actual);
        }

    }
}
