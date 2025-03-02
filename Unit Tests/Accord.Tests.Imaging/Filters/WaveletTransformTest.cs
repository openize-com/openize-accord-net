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

namespace Openize.Accord.Tests.Imaging.Filters
{
    using System.Drawing;
    using Openize.Accord.Math.Matrix;
    using NUnit.Framework;
    using Image = Openize.Accord.Imaging.AForge.Imaging.Image;
    using Openize.Accord.Tests.Imaging.Properties;
    using Openize.Accord.Imaging.AForge.Imaging.Filters.Transform;
    using Openize.Accord.Imaging.Converters;
    using Openize.Accord.Imaging.Filters;
    using Openize.Accord.Math.Wavelets;
    using Openize.Accord.Math.Wavelets.Base;

    [TestFixture]
    public class WaveletTransformTest
    {

        [Test]
        public void Example1()
        {
            Bitmap image = Image.Clone(Resources.lena512);

            // Create a new Haar Wavelet transform filter
            var wavelet = new WaveletTransform(new Haar(1));

            // Apply the Wavelet transformation
            Bitmap result = wavelet.Apply(image);

            // Show on the screen
            //ImageBox.Show(result);
            Assert.IsNotNull(result);

            // Extract only one of the resulting images
            var crop = new Crop(new Rectangle(0, 0,
                image.Width / 2, image.Height / 2));

            Bitmap quarter = crop.Apply(result);

            // Show on the screen
            //ImageBox.Show(quarter);
            Assert.IsNotNull(quarter);
        }

        [Test]
        public void WaveletTransformConstructorTest()
        {
            // Start with a grayscale image
            Bitmap src = Image.Clone(Resources.lena512);

            // Create a wavelet filter            
            IWavelet wavelet = new Haar(2);
            WaveletTransform target = new WaveletTransform(wavelet);

            // Apply the transformation
            Bitmap dst = target.Apply(src);

            // Revert the transformation
            target.Backward = true;
            Bitmap org = target.Apply(dst);

            double[,] actual; new ImageToMatrix().Convert(org, out actual);
            double[,] expected; new ImageToMatrix().Convert(src, out expected);

            Assert.IsTrue(actual.IsEqual(expected, atol: 0.102));
        }
    }
}
