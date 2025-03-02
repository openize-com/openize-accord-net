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
    using Openize.Accord.Imaging.Converters;
    using Openize.Accord.Imaging.Filters;

    [TestFixture]
    public class LogExpTest
    {

        [Test]
        public void ApplyTest()
        {
            Bitmap input = Image.Clone(Resources.lena_color);

            Logarithm log = new Logarithm();
            Bitmap output = log.Apply(input);

            Exponential exp = new Exponential();
            Bitmap reconstruction = exp.Apply(output);

            //ImageBox.Show("input", input);
            //ImageBox.Show("output", output);
            //ImageBox.Show("reconstruction", reconstruction);
        }

        [Test]
        public void ProcessImageTest()
        {
            double[,] diag = Matrix.Magic(5);

            Bitmap input;
            new MatrixToImage(0, 100).Convert(diag, out input);

            Logarithm log = new Logarithm();
            Bitmap output = log.Apply(input);

            Exponential exp = new Exponential();
            Bitmap reconstruction = exp.Apply(output);


            double[,] actual;
            new ImageToMatrix(0, 100).Convert(reconstruction, out actual);
            
            double[,] expected = diag;

            Assert.IsTrue(expected.IsEqual(actual, 1.5));
        }

    }
}
