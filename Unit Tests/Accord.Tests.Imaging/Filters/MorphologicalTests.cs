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
    using NUnit.Framework;
    using Openize.Accord.Imaging.Accord.DataSets;
    using Openize.Accord.Imaging.AForge.Imaging.Filters.Morphology.Specific_Optimizations;
    using Openize.Accord.Imaging.Filters;

    [TestFixture]
    [Ignore("didn't work initially")]
    public class MorphologicalTest
    {

        [Test]
        [Category("Slow")]
        public void BinaryDilation3x3Test1()
        {
            string basePath = NUnit.Framework.TestContext.CurrentContext.TestDirectory;

            #region doc_binary_dilation_3x3
            // Let's start with one of the default 
            // color test images in the framework:
            var test = new TestImages(basePath);

            // Let's get Lena's picture
            Bitmap bmp = test["lena.bmp"];

            // And transform it to a binary mask 
            // using Niblack's threshold class
            var niblack = new NiblackThreshold();
            Bitmap binary = niblack.Apply(bmp);

            // The result can be seen below:
            // ImageBox.Show(binary);

            // Now, let's finally apply the dilation 
            // filter to the binarized image below:
            var dil3x3 = new BinaryDilation3x3();
            Bitmap result = dil3x3.Apply(binary);

            // The result can be seen below:
            // ImageBox.Show(result);
            #endregion

            //result.Save(@"C:\Projects\morpho-dilation3x3-result.png");
            //binary.Save(@"C:\Projects\morpho-dilation3x3-binary.png");
        }
    }
}
