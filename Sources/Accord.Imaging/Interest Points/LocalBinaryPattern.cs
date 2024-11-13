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

namespace FileFormat.Accord.Imaging.Interest_Points
{
    using System;
    using System.Collections.Generic;
    using System.Drawing;
    using System.Drawing.Imaging;
    using AForge.Imaging;
    using AForge.Imaging.Filters.Color_Filters;
    using Base;
    using global::Accord.Math;
    using Math;
    using Math.Core;
    using Math.Matrix;

    /// <summary>
    ///   Local Binary Patterns.
    /// </summary>
    /// 
    /// <remarks>
    /// <para>
    ///    Local binary patterns (LBP) is a type of feature used for classification
    ///    in computer vision. LBP is the particular case of the Texture Spectrum 
    ///    model proposed in 1990. LBP was first described in 1994. It has since 
    ///    been found to be a powerful feature for texture classification; it has
    ///    further been determined that when LBP is combined with the Histogram of
    ///    oriented gradients (HOG) classifier, it improves the detection performance
    ///    considerably on some datasets. </para>
    /// 
    /// <para>
    ///   References:
    ///   <list type="bullet">
    ///     <item><description>
    ///       Wikipedia Contributors, "Local Binary Patterns". Available at
    ///       http://en.wikipedia.org/wiki/Local_binary_patterns </description></item>
    ///   </list>
    /// </para>
    /// </remarks>
    /// 
    /// <example>
    /// <para>
    ///   The first example shows how to extract LBP descriptors given an image.</para>
    ///   <code source="Unit Tests\Accord.Tests.Imaging\LocalBinaryPatternsTest.cs" region="doc_apply" />
    ///   <para><b>Input image:</b></para>
    ///   <img src="..\images\imaging\wood_texture.jpg" width="320" height="240" />
    ///   
    /// <para>
    ///   The second example shows how to use the LBP feature extractor as part of a
    ///   Bag-of-Words model in order to perform texture image classification:</para>
    ///   <code source="Unit Tests\Accord.Tests.Vision\Imaging\BagOfVisualWordsTest.cs" region="doc_feature_lbp" />
    /// </example>
    /// 
    [Serializable]
    public class LocalBinaryPattern : BaseFeatureExtractor<FeatureDescriptor>
    {
        const int numberOfBins = 256;

        int cellSize = 6;  // size of the cell, in number of pixels
        int blockSize = 3; // size of the block, in number of cells
        bool normalize = true;

        double epsilon = 1e-10;


        int[,] patterns;
        int[,][] histograms;


        /// <summary>
        ///   Gets the size of a cell, in pixels. Default is 6.
        /// </summary>
        /// 
        public int CellSize
        {
            get { return this.cellSize; }
            set { this.cellSize = value; }
        }

        /// <summary>
        ///   Gets the size of a block, in pixels. Default is 3.
        /// </summary>
        /// 
        public int BlockSize
        {
            get { return this.blockSize; }
            set { this.blockSize = value; }
        }

        /// <summary>
        ///   Gets the set of local binary patterns computed for each
        ///   pixel in the last call to to <see cref="BaseFeatureExtractor{TFeature}.Transform(Bitmap)"/>.
        /// </summary>
        /// 
        public int[,] Patterns { get { return this.patterns; } }

        /// <summary>
        ///   Gets the histogram computed at each cell.
        /// </summary>
        /// 
        public int[,][] Histograms { get { return this.histograms; } }

        /// <summary>
        ///   Gets or sets whether to normalize final 
        ///   histogram feature vectors. Default is true.
        /// </summary>
        /// 
        public bool Normalize
        {
            get { return this.normalize; }
            set { this.normalize = value; }
        }


        /// <summary>
        ///   Initializes a new instance of the <see cref="LocalBinaryPattern"/> class.
        /// </summary>
        /// 
        /// <param name="blockSize">
        ///   The size of a block, measured in cells. Default is 3.</param>
        /// <param name="cellSize">
        ///   The size of a cell, measured in pixels. If set to zero, the entire
        ///   image will be used at once, forming a single block. Default is 6.</param>
        /// <param name="normalize">
        ///   Whether to normalize generated histograms. Default is true.</param>
        /// 
        public LocalBinaryPattern(int blockSize = 3, int cellSize = 6, bool normalize = true)
        {
            this.cellSize = cellSize;
            this.blockSize = blockSize;
            this.normalize = normalize;

            base.SupportedFormats.UnionWith(new[] {
                PixelFormat.Format8bppIndexed,
                PixelFormat.Format24bppRgb,
                PixelFormat.Format32bppRgb,
                PixelFormat.Format32bppArgb });
        }

        /// <summary>
        ///   This method should be implemented by inheriting classes to implement the 
        ///   actual feature extraction, transforming the input image into a list of features.
        /// </summary>
        /// 
        protected override IEnumerable<FeatureDescriptor> InnerTransform(UnmanagedImage image)
        {
            // make sure we have grayscale image
            UnmanagedImage grayImage = null;

            if (image.PixelFormat == PixelFormat.Format8bppIndexed)
            {
                grayImage = image;
            }
            else
            {
                // create temporary grayscale image
                grayImage = Grayscale.CommonAlgorithms.BT709.Apply(image);
            }


            // get source image size
            int width = grayImage.Width;
            int height = grayImage.Height;
            int stride = grayImage.Stride;
            int offset = stride - width;

            // 1. Calculate 8-pixel neighborhood binary patterns 
            if (this.patterns == null || height > this.patterns.GetLength(0) || width > this.patterns.GetLength(1))
            {
                this.patterns = new int[height, width];
            }
            else
            {
                System.Diagnostics.Debug.Write(String.Format("Reusing storage for patterns. " +
                    "Need ({0}, {1}), have ({1}, {2})", height, width, this.patterns.Rows(), this.patterns.Columns()));
            }

            unsafe
            {
                fixed (int* ptrPatterns = this.patterns)
                {
                    // Begin skipping first line
                    byte* src = (byte*)grayImage.ImageData.ToPointer() + stride;
                    int* neighbors = ptrPatterns + width;

                    // for each line
                    for (int y = 1; y < height - 1; y++)
                    {
                        // skip first column
                        neighbors++; src++;

                        // for each inner pixel in line (skipping first and last)
                        for (int x = 1; x < width - 1; x++, src++, neighbors++)
                        {
                            // Retrieve the pixel neighborhood
                            byte a11 = src[+stride + 1], a12 = src[+1], a13 = src[-stride + 1];
                            byte a21 = src[+stride + 0], a22 = src[0], a23 = src[-stride + 0];
                            byte a31 = src[+stride - 1], a32 = src[-1], a33 = src[-stride - 1];

                            int sum = 0;
                            if (a22 < a11) sum += 1 << 0;
                            if (a22 < a12) sum += 1 << 1;
                            if (a22 < a13) sum += 1 << 2;
                            if (a22 < a21) sum += 1 << 3;
                            if (a22 < a23) sum += 1 << 4;
                            if (a22 < a31) sum += 1 << 5;
                            if (a22 < a32) sum += 1 << 6;
                            if (a22 < a33) sum += 1 << 7;

                            *neighbors = sum;
                        }

                        // Skip last column
                        neighbors++; src += offset + 1;
                    }
                }
            }

            // Free some resources which wont be needed anymore
            if (image.PixelFormat != PixelFormat.Format8bppIndexed)
                grayImage.Dispose();


            // 2. Compute cell histograms
            int cellCountX;
            int cellCountY;

            if (this.cellSize > 0)
            {
                cellCountX = (int)Math.Floor(width / (double)this.cellSize);
                cellCountY = (int)Math.Floor(height / (double)this.cellSize);

                if (this.histograms == null || cellCountX > this.histograms.Rows() || cellCountY > this.histograms.Columns())
                {
                    this.histograms = new int[cellCountX, cellCountY][];
                    for (int i = 0; i < cellCountX; i++)
                        for (int j = 0; j < cellCountY; j++)
                            this.histograms[i, j] = new int[numberOfBins];
                }
                else
                {
                    System.Diagnostics.Debug.Write(String.Format("Reusing storage for histograms. " +
                        "Need ({0}, {1}), have ({1}, {2})", cellCountX, cellCountY, this.histograms.Rows(), this.histograms.Columns()));
                }

                // For each cell
                for (int i = 0; i < cellCountX; i++)
                {
                    for (int j = 0; j < cellCountY; j++)
                    {
                        // Compute the histogram
                        int[] histogram = this.histograms[i, j];

                        int startCellX = i * this.cellSize;
                        int startCellY = j * this.cellSize;

                        // for each pixel in the cell
                        for (int x = 0; x < this.cellSize; x++)
                            for (int y = 0; y < this.cellSize; y++)
                                histogram[this.patterns[startCellY + y, startCellX + x]]++;
                    }
                }
            }
            else
            {
                cellCountX = 1;
                cellCountY = 1;

                if (this.histograms == null)
                {
                    this.histograms = new int[,][] { { new int[numberOfBins] } };
                }
                else
                {
                    System.Diagnostics.Debug.Write(String.Format("Reusing storage for histograms. " +
                        "Need ({0}, {1}), have ({1}, {2})", cellCountX, cellCountY, this.histograms.Rows(), this.histograms.Columns()));
                }

                int[] histogram = this.histograms[0, 0];

                for (int i = 0; i < height; i++)
                    for (int j = 0; j < width; j++)
                        histogram[this.patterns[i, j]]++;
            }

            // 3. Group the cells into larger, normalized blocks
            int blocksCountX;
            int blocksCountY;

            if (this.blockSize > 0)
            {
                blocksCountX = (int)Math.Floor(cellCountX / (double)this.blockSize);
                blocksCountY = (int)Math.Floor(cellCountY / (double)this.blockSize);
            }
            else
            {
                this.blockSize = blocksCountX = blocksCountY = 1;
            }


            var blocks = new List<FeatureDescriptor>();

            for (int i = 0; i < blocksCountX; i++)
            {
                for (int j = 0; j < blocksCountY; j++)
                {
                    double[] block = new double[this.blockSize * this.blockSize * numberOfBins];

                    int startBlockX = i * this.blockSize;
                    int startBlockY = j * this.blockSize;
                    int c = 0;

                    // for each cell in the block
                    for (int x = 0; x < this.blockSize; x++)
                    {
                        for (int y = 0; y < this.blockSize; y++)
                        {
                            int[] histogram = this.histograms[startBlockX + x, startBlockY + y];

                            // Copy all histograms to the block vector
                            for (int k = 0; k < histogram.Length; k++)
                                block[c++] = histogram[k];
                        }
                    }

                    // TODO: Remove this block and instead propose a general architecture 
                    //       for applying normalizations to descriptor blocks
                    if (this.normalize)
                        block.Divide(block.Euclidean() + this.epsilon, result: block);

                    blocks.Add(block);
                }
            }

            return blocks;
        }

        /// <summary>
        /// Creates a new object that is a copy of the current instance.
        /// </summary>
        /// 
        protected override object Clone(ISet<PixelFormat> supportedFormats)
        {
            var clone = new LocalBinaryPattern(this.blockSize, this.cellSize, this.normalize);
            clone.epsilon = this.epsilon;
            clone.SupportedFormats = supportedFormats;
            return clone;
        }

    }
}
