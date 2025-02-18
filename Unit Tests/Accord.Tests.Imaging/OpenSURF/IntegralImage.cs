namespace Openize.Accord.Tests.Imaging.OpenSURF
{
    using System;
    using System.Drawing;

    public class IntegralImage
    {
        const float cR = .2989f;
        const float cG = .5870f;
        const float cB = .1140f;

        internal float[,] Matrix;
        public int Width, Height;

        public float this[int y, int x]
        {
            get { return this.Matrix[y, x]; }
            set { this.Matrix[y, x] = value; }
        }

        private IntegralImage(int width, int height)
        {
            this.Width = width;
            this.Height = height;

            this.Matrix = new float[height, width];
        }

        public static IntegralImage FromImage(Bitmap image)
        {
            IntegralImage pic = new IntegralImage(image.Width, image.Height);

            float rowsum = 0;
            for (int x = 0; x < image.Width; x++)
            {
                Color c = image.GetPixel(x, 0);
                rowsum += (cR * c.R + cG * c.G + cB * c.B) / 255f;
                pic[0, x] = rowsum;
            }


            for (int y = 1; y < image.Height; y++)
            {
                rowsum = 0;
                for (int x = 0; x < image.Width; x++)
                {
                    Color c = image.GetPixel(x, y);
                    rowsum += (cR * c.R + cG * c.G + cB * c.B) / 255f;

                    // integral image is rowsum + value above        
                    pic[y, x] = rowsum + pic[y - 1, x];
                }
            }

            return pic;
        }


        public float BoxIntegral(int row, int col, int rows, int cols)
        {
            // The subtraction by one for row/col is because row/col is inclusive.
            int r1 = Math.Min(row, this.Height) - 1;
            int c1 = Math.Min(col, this.Width) - 1;
            int r2 = Math.Min(row + rows, this.Height) - 1;
            int c2 = Math.Min(col + cols, this.Width) - 1;

            float A = 0, B = 0, C = 0, D = 0;
            if (r1 >= 0 && c1 >= 0) A = this.Matrix[r1, c1];
            if (r1 >= 0 && c2 >= 0) B = this.Matrix[r1, c2];
            if (r2 >= 0 && c1 >= 0) C = this.Matrix[r2, c1];
            if (r2 >= 0 && c2 >= 0) D = this.Matrix[r2, c2];

            return Math.Max(0, A - B - C + D);
        }

        /// <summary>
        /// Get Haar Wavelet X response
        /// </summary>
        /// <param name="row"></param>
        /// <param name="column"></param>
        /// <param name="size"></param>
        /// <returns></returns>
        public float HaarX(int row, int column, int size)
        {
            return this.BoxIntegral(row - size / 2, column, size, size / 2)
              - 1 * this.BoxIntegral(row - size / 2, column - size / 2, size, size / 2);
        }

        /// <summary>
        /// Get Haar Wavelet Y response
        /// </summary>
        /// <param name="row"></param>
        /// <param name="column"></param>
        /// <param name="size"></param>
        /// <returns></returns>
        public float HaarY(int row, int column, int size)
        {
            return this.BoxIntegral(row, column - size / 2, size / 2, size)
              - 1 * this.BoxIntegral(row - size / 2, column - size / 2, size / 2, size);
        }
    }
}
