// Accord Math Library
// The Accord.NET Framework
// http://accord-framework.net
//
// Copyright © César Souza, 2009-2017
// cesarsouza at gmail.com
//
// Original work copyright © Lutz Roeder, 2000
//  Adapted from Mapack for .NET, September 2000
//  Adapted from Mapack for COM and Jama routines
//  http://www.aisto.com/roeder/dotnet
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

namespace Openize.Accord.Math.Decompositions
{
    using System;
    using Openize.Accord.Math.Matrix;
    using Openize.Accord.Math;

    /// <summary>
    ///     Determines the eigenvalues and eigenvectors of a real square matrix.
    /// </summary>
    ///
    /// <remarks>
    ///   <para>
    ///     In the mathematical discipline of linear algebra, eigendecomposition
    ///     or sometimes spectral decomposition is the factorization of a matrix
    ///     into a canonical form, whereby the matrix is represented in terms of
    ///     its eigenvalues and eigenvectors.</para>
    ///   <para>
    ///     If <c>A</c> is symmetric, then <c>A = V * D * V'</c> and <c>A = V * V'</c>
    ///     where the eigenvalue matrix <c>D</c> is diagonal and the eigenvector matrix <c>V</c> is orthogonal.
    ///     If <c>A</c> is not symmetric, the eigenvalue matrix <c>D</c> is block diagonal
    ///     with the real eigenvalues in 1-by-1 blocks and any complex eigenvalues,
    ///     <c>lambda + i*mu</c>, in 2-by-2 blocks, <c>[lambda, mu; -mu, lambda]</c>.
    ///     The columns of <c>V</c> represent the eigenvectors in the sense that <c>A * V = V * D</c>.
    ///     The matrix V may be badly conditioned, or even singular, so the validity of the equation
    ///     <c>A = V * D * inverse(V)</c> depends upon the condition of <c>V</c>.
    ///   </para>
    /// </remarks>
    /// 
    public sealed class JaggedEigenvalueDecomposition : ICloneable
    {
        private int n;              // matrix dimension
        private Double[] d, e;      // storage of eigenvalues.
        private Double[][] V;       // storage of eigenvectors.
        private Double[][] H;       // storage of nonsymmetric Hessenberg form.
        private Double[] ort;       // storage for nonsymmetric algorithm.
        private bool symmetric;

		private int? rank;
		private Double[][] diagonalMatrix;

        private const Double eps = 2 * Constants.DoubleEpsilon;


        /// <summary>
        ///   Returns the effective numerical matrix rank.
        /// </summary>
        ///
        /// <value>Number of non-negligible eigen values.</value>
        ///
        public int Rank
        {
            get
            {
				if (this.rank.HasValue)
					return this.rank.Value;

                Double tol = this.n * this.d[0] * eps;

                int r = 0;
                for (int i = 0; i < this.d.Length; i++)
                    if (this.d[i] > tol) r++;

                return (int)(this.rank = r);
            }
        }

        /// <summary>
        ///   Construct an eigenvalue decomposition.</summary>
        /// <param name="value">
        ///   The matrix to be decomposed.</param>
        /// <param name="inPlace">
        ///   Pass <see langword="true"/> to perform the decomposition in place. The matrix
        ///   <paramref name="value"/> will be destroyed in the process, resulting in less
        ///   memory comsumption.</param>
        /// <param name="sort">
        ///   Pass <see langword="true"/> to sort the eigenvalues and eigenvectors at the end
        ///   of the decomposition.</param>
        ///
        public JaggedEigenvalueDecomposition(Double[][] value, bool inPlace = false, bool sort = false)
            : this(value, assumeSymmetric: value.IsSymmetric(), inPlace: inPlace, sort: sort)
        {
        }
        
        /// <summary>
        ///   Construct an eigenvalue decomposition.</summary>
        /// <param name="value">
        ///   The matrix to be decomposed.</param>
        /// <param name="assumeSymmetric">
        ///   Defines if the matrix should be assumed as being symmetric
        ///   regardless if it is or not. Default is <see langword="false"/>.</param>
        /// <param name="inPlace">
        ///   Pass <see langword="true"/> to perform the decomposition in place. The matrix
        ///   <paramref name="value"/> will be destroyed in the process, resulting in less
        ///   memory comsumption.</param>
        /// <param name="sort">
        ///   Pass <see langword="true"/> to sort the eigenvalues and eigenvectors at the end
        ///   of the decomposition.</param>
        ///
        public JaggedEigenvalueDecomposition(Double[][] value, bool assumeSymmetric, 
            bool inPlace = false, bool sort = false)
        {
            if (value == null)
                throw new ArgumentNullException("value", "Matrix cannot be null.");

            if (value.Length != value[0].Length)
                throw new ArgumentException("Matrix is not a square matrix.", "value");

            this.n = value.Length;
            this.V = new Double[this.n][];
            for (int i = 0; i < this.n; i++)
              this.V[i] = new Double[this.n];
            this.d = new Double[this.n];
            this.e = new Double[this.n];


            this.symmetric = assumeSymmetric;

            if (this.symmetric)
            {
                this.V = inPlace ? value : (Double[][])value.MemberwiseClone();

                // Tridiagonalize.
                this.tred2();

                // Diagonalize.
                this.tql2();
            }
            else
            {
                this.H = inPlace ? value : (Double[][])value.MemberwiseClone();

                this.ort = new Double[this.n];

                // Reduce to Hessenberg form.
                this.orthes();

                // Reduce Hessenberg to real Schur form.
                this.hqr2();
            }

            if (sort)
            {
                // Sort eigenvalues and vectors in descending order
                var idx = Openize.Accord.Math.Vector.Vector.Range(this.n);
                Array.Sort(idx, (i, j) => 
                {
                    if (Math.Abs(this.d[i]) == Math.Abs(this.d[j]))
                        return -Math.Abs(this.e[i]).CompareTo(Math.Abs(this.e[j]));
                    return -Math.Abs(this.d[i]).CompareTo(Math.Abs(this.d[j]));
                });

                this.d = this.d.Get(idx);
                this.e = this.e.Get(idx);
                this.V = this.V.Get(null, idx);
            }
        }


        /// <summary>Returns the real parts of the eigenvalues.</summary>
        public Double[] RealEigenvalues
        {
            get { return this.d; }
        }

        /// <summary>Returns the imaginary parts of the eigenvalues.</summary>    
        public Double[] ImaginaryEigenvalues
        {
            get { return this.e; }
        }

        /// <summary>Returns the eigenvector matrix.</summary>
        public Double[][] Eigenvectors
        {
            get { return this.V; }
        }

        /// <summary>Returns the block diagonal eigenvalue matrix.</summary>
        public Double[][] DiagonalMatrix
        {
            get
            {
				if (this.diagonalMatrix != null)
					return this.diagonalMatrix;

                var x = new Double[this.n][];
                for (int i = 0; i < this.n; i++)
                  x[i] = new Double[this.n];

                for (int i = 0; i < this.n; i++)
                {
                    for (int j = 0; j < this.n; j++)
                        x[i][j] = 0;

                    x[i][i] = this.d[i];
                    if (this.e[i] > 0)
                    {
                        x[i][i + 1] = this.e[i];
                    }
                    else if (this.e[i] < 0)
                    {
                        x[i][i - 1] = this.e[i];
                    }
                }

                return this.diagonalMatrix = x;
            }
        }


        #region Private methods
        private void tred2()
        {
            // Symmetric Householder reduction to tridiagonal form.
            // This is derived from the Algol procedures tred2 by Bowdler, Martin, Reinsch, and Wilkinson, 
            // Handbook for Auto. Comp., Vol.ii-Linear Algebra, and the corresponding Fortran subroutine in EISPACK.
            for (int j = 0; j < this.n; j++)
                this.d[j] = this.V[this.n - 1][j];

            // Householder reduction to tridiagonal form.
            for (int i = this.n - 1; i > 0; i--)
            {
                // Scale to avoid under/overflow.
                Double scale = 0;
                Double h = 0;
                for (int k = 0; k < i; k++)
                    scale = scale + System.Math.Abs(this.d[k]);

                if (scale == 0)
                {
                    this.e[i] = this.d[i - 1];
                    for (int j = 0; j < i; j++)
                    {
                        this.d[j] = this.V[i - 1][j];
                        this.V[i][j] = 0;
                        this.V[j][i] = 0;
                    }
                }
                else
                {
                    // Generate Householder vector.
                    for (int k = 0; k < i; k++)
                    {
                        this.d[k] /= scale;
                        h += this.d[k] * this.d[k];
                    }

                    Double f = this.d[i - 1];
                    Double g = (Double)System.Math.Sqrt(h);
                    if (f > 0) g = -g;

                    this.e[i] = scale * g;
                    h = h - f * g;
                    this.d[i - 1] = f - g;
                    for (int j = 0; j < i; j++)
                        this.e[j] = 0;

                    // Apply similarity transformation to remaining columns.
                    for (int j = 0; j < i; j++)
                    {
                        f = this.d[j];
                        this.V[j][i] = f;
                        g = this.e[j] + this.V[j][j] * f;
                        for (int k = j + 1; k <= i - 1; k++)
                        {
                            g += this.V[k][j] * this.d[k];
                            this.e[k] += this.V[k][j] * f;
                        }
                        this.e[j] = g;
                    }

                    f = 0;
                    for (int j = 0; j < i; j++)
                    {
                        this.e[j] /= h;
                        f += this.e[j] * this.d[j];
                    }

                    Double hh = f / (h + h);
                    for (int j = 0; j < i; j++)
                        this.e[j] -= hh * this.d[j];

                    for (int j = 0; j < i; j++)
                    {
                        f = this.d[j];
                        g = this.e[j];
                        for (int k = j; k <= i - 1; k++)
                            this.V[k][j] -= (f * this.e[k] + g * this.d[k]);

                        this.d[j] = this.V[i - 1][j];
                        this.V[i][j] = 0;
                    }
                }
                this.d[i] = h;
            }

            // Accumulate transformations.
            for (int i = 0; i < this.n - 1; i++)
            {
                this.V[this.n - 1][i] = this.V[i][i];
                this.V[i][i] = 1;
                Double h = this.d[i + 1];
                if (h != 0)
                {
                    for (int k = 0; k <= i; k++)
                        this.d[k] = this.V[k][i + 1] / h;

                    for (int j = 0; j <= i; j++)
                    {
                        Double g = 0;
                        for (int k = 0; k <= i; k++)
                            g += this.V[k][i + 1] * this.V[k][j];
                        for (int k = 0; k <= i; k++)
                            this.V[k][j] -= g * this.d[k];
                    }
                }

                for (int k = 0; k <= i; k++)
                    this.V[k][i + 1] = 0;
            }

            for (int j = 0; j < this.n; j++)
            {
                this.d[j] = this.V[this.n - 1][j];
                this.V[this.n - 1][j] = 0;
            }

            this.V[this.n - 1][this.n - 1] = 1;
            this.e[0] = 0;
        }

        private void tql2()
        {
            // Symmetric tridiagonal QL algorithm.
            // This is derived from the Algol procedures tql2, by Bowdler, Martin, Reinsch, and Wilkinson, 
            // Handbook for Auto. Comp., Vol.ii-Linear Algebra, and the corresponding Fortran subroutine in EISPACK.
            for (int i = 1; i < this.n; i++)
                this.e[i - 1] = this.e[i];

            this.e[this.n - 1] = 0;

            Double f = 0;
            Double tst1 = 0;
            Double eps = 2 * Constants.DoubleEpsilon;

            for (int l = 0; l < this.n; l++)
            {
                // Find small subdiagonal element.
                tst1 = System.Math.Max(tst1, System.Math.Abs(this.d[l]) + System.Math.Abs(this.e[l]));
                int m = l;
                while (m < this.n)
                {
                    if (System.Math.Abs(this.e[m]) <= eps * tst1)
                        break;
                    m++;
                }

                // If m == l, d[l] is an eigenvalue, otherwise, iterate.
                if (m > l)
                {
                    int iter = 0;
                    do
                    {
                        iter = iter + 1;  // (Could check iteration count here.)

                        // Compute implicit shift
                        Double g = this.d[l];
                        Double p = (this.d[l + 1] - g) / (2 * this.e[l]);
                        Double r = global::Openize.Accord.Math.Tools.Hypotenuse(p, 1);
                        if (p < 0)
                        {
                            r = -r;
                        }

                        this.d[l] = this.e[l] / (p + r);
                        this.d[l + 1] = this.e[l] * (p + r);
                        Double dl1 = this.d[l + 1];
                        Double h = g - this.d[l];
                        for (int i = l + 2; i < this.n; i++)
                        {
                            this.d[i] -= h;
                        }

                        f = f + h;

                        // Implicit QL transformation.
                        p = this.d[m];
                        Double c = 1;
                        Double c2 = c;
                        Double c3 = c;
                        Double el1 = this.e[l + 1];
                        Double s = 0;
                        Double s2 = 0;
                        for (int i = m - 1; i >= l; i--)
                        {
                            c3 = c2;
                            c2 = c;
                            s2 = s;
                            g = c * this.e[i];
                            h = c * p;
                            r = global::Openize.Accord.Math.Tools.Hypotenuse(p, this.e[i]);
                            this.e[i + 1] = s * r;
                            s = this.e[i] / r;
                            c = p / r;
                            p = c * this.d[i] - s * g;
                            this.d[i + 1] = h + s * (c * g + s * this.d[i]);

                            // Accumulate transformation.
                            for (int k = 0; k < this.n; k++)
                            {
                                h = this.V[k][i + 1];
                                this.V[k][i + 1] = s * this.V[k][i] + c * h;
                                this.V[k][i] = c * this.V[k][i] - s * h;
                            }
                        }

                        p = -s * s2 * c3 * el1 * this.e[l] / dl1;
                        this.e[l] = s * p;
                        this.d[l] = c * p;

                        // Check for convergence.
                    }
                    while (System.Math.Abs(this.e[l]) > eps * tst1);
                }
                this.d[l] = this.d[l] + f;
                this.e[l] = 0;
            }

            // Sort eigenvalues and corresponding vectors.
            for (int i = 0; i < this.n - 1; i++)
            {
                int k = i;
                Double p = this.d[i];
                for (int j = i + 1; j < this.n; j++)
                {
                    if (this.d[j] < p)
                    {
                        k = j;
                        p = this.d[j];
                    }
                }

                if (k != i)
                {
                    this.d[k] = this.d[i];
                    this.d[i] = p;
                    for (int j = 0; j < this.n; j++)
                    {
                        p = this.V[j][i];
                        this.V[j][i] = this.V[j][k];
                        this.V[j][k] = p;
                    }
                }
            }
        }

        private void orthes()
        {
            // Nonsymmetric reduction to Hessenberg form.
            // This is derived from the Algol procedures orthes and ortran, by Martin and Wilkinson, 
            // Handbook for Auto. Comp., Vol.ii-Linear Algebra, and the corresponding Fortran subroutines in EISPACK.
            int low = 0;
            int high = this.n - 1;

            for (int m = low + 1; m <= high - 1; m++)
            {
                // Scale column.

                Double scale = 0;
                for (int i = m; i <= high; i++)
                    scale = scale + System.Math.Abs(this.H[i][m - 1]);

                if (scale != 0)
                {
                    // Compute Householder transformation.
                    Double h = 0;
                    for (int i = high; i >= m; i--)
                    {
                        this.ort[i] = this.H[i][m - 1] / scale;
                        h += this.ort[i] * this.ort[i];
                    }

                    Double g = (Double)System.Math.Sqrt(h);
                    if (this.ort[m] > 0) g = -g;

                    h = h - this.ort[m] * g;
                    this.ort[m] = this.ort[m] - g;

                    // Apply Householder similarity transformation
                    // H = (I - u * u' / h) * H * (I - u * u') / h)
                    for (int j = m; j < this.n; j++)
                    {
                        Double f = 0;
                        for (int i = high; i >= m; i--)
                            f += this.ort[i] * this.H[i][j];

                        f = f / h;
                        for (int i = m; i <= high; i++)
                            this.H[i][j] -= f * this.ort[i];
                    }

                    for (int i = 0; i <= high; i++)
                    {
                        Double f = 0;
                        for (int j = high; j >= m; j--)
                            f += this.ort[j] * this.H[i][j];

                        f = f / h;
                        for (int j = m; j <= high; j++)
                            this.H[i][j] -= f * this.ort[j];
                    }

                    this.ort[m] = scale * this.ort[m];
                    this.H[m][m - 1] = scale * g;
                }
            }

            // Accumulate transformations (Algol's ortran).
            for (int i = 0; i < this.n; i++)
                for (int j = 0; j < this.n; j++)
                    this.V[i][j] = (i == j ? 1 : 0);

            for (int m = high - 1; m >= low + 1; m--)
            {
                if (this.H[m][m - 1] != 0)
                {
                    for (int i = m + 1; i <= high; i++)
                        this.ort[i] = this.H[i][m - 1];

                    for (int j = m; j <= high; j++)
                    {
                        Double g = 0;
                        for (int i = m; i <= high; i++)
                            g += this.ort[i] * this.V[i][j];

                        // Double division avoids possible underflow.
                        g = (g / this.ort[m]) / this.H[m][m - 1];
                        for (int i = m; i <= high; i++)
                            this.V[i][j] += g * this.ort[i];
                    }
                }
            }
        }

        private static void cdiv(Double xr, Double xi, Double yr, Double yi,
            out Double cdivr, out Double cdivi)
        {
            // Complex scalar division.
            Double r;
            Double d;
            if (System.Math.Abs(yr) > System.Math.Abs(yi))
            {
                r = yi / yr;
                d = yr + r * yi;
                cdivr = (xr + r * xi) / d;
                cdivi = (xi - r * xr) / d;
            }
            else
            {
                r = yr / yi;
                d = yi + r * yr;
                cdivr = (r * xr + xi) / d;
                cdivi = (r * xi - xr) / d;
            }
        }

        private void hqr2()
        {
            // Nonsymmetric reduction from Hessenberg to real Schur form.   
            // This is derived from the Algol procedure hqr2, by Martin and Wilkinson, Handbook for Auto. Comp.,
            // Vol.ii-Linear Algebra, and the corresponding  Fortran subroutine in EISPACK.
            int nn = this.n;
            int n = nn - 1;
            int low = 0;
            int high = nn - 1;
            Double eps = 2 * Constants.DoubleEpsilon;
            Double exshift = 0;
            Double p = 0;
            Double q = 0;
            Double r = 0;
            Double s = 0;
            Double z = 0;
            Double t;
            Double w;
            Double x;
            Double y;

            // Store roots isolated by balanc and compute matrix norm
            Double norm = 0;
            for (int i = 0; i < nn; i++)
            {
                if (i < low | i > high)
                {
                    this.d[i] = this.H[i][i];
                    this.e[i] = 0;
                }

                for (int j = System.Math.Max(i - 1, 0); j < nn; j++)
                    norm = norm + System.Math.Abs(this.H[i][j]);
            }

            // Outer loop over eigenvalue index
            int iter = 0;
            while (n >= low)
            {
                // Look for single small sub-diagonal element
                int l = n;
                while (l > low)
                {
                    s = System.Math.Abs(this.H[l - 1][l - 1]) + System.Math.Abs(this.H[l][l]);

                    if (s == 0) 
                        s = norm;

                    if (Double.IsNaN(s))
                        break;

                    if (System.Math.Abs(this.H[l][l - 1]) < eps * s)
                        break;

                    l--;
                }

                // Check for convergence
                if (l == n)
                {
                    // One root found
                    this.H[n][n] = this.H[n][n] + exshift;
                    this.d[n] = this.H[n][n];
                    this.e[n] = 0;
                    n--;
                    iter = 0;
                }
                else if (l == n - 1)
                {
                    // Two roots found
                    w = this.H[n][n - 1] * this.H[n - 1][n];
                    p = (this.H[n - 1][n - 1] - this.H[n][n]) / 2;
                    q = p * p + w;
                    z = (Double)System.Math.Sqrt(System.Math.Abs(q));
                    this.H[n][n] = this.H[n][n] + exshift;
                    this.H[n - 1][n - 1] = this.H[n - 1][n - 1] + exshift;
                    x = this.H[n][n];

                    if (q >= 0)
                    {
                        // Real pair
                        z = (p >= 0) ? (p + z) : (p - z);
                        this.d[n - 1] = x + z;
                        this.d[n] = this.d[n - 1];
                        if (z != 0)
                            this.d[n] = x - w / z;
                        this.e[n - 1] = 0;
                        this.e[n] = 0;
                        x = this.H[n][n - 1];
                        s = System.Math.Abs(x) + System.Math.Abs(z);
                        p = x / s;
                        q = z / s;
                        r = (Double)System.Math.Sqrt(p * p + q * q);
                        p = p / r;
                        q = q / r;

                        // Row modification
                        for (int j = n - 1; j < nn; j++)
                        {
                            z = this.H[n - 1][j];
                            this.H[n - 1][j] = q * z + p * this.H[n][j];
                            this.H[n][j] = q * this.H[n][j] - p * z;
                        }

                        // Column modification
                        for (int i = 0; i <= n; i++)
                        {
                            z = this.H[i][n - 1];
                            this.H[i][n - 1] = q * z + p * this.H[i][n];
                            this.H[i][n] = q * this.H[i][n] - p * z;
                        }

                        // Accumulate transformations
                        for (int i = low; i <= high; i++)
                        {
                            z = this.V[i][n - 1];
                            this.V[i][n - 1] = q * z + p * this.V[i][n];
                            this.V[i][n] = q * this.V[i][n] - p * z;
                        }
                    }
                    else
                    {
                        // Complex pair
                        this.d[n - 1] = x + p;
                        this.d[n] = x + p;
                        this.e[n - 1] = z;
                        this.e[n] = -z;
                    }

                    n = n - 2;
                    iter = 0;
                }
                else
                {
                    // No convergence yet     

                    // Form shift
                    x = this.H[n][n];
                    y = 0;
                    w = 0;
                    if (l < n)
                    {
                        y = this.H[n - 1][n - 1];
                        w = this.H[n][n - 1] * this.H[n - 1][n];
                    }

                    // Wilkinson's original ad hoc shift
                    if (iter == 10)
                    {
                        exshift += x;
                        for (int i = low; i <= n; i++)
                            this.H[i][i] -= x;

                        s = System.Math.Abs(this.H[n][n - 1]) + System.Math.Abs(this.H[n - 1][n - 2]);
                        x = y = (Double)0.75 * s;
                        w = (Double)(-0.4375) * s * s;
                    }

                    // MATLAB's new ad hoc shift
                    if (iter == 30)
                    {
                        s = (y - x) / 2;
                        s = s * s + w;
                        if (s > 0)
                        {
                            s = (Double)System.Math.Sqrt(s);
                            if (y < x) s = -s;
                            s = x - w / ((y - x) / 2 + s);
                            for (int i = low; i <= n; i++)
                                this.H[i][i] -= s;
                            exshift += s;
                            x = y = w = (Double)0.964;
                        }
                    }

                    iter = iter + 1;

                    // Look for two consecutive small sub-diagonal elements
                    int m = n - 2;
                    while (m >= l)
                    {
                        z = this.H[m][m];
                        r = x - z;
                        s = y - z;
                        p = (r * s - w) / this.H[m + 1][m] + this.H[m][m + 1];
                        q = this.H[m + 1][m + 1] - z - r - s;
                        r = this.H[m + 2][m + 1];
                        s = System.Math.Abs(p) + System.Math.Abs(q) + System.Math.Abs(r);
                        p = p / s;
                        q = q / s;
                        r = r / s;
                        if (m == l)
                            break;
                        if (System.Math.Abs(this.H[m][m - 1]) * (System.Math.Abs(q) + System.Math.Abs(r)) < eps * (System.Math.Abs(p) * (System.Math.Abs(this.H[m - 1][m - 1]) + System.Math.Abs(z) + System.Math.Abs(this.H[m + 1][m + 1]))))
                            break;
                        m--;
                    }

                    for (int i = m + 2; i <= n; i++)
                    {
                        this.H[i][i - 2] = 0;
                        if (i > m + 2)
                            this.H[i][i - 3] = 0;
                    }

                    // Double QR step involving rows l:n and columns m:n
                    for (int k = m; k <= n - 1; k++)
                    {
                        bool notlast = (k != n - 1);
                        if (k != m)
                        {
                            p = this.H[k][k - 1];
                            q = this.H[k + 1][k - 1];
                            r = (notlast ? this.H[k + 2][k - 1] : 0);
                            x = System.Math.Abs(p) + System.Math.Abs(q) + System.Math.Abs(r);
                            if (x != 0)
                            {
                                p = p / x;
                                q = q / x;
                                r = r / x;
                            }
                        }

                        if (x == 0) 
                            break;

                        s = (Double)System.Math.Sqrt(p * p + q * q + r * r);
                        if (p < 0) s = -s;

                        if (s != 0)
                        {
                            if (k != m)
                                this.H[k][k - 1] = -s * x;
                            else
                                if (l != m)
                                    this.H[k][k - 1] = -this.H[k][k - 1];

                            p = p + s;
                            x = p / s;
                            y = q / s;
                            z = r / s;
                            q = q / p;
                            r = r / p;

                            // Row modification
                            for (int j = k; j < nn; j++)
                            {
                                p = this.H[k][j] + q * this.H[k + 1][j];
                                if (notlast)
                                {
                                    p = p + r * this.H[k + 2][j];
                                    this.H[k + 2][j] = this.H[k + 2][j] - p * z;
                                }

                                this.H[k][j] = this.H[k][j] - p * x;
                                this.H[k + 1][j] = this.H[k + 1][j] - p * y;
                            }

                            // Column modification
                            for (int i = 0; i <= System.Math.Min(n, k + 3); i++)
                            {
                                p = x * this.H[i][k] + y * this.H[i][k + 1];
                                if (notlast)
                                {
                                    p = p + z * this.H[i][k + 2];
                                    this.H[i][k + 2] = this.H[i][k + 2] - p * r;
                                }

                                this.H[i][k] = this.H[i][k] - p;
                                this.H[i][k + 1] = this.H[i][k + 1] - p * q;
                            }

                            // Accumulate transformations
                            for (int i = low; i <= high; i++)
                            {
                                p = x * this.V[i][k] + y * this.V[i][k + 1];
                                if (notlast)
                                {
                                    p = p + z * this.V[i][k + 2];
                                    this.V[i][k + 2] = this.V[i][k + 2] - p * r;
                                }

                                this.V[i][k] = this.V[i][k] - p;
                                this.V[i][k + 1] = this.V[i][k + 1] - p * q;
                            }
                        }
                    }
                }
            }

            // Backsubstitute to find vectors of upper triangular form
            if (norm == 0)
            {
                return;
            }

            for (n = nn - 1; n >= 0; n--)
            {
                p = this.d[n];
                q = this.e[n];

                // Real vector
                if (q == 0)
                {
                    int l = n;
                    this.H[n][n] = 1;
                    for (int i = n - 1; i >= 0; i--)
                    {
                        w = this.H[i][i] - p;
                        r = 0;
                        for (int j = l; j <= n; j++)
                            r = r + this.H[i][j] * this.H[j][n];

                        if (this.e[i] < 0)
                        {
                            z = w;
                            s = r;
                        }
                        else
                        {
                            l = i;
                            if (this.e[i] == 0)
                            {
                                this.H[i][n] = (w != 0) ? (-r / w) : (-r / (eps * norm));
                            }
                            else
                            {
                                // Solve real equations
                                x = this.H[i][i + 1];
                                y = this.H[i + 1][i];
                                q = (this.d[i] - p) * (this.d[i] - p) + this.e[i] * this.e[i];
                                t = (x * s - z * r) / q;
                                this.H[i][n] = t;
                                this.H[i + 1][n] = (System.Math.Abs(x) > System.Math.Abs(z)) ? ((-r - w * t) / x) : ((-s - y * t) / z);
                            }

                            // Overflow control
                            t = System.Math.Abs(this.H[i][n]);
                            if ((eps * t) * t > 1)
                                for (int j = i; j <= n; j++)
                                    this.H[j][n] = this.H[j][n] / t;
                        }
                    }
                }
                else if (q < 0)
                {
                    // Complex vector
                    int l = n - 1;

                    // Last vector component imaginary so matrix is triangular
                    if (System.Math.Abs(this.H[n][n - 1]) > System.Math.Abs(this.H[n - 1][n]))
                    {
                        this.H[n - 1][n - 1] = q / this.H[n][n - 1];
                        this.H[n - 1][n] = -(this.H[n][n] - p) / this.H[n][n - 1];
                    }
                    else
                    {
                        cdiv(0, -this.H[n - 1][n], this.H[n - 1][n - 1] - p, q, out this.H[n - 1][n - 1], out this.H[n - 1][n]);
                    }

                    this.H[n][n - 1] = 0;
                    this.H[n][n] = 1;
                    for (int i = n - 2; i >= 0; i--)
                    {
                        Double ra, sa, vr, vi;
                        ra = 0;
                        sa = 0;
                        for (int j = l; j <= n; j++)
                        {
                            ra = ra + this.H[i][j] * this.H[j][n - 1];
                            sa = sa + this.H[i][j] * this.H[j][n];
                        }

                        w = this.H[i][i] - p;

                        if (this.e[i] < 0)
                        {
                            z = w;
                            r = ra;
                            s = sa;
                        }
                        else
                        {
                            l = i;
                            if (this.e[i] == 0)
                            {
                                cdiv(-ra, -sa, w, q, out this.H[i][n - 1], out this.H[i][n]);
                            }
                            else
                            {
                                // Solve complex equations
                                x = this.H[i][i + 1];
                                y = this.H[i + 1][i];
                                vr = (this.d[i] - p) * (this.d[i] - p) + this.e[i] * this.e[i] - q * q;
                                vi = (this.d[i] - p) * 2 * q;
                                if (vr == 0 & vi == 0)
                                    vr = eps * norm * (System.Math.Abs(w) + System.Math.Abs(q) + System.Math.Abs(x) + System.Math.Abs(y) + System.Math.Abs(z));
                                cdiv(x * r - z * ra + q * sa, x * s - z * sa - q * ra, vr, vi, out this.H[i][n - 1], out this.H[i][n]);
                                if (System.Math.Abs(x) > (System.Math.Abs(z) + System.Math.Abs(q)))
                                {
                                    this.H[i + 1][n - 1] = (-ra - w * this.H[i][n - 1] + q * this.H[i][n]) / x;
                                    this.H[i + 1][n] = (-sa - w * this.H[i][n] - q * this.H[i][n - 1]) / x;
                                }
                                else
                                {
                                    cdiv(-r - y * this.H[i][n - 1], -s - y * this.H[i][n], z, q, out this.H[i + 1][n - 1], out this.H[i + 1][n]);
                                }
                            }

                            // Overflow control
                            t = System.Math.Max(System.Math.Abs(this.H[i][n - 1]), System.Math.Abs(this.H[i][n]));
                            if ((eps * t) * t > 1)
                            {
                                for (int j = i; j <= n; j++)
                                {
                                    this.H[j][n - 1] = this.H[j][n - 1] / t;
                                    this.H[j][n] = this.H[j][n] / t;
                                }
                            }
                        }
                    }
                }
            }

            // Vectors of isolated roots
            for (int i = 0; i < nn; i++)
                if (i < low | i > high)
                    for (int j = i; j < nn; j++)
                        this.V[i][j] = this.H[i][j];

            // Back transformation to get eigenvectors of original matrix
            for (int j = nn - 1; j >= low; j--)
            {
                for (int i = low; i <= high; i++)
                {
                    z = 0;
                    for (int k = low; k <= System.Math.Min(j, high); k++)
                        z = z + this.V[i][k] * this.H[k][j];
                    this.V[i][j] = z;
                }
            }
        }
        #endregion

        /// <summary>
        ///   Reverses the decomposition, reconstructing the original matrix <c>X</c>.
        /// </summary>
        /// 
        public Double[][] Reverse()
        {
            return this.V.DotWithDiagonal(this.d).Divide(this.V);
        }


        #region ICloneable Members

        private JaggedEigenvalueDecomposition()
        {
        }

        /// <summary>
        ///   Creates a new object that is a copy of the current instance.
        /// </summary>
        /// <returns>
        ///   A new object that is a copy of this instance.
        /// </returns>
        public object Clone()
        {
            var clone = new JaggedEigenvalueDecomposition();
            clone.d = (Double[])this.d.Clone();
            clone.e = (Double[])this.e.Clone();
            clone.H = this.H.MemberwiseClone();
            clone.n = this.n;
            clone.ort = (Double[])this.ort;
            clone.symmetric = this.symmetric;
            clone.V = this.V.MemberwiseClone();
            return clone;
        }

        #endregion

    }
}

