// Accord Math Library
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

namespace Openize.Accord.Math.Decompositions
{
    using System;
    using Openize.Accord.Math.Matrix;
    using Openize.Accord.Math.Vector;
    using Openize.Accord.Math;

    /// <summary>
    ///   Determines the Generalized eigenvalues and eigenvectors of two real square matrices.
    /// </summary>
    /// <remarks>
    ///   <para>
    ///     A generalized eigenvalue problem is the problem of finding a vector <c>v</c> that
    ///     obeys <c>A * v = λ * B * v</c> where <c>A</c> and <c>B</c> are matrices. If <c>v</c>
    ///     obeys this equation, with some <c>λ</c>, then we call <c>v</c> the generalized eigenvector
    ///     of <c>A</c> and <c>B</c>, and <c>λ</c> is called the generalized eigenvalue of <c>A</c>
    ///     and <c>B</c> which corresponds to the generalized eigenvector <c>v</c>. The possible
    ///     values of <c>λ</c>, must obey the identity <c>det(A - λ*B) = 0</c>.</para>
    ///   <para>
    ///     Part of this code has been adapted from the original EISPACK routines in Fortran.</para>
    ///  
    ///   <para>
    ///     References:
    ///     <list type="bullet">
    ///       <item><description>
    ///         <a href="http://en.wikipedia.org/wiki/Generalized_eigenvalue_problem#Generalized_eigenvalue_problem">
    ///         http://en.wikipedia.org/wiki/Generalized_eigenvalue_problem#Generalized_eigenvalue_problem</a>
    ///       </description></item>
    ///       <item><description>
    ///         <a href="http://www.netlib.org/eispack/">
    ///         http://www.netlib.org/eispack/</a>
    ///       </description></item>
    ///     </list>
    ///   </para>
    /// </remarks>
    /// 
    /// <example>
    /// // Suppose we have the following 
    /// // matrices A and B shown below:
    /// 
    /// double[][] A = 
    /// {
    ///     new double[] { 1, 2, 3},
    ///     new double[] { 8, 1, 4},
    ///     new double[] { 3, 2, 3}
    /// };
    /// 
    ///  double[][] B = 
    ///  {
    ///     new double[] { 5, 1, 1},
    ///     new double[] { 1, 5, 1},
    ///     new double[] { 1, 1, 5}
    /// };
    /// 
    /// // Now, suppose we would like to find values for λ 
    /// // that are solutions for the equation det(A - λB) = 0
    /// 
    /// // For this, we can use a Generalized Eigendecomposition
    /// var gevd = new JaggedGeneralizedEigenvalueDecomposition(A, B);
    /// 
    /// // Now, if A and B are Hermitian and B is positive
    /// // -definite, then the eigenvalues λ will be real:
    /// double[] lambda = gevd.RealEigenvalues;
    /// 
    /// // Lets check if they are indeed a solution:
    /// for (int i = 0; i &lt; lambda.Length; i++)
    /// {
    ///     // Compute the determinant equation show above
    ///     double det = Matrix.Determinant(A.Subtract(lambda[i].Multiply(B))); // almost zero
    /// }
    /// </example>
    /// 
    public sealed class JaggedGeneralizedEigenvalueDecomposition : ICloneable
    {
        private int n;
        private double[] ar;
        private double[] ai;
        private double[] beta;
        private double[][] Z;

        private const Double eps = 2 * Constants.DoubleEpsilon;



        /// <summary>
        ///   Constructs a new generalized eigenvalue decomposition.
        /// </summary>
        /// 
        /// <param name="a">The first matrix of the (A,B) matrix pencil.</param>
        /// <param name="b">The second matrix of the (A,B) matrix pencil.</param>
        /// <param name="sort">
        ///   Pass <see langword="true"/> to sort the eigenvalues and eigenvectors at the end
        ///   of the decomposition.</param>
        ///
        public JaggedGeneralizedEigenvalueDecomposition(double[][] a, double[][] b, bool sort = false)
        {
            if (a == null)
                throw new ArgumentNullException("a", "Matrix A cannot be null.");

            if (b == null)
                throw new ArgumentNullException("b", "Matrix B cannot be null.");

            if (a.Rows() != a.Columns())
                throw new ArgumentException("Matrix is not a square matrix.", "a");

            if (b.Rows() != b.Columns())
                throw new ArgumentException("Matrix is not a square matrix.", "b");

            if (a.Rows() != b.Rows() || a.Columns() != b.Columns())
                throw new ArgumentException("Matrix dimensions do not match", "b");

            this.n = a.Rows();
            this.ar = new double[this.n];
            this.ai = new double[this.n];
            this.beta = new double[this.n];
            this.Z = Jagged.Zeros(this.n, this.n);
            var A = Vector.Copy(a);
            var B = Vector.Copy(b);

            bool matz = true;
            int ierr = 0;


            // reduces A to upper Hessenberg form and B to upper
            // triangular form using orthogonal transformations
            qzhes(this.n, A, B, matz, this.Z);

            // reduces the Hessenberg matrix A to quasi-triangular form
            // using orthogonal transformations while maintaining the
            // triangular form of the B matrix.
            qzit(this.n, A, B, Constants.DoubleEpsilon, matz, this.Z, ref ierr);

            // reduces the quasi-triangular matrix further, so that any
            // remaining 2-by-2 blocks correspond to pairs of complex
            // eigenvalues, and returns quantities whose ratios give the
            // generalized eigenvalues.
            qzval(this.n, A, B, this.ar, this.ai, this.beta, matz, this.Z);

            // computes the eigenvectors of the triangular problem and
            // transforms the results back to the original coordinate system.
            qzvec(this.n, A, B, this.ar, this.ai, this.beta, this.Z);

            if (sort)
            {
                // Sort eigenvalues and vectors in descending order
                var idx = Vector.Range(this.n);
                Array.Sort(idx, (i, j) =>
                {
                    if (Math.Abs(this.ar[i]) == Math.Abs(this.ar[j]))
                        return -Math.Abs(this.ai[i]).CompareTo(Math.Abs(this.ai[j]));
                    return -Math.Abs(this.ar[i]).CompareTo(Math.Abs(this.ar[j]));
                });

                this.ar = this.ar.Get(idx);
                this.ai = this.ai.Get(idx);
                this.beta = this.beta.Get(idx);
                this.Z = this.Z.Get(null, idx);
            }
        }

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
                Double tol = this.n * this.ar[0] * eps;

                int r = 0;
                for (int i = 0; i < this.ar.Length; i++)
                    if (this.ar[i] > tol) r++;

                return r;
            }
        }

        /// <summary>Returns the real parts of the alpha values.</summary>
        public double[] RealAlphas
        {
            get { return this.ar; }
        }

        /// <summary>Returns the imaginary parts of the alpha values.</summary>
        public double[] ImaginaryAlphas
        {
            get { return this.ai; }
        }

        /// <summary>Returns the beta values.</summary>
        public double[] Betas
        {
            get { return this.beta; }
        }

        /// <summary>
        ///   Returns true if matrix B is singular.
        /// </summary>
        /// <remarks>
        ///   This method checks if any of the generated betas is zero. It
        ///   does not says that the problem is singular, but only that one
        ///   of the matrices of the pencil (A,B) is singular.
        /// </remarks>
        public bool IsSingular
        {
            get
            {
                for (int i = 0; i < this.n; i++)
                    if (this.beta[i] == 0)
                        return true;
                return false;
            }
        }

        /// <summary>
        ///   Returns true if the eigenvalue problem is degenerate (ill-posed).
        /// </summary>
        public bool IsDegenerate
        {
            get
            {
                for (int i = 0; i < this.n; i++)
                    if (this.beta[i] == 0 && this.ar[i] == 0)
                        return true;
                return false;
            }
        }

        /// <summary>Returns the real parts of the eigenvalues.</summary>
        /// <remarks>
        ///   The eigenvalues are computed using the ratio alpha[i]/beta[i],
        ///   which can lead to valid, but infinite eigenvalues.
        /// </remarks>
        public double[] RealEigenvalues
        {
            get
            {
                // ((alfr+i*alfi)/beta)
                double[] eval = new double[this.n];
                for (int i = 0; i < this.n; i++)
                    eval[i] = this.ar[i] / this.beta[i];
                return eval;
            }
        }

        /// <summary>Returns the imaginary parts of the eigenvalues.</summary>	
        /// <remarks>
        ///   The eigenvalues are computed using the ratio alpha[i]/beta[i],
        ///   which can lead to valid, but infinite eigenvalues.
        /// </remarks>
        public double[] ImaginaryEigenvalues
        {
            get
            {
                // ((alfr+i*alfi)/beta)
                double[] eval = new double[this.n];
                for (int i = 0; i < this.n; i++)
                    eval[i] = this.ai[i] / this.beta[i];
                return eval;
            }
        }

        /// <summary>Returns the eigenvector matrix.</summary>
        public double[][] Eigenvectors
        {
            get
            {
                return this.Z;
            }
        }

        /// <summary>Returns the block diagonal eigenvalue matrix.</summary>
        public double[][] DiagonalMatrix
        {
            get
            {
                double[][] x = Jagged.Zeros(this.n, this.n);

                for (int i = 0; i < this.n; i++)
                {
                    for (int j = 0; j < this.n; j++)
                        x[i][j] = 0.0;

                    x[i][i] = this.ar[i] / this.beta[i];
                    if (this.ai[i] > 0)
                        x[i][i + 1] = this.ai[i] / this.beta[i];
                    else if (this.ai[i] < 0)
                        x[i][i - 1] = this.ai[i] / this.beta[i];
                }

                return x;
            }
        }



        #region EISPACK Routines
        /// <summary>
        ///   Adaptation of the original Fortran QZHES routine from EISPACK.
        /// </summary>
        /// <remarks>
        ///   This subroutine is the first step of the qz algorithm
        ///   for solving generalized matrix eigenvalue problems,
        ///   Siam J. Numer. anal. 10, 241-256(1973) by Moler and Stewart.
        ///
        ///   This subroutine accepts a pair of real general matrices and
        ///   reduces one of them to upper Hessenberg form and the other
        ///   to upper triangular form using orthogonal transformations.
        ///   it is usually followed by  qzit,  qzval  and, possibly,  qzvec.
        ///   
        ///   For the full documentation, please check the original function.
        /// </remarks>
        private static int qzhes(int n, double[][] a, double[][] b, bool matz, double[][] z)
        {
            int i, j, k, l;
            double r, s, t;
            int l1;
            double u1, u2, v1, v2;
            int lb, nk1;
            double rho;


            if (matz)
            {
                // If we are interested in computing the
                //  eigenvectors, set Z to identity(n,n)
                for (j = 0; j < n; ++j)
                {
                    for (i = 0; i < n; ++i)
                        z[i][j] = 0.0;
                    z[j][j] = 1.0;
                }
            }

            // Reduce b to upper triangular form
            if (n <= 1) return 0;
            for (l = 0; l < n - 1; ++l)
            {
                l1 = l + 1;
                s = 0.0;

                for (i = l1; i < n; ++i)
                    s += (System.Math.Abs(b[i][l]));

                if (s == 0.0) continue;
                s += (System.Math.Abs(b[l][l]));
                r = 0.0;

                for (i = l; i < n; ++i)
                {
                    // Computing 2nd power
                    b[i][l] /= s;
                    r += b[i][l] * b[i][l];
                }

                r = Special.Sign(System.Math.Sqrt(r), b[l][l]);
                b[l][l] += r;
                rho = r * b[l][l];

                for (j = l1; j < n; ++j)
                {
                    t = 0.0;
                    for (i = l; i < n; ++i)
                        t += b[i][l] * b[i][j];
                    t = -t / rho;
                    for (i = l; i < n; ++i)
                        b[i][j] += t * b[i][l];
                }

                for (j = 0; j < n; ++j)
                {
                    t = 0.0;
                    for (i = l; i < n; ++i)
                        t += b[i][l] * a[i][j];
                    t = -t / rho;
                    for (i = l; i < n; ++i)
                        a[i][j] += t * b[i][l];
                }

                b[l][l] = -s * r;
                for (i = l1; i < n; ++i)
                    b[i][l] = 0.0;
            }

            // Reduce a to upper Hessenberg form, while keeping b triangular
            if (n == 2) 
                return 0;

            for (k = 0; k < n - 2; ++k)
            {
                nk1 = n - 2 - k;

                // for l=n-1 step -1 until k+1 do
                for (lb = 0; lb < nk1; ++lb)
                {
                    l = n - lb - 2;
                    l1 = l + 1;

                    // Zero a(l+1,k)
                    s = (System.Math.Abs(a[l][k])) + (System.Math.Abs(a[l1][k]));

                    if (s == 0.0) continue;
                    u1 = a[l][k] / s;
                    u2 = a[l1][k] / s;
                    r = Special.Sign(Math.Sqrt(u1 * u1 + u2 * u2), u1);
                    v1 = -(u1 + r) / r;
                    v2 = -u2 / r;
                    u2 = v2 / v1;

                    for (j = k; j < n; ++j)
                    {
                        t = a[l][j] + u2 * a[l1][j];
                        a[l][j] += t * v1;
                        a[l1][j] += t * v2;
                    }

                    a[l1][k] = 0.0;

                    for (j = l; j < n; ++j)
                    {
                        t = b[l][j] + u2 * b[l1][j];
                        b[l][j] += t * v1;
                        b[l1][j] += t * v2;
                    }

                    // Zero b(l+1,l)
                    s = (System.Math.Abs(b[l1][l1])) + (System.Math.Abs(b[l1][l]));

                    if (s == 0.0) continue;
                    u1 = b[l1][l1] / s;
                    u2 = b[l1][l] / s;
                    r = Special.Sign(Math.Sqrt(u1 * u1 + u2 * u2), u1);
                    v1 = -(u1 + r) / r;
                    v2 = -u2 / r;
                    u2 = v2 / v1;

                    for (i = 0; i <= l1; ++i)
                    {
                        t = b[i][l1] + u2 * b[i][l];
                        b[i][l1] += t * v1;
                        b[i][l] += t * v2;
                    }

                    b[l1][l] = 0.0;

                    for (i = 0; i < n; ++i)
                    {
                        t = a[i][l1] + u2 * a[i][l];
                        a[i][l1] += t * v1;
                        a[i][l] += t * v2;
                    }

                    if (matz)
                    {
                        for (i = 0; i < n; ++i)
                        {
                            t = z[i][l1] + u2 * z[i][l];
                            z[i][l1] += t * v1;
                            z[i][l] += t * v2;
                        }
                    }
                }
            }

            return 0;
        }

        /// <summary>
        ///   Adaptation of the original Fortran QZIT routine from EISPACK.
        /// </summary>
        /// <remarks>
        ///   This subroutine is the second step of the qz algorithm
        ///   for solving generalized matrix eigenvalue problems,
        ///   Siam J. Numer. anal. 10, 241-256(1973) by Moler and Stewart,
        ///   as modified in technical note nasa tn d-7305(1973) by ward.
        ///   
        ///   This subroutine accepts a pair of real matrices, one of them
        ///   in upper Hessenberg form and the other in upper triangular form.
        ///   it reduces the Hessenberg matrix to quasi-triangular form using
        ///   orthogonal transformations while maintaining the triangular form
        ///   of the other matrix.  it is usually preceded by  qzhes  and
        ///   followed by  qzval  and, possibly,  qzvec.
        ///   
        ///   For the full documentation, please check the original function.
        /// </remarks>
        private static int qzit(int n, double[][] a, double[][] b, double eps1, bool matz, double[][] z, ref int ierr)
        {
            int i, j, k, l = 0;
            double r, s, t, a1, a2, a3 = 0;
            int k1, k2, l1, ll;
            double u1, u2, u3;
            double v1, v2, v3;
            double a11, a12, a21, a22, a33, a34, a43, a44;
            double b11, b12, b22, b33, b34, b44;
            int na, en, ld;
            double ep;
            double sh = 0;
            int km1, lm1 = 0;
            double ani, bni;
            int ish, itn, its, enm2, lor1;
            double epsa, epsb, anorm = 0, bnorm = 0;
            int enorn;
            bool notlas;


            ierr = 0;

            #region Compute epsa and epsb
            for (i = 0; i < n; ++i)
            {
                ani = 0.0;
                bni = 0.0;

                if (i != 0)
                    ani = (Math.Abs(a[i][(i - 1)]));

                for (j = i; j < n; ++j)
                {
                    ani += Math.Abs(a[i][j]);
                    bni += Math.Abs(b[i][j]);
                }

                if (ani > anorm) anorm = ani;
                if (bni > bnorm) bnorm = bni;
            }

            if (anorm == 0.0) anorm = 1.0;
            if (bnorm == 0.0) bnorm = 1.0;

            ep = eps1;
            if (ep == 0.0)
            {
                // Use round-off level if eps1 is zero
                ep = Special.Epslon(1.0);
            }

            epsa = ep * anorm;
            epsb = ep * bnorm;
            #endregion


            // Reduce a to quasi-triangular form][while keeping b triangular
            lor1 = 0;
            enorn = n;
            en = n - 1;
            itn = n * 30;

        // Begin QZ step
        L60:
            if (en <= 1) goto L1001;
            if (!matz) enorn = en + 1;

            its = 0;
            na = en - 1;
            enm2 = na;

        L70:
            ish = 2;
            // Check for convergence or reducibility.
            for (ll = 0; ll <= en; ++ll)
            {
                lm1 = en - ll - 1;
                l = lm1 + 1;

                if (l + 1 == 1)
                    goto L95;

                if ((Math.Abs(a[l][lm1])) <= epsa)
                    break;
            }

        L90:
            a[l][lm1] = 0.0;
            if (l < na) goto L95;

            // 1-by-1 or 2-by-2 block isolated
            en = lm1;
            goto L60;

        // Check for small top of b 
        L95:
            ld = l;

        L100:
            l1 = l + 1;
            b11 = b[l][l];

            if (Math.Abs(b11) > epsb) goto L120;

            b[l][l] = 0.0;
            s = (Math.Abs(a[l][l]) + Math.Abs(a[l1][l]));
            u1 = a[l][l] / s;
            u2 = a[l1][l] / s;
            r = Special.Sign(Math.Sqrt(u1 * u1 + u2 * u2), u1);
            v1 = -(u1 + r) / r;
            v2 = -u2 / r;
            u2 = v2 / v1;

            for (j = l; j < enorn; ++j)
            {
                t = a[l][j] + u2 * a[l1][j];
                a[l][j] += t * v1;
                a[l1][j] += t * v2;

                t = b[l][j] + u2 * b[l1][j];
                b[l][j] += t * v1;
                b[l1][j] += t * v2;
            }

            if (l != 0)
                a[l][lm1] = -a[l][lm1];

            lm1 = l;
            l = l1;
            goto L90;

        L120:
            a11 = a[l][l] / b11;
            a21 = a[l1][l] / b11;
            if (ish == 1) goto L140;

            // Iteration strategy
            if (itn == 0) goto L1000;
            if (its == 10) goto L155;

            // Determine type of shift
            b22 = b[l1][l1];
            if (Math.Abs(b22) < epsb) b22 = epsb;
            b33 = b[na][na];
            if (Math.Abs(b33) < epsb) b33 = epsb;
            b44 = b[en][en];
            if (Math.Abs(b44) < epsb) b44 = epsb;
            a33 = a[na][na] / b33;
            a34 = a[na][en] / b44;
            a43 = a[en][na] / b33;
            a44 = a[en][en] / b44;
            b34 = b[na][en] / b44;
            t = (a43 * b34 - a33 - a44) * .5;
            r = t * t + a34 * a43 - a33 * a44;
            if (r < 0.0) goto L150;

            // Determine single shift zero-th column of a
            ish = 1;
            r = Math.Sqrt(r);
            sh = -t + r;
            s = -t - r;
            if (Math.Abs(s - a44) < Math.Abs(sh - a44))
                sh = s;

            // Look for two consecutive small sub-diagonal elements of a.
            for (ll = ld; ll + 1 <= enm2; ++ll)
            {
                l = enm2 + ld - ll - 1;

                if (l == ld)
                    goto L140;

                lm1 = l - 1;
                l1 = l + 1;
                t = a[l + 1][l + 1];

                if (Math.Abs(b[l][l]) > epsb)
                    t -= sh * b[l][l];

                if (Math.Abs(a[l][lm1]) <= (Math.Abs(t / a[l1][l])) * epsa)
                    goto L100;
            }

        L140:
            a1 = a11 - sh;
            a2 = a21;
            if (l != ld)
                a[l][lm1] = -a[l][lm1];
            goto L160;

        // Determine double shift zero-th column of a
        L150:
            a12 = a[l][l1] / b22;
            a22 = a[l1][l1] / b22;
            b12 = b[l][l1] / b22;
            a1 = ((a33 - a11) * (a44 - a11) - a34 * a43 + a43 * b34 * a11) / a21 + a12 - a11 * b12;
            a2 = a22 - a11 - a21 * b12 - (a33 - a11) - (a44 - a11) + a43 * b34;
            a3 = a[l1 + 1][l1] / b22;
            goto L160;

        // Ad hoc shift
        L155:
            a1 = 0.0;
            a2 = 1.0;
            a3 = 1.1605;

        L160:
            ++its;
            --itn;

            if (!matz) lor1 = ld;

            // Main loop
            for (k = l; k <= na; ++k)
            {
                notlas = k != na && ish == 2;
                k1 = k + 1;
                k2 = k + 2;

                km1 = Math.Max(k, l + 1) - 1; // Computing MAX
                ll = Math.Min(en, k1 + ish);  // Computing MIN

                if (notlas) goto L190;

                // Zero a(k+1,k-1)
                if (k == l) goto L170;
                a1 = a[k][km1];
                a2 = a[k1][km1];

            L170:
                s = Math.Abs(a1) + Math.Abs(a2);
                if (s == 0.0) goto L70;
                u1 = a1 / s;
                u2 = a2 / s;
                r = Special.Sign(Math.Sqrt(u1 * u1 + u2 * u2), u1);
                v1 = -(u1 + r) / r;
                v2 = -u2 / r;
                u2 = v2 / v1;

                for (j = km1; j < enorn; ++j)
                {
                    t = a[k][j] + u2 * a[k1][j];
                    a[k][j] += t * v1;
                    a[k1][j] += t * v2;

                    t = b[k][j] + u2 * b[k1][j];
                    b[k][j] += t * v1;
                    b[k1][j] += t * v2;
                }

                if (k != l)
                    a[k1][km1] = 0.0;
                goto L240;

                // Zero a(k+1,k-1) and a(k+2,k-1)
            L190:
                if (k == l) goto L200;
                a1 = a[k][km1];
                a2 = a[k1][km1];
                a3 = a[k2][km1];

            L200:
                s = Math.Abs(a1) + Math.Abs(a2) + Math.Abs(a3);
                if (s == 0.0) goto L260;
                u1 = a1 / s;
                u2 = a2 / s;
                u3 = a3 / s;
                r = Special.Sign(Math.Sqrt(u1 * u1 + u2 * u2 + u3 * u3), u1);
                v1 = -(u1 + r) / r;
                v2 = -u2 / r;
                v3 = -u3 / r;
                u2 = v2 / v1;
                u3 = v3 / v1;

                for (j = km1; j < enorn; ++j)
                {
                    t = a[k][j] + u2 * a[k1][j] + u3 * a[k2][j];
                    a[k][j] += t * v1;
                    a[k1][j] += t * v2;
                    a[k2][j] += t * v3;

                    t = b[k][j] + u2 * b[k1][j] + u3 * b[k2][j];
                    b[k][j] += t * v1;
                    b[k1][j] += t * v2;
                    b[k2][j] += t * v3;
                }

                if (k == l) goto L220;
                a[k1][km1] = 0.0;
                a[k2][km1] = 0.0;

            // Zero b(k+2,k+1) and b(k+2,k)
            L220:
                s = (Math.Abs(b[k2][k2])) + (Math.Abs(b[k2][k1])) + (Math.Abs(b[k2][k]));
                if (s == 0.0) goto L240;
                u1 = b[k2][k2] / s;
                u2 = b[k2][k1] / s;
                u3 = b[k2][k] / s;
                r = Special.Sign(Math.Sqrt(u1 * u1 + u2 * u2 + u3 * u3), u1);
                v1 = -(u1 + r) / r;
                v2 = -u2 / r;
                v3 = -u3 / r;
                u2 = v2 / v1;
                u3 = v3 / v1;

                for (i = lor1; i < ll + 1; ++i)
                {
                    t = a[i][k2] + u2 * a[i][k1] + u3 * a[i][k];
                    a[i][k2] += t * v1;
                    a[i][k1] += t * v2;
                    a[i][k] += t * v3;

                    t = b[i][k2] + u2 * b[i][k1] + u3 * b[i][k];
                    b[i][k2] += t * v1;
                    b[i][k1] += t * v2;
                    b[i][k] += t * v3;
                }

                b[k2][k] = 0.0;
                b[k2][k1] = 0.0;

                if (matz)
                {
                    for (i = 0; i < n; ++i)
                    {
                        t = z[i][k2] + u2 * z[i][k1] + u3 * z[i][k];
                        z[i][k2] += t * v1;
                        z[i][k1] += t * v2;
                        z[i][k] += t * v3;
                    }
                }

            // Zero b(k+1,k)
            L240:
                s = (Math.Abs(b[k1][k1])) + (Math.Abs(b[k1][k]));
                if (s == 0.0) goto L260;
                u1 = b[k1][k1] / s;
                u2 = b[k1][k] / s;
                r = Special.Sign(Math.Sqrt(u1 * u1 + u2 * u2), u1);
                v1 = -(u1 + r) / r;
                v2 = -u2 / r;
                u2 = v2 / v1;

                for (i = lor1; i < ll + 1; ++i)
                {
                    t = a[i][k1] + u2 * a[i][k];
                    a[i][k1] += t * v1;
                    a[i][k] += t * v2;

                    t = b[i][k1] + u2 * b[i][k];
                    b[i][k1] += t * v1;
                    b[i][k] += t * v2;
                }

                b[k1][k] = 0.0;

                if (matz)
                {
                    for (i = 0; i < n; ++i)
                    {
                        t = z[i][k1] + u2 * z[i][k];
                        z[i][k1] += t * v1;
                        z[i][k] += t * v2;
                    }
                }

            L260:
                ;
            }

            goto L70; // End QZ step

        // Set error -- all eigenvalues have not converged after 30*n iterations
        L1000:
            ierr = en + 1;

        // Save epsb for use by qzval and qzvec
        L1001:
            if (n > 1)
                b[n - 1][0] = epsb;
            return 0;
        }

        /// <summary>
        ///   Adaptation of the original Fortran QZVAL routine from EISPACK.
        /// </summary>
        /// <remarks>
        ///   This subroutine is the third step of the qz algorithm
        ///   for solving generalized matrix eigenvalue problems,
        ///   Siam J. Numer. anal. 10, 241-256(1973) by Moler and Stewart.
        ///   
        ///   This subroutine accepts a pair of real matrices, one of them
        ///   in quasi-triangular form and the other in upper triangular form.
        ///   it reduces the quasi-triangular matrix further, so that any
        ///   remaining 2-by-2 blocks correspond to pairs of complex
        ///   Eigenvalues, and returns quantities whose ratios give the
        ///   generalized eigenvalues.  it is usually preceded by  qzhes
        ///   and  qzit  and may be followed by  qzvec.
        ///   
        ///   For the full documentation, please check the original function.
        /// </remarks>
        private static int qzval(int n, double[][] a, double[][] b, double[] alfr, double[] alfi, double[] beta, bool matz, double[][] z)
        {
            int i, j;
            int na, en, nn;
            double c, d, e = 0;
            double r, s, t;
            double a1, a2, u1, u2, v1, v2;
            double a11, a12, a21, a22;
            double b11, b12, b22;
            double di, ei;
            double an = 0, bn;
            double cq, dr;
            double cz, ti, tr;
            double a1i, a2i, a11i, a12i, a22i, a11r, a12r, a22r;
            double sqi, ssi, sqr, szi, ssr, szr;

            double epsb = b[n - 1][0];
            int isw = 1;


            // Find eigenvalues of quasi-triangular matrices.
            for (nn = 0; nn < n; ++nn)
            {
                en = n - nn - 1;
                na = en - 1;

                if (isw == 2) goto L505;
                if (en == 0) goto L410;
                if (a[en][na] != 0.0) goto L420;

            // 1-by-1 block][one real root
            L410:
                alfr[en] = a[en][en];
                if (b[en][en] < 0.0)
                {
                    alfr[en] = -alfr[en];
                }
                beta[en] = (Math.Abs(b[en][en]));
                alfi[en] = 0.0;
                goto L510;

            // 2-by-2 block
            L420:
                if (Math.Abs(b[na][na]) <= epsb) goto L455;
                if (Math.Abs(b[en][en]) > epsb) goto L430;
                a1 = a[en][en];
                a2 = a[en][na];
                bn = 0.0;
                goto L435;

            L430:
                an = Math.Abs(a[na][na]) + Math.Abs(a[na][en]) + Math.Abs(a[en][na]) + Math.Abs(a[en][en]);
                bn = Math.Abs(b[na][na]) + Math.Abs(b[na][en]) + Math.Abs(b[en][en]);
                a11 = a[na][na] / an;
                a12 = a[na][en] / an;
                a21 = a[en][na] / an;
                a22 = a[en][en] / an;
                b11 = b[na][na] / bn;
                b12 = b[na][en] / bn;
                b22 = b[en][en] / bn;
                e = a11 / b11;
                ei = a22 / b22;
                s = a21 / (b11 * b22);
                t = (a22 - e * b22) / b22;

                if (Math.Abs(e) <= Math.Abs(ei))
                    goto L431;

                e = ei;
                t = (a11 - e * b11) / b11;

            L431:
                c = (t - s * b12) * .5;
                d = c * c + s * (a12 - e * b12);
                if (d < 0.0) goto L480;

                // Two real roots. Zero both a(en,na) and b(en,na)
                e += c + Special.Sign(Math.Sqrt(d), c);
                a11 -= e * b11;
                a12 -= e * b12;
                a22 -= e * b22;

                if (Math.Abs(a11) + Math.Abs(a12) < Math.Abs(a21) + Math.Abs(a22))
                    goto L432;

                a1 = a12;
                a2 = a11;
                goto L435;

            L432:
                a1 = a22;
                a2 = a21;

            // Choose and apply real z
            L435:
                s = Math.Abs(a1) + Math.Abs(a2);
                u1 = a1 / s;
                u2 = a2 / s;
                r = Special.Sign(Math.Sqrt(u1 * u1 + u2 * u2), u1);
                v1 = -(u1 + r) / r;
                v2 = -u2 / r;
                u2 = v2 / v1;

                for (i = 0; i <= en; ++i)
                {
                    t = a[i][en] + u2 * a[i][na];
                    a[i][en] += t * v1;
                    a[i][na] += t * v2;

                    t = b[i][en] + u2 * b[i][na];
                    b[i][en] += t * v1;
                    b[i][na] += t * v2;
                }

                if (matz)
                {
                    for (i = 0; i < n; ++i)
                    {
                        t = z[i][en] + u2 * z[i][na];
                        z[i][en] += t * v1;
                        z[i][na] += t * v2;
                    }
                }

                if (bn == 0.0) goto L475;
                if (an < System.Math.Abs(e) * bn) goto L455;
                a1 = b[na][na];
                a2 = b[en][na];
                goto L460;

            L455:
                a1 = a[na][na];
                a2 = a[en][na];

            // Choose and apply real q
            L460:
                s = System.Math.Abs(a1) + System.Math.Abs(a2);
                if (s == 0.0) goto L475;
                u1 = a1 / s;
                u2 = a2 / s;
                r = Special.Sign(Math.Sqrt(u1 * u1 + u2 * u2), u1);
                v1 = -(u1 + r) / r;
                v2 = -u2 / r;
                u2 = v2 / v1;

                for (j = na; j < n; ++j)
                {
                    t = a[na][j] + u2 * a[en][j];
                    a[na][j] += t * v1;
                    a[en][j] += t * v2;

                    t = b[na][j] + u2 * b[en][j];
                    b[na][j] += t * v1;
                    b[en][j] += t * v2;
                }

            L475:
                a[en][na] = 0.0;
                b[en][na] = 0.0;
                alfr[na] = a[na][na];
                alfr[en] = a[en][en];

                if (b[na][na] < 0.0)
                    alfr[na] = -alfr[na];

                if (b[en][en] < 0.0)
                    alfr[en] = -alfr[en];

                beta[na] = (System.Math.Abs(b[na][na]));
                beta[en] = (System.Math.Abs(b[en][en]));
                alfi[en] = 0.0;
                alfi[na] = 0.0;
                goto L505;

                // Two complex roots
            L480:
                e += c;
                ei = System.Math.Sqrt(-d);
                a11r = a11 - e * b11;
                a11i = ei * b11;
                a12r = a12 - e * b12;
                a12i = ei * b12;
                a22r = a22 - e * b22;
                a22i = ei * b22;

                if (System.Math.Abs(a11r) + System.Math.Abs(a11i) +
                    System.Math.Abs(a12r) + System.Math.Abs(a12i) <
                    System.Math.Abs(a21) + System.Math.Abs(a22r)
                    + System.Math.Abs(a22i))
                    goto L482;

                a1 = a12r;
                a1i = a12i;
                a2 = -a11r;
                a2i = -a11i;
                goto L485;

            L482:
                a1 = a22r;
                a1i = a22i;
                a2 = -a21;
                a2i = 0.0;

            // Choose complex z
            L485:
                cz = System.Math.Sqrt(a1 * a1 + a1i * a1i);
                if (cz == 0.0) goto L487;
                szr = (a1 * a2 + a1i * a2i) / cz;
                szi = (a1 * a2i - a1i * a2) / cz;
                r = System.Math.Sqrt(cz * cz + szr * szr + szi * szi);
                cz /= r;
                szr /= r;
                szi /= r;
                goto L490;

            L487:
                szr = 1.0;
                szi = 0.0;

            L490:
                if (an < (System.Math.Abs(e) + ei) * bn) goto L492;
                a1 = cz * b11 + szr * b12;
                a1i = szi * b12;
                a2 = szr * b22;
                a2i = szi * b22;
                goto L495;

            L492:
                a1 = cz * a11 + szr * a12;
                a1i = szi * a12;
                a2 = cz * a21 + szr * a22;
                a2i = szi * a22;

            // Choose complex q
            L495:
                cq = System.Math.Sqrt(a1 * a1 + a1i * a1i);
                if (cq == 0.0) goto L497;
                sqr = (a1 * a2 + a1i * a2i) / cq;
                sqi = (a1 * a2i - a1i * a2) / cq;
                r = System.Math.Sqrt(cq * cq + sqr * sqr + sqi * sqi);
                cq /= r;
                sqr /= r;
                sqi /= r;
                goto L500;

            L497:
                sqr = 1.0;
                sqi = 0.0;

            // Compute diagonal elements that would result if transformations were applied
            L500:
                ssr = sqr * szr + sqi * szi;
                ssi = sqr * szi - sqi * szr;
                i = 0;
                tr = cq * cz * a11 + cq * szr * a12 + sqr * cz * a21 + ssr * a22;
                ti = cq * szi * a12 - sqi * cz * a21 + ssi * a22;
                dr = cq * cz * b11 + cq * szr * b12 + ssr * b22;
                di = cq * szi * b12 + ssi * b22;
                goto L503;

            L502:
                i = 1;
                tr = ssr * a11 - sqr * cz * a12 - cq * szr * a21 + cq * cz * a22;
                ti = -ssi * a11 - sqi * cz * a12 + cq * szi * a21;
                dr = ssr * b11 - sqr * cz * b12 + cq * cz * b22;
                di = -ssi * b11 - sqi * cz * b12;

            L503:
                t = ti * dr - tr * di;
                j = na;

                if (t < 0.0)
                    j = en;

                r = Math.Sqrt(dr * dr + di * di);
                beta[j] = bn * r;
                alfr[j] = an * (tr * dr + ti * di) / r;
                alfi[j] = an * t / r;
                if (i == 0) goto L502;

            L505:
                isw = 3 - isw;

            L510:
                ;
            }

            b[n - 1][0] = epsb;

            return 0;
        }

        /// <summary>
        ///   Adaptation of the original Fortran QZVEC routine from EISPACK.
        /// </summary>
        /// <remarks>
        ///   This subroutine is the optional fourth step of the qz algorithm
        ///   for solving generalized matrix eigenvalue problems,
        ///   Siam J. Numer. anal. 10, 241-256(1973) by Moler and Stewart.
        ///   
        ///   This subroutine accepts a pair of real matrices, one of them in
        ///   quasi-triangular form (in which each 2-by-2 block corresponds to
        ///   a pair of complex eigenvalues) and the other in upper triangular
        ///   form.  It computes the eigenvectors of the triangular problem and
        ///   transforms the results back to the original coordinate system.
        ///   it is usually preceded by  qzhes,  qzit, and  qzval.
        ///   
        ///   For the full documentation, please check the original function.
        /// </remarks>
        private static int qzvec(int n, double[][] a, double[][] b, double[] alfr, double[] alfi, double[] beta, double[][] z)
        {
            int i, j, k, m;
            int na, ii, en, jj, nn, enm2;
            double d, q;
            double r = 0, s = 0, t, w, x = 0, y, t1, t2, w1, x1 = 0, z1 = 0, di;
            double ra, dr, sa;
            double ti, rr, tr, zz = 0;
            double alfm, almi, betm, almr;

            double epsb = b[n - 1][0];
            int isw = 1;


            // for en=n step -1 until 1 do --
            for (nn = 0; nn < n; ++nn)
            {
                en = n - nn - 1;
                na = en - 1;
                if (isw == 2) goto L795;
                if (alfi[en] != 0.0) goto L710;

                // Real vector
                m = en;
                b[en][en] = 1.0;
                if (na == -1) goto L800;
                alfm = alfr[m];
                betm = beta[m];

                // for i=en-1 step -1 until 1 do --
                for (ii = 0; ii <= na; ++ii)
                {
                    i = en - ii - 1;
                    w = betm * a[i][i] - alfm * b[i][i];
                    r = 0.0;

                    for (j = m; j <= en; ++j)
                        r += (betm * a[i][j] - alfm * b[i][j]) * b[j][en];

                    if (i == 0 || isw == 2)
                        goto L630;

                    if (betm * a[i][i - 1] == 0.0)
                        goto L630;

                    zz = w;
                    s = r;
                    goto L690;

                L630:
                    m = i;
                    if (isw == 2) goto L640;

                    // Real 1-by-1 block
                    t = w;
                    if (w == 0.0)
                        t = epsb;
                    b[i][en] = -r / t;
                    goto L700;

                // Real 2-by-2 block
                L640:
                    x = betm * a[i][i + 1] - alfm * b[i][i + 1];
                    y = betm * a[i + 1][i];
                    q = w * zz - x * y;
                    t = (x * s - zz * r) / q;
                    b[i][en] = t;
                    if (Math.Abs(x) <= Math.Abs(zz)) goto L650;
                    b[i + 1][en] = (-r - w * t) / x;
                    goto L690;

                L650:
                    b[i + 1][en] = (-s - y * t) / zz;

                L690:
                    isw = 3 - isw;

                L700:
                    ;
                }
                // End real vector
                goto L800;

            // Complex vector
            L710:
                m = na;
                almr = alfr[m];
                almi = alfi[m];
                betm = beta[m];

                // last vector component chosen imaginary so that eigenvector matrix is triangular
                y = betm * a[en][na];
                b[na][na] = -almi * b[en][en] / y;
                b[na][en] = (almr * b[en][en] - betm * a[en][en]) / y;
                b[en][na] = 0.0;
                b[en][en] = 1.0;
                enm2 = na;
                if (enm2 == 0) goto L795;

                // for i=en-2 step -1 until 1 do --
                for (ii = 0; ii < enm2; ++ii)
                {
                    i = na - ii - 1;
                    w = betm * a[i][i] - almr * b[i][i];
                    w1 = -almi * b[i][i];
                    ra = 0.0;
                    sa = 0.0;

                    for (j = m; j <= en; ++j)
                    {
                        x = betm * a[i][j] - almr * b[i][j];
                        x1 = -almi * b[i][j];
                        ra = ra + x * b[j][na] - x1 * b[j][en];
                        sa = sa + x * b[j][en] + x1 * b[j][na];
                    }

                    if (i == 0 || isw == 2) goto L770;
                    if (betm * a[i][i - 1] == 0.0) goto L770;

                    zz = w;
                    z1 = w1;
                    r = ra;
                    s = sa;
                    isw = 2;
                    goto L790;

                L770:
                    m = i;
                    if (isw == 2) goto L780;

                    // Complex 1-by-1 block 
                    tr = -ra;
                    ti = -sa;

                L773:
                    dr = w;
                    di = w1;

                    // Complex divide (t1,t2) = (tr,ti) / (dr,di)
                L775:
                    if (Math.Abs(di) > Math.Abs(dr)) goto L777;
                    rr = di / dr;
                    d = dr + di * rr;
                    t1 = (tr + ti * rr) / d;
                    t2 = (ti - tr * rr) / d;

                    switch (isw)
                    {
                        case 1: goto L787;
                        case 2: goto L782;
                    }

                L777:
                    rr = dr / di;
                    d = dr * rr + di;
                    t1 = (tr * rr + ti) / d;
                    t2 = (ti * rr - tr) / d;
                    switch (isw)
                    {
                        case 1: goto L787;
                        case 2: goto L782;
                    }

                   // Complex 2-by-2 block 
                L780:
                    x = betm * a[i][i + 1] - almr * b[i][i + 1];
                    x1 = -almi * b[i][i + 1];
                    y = betm * a[i + 1][i];
                    tr = y * ra - w * r + w1 * s;
                    ti = y * sa - w * s - w1 * r;
                    dr = w * zz - w1 * z1 - x * y;
                    di = w * z1 + w1 * zz - x1 * y;
                    if (dr == 0.0 && di == 0.0)
                        dr = epsb;
                    goto L775;

                L782:
                    b[i + 1][na] = t1;
                    b[i + 1][en] = t2;
                    isw = 1;
                    if (Math.Abs(y) > Math.Abs(w) + Math.Abs(w1))
                        goto L785;
                    tr = -ra - x * b[(i + 1)][na] + x1 * b[(i + 1)][en];
                    ti = -sa - x * b[(i + 1)][en] - x1 * b[(i + 1)][na];
                    goto L773;

                L785:
                    t1 = (-r - zz * b[(i + 1)][na] + z1 * b[(i + 1)][en]) / y;
                    t2 = (-s - zz * b[(i + 1)][en] - z1 * b[(i + 1)][na]) / y;

                L787:
                    b[i][na] = t1;
                    b[i][en] = t2;

                L790:
                    ;
                }

                // End complex vector
            L795:
                isw = 3 - isw;

            L800:
                ;
            }

            // End back substitution. Transform to original coordinate system.
            for (jj = 0; jj < n; ++jj)
            {
                j = n - jj - 1;

                for (i = 0; i < n; ++i)
                {
                    zz = 0.0;
                    for (k = 0; k <= j; ++k)
                        zz += z[i][k] * b[k][j];
                    z[i][j] = zz;
                }
            }

            // Normalize so that modulus of largest component of each vector is 1.
            // (isw is 1 initially from before)
            for (j = 0; j < n; ++j)
            {
                d = 0.0;
                if (isw == 2) goto L920;
                if (alfi[j] != 0.0) goto L945;

                for (i = 0; i < n; ++i)
                {
                    if ((Math.Abs(z[i][j])) > d)
                        d = (Math.Abs(z[i][j]));
                }

                for (i = 0; i < n; ++i)
                    z[i][j] /= d;

                goto L950;

            L920:
                for (i = 0; i < n; ++i)
                {
                    r = System.Math.Abs(z[i][j - 1]) + System.Math.Abs(z[i][j]);
                    if (r != 0.0)
                    {
                        // Computing 2nd power
                        double u1 = z[i][j - 1] / r;
                        double u2 = z[i][j] / r;
                        r *= Math.Sqrt(u1 * u1 + u2 * u2);
                    }
                    if (r > d)
                        d = r;
                }

                for (i = 0; i < n; ++i)
                {
                    z[i][j - 1] /= d;
                    z[i][j] /= d;
                }

            L945:
                isw = 3 - isw;

            L950:
                ;
            }

            return 0;
        }

        #endregion



        #region ICloneable Members

        private JaggedGeneralizedEigenvalueDecomposition()
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
            var clone = new JaggedGeneralizedEigenvalueDecomposition();
            clone.ai = this.ai.Copy();
            clone.ar = this.ar.Copy();
            clone.beta = this.beta.Copy();
            clone.n = this.n;
            clone.Z = Vector.Copy(this.Z);
            return clone;
        }

        #endregion

    }
}

