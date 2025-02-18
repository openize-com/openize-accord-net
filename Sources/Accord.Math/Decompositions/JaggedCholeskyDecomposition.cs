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
    using System.Threading.Tasks;
    using Base;
    using Openize.Accord.Math.Matrix;
    using Openize.Accord.Math.Vector;
    using Openize.Accord.Core.Exceptions;
    using Vector = Openize.Accord.Math.Vector.Vector;

    /// <summary>
    ///   Cholesky Decomposition of a symmetric, positive definite matrix.
    /// </summary>
    ///
    /// <remarks>
    ///   <para>
    ///     For a symmetric, positive definite matrix <c>A</c>, the Cholesky decomposition is a
    ///     lower triangular matrix <c>L</c> so that <c>A = L * L'</c>.
    ///     If the matrix is not positive definite, the constructor returns a partial 
    ///     decomposition and sets two internal variables that can be queried using the
    ///     <see cref="IsUndefined"/> and <see cref="IsPositiveDefinite"/> properties.</para>
    ///   <para>
    ///     Any square matrix A with non-zero pivots can be written as the product of a
    ///     lower triangular matrix L and an upper triangular matrix U; this is called
    ///     the LU decomposition. However, if A is symmetric and positive definite, we
    ///     can choose the factors such that U is the transpose of L, and this is called
    ///     the Cholesky decomposition. Both the LU and the Cholesky decomposition are
    ///     used to solve systems of linear equations.</para>
    ///   <para>
    ///     When it is applicable, the Cholesky decomposition is twice as efficient
    ///     as the LU decomposition.</para>
    ///    </remarks>
    ///    
    [Serializable]
    public sealed class JaggedCholeskyDecomposition : ICloneable, ISolverArrayDecomposition<Double>
    {

        private Double[][] L;
        private Double[] D;
        private int n;

        private bool positiveDefinite;
        private bool undefined;
        private bool robust;
        private bool destroyed;

        // cache for lazy evaluation
        private Double[][] leftTriangularFactor;
        private Double[][] diagonalMatrix;
        private Double? determinant;
        private double? lndeterminant;
        private bool? nonsingular;


        /// <summary>
        ///   Constructs a new Cholesky Decomposition.
        /// </summary>
        /// 
        /// <param name="value">
        ///   The symmetric matrix, given in upper triangular form, to be decomposed.</param>
        /// <param name="robust">
        ///   True to perform a square-root free LDLt decomposition, false otherwise.</param>
        /// <param name="inPlace">
        ///   True to perform the decomposition in place, storing the factorization in the
        ///   lower triangular part of the given matrix.</param>
        /// <param name="valueType">
        ///   How to interpret the matrix given to be decomposed. Using this parameter, a lower or
        ///   upper-triangular matrix can be interpreted as a symmetric matrix by assuming both lower
        ///   and upper parts contain the same elements. Use this parameter in conjunction with inPlace
        ///   to save memory by storing the original matrix and its decomposition at the same memory
        ///   location (lower part will contain the decomposition's L matrix, upper part will contains 
        ///   the original matrix).</param>
        /// 
        public JaggedCholeskyDecomposition(Double[][] value, bool robust = false, 
            bool inPlace = false, MatrixType valueType = MatrixType.UpperTriangular)
        {
            if (value.Rows() != value.Columns())
                throw new DimensionMismatchException("value", "Matrix is not square.");

            if (!inPlace)
                value = value.Copy();

            this.n = value.Rows();
            this.L = value.ToUpperTriangular(valueType, result: value);
            this.robust = robust;

            if (robust)
            {
                this.LDLt(); // Compute square-root free decomposition
            }
            else
            {
                this.LLt(); // Compute standard Cholesky decomposition
            }
        }

        /// <summary>
        ///   Gets whether the decomposed matrix was positive definite.
        /// </summary>
        ///
        public bool IsPositiveDefinite
        {
            get { return this.positiveDefinite && !this.undefined; }
        }

        /// <summary>
        ///   Gets a value indicating whether the LDLt factorization
        ///   has been computed successfully or if it is undefined.
        /// </summary>
        /// 
        /// <value>
        ///     <c>true</c> if the factorization is not defined; otherwise, <c>false</c>.
        /// </value>
        /// 
        public bool IsUndefined
        {
            get { return this.undefined; }
        }

        /// <summary>
        ///   Gets the left (lower) triangular factor
        ///   <c>L</c> so that <c>A = L * D * L'</c>.
        /// </summary>
        /// 
        public Double[][] LeftTriangularFactor
        {
            get
            {
                if (this.leftTriangularFactor == null)
                {
                    if (this.destroyed)
                        throw new InvalidOperationException("The decomposition has been destroyed.");
                        
                    if (this.undefined)
                        throw new InvalidOperationException("The decomposition is undefined (zero in diagonal).");

                    this.leftTriangularFactor = this.L.GetLowerTriangle();
                }

                return this.leftTriangularFactor;
            }
        }

        /// <summary>
        ///   Gets the block diagonal matrix of diagonal elements in a LDLt decomposition.
        /// </summary>        
        ///   
        public Double[][] DiagonalMatrix
        {
            get
            {
                if (this.diagonalMatrix == null)
                {
                    if (this.destroyed)
                        throw new InvalidOperationException("The decomposition has been destroyed.");
                        
                    this.diagonalMatrix = Jagged.Diagonal(this.D);
                }

                return this.diagonalMatrix;
            }
        }

        /// <summary>
        ///   Gets the one-dimensional array of diagonal elements in a LDLt decomposition.
        /// </summary>        
        /// 
        public Double[] Diagonal
        {
            get { return this.D; }
        }

        /// <summary>
        ///   Gets the determinant of the decomposed matrix.
        /// </summary>
        /// 
        public Double Determinant
        {
            get
            {
                if (!this.determinant.HasValue)
                {
                    if (this.destroyed)
                        throw new InvalidOperationException("The decomposition has been destroyed.");

                    if (this.undefined)
                        throw new InvalidOperationException("The decomposition is undefined (zero in diagonal).");

                    Double detL = 1, detD = 1;
                    for (int i = 0; i < this.n; i++)
                        detL *= this.L[i][i];

                    if (this.D != null)
                    {
                        for (int i = 0; i < this.n; i++)
                            detD *= this.D[i];
                    }

                    this.determinant = detL * detL * detD;
                }

                return this.determinant.Value;
            }
        }

        /// <summary>
        ///   If the matrix is positive-definite, gets the
        ///   log-determinant of the decomposed matrix.
        /// </summary>
        /// 
        public double LogDeterminant
        {
            get
            {
                if (!this.lndeterminant.HasValue)
                {
                    if (this.destroyed)
                        throw new InvalidOperationException("The decomposition has been destroyed.");
                    
                    if (this.undefined)
                        throw new InvalidOperationException("The decomposition is undefined (zero in diagonal).");

                    double detL = 0, detD = 0;
                    for (int i = 0; i < this.n; i++)
                        detL += Math.Log((double)this.L[i][i]);

                    if (this.D != null)
                    {
                        for (int i = 0; i < this.D.Length; i++)
                            detD += Math.Log((double)this.D[i]);
                    }

                    this.lndeterminant = detL + detL + detD;
                }

                return this.lndeterminant.Value;
            }
        }

        /// <summary>
        ///   Gets a value indicating whether the decomposed
        ///   matrix is non-singular (i.e. invertible).
        /// </summary>
        /// 
        public bool Nonsingular
        {
            get
            {
                if (!this.nonsingular.HasValue)
                {
                    if (this.destroyed)
                        throw new InvalidOperationException("The decomposition has been destroyed.");
                        
                    bool nonSingular = true;
                    for (int i = 0; i < this.n && nonSingular; i++)
                        if (this.L[i][i] == 0 || this.D[i] == 0) nonSingular = false;

                    this.nonsingular = nonSingular;
                }

                return this.nonsingular.Value;
            }
        }


        private unsafe void LLt()
        {
            this.D = Vector.Ones<Double>(this.n);

            this.positiveDefinite = true;

            for (int j = 0; j < this.n; j++)
            {
                Double s = 0;
                for (int k = 0; k < j; k++)
                {
                    Double t = this.L[k][j];
                    for (int i = 0; i < k; i++)
                        t -= this.L[j][i] * this.L[k][i];
                    t = t / this.L[k][k];

                    this.L[j][k] = t;
                    s += t * t;
                }

                s = this.L[j][j] - s;

                // Use a tolerance for positive-definiteness
                this.positiveDefinite &= (s > (Double)1e-14 * Math.Abs(this.L[j][j]));

                this.L[j][j] = (Double)Math.Sqrt((double)s);
            }
        }


        private unsafe void LDLt()
        {
            this.D = new Double[this.n];

            this.positiveDefinite = true;

            Double[] v = new Double[this.n];
            for (int i = 0; i < this.L.Length; i++)
            {
                for (int j = 0; j < i; j++)
                    v[j] = this.L[i][j] * this.D[j];

                Double d = 0;
                for (int k = 0; k < i; k++)
                    d += this.L[i][k] * v[k];

                d = this.D[i] = v[i] = this.L[i][i] - d;

                // Use a tolerance for positive-definiteness
                this.positiveDefinite &= (v[i] > (Double)1e-14 * Math.Abs(this.L[i][i]));

                // If one of the diagonal elements is zero, the 
                // decomposition (without pivoting) is undefined.
                if (v[i] == 0) { this.undefined = true; return; }
                
                Parallel.For(i + 1, this.L.Length, k =>
                {
                     Double s = 0;
                     for (int j = 0; j < i; j++)
                         s += this.L[k][j] * v[j];

                     this.L[k][i] = (this.L[i][k] - s) / d;
                });
            }

            for (int i = 0; i < this.L.Length; i++)
                this.L[i][i] = 1;
        }

        /// <summary>
        ///   Solves a set of equation systems of type <c>A * X = B</c>.
        /// </summary>
        /// 
        /// <param name="value">Right hand side matrix with as many rows as <c>A</c> and any number of columns.</param>
        /// <returns>Matrix <c>X</c> so that <c>L * L' * X = B</c>.</returns>
        /// <exception cref="T:System.ArgumentException">Matrix dimensions do not match.</exception>
        /// <exception cref="T:System.NonSymmetricMatrixException">Matrix is not symmetric.</exception>
        /// <exception cref="T:System.NonPositiveDefiniteMatrixException">Matrix is not positive-definite.</exception>
        /// 
        public Double[][] Solve(Double[][] value)
        {
            return this.Solve(value, false);
        }

        /// <summary>
        ///   Solves a set of equation systems of type <c>A * X = B</c>.
        /// </summary>
        /// 
        /// <param name="value">Right hand side matrix with as many rows as <c>A</c> and any number of columns.</param>
        /// <returns>Matrix <c>X</c> so that <c>L * L' * X = B</c>.</returns>
        /// <exception cref="T:System.ArgumentException">Matrix dimensions do not match.</exception>
        /// <exception cref="T:System.NonSymmetricMatrixException">Matrix is not symmetric.</exception>
        /// <exception cref="T:System.NonPositiveDefiniteMatrixException">Matrix is not positive-definite.</exception>
        /// <param name="inPlace">True to compute the solving in place, false otherwise.</param>
        /// 
        public Double[][] Solve(Double[][] value, bool inPlace)
        {
            if (value == null)
                throw new ArgumentNullException("value");

            if (value.Length != this.n)
                throw new ArgumentException("Argument matrix should have the same number of rows as the decomposed matrix.", "value");

            if (!this.robust && !this.positiveDefinite)
                throw new NonPositiveDefiniteMatrixException("Decomposed matrix is not positive definite.");

            if (this.destroyed)
                throw new InvalidOperationException("The decomposition has been destroyed.");

            if (this.undefined)
                throw new InvalidOperationException("The decomposition is undefined (zero in diagonal).");

            int count = value[0].Length;
            var B = inPlace ? value : value.MemberwiseClone();

            // Solve L*Y = B;
            for (int k = 0; k < this.n; k++)
            {
                for (int j = 0; j < B[k].Length; j++)
                {
                    for (int i = 0; i < k; i++)
                        B[k][j] -= B[i][j] * this.L[k][i];
                    B[k][j] /= this.L[k][k];
                }
            }

            if (this.robust)
            {
                for (int k = 0; k < this.D.Length; k++)
                    for (int j = 0; j < B[k].Length; j++)
                        B[k][j] /= this.D[k];
            }

            // Solve L'*X = Y;
            for (int k = this.n - 1; k >= 0; k--)
            {
                for (int j = 0; j < B[k].Length; j++)
                {
                    for (int i = k + 1; i < this.n; i++)
                        B[k][j] -= B[i][j] * this.L[i][k];

                    B[k][j] /= this.L[k][k];
                }
            }

            return B;
        }

        /// <summary>
        ///   Solves a set of equation systems of type <c>A * x = b</c>.
        /// </summary>
        /// 
        /// <param name="value">Right hand side column vector with as many rows as <c>A</c>.</param>
        /// <returns>Vector <c>x</c> so that <c>L * L' * x = b</c>.</returns>
        /// <exception cref="T:System.ArgumentException">Matrix dimensions do not match.</exception>
        /// <exception cref="T:System.NonSymmetricMatrixException">Matrix is not symmetric.</exception>
        /// <exception cref="T:System.NonPositiveDefiniteMatrixException">Matrix is not positive-definite.</exception>
        /// 
        public Double[] Solve(Double[] value)
        {
            return this.Solve(value, false);
        }

        /// <summary>
        ///   Solves a set of equation systems of type <c>A * x = b</c>.
        /// </summary>
        /// 
        /// <param name="value">Right hand side column vector with as many rows as <c>A</c>.</param>
        /// <returns>Vector <c>x</c> so that <c>L * L' * x = b</c>.</returns>
        /// <exception cref="T:System.ArgumentException">Matrix dimensions do not match.</exception>
        /// <exception cref="T:System.NonSymmetricMatrixException">Matrix is not symmetric.</exception>
        /// <exception cref="T:System.NonPositiveDefiniteMatrixException">Matrix is not positive-definite.</exception>
        /// <param name="inPlace">True to compute the solving in place, false otherwise.</param>
        /// 
        public Double[] Solve(Double[] value, bool inPlace)
        {
            if (value == null)
                throw new ArgumentNullException("value");

            if (value.Length != this.n)
                throw new ArgumentException("Argument vector should have the same length as rows in the decomposed matrix.", "value");

            if (!this.robust && !this.positiveDefinite)
                throw new NonPositiveDefiniteMatrixException("Decomposed matrix is not positive definite.");
                
            if (this.destroyed)
                throw new InvalidOperationException("The decomposition has been destroyed.");

            if (this.undefined)
                throw new InvalidOperationException("The decomposition is undefined (zero in diagonal).");

            var B = inPlace ? value : value.Copy();

            // Solve L*Y = B;
            for (int k = 0; k < this.n; k++)
            {
                for (int i = 0; i < k; i++)
                    B[k] -= B[i] * this.L[k][i];
                B[k] /= this.L[k][k];
            }

            if (this.robust)
            {
                for (int k = 0; k < this.D.Length; k++)
                    B[k] /= this.D[k];
            }

            // Solve L'*X = Y;
            for (int k = this.n - 1; k >= 0; k--)
            {
                for (int i = k + 1; i < this.n; i++)
                    B[k] -= B[i] * this.L[i][k];
                B[k] /= this.L[k][k];
            }
            
            return B;
        }

		/// <summary>
        ///   Solves a set of equation systems of type <c>A * X = B</c> where B is a diagonal matrix.
        /// </summary>
        /// <param name="diagonal">Diagonal fo the right hand side matrix with as many rows as <c>A</c>.</param>
        /// <returns>Matrix <c>X</c> so that <c>L * U * X = B</c>.</returns>
        /// 
        public Double[][] SolveForDiagonal(Double[] diagonal)
        {
            if (diagonal == null)
                throw new ArgumentNullException("diagonal");

            return this.Solve(Jagged.Diagonal(diagonal));
        }

        /// <summary>
        ///   Solves a set of equation systems of type <c>A * X = I</c>.
        /// </summary>
        /// 
        public Double[][] Inverse()
        {
            return this.Solve(Jagged.Identity<Double>(this.n));
        }

        /// <summary>
        ///   Computes the diagonal of the inverse of the decomposed matrix.
        /// </summary>
        /// 
        public Double[] InverseDiagonal(bool destroy = false)
        {
            return this.InverseDiagonal(new Double[this.n], destroy);
        }
        
        /// <summary>
        ///   Computes the diagonal of the inverse of the decomposed matrix.
        /// </summary>
        ///
        /// <param name="destroy">True to conserve memory by reusing the
        ///    same space used to hold the decomposition, thus destroying
        ///    it in the process. Pass false otherwise.</param>
        /// <param name="result">The array to hold the result of the 
        ///    computation. Should be of same length as the the diagonal
        ///    of the original matrix.</param>
        /// 
        public Double[] InverseDiagonal(Double[] result, bool destroy = false)
        {
            if (!this.robust && !this.positiveDefinite)
                throw new NonPositiveDefiniteMatrixException("Matrix is not positive definite.");

            if (this.destroyed)
                throw new InvalidOperationException("The decomposition has been destroyed.");

            if (this.undefined)
                throw new InvalidOperationException("The decomposition is undefined (zero in diagonal).");

            Double[][] S;

            if (destroy)
            {
                S = this.L; this.destroyed = true;
            }
            else
            {
                S = Jagged.Zeros<Double>(this.n, this.n);
            }

            // References:
            // http://books.google.com/books?id=myzIPBwyBbcC&pg=PA119

            // Compute the inverse S of the lower triangular matrix L
            // and store in place of the upper triangular part of S.

            for (int j = this.n - 1; j >= 0; j--)
            {
                S[j][j] = 1 / this.L[j][j];
                for (int i = j - 1; i >= 0; i--)
                {
                    Double sum = 0;
                    for (int k = i + 1; k <= j; k++)
                        sum += this.L[k][i] * S[k][j];
                    S[i][j] = -sum / this.L[i][i];
                }
            }

            // Compute the 2-norm squared of the rows
            // of the upper (right) triangular matrix S.
            if (this.robust)
            {
                for (int i = 0; i < S.Length; i++)
                {
                    Double sum = 0;
                    for (int j = i; j < S[i].Length; j++)
                        sum += S[i][j] * S[i][j] / this.D[j];
                    result[i] = sum;
                }
            }
            else
            {
                for (int i = 0; i < S.Length; i++)
                {
                    Double sum = 0;
                    for (int j = i; j < S[i].Length; j++)
                        sum += S[i][j] * S[i][j];
                    result[i] = sum;
                }
            }

            return result;
        }

        /// <summary>
        ///   Computes the trace of the inverse of the decomposed matrix.
        /// </summary>
        ///
        /// <param name="destroy">True to conserve memory by reusing the
        ///    same space used to hold the decomposition, thus destroying
        ///    it in the process. Pass false otherwise.</param>
        /// 
        public Double InverseTrace(bool destroy = false)
        {
            if (!this.robust && !this.positiveDefinite)
                throw new NonPositiveDefiniteMatrixException("Matrix is not positive definite.");

            if (this.destroyed)
                throw new InvalidOperationException("The decomposition has been destroyed.");
                
            if (this.undefined)
                throw new InvalidOperationException("The decomposition is undefined (zero in diagonal).");

            Double[][] S;

            if (destroy)
            {
                S = this.L; this.destroyed = true;
            }
            else
            {
                S = Jagged.Zeros<Double>(this.n, this.n);
            }

            // References:
            // http://books.google.com/books?id=myzIPBwyBbcC&pg=PA119

            // Compute the inverse S of the lower triangular matrix L
            // and store in place of the upper triangular part of S.

            for (int j = this.n - 1; j >= 0; j--)
            {
                S[j][j] = 1 / this.L[j][j];
                for (int i = j - 1; i >= 0; i--)
                {
                    Double sum = 0;
                    for (int k = i + 1; k <= j; k++)
                        sum += this.L[k][i] * S[k][j];
                    S[i][j] = -sum / this.L[i][i];
                }
            }

            // Compute the 2-norm squared of the rows
            // of the upper (right) triangular matrix S.
            Double trace = 0;
            
            if (this.robust)
            {
                for (int i = 0; i < S.Length; i++)
                    for (int j = i; j < S[i].Length; j++)
                        trace += S[i][j] * S[i][j] / this.D[j];
            }
            else
            {
                for (int i = 0; i < S.Length; i++)
                    for (int j = i; j < S[i].Length; j++)
                        trace += S[i][j] * S[i][j];
            }
            
            return trace;
        }

        /// <summary>
        ///   Reverses the decomposition, reconstructing the original matrix <c>X</c>.
        /// </summary>
        /// 
        public Double[][] Reverse()
        {
            if (this.destroyed)
                throw new InvalidOperationException("The decomposition has been destroyed.");
                
            if (this.undefined)
                throw new InvalidOperationException("The decomposition is undefined (zero in diagonal).");

            if (this.robust)
                return this.LeftTriangularFactor.Dot(this.DiagonalMatrix).DotWithTransposed(this.LeftTriangularFactor);
            return this.LeftTriangularFactor.DotWithTransposed(this.LeftTriangularFactor);
        }

        /// <summary>
        ///   Computes <c>(Xt * X)^1</c> (the inverse of the covariance matrix). This
        ///   matrix can be used to determine standard errors for the coefficients when
        ///   solving a linear set of equations through any of the <see cref="Solve(Double[][])"/>
        ///   methods.
        /// </summary>
        /// 
        public Double[][] GetInformationMatrix()
        {
            var X = this.Reverse();
            return X.TransposeAndDot(X).Inverse();
        }

        /// <summary>
        ///   Creates a new Cholesky decomposition directly from
        ///   an already computed left triangular matrix <c>L</c>.
        /// </summary>
        /// <param name="leftTriangular">The left triangular matrix from a Cholesky decomposition.</param>
        /// 
        public static JaggedCholeskyDecomposition FromLeftTriangularMatrix(Double[][] leftTriangular)
        {
            var chol = new JaggedCholeskyDecomposition();
            chol.n = leftTriangular.Length;
            chol.L = leftTriangular;
            chol.positiveDefinite = true;
            chol.robust = false;
            chol.D = Vector.Ones<Double>(chol.n);
            return chol;
        }


        #region ICloneable Members

        private JaggedCholeskyDecomposition()
        {
        }

        /// <summary>
        ///   Creates a new object that is a copy of the current instance.
        /// </summary>
        /// <returns>
        ///   A new object that is a copy of this instance.
        /// </returns>
        /// 
        public object Clone()
        {
            var clone = new JaggedCholeskyDecomposition();
            clone.L = this.L.MemberwiseClone();
            clone.D = (Double[])this.D.Clone();
            clone.destroyed = this.destroyed;
            clone.n = this.n;
            clone.undefined = this.undefined;
            clone.robust = this.robust;
            clone.positiveDefinite = this.positiveDefinite;
            return clone;
        }

        #endregion

    }
}

