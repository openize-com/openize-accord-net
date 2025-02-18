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
    using Openize.Accord.Math.Core;
    using Openize.Accord.Math.Matrix;
    using Openize.Accord.Core.Exceptions;
    using Openize.Accord.Math;

    /// <summary>
    ///   Gram-Schmidt Orthogonalization.
    /// </summary>
    /// 
    public class GramSchmidtOrthogonalization
    {
        private double[,] q;
        private double[,] r;


        /// <summary>
        ///   Initializes a new instance of the <see cref="GramSchmidtOrthogonalization"/> class.
        /// </summary>
        /// 
        /// <param name="value">The matrix <c>A</c> to be decomposed.</param>
        /// 
        public GramSchmidtOrthogonalization(double[,] value)
            : this(value, true)
        {
        }

        /// <summary>
        ///   Initializes a new instance of the <see cref="GramSchmidtOrthogonalization"/> class.
        /// </summary>
        /// 
        /// <param name="value">The matrix <c>A</c> to be decomposed.</param>
        /// <param name="modified">True to use modified Gram-Schmidt; false
        ///   otherwise. Default is true (and is the recommended setup).</param>
        /// 
        public GramSchmidtOrthogonalization(double[,] value, bool modified)
        {
            if (value.GetLength(0) != value.GetLength(1))
                throw new DimensionMismatchException("value", "Matrix must be square.");

            int size = value.GetLength(0);

            this.q = new double[size, size];
            this.r = new double[size, size];

            if (modified)
            {
                for (int j = 0; j < size; j++)
                {
                    double[] v = value.GetColumn(j);

                    for (int i = 0; i < j; i++)
                    {
                        this.r[i, j] = this.q.GetColumn(i).Dot(v);
                        var t = this.r[i, j].Multiply(this.q.GetColumn(i));
                        v.Subtract(t, result: v);
                    }

                    this.r[j, j] = Norm.Euclidean(v);
                    this.q.SetColumn(j, v.Divide(this.r[j, j]));
                }
            }

            else
            {
                for (int j = 0; j < size; j++)
                {
                    double[] v = value.GetColumn(j);
                    double[] a = value.GetColumn(j);

                    for (int i = 0; i < j; i++)
                    {
                        this.r[i, j] = this.q.GetColumn(j).Dot(a);
                        v = v.Subtract(this.r[i, j].Multiply(this.q.GetColumn(i)));
                    }

                    this.r[j, j] = Norm.Euclidean(v);
                    this.q.SetColumn(j, v.Divide(this.r[j, j]));
                }
            }
        }


        /// <summary>
        ///   Returns the orthogonal factor matrix <c>Q</c>.
        /// </summary>
        /// 
        public double[,] OrthogonalFactor
        {
            get { return this.q; }
        }

        /// <summary>
        ///   Returns the upper triangular factor matrix <c>R</c>.
        /// </summary>
        /// 
        public double[,] UpperTriangularFactor
        {
            get { return this.r; }
        }
    }
}
