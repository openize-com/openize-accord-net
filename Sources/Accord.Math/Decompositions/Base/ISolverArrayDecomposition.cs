﻿// Accord Math Library
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

namespace Openize.Accord.Math.Decompositions.Base
{

    /// <summary>
    ///   Common interface for matrix decompositions which
    ///   can be used to solve linear systems of equations
    ///   involving jagged array matrices.
    /// </summary>
    /// 
    public interface ISolverArrayDecomposition<T> where T : struct
    {

        /// <summary>
        ///   Solves a set of equation systems of type <c>A * X = B</c>.
        /// </summary>
        /// 
        T[][] Solve(T[][] value);

        /// <summary>
        ///   Solves a set of equation systems of type <c>A * X = B</c>.
        /// </summary>
        /// 
        T[] Solve(T[] value);

        /// <summary>
        ///   Solves a set of equation systems of type <c>A * X = B</c> where B is a diagonal matrix.
        /// </summary>
        /// 
        T[][] SolveForDiagonal(T[] diagonal);

        /// <summary>
        ///   Solves a set of equation systems of type <c>A * X = I</c>.
        /// </summary>
        /// 
        T[][] Inverse();

        /// <summary>
        ///   Computes <c>(Xt * X)^1</c> (the inverse of the covariance matrix). This
        ///   matrix can be used to determine standard errors for the coefficients when
        ///   solving a linear set of equations through any of the <see cref="Solve(T[][])"/>
        ///   methods.
        /// </summary>
        /// 
        T[][] GetInformationMatrix();

        /// <summary>
        ///   Reverses the decomposition, reconstructing the original matrix <c>X</c>.
        /// </summary>
        /// 
        T[][] Reverse();
    }
}
