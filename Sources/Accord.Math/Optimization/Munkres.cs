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
//
// This class contains the original code written by Robert A. Pilgrim, originally 
// shared under the permissive MIT license. The original license text is reproduced
// below:
//  
//    The MIT License (MIT)
//  
//    Copyright (c) 2000 Robert A. Pilgrim
//                       Murray State University
//                       Dept. of Computer Science & Information Systems
//                       Murray,Kentucky
//  
//    Permission is hereby granted, free of charge, to any person obtaining a copy
//    of this software and associated documentation files (the "Software"), to deal
//    in the Software without restriction, including without limitation the rights
//    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//    copies of the Software, and to permit persons to whom the Software is
//    furnished to do so, subject to the following conditions:
//  
//    The above copyright notice and this permission notice shall be included in all
//    copies or substantial portions of the Software.
//  
//    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//   SOFTWARE.
//  
// This class is also based on the original code written by Yi Cao, originally shared 
// under the permissive BSD license. The original license text is reproduced below:
//
//   Copyright(c) 2009, Yi Cao
//   All rights reserved.
//   
//   Redistribution and use in source and binary forms, with or without
//   modification, are permitted provided that the following conditions are
//   met:
//   
//   * Redistributions of source code must retain the above copyright
//   notice, this list of conditions and the following disclaimer.
//   * Redistributions in binary form must reproduce the above copyright
//   notice, this list of conditions and the following disclaimer in 
//   the documentation and/or other materials provided with the distribution
//   
//   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
//   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
//   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
//   ARE DISCLAIMED.IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
//   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
//   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
//   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
//   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
//   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
//   POSSIBILITY OF SUCH DAMAGE.
//

namespace Openize.Accord.Math.Optimization
{
    using System;
    using System.Collections.Generic;
    using System.Linq;
    using Core;
    using Openize.Accord.Math.Core;
    using Openize.Accord.Math.Matrix;
    using Openize.Accord.Math.Optimization.Base;
    using Elementwise = Openize.Accord.Math.Core.Elementwise;
    using Vector = Openize.Accord.Math.Vector.Vector;

    /// <summary>
    ///   Hungarian method for solving the assignment problem, also known
    ///   as the Kuhn–Munkres algorithm or Munkres assignment algorithm. 
    /// </summary>
    /// 
    ///   
    /// <remarks>
    /// <para>
    ///   The Hungarian method is a combinatorial optimization algorithm that solves the assignment
    ///   problem in polynomial time and which anticipated later primal-dual methods. It was developed
    ///   and published in 1955 by Harold Kuhn, who gave the name "Hungarian method" because the algorithm 
    ///   was largely based on the earlier works of two Hungarian mathematicians: Dénes Kőnig and Jenő
    ///   Egerváry.</para>
    /// <para>
    ///   James Munkres reviewed the algorithm in 1957 and observed that it is (strongly) polynomial. 
    ///   Since then the algorithm has been known also as the Kuhn–Munkres algorithm or Munkres assignment
    ///   algorithm.The time complexity of the original algorithm was O(n^4), however Edmonds and Karp, and 
    ///   independently Tomizawa noticed that it can be modified to achieve an O(n^3) running time. Ford and 
    ///   Fulkerson extended the method to general transportation problems. In 2006, it was discovered that
    ///   Carl Gustav Jacobi had solved the assignment problem in the 19th century, and the solution had been
    ///   published posthumously in 1890 in Latin.</para>
    /// 
    /// <para>
    ///   This code has been based on the original MIT-licensed code by R.A. Pilgrim, available in 
    ///   http://csclab.murraystate.edu/~bob.pilgrim/445/munkres.html, and on the BSD-licensed code
    ///   by Yi Cao, available in http://fr.mathworks.com/matlabcentral/fileexchange/20652-hungarian-algorithm-for-linear-assignment-problems--v2-3-
    /// </para>
    ///   
    /// <para>
    ///   References:
    ///   <list type="bullet">
    ///     <item><description><a href="http://fr.mathworks.com/matlabcentral/fileexchange/20652-hungarian-algorithm-for-linear-assignment-problems--v2-3-">
    ///       Yi Cao (2011). Hungarian Algorithm for Linear Assignment Problems (V2.3). Available in http://fr.mathworks.com/matlabcentral/fileexchange/20652-hungarian-algorithm-for-linear-assignment-problems--v2-3- </a></description></item>
    ///     <item><description><a href="http://csclab.murraystate.edu/~bob.pilgrim/445/munkres.html">
    ///       R. A. Pilgrim (2000). Munkres' Assignment Algorithm Modified for 
    ///       Rectangular Matrices. Available in http://csclab.murraystate.edu/~bob.pilgrim/445/munkres.html </a></description></item>
    ///     <item><description><a href="https://en.wikipedia.org/wiki/Hungarian_algorithm">
    ///       Wikipedia contributors. "Hungarian algorithm." Wikipedia, The Free Encyclopedia.
    ///       Wikipedia, The Free Encyclopedia, 23 Jan. 2016. </a></description></item>
    ///   </list>
    /// </para>   
    /// </remarks>
    /// 
    /// <example>
    /// <code source="Unit Tests\Accord.Tests.Math\Optimization\MunkresTest.cs" region="doc_example" />
    /// </example>
    /// 
    /// <seealso cref="IOptimizationMethod" />
    /// 
    public class Munkres : IOptimizationMethod
    {
        private double[][] originalMatrix;
        private double[][] costMatrix;
        private double[][] validCost;
        private bool[][] stars;
        private bool[] rowCover;
        private bool[] colCover;
        internal int[] starZ;
        internal int[] primeZ;
        private bool[][] validMap;

        private double tolerance = 1e-10;

        private double[] minRow;
        private double[] minCol;

        internal bool[] validRow;
        internal bool[] validCol;

        private int path_row_0;
        private int path_col_0;

        private int nRows;
        private int nCols;
        private int n;

        /// <summary>
        ///   Gets the minimum values across the cost matrix's rows.
        /// </summary>
        /// 
        public double[] MinRow { get { return this.minRow; } }

        /// <summary>
        ///   Gets the minimum values across the cost matrix's columns.
        /// </summary>
        /// 
        public double[] MinCol { get { return this.minCol; } }

        /// <summary>
        ///   Gets a boolean mask indicating which rows contain at least one valid element.
        /// </summary>
        /// 
        public bool[] ValidRow { get { return this.validRow; } }

        /// <summary>
        ///   Gets a boolean mask indicating which columns contain at least one valid element.
        /// </summary>
        /// 
        public bool[] ValidCol { get { return this.validCol; } }

        /// <summary>
        ///   Gets or sets the tolerance value used when performing cost 
        ///   comparisons. Default is 1e-10. If the algorithm takes too
        ///   much time to finish, try decreasing this value.
        /// </summary>
        /// 
        public double Tolerance
        {
            get { return this.tolerance; }
            set { this.tolerance = value; }
        }

        /// <summary>
        ///   Gets or sets the cost matrix for this assignment algorithm. This is
        ///   a (W x T) matrix where N corresponds to the <see cref="NumberOfWorkers"/>
        ///   and T to the <see cref="NumberOfTasks"/>.
        /// </summary>
        /// 
        /// <value>The cost matrix.</value>
        /// 
        public double[][] CostMatrix { get { return this.originalMatrix; } }

        /// <summary>
        /// Gets or sets the number of variables in this optimization problem
        /// (<see cref="NumberOfTasks"/> * <see cref="NumberOfWorkers"/>).
        /// </summary>
        /// 
        public int NumberOfVariables
        {
            get { return this.NumberOfTasks * this.NumberOfWorkers; }
            set { throw new InvalidOperationException(); }
        }

        /// <summary>
        ///   Gets the number of variables (free parameters)
        ///   in the optimization problem. In the assigment
        ///   problem, this gives the number of jobs (or tasks)
        ///   to be performed.
        /// </summary>
        /// 
        /// <value>The number of tasks in the assignment problem.</value>
        /// 
        public int NumberOfTasks
        {
            get { return this.originalMatrix.Columns(); }
        }

        /// <summary>
        ///   Gets or sets the number of workers in the assignment algorithm.
        ///   The workers are the entites that can be assigned jobs according
        ///   to the costs in the <see cref="CostMatrix"/>.
        /// </summary>
        /// 
        /// <value>The number of workers.</value>
        /// 
        public int NumberOfWorkers
        {
            get { return this.originalMatrix.Rows(); }
        }

        /// <summary>
        /// Gets the current solution found, the values of
        /// the parameters which optimizes the function.
        /// </summary>
        /// 
        /// <value>The solution.</value>
        /// 
        public double[] Solution { get; set; }

        /// <summary>
        ///   Gets the output of the function at the current <see cref="Solution" />.
        /// </summary>
        /// 
        /// <value>The value.</value>
        /// 
        public double Value { get; protected set; }


        /// <summary>
        /// Initializes a new instance of the <see cref="Munkres"/> class.
        /// </summary>
        /// 
        /// <param name="numberOfJobs">The number of jobs (tasks) that can be assigned.</param>
        /// <param name="numberOfWorkers">The number of workers that can receive an assignment.</param>
        /// 
        public Munkres(int numberOfJobs, int numberOfWorkers)
        {
            this.init(Jagged.Zeros(numberOfWorkers, numberOfJobs));
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="Munkres"/> class.
        /// </summary>
        /// 
        /// <param name="costMatrix">The cost matrix where each row represents
        ///   a worker, each column represents a task, and each individual element
        ///   represents how much it costs for a particular worker to receive (be
        ///   assigned) a particular task.</param>
        /// 
        public Munkres(double[][] costMatrix)
        {
            this.init(costMatrix);
        }

        private void init(double[][] costMatrix)
        {
            this.originalMatrix = costMatrix;
            this.Solution = new double[this.NumberOfWorkers];
        }

        /// <summary>
        ///   Finds the minimum value of a function. The solution vector
        ///   will be made available at the <see cref="Solution" /> property.
        /// </summary>
        /// 
        /// <returns>Returns <c>true</c> if the method converged to a <see cref="Solution" />.
        /// In this case, the found value will also be available at the <see cref="Value" />
        /// property.</returns>
        /// 
        public bool Minimize()
        {
            return this.run(Matrix.Copy(this.originalMatrix));
        }


        /// <summary>
        ///   Finds the maximum value of a function. The solution vector
        ///   will be made available at the <see cref="Solution" /> property.
        /// </summary>
        /// 
        /// <returns>Returns <c>true</c> if the method converged to a <see cref="Solution" />.
        /// In this case, the found value will also be available at the <see cref="Value" />
        /// property.</returns>
        /// 
        public bool Maximize()
        {
            return this.run(this.originalMatrix.Multiply(-1));
        }

        private bool run(double[][] m)
        {
            this.costMatrix = m;

            int step = 0;

            while (step >= 0)
            {
                step = this.RunStep(step);
            }

            return true;
        }


        internal int RunStep(int step)
        {
            switch (step)
            {
                case 0:
                    return this.step_zero();
                case 1:
                    return this.step_one();
                case 2:
                    return this.step_two();
                case 3:
                    return this.step_three();
                case 4:
                    return this.step_four();
                case 5:
                    return this.step_five();
                case 6:
                    return this.step_six();
                case 7:
                    return this.step_seven();
                default:
                    throw new InvalidOperationException();
            }
        }


        /// <summary>
        ///  Preprocesses the cost matrix to remove infinities and invalid values.
        /// </summary>
        /// 
        /// <returns>Go to step 1.</returns>
        /// 
        private int step_zero()
        {
            this.validRow = new bool[this.NumberOfWorkers];
            this.validCol = new bool[this.NumberOfTasks];
            this.validMap = Jagged.Create<bool>(this.NumberOfWorkers, this.NumberOfTasks);


            double sum = 0;
            double max = Double.NegativeInfinity;

            for (int i = 0; i < this.validRow.Length; i++)
            {
                for (int j = 0; j < this.validCol.Length; j++)
                {
                    double v = Math.Abs(this.costMatrix[i][j]);

                    if (!Double.IsInfinity(v) && !Double.IsNaN(v))
                    {
                        this.validMap[i][j] = this.validRow[i] = this.validCol[j] = true;

                        sum += v;
                        if (v > max)
                            max = v;
                    }
                }
            }

            double bigM = Math.Pow(10, Math.Ceiling(Math.Log10(sum)) + 1);
            for (int i = 0; i < this.costMatrix.Length; i++)
                for (int j = 0; j < this.costMatrix[i].Length; j++)
                    if (!this.validMap[i][j])
                        this.costMatrix[i][j] = Math.Sign(this.costMatrix[i][j]) * bigM;

            this.nRows = Enumerable.Count(this.validRow, x => x);
            this.nCols = Enumerable.Count(this.validCol, x => x);
            this.n = Math.Max(this.nRows, this.nCols);

            if (this.n == 0)
                throw new InvalidOperationException("There are no valid values in the cost matrix.");

            this.validCost = this.costMatrix.Get(this.validRow, this.validCol, result: Jagged.Create(this.n, this.n, 10.0 * max));

            this.rowCover = new bool[this.n];
            this.colCover = new bool[this.n];
            this.stars = Jagged.Create(this.n, this.n, false);

            return 1;
        }

        /// <summary>
        ///  For each row of the cost matrix, find the smallest element
        ///  and subtract it from every element in its row.
        /// </summary>
        /// 
        /// <returns>Go to step 2.</returns>
        /// 
        private int step_one()
        {
            this.minRow = Matrix.Min(this.validCost, dimension: 1);
            this.minCol = Matrix.Min(this.validCost.Subtract(this.MinRow, dimension: (VectorType)1), dimension: 0);

            return 2;
        }

        /// <summary>
        ///   Find a zero (Z) in the resulting matrix. If there is no starred 
        ///   zero in its row or column, star Z. Repeat for each element in the 
        ///   matrix.
        /// </summary>
        /// 
        /// <returns>Go to step 3.</returns>
        /// 
        private int step_two()
        {
            double[][] min;
            this.stars = findZeros(this.validCost, this.minRow, this.minCol, out min);

            this.starZ = Vector.Create(this.n, -1);
            this.primeZ = Vector.Create(this.n, -1);

            for (int j = 0; j < this.stars[0].Length; j++)
            {
                for (int i = 0; i < this.stars.Length; i++)
                {
                    if (this.stars[i][j])
                    {
                        this.starZ[i] = j;
                        this.stars.SetColumn(j, false);
                        this.stars.SetRow(i, false);
                    }
                }
            }

            return 3;
        }

        internal static bool[][] findZeros(double[][] C, double[] minRow, double[] minCol, out double[][] min)
        {
            min = Jagged.CreateAs(C);
            for (int r = 0; r < minRow.Length; r++)
                for (int c = 0; c < minCol.Length; c++)
                    min[r][c] = minRow[r] + minCol[c];

            return Elementwise.Equals(C, min);
        }

        /// <summary>
        ///   Cover each column containing a starred zero. If K columns are covered, 
        ///   the starred zeros describe a complete set of unique assignments. In this 
        ///   case, go to DONE, otherwise, go to Step 4.
        /// </summary>
        /// 
        /// <returns>If K columns are covered, returns 7. Otherwise, returns 4.</returns>
        /// 
        private int step_three()
        {
            this.colCover.Clear();
            this.rowCover.Clear();
            this.primeZ.Clear();

            bool done = true;
            int count = 0;
            for (int i = 0; i < this.starZ.Length; i++)
            {
                int j = this.starZ[i];
                if (j >= 0)
                {
                    this.colCover[j] = true;
                    count++;
                }
                else
                {
                    done = false;
                }
            }

            if (done)
                return 7; // done
            return 4;
        }


        /// <summary>
        ///   Find a noncovered zero and prime it. If there is no starred zero 
        ///   in the row containing this primed zero, Go to Step 5. Otherwise, 
        ///   cover this row and uncover the column containing the starred zero. 
        ///   Continue in this manner until there are no uncovered zeros left. 
        ///   Save the smallest uncovered value and Go to Step 6.        
        /// </summary>
        /// 
        /// <returns>Goes to step 5 or 6.</returns>
        /// 
        private int step_four()
        {
            List<Tuple<int, int>> zeros = this.find_all_zeros();

            while (zeros.Count > 0)
            {
                this.path_row_0 = zeros[0].Item1;
                this.path_col_0 = zeros[0].Item2;
                this.primeZ[this.path_row_0] = this.path_col_0;

                int stz = this.starZ[this.path_row_0];
                if (stz == -1)
                    return 5;

                this.rowCover[this.path_row_0] = true;
                this.colCover[stz] = false;

                zeros.RemoveAll(x => x.Item1 == this.path_row_0);

                // Update
                for (int r = 0; r < Math.Min(this.costMatrix.Length, this.rowCover.Length); r++)
                {
                    if (this.rowCover[r])
                        continue;

                    double a = this.costMatrix[r][stz];
                    double b = this.MinRow[r] + this.MinCol[stz];

                    if (a.IsEqual(b, rtol: this.tolerance))
                        zeros.Add(Tuple.Create(r, stz));
                }
            }

            return 6;
        }

        private List<Tuple<int, int>> find_all_zeros()
        {
            var zeros = new List<Tuple<int, int>>();
            for (int c = 0; c < this.colCover.Length; c++)
            {
                if (this.colCover[c])
                    continue;

                for (int r = 0; r < this.rowCover.Length; r++)
                {
                    if (this.rowCover[r])
                        continue;

                    double a = this.validCost[r][c];
                    double b = this.MinCol[c] + this.MinRow[r];

                    if (a.IsEqual(b, rtol: this.tolerance))
                        zeros.Add(Tuple.Create(r, c));
                }
            }

            return zeros;
        }

        /// <summary>
        ///   Construct a series of alternating primed and starred zeros as follows.  
        ///   Let Z0 represent the uncovered primed zero found in Step 4. Let Z1 denote 
        ///   the starred zero in the column of Z0 (if any). Let Z2 denote the primed zero 
        ///   in the row of Z1 (there will always be one). Continue until the series 
        ///   terminates at a primed zero that has no starred zero in its column.  
        ///   Unstar each starred zero of the series, star each primed zero of the series, 
        ///   erase all primes and uncover every line in the matrix.  
        /// </summary>
        /// 
        /// <returns>Return to Step 3.</returns>
        /// 
        private int step_five()
        {
            int rowZ1 = Array.IndexOf(this.starZ, this.path_col_0);
            this.starZ[this.path_row_0] = this.path_col_0;

            while (rowZ1 >= 0)
            {
                this.starZ[rowZ1] = -1;
                this.path_col_0 = this.primeZ[rowZ1];
                this.path_row_0 = rowZ1;
                rowZ1 = Array.IndexOf(this.starZ, this.path_col_0);
                this.starZ[this.path_row_0] = this.path_col_0;
            }

            return 3;
        }


        /// <summary>
        ///   Add the value found in Step 4 to every element of each covered row, and subtract 
        ///   it from every element of each uncovered column.  
        /// </summary>
        /// 
        /// <returns>Return to step 4.</returns>
        /// 
        private int step_six()
        {
            double minval = this.find_smallest();

            for (int r = 0; r < this.rowCover.Length; r++)
                if (this.rowCover[r])
                    this.minRow[r] -= minval;

            for (int c = 0; c < this.colCover.Length; c++)
                if (!this.colCover[c])
                    this.minCol[c] += minval;

            return 4;
        }

        //methods to support step 6
        private double find_smallest()
        {
            double minval = Double.PositiveInfinity;
            for (int c = 0; c < this.colCover.Length; c++)
            {
                if (this.colCover[c])
                    continue;

                for (int r = 0; r < this.rowCover.Length; r++)
                {
                    if (this.rowCover[r])
                        continue;

                    double v = this.validCost[r][c] - (this.MinRow[r] + this.MinCol[c]);

                    if (v < minval)
                        minval = v;
                }
            }

            return minval;
        }

        private int step_seven()
        {
            // DONE: Assignment pairs are indicated by the positions of the starred zeros in the 
            // cost matrix.If C(i, j) is a starred zero, then the element associated with row i 
            // is assigned to the element associated with column j.
            //
            //                     (http://csclab.murraystate.edu/~bob.pilgrim/445/munkres.html)
            double value = 0;

            for (int i = 0; i < this.nRows; i++)
            {
                this.Solution[i] = Double.NaN;
                int j = this.starZ[i];

                if (j >= 0 && j < this.nCols && this.validMap[i][j])
                {
                    this.Solution[i] = j;
                    value += this.validCost[i][j];
                }
            }

            this.Value = value;
            return -1;
        }


    }
}
