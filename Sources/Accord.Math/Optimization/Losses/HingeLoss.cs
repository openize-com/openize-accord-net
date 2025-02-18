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

namespace Openize.Accord.Math.Optimization.Losses
{
    using System;
    using Openize.Accord.Math.Matrix;
    using Openize.Accord.Math.Accord.Statistics;

    /// <summary>
    ///   Hinge loss.
    /// </summary>
    /// 
    [Serializable]
    public struct HingeLoss : ILoss<double[]>,
        IDifferentiableLoss<bool, double, double>,
        IDifferentiableLoss<double, double, double>
    {

        private bool[][] expected;

        /// <summary>
        ///   Gets the expected outputs (the ground truth).
        /// </summary>
        /// 
        public bool[][] Expected
        {
            get { return this.expected; }
            set { this.expected = value; }
        }

        /// <summary>
        ///   Initializes a new instance of the <see cref="HingeLoss"/> class.
        /// </summary>
        /// 
        /// <param name="expected">The expected outputs (ground truth).</param>
        /// 
        public HingeLoss(double[][] expected)
        {
            this.expected = Classes.Decide(expected);
        }

        /// <summary>
        ///   Initializes a new instance of the <see cref="HingeLoss"/> class.
        /// </summary>
        /// 
        /// <param name="expected">The expected outputs (ground truth).</param>
        /// 
        public HingeLoss(double[] expected)
        {
            this.expected = Classes.Decide(Jagged.ColumnVector(expected));
        }

        /// <summary>
        ///   Initializes a new instance of the <see cref="HingeLoss"/> class.
        /// </summary>
        /// 
        /// <param name="expected">The expected outputs (ground truth).</param>
        /// 

        public HingeLoss(int[] expected)
        {
            if (Classes.IsMinusOnePlusOne(expected))
                expected = expected.ToZeroOne();

            this.expected = Jagged.OneHot<bool>(expected);
        }

        /// <summary>
        ///   Initializes a new instance of the <see cref="HingeLoss"/> class.
        /// </summary>
        /// 
        /// <param name="expected">The expected outputs (ground truth).</param>
        /// 
        public HingeLoss(bool[] expected)
        {
            this.expected = Jagged.ColumnVector(expected);
        }

        /// <summary>
        ///   Computes the loss between the expected values (ground truth)
        ///   and the given actual values that have been predicted.
        /// </summary>
        /// 
        /// <param name="actual">The actual values that have been predicted.</param>
        /// 
        /// <returns>
        ///   The loss value between the expected values and
        ///   the actual predicted values.
        /// </returns>
        /// 
        public double Loss(double[][] actual)
        {
            double error = 0;
            for (int i = 0; i < this.Expected.Length; i++)
                for (int j = 0; j < this.Expected[i].Length; j++)
                    error += this.Loss(this.Expected[i][j], actual[i][j]);
            return error;
        }

        /// <summary>
        ///   Computes the loss between the expected values (ground truth)
        ///   and the given actual values that have been predicted.
        /// </summary>
        /// 
        /// <param name="actual">The actual values that have been predicted.</param>
        /// 
        /// <returns>
        ///   The loss value between the expected values and
        ///   the actual predicted values.
        /// </returns>
        /// 
        public double Loss(double[] actual)
        {
            double error = 0;
            for (int i = 0; i < this.Expected.Length; i++)
                error += this.Loss(this.Expected[i][0], actual[i]);
            return error;
        }


        /// <summary>
        ///   Computes the derivative of the loss between the expected values (ground truth) 
        ///   and the given actual values that have been predicted.
        /// </summary>
        /// 
        /// <param name="actual">The actual values that have been predicted.</param>
        /// <param name="expected">The expected values that should have been predicted.</param>
        /// 
        /// <returns>The loss value between the expected values and
        ///   the actual predicted values.</returns>
        /// 
        public double Loss(bool expected, double actual)
        {
            if (expected)
            {
                if (actual > 1)
                    return 0;
                return 1 - actual;
            }
            else
            {
                if (actual < -1)
                    return 0;
                return 1 + actual;
            }
        }

        /// <summary>
        ///   Computes the derivative of the loss between the expected values (ground truth) 
        ///   and the given actual values that have been predicted.
        /// </summary>
        /// 
        /// <param name="actual">The actual values that have been predicted.</param>
        /// <param name="expected">The expected values that should have been predicted.</param>
        /// 
        /// <returns>The loss value between the expected values and
        ///   the actual predicted values.</returns>
        /// 
        public double Derivative(bool expected, double actual)
        {
            if (expected)
            {
                if (actual > 1)
                    return 0;
                return actual;
            }
            else
            {
                if (actual < -1)
                    return 0;
                return actual;
            }
        }

        /// <summary>
        ///   Computes the derivative of the loss between the expected values (ground truth) 
        ///   and the given actual values that have been predicted.
        /// </summary>
        /// 
        /// <param name="actual">The actual values that have been predicted.</param>
        /// <param name="expected">The expected values that should have been predicted.</param>
        /// 
        /// <returns>The loss value between the expected values and
        ///   the actual predicted values.</returns>
        /// 
        public double Loss(double expected, double actual)
        {
            return this.Loss(Classes.Decide(expected), actual);
        }

        /// <summary>
        ///   Computes the derivative of the loss between the expected values (ground truth) 
        ///   and the given actual values that have been predicted.
        /// </summary>
        /// 
        /// <param name="actual">The actual values that have been predicted.</param>
        /// <param name="expected">The expected values that should have been predicted.</param>
        /// 
        /// <returns>The loss value between the expected values and
        ///   the actual predicted values.</returns>
        /// 
        public double Derivative(double expected, double actual)
        {
            return this.Derivative(Classes.Decide(expected), actual);
        }

    }
}
