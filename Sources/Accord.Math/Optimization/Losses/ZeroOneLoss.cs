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
    ///   Accuracy loss, also known as zero-one-loss. This class
    ///   provides exactly the same functionality as <see cref="ZeroOneLoss"/>
    ///   but has a more intuitive name. Both classes are interchangeable.
    /// </summary>
    /// 
    [Serializable]
    public class AccuracyLoss : ZeroOneLoss
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="AccuracyLoss"/> class.
        /// </summary>
        /// <param name="expected">The expected outputs (ground truth).</param>
        public AccuracyLoss(double[][] expected)
            : base(expected)
        {
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="AccuracyLoss"/> class.
        /// </summary>
        /// 
        /// <param name="expected">The expected outputs (ground truth).</param>
        /// 
        public AccuracyLoss(double[] expected)
            : base(expected)
        {
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="AccuracyLoss"/> class.
        /// </summary>
        /// 
        /// <param name="expected">The expected outputs (ground truth).</param>
        /// 
        public AccuracyLoss(int[] expected)
            : base(expected)
        {
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="ZeroOneLoss"/> class.
        /// </summary>
        /// 
        /// <param name="expected">The expected outputs (ground truth).</param>
        /// 
        public AccuracyLoss(bool[] expected) 
            : base(expected)
        {
        }
    }

    /// <summary>
    ///   Accuracy loss, also known as zero-one-loss.
    /// </summary>
    /// 
    [Serializable]
    public class ZeroOneLoss : LossBase<int[]>, ILoss<bool[]>,
        ILoss<double[][]>, ILoss<double[]>
    {
        private bool mean = true;

        /// <summary>
        ///   Gets or sets a value indicating whether the average 
        ///   accuracy loss should be computed. Default is true.
        /// </summary>
        /// 
        /// <value>
        ///   <c>true</c> if the average accuracy loss should be computed; otherwise, <c>false</c>.
        /// </value>
        /// 
        public bool Mean
        {
            get { return this.mean; }
            set { this.mean = value; }
        }
        
        /// <summary>
        /// Gets or sets the number of classes.
        /// </summary>
        /// <value>The number of classes.</value>
        public int NumberOfClasses { get; set; }

        /// <summary>
        /// This flag indicates whether the expected class labels are binary.
        /// </summary>
        /// 
        public bool IsBinary
        {
            get { return this.NumberOfClasses == 2; }
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="ZeroOneLoss"/> class.
        /// </summary>
        /// <param name="expected">The expected outputs (ground truth).</param>
        public ZeroOneLoss(double[][] expected) 
            : this(expected.ArgMax(dimension: 0))
        {
        }
        
        /// <summary>
        /// Initializes a new instance of the <see cref="ZeroOneLoss"/> class.
        /// </summary>
        /// <param name="classes">The number of classes.</param>
        /// <param name="expected">The expected outputs (ground truth).</param>
        public ZeroOneLoss(int classes, double[][] expected) 
            : this(classes, expected.ArgMax(dimension: 0))
        {
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="ZeroOneLoss"/> class.
        /// </summary>
        /// 
        /// <param name="expected">The expected outputs (ground truth).</param>
        /// 
        public ZeroOneLoss(double[] expected) 
            : this(expected.ToZeroOne())
        {
        }
        
        /// <summary>
        /// Initializes a new instance of the <see cref="ZeroOneLoss"/> class.
        /// </summary>
        /// <param name="classes">The number of classes.</param>
        /// <param name="expected">The expected outputs (ground truth).</param>
        /// 
        public ZeroOneLoss(int classes, double[] expected) 
            : this(classes, expected.ToZeroOne())
        {
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="ZeroOneLoss"/> class.
        /// </summary>
        /// 
        /// <param name="expected">The expected outputs (ground truth).</param>
        /// 
        public ZeroOneLoss(int[] expected) 
            : this(expected.Max() + 1, expected)
        {
        }
        
        /// <summary>
        /// Initializes a new instance of the <see cref="ZeroOneLoss"/> class.
        /// </summary>
        /// <param name="classes">The number of classes.</param>
        /// <param name="expected">The expected outputs (ground truth).</param>
        /// 
        public ZeroOneLoss(int classes, int[] expected)
        {
            this.NumberOfClasses = classes;
            this.Expected = this.NumberOfClasses == 2 && expected.IsMinusOnePlusOne() ? 
                expected.ToZeroOne() : 
                expected;
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="ZeroOneLoss"/> class.
        /// </summary>
        /// 
        /// <param name="expected">The expected outputs (ground truth).</param>
        /// 
        public ZeroOneLoss(bool[] expected)
            : this(2, expected.ToInt32())
        {
        }

        /// <summary>
        /// Computes the loss between the expected values (ground truth)
        /// and the given actual values that have been predicted.
        /// </summary>
        /// <param name="actual">The actual values that have been predicted.</param>
        /// <returns>
        /// The loss value between the expected values and
        /// the actual predicted values.
        /// </returns>
        public double Loss(double[][] actual)
        {
            return this.Loss(actual.ArgMax(dimension: 0));
        }

        /// <summary>
        /// Computes the loss between the expected values (ground truth)
        /// and the given actual values that have been predicted.
        /// </summary>
        /// <param name="actual">The actual values that have been predicted.</param>
        /// <returns>
        /// The loss value between the expected values and
        /// the actual predicted values.
        /// </returns>
        public override double Loss(int[] actual)
        {
            this.NumberOfClasses = Math.Max(this.NumberOfClasses, actual.Max() + 1);

            if (this.NumberOfClasses == 2)
            {
                actual = actual.ToZeroOne();
            }
            
            int error = 0;
            for (int i = 0; i < this.Expected.Length; i++)
                if (this.Expected[i] != actual[i])
                    error++;

            if (this.mean)
                return error / (double)this.Expected.Length;
            return error;
        }

        /// <summary>
        /// Computes the loss between the expected values (ground truth)
        /// and the given actual values that have been predicted.
        /// </summary>
        /// <param name="actual">The actual values that have been predicted.</param>
        /// <returns>
        /// The loss value between the expected values and
        /// the actual predicted values.
        /// </returns>
        public double Loss(bool[] actual)
        {
            return this.Loss(actual.ToZeroOne());
        }

        /// <summary>
        /// Computes the loss between the expected values (ground truth)
        /// and the given actual values that have been predicted.
        /// </summary>
        /// <param name="actual">The actual values that have been predicted.</param>
        /// <returns>
        /// The loss value between the expected values and
        /// the actual predicted values.
        /// </returns>
        public double Loss(double[] actual)
        {
            return this.Loss(Classes.Decide(actual));
        }
    }
}
