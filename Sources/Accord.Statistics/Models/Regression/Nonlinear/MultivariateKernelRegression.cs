﻿// Accord Statistics Library
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

namespace Openize.Accord.Statistics.Models.Regression.Nonlinear
{
    using System;
    using Accord.MachineLearning.Classifiers;
    using Openize.Accord.Math.Core;
    using Openize.Accord.Math.Matrix;
    using Math.Core;
    using Openize.Accord.Statistics.Analysis;
    using Openize.Accord.Statistics.Kernels.Base;

    /// <summary>
    ///   Multivariate non-linear regression using Kernels.
    /// </summary>
    /// 
    /// <typeparam name="TKernel">The kernel function.</typeparam>
    /// 
    /// <remarks>
    ///   Note: For the moment, the only purpose of this class is to provide a IClassifier
    ///   implementation for models generated by <see cref="KernelDiscriminantAnalysis"/>.
    ///   If you would like to train Support Vector Machines for Regression problems, please
    ///   take a look at the documentation for <c>FanChenLinSupportVectorRegression</c> class.
    /// </remarks>
    /// 
    /// <example>
    /// <code source="Unit Tests\Accord.Tests.Statistics\Analysis\KernelDiscriminantAnalysisTest.cs" region="doc_learn" />
    /// </example>
    /// 
    [Serializable]
    public class MultivariateKernelRegression<TKernel> : MultipleTransformBase<double[], double>
        where TKernel : IKernel<double[]>
    {
        /// <summary>
        ///   Gets or sets the kernel function.
        /// </summary>
        /// 
        public TKernel Kernel { get; set; }

        /// <summary>
        ///   Gets or sets the original input data that is needed to 
        ///   compute the kernel (Gram) matrices for the regression.
        /// </summary>
        /// 
        public double[][] BasisVectors { get; set; }

        /// <summary>
        ///   Gets or sets the linear weights of the regression model. The
        ///   intercept term is not stored in this vector, but is instead
        ///   available through the <see cref="Intercept"/> property.
        /// </summary>
        /// 
        public double[][] Weights { get; set; }

        /// <summary>
        ///   Gets or sets the intercept value for the regression.
        /// </summary>
        /// 
        [Obsolete()]
        public double[] Intercept { get; set; }

        /// <summary>
        ///   Gets or sets the mean values (to be subtracted from samples).
        /// </summary>
        /// 
        public double[] Means { get; set; }

        /// <summary>
        ///   Gets or sets the standard deviations (to be divided from samples).
        /// </summary>
        /// 
        public double[] StandardDeviations { get; set; }

        /// <summary>
        ///   Gets or sets the means of the data in feature space (to center samples).
        /// </summary>
        /// 
        public double[] FeatureMeans { get; set; }

        /// <summary>
        ///   Gets or sets the grand mean of the data in feature space (to center samples).
        /// </summary>
        /// 
        public double FeatureGrandMean { get; set; }

        /// <summary>
        /// Applies the transformation to an input, producing an associated output.
        /// </summary>
        /// <param name="input">The input data to which the transformation should be applied.</param>
        /// <param name="result">The location where the output should be stored.</param>
        /// <returns>
        /// The output generated by applying this transformation to the given input.
        /// </returns>
        public override double[][] Transform(double[][] input, double[][] result)
        {
            if (this.Means != null)
                input = input.Subtract(this.Means, dimension: (VectorType)0);
            if (this.StandardDeviations != null)
                input = input.Divide(this.StandardDeviations, dimension: (VectorType)0);

            // Create the Kernel matrix
            var newK = this.Kernel.ToJagged2(x: input, y: this.BasisVectors);

            if (this.FeatureMeans != null)
                global::Openize.Accord.Statistics.Kernels.Base.Kernel.Center(newK, this.FeatureMeans, this.FeatureGrandMean, result: newK);

            // Project into the kernel principal components
            return Matrix.DotWithTransposed(newK, this.Weights, result: result);
        }
    }


    /// <summary>
    ///   Multivariate non-linear regression using Kernels.
    /// </summary>
    /// 
    [Serializable]
    public class MultivariateKernelRegression : MultivariateKernelRegression<IKernel>
    {

    }
}
