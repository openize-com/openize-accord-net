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

namespace Openize.Accord.Statistics.Distributions.Fitting
{
    using System;
    using Multivariate.Continuous;
    using Openize.Accord.Statistics.Distributions.Fitting.Base;

    /// <summary>
    ///   Smoothing rule function definition for <see cref="MultivariateEmpiricalDistribution">
    ///   Empirical distributions</see>. 
    /// </summary>
    /// 
    /// <param name="observations">The observations for the empirical distribution.</param>
    /// <param name="weights">The fractional importance for each sample. Those values must sum up to one.</param>
    /// <param name="repeats">The number of times each sample should be repeated.</param>
    /// 
    /// <returns>An estimative of the smoothing parameter.</returns>
    /// 
    public delegate double[,] MultivariateSmoothingRule(double[][] observations, double[] weights = null, int[] repeats = null);

    /// <summary>
    ///   Estimation options for <see cref="Accord.Statistics.Distributions.Multivariate
    ///   .MultivariateEmpiricalDistribution">Multivariate Empirical distributions</see>.
    /// </summary>
    /// 
    [Serializable]
    public class MultivariateEmpiricalOptions : IFittingOptions
    {
        /// <summary>
        ///   Gets or sets the smoothing rule used to compute the smoothing
        ///   parameter in the <see cref="MultivariateEmpiricalDistribution"/>.
        ///   Default is to use <see cref="MultivariateEmpiricalDistribution.SilvermanRule(double[][], double[], int[])">
        ///   Silverman's rule</see>.
        /// </summary>
        /// 
        public MultivariateSmoothingRule SmoothingRule { get; set; }

        /// <summary>
        ///   Gets or sets whether the empirical distribution should be take the
        ///   observation and weight vectors directly instead of making a copy
        ///   beforehand.
        /// </summary>
        /// 
        public bool InPlace { get; set; }

        /// <summary>
        ///   Initializes a new instance of the <see cref="NormalOptions"/> class.
        /// </summary>
        /// 
        public MultivariateEmpiricalOptions()
        {
            this.SmoothingRule = MultivariateEmpiricalDistribution.SilvermanRule;
            this.InPlace = false;
        }

        /// <summary>
        /// Creates a new object that is a copy of the current instance.
        /// </summary>
        /// <returns>
        /// A new object that is a copy of this instance.
        /// </returns>
        public object Clone()
        {
            return this.MemberwiseClone();
        }
    }
}
