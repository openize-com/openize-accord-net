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
    using Openize.Accord.Core;
    using Openize.Accord.Statistics.Distributions.Fitting.Base;
    using Univariate.Continuous;

    /// <summary>
    ///   Options for Survival distributions.
    /// </summary>
    /// 
    public class SurvivalOptions : IFittingOptions
    {
        /// <summary>
        ///   Default survival estimation method. Returns <see cref="SurvivalEstimator.FlemingHarrington"/>.
        /// </summary>
        /// 
        public const SurvivalEstimator DefaultSurvival = SurvivalEstimator.FlemingHarrington;

        /// <summary>
        ///   Gets or sets the values for
        ///   the right-censoring variable.
        /// </summary>
        /// 
        public SurvivalOutcome[] Outcome { get; set; }

        /// <summary>
        ///   Initializes a new instance of the <see cref="SurvivalOptions"/> class.
        /// </summary>
        /// 
        public SurvivalOptions()
        {
        }

        /// <summary>
        /// Creates a new object that is a copy of the current instance.
        /// </summary>
        /// <returns>
        /// A new object that is a copy of this instance.
        /// </returns>
        public virtual object Clone()
        {
            return this.MemberwiseClone();
        }

    }

    /// <summary>
    ///   Options for Empirical Hazard distributions.
    /// </summary>
    /// 
    public class EmpiricalHazardOptions : SurvivalOptions
    {
        /// <summary>
        ///   Default hazard estimator. Returns <see cref="HazardEstimator.BreslowNelsonAalen"/>.
        /// </summary>
        /// 
        public const HazardEstimator DefaultEstimator = HazardEstimator.BreslowNelsonAalen;

        /// <summary>
        ///   Default tie handling method. Returns <see cref="HazardTiesMethod.Efron"/>.
        /// </summary>
        /// 
        public const HazardTiesMethod DefaultTies = HazardTiesMethod.Efron;


        /// <summary>
        ///   Gets or sets the estimator to be used. Default is <see cref="HazardEstimator.BreslowNelsonAalen"/>.
        /// </summary>
        /// 
        public HazardEstimator Estimator { get; set; }

        /// <summary>
        ///   Gets or sets the tie handling method to be used. Default is <see cref="HazardTiesMethod.Efron"/>.
        /// </summary>
        /// 
        public HazardTiesMethod Ties { get; set; }


        /// <summary>
        ///   Initializes a new instance of the <see cref="EmpiricalHazardOptions"/> class.
        /// </summary>
        /// 
        public EmpiricalHazardOptions()
        {
            this.Estimator = DefaultEstimator;
            this.Ties = DefaultTies;
        }

        /// <summary>
        ///   Initializes a new instance of the <see cref="EmpiricalHazardOptions"/> class.
        /// </summary>
        /// 
        public EmpiricalHazardOptions(HazardEstimator estimator, SurvivalOutcome[] output)
        {
            this.Estimator = estimator;
            this.Outcome = output;
            this.Ties = DefaultTies;
        }

        /// <summary>
        ///   Initializes a new instance of the <see cref="EmpiricalHazardOptions"/> class.
        /// </summary>
        /// 
        public EmpiricalHazardOptions(HazardEstimator estimator, int[] output)
        {
            this.Estimator = estimator;
            this.Outcome = output.To<SurvivalOutcome[]>();
            this.Ties = HazardTiesMethod.Breslow;
        }

        /// <summary>
        ///   Initializes a new instance of the <see cref="EmpiricalHazardOptions"/> class.
        /// </summary>
        /// 
        public EmpiricalHazardOptions(HazardEstimator estimator, HazardTiesMethod ties, SurvivalOutcome[] outcome)
        {
            this.Estimator = estimator;
            this.Outcome = outcome;
            this.Ties = ties;
        }

        /// <summary>
        ///   Initializes a new instance of the <see cref="EmpiricalHazardOptions"/> class.
        /// </summary>
        /// 
        public EmpiricalHazardOptions(HazardEstimator estimator, HazardTiesMethod ties, int[] output)
        {
            this.Estimator = estimator;
            this.Outcome = output.To<SurvivalOutcome[]>();
            this.Ties = ties;
        }

        /// <summary>
        /// Creates a new object that is a copy of the current instance.
        /// </summary>
        /// <returns>
        /// A new object that is a copy of this instance.
        /// </returns>
        public override object Clone()
        {
            return this.MemberwiseClone();
        }
    }
}
