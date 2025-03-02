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

namespace Openize.Accord.Statistics.Testing.Power.Base
{
    using System;
    using Openize.Accord.Math.Optimization;
    using Testing.Base;

    /// <summary>
    ///   Base class for one sample power analysis methods. 
    ///   This class cannot be instantiated.
    /// </summary>
    /// 
    [Serializable]
    public abstract class BaseOneSamplePowerAnalysis : IPowerAnalysis
    {

        /// <summary>
        ///   Gets the test type.
        /// </summary>
        /// 
        public DistributionTail Tail { get; set; }

        /// <summary>
        ///   Gets or sets the power of the test, also
        ///   known as the (1-Beta error rate).
        /// </summary>
        /// 
        public double Power { get; set; }

        /// <summary>
        ///   Gets or sets the significance level
        ///   for the test. Also known as alpha.
        /// </summary>
        /// 
        public double Size { get; set; }

        /// <summary>
        ///   Gets or sets the number of samples 
        ///   considered in the test.
        /// </summary>
        /// 
        public double Samples { get; set; }

        /// <summary>
        ///   Gets or sets the effect size of the test.
        /// </summary>
        /// 
        public double Effect { get; set; }

        /// <summary>
        ///   Constructs a new power analysis for a one-sample test.
        /// </summary>
        /// 
        protected BaseOneSamplePowerAnalysis(DistributionTail tail)
        {
            this.Tail = tail;
            this.Size = 0.05;
            this.Power = 0.8;
            this.Samples = 2;
        }

        /// <summary>
        ///   Computes the power for a test with givens values of 
        ///   <see cref="Effect">effect size</see> and <see cref="Samples">
        ///   number of samples</see> under <see cref="Size"/>.
        /// </summary>
        /// 
        /// <returns>The power for the test under the given conditions.</returns>
        /// 
        public abstract void ComputePower();

        /// <summary>
        ///   Computes the minimum detectable effect size for the test
        ///   considering the power given in <see cref="Power"/>, the
        ///   number of samples in <see cref="Samples"/> and the significance
        ///   level <see cref="Size"/>.
        /// </summary>
        /// 
        /// <returns>The minimum detectable <see cref="Effect">effect
        /// size</see> for the test under the given conditions.</returns>
        /// 
        public virtual void ComputeEffect()
        {
            double requiredPower = this.Power;

            // Attempt to locate the optimal sample size
            // to attain the required power by locating
            // a zero in the difference function:

            double sol = BrentSearch.FindRoot(eff =>
            {
                this.Effect = eff;
                this.ComputePower();

                return requiredPower - this.Power;
            },

            lowerBound: 1e-5,
            upperBound: 1e+5);

            
            // Check it
            this.Effect = sol;
            this.ComputePower();

            double newPower = this.Power;
            this.Power = requiredPower;

            if (Math.Abs(requiredPower - newPower) > 1e-5)
                this.Effect = Double.NaN;
        }

        /// <summary>
        ///   Computes recommended sample size for the test to attain
        ///   the power indicated in <see cref="Power"/> considering
        ///   values of <see cref="Effect"/> and <see cref="Size"/>.
        /// </summary>
        /// 
        /// <returns>Recommended sample size for attaining the given
        /// <see cref="Power"/> for size effect <see cref="Effect"/>
        /// under the given <see cref="Size"/>.</returns>
        /// 
        public virtual void ComputeSamples()
        {
            double requiredPower = this.Power;

            // Attempt to locate the optimal sample size
            // to attain the required power by locating
            // a zero in the difference function:
            var search = new BrentSearch(n =>
            {
                this.Samples = n;
                this.ComputePower();

                return this.Power;
            },

            lowerBound: 2,
            upperBound: 1e+4);

            this.Samples = search.Find(requiredPower) ? search.Solution : double.NaN;
            this.Power = requiredPower;
        }

        /// <summary>
        ///   Creates a new object that is a copy of the current instance.
        /// </summary>
        /// 
        /// <returns>
        ///   A new object that is a copy of this instance.
        /// </returns>
        /// 
        public object Clone()
        {
            return this.MemberwiseClone();
        }

        /// <summary>
        ///   Converts the numeric power of this test to its equivalent string representation.
        /// </summary>
        /// 
        public string ToString(string format, IFormatProvider formatProvider)
        {
            return this.Power.ToString(format, formatProvider);
        }

        /// <summary>
        ///   Converts the numeric power of this test to its equivalent string representation.
        /// </summary>
        /// 
        public override string ToString()
        {
            return this.Power.ToString(System.Globalization.CultureInfo.CurrentCulture);
        }

        /// <summary>
        ///   Gets the minimum difference in the experiment units
        ///   to which it is possible to detect a difference.
        /// </summary>
        /// 
        /// <param name="standardDeviation">The standard deviation for the samples.</param>
        /// 
        /// <returns>The minimum difference in means which can be detected by the test.</returns>
        /// 
        public double GetDiferentiableUnits(double standardDeviation)
        {
            return this.Effect * standardDeviation;
        }

    }
}
