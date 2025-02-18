// Accord Statistics Library
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

#pragma warning disable 612, 618

namespace Openize.Accord.Statistics.Models.Fields.Learning.Hidden
{
    using System;
    using Accord.MachineLearning.Learning;
    using Openize.Accord.Core.MachineLearning;
    using Openize.Accord.Math.Accord.Statistics;
    using Openize.Accord.Math.Optimization.Unconstrained;

    /// <summary>
    ///   Quasi-Newton (L-BFGS) learning algorithm for <see cref="HiddenConditionalRandomField{T}">
    ///   Hidden Conditional Hidden Fields</see>.
    /// </summary>
    /// 
    /// <typeparam name="T">The type of the observations.</typeparam>
    /// 
    /// <example>
    ///   <code source="Unit Tests\Accord.Tests.Statistics\Models\Fields\HiddenConditionalRandomFieldTest.cs" region="doc_learn_1" />
    ///   <code source="Unit Tests\Accord.Tests.Statistics\Models\Fields\HiddenConditionalRandomFieldTest.cs" region="doc_learn_2" />
    ///   <code source="Unit Tests\Accord.Tests.Statistics\Models\Fields\HiddenConditionalRandomFieldTest.cs" region="doc_learn_3" />
    ///   
    ///   <para>
    ///   The next example shows how to use the learning algorithms in a real-world dataset,
    ///   including training and testing in separate sets and evaluating its performance:</para>
    ///   <code source="Unit Tests\Accord.Tests.Statistics\Models\Fields\Learning\NormalQuasiNewtonHiddenLearningTest.cs" region="doc_learn_pendigits" />
    /// </example>
    /// 
    /// <seealso cref="HiddenGradientDescentLearning{T}"/>
    /// <seealso cref="HiddenResilientGradientLearning{T}"/>
    /// 
    public class HiddenQuasiNewtonLearning<T> : BaseHiddenGradientOptimizationLearning<T, BoundedBroydenFletcherGoldfarbShanno>,
        ISupervisedLearning<HiddenConditionalRandomField<T>, T[], int>, IParallel,
        IHiddenConditionalRandomFieldLearning<T>, IConvergenceLearning, IDisposable
    {
        int IConvergenceLearning.Iterations
        {
            get
            {
                if (this.Optimizer == null)
                    return 0;
                return this.Optimizer.Iterations;
            }
            set { throw new NotImplementedException(); }
        }

        /// <summary>
        /// Gets the current iteration number.
        /// </summary>
        /// <value>The current iteration.</value>
        public int CurrentIteration
        {
            get
            {
                if (this.Optimizer != null)
                    return this.Optimizer.Iterations;
                return 0;
            }
        }

        /// <summary>
        ///   Constructs a new L-BFGS learning algorithm.
        /// </summary>
        /// 
        public HiddenQuasiNewtonLearning()
        {
        }

        /// <summary>
        ///   Constructs a new L-BFGS learning algorithm.
        /// </summary>
        /// 
        public HiddenQuasiNewtonLearning(HiddenConditionalRandomField<T> model)
        {
            this.Model = model;
        }

        /// <summary>
        /// Inheritors of this class should create the optimization algorithm in this
        /// method, using the current <see cref="P:Openize.Accord.Statistics.Models.Fields.Learning.Hidden.BaseHiddenGradientOptimizationLearning`2.MaxIterations" /> and <see cref="P:Openize.Accord.Statistics.Models.Fields.Learning.Hidden.BaseHiddenGradientOptimizationLearning`2.Tolerance" />
        /// settings.
        /// </summary>
        /// <returns>BoundedBroydenFletcherGoldfarbShanno.</returns>
        protected override BoundedBroydenFletcherGoldfarbShanno CreateOptimizer()
        {
            var lbfgs = new BoundedBroydenFletcherGoldfarbShanno(this.Model.Function.Weights.Length)
            {
                FunctionTolerance = this.Tolerance,
                MaxIterations = this.MaxIterations,
            };

            for (int i = 0; i < lbfgs.UpperBounds.Length; i++)
            {
                lbfgs.UpperBounds[i] = 1e10;
                lbfgs.LowerBounds[i] = -1e100;
            }

            return lbfgs;
        }

    }
}
