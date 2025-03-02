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

namespace Openize.Accord.Statistics.Accord.MachineLearning.Classifiers
{
    using System;
    using Openize.Accord.Math.Matrix;
    using Openize.Accord.Core.MachineLearning;
    using Openize.Accord.Core.MachineLearning.Classifiers;

    /// <summary>
    ///   Base class for multi-class and multi-label classifiers.
    /// </summary>
    /// 
    /// <typeparam name="TInput">The data type for the input data. Default is double[].</typeparam>
    /// <typeparam name="TClasses">The data type for the classes. Default is int.</typeparam>
    /// 
    [Serializable]
    public abstract class ClassifierBase<TInput, TClasses> :
        TransformBase<TInput, TClasses>, 
        IClassifier<TInput, TClasses>
    {
        /// <summary>
        /// Gets the number of classes expected and recognized by the classifier.
        /// </summary>
        /// <value>The number of classes.</value>
        public virtual int NumberOfClasses { get; set; }

        /// <summary>
        ///   Computes a class-label decision for a given <paramref name="input"/>.
        /// </summary>
        /// 
        /// <param name="input">The input vector that should be classified into
        ///   one of the <see cref="ITransform.NumberOfOutputs"/> possible classes.</param>
        /// 
        /// <returns>A class-label that best described <paramref name="input"/> according
        /// to this classifier.</returns>
        /// 
        public abstract TClasses Decide(TInput input);

        /// <summary>
        ///   Computes class-label decisions for a given set of <paramref name="input"/> vectors.
        /// </summary>
        /// 
        /// <param name="input">The input vectors that should be classified into
        ///   one of the <see cref="ITransform.NumberOfOutputs"/> possible classes.</param>
        /// 
        /// <returns>The class-labels that best described each <paramref name="input"/>
        /// vector according to this classifier.</returns>
        /// 
        public TClasses[] Decide(TInput[] input)
        {
            return this.Decide(input, new TClasses[input.Length]);
        }


        /// <summary>
        ///  Computes a class-label decision for a given <paramref name="input" />.
        /// 
        /// </summary>
        /// 
        /// <param name="input">The input vector that should be classified into
        ///   one of the <see cref="ITransform.NumberOfOutputs" /> possible classes.</param>
        /// <param name="result">The location where to store the class-labels.</param>   
        /// 
        /// <returns>
        ///   A class-label that best described <paramref name="input" /> according
        ///   to this classifier.
        /// </returns>
        /// 
        public virtual TClasses[] Decide(TInput[] input, TClasses[] result)
        {
            for (int i = 0; i < input.Length; i++)
                result[i] = this.Decide(input[i]);
            return result;
        }

        /// <summary>
        /// Applies the transformation to an input, producing an associated output.
        /// </summary>
        /// <param name="input">The input data to which the transformation should be applied.</param>
        /// <returns>
        /// The output generated by applying this transformation to the given input.
        /// </returns>
        public override TClasses Transform(TInput input)
        {
            return this.Decide(input);
        }

        /// <summary>
        /// Applies the transformation to an input, producing an associated output.
        /// </summary>
        /// <param name="input">The input data to which the transformation should be applied.</param>
        /// <param name="result">The location where to store the class-labels.</param>
        /// <returns>
        /// The output generated by applying this transformation to the given input.
        /// </returns>
        public override TClasses[] Transform(TInput[] input, TClasses[] result)
        {
            return this.Decide(input, result);
        }




        internal T[] create<T>(TInput input)
        {
            return new T[this.NumberOfClasses];
        }

        internal T[][] create<T>(TInput[] input)
        {
            return Jagged.Create<T>(input.Length, this.NumberOfClasses);
        }

        internal T[] createOrReuse<T>(TInput input, T[] decision)
        {
            if (decision == null)
                decision = new T[this.NumberOfClasses];
            return decision;
        }

        internal T[][] createOrReuse<T>(TInput[] input, T[][] decision)
        {
            if (decision == null)
                decision = Jagged.Create<T>(input.Length, this.NumberOfClasses);
            return decision;
        }

        internal T[] createOrReuse<T>(TInput[] input, T[] decision)
        {
            if (decision == null)
                decision = new T[input.Length];
            return decision;
        }

    }
}
