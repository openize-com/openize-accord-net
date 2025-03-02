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

namespace Openize.Accord.Statistics.Accord.MachineLearning.Classifiers.Multiple.Multilabel
{
    using System;
    using Openize.Accord.Math.Matrix;
    using Openize.Accord.Core.MachineLearning;
    using Openize.Accord.Core.MachineLearning.Classifiers;
    using Openize.Accord.Core.MachineLearning.Classifiers.Multilabel;
    using Vector = Openize.Accord.Math.Vector.Vector;

    /// <summary>
    ///   Base class for <see cref="IMultilabelClassifier{TInput}">
    ///   multi-label classifiers</see>.
    /// </summary>
    /// 
    /// <typeparam name="TInput">The data type for the input data. Default is double[].</typeparam>
    /// 
    [Serializable]
    public abstract class MultilabelClassifierBase<TInput> :
        ClassifierBase<TInput, bool[]>,
        IMultilabelClassifier<TInput>
    {

        // Main, overridable methods

        /// <summary>
        ///   Computes whether a class label applies to an <paramref name="input"/> vector.
        /// </summary>
        /// 
        /// <param name="input">The input vectors that should be classified as
        ///   any of the <see cref="ITransform.NumberOfOutputs"/> possible classes.</param>
        /// <param name="classIndex">The class label index to be tested.</param>
        /// 
        /// <returns>A boolean value indicating whether the given <paramref name="classIndex">
        /// class label</paramref> applies to the <paramref name="input"/> vector.</returns>
        /// 
        public abstract bool Decide(TInput input, int classIndex);

        /// <summary>
        /// Computes a class-label decision for a given <paramref name="input" />.
        /// </summary>
        /// <param name="input">The input vector that should be classified into
        /// one of the <see cref="ITransform.NumberOfOutputs" /> possible classes.</param>
        /// <returns>
        /// A class-label that best described <paramref name="input" /> according
        /// to this classifier.
        /// </returns>
        public override bool[] Decide(TInput input)
        {
            return this.Decide(input, new bool[this.NumberOfOutputs]);
        }

        /// <summary>
        /// Computes class-label decisions for the given <paramref name="input" />.
        /// </summary>
        /// <param name="input">The input vectors that should be classified as
        /// any of the <see cref="ITransform.NumberOfOutputs" /> possible classes.</param>
        /// <param name="result">The location where to store the class-labels.</param>
        /// <returns>
        /// A set of class-labels that best describe the <paramref name="input" />
        /// vectors according to this classifier.
        /// </returns>
        public virtual bool[] Decide(TInput input, bool[] result)
        {
            for (int i = 0; i < result.Length; i++)
                result[i] = this.Decide(input, i);
            return result;
        }



        // Input

        int[] IClassifier<TInput, int[]>.Decide(TInput input)
        {
            return this.Decide(input, new int[this.NumberOfOutputs]);
        }

        double[] IClassifier<TInput, double[]>.Decide(TInput input)
        {
            return this.Decide(input, new double[this.NumberOfOutputs]);
        }



        // Input[]

        double[][] IClassifier<TInput, double[]>.Decide(TInput[] input)
        {
            return this.Decide(input, this.create<double>(input));
        }

        int[][] IClassifier<TInput, int[]>.Decide(TInput[] input)
        {
            return this.Decide(input, this.create<int>(input));
        }




        // Input, result

        /// <summary>
        /// Computes class-label decisions for the given <paramref name="input" />.
        /// </summary>
        /// <param name="input">The input vectors that should be classified as
        /// any of the <see cref="ITransform.NumberOfOutputs" /> possible classes.</param>
        /// <param name="result">The location where to store the class-labels.</param>
        /// <returns>
        /// A set of class-labels that best describe the <paramref name="input" />
        /// vectors according to this classifier.
        /// </returns>
        public int[] Decide(TInput input, int[] result)
        {
            return Vector.KHot<int>(this.Decide(input), result);
        }


        /// <summary>
        /// Computes class-label decisions for the given <paramref name="input" />.
        /// </summary>
        /// <param name="input">The input vectors that should be classified as
        /// any of the <see cref="ITransform.NumberOfOutputs" /> possible classes.</param>
        /// <param name="result">The location where to store the class-labels.</param>
        /// <returns>
        /// A set of class-labels that best describe the <paramref name="input" />
        /// vectors according to this classifier.
        /// </returns>
        public double[] Decide(TInput input, double[] result)
        {
            return Vector.KHot<double>(this.Decide(input), result);
        }



        // Input[], result[]

        /// <summary>
        /// Computes a class-label decision for a given <paramref name="input" />.
        /// </summary>
        /// <param name="input">The input vector that should be classified into
        /// one of the <see cref="ITransform.NumberOfOutputs" /> possible classes.</param>
        /// <param name="result">The location where to store the class-labels.</param>
        /// <returns>
        /// A class-label that best described <paramref name="input" /> according
        /// to this classifier.
        /// </returns>
        public int[][] Decide(TInput[] input, int[][] result)
        {
            return Jagged.KHot<int>(this.Decide(input), result);
        }

        /// <summary>
        /// Computes a class-label decision for a given <paramref name="input" />.
        /// </summary>
        /// <param name="input">The input vector that should be classified into
        /// one of the <see cref="ITransform.NumberOfOutputs" /> possible classes.</param>
        /// <param name="result">The location where to store the class-labels.</param>
        /// <returns>
        /// A class-label that best described <paramref name="input" /> according
        /// to this classifier.
        /// </returns>
        public double[][] Decide(TInput[] input, double[][] result)
        {
            return Jagged.KHot<double>(this.Decide(input), result);
        }

        


        // Transform

        double[] ICovariantTransform<TInput, double[]>.Transform(TInput input)
        {
            return this.Transform(input, this.create<double>(input));
        }

        double[][] ICovariantTransform<TInput, double[]>.Transform(TInput[] input)
        {
            return this.Transform(input, this.create<double>(input));
        }

        bool[] ICovariantTransform<TInput, bool[]>.Transform(TInput input)
        {
            return this.Transform(input, this.create<bool>(input));
        }

        bool[][] ICovariantTransform<TInput, bool[]>.Transform(TInput[] input)
        {
            return this.Transform(input, this.create<bool>(input));
        }

        int[] ICovariantTransform<TInput, int[]>.Transform(TInput input)
        {
            return this.Transform(input, this.create<int>(input));
        }

        int[][] ICovariantTransform<TInput, int[]>.Transform(TInput[] input)
        {
            return this.Transform(input, this.create<int>(input));
        }



        /// <summary>
        /// Applies the transformation to an input, producing an associated output.
        /// </summary>
        /// <param name="input">The input data to which the transformation should be applied.</param>
        /// <param name="result">A location to store the output, avoiding unnecessary memory allocations.</param>
        /// <returns>
        /// The output generated by applying this transformation to the given input.
        /// </returns>
        virtual public bool[] Transform(TInput input, bool[] result)
        {
            return this.Decide(input, result);
        }

        /// <summary>
        /// Applies the transformation to an input, producing an associated output.
        /// </summary>
        /// <param name="input">The input data to which the transformation should be applied.</param>
        /// <param name="result">A location to store the output, avoiding unnecessary memory allocations.</param>
        /// <returns>
        /// The output generated by applying this transformation to the given input.
        /// </returns>
        virtual public int[] Transform(TInput input, int[] result)
        {
            return this.Decide(input, result);
        }

        /// <summary>
        /// Applies the transformation to an input, producing an associated output.
        /// </summary>
        /// <param name="input">The input data to which the transformation should be applied.</param>
        /// <param name="result">A location to store the output, avoiding unnecessary memory allocations.</param>
        /// <returns>
        /// The output generated by applying this transformation to the given input.
        /// </returns>
        virtual public double[] Transform(TInput input, double[] result)
        {
            return this.Decide(input, result);
        }

        /// <summary>
        /// Applies the transformation to an input, producing an associated output.
        /// </summary>
        /// <param name="input">The input data to which the transformation should be applied.</param>
        /// <param name="result">A location to store the output, avoiding unnecessary memory allocations.</param>
        /// <returns>
        /// The output generated by applying this transformation to the given input.
        /// </returns>
        public override bool[][] Transform(TInput[] input, bool[][] result)
        {
            return this.Decide(input, result);
        }

        /// <summary>
        /// Applies the transformation to an input, producing an associated output.
        /// </summary>
        /// <param name="input">The input data to which the transformation should be applied.</param>
        /// <param name="result">A location to store the output, avoiding unnecessary memory allocations.</param>
        /// <returns>
        /// The output generated by applying this transformation to the given input.
        /// </returns>
        virtual public int[][] Transform(TInput[] input, int[][] result)
        {
            return this.Decide(input, result);
        }

        /// <summary>
        /// Applies the transformation to an input, producing an associated output.
        /// </summary>
        /// <param name="input">The input data to which the transformation should be applied.</param>
        /// <param name="result">A location to store the output, avoiding unnecessary memory allocations.</param>
        /// <returns>
        /// The output generated by applying this transformation to the given input.
        /// </returns>
        virtual public double[][] Transform(TInput[] input, double[][] result)
        {
            return this.Decide(input, result);
        }
    }
}
