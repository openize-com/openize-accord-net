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

namespace Openize.Accord.Statistics.Accord.MachineLearning.Classifiers.Binary
{
    using System;
    using Openize.Accord.Math.Matrix;
    using Openize.Accord.Core.MachineLearning;
    using Openize.Accord.Core.MachineLearning.Classifiers;
    using Openize.Accord.Core.MachineLearning.Classifiers.Multiclass;
    using Openize.Accord.Core.MachineLearning.Classifiers.Multiclass.Binary;
    using Openize.Accord.Core.MachineLearning.Classifiers.Multilabel;
    using Openize.Accord.Math.Accord.Statistics;
    using Vector = Openize.Accord.Math.Vector.Vector;

    /// <summary>
    /// Base class for binary classifiers.
    /// </summary>
    /// <typeparam name="TInput">The data type for the input data. Default is double[].</typeparam>
    [Serializable]
    public abstract class BinaryClassifierBase<TInput> :
        ClassifierBase<TInput, bool>,
        IBinaryClassifier<TInput>
    {

        /// <summary>
        /// Initializes a new instance of the <see cref="BinaryClassifierBase{TInput}"/> class.
        /// </summary>
        public BinaryClassifierBase()
        {
            this.NumberOfOutputs = 1;
            this.NumberOfClasses = 2;
        }

        // Input

        double IClassifier<TInput, double>.Decide(TInput input)
        {
            return this.Decide(input).ToZeroOne();
        }

        int IClassifier<TInput, int>.Decide(TInput input)
        {
            return this.Decide(input).ToZeroOne();
        }
        int IMulticlassClassifier<TInput>.Decide(TInput input)
        {
            return this.ToMulticlass<int>().Decide(input);
        }

        int[] IMulticlassClassifier<TInput>.Decide(TInput[] input)
        {
            return this.ToMulticlass<int>().Decide(input);
        }
        bool[] IClassifier<TInput, bool[]>.Decide(TInput input)
        {
            return this.Decide(input, new bool[this.NumberOfOutputs]);
        }

        int[] IClassifier<TInput, int[]>.Decide(TInput input)
        {
            return this.ToMulticlass().Decide(input, new int[this.NumberOfOutputs]);
        }

        double[] IClassifier<TInput, double[]>.Decide(TInput input)
        {
            return this.ToMulticlass().Decide(input, new double[this.NumberOfOutputs]);
        }



        // Input[]

        double[] IClassifier<TInput, double>.Decide(TInput[] input)
        {
            return this.ToMulticlass().Decide(input, new double[input.Length]);
        }

        int[] IClassifier<TInput, int>.Decide(TInput[] input)
        {
            return this.ToMulticlass().Decide(input, new int[input.Length]);
        }

        bool[][] IClassifier<TInput, bool[]>.Decide(TInput[] input)
        {
            return this.ToMulticlass().Decide(input, this.create<bool>(input));
        }

        double[][] IClassifier<TInput, double[]>.Decide(TInput[] input)
        {
            return this.ToMulticlass().Decide(input, this.create<double>(input));
        }

        int[][] IClassifier<TInput, int[]>.Decide(TInput[] input)
        {
            return this.ToMulticlass().Decide(input, this.create<int>(input));
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
        public bool[] Decide(TInput input, bool[] result)
        {
            return Vector.OneHot<bool>(this.Decide(input), result);
        }

        int[] IMultilabelClassifier<TInput, int[]>.Decide(TInput input, int[] result)
        {
            return Vector.OneHot<int>(this.Decide(input), result);
        }

        double[] IMultilabelClassifier<TInput, double[]>.Decide(TInput input, double[] result)
        {
            return Vector.OneHot<double>(this.Decide(input), result);
        }



        // Input[], result[]

        int[] IClassifier<TInput, int>.Decide(TInput[] input, int[] result)
        {
            for (int i = 0; i < input.Length; i++)
                result[i] = this.Decide(input[i]).ToZeroOne();
            return result;
        }

        double[] IClassifier<TInput, double>.Decide(TInput[] input, double[] result)
        {
            for (int i = 0; i < input.Length; i++)
                result[i] = this.Decide(input[i]).ToZeroOne();
            return result;
        }

        
        bool[][] IClassifier<TInput, bool[]>.Decide(TInput[] input, bool[][] result)
        {
            return Jagged.OneHot<bool>(this.Decide(input), result);
        }


        int[][] IClassifier<TInput, int[]>.Decide(TInput[] input, int[][] result)
        {
            return Jagged.OneHot<int>(this.Decide(input), result);
        }

        double[][] IClassifier<TInput, double[]>.Decide(TInput[] input, double[][] result)
        {
            return Jagged.OneHot<double>(this.Decide(input), result);
        }






        // Transform

        int ICovariantTransform<TInput, int>.Transform(TInput input)
        {
            return ((IClassifier<TInput, int>)this).Decide(input);
        }

        
        double ICovariantTransform<TInput, double>.Transform(TInput input)
        {
            return ((IClassifier<TInput, double>)this).Decide(input);
        }

        
        double[] ICovariantTransform<TInput, double>.Transform(TInput[] input)
        {
            return this.Transform(input, new double[input.Length]);
        }

        
        double[] ICovariantTransform<TInput, double[]>.Transform(TInput input)
        {
            return this.Transform(input, this.create<double>(input));
        }


        int[] ICovariantTransform<TInput, int>.Transform(TInput[] input)
        {
            return this.Transform(input, new int[input.Length]);
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
            return this.ToMulticlass().Decide(input, result);
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
            return this.ToMulticlass().Decide(input, result);
        }

        /// <summary>
        /// Applies the transformation to an input, producing an associated output.
        /// </summary>
        /// <param name="input">The input data to which the transformation should be applied.</param>
        /// <param name="result">A location to store the output, avoiding unnecessary memory allocations.</param>
        /// <returns>
        /// The output generated by applying this transformation to the given input.
        /// </returns>
        virtual public int[] Transform(TInput[] input, int[] result)
        {
            return this.ToMulticlass().Decide(input, result);
        }

        /// <summary>
        /// Applies the transformation to an input, producing an associated output.
        /// </summary>
        /// <param name="input">The input data to which the transformation should be applied.</param>
        /// <param name="result">A location to store the output, avoiding unnecessary memory allocations.</param>
        /// <returns>
        /// The output generated by applying this transformation to the given input.
        /// </returns>

        virtual public double[] Transform(TInput[] input, double[] result)
        {
            return this.ToMulticlass().Decide(input, result);
        }

        /// <summary>
        /// Applies the transformation to an input, producing an associated output.
        /// </summary>
        /// <param name="input">The input data to which the transformation should be applied.</param>
        /// <param name="result">A location to store the output, avoiding unnecessary memory allocations.</param>
        /// <returns>
        /// The output generated by applying this transformation to the given input.
        /// </returns>

        virtual public bool[][] Transform(TInput[] input, bool[][] result)
        {
            return this.ToMulticlass().Decide(input, result);
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
            return this.ToMulticlass().Decide(input, result);
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
            return this.ToMulticlass().Decide(input, result);
        }



        /// <summary>
        /// Views this instance as a multi-class classifier,
        /// giving access to more advanced methods, such as the prediction
        /// of integer labels.
        /// </summary>
        /// <returns>
        /// This instance seen as an <see cref="IMulticlassClassifier{TInput}" />.
        /// </returns>
        public IMulticlassClassifier<TInput> ToMulticlass()
        {
            return (IMulticlassClassifier<TInput>)this;
        }

        /// <summary>
        /// Views this instance as a multi-class classifier,
        /// giving access to more advanced methods, such as the prediction
        /// of integer labels.
        /// </summary>
        /// <returns>
        /// This instance seen as an <see cref="IMulticlassClassifier{TInput}" />.
        /// </returns>
        public IMulticlassClassifier<TInput, T> ToMulticlass<T>()
        {
            return (IMulticlassClassifier<TInput, T>)this;
        }

        /// <summary>
        /// Views this instance as a multi-label classifier,
        /// giving access to more advanced methods, such as the prediction
        /// of one-hot vectors.
        /// </summary>
        /// <returns>
        /// This instance seen as an <see cref="IMultilabelLikelihoodClassifier{TInput}" />.
        /// </returns>
        public IMultilabelClassifier<TInput> ToMultilabel()
        {
            return (IMultilabelClassifier<TInput>)this;
        }
    }

}
