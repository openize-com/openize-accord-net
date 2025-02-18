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

namespace Openize.Accord.Statistics.Accord.MachineLearning.Classifiers.Multiple.Multiclass
{
    using System;
    using Openize.Accord.Math.Matrix;
    using Openize.Accord.Core.MachineLearning;
    using Openize.Accord.Core.MachineLearning.Classifiers;
    using Openize.Accord.Core.MachineLearning.Classifiers.Multiclass;
    using Openize.Accord.Core.MachineLearning.Classifiers.Multilabel;
    using Vector = Openize.Accord.Math.Vector.Vector;

    /// <summary>
    /// Base class for multi-class classifiers.
    /// </summary>
    /// 
    [Serializable]
    public abstract class MulticlassClassifierBase :
        MulticlassClassifierBase<double[]>,
        IMulticlassClassifier<int[]>,
        IMulticlassClassifier<float[]>
    {

        /// <summary>
        /// Computes a class-label decision for a given <paramref name="input" />.
        /// </summary>
        /// <param name="input">The input vector that should be classified into
        /// one of the <see cref="ITransform.NumberOfOutputs" /> possible classes.</param>
        /// <returns>
        /// A class-label that best described <paramref name="input" /> according
        /// to this classifier.
        /// </returns>
        public virtual int Decide(int[] input)
        {
            return this.Decide(input.ToDouble());
        }

        /// <summary>
        /// Computes a class-label decision for a given <paramref name="input" />.
        /// </summary>
        /// <param name="input">The input vector that should be classified into
        /// one of the <see cref="ITransform.NumberOfOutputs" /> possible classes.</param>
        /// <returns>
        /// A class-label that best described <paramref name="input" /> according
        /// to this classifier.
        /// </returns>
        public virtual int Decide(float[] input)
        {
            return this.Decide(input.ToDouble());
        }



        /// <summary>
        /// Computes a class-label decision for a given <paramref name="input" />.
        /// </summary>
        /// <param name="input">The input vector that should be classified into
        /// one of the <see cref="ITransform.NumberOfOutputs" /> possible classes.</param>
        /// <returns>
        /// A class-label that best described <paramref name="input" /> according
        /// to this classifier.
        /// </returns>
        public int[] Decide(int[][] input)
        {
            return this.Decide(input, new int[input.Length]);
        }

        /// <summary>
        /// Computes a class-label decision for a given <paramref name="input" />.
        /// </summary>
        /// <param name="input">The input vector that should be classified into
        /// one of the <see cref="ITransform.NumberOfOutputs" /> possible classes.</param>
        /// <returns>
        /// A class-label that best described <paramref name="input" /> according
        /// to this classifier.
        /// </returns>
        public int[] Decide(float[][] input)
        {
            return this.Decide(input, new int[input.Length]);
        }


        double IClassifier<int[], double>.Decide(int[] input)
        {
            return this.Decide(input);
        }

        double IClassifier<float[], double>.Decide(float[] input)
        {
            return this.Decide(input);
        }


        double[] IClassifier<int[], double>.Decide(int[][] input)
        {
            return this.Decide(input, new double[input.Length]);
        }

        double[] IClassifier<float[], double>.Decide(float[][] input)
        {
            return this.Decide(input, new double[input.Length]);
        }


        double[][] IClassifier<int[], double[]>.Decide(int[][] input)
        {
            return this.Decide(input, Jagged.Create<double>(input.Length, this.NumberOfOutputs));
        }

        double[][] IClassifier<float[], double[]>.Decide(float[][] input)
        {
            return this.Decide(input, Jagged.Create<double>(input.Length, this.NumberOfOutputs));
        }

        bool[][] IClassifier<int[], bool[]>.Decide(int[][] input)
        {
            return this.Decide(input, Jagged.Create<bool>(input.Length, this.NumberOfOutputs));
        }

        bool[][] IClassifier<float[], bool[]>.Decide(float[][] input)
        {
            return this.Decide(input, Jagged.Create<bool>(input.Length, this.NumberOfOutputs));
        }

        double[] IClassifier<int[], double[]>.Decide(int[] input)
        {
            return this.Decide(input, new double[this.NumberOfOutputs]);
        }

        double[] IClassifier<float[], double[]>.Decide(float[] input)
        {
            return this.Decide(input, new double[this.NumberOfOutputs]);
        }

        bool[] IClassifier<int[], bool[]>.Decide(int[] input)
        {
            return this.Decide(input, new bool[this.NumberOfOutputs]);
        }

        bool[] IClassifier<float[], bool[]>.Decide(float[] input)
        {
            return this.Decide(input, new bool[this.NumberOfOutputs]);
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
        public int[] Decide(int[][] input, int[] result)
        {
            var cast = new double[this.NumberOfInputs];
            for (int i = 0; i < input.Length; i++)
            {
                input[i].ToDouble(result: cast);
                result[i] = this.Decide(cast);
            }
            return result;
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
        public double[] Decide(float[] input, double[] result)
        {
            return Vector.OneHot(this.Decide(input), result);
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
        public double[] Decide(int[] input, double[] result)
        {
            return Vector.OneHot(this.Decide(input), result);
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
        public bool[] Decide(float[] input, bool[] result)
        {
            return Vector.OneHot(this.Decide(input), result);
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
        public bool[] Decide(int[] input, bool[] result)
        {
            return Vector.OneHot(this.Decide(input), result);
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
        public int[] Decide(float[][] input, int[] result)
        {
            var cast = new double[this.NumberOfInputs];
            for (int i = 0; i < input.Length; i++)
            {
                input[i].ToDouble(result: cast);
                result[i] = this.Decide(cast);
            }
            return result;
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
        public int[][] Decide(float[][] input, int[][] result)
        {
            var cast = new double[this.NumberOfInputs];
            for (int i = 0; i < input.Length; i++)
            {
                input[i].ToDouble(result: cast);
                Vector.OneHot(this.Decide(cast), result[i]);
            }
            return result;
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
        public double[] Decide(int[][] input, double[] result)
        {
            var cast = new double[this.NumberOfInputs];
            for (int i = 0; i < input.Length; i++)
            {
                input[i].ToDouble(result: cast);
                result[i] = this.Decide(cast);
            }
            return result;
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
        public double[] Decide(float[][] input, double[] result)
        {
            var cast = new double[this.NumberOfInputs];
            for (int i = 0; i < input.Length; i++)
            {
                input[i].ToDouble(result: cast);
                result[i] = this.Decide(cast);
            }
            return result;
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
        public int[][] Decide(int[][] input, int[][] result)
        {
            var cast = new double[this.NumberOfInputs];
            for (int i = 0; i < input.Length; i++)
            {
                input[i].ToDouble(result: cast);
                Vector.OneHot(this.Decide(cast), result[i]);
            }
            return result;
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
        public double[][] Decide(float[][] input, double[][] result)
        {
            var cast = new double[this.NumberOfInputs];
            for (int i = 0; i < input.Length; i++)
            {
                input[i].ToDouble(result: cast);
                Vector.OneHot(this.Decide(cast), result[i]);
            }
            return result;
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
        public double[][] Decide(int[][] input, double[][] result)
        {
            var cast = new double[this.NumberOfInputs];
            for (int i = 0; i < input.Length; i++)
            {
                input[i].ToDouble(result: cast);
                Vector.OneHot(this.Decide(cast), result[i]);
            }
            return result;
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
        public bool[][] Decide(float[][] input, bool[][] result)
        {
            var cast = new double[this.NumberOfInputs];
            for (int i = 0; i < input.Length; i++)
            {
                input[i].ToDouble(result: cast);
                Vector.OneHot(this.Decide(cast), result[i]);
            }
            return result;
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
        public bool[][] Decide(int[][] input, bool[][] result)
        {
            var cast = new double[this.NumberOfInputs];
            for (int i = 0; i < input.Length; i++)
            {
                input[i].ToDouble(result: cast);
                Vector.OneHot(this.Decide(cast), result[i]);
            }
            return result;
        }



        int[] IClassifier<float[], int[]>.Decide(float[] input)
        {
            var t = (IMultilabelClassifier<float[]>)this;
            return t.Decide(input, new int[this.NumberOfOutputs]);
        }

        int[] IClassifier<int[], int[]>.Decide(int[] input)
        {
            var t = (IMultilabelClassifier<double[]>)this;
            return t.Decide(input.ToDouble(), new int[this.NumberOfOutputs]);
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
        public int[] Decide(float[] input, int[] result)
        {
            return this.Decide(input.ToDouble(), result);
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
        public int[] Decide(int[] input, int[] result)
        {
            return this.Decide(input.ToDouble(), result);
        }



        int[][] IClassifier<float[], int[]>.Decide(float[][] input)
        {
            return this.Decide(input, Jagged.Create<int>(input.Length, this.NumberOfOutputs));
        }

        int[][] IClassifier<int[], int[]>.Decide(int[][] input)
        {
            return this.Decide(input, Jagged.Create<int>(input.Length, this.NumberOfOutputs));
        }








        /// <summary>
        /// Applies the transformation to an input, producing an associated output.
        /// </summary>
        /// <param name="input">The input data to which the transformation should be applied.</param>
        /// <returns>
        /// The output generated by applying this transformation to the given input.
        /// </returns>
        public int Transform(int[] input)
        {
            return this.Decide(input);
        }

        /// <summary>
        /// Applies the transformation to an input, producing an associated output.
        /// </summary>
        /// <param name="input">The input data to which the transformation should be applied.</param>
        /// <param name="result">A location to store the output, avoiding unnecessary memory allocations.</param>
        /// <returns>
        /// The output generated by applying this transformation to the given input.
        /// </returns>
        public double[] Transform(int[] input, double[] result)
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
        public double[][] Transform(int[][] input, double[][] result)
        {
            return this.Decide(input, result);
        }








        /// <summary>
        /// Applies the transformation to an input, producing an associated output.
        /// </summary>
        /// <param name="input">The input data to which the transformation should be applied.</param>
        /// <returns>
        /// The output generated by applying this transformation to the given input.
        /// </returns>
        public int Transform(float[] input)
        {
            return this.Decide(input);
        }

        /// <summary>
        /// Applies the transformation to an input, producing an associated output.
        /// </summary>
        /// <param name="input">The input data to which the transformation should be applied.</param>
        /// <param name="result">A location to store the output, avoiding unnecessary memory allocations.</param>
        /// <returns>
        /// The output generated by applying this transformation to the given input.
        /// </returns>
        public int Transform(float[] input, int result)
        {
            return this.Decide(input);
        }

        /// <summary>
        /// Applies the transformation to an input, producing an associated output.
        /// </summary>
        /// <param name="input">The input data to which the transformation should be applied.</param>
        /// <returns>
        /// The output generated by applying this transformation to the given input.
        /// </returns>
        public int[] Transform(float[][] input)
        {
            return this.Decide(input);
        }

        /// <summary>
        /// Applies the transformation to an input, producing an associated output.
        /// </summary>
        /// <param name="input">The input data to which the transformation should be applied.</param>
        /// <param name="result">A location to store the output, avoiding unnecessary memory allocations.</param>
        /// <returns>
        /// The output generated by applying this transformation to the given input.
        /// </returns>
        public int[] Transform(float[][] input, int[] result)
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
        public double[] Transform(float[] input, double[] result)
        {
            return Vector.OneHot(this.Decide(input), result);
        }

        /// <summary>
        /// Applies the transformation to an input, producing an associated output.
        /// </summary>
        /// <param name="input">The input data to which the transformation should be applied.</param>
        /// <param name="result">A location to store the output, avoiding unnecessary memory allocations.</param>
        /// <returns>
        /// The output generated by applying this transformation to the given input.
        /// </returns>
        public double[][] Transform(float[][] input, double[][] result)
        {
            return Jagged.OneHot(this.Decide(input), result);
        }

        /// <summary>
        /// Applies the transformation to an input, producing an associated output.
        /// </summary>
        /// <param name="input">The input data to which the transformation should be applied.</param>
        /// <param name="result">A location to store the output, avoiding unnecessary memory allocations.</param>
        /// <returns>
        /// The output generated by applying this transformation to the given input.
        /// </returns>
        public bool[] Transform(float[] input, bool[] result)
        {
            return Vector.OneHot<bool>(this.Decide(input), result);
        }

        /// <summary>
        /// Applies the transformation to an input, producing an associated output.
        /// </summary>
        /// <param name="input">The input data to which the transformation should be applied.</param>
        /// <param name="result">A location to store the output, avoiding unnecessary memory allocations.</param>
        /// <returns>
        /// The output generated by applying this transformation to the given input.
        /// </returns>
        public bool[][] Transform(float[][] input, bool[][] result)
        {
            return Jagged.OneHot<bool>(this.Decide(input), result);
        }

        /// <summary>
        /// Applies the transformation to an input, producing an associated output.
        /// </summary>
        /// <param name="input">The input data to which the transformation should be applied.</param>
        /// <param name="result">A location to store the output, avoiding unnecessary memory allocations.</param>
        /// <returns>
        /// The output generated by applying this transformation to the given input.
        /// </returns>
        public int[] Transform(float[] input, int[] result)
        {
            return Vector.OneHot<int>(this.Decide(input), result);
        }

        /// <summary>
        /// Applies the transformation to an input, producing an associated output.
        /// </summary>
        /// <param name="input">The input data to which the transformation should be applied.</param>
        /// <param name="result">A location to store the output, avoiding unnecessary memory allocations.</param>
        /// <returns>
        /// The output generated by applying this transformation to the given input.
        /// </returns>
        public int[][] Transform(float[][] input, int[][] result)
        {
            return Jagged.OneHot<int>(this.Decide(input), result);
        }

        /// <summary>
        /// Applies the transformation to an input, producing an associated output.
        /// </summary>
        /// <param name="input">The input data to which the transformation should be applied.</param>
        /// <param name="result">A location to store the output, avoiding unnecessary memory allocations.</param>
        /// <returns>
        /// The output generated by applying this transformation to the given input.
        /// </returns>
        public int[][] Transform(int[][] input, int[][] result)
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
        public int[] Transform(int[] input, int[] result)
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
        public int[] Transform(int[][] input, int[] result)
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
        public bool[][] Transform(int[][] input, bool[][] result)
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
        public double[] Transform(int[][] input, double[] result)
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
        public double[] Transform(float[][] input, double[] result)
        {
            return this.Decide(input, result);
        }




        double[] ICovariantTransform<int[], double[]>.Transform(int[] input)
        {
            return Vector.OneHot(this.Decide(input), this.NumberOfOutputs);
        }

        double[][] ICovariantTransform<int[], double[]>.Transform(int[][] input)
        {
            return Jagged.OneHot(this.Decide(input), this.NumberOfOutputs);
        }

        double[] ICovariantTransform<float[], double[]>.Transform(float[] input)
        {
            return Vector.OneHot(this.Decide(input), this.NumberOfOutputs);
        }

        double[][] ICovariantTransform<float[], double[]>.Transform(float[][] input)
        {
            return Jagged.OneHot(this.Decide(input), this.NumberOfOutputs);
        }

        bool[] ICovariantTransform<float[], bool[]>.Transform(float[] input)
        {
            return Vector.OneHot<bool>(this.Decide(input), this.NumberOfOutputs);
        }

        bool[][] ICovariantTransform<float[], bool[]>.Transform(float[][] input)
        {
            return Jagged.OneHot<bool>(this.Decide(input), this.NumberOfOutputs);
        }

        int[] ICovariantTransform<float[], int[]>.Transform(float[] input)
        {
            return Vector.OneHot<int>(this.Decide(input), this.NumberOfOutputs);
        }

        int[][] ICovariantTransform<float[], int[]>.Transform(float[][] input)
        {
            return Jagged.OneHot<int>(this.Decide(input), this.NumberOfOutputs);
        }

        int[] ICovariantTransform<int[], int[]>.Transform(int[] input)
        {
            return Vector.OneHot<int>(this.Decide(input), this.NumberOfOutputs);
        }

        int[][] ICovariantTransform<int[], int[]>.Transform(int[][] input)
        {
            return Jagged.OneHot<int>(this.Decide(input), this.NumberOfOutputs);
        }

        bool[] ICovariantTransform<int[], bool[]>.Transform(int[] input)
        {
            return Vector.OneHot<bool>(this.Decide(input), this.NumberOfOutputs);
        }

        bool[][] ICovariantTransform<int[], bool[]>.Transform(int[][] input)
        {
            return Jagged.OneHot<bool>(this.Decide(input), this.NumberOfOutputs);
        }




        int ICovariantTransform<int[], int>.Transform(int[] input)
        {
            return this.Decide(input);
        }

        int[] ICovariantTransform<int[], int>.Transform(int[][] input)
        {
            return this.Decide(input);
        }



        double ICovariantTransform<int[], double>.Transform(int[] input)
        {
            return this.Decide(input);
        }

        double[] ICovariantTransform<int[], double>.Transform(int[][] input)
        {
            var result = new double[input.Length];
            return this.Decide(input, result);
        }


        double ICovariantTransform<float[], double>.Transform(float[] input)
        {
            return this.Decide(input);
        }

        double[] ICovariantTransform<float[], double>.Transform(float[][] input)
        {
            var result = new double[input.Length];
            return this.Decide(input, result);
        }


        IMultilabelClassifier<int[]> IMulticlassClassifier<int[]>.ToMultilabel()
        {
            return (IMultilabelClassifier<int[]>)this;
        }

        IMultilabelClassifier<float[]> IMulticlassClassifier<float[]>.ToMultilabel()
        {
            return (IMultilabelClassifier<float[]>)this;
        }
    }


}
