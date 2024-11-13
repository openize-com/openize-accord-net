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

namespace FileFormat.Accord.Statistics.Accord.MachineLearning.Classifiers.Binary
{
    using System;
    using FileFormat.Accord.Core.MachineLearning;
    using FileFormat.Accord.Core.MachineLearning.Classifiers;
    using FileFormat.Accord.Core.MachineLearning.Classifiers.Multiclass;
    using FileFormat.Accord.Core.MachineLearning.Classifiers.Multiclass.Binary;
    using FileFormat.Accord.Core.MachineLearning.Classifiers.Multilabel;
    using global::Accord.Math;
    using Math.Accord.Statistics;
    using Math.Core;
    using Math.Matrix;
    using Vector = Math.Vector.Vector;

    /// <summary>
    ///   Base class for <see cref="IBinaryLikelihoodClassifier{TInput}">
    ///   generative binary classifiers</see>.
    /// </summary>
    /// 
    /// <typeparam name="TInput">The data type for the input data. Default is double[].</typeparam>
    /// 
    [Serializable]
    public abstract class BinaryLikelihoodClassifierBase<TInput> :
        BinaryScoreClassifierBase<TInput>,
        IBinaryLikelihoodClassifier<TInput>
    {

        internal const double PROBABILITY_DECISION_THRESHOLD = 0.5;
        internal readonly double LOGLIKELIHOOD_DECISION_THRESHOLD = Math.Log(0.5);



        // Main, overridable methods

        /// <summary>
        /// Predicts a class label vector for the given input vectors, returning the
        /// log-likelihood that the input vector belongs to its predicted class.
        /// </summary>
        /// <param name="input">The input vector.</param>
        /// <param name="result">An array where the log-likelihoods will be stored,
        /// avoiding unnecessary memory allocations.</param>
        /// <returns>System.Double[].</returns>
        public abstract double[] LogLikelihood(TInput[] input, double[] result);


        /// <summary>
        /// Computes a numerical score measuring the association between
        /// the given <paramref name="input" /> vector and each class.
        /// </summary>
        /// <param name="input">The input vector.</param>
        /// <param name="result">An array where the result will be stored,
        /// avoiding unnecessary memory allocations.</param>
        /// <returns>System.Double[].</returns>
        public override double[] Score(TInput[] input, double[] result)
        {
            return this.LogLikelihood(input, result);
        }





        #region LogLikelihood

        // Input, classIndex
        double IMultilabelLikelihoodClassifier<TInput>.LogLikelihood(TInput input, int classIndex)
        {
            return this.ToMultilabel().LogLikelihoods(input)[classIndex];
        }

        double[] IMultilabelLikelihoodClassifier<TInput>.LogLikelihood(TInput[] input, int classIndex)
        {
            return ((IMulticlassLikelihoodClassifier<TInput>)this).LogLikelihood(input, classIndex, new double[input.Length]);
        }

        double[] IMultilabelLikelihoodClassifier<TInput>.LogLikelihood(TInput[] input, int classIndex, double[] result)
        {
            var temp = new double[this.NumberOfClasses];
            for (int i = 0; i < input.Length; i++)
            {
                this.ToMultilabel().LogLikelihoods(input[i], result: temp);
                result[i] = temp[classIndex];
            }

            return result;
        }


        // Input, classIndex[]
        double[] IMultilabelLikelihoodClassifier<TInput>.LogLikelihood(TInput[] input, int[] classIndex)
        {
            return ((IMulticlassLikelihoodClassifier<TInput>)this).LogLikelihood(input, classIndex, new double[input.Length]);
        }

        double[] IMultilabelLikelihoodClassifier<TInput>.LogLikelihood(TInput[] input, int[] classIndex, double[] result)
        {
            var temp = new double[this.NumberOfClasses];
            for (int i = 0; i < input.Length; i++)
            {
                this.ToMultilabel().LogLikelihoods(input[i], result: temp);
                result[i] = temp[classIndex[i]];
            }

            return result;
        }


        // Input

        /// <summary>
        /// Predicts a class label vector for the given input vector, returning the
        /// log-likelihood that the input vector belongs to its predicted class.
        /// </summary>
        /// <param name="input">The input vector.</param>
        public double LogLikelihood(TInput input)
        {
            return this.LogLikelihood(new TInput[] { input })[0];
        }

        /// <summary>
        /// Computes the log-likelihood that the given input
        /// vector belongs to each of the possible classes.
        /// </summary>
        /// <param name="input">The input vector.</param>
        /// <returns>System.Double[].</returns>
        public double[] LogLikelihoods(TInput input)
        {
            return this.ToMultilabel().LogLikelihoods(input, new double[this.NumberOfClasses]);
        }

        /// <summary>
        /// Predicts a class label vector for the given input vector, returning the
        /// log-likelihood that the input vector belongs to its predicted class.
        /// </summary>
        /// <param name="input">The input vector.</param>
        public double[] LogLikelihood(TInput[] input)
        {
            return this.LogLikelihood(input, new double[input.Length]);
        }

        /// <summary>
        /// Computes the log-likelihoods that the given input
        /// vectors belongs to each of the possible classes.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <returns>System.Double[][].</returns>
        public double[][] LogLikelihoods(TInput[] input)
        {
            return this.ToMultilabel().LogLikelihoods(input, this.create<double>(input));
        }


        // Input, result

        /// <summary>
        /// Computes the log-likelihood that the given input
        /// vector belongs to each of the possible classes.
        /// </summary>
        /// <param name="input">The input vector.</param>
        /// <param name="result">An array where the log-likelihoods will be stored,
        /// avoiding unnecessary memory allocations.</param>
        /// <returns>System.Double[].</returns>
        public double[] LogLikelihoods(TInput input, double[] result)
        {
            return this.LogLikelihoods(new[] { input }, new[] { result })[0];
        }


        /// <summary>
        /// Computes the log-likelihoods that the given input
        /// vectors belongs to each of the possible classes.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <param name="result">An array where the log-likelihoods will be stored,
        /// avoiding unnecessary memory allocations.</param>
        /// <returns>System.Double[][].</returns>
        /// <exception cref="Exception"></exception>
        public double[][] LogLikelihoods(TInput[] input, double[][] result)
        {
            TInput[] p = new TInput[1];
            double[] o = new double[1];

            for (int i = 0; i < input.Length; i++)
            {
                p[0] = input[i];
                double d = this.LogLikelihood(p, result: o)[0];
                result[i][CLASS_NEGATIVE] = Math.Log(1.0 - Math.Exp(d));
                result[i][CLASS_POSITIVE] = d;
#if DEBUG
                double sum = result[i].Exp().Sum();
                if (!sum.IsEqual(1, rtol: 1e-5))
                    throw new Exception();
#endif
            }

            return result;
        }



        // Input, decision

        /// <summary>
        ///   Predicts a class label vector for the given input vector, returning the
        ///   log-likelihood that the input vector belongs to its predicted class.
        /// </summary>
        ///
        /// <param name="input">The input vector.</param>
        /// <param name="decision">The class label predicted by the classifier.</param>
        ///
        public double LogLikelihood(TInput input, out bool decision)
        {
            bool[] d = new bool[1];
            double s = this.LogLikelihood(new TInput[] { input }, ref d)[0];
            decision = d[0];
            return s;
        }

        /// <summary>
        /// Predicts a class label vector for the given input vector, returning the
        /// log-likelihood that the input vector belongs to its predicted class.
        /// </summary>
        /// <param name="input">The input vector.</param>
        /// <param name="decision">The class label predicted by the classifier.</param>
        /// <returns>System.Double.</returns>
        public double LogLikelihood(TInput input, out double decision)
        {
            double[] d = new double[1];
            double s = this.LogLikelihood(new TInput[] { input }, ref d)[0];
            decision = d[0];
            return s;
        }

        /// <summary>
        /// Predicts a class label vector for the given input vector, returning the
        /// log-likelihood that the input vector belongs to its predicted class.
        /// </summary>
        /// <param name="input">The input vector.</param>
        /// <param name="decision">The class label predicted by the classifier.</param>
        /// <returns>System.Double.</returns>
        public double LogLikelihood(TInput input, out int decision)
        {
            int[] d = new int[1];
            double s = this.LogLikelihood(new TInput[] { input }, ref d)[0];
            decision = d[0];
            return s;
        }



        /// <summary>
        /// Predicts a class label vector for the given input vector, returning the
        /// log-likelihoods of the input vector belonging to each possible class.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <param name="decision">The class labels associated with each input
        /// vector, as predicted by the classifier. If passed as null, the classifier
        /// will create a new array.</param>
        /// <returns>System.Double[].</returns>
        public double[] LogLikelihoods(TInput input, out bool decision)
        {
            return this.LogLikelihoods(input, out decision, new double[this.NumberOfClasses]);
        }

        /// <summary>
        /// Predicts a class label vector for the given input vector, returning the
        /// log-likelihoods of the input vector belonging to each possible class.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <param name="decision">The class labels associated with each input
        /// vector, as predicted by the classifier. If passed as null, the classifier
        /// will create a new array.</param>
        /// <returns>System.Double[].</returns>
        public double[] LogLikelihoods(TInput input, out int decision)
        {
            return this.ToMultilabel().LogLikelihoods(input, out decision, new double[this.NumberOfClasses]);
        }

        /// <summary>
        /// Predicts a class label vector for the given input vector, returning the
        /// log-likelihoods of the input vector belonging to each possible class.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <param name="decision">The class labels associated with each input
        /// vector, as predicted by the classifier. If passed as null, the classifier
        /// will create a new array.</param>
        /// <returns>System.Double[].</returns>
        public double[] LogLikelihoods(TInput input, out double decision)
        {
            return this.ToMultilabel().LogLikelihoods(input, out decision, new double[this.NumberOfClasses]);
        }

        /// <summary>
        /// Predicts a class label vector for the given input vector, returning the
        /// log-likelihoods of the input vector belonging to each possible class.
        /// </summary>
        /// <param name="input">The input vector.</param>
        /// <param name="decision">The class label predicted by the classifier.</param>
        public double[] LogLikelihoods(TInput input, ref bool[] decision)
        {
            return this.LogLikelihoods(input, ref decision, new double[this.NumberOfClasses]);
        }

        /// <summary>
        /// Predicts a class label vector for the given input vector, returning the
        /// log-likelihoods of the input vector belonging to each possible class.
        /// </summary>
        /// <param name="input">The input vector.</param>
        /// <param name="decision">The class label predicted by the classifier.</param>
        /// <returns>System.Double[].</returns>
        public double[] LogLikelihoods(TInput input, ref int[] decision)
        {
            return this.ToMultilabel().LogLikelihoods(input, ref decision, new double[this.NumberOfClasses]);
        }

        /// <summary>
        /// Predicts a class label vector for the given input vector, returning the
        /// log-likelihoods of the input vector belonging to each possible class.
        /// </summary>
        /// <param name="input">The input vector.</param>
        /// <param name="decision">The class label predicted by the classifier.</param>
        /// <returns>System.Double[].</returns>
        public double[] LogLikelihoods(TInput input, ref double[] decision)
        {
            return this.ToMultilabel().LogLikelihoods(input, ref decision, new double[this.NumberOfClasses]);
        }





        // Input, decision, result


        /// <summary>
        /// Predicts a class label vector for the given input vector, returning the
        /// log-likelihoods of the input vector belonging to each possible class.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <param name="decision">The class labels associated with each input
        /// vector, as predicted by the classifier. If passed as null, the classifier
        /// will create a new array.</param>
        /// <param name="result">An array where the log-likelihoods will be stored,
        ///   avoiding unnecessary memory allocations.</param>
        public double[] LogLikelihoods(TInput input, out bool decision, double[] result)
        {
            bool[] d = new bool[1];
            double[][] r = new[] { result };
            double[] s = this.LogLikelihoods(new TInput[] { input }, ref d, r)[0];
            decision = d[0];
            return s;
        }


        /// <summary>
        /// Predicts a class label vector for the given input vector, returning the
        /// log-likelihoods of the input vector belonging to each possible class.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <param name="decision">The class labels associated with each input
        /// vector, as predicted by the classifier. If passed as null, the classifier
        /// will create a new array.</param>
        /// <param name="result">An array where the distances will be stored,
        /// avoiding unnecessary memory allocations.</param>
        /// <returns>System.Double[].</returns>
        public double[] LogLikelihoods(TInput input, out double decision, double[] result)
        {
            double[] d = new double[1];
            double[][] r = new[] { result };
            double[] s = this.LogLikelihoods(new TInput[] { input }, ref d, r)[0];
            decision = d[0];
            return s;
        }


        /// <summary>
        /// Predicts a class label vector for the given input vector, returning the
        /// log-likelihoods of the input vector belonging to each possible class.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <param name="decision">The class labels associated with each input
        /// vector, as predicted by the classifier. If passed as null, the classifier
        /// will create a new array.</param>
        /// <param name="result">An array where the distances will be stored,
        /// avoiding unnecessary memory allocations.</param>
        /// <returns>System.Double[].</returns>
        public double[] LogLikelihoods(TInput input, out int decision, double[] result)
        {
            int[] d = new int[1];
            double[][] r = new[] { result };
            double[] s = this.LogLikelihoods(new TInput[] { input }, ref d, r)[0];
            decision = d[0];
            return s;
        }

        /// <summary>
        /// Predicts a class label vector for the given input vector, returning the
        /// log-likelihoods of the input vector belonging to each possible class.
        /// </summary>
        /// <param name="input">The input vector.</param>
        /// <param name="decision">The class label predicted by the classifier.</param>
        /// <param name="result">An array where the log-likelihoods will be stored,
        ///   avoiding unnecessary memory allocations.</param>
        public double[] LogLikelihoods(TInput input, ref bool[] decision, double[] result)
        {
            bool[] d = new bool[1];
            double[][] r = new[] { result };
            double[] s = this.LogLikelihoods(new TInput[] { input }, ref d, r)[0];
            decision = Vector.OneHot(d[0], this.createOrReuse(input, decision));
            return s;
        }

        /// <summary>
        /// Predicts a class label vector for the given input vector, returning the
        /// log-likelihoods of the input vector belonging to each possible class.
        /// </summary>
        /// <param name="input">The input vector.</param>
        /// <param name="decision">The class label predicted by the classifier.</param>
        /// <param name="result">An array where the distances will be stored,
        /// avoiding unnecessary memory allocations.</param>
        /// <returns>System.Double[].</returns>
        public double[] LogLikelihoods(TInput input, ref int[] decision, double[] result)
        {
            bool[] d = new bool[1];
            double[][] r = new[] { result };
            double[] s = this.LogLikelihoods(new TInput[] { input }, ref d, r)[0];
            decision = Vector.OneHot(d[0], this.createOrReuse(input, decision));
            return s;
        }

        /// <summary>
        /// Predicts a class label vector for the given input vector, returning the
        /// log-likelihoods of the input vector belonging to each possible class.
        /// </summary>
        /// <param name="input">The input vector.</param>
        /// <param name="decision">The class label predicted by the classifier.</param>
        /// <param name="result">An array where the distances will be stored,
        /// avoiding unnecessary memory allocations.</param>
        /// <returns>System.Double[].</returns>
        public double[] LogLikelihoods(TInput input, ref double[] decision, double[] result)
        {
            bool[] d = new bool[1];
            double[][] r = new[] { result };
            double[] s = this.LogLikelihoods(new TInput[] { input }, ref d, r)[0];
            decision = Vector.OneHot(d[0], this.createOrReuse(input, decision));
            return s;
        }


        // Input[], decision[]

        /// <summary>
        /// Predicts a class label for each input vector, returning the
        /// log-likelihood that each vector belongs to its predicted class.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <param name="decision">The class labels associated with each input
        /// vector, as predicted by the classifier. If passed as null, the classifier
        /// will create a new array.</param>
        /// <returns>System.Double[].</returns>
        public double[] LogLikelihood(TInput[] input, ref int[] decision)
        {
            return this.ToMulticlass().LogLikelihood(input, ref decision, new double[input.Length]);
        }

        /// <summary>
        /// Predicts a class label for each input vector, returning the
        /// log-likelihood that each vector belongs to its predicted class.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <param name="decision">The class labels associated with each input
        /// vector, as predicted by the classifier. If passed as null, the classifier
        /// will create a new array.</param>
        /// <returns>System.Double[].</returns>
        public double[] LogLikelihood(TInput[] input, ref double[] decision)
        {
            return ((IMulticlassLikelihoodClassifier<TInput>)this).LogLikelihood(input, ref decision, new double[input.Length]);
        }

        /// <summary>
        /// Predicts a class label for each input vector, returning the
        /// log-likelihood that each vector belongs to its predicted class.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <param name="decision">The class labels associated with each input
        /// vector, as predicted by the classifier. If passed as null, the classifier
        /// will create a new array.</param>
        public double[] LogLikelihood(TInput[] input, ref bool[] decision)
        {
            return this.LogLikelihood(input, ref decision, new double[input.Length]);
        }







        /// <summary>
        /// Predicts a class label vector for each input vector, returning the
        /// log-likelihoods of the input vector belonging to each possible class.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <param name="decision">The labels predicted by the classifier.</param>
        /// <returns>System.Double[][].</returns>
        public double[][] LogLikelihoods(TInput[] input, ref int[] decision)
        {
            return this.ToMultilabel().LogLikelihoods(input, ref decision, this.create<double>(input));
        }


        /// <summary>
        /// Predicts a class label vector for each input vector, returning the
        /// log-likelihoods of the input vector belonging to each possible class.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <param name="decision">The labels predicted by the classifier.</param>
        /// <returns>System.Double[][].</returns>
        public double[][] LogLikelihoods(TInput[] input, ref double[] decision)
        {
            return this.ToMultilabel().LogLikelihoods(input, ref decision, this.create<double>(input));
        }


        /// <summary>
        /// Predicts a class label vector for each input vector, returning the
        /// log-likelihoods of the input vector belonging to each possible class.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <param name="decision">The labels predicted by the classifier.</param>
        public double[][] LogLikelihoods(TInput[] input, ref bool[] decision)
        {
            return this.LogLikelihoods(input, ref decision, this.create<double>(input));
        }


        double[][] IMultilabelLikelihoodClassifierBase<TInput, int[]>.LogLikelihoods(TInput[] input, ref int[][] decision)
        {
            return this.ToMultilabel().LogLikelihoods(input, ref decision, this.create<double>(input));
        }


        /// <summary>
        /// Predicts a class label vector for each input vector, returning the
        /// log-likelihoods of the input vector belonging to each possible class.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <param name="decision">The labels predicted by the classifier.</param>
        public double[][] LogLikelihoods(TInput[] input, ref bool[][] decision)
        {
            return this.ToMultilabel().LogLikelihoods(input, ref decision, this.create<double>(input));
        }


        /// <summary>
        /// Predicts a class label vector for each input vector, returning the
        /// log-likelihoods of the input vector belonging to each possible class.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <param name="decision">The labels predicted by the classifier.</param>
        /// <returns>System.Double[][].</returns>
        public double[][] LogLikelihoods(TInput[] input, ref double[][] decision)
        {
            return this.ToMultilabel().LogLikelihoods(input, ref decision, this.create<double>(input));
        }




        // Input[], decision[], result[]        

        /// <summary>
        /// Predicts a class label for each input vector, returning the
        /// log-likelihood that each vector belongs to its predicted class.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <param name="decision">The class labels associated with each input
        /// vector, as predicted by the classifier. If passed as null, the classifier
        /// will create a new array.</param>
        /// <param name="result">An array where the log-likelihoods will be stored,
        ///   avoiding unnecessary memory allocations.</param>
        public double[] LogLikelihood(TInput[] input, ref bool[] decision, double[] result)
        {
            this.LogLikelihood(input, result);
            decision = this.createOrReuse(input, decision);
            for (int i = 0; i < input.Length; i++)
                decision[i] = Classes.Decide(result[i] - this.LOGLIKELIHOOD_DECISION_THRESHOLD);
            return result;
        }

        /// <summary>
        /// Predicts a class label for each input vector, returning the
        /// log-likelihood that each vector belongs to its predicted class.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <param name="decision">The class labels associated with each input
        /// vector, as predicted by the classifier. If passed as null, the classifier
        /// will create a new array.</param>
        /// <param name="result">An array where the log-likelihoods will be stored,
        /// avoiding unnecessary memory allocations.</param>
        /// <returns>System.Double[].</returns>
        public double[] LogLikelihood(TInput[] input, ref int[] decision, double[] result)
        {
            this.LogLikelihood(input, result);
            decision = this.createOrReuse(input, decision);
            for (int i = 0; i < input.Length; i++)
                decision[i] = Classes.ToZeroOne(result[i] - this.LOGLIKELIHOOD_DECISION_THRESHOLD);
            return result;
        }

        /// <summary>
        /// Predicts a class label for each input vector, returning the
        /// log-likelihood that each vector belongs to its predicted class.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <param name="decision">The class labels associated with each input
        /// vector, as predicted by the classifier. If passed as null, the classifier
        /// will create a new array.</param>
        /// <param name="result">An array where the log-likelihoods will be stored,
        /// avoiding unnecessary memory allocations.</param>
        /// <returns>System.Double[].</returns>
        public double[] LogLikelihood(TInput[] input, ref double[] decision, double[] result)
        {
            this.LogLikelihood(input, result);
            decision = this.createOrReuse(input, decision);
            for (int i = 0; i < input.Length; i++)
                decision[i] = Classes.ToZeroOne(result[i] - this.LOGLIKELIHOOD_DECISION_THRESHOLD);
            return result;
        }

        /// <summary>
        /// Predicts a class label vector for each input vector, returning the
        /// log-likelihoods of the input vector belonging to each possible class.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <param name="decision">The labels predicted by the classifier.</param>
        /// <param name="result">An array where the log-likelihoods will be stored,
        ///   avoiding unnecessary memory allocations.</param>
        public double[][] LogLikelihoods(TInput[] input, ref bool[][] decision, double[][] result)
        {
            this.ToMultilabel().LogLikelihoods(input, result);
            decision = this.createOrReuse(input, decision);
            for (int i = 0; i < result.Length; i++)
                Vector.OneHot(Classes.ToZeroOne(result[i][CLASS_POSITIVE] - this.LOGLIKELIHOOD_DECISION_THRESHOLD), result: decision[i]);
            return result;
        }

        /// <summary>
        /// Predicts a class label vector for each input vector, returning the
        /// log-likelihoods of the input vector belonging to each possible class.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <param name="decision">The labels predicted by the classifier.</param>
        /// <param name="result">An array where the log-likelihoods will be stored,
        /// avoiding unnecessary memory allocations.</param>
        /// <returns>System.Double[][].</returns>
        public double[][] LogLikelihoods(TInput[] input, ref int[][] decision, double[][] result)
        {
            this.ToMultilabel().LogLikelihoods(input, result);
            decision = this.createOrReuse(input, decision);
            for (int i = 0; i < result.Length; i++)
                Vector.OneHot(Classes.ToZeroOne(result[i][CLASS_POSITIVE] - this.LOGLIKELIHOOD_DECISION_THRESHOLD), result: decision[i]);
            return result;
        }

        /// <summary>
        /// Predicts a class label vector for each input vector, returning the
        /// log-likelihoods of the input vector belonging to each possible class.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <param name="decision">The labels predicted by the classifier.</param>
        /// <param name="result">An array where the log-likelihoods will be stored,
        /// avoiding unnecessary memory allocations.</param>
        /// <returns>System.Double[][].</returns>
        public double[][] LogLikelihoods(TInput[] input, ref double[][] decision, double[][] result)
        {
            this.ToMultilabel().LogLikelihoods(input, result);
            decision = this.createOrReuse(input, decision);
            for (int i = 0; i < result.Length; i++)
                Vector.OneHot(Classes.ToZeroOne(result[i][CLASS_POSITIVE] - this.LOGLIKELIHOOD_DECISION_THRESHOLD), result: decision[i]);
            return result;
        }


        /// <summary>
        /// Predicts a class label vector for each input vector, returning the
        /// log-likelihoods of the input vector belonging to each possible class.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <param name="decision">The labels predicted by the classifier.</param>
        /// <param name="result">An array where the log-likelihoods will be stored,
        ///   avoiding unnecessary memory allocations.</param>
        public double[][] LogLikelihoods(TInput[] input, ref bool[] decision, double[][] result)
        {
            this.ToMultilabel().LogLikelihoods(input, result);
            decision = this.createOrReuse(input, decision);
            for (int i = 0; i < result.Length; i++)
                decision[i] = Classes.Decide(result[i][CLASS_POSITIVE] - this.LOGLIKELIHOOD_DECISION_THRESHOLD);
            return result;
        }




        /// <summary>
        /// Predicts a class label vector for each input vector, returning the
        /// log-likelihoods of the input vector belonging to each possible class.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <param name="decision">The labels predicted by the classifier.</param>
        /// <param name="result">An array where the log-likelihoods will be stored,
        /// avoiding unnecessary memory allocations.</param>
        /// <returns>System.Double[][].</returns>
        public double[][] LogLikelihoods(TInput[] input, ref int[] decision, double[][] result)
        {
            this.ToMultilabel().LogLikelihoods(input, result);
            decision = this.createOrReuse(input, decision);
            for (int i = 0; i < result.Length; i++)
                decision[i] = Classes.ToZeroOne(result[i][CLASS_POSITIVE] - this.LOGLIKELIHOOD_DECISION_THRESHOLD);
            return result;
        }


        /// <summary>
        /// Predicts a class label vector for each input vector, returning the
        /// log-likelihoods of the input vector belonging to each possible class.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <param name="decision">The labels predicted by the classifier.</param>
        /// <param name="result">An array where the log-likelihoods will be stored,
        /// avoiding unnecessary memory allocations.</param>
        /// <returns>System.Double[][].</returns>
        public double[][] LogLikelihoods(TInput[] input, ref double[] decision, double[][] result)
        {
            this.ToMultilabel().LogLikelihoods(input, result);
            decision = this.createOrReuse(input, decision);
            for (int i = 0; i < result.Length; i++)
                decision[i] = Classes.ToZeroOne(result[i][CLASS_POSITIVE] - this.LOGLIKELIHOOD_DECISION_THRESHOLD);
            return result;
        }

        #endregion





        #region Probability


        /// <summary>
        /// Predicts a class label for the given input vector, returning the
        /// probability that the input vector belongs to its predicted class.
        /// </summary>
        /// <param name="input">The input vector.</param>
        public double Probability(TInput input)
        {
            return this.Probability(new TInput[] { input })[0];
        }


        // Input, classIndex
        double IMultilabelLikelihoodClassifier<TInput>.Probability(TInput input, int classIndex)
        {
            return this.ToMultilabel().Probabilities(input)[classIndex];
        }

        double[] IMultilabelLikelihoodClassifier<TInput>.Probability(TInput[] input, int classIndex)
        {
            return ((IMulticlassLikelihoodClassifier<TInput>)this).Probability(input, classIndex, new double[input.Length]);
        }


        double[] IMultilabelLikelihoodClassifier<TInput>.Probability(TInput[] input, int classIndex, double[] result)
        {
            var temp = new double[this.NumberOfClasses];
            for (int i = 0; i < input.Length; i++)
                result[i] = this.ToMultilabel().Probabilities(input[i], result: temp)[classIndex];
            return result;
        }


        // Input, classIndex[]
        double[] IMultilabelLikelihoodClassifier<TInput>.Probability(TInput[] input, int[] classIndex)
        {
            return ((IMulticlassLikelihoodClassifier<TInput>)this).Probability(input, classIndex, new double[input.Length]);
        }

        double[] IMultilabelLikelihoodClassifier<TInput>.Probability(TInput[] input, int[] classIndex, double[] result)
        {
            var temp = new double[this.NumberOfClasses];
            for (int i = 0; i < input.Length; i++)
                result[i] = this.ToMultilabel().Probabilities(input[i], result: temp)[classIndex[i]];
            return result;
        }


        // Input

        /// <summary>
        /// Computes the probabilities that the given input
        /// vector belongs to each of the possible classes.
        /// </summary>
        /// <param name="input">The input vector.</param>
        /// <returns>System.Double[].</returns>
        public double[] Probabilities(TInput input)
        {
            return this.ToMultilabel().Probabilities(input, new double[this.NumberOfClasses]);
        }

        /// <summary>
        /// Predicts a class label for the given input vector, returning the
        /// probability that the input vector belongs to its predicted class.
        /// </summary>
        /// <param name="input">The input vector.</param>
        public double[] Probability(TInput[] input)
        {
            return this.Probability(input, new double[input.Length]);
        }

        /// <summary>
        /// Computes the probabilities that the given input
        /// vectors belongs to each of the possible classes.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <returns>System.Double[][].</returns>
        public double[][] Probabilities(TInput[] input)
        {
            return this.ToMultilabel().Probabilities(input, this.create<double>(input));
        }




        // Input, result

        /// <summary>
        /// Computes the probabilities that the given input
        /// vector belongs to each of the possible classes.
        /// </summary>
        /// <param name="input">The input vector.</param>
        /// <param name="result">An array where the probabilities will be stored,
        /// avoiding unnecessary memory allocations.</param>
        /// <returns>System.Double[].</returns>
        public double[] Probabilities(TInput input, double[] result)
        {
            return this.Probabilities(new TInput[] { input }, result: new[] { result })[0];
        }

        /// <summary>
        /// Predicts a class label for the given input vector, returning the
        /// probability that the input vector belongs to its predicted class.
        /// </summary>
        /// <param name="input">The input vector.</param>
        /// <param name="result">An array where the probabilities will be stored,
        ///   avoiding unnecessary memory allocations.</param>
        public virtual double[] Probability(TInput[] input, double[] result)
        {
            this.LogLikelihood(input, result: result);
            Elementwise.Exp(result, result: result);
            return result;
        }

        /// <summary>
        /// Computes the probabilities that the given input
        /// vectors belongs to each of the possible classes.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <param name="result">An array where the probabilities will be stored,
        /// avoiding unnecessary memory allocations.</param>
        /// <returns>System.Double[][].</returns>
        public double[][] Probabilities(TInput[] input, double[][] result)
        {
            double[] r = this.Probability(input, new double[input.Length]);

            for (int i = 0; i < input.Length; i++)
            {
                result[i][CLASS_NEGATIVE] = 1.0 - r[i];
                result[i][CLASS_POSITIVE] = r[i];
            }

            return result;
        }



        // Input, decision


        /// <summary>
        ///   Predicts a class label for the given input vector, returning the
        ///   probability that the input vector belongs to its predicted class.
        /// </summary>
        /// 
        /// <param name="input">The input vector.</param>
        /// <param name="decision">The class label predicted by the classifier.</param>
        ///
        public double Probability(TInput input, out bool decision)
        {
            bool[] d = new bool[1];
            double s = this.Probability(new TInput[] { input }, ref d)[0];
            decision = d[0];
            return s;
        }

        /// <summary>
        /// Predicts a class label for the given input vector, returning the
        /// probability that the input vector belongs to its predicted class.
        /// </summary>
        /// <param name="input">The input vector.</param>
        /// <param name="decision">The class label predicted by the classifier.</param>
        /// <returns>System.Double.</returns>
        public double Probability(TInput input, out double decision)
        {
            double[] d = new double[1];
            double s = this.Probability(new TInput[] { input }, ref d)[0];
            decision = d[0];
            return s;
        }

        /// <summary>
        /// Predicts a class label for the given input vector, returning the
        /// probability that the input vector belongs to its predicted class.
        /// </summary>
        /// <param name="input">The input vector.</param>
        /// <param name="decision">The class label predicted by the classifier.</param>
        /// <returns>System.Double.</returns>
        public double Probability(TInput input, out int decision)
        {
            int[] d = new int[1];
            double s = this.Probability(new TInput[] { input }, ref d)[0];
            decision = d[0];
            return s;
        }



        /// <summary>
        /// Predicts a class label vector for the given input vector, returning the
        /// probabilities of the input vector belonging to each possible class.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <param name="decision">The class labels associated with each input
        /// vector, as predicted by the classifier. If passed as null, the classifier
        /// will create a new array.</param>
        public double[] Probabilities(TInput input, out bool decision)
        {
            return this.Probabilities(input, out decision, new double[this.NumberOfClasses]);
        }

        /// <summary>
        /// Predicts a class label vector for the given input vector, returning the
        /// probabilities of the input vector belonging to each possible class.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <param name="decision">The class labels associated with each input
        /// vector, as predicted by the classifier. If passed as null, the classifier
        /// will create a new array.</param>
        /// <returns>System.Double[].</returns>
        public double[] Probabilities(TInput input, out int decision)
        {
            return this.ToMultilabel().Probabilities(input, out decision, new double[this.NumberOfClasses]);
        }

        /// <summary>
        /// Predicts a class label vector for the given input vector, returning the
        /// probabilities of the input vector belonging to each possible class.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <param name="decision">The class labels associated with each input
        /// vector, as predicted by the classifier. If passed as null, the classifier
        /// will create a new array.</param>
        /// <returns>System.Double[].</returns>
        public double[] Probabilities(TInput input, out double decision)
        {
            return this.ToMultilabel().Probabilities(input, out decision, new double[this.NumberOfClasses]);
        }

        /// <summary>
        /// Predicts a class label vector for the given input vector, returning the
        /// probabilities of the input vector belonging to each possible class.
        /// </summary>
        /// <param name="input">The input vector.</param>
        /// <param name="decision">The class label predicted by the classifier.</param>
        public double[] Probabilities(TInput input, ref bool[] decision)
        {
            return this.ToMultilabel().Probabilities(input, ref decision, new double[this.NumberOfClasses]);
        }

        /// <summary>
        /// Predicts a class label vector for the given input vector, returning the
        /// probabilities of the input vector belonging to each possible class.
        /// </summary>
        /// <param name="input">The input vector.</param>
        /// <param name="decision">The class label predicted by the classifier.</param>
        /// <returns>System.Double[].</returns>
        public double[] Probabilities(TInput input, ref int[] decision)
        {
            return this.ToMultilabel().Probabilities(input, ref decision, new double[this.NumberOfClasses]);
        }

        /// <summary>
        /// Probabilitieses the specified input.
        /// </summary>
        /// <param name="input">The input.</param>
        /// <param name="decision">The decision.</param>
        /// <returns>System.Double[].</returns>
        public double[] Probabilities(TInput input, ref double[] decision)
        {
            return this.ToMultilabel().Probabilities(input, ref decision, new double[this.NumberOfClasses]);
        }




        // Input, decision, result

        /// <summary>
        /// Predicts a class label vector for the given input vector, returning the
        /// probabilities of the input vector belonging to each possible class.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <param name="decision">The class labels associated with each input
        /// vector, as predicted by the classifier. If passed as null, the classifier
        /// will create a new array.</param>
        /// <param name="result">An array where the probabilities will be stored,
        ///   avoiding unnecessary memory allocations.</param>
        public double[] Probabilities(TInput input, out bool decision, double[] result)
        {
            bool[] d = new bool[1];
            double[][] r = new[] { result };
            double[] s = this.LogLikelihoods(new TInput[] { input }, ref d, r)[0];
            decision = d[0];
            return s;
        }

        /// <summary>
        /// Predicts a class label vector for the given input vector, returning the
        /// probabilities of the input vector belonging to each possible class.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <param name="decision">The class labels associated with each input
        /// vector, as predicted by the classifier. If passed as null, the classifier
        /// will create a new array.</param>
        /// <param name="result">An array where the distances will be stored,
        /// avoiding unnecessary memory allocations.</param>
        /// <returns>System.Double[].</returns>
        public double[] Probabilities(TInput input, out double decision, double[] result)
        {
            double[] d = new double[1];
            double[][] r = new[] { result };
            double[] s = this.LogLikelihoods(new TInput[] { input }, ref d, r)[0];
            decision = d[0];
            return s;
        }

        /// <summary>
        /// Predicts a class label vector for the given input vector, returning the
        /// probabilities of the input vector belonging to each possible class.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <param name="decision">The class labels associated with each input
        /// vector, as predicted by the classifier. If passed as null, the classifier
        /// will create a new array.</param>
        /// <param name="result">An array where the distances will be stored,
        /// avoiding unnecessary memory allocations.</param>
        /// <returns>System.Double[].</returns>
        public double[] Probabilities(TInput input, out int decision, double[] result)
        {
            int[] d = new int[1];
            double[][] r = new[] { result };
            double[] s = this.LogLikelihoods(new TInput[] { input }, ref d, r)[0];
            decision = d[0];
            return s;
        }

        /// <summary>
        /// Predicts a class label vector for the given input vector, returning the
        /// probabilities of the input vector belonging to each possible class.
        /// </summary>
        /// <param name="input">The input vector.</param>
        /// <param name="decision">The class label predicted by the classifier.</param>
        /// <param name="result">An array where the distances will be stored,
        /// avoiding unnecessary memory allocations.</param>
        /// <returns>System.Double[].</returns>
        public double[] Probabilities(TInput input, ref bool[] decision, double[] result)
        {
            bool[] d = new bool[1];
            double[][] r = new[] { result };
            double[] s = this.LogLikelihoods(new TInput[] { input }, ref d, r)[0];
            decision = Vector.OneHot(d[0], this.createOrReuse(input, decision));
            return s;
        }

        /// <summary>
        /// Predicts a class label vector for the given input vector, returning the
        /// probabilities of the input vector belonging to each possible class.
        /// </summary>
        /// <param name="input">The input vector.</param>
        /// <param name="decision">The class label predicted by the classifier.</param>
        /// <param name="result">An array where the distances will be stored,
        /// avoiding unnecessary memory allocations.</param>
        /// <returns>System.Double[].</returns>
        public double[] Probabilities(TInput input, ref int[] decision, double[] result)
        {
            bool[] d = new bool[1];
            double[][] r = new[] { result };
            double[] s = this.LogLikelihoods(new TInput[] { input }, ref d, r)[0];
            decision = Vector.OneHot(d[0], this.createOrReuse(input, decision));
            return s;
        }

        /// <summary>
        /// Predicts a class label vector for the given input vector, returning the
        /// probabilities of the input vector belonging to each possible class.
        /// </summary>
        /// <param name="input">The input vector.</param>
        /// <param name="decision">The class label predicted by the classifier.</param>
        /// <param name="result">An array where the distances will be stored,
        /// avoiding unnecessary memory allocations.</param>
        /// <returns>System.Double[].</returns>
        public double[] Probabilities(TInput input, ref double[] decision, double[] result)
        {
            bool[] d = new bool[1];
            double[][] r = new[] { result };
            double[] s = this.LogLikelihoods(new TInput[] { input }, ref d, r)[0];
            decision = Vector.OneHot(d[0], this.createOrReuse(input, decision));
            return s;
        }


        // Input[], decision[]

        /// <summary>
        /// Predicts a class label for each input vector, returning the
        /// probability that each vector belongs to its predicted class.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <param name="decision">The class labels associated with each input
        /// vector, as predicted by the classifier. If passed as null, the classifier
        /// will create a new array.</param>
        /// <returns>System.Double[].</returns>
        public double[] Probability(TInput[] input, ref int[] decision)
        {
            return this.ToMulticlass().Probability(input, ref decision, new double[input.Length]);
        }

        /// <summary>
        /// Predicts a class label for each input vector, returning the
        /// probability that each vector belongs to its predicted class.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <param name="decision">The class labels associated with each input
        /// vector, as predicted by the classifier. If passed as null, the classifier
        /// will create a new array.</param>
        /// <returns>System.Double[].</returns>
        public double[] Probability(TInput[] input, ref double[] decision)
        {
            return ((IMulticlassLikelihoodClassifier<TInput>)this).Probability(input, ref decision, new double[input.Length]);
        }

        /// <summary>
        /// Predicts a class label for each input vector, returning the
        /// probability that each vector belongs to its predicted class.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <param name="decision">The class labels associated with each input
        /// vector, as predicted by the classifier. If passed as null, the classifier
        /// will create a new array.</param>
        /// <returns></returns>
        public double[] Probability(TInput[] input, ref bool[] decision)
        {
            return this.Probability(input, ref decision, new double[input.Length]);
        }






        /// <summary>
        /// Predicts a class label vector for each input vector, returning the
        /// probabilities of the input vector belonging to each possible class.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <param name="decision">The labels predicted by the classifier.</param>
        /// <returns>System.Double[][].</returns>
        public double[][] Probabilities(TInput[] input, ref int[] decision)
        {
            return this.ToMultilabel().Probabilities(input, ref decision, this.create<double>(input));
        }

        /// <summary>
        /// Predicts a class label vector for each input vector, returning the
        /// probabilities of the input vector belonging to each possible class.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <param name="decision">The labels predicted by the classifier.</param>
        /// <returns>System.Double[][].</returns>
        public double[][] Probabilities(TInput[] input, ref double[] decision)
        {
            return this.ToMultilabel().Probabilities(input, ref decision, this.create<double>(input));
        }

        /// <summary>
        /// Predicts a class label vector for each input vector, returning the
        /// probabilities of the input vector belonging to each possible class.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <param name="decision">The labels predicted by the classifier.</param>
        public double[][] Probabilities(TInput[] input, ref bool[] decision)
        {
            return this.Probabilities(input, ref decision, this.create<double>(input));
        }

        /// <summary>
        /// Predicts a class label vector for each input vector, returning the
        /// probabilities of the input vector belonging to each possible class.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <param name="decision">The labels predicted by the classifier.</param>
        /// <returns>System.Double[][].</returns>
        public double[][] Probabilities(TInput[] input, ref int[][] decision)
        {
            return this.ToMultilabel().Probabilities(input, ref decision, this.create<double>(input));
        }

        /// <summary>
        /// Predicts a class label vector for each input vector, returning the
        /// probabilities of the input vector belonging to each possible class.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <param name="decision">The labels predicted by the classifier.</param>
        /// <returns>System.Double[][].</returns>
        public double[][] Probabilities(TInput[] input, ref bool[][] decision)
        {
            return this.ToMultilabel().Probabilities(input, ref decision, this.create<double>(input));
        }

        /// <summary>
        /// Predicts a class label vector for each input vector, returning the
        /// probabilities of the input vector belonging to each possible class.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <param name="decision">The labels predicted by the classifier.</param>
        /// <returns>System.Double[][].</returns>
        public double[][] Probabilities(TInput[] input, ref double[][] decision)
        {
            return this.ToMultilabel().Probabilities(input, ref decision, this.create<double>(input));
        }




        // Input[], decision[], result[]        

        /// <summary>
        /// Predicts a class label for each input vector, returning the
        /// probability that each vector belongs to its predicted class.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <param name="decision">The class labels associated with each input
        /// vector, as predicted by the classifier. If passed as null, the classifier
        /// will create a new array.</param>
        /// <param name="result">An array where the probabilities will be stored,
        ///   avoiding unnecessary memory allocations.</param>
        public double[] Probability(TInput[] input, ref bool[] decision, double[] result)
        {
            this.Probability(input, result);
            decision = this.createOrReuse(input, decision);
            for (int i = 0; i < input.Length; i++)
                decision[i] = Classes.Decide(result[i] - PROBABILITY_DECISION_THRESHOLD);
            return result;
        }

        /// <summary>
        /// Predicts a class label for each input vector, returning the
        /// probability that each vector belongs to its predicted class.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <param name="decision">The class labels associated with each input
        /// vector, as predicted by the classifier. If passed as null, the classifier
        /// will create a new array.</param>
        /// <param name="result">An array where the probabilities will be stored,
        /// avoiding unnecessary memory allocations.</param>
        /// <returns>System.Double[].</returns>
        public double[] Probability(TInput[] input, ref int[] decision, double[] result)
        {
            this.Probability(input, result);
            decision = this.createOrReuse(input, decision);
            for (int i = 0; i < input.Length; i++)
                decision[i] = Classes.ToZeroOne(result[i] - PROBABILITY_DECISION_THRESHOLD);
            return result;
        }

        /// <summary>
        /// Predicts a class label for each input vector, returning the
        /// probability that each vector belongs to its predicted class.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <param name="decision">The class labels associated with each input
        /// vector, as predicted by the classifier. If passed as null, the classifier
        /// will create a new array.</param>
        /// <param name="result">An array where the probabilities will be stored,
        /// avoiding unnecessary memory allocations.</param>
        /// <returns>System.Double[].</returns>
        public double[] Probability(TInput[] input, ref double[] decision, double[] result)
        {
            this.Probability(input, result);
            decision = this.createOrReuse(input, decision);
            for (int i = 0; i < input.Length; i++)
                decision[i] = Classes.ToZeroOne(result[i] - PROBABILITY_DECISION_THRESHOLD);
            return result;
        }

        /// <summary>
        /// Predicts a class label vector for each input vector, returning the
        /// probabilities of the input vector belonging to each possible class.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <param name="decision">The labels predicted by the classifier.</param>
        /// <param name="result">An array where the probabilities will be stored,
        /// avoiding unnecessary memory allocations.</param>
        /// <returns>System.Double[][].</returns>
        public double[][] Probabilities(TInput[] input, ref bool[][] decision, double[][] result)
        {
            this.ToMultilabel().LogLikelihoods(input, result);
            decision = this.createOrReuse(input, decision);
            for (int i = 0; i < result.Length; i++)
                Vector.OneHot(Classes.ToZeroOne(result[i][CLASS_POSITIVE] - PROBABILITY_DECISION_THRESHOLD), result: decision[i]);
            return result;
        }

        /// <summary>
        /// Predicts a class label vector for each input vector, returning the
        /// probabilities of the input vector belonging to each possible class.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <param name="decision">The labels predicted by the classifier.</param>
        /// <param name="result">An array where the probabilities will be stored,
        /// avoiding unnecessary memory allocations.</param>
        /// <returns>System.Double[][].</returns>
        public double[][] Probabilities(TInput[] input, ref int[][] decision, double[][] result)
        {
            this.ToMultilabel().Probabilities(input, result);
            decision = this.createOrReuse(input, decision);
            for (int i = 0; i < result.Length; i++)
                Vector.OneHot(Classes.ToZeroOne(result[i][CLASS_POSITIVE] - PROBABILITY_DECISION_THRESHOLD), result: decision[i]);
            return result;
        }

        /// <summary>
        /// Predicts a class label vector for each input vector, returning the
        /// probabilities of the input vector belonging to each possible class.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <param name="decision">The labels predicted by the classifier.</param>
        /// <param name="result">An array where the probabilities will be stored,
        /// avoiding unnecessary memory allocations.</param>
        /// <returns>System.Double[][].</returns>
        public double[][] Probabilities(TInput[] input, ref double[][] decision, double[][] result)
        {
            this.ToMultilabel().Probabilities(input, result);
            decision = this.createOrReuse(input, decision);
            for (int i = 0; i < result.Length; i++)
                Vector.OneHot(Classes.ToZeroOne(result[i][CLASS_POSITIVE] - PROBABILITY_DECISION_THRESHOLD), result: decision[i]);
            return result;
        }


        /// <summary>
        /// Predicts a class label vector for each input vector, returning the
        /// probabilities of the input vector belonging to each possible class.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <param name="decision">The labels predicted by the classifier.</param>
        /// <param name="result">An array where the probabilities will be stored,
        ///   avoiding unnecessary memory allocations.</param>
        public double[][] Probabilities(TInput[] input, ref bool[] decision, double[][] result)
        {
            this.ToMultilabel().Probabilities(input, result);
            decision = this.createOrReuse(input, decision);
            for (int i = 0; i < result.Length; i++)
                decision[i] = Classes.Decide(result[i][CLASS_POSITIVE] - PROBABILITY_DECISION_THRESHOLD);
            return result;
        }

        /// <summary>
        /// Predicts a class label vector for each input vector, returning the
        /// probabilities of the input vector belonging to each possible class.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <param name="decision">The labels predicted by the classifier.</param>
        /// <param name="result">An array where the probabilities will be stored,
        /// avoiding unnecessary memory allocations.</param>
        /// <returns>System.Double[][].</returns>
        public double[][] Probabilities(TInput[] input, ref int[] decision, double[][] result)
        {
            this.ToMultilabel().Probabilities(input, result);
            decision = this.createOrReuse(input, decision);
            for (int i = 0; i < result.Length; i++)
                decision[i] = Classes.ToZeroOne(result[i][CLASS_POSITIVE] - PROBABILITY_DECISION_THRESHOLD);
            return result;
        }

        /// <summary>
        /// Predicts a class label vector for each input vector, returning the
        /// probabilities of the input vector belonging to each possible class.
        /// </summary>
        /// <param name="input">A set of input vectors.</param>
        /// <param name="decision">The labels predicted by the classifier.</param>
        /// <param name="result">An array where the probabilities will be stored,
        /// avoiding unnecessary memory allocations.</param>
        /// <returns>System.Double[][].</returns>
        public double[][] Probabilities(TInput[] input, ref double[] decision, double[][] result)
        {
            this.ToMultilabel().Probabilities(input, result);
            decision = this.createOrReuse(input, decision);
            for (int i = 0; i < result.Length; i++)
                decision[i] = Classes.ToZeroOne(result[i][CLASS_POSITIVE] - PROBABILITY_DECISION_THRESHOLD);
            return result;
        }

        #endregion





        // Transform

        double ICovariantTransform<TInput, double>.Transform(TInput input)
        {
            return this.Probability(input);
        }

        /// <summary>
        /// Applies the transformation to a set of input vectors,
        /// producing an associated set of output vectors.
        /// </summary>
        /// <param name="input">The input data to which
        /// the transformation should be applied.</param>
        /// <param name="result">The location to where to store the
        /// result of this transformation.</param>
        /// <returns>
        /// The output generated by applying this
        /// transformation to the given input.
        /// </returns>
        public override double[] Transform(TInput input, double[] result)
        {
            return this.ToMultilabel().Probabilities(input, result);
        }

        /// <summary>
        /// Applies the transformation to a set of input vectors,
        /// producing an associated set of output vectors.
        /// </summary>
        /// <param name="input">The input data to which
        /// the transformation should be applied.</param>
        /// <param name="result">The location to where to store the
        /// result of this transformation.</param>
        /// <returns>
        /// The output generated by applying this
        /// transformation to the given input.
        /// </returns>
        public override double[][] Transform(TInput[] input, double[][] result)
        {
            return this.ToMultilabel().Probabilities(input, result);
        }

        /// <summary>
        /// Applies the transformation to a set of input vectors,
        /// producing an associated set of output vectors.
        /// </summary>
        /// <param name="input">The input data to which
        /// the transformation should be applied.</param>
        /// <param name="result">The location to where to store the
        /// result of this transformation.</param>
        /// <returns>
        /// The output generated by applying this
        /// transformation to the given input.
        /// </returns>
        public override double[] Transform(TInput[] input, double[] result)
        {
            return this.Probability(input, result);
        }


        int IMulticlassLikelihoodClassifier<TInput>.Decide(TInput input)
        {
            return ((IClassifier<TInput, int>)this).Decide(input);
        }

        int[] IMulticlassLikelihoodClassifier<TInput>.Decide(TInput[] input)
        {
            return ((IClassifier<TInput, int>)this).Decide(input);
        }



        /// <summary>
        /// Views this instance as a multi-class generative classifier,
        /// giving access to more advanced methods, such as the prediction
        /// of integer labels.
        /// </summary>
        /// <returns>
        /// This instance seen as an <see cref="IMulticlassLikelihoodClassifier{TInput}" />.
        /// </returns>
        new public IMulticlassLikelihoodClassifier<TInput> ToMulticlass()
        {
            return (IMulticlassLikelihoodClassifier<TInput>)this;
        }

        /// <summary>
        /// Views this instance as a multi-class generative classifier,
        /// giving access to more advanced methods, such as the prediction
        /// of integer labels.
        /// </summary>
        /// <returns>
        /// This instance seen as an <see cref="IMulticlassLikelihoodClassifier{TInput}" />.
        /// </returns>
        new public IMulticlassLikelihoodClassifier<TInput, T> ToMulticlass<T>()
        {
            return (IMulticlassLikelihoodClassifier<TInput, T>)this;
        }

        /// <summary>
        /// Views this instance as a multi-label generative classifier,
        /// giving access to more advanced methods, such as the prediction
        /// of one-hot vectors.
        /// </summary>
        /// <returns>
        /// This instance seen as an <see cref="IMultilabelLikelihoodClassifier{TInput}" />.
        /// </returns>
        new public IMultilabelLikelihoodClassifier<TInput> ToMultilabel()
        {
            return (IMultilabelLikelihoodClassifier<TInput>)this;
        }

    }
}
