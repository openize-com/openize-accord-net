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

namespace Openize.Accord.Statistics.Filters
{
    using System;
    using System.Collections.Generic;
    using System.Data;
    using System.Linq.Expressions;
    using Accord.MachineLearning.Learning;
    using Base;
    using Openize.Accord.Core.Collections;
    using Openize.Accord.Core.MachineLearning;

    public partial class Discretization<TInput, TOutput>
    {
        /// <summary>
        ///   Options for the discretization filter.
        /// </summary>
        /// 
        [Serializable]
        public class Options : ColumnOptionsBase<Discretization<TInput, TOutput>>,
            ITransform<TInput, TOutput>,
            IUnsupervisedLearning<Options, TInput, TOutput>
        {
            /// <summary>
            ///   Gets the map between matching rules and the output 
            ///   that should be produced/inserted when they match.
            /// </summary>
            /// 
            public TwoWayDictionary<Expression<Func<TInput, bool>>, Expression<Func<TInput, TOutput>>> Mapping { get; private set; }

            [NonSerialized]
            private Dictionary<Expression<Func<TInput, bool>>, Func<TInput, bool>> keyCache = new Dictionary<Expression<Func<TInput, bool>>, Func<TInput, bool>>();

            [NonSerialized]
            private Dictionary<Expression<Func<TInput, TOutput>>, Func<TInput, TOutput>> valueCache = new Dictionary<Expression<Func<TInput, TOutput>>, Func<TInput, TOutput>>();

            /// <summary>
            ///   Gets the number of symbols used to code this variable.
            /// </summary>
            /// 
            public int NumberOfSymbols
            {
                get { return this.Mapping.Count; }
            }

            /// <summary>
            ///   Gets the number of inputs accepted by the model (value will be 1).
            /// </summary>
            /// 
            /// <value>The number of inputs.</value>
            /// 
            public int NumberOfInputs { get; set; }

            /// <summary>
            /// Gets the number of outputs generated by the model (value will be 1).
            /// </summary>
            /// 
            /// <value>The number of outputs.</value>
            /// 
            public int NumberOfOutputs { get; set; }

            /// <summary>
            /// Applies the transformation to an input, producing an associated output.
            /// </summary>
            /// <param name="input">The input data to which the transformation should be applied.</param>
            /// <returns>The output generated by applying this transformation to the given input.</returns>
            public TOutput Transform(TInput input)
            {
                foreach (var pair in this.Mapping)
                {
                    Func<TInput, bool> func;
                    if (!this.keyCache.TryGetValue(pair.Key, out func))
                        this.keyCache[pair.Key] = func = pair.Key.Compile();

                    if (func(input))
                    {
                        Func<TInput, TOutput> transform;
                        if (!this.valueCache.TryGetValue(pair.Value, out transform))
                            this.valueCache[pair.Value] = transform = pair.Value.Compile();

                        return transform(input);
                    }
                }

                throw new Exception("The value could not be transformed using any of the rules. " +
                    "Please verify if the provided set of rules can cover all possible inputs.");
            }

            /// <summary>
            /// Applies the transformation to a set of input vectors,
            /// producing an associated set of output vectors.
            /// </summary>
            /// <param name="input">The input data to which
            /// the transformation should be applied.</param>
            /// <returns>The output generated by applying this
            /// transformation to the given input.</returns>
            public TOutput[] Transform(TInput[] input)
            {
                return this.Transform(input, new TOutput[input.Length]);
            }


            /// <summary>
            /// Applies the transformation to a set of input vectors,
            /// producing an associated set of output vectors.
            /// </summary>
            /// <param name="input">The input data to which
            /// the transformation should be applied.</param>
            /// <param name="result">The location to where to store the
            /// result of this transformation.</param>
            /// <returns>The output generated by applying this
            /// transformation to the given input.</returns>
            public TOutput[] Transform(TInput[] input, TOutput[] result)
            {
                for (int i = 0; i < input.Length; i++)
                    result[i] = this.Transform(input[i]);
                return result;
            }

#if !NETSTANDARD1_4
            /// <summary>
            /// Applies the transformation to an input, producing an associated output.
            /// </summary>
            /// <param name="input">The input data to which the transformation should be applied.</param>
            /// <returns>The output generated by applying this transformation to the given input.</returns>
            public TOutput Transform(DataRow input)
            {
                return this.Transform((TInput)input[this.ColumnName]);
            }

            /// <summary>
            /// Applies the transformation to an input, producing an associated output.
            /// </summary>
            /// <param name="input">The input data to which the transformation should be applied.</param>
            /// <returns>The output generated by applying this transformation to the given input.</returns>
            public TOutput[] Transform(DataRow[] input)
            {
                return this.Transform(input, new TOutput[input.Length]);
            }

            /// <summary>
            /// Applies the transformation to a set of input vectors,
            /// producing an associated set of output vectors.
            /// </summary>
            /// <param name="input">The input data to which
            /// the transformation should be applied.</param>
            /// <param name="result">The location to where to store the
            /// result of this transformation.</param>
            /// <returns>The output generated by applying this
            /// transformation to the given input.</returns>
            public TOutput[] Transform(DataRow[] input, TOutput[] result)
            {
                for (int i = 0; i < input.Length; i++)
                    result[i] = this.Transform(input[i]);
                return result;
            }
#endif

            /// <summary>
            /// Learns a model that can map the given inputs to the desired outputs.
            /// </summary>
            /// <param name="x">The model inputs.</param>
            /// <param name="weights">The weight of importance for each input sample.</param>
            /// <returns>A model that has learned how to produce suitable outputs
            /// given the input data <paramref name="x" />.</returns>
            /// <exception cref="System.ArgumentException">Weights are not supported and should be null.</exception>
            public Options Learn(TInput[] x, double[] weights = null)
            {
                if (weights != null)
                    throw new ArgumentException(global::Accord.Properties.Resources.NotSupportedWeights, "weights");

                this.NumberOfInputs = 1;
                this.NumberOfOutputs = 1;

                return this;
            }

#if !NETSTANDARD1_4
            /// <summary>
            /// Learns a model that can map the given inputs to the desired outputs.
            /// </summary>
            /// <param name="x">The model inputs.</param>
            /// <param name="weights">The weight of importance for each input sample.</param>
            /// <returns>A model that has learned how to produce suitable outputs
            /// given the input data <paramref name="x" />.</returns>
            /// <exception cref="System.ArgumentException">Weights are not supported and should be null.</exception>
            public Options Learn(DataTable x, double[] weights = null)
            {
                if (weights != null)
                    throw new ArgumentException("Weights are not supported and should be null.");

                this.NumberOfInputs = 1;
                this.NumberOfOutputs = 1;

                return this;
            }

            /// <summary>
            /// Learns a model that can map the given inputs to the desired outputs.
            /// </summary>
            /// <param name="x">The model inputs.</param>
            /// <param name="weights">The weight of importance for each input sample.</param>
            /// <returns>A model that has learned how to produce suitable outputs
            /// given the input data <paramref name="x" />.</returns>
            /// <exception cref="System.ArgumentException">Weights are not supported and should be null.</exception>
            public Options Learn(DataRow[] x, double[] weights = null)
            {
                if (weights != null)
                    throw new ArgumentException(global::Accord.Properties.Resources.NotSupportedWeights, "weights");

                this.NumberOfInputs = 1;
                this.NumberOfOutputs = 1;

                return this;
            }
#endif

            /// <summary>
            /// Learns a model that can map the given inputs to the desired outputs.
            /// </summary>
            /// <param name="x">The model inputs.</param>
            /// <param name="weights">The weight of importance for each input sample.</param>
            /// <returns>A model that has learned how to produce suitable outputs
            /// given the input data <paramref name="x" />.</returns>
            /// <exception cref="System.ArgumentException">Weights are not supported and should be null.</exception>
            public Options Learn(object[] x, double[] weights = null)
            {
                if (weights != null)
                    throw new ArgumentException("Weights are not supported and should be null.");

                this.NumberOfInputs = 1;
                this.NumberOfOutputs = 1;

                return this;
            }

            /// <summary>
            /// Computes a class-label decision for a given <paramref name="input" />.
            /// </summary>
            /// <param name="input">The input vector that should be classified into
            /// one of the <see cref="P:Accord.MachineLearning.ITransform.NumberOfOutputs" /> possible classes.</param>
            /// <returns>A class-label that best described <paramref name="input" /> according
            /// to this classifier.</returns>
            public TOutput Decide(TInput input)
            {
                return this.Transform(input);
            }

            /// <summary>
            /// Computes class-label decisions for each vector in the given <paramref name="input" />.
            /// </summary>
            /// <param name="input">The input vectors that should be classified into
            /// one of the <see cref="P:Accord.MachineLearning.ITransform.NumberOfOutputs" /> possible classes.</param>
            /// <returns>The class-labels that best describe each <paramref name="input" />
            /// vectors according to this classifier.</returns>
            public TOutput[] Decide(TInput[] input)
            {
                return this.Transform(input);
            }

            /// <summary>
            /// Computes class-label decisions for each vector in the given <paramref name="input" />.
            /// </summary>
            /// <param name="input">The input vectors that should be classified into
            /// one of the <see cref="P:Accord.MachineLearning.ITransform.NumberOfOutputs" /> possible classes.</param>
            /// <param name="result">The location where to store the class-labels.</param>
            /// <returns>The class-labels that best describe each <paramref name="input" />
            /// vectors according to this classifier.</returns>
            public TOutput[] Decide(TInput[] input, TOutput[] result)
            {
                return this.Transform(input, result);
            }

            /// <summary>
            ///   Constructs a new Options object.
            /// </summary>
            /// 
            public Options()
                : this("New column")
            {
            }

            /// <summary>
            ///   Constructs a new Options object for the given column.
            /// </summary>
            /// 
            /// <param name="name">
            ///   The name of the column to create this options for.
            /// </param>
            /// 
            public Options(String name)
                : this(name, new Dictionary<Expression<Func<TInput, bool>>, Expression<Func<TInput, TOutput>>>())
            {
            }

            /// <summary>
            ///   Constructs a new Options object for the given column.
            /// </summary>
            /// 
            /// <param name="name">
            ///   The name of the column to create this options for.
            /// </param>
            /// 
            /// <param name="map">The initial mapping for this column.</param>
            /// 
            public Options(String name, Dictionary<Expression<Func<TInput, bool>>, Expression<Func<TInput, TOutput>>> map)
                : base(name)
            {
                this.Mapping = new TwoWayDictionary<Expression<Func<TInput, bool>>, Expression<Func<TInput, TOutput>>>(map);
                this.NumberOfInputs = 1;
                this.NumberOfOutputs = 1;
            }
        }
    }
}
