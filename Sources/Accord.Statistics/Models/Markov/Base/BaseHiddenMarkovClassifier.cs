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

namespace FileFormat.Accord.Statistics.Models.Markov.Base
{
    using System;
    using System.Collections.Generic;
    using System.Runtime.Serialization;
    using System.Threading.Tasks;
    using global::Accord.Math;
    using Math;
    using Math.Matrix;

    /// <summary>
    ///   Base class for (HMM) Sequence Classifiers. 
    ///   This class cannot be instantiated.
    /// </summary>
    /// 
    [Obsolete("This class will be removed.")]
    [Serializable]
    public abstract class BaseHiddenMarkovClassifier<TModel> : IEnumerable<TModel>
        where TModel : IHiddenMarkovModel
    {

        private TModel[] models;
        private double[] classPriors;

        // Threshold (rejection) model
        private TModel threshold;
        private double weight = 1;


        /// <summary>
        ///   Initializes a new instance of the <see cref="BaseHiddenMarkovClassifier&lt;T&gt;"/> class.
        /// </summary>
        /// <param name="classes">The number of classes in the classification problem.</param>
        /// 
        protected BaseHiddenMarkovClassifier(int classes)
        {
            this.models = new TModel[classes];

            this.classPriors = new double[classes];
            for (int i = 0; i < this.classPriors.Length; i++)
                this.classPriors[i] = 1.0 / this.classPriors.Length;
        }

        /// <summary>
        ///   Initializes a new instance of the <see cref="BaseHiddenMarkovClassifier&lt;T&gt;"/> class.
        /// </summary>
        /// <param name="models">The models specializing in each of the classes of the classification problem.</param>
        /// 
        protected BaseHiddenMarkovClassifier(TModel[] models)
        {
            this.models = models;

            this.classPriors = new double[models.Length];
            for (int i = 0; i < this.classPriors.Length; i++)
                this.classPriors[i] = 1.0 / this.classPriors.Length;
        }

        /// <summary>
        ///   Gets or sets the threshold model.
        /// </summary>
        /// 
        /// <remarks>
        /// <para>
        ///   For gesture spotting, Lee and Kim introduced a threshold model which is
        ///   composed of parts of the models in a hidden Markov sequence classifier.</para>
        /// <para>
        ///   The threshold model acts as a baseline for decision rejection. If none of
        ///   the classifiers is able to produce a higher likelihood than the threshold
        ///   model, the decision is rejected.</para>
        /// <para>
        ///   In the original Lee and Kim publication, the threshold model is constructed
        ///   by creating a fully connected ergodic model by removing all outgoing transitions
        ///   of states in all gesture models and fully connecting those states.</para>
        /// <para>
        ///   References:
        ///   <list type="bullet">
        ///     <item><description>
        ///        H. Lee, J. Kim, An HMM-based threshold model approach for gesture
        ///        recognition, IEEE Trans. Pattern Anal. Mach. Intell. 21 (10) (1999)
        ///        961–973.</description></item>
        ///   </list></para>
        /// </remarks>
        /// 
        public TModel Threshold
        {
            get { return this.threshold; }
            set { this.threshold = value; }
        }

        /// <summary>
        ///   Gets or sets a value governing the rejection given by
        ///   a threshold model (if present). Increasing this value
        ///   will result in higher rejection rates. Default is 1.
        /// </summary>
        /// 
        public double Sensitivity
        {
            get { return this.weight; }
            set { this.weight = value; }
        }


        /// <summary>
        ///   Gets the collection of models specialized in each 
        ///   class of the sequence classification problem.
        /// </summary>
        /// 
        public TModel[] Models
        {
            get { return this.models; }
        }

        /// <summary>
        ///   Gets the <see cref="IHiddenMarkovModel">Hidden Markov
        ///   Model</see> implementation responsible for recognizing
        ///   each of the classes given the desired class label.
        /// </summary>
        /// <param name="label">The class label of the model to get.</param>
        /// 
        public TModel this[int label]
        {
            get { return this.models[label]; }
        }

        /// <summary>
        ///   Gets the number of classes which can be recognized by this classifier.
        /// </summary>
        /// 
        public int Classes
        {
            get { return this.models.Length; }
        }

        /// <summary>
        ///   Gets the prior distribution assumed for the classes.
        /// </summary>
        /// 
        public double[] Priors
        {
            get { return this.classPriors; }
        }

        /// <summary>
        ///   Computes the most likely class for a given sequence.
        /// </summary>
        /// 
        /// <param name="sequence">The sequence of observations.</param>
        /// 
        /// <returns>Return the label of the given sequence, or -1 if it has
        /// been rejected by the <see cref="Threshold">threshold model</see>.</returns>
        /// 
        protected int Compute(Array sequence)
        {
            double[] likelihoods;
            return this.Compute(sequence, out likelihoods);
        }

        /// <summary>
        ///   Computes the most likely class for a given sequence.
        /// </summary>
        /// 
        /// <param name="sequence">The sequence of observations.</param>
        /// <param name="response">The probability of the assigned class.</param>
        /// 
        /// <returns>Return the label of the given sequence, or -1 if it has
        /// been rejected by the <see cref="Threshold">threshold model</see>.</returns>
        /// 
        protected int Compute(Array sequence, out double response)
        {
            double[] likelihoods;
            int label = this.Compute(sequence, out likelihoods);

            if (label == -1)
            {
                // The sequence has been rejected
                response = 1.0 - likelihoods.Sum();
            }
            else
            {
                // Successful recognition
                response = likelihoods[label];
            }

            return label;
        }


        /// <summary>
        ///   Computes the most likely class for a given sequence.
        /// </summary>
        /// 
        /// <param name="sequence">The sequence of observations.</param>
        /// <param name="responsibilities">The probabilities for each class.</param>
        /// 
        /// <returns>Return the label of the given sequence, or -1 if it has
        /// been rejected by the <see cref="Threshold">threshold model</see>.</returns>
        /// 
        protected int Compute(Array sequence, out double[] responsibilities)
        {
            double maxValue, thresholdValue;
            int imax = this.compute(sequence, out responsibilities, out thresholdValue, out maxValue);

            // Convert to probabilities
            for (int i = 0; i < responsibilities.Length; i++)
                responsibilities[i] = Math.Exp(responsibilities[i]);

            // return the index of the most likely model
            // or -1 if the sequence has been rejected.
            //
            return (thresholdValue > maxValue) ? -1 : imax;
        }

        private int compute(Array sequence, out double[] logLikelihoods, 
            out double rejectionValue, out double maxValue)
        {
            logLikelihoods = new double[this.models.Length];

            double[] responses = logLikelihoods;
            double rejection = Double.NegativeInfinity;

            // For every model in the set (including the threshold model)
#if SERIAL
            for (int i = 0; i < models.Length+1; i++)
#else
            Parallel.For(0, this.models.Length + 1, i =>
#endif
            {
                if (i < this.models.Length)
                {
                    // Evaluate the probability of the sequence
                    responses[i] = this.models[i].Evaluate(sequence) + Math.Log(this.classPriors[i]);
                }
                else if (this.threshold != null)
                {
                    // Evaluate the current rejection threshold 
                    rejection = this.threshold.Evaluate(sequence) + Math.Log(this.weight);
                }
            }
#if !SERIAL
);
#endif

            logLikelihoods = responses;
            rejectionValue = rejection;

            // Compute posterior likelihoods
            double lnsum = Double.NegativeInfinity;
            for (int i = 0; i < this.classPriors.Length; i++)
                lnsum = Special.LogSum(lnsum, logLikelihoods[i]);

            // Compute threshold model posterior likelihood
            if (this.threshold != null)
                lnsum = Special.LogSum(lnsum, rejection);

            int imax; 

            // Get the index of the most likely model
            maxValue = logLikelihoods.Max(out imax);

            // Normalize if different from zero
            if (lnsum != Double.NegativeInfinity)
            {
                for (int i = 0; i < logLikelihoods.Length; i++)
                    logLikelihoods[i] -= lnsum;
            }

            return imax;
        }

        /// <summary>
        ///   Computes the log-likelihood that a sequence
        ///   belongs to a given class according to this
        ///   classifier.
        /// </summary>
        /// <param name="sequence">The sequence of observations.</param>
        /// <param name="output">The output class label.</param>
        /// 
        /// <returns>The log-likelihood of the sequence belonging to the given class.</returns>
        /// 
        protected double LogLikelihood(Array sequence, int output)
        {
            double[] logLikelihoods;

            double rejectionValue, maxValue;
            this.compute(sequence, out logLikelihoods, out rejectionValue, out maxValue);

            return logLikelihoods[output];
        }

        /// <summary>
        ///   Computes the log-likelihood that a sequence
        ///   belongs any of the classes in the classifier.
        /// </summary>
        /// <param name="sequence">The sequence of observations.</param>
        /// 
        /// <returns>The log-likelihood of the sequence belonging to the classifier.</returns>
        /// 
        protected double LogLikelihood(Array sequence)
        {
            double sum = Double.NegativeInfinity;

            for (int i = 0; i < this.models.Length; i++)
            {
                double prior = Math.Log(this.classPriors[i]);
                double model = this.models[i].Evaluate(sequence);
                double result = Special.LogSum(prior, model);

                sum = Special.LogSum(sum, result);
            }

            return sum;
        }


        /// <summary>
        ///   Computes the log-likelihood of a set of sequences
        ///   belonging to their given respective classes according
        ///   to this classifier.
        /// </summary>
        /// <param name="sequences">A set of sequences of observations.</param>
        /// <param name="outputs">The output class label for each sequence.</param>
        /// 
        /// <returns>The log-likelihood of the sequences belonging to the given classes.</returns>
        /// 
        protected double LogLikelihood(Array[] sequences, int[] outputs)
        {
            double sum = 0;
            for (int i = 0; i < sequences.Length; i++)
            {
                double[] logLikelihoods;

                double rejectionValue, maxValue;
                this.compute(sequences[i], out logLikelihoods, out rejectionValue, out maxValue);

                sum += logLikelihoods[outputs[i]];
            }

            return sum;
        }


        [OnDeserializing]
        private void setSerializationDefaults(StreamingContext sc)
        {
            // In case this was an older model which didn't include a
            // configurable rejection threshold, initialize it with 1.

            this.weight = 1;
        }

        /// <summary>
        ///   Returns an enumerator that iterates through the models in the classifier.
        /// </summary>
        /// 
        /// <returns>
        ///   A <see cref="T:System.Collections.Generic.IEnumerator`1"/> that 
        ///   can be used to iterate through the collection.
        /// </returns>
        /// 
        public IEnumerator<TModel> GetEnumerator()
        {
            foreach (var model in this.models)
                yield return model;
        }

        /// <summary>
        ///   Returns an enumerator that iterates through the models in the classifier.
        /// </summary>
        /// 
        /// <returns>
        ///   A <see cref="T:System.Collections.Generic.IEnumerator`1"/> that 
        ///   can be used to iterate through the collection.
        /// </returns>
        /// 
        System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator()
        {
            foreach (var model in this.models)
                yield return model;
        }

    }

}
