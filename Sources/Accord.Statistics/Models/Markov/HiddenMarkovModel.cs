// Accord Statistics Library
// The Accord.NET Framework
// http://accord-framework.net
//
// Copyright © César Souza, 2009-2017
// Copyright © Guilherme Pedroso, 2009
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

namespace Openize.Accord.Statistics.Models.Markov
{
#pragma warning disable 612, 618
    using System;
    using System.IO;
    using System.Runtime.Serialization.Formatters.Binary;
    using Base;
    using Distributions;
    using Openize.Accord.Math;
    using Openize.Accord.Math.Core;
    using Openize.Accord.Math.Matrix;
    using Openize.Accord.Math.Vector;
    using Math;
    using Openize.Accord.Statistics.Distributions.Univariate.Discrete;
    using Topology;

    /// <summary>
    ///   Discrete-density Hidden Markov Model.
    /// </summary>
    /// 
    /// <remarks>
    /// <para>
    ///   Hidden Markov Models (HMM) are stochastic methods to model temporal and sequence
    ///   data. They are especially known for their application in temporal pattern recognition
    ///   such as speech, handwriting, gesture recognition, part-of-speech tagging, musical
    ///   score following, partial discharges and bioinformatics.</para>
    ///   
    /// <para>
    ///   This page refers to the discrete-density version of the model. For arbitrary
    ///   density (probability distribution) definitions, please see 
    ///   <see cref="HiddenMarkovModel{TDistribution}"/>.
    /// </para>
    ///   
    /// <para>
    ///   Dynamical systems of discrete nature assumed to be governed by a Markov chain emits
    ///   a sequence of observable outputs. Under the Markov assumption, it is also assumed that
    ///   the latest output depends only on the current state of the system. Such states are often
    ///   not known from the observer when only the output values are observable.</para>
    ///   
    ///  <para>
    ///   Assuming the Markov probability, the probability of any sequence of observations 
    ///   occurring when following a given sequence of states can be stated as</para>
    ///      
    ///  <p align="center">
    ///      <img src="..\images\hmm\hmm-joint-probability.png" width="383" height="133" /></p>
    ///      
    ///  <para>
    ///  in which the probabilities <c>p(y<sub>t</sub>|y<sub>t-1</sub>)</c> can be read as the 
    ///  probability of being currently in state <c>y<sub>t</sub></c> given we just were in the
    ///  state<c> y<sub>t-1</sub></c> at the previous instant <c>t-1</c>, and the probability
    ///  <c> p(x<sub>t</sub>|y<sub>t</sub>)</c> can be understood as the probability of observing 
    ///  <c><strong>x<sub>t</sub></strong></c> at instant t given we are currently in the state 
    ///  <c>y<sub>t</sub></c>. To compute those probabilities, we simple use two matrices <strong>
    ///  <c><strong>A</strong></c></strong> and <strong><c><strong>B</strong></c></strong>. 
    ///  The matrix <strong><c><strong>A</strong></c></strong> is the matrix of state probabilities:
    ///  it gives the probabilities <c>p(y<sub>t</sub>|y<sub>t-1</sub>)</c> of jumping from one state
    ///  to the other, and the matrix B is the matrix of observation probabilities, which gives the
    ///  distribution density <c>p(<strong>x<sub>t</sub></strong>|y<sub>t</sub>)</c> associated 
    ///  a given state <c>y<sub>t</sub></c>. In the discrete case, <c><strong><c><strong>
    ///  B</strong></c></strong></c> is really a matrix. In the continuous case, <c><strong>
    ///  B</strong></c> is a vector of probability distributions. The overall model definition
    ///  can then be stated by the tuple</para>
    ///  
    ///  <p align="center">
    ///      <img src="..\images\hmm\hmm-tuple.png" width="159" height="42" /></p>
    ///      
    /// <para>
    ///  in which <em><c><em>n</em></c></em> is an integer representing the total number 
    ///  of states in the system, <strong><c><strong>A</strong></c></strong> is a matrix 
    ///  of transition probabilities, <strong><c><strong>B</strong></c></strong> is either
    ///  a matrix of observation probabilities (in the discrete case) or a vector of probability
    ///  distributions (in the general case) and <c><strong>p</strong></c> is a vector of 
    ///  initial state probabilities determining the probability of starting in each of the 
    ///  possible states in the model.</para>
    ///   
    /// <para>
    ///   Hidden Markov Models attempt to model such systems and allow, among other things,
    ///   <list type="number">
    ///     <item><description>
    ///       To infer the most likely sequence of states that produced a given output sequence,</description></item>
    ///     <item><description>
    ///       Infer which will be the most likely next state (and thus predicting the next output),</description></item>
    ///     <item><description>
    ///       Calculate the probability that a given sequence of outputs originated from the system
    ///       (allowing the use of hidden Markov models for sequence classification).</description></item>
    ///     </list></para>
    ///     
    /// <para>     
    ///   The “hidden” in Hidden Markov Models comes from the fact that the observer does not
    ///   know in which state the system may be in, but has only a probabilistic insight on where
    ///   it should be.</para>
    ///   
    /// <para>
    ///   To learn a Markov model, you can find a list of both <see cref="ISupervisedLearning">
    ///   supervised</see> and <see cref="IUnsupervisedLearning">unsupervised</see> learning 
    ///   algorithms in the <see cref="Accord.Statistics.Models.Markov.Learning"/> namespace.</para>
    ///   
    /// <para>
    ///   References:
    ///   <list type="bullet">
    ///     <item><description>
    ///       Wikipedia contributors. "Linear regression." Wikipedia, the Free Encyclopedia.
    ///       Available at: http://en.wikipedia.org/wiki/Hidden_Markov_model </description></item>
    ///     <item><description>
    ///       Nikolai Shokhirev, Hidden Markov Models. Personal website. Available at:
    ///       http://www.shokhirev.com/nikolai/abc/alg/hmm/hmm.html </description></item>
    ///     <item><description>
    ///       X. Huang, A. Acero, H. Hon. "Spoken Language Processing." pp 396-397. 
    ///       Prentice Hall, 2001.</description></item>
    ///     <item><description>
    ///       Dawei Shen. Some mathematics for HMMs, 2008. Available at:
    ///       http://courses.media.mit.edu/2010fall/mas622j/ProblemSets/ps4/tutorial.pdf </description></item>
    ///   </list></para>
    /// </remarks>
    /// 
    /// <example>
    ///   <para>The example below reproduces the same example given in the Wikipedia
    ///   entry for the Viterbi algorithm (http://en.wikipedia.org/wiki/Viterbi_algorithm).
    ///   In this example, the model's parameters are initialized manually. However, it is
    ///   possible to learn those automatically using <see cref="BaumWelchLearning"/>.</para>
    ///   
    /// <code source="Unit Tests\Accord.Tests.Statistics\Models\Markov\HiddenMarkovModelTest.cs" region="doc_decode" />
    /// 
    /// <para>If you would like to learn the a hidden Markov model straight
    ///   from a dataset, you can use:</para>
    ///   
    /// <code source="Unit Tests\Accord.Tests.Statistics\Models\Markov\HiddenMarkovModelTest.cs" region="doc_learn"/>
    /// 
    /// <para>Hidden Markov Models are generative models, and as such, can be used
    ///   to generate new samples following the structure that they have learned from
    ///   the data</para>
    /// 
    /// <code source="Unit Tests\Accord.Tests.Statistics\Models\Markov\HiddenMarkovModelTest.cs" region="doc_generate"/>
    ///
    /// <para>Hidden Markov Models can also be used to predict the next observation in a sequence. This can be done by
    ///   inspecting the forward matrix of probabilities for the sequence and inspecting all states and possible symbols
    ///   to find which state-observation combination would be the most likely after the current ones. This limits the 
    ///   applicability of this model to only very short-term predictions (i.e. most likely, only the most immediate next 
    ///   observation). </para>
    /// 
    /// <code source="Unit Tests\Accord.Tests.Statistics\Models\Markov\HiddenMarkovModelTest.cs" region="doc_predict"/>
    /// 
    /// <para>For more examples on how to learn discrete models, please see the <see cref="BaumWelchLearning"/> documentation
    ///   page. For continuous models (models that can model more than just integer labels),
    ///   please see <see cref="BaumWelchLearning{TDistribution, TObservation}"/>.</para>
    /// </example>
    /// 
    /// <seealso cref="BaumWelchLearning">Baum-Welch, one of the most famous 
    ///   learning algorithms for Hidden Markov Models.</seealso>
    /// <seealso cref="HiddenMarkovModel{TDistribution}">Arbitrary-density 
    ///   Hidden Markov Model.</seealso>
    /// <seealso cref="Accord.Statistics.Models.Markov.Learning"/>
    /// 
    [Serializable]
    public class HiddenMarkovModel : HiddenMarkovModel<GeneralDiscreteDistribution, int>, IHiddenMarkovModel, ICloneable
    {

        private double[][] logB; // emission probabilities

        // The other parameters are defined in HiddenMarkovModelBase
        // private double[,] A; // Transition probabilities
        // private double[] pi; // Initial state probabilities


        // Size of vocabulary
        private int symbols;

        /// <summary>
        ///   Please use <see cref="HiddenMarkovModel{T, U}.LogTransitions"/> instead.
        /// </summary>
        /// 
        [Obsolete("Please use LogTransitions instead.")]
        public double[,] Transitions
        {
            get
            {
                return base.LogTransitions.ToMatrix();
            }
        }

        /// <summary>
        ///   Please use <see cref="HiddenMarkovModel{T, U}.LogInitial"/> instead.
        /// </summary>
        /// 
        [Obsolete("Please use LogInitial instead.")]
        public new double[] Probabilities { get { return base.LogInitial; } }


        /// <summary>
        ///   Please use <see cref="NumberOfSymbols"/> instead.
        /// </summary>
        /// 
        [Obsolete("Please use NumberOfSymbols instead.")]
        public int Symbols
        {
            get { return this.symbols; }
        }

        /// <summary>
        ///   Gets the number of symbols in this model's alphabet.
        /// </summary>
        /// 
        public int NumberOfSymbols
        {
            get { return this.symbols; }
        }

        /// <summary>
        ///   Please use <see cref="LogEmissions"/> instead.
        /// </summary>
        /// 
        [Obsolete("Please use LogEmissions instead.")]
        public new double[,] Emissions
        {
            get { return this.LogEmissions.ToMatrix(); }
        }

        /// <summary>
        ///   Gets the log-emission matrix <c>log(B)</c> for this model.
        /// </summary>
        /// 
        public double[][] LogEmissions
        {
            get
            {
                if (this.logB == null)
                {
                    this.logB = new double[base.Emissions.Length][];
                    for (int i = 0; i < base.Emissions.Length; i++)
                        this.logB[i] = base.Emissions[i].Frequencies;
                }

                return this.logB;
            }
        }


        /// <summary>
        ///   Constructs a new Hidden Markov Model.
        /// </summary>
        /// 
        /// <param name="topology">
        ///   A <see cref="Topology"/> object specifying the initial values of the matrix of transition 
        ///   probabilities <c>A</c> and initial state probabilities <c>pi</c> to be used by this model.
        /// </param>
        /// <param name="emissions">The emissions matrix B for this model.</param>
        /// <param name="logarithm">Set to true if the matrices are given with logarithms of the
        /// intended probabilities; set to false otherwise. Default is false.</param>
        /// 
        public HiddenMarkovModel(ITopology topology, double[,] emissions, bool logarithm = false)
            : base(topology, GeneralDiscreteDistribution.FromMatrix((logarithm) ? emissions : emissions.Log(), true))
        {
            this.symbols = base.Emissions[0].Length;
        }

        /// <summary>
        ///   Constructs a new Hidden Markov Model.
        /// </summary>
        /// 
        /// <param name="topology">
        ///   A <see cref="Topology"/> object specifying the initial values of the matrix of transition 
        ///   probabilities <c>A</c> and initial state probabilities <c>pi</c> to be used by this model.
        /// </param>
        /// <param name="emissions">The emissions matrix B for this model.</param>
        /// <param name="logarithm">Set to true if the matrices are given with logarithms of the
        /// intended probabilities; set to false otherwise. Default is false.</param>
        /// 
        public HiddenMarkovModel(ITopology topology, double[][] emissions, bool logarithm = false)
            : base(topology, GeneralDiscreteDistribution.FromMatrix((logarithm) ? emissions : emissions.Log(), true))
        {
            this.symbols = base.Emissions[0].Length;
        }

        /// <summary>
        ///   Constructs a new Hidden Markov Model.
        /// </summary>
        /// 
        /// <param name="topology">
        ///   A <see cref="Topology"/> object specifying the initial values of the matrix of transition 
        ///   probabilities <c>A</c> and initial state probabilities <c>pi</c> to be used by this model.
        /// </param>
        /// <param name="symbols">The number of output symbols used for this model.</param>
        /// 
        public HiddenMarkovModel(ITopology topology, int symbols)
            : this(topology, symbols, false)
        {
        }

        /// <summary>
        ///   Constructs a new Hidden Markov Model.
        /// </summary>
        /// 
        /// <param name="topology">
        ///   A <see cref="Topology"/> object specifying the initial values of the matrix of transition 
        ///   probabilities <c>A</c> and initial state probabilities <c>pi</c> to be used by this model.
        /// </param>
        /// <param name="symbols">The number of output symbols used for this model.</param>
        /// <param name="random">Whether to initialize emissions with random probabilities
        ///   or uniformly with <c>1 / number of symbols</c>. Default is false (default is
        ///   to use <c>1/symbols</c>).</param>
        /// 
        public HiddenMarkovModel(ITopology topology, int symbols, bool random)
            : base(topology)
        {
            if (symbols <= 0)
                throw new ArgumentOutOfRangeException("symbols");

            this.symbols = symbols;

            double[][] b;
            if (random)
                b = Jagged.Random(this.States, symbols);
            else
                b = Jagged.Ones(this.States, symbols);

            base.Emissions = GeneralDiscreteDistribution.FromMatrix(b.Log(), true);
        }

        /// <summary>
        ///   Constructs a new Hidden Markov Model.
        /// </summary>
        /// 
        /// <param name="transitions">The transitions matrix A for this model.</param>
        /// <param name="emissions">The emissions matrix B for this model.</param>
        /// <param name="initial">The initial state probabilities for this model.</param>
        /// <param name="logarithm">Set to true if the matrices are given with logarithms of the
        /// intended probabilities; set to false otherwise. Default is false.</param>
        /// 
        public HiddenMarkovModel(double[,] transitions, double[,] emissions, double[] initial, bool logarithm = false)
            : this(new Custom(transitions, initial, logarithm), emissions, logarithm)
        {
        }

        /// <summary>
        ///   Constructs a new Hidden Markov Model.
        /// </summary>
        /// 
        /// <param name="transitions">The transitions matrix A for this model.</param>
        /// <param name="emissions">The emissions matrix B for this model.</param>
        /// <param name="initial">The initial state probabilities for this model.</param>
        /// <param name="logarithm">Set to true if the matrices are given with logarithms of the
        /// intended probabilities; set to false otherwise. Default is false.</param>
        /// 
        public HiddenMarkovModel(double[][] transitions, double[][] emissions, double[] initial, bool logarithm = false)
            : this(new Custom(transitions, initial, logarithm), emissions, logarithm)
        {
        }

        /// <summary>
        ///   Constructs a new Hidden Markov Model.
        /// </summary>
        /// 
        /// <param name="states">The number of states for this model.</param>
        /// <param name="symbols">The number of output symbols used for this model.</param>
        /// 
        public HiddenMarkovModel(int states, int symbols)
            : this(new Ergodic(states), symbols)
        {
        }

        /// <summary>
        ///   Constructs a new Hidden Markov Model.
        /// </summary>
        /// 
        /// <param name="states">The number of states for this model.</param>
        /// <param name="symbols">The number of output symbols used for this model.</param>
        /// <param name="random">Whether to initialize the model transitions and emissions 
        ///   with random probabilities or uniformly with <c>1 / number of states</c> (for
        ///   transitions) and <c>1 / number of symbols</c> (for emissions). Default is false.</param>
        /// 
        public HiddenMarkovModel(int states, int symbols, bool random)
            : this(new Ergodic(states, random), symbols, random)
        {
        }




        /// <summary>
        ///   Predicts next observations occurring after a given observation sequence.
        /// </summary>
        /// 
        /// <param name="observations">A sequence of observations. Predictions will be made regarding 
        ///   the next observations that should be coming after the last observation in this sequence.</param>
        /// <param name="next">The number of observations to be predicted. Default is 1.</param>
        /// 
        /// <example>
        /// <code source="Unit Tests\Accord.Tests.Statistics\Models\Markov\HiddenMarkovModelTest.cs" region="doc_predict"/>
        /// </example>
        /// 
        public override int[] Predict(int[] observations, int next)
        {
            double logLikelihood;
            double[][] logLikelihoods;
            return this.Predict(observations, next, out logLikelihood, out logLikelihoods);
        }

        /// <summary>
        ///   Predicts next observations occurring after a given observation sequence.
        /// </summary>
        /// 
        /// <param name="observations">A sequence of observations. Predictions will be made regarding 
        ///   the next observations that should be coming after the last observation in this sequence.</param>
        /// <param name="next">The number of observations to be predicted. Default is 1.</param>
        /// <param name="logLikelihoods">The log-likelihood of the different symbols for each predicted
        ///   next observations. In order to convert those values to probabilities, exponentiate the
        ///   values in the vectors (using the Exp function) and divide each value by their vector's sum.</param>
        /// 
        /// <example>
        /// <code source="Unit Tests\Accord.Tests.Statistics\Models\Markov\HiddenMarkovModelTest.cs" region="doc_predict"/>
        /// </example>
        /// 
        public int[] Predict(int[] observations, int next, out double[][] logLikelihoods)
        {
            double logLikelihood;
            return this.Predict(observations, next, out logLikelihood, out logLikelihoods);
        }

        /// <summary>
        ///   Predicts the next observation occurring after a given observation sequence.
        /// </summary>
        /// 
        /// <param name="observations">A sequence of observations. Predictions will be made regarding 
        ///   the next observation that should be coming after the last observation in this sequence.</param>
        /// <param name="logLikelihoods">The log-likelihood of the different symbols for the next observation.
        ///   In order to convert those values to probabilities, exponentiate the values in the vector (using
        ///   the Exp function) and divide each value by the vector sum. This will give the probability of each
        ///   next possible symbol to be the next observation in the sequence.</param>
        /// 
        /// <example>
        /// <code source="Unit Tests\Accord.Tests.Statistics\Models\Markov\HiddenMarkovModelTest.cs" region="doc_predict"/>
        /// </example>
        /// 
        public int Predict(int[] observations, out double[] logLikelihoods)
        {
            double[][] ll;
            double logLikelihood;
            int prediction = this.Predict(observations, 1, out logLikelihood, out ll)[0];
            logLikelihoods = ll[0];
            return prediction;
        }

        /// <summary>
        ///   Predicts the next observations occurring after a given observation sequence.
        /// </summary>
        /// 
        /// <param name="observations">A sequence of observations. Predictions will be made regarding 
        ///   the next observations that should be coming after the last observation in this sequence.</param>
        /// <param name="next">The number of observations to be predicted. Default is 1.</param>
        /// <param name="logLikelihoods">The log-likelihood of the different symbols for each predicted
        ///   next observations. In order to convert those values to probabilities, exponentiate the
        ///   values in the vectors (using the Exp function) and divide each value by their vector's sum.</param>
        /// <param name="logLikelihood">The log-likelihood of the given sequence, plus the predicted
        ///   next observations. Exponentiate this value (use the System.Math.Exp function) to obtain
        ///   a <c>likelihood</c> value.</param>
        ///   
        /// <example>
        /// <code source="Unit Tests\Accord.Tests.Statistics\Models\Markov\HiddenMarkovModelTest.cs" region="doc_predict"/>
        /// </example>
        /// 
        public int[] Predict(int[] observations, int next, out double logLikelihood, out double[][] logLikelihoods)
        {
            int states = this.States;
            int T = next;
            double[][] lnA = this.LogTransitions;

            int[] prediction = new int[next];
            logLikelihoods = new double[next][];

            try
            {
                // Compute forward probabilities for the given observation sequence.
                double[,] lnFw0 = ForwardBackwardAlgorithm.LogForward(this, observations, out logLikelihood);

                // Create a matrix to store the future probabilities for the prediction
                // sequence and copy the latest forward probabilities on its first row.
                double[,] lnFwd = new double[T + 1, states];


                // 1. Initialization
                for (int i = 0; i < states; i++)
                    lnFwd[0, i] = lnFw0[observations.Length - 1, i];

                // 2. Induction
                for (int t = 0; t < T; t++)
                {
                    double[] weights = new double[this.symbols];
                    for (int s = 0; s < this.symbols; s++)
                    {
                        weights[s] = Double.NegativeInfinity;

                        for (int i = 0; i < states; i++)
                        {
                            double sum = Double.NegativeInfinity;
                            for (int j = 0; j < states; j++)
                                sum = Special.LogSum(sum, lnFwd[t, j] + lnA[j][i]);
                            lnFwd[t + 1, i] = sum + this.logB[i][s];

                            weights[s] = Special.LogSum(weights[s], lnFwd[t + 1, i]);
                        }
                    }

                    double sumWeight = Double.NegativeInfinity;
                    for (int i = 0; i < weights.Length; i++)
                        sumWeight = Special.LogSum(sumWeight, weights[i]);
                    for (int i = 0; i < weights.Length; i++)
                        weights[i] -= sumWeight;


                    // Select most probable symbol
                    double maxWeight = weights[0];
                    prediction[t] = 0;
                    for (int i = 1; i < weights.Length; i++)
                    {
                        if (weights[i] > maxWeight)
                        {
                            maxWeight = weights[i];
                            prediction[t] = i;
                        }
                    }

                    // Recompute log-likelihood
                    logLikelihoods[t] = weights;
                    logLikelihood = maxWeight;
                }


                return prediction;
            }
            catch (IndexOutOfRangeException ex)
            {
                this.checkObservations(ex, observations);
                throw;
            }
        }


        /// <summary>
        ///   Converts this <see cref="HiddenMarkovModel">Discrete density Hidden Markov Model</see>
        ///   into a <see cref="HiddenMarkovModel{TDistribution}">arbitrary density model</see>.
        /// </summary>
        /// 
        [Obsolete("Please use ToGenericModel instead.")]
        public HiddenMarkovModel<GeneralDiscreteDistribution> ToContinuousModel()
        {
            var transitions = (double[,])this.LogTransitions.Clone();
            var probabilities = (double[])this.LogInitial.Clone();

            var emissions = new GeneralDiscreteDistribution[this.States];
            for (int i = 0; i < emissions.Length; i++)
                emissions[i] = new GeneralDiscreteDistribution(Matrix.GetRow(this.Emissions, i));

            return new HiddenMarkovModel<GeneralDiscreteDistribution>(transitions, emissions, probabilities);
        }

        /// <summary>
        ///   Converts this <see cref="HiddenMarkovModel">Discrete density Hidden Markov Model</see>
        ///   into a <see cref="HiddenMarkovModel{TDistribution, TObservation}">arbitrary density model</see>.
        /// </summary>
        /// 
        public HiddenMarkovModel<GeneralDiscreteDistribution, int> ToGenericModel()
        {
            return (HiddenMarkovModel<GeneralDiscreteDistribution, int>)this;
        }

        /// <summary>
        ///   Converts this <see cref="HiddenMarkovModel">Discrete density Hidden Markov Model</see>
        ///   to a <see cref="HiddenMarkovModel{TDistribution}">Continuous density model</see>.
        /// </summary>
        public static explicit operator HiddenMarkovModel<GeneralDiscreteDistribution>(HiddenMarkovModel model)
        {
            return model.ToContinuousModel();
        }




        /// <summary>
        ///   Constructs a new discrete-density Hidden Markov Model.
        /// </summary>
        /// 
        /// <param name="transitions">The transitions matrix A for this model.</param>
        /// <param name="emissions">The emissions matrix B for this model.</param>
        /// <param name="probabilities">The initial state probabilities for this model.</param>
        /// <param name="logarithm">Set to true if the matrices are given with logarithms of the
        /// intended probabilities; set to false otherwise. Default is false.</param>
        /// 
        [Obsolete("Please use CreateDiscrete() instead.")]
        public static HiddenMarkovModel<GeneralDiscreteDistribution> CreateGeneric(double[,] transitions,
            double[,] emissions, double[] probabilities, bool logarithm = false)
        {
            ITopology topology = new Custom(transitions, probabilities, logarithm);

            if (emissions == null)
            {
                throw new ArgumentNullException("emissions");
            }

            if (emissions.GetLength(0) != topology.States)
            {
                throw new ArgumentException(
                    "The emission matrix should have the same number of rows as the number of states in the model.",
                    "emissions");
            }


            // Initialize B using a discrete distribution
            var B = new GeneralDiscreteDistribution[topology.States];
            for (int i = 0; i < B.Length; i++)
                B[i] = new GeneralDiscreteDistribution(Matrix.GetRow(emissions, i));

            return new HiddenMarkovModel<GeneralDiscreteDistribution>(topology, B);
        }

        /// <summary>
        ///   Constructs a new Hidden Markov Model with discrete state probabilities.
        /// </summary>
        /// 
        /// <param name="topology">
        ///   A <see cref="Topology"/> object specifying the initial values of the matrix of transition 
        ///   probabilities <c>A</c> and initial state probabilities <c>pi</c> to be used by this model.
        /// </param>
        /// <param name="symbols">The number of output symbols used for this model.</param>
        /// 
        [Obsolete("Please use CreateDiscrete() instead.")]
        public static HiddenMarkovModel<GeneralDiscreteDistribution> CreateGeneric(ITopology topology, int symbols)
        {
            return CreateGeneric(topology, symbols, false);
        }

        /// <summary>
        ///   Constructs a new Hidden Markov Model with discrete state probabilities.
        /// </summary>
        /// 
        /// <param name="topology">
        ///   A <see cref="Topology"/> object specifying the initial values of the matrix of transition 
        ///   probabilities <c>A</c> and initial state probabilities <c>pi</c> to be used by this model.
        /// </param>
        /// <param name="symbols">The number of output symbols used for this model.</param>
        /// <param name="random">Whether to initialize emissions with random probabilities
        ///   or uniformly with <c>1 / number of symbols</c>. Default is false (default is
        ///   to use <c>1/symbols</c>).</param>
        /// 
        [Obsolete("Please use CreateDiscrete() instead.")]
        public static HiddenMarkovModel<GeneralDiscreteDistribution> CreateGeneric(ITopology topology, int symbols, bool random)
        {
            if (symbols <= 0)
            {
                throw new ArgumentOutOfRangeException("symbols",
                    "Number of symbols should be higher than zero.");
            }

            double[,] A;
            double[] pi;
            topology.Create(true, out A, out pi);

            // Initialize B with a uniform discrete distribution
            var B = new GeneralDiscreteDistribution[topology.States];

            if (random)
            {
                for (int i = 0; i < B.Length; i++)
                {
                    double[] probabilities = new double[symbols];

                    double sum = 0;
                    for (int j = 0; j < probabilities.Length; j++)
                        sum += probabilities[j] = global::Openize.Accord.Math.Random.Generator.Random.NextDouble();

                    for (int j = 0; j < probabilities.Length; j++)
                        probabilities[j] /= sum;

                    B[i] = new GeneralDiscreteDistribution(true, probabilities.Log());
                }
            }
            else
            {
                for (int i = 0; i < B.Length; i++)
                    B[i] = new GeneralDiscreteDistribution(logarithm: true, symbols: symbols);
            }


            return new HiddenMarkovModel<GeneralDiscreteDistribution>(A, B, pi, logarithm: true);
        }

        /// <summary>
        ///   Constructs a new Hidden Markov Model with discrete state probabilities.
        /// </summary>
        /// 
        /// <param name="states">The number of states for this model.</param>
        /// <param name="symbols">The number of output symbols used for this model.</param>
        /// 
        [Obsolete("Please use CreateDiscrete() instead.")]
        public static HiddenMarkovModel<GeneralDiscreteDistribution> CreateGeneric(int states, int symbols)
        {
            return CreateGeneric(new Ergodic(states), symbols);
        }

        /// <summary>
        ///   Constructs a new Hidden Markov Model with discrete state probabilities.
        /// </summary>
        /// 
        /// <param name="states">The number of states for this model.</param>
        /// <param name="symbols">The number of output symbols used for this model.</param>
        /// <param name="random">Whether to initialize emissions with random probabilities
        ///   or uniformly with <c>1 / number of symbols</c>. Default is false (default is
        ///   to use <c>1/symbols</c>).</param>
        /// 
        [Obsolete("Please use CreateDiscrete() instead.")]
        public static HiddenMarkovModel<GeneralDiscreteDistribution> CreateGeneric(int states, int symbols, bool random)
        {
            return CreateGeneric(new Ergodic(states, random), symbols, random);
        }


        /// <summary>
        ///   Creates a discrete hidden Markov model using the generic interface.
        /// </summary>
        /// 
        /// <param name="transitions">The transitions matrix A for this model.</param>
        /// <param name="emissions">The emissions matrix B for this model.</param>
        /// <param name="probabilities">The initial state probabilities for this model.</param>
        /// <param name="logarithm">Set to true if the matrices are given with logarithms of the
        /// intended probabilities; set to false otherwise. Default is false.</param>
        /// 
        public static HiddenMarkovModel<GeneralDiscreteDistribution, int> CreateDiscrete(double[,] transitions,
            double[,] emissions, double[] probabilities, bool logarithm = false)
        {
            ITopology topology = new Custom(transitions, probabilities, logarithm);

            if (emissions == null)
            {
                throw new ArgumentNullException("emissions");
            }

            if (emissions.GetLength(0) != topology.States)
            {
                throw new ArgumentException(
                    "The emission matrix should have the same number of rows as the number of states in the model.",
                    "emissions");
            }


            // Initialize B using a discrete distribution
            var B = new GeneralDiscreteDistribution[topology.States];
            for (int i = 0; i < B.Length; i++)
                B[i] = new GeneralDiscreteDistribution(Matrix.GetRow(emissions, i));

            return new HiddenMarkovModel<GeneralDiscreteDistribution, int>(topology, B);
        }

        /// <summary>
        ///   Creates a discrete hidden Markov model using the generic interface.
        /// </summary>
        /// 
        /// <param name="topology">
        ///   A <see cref="Topology"/> object specifying the initial values of the matrix of transition 
        ///   probabilities <c>A</c> and initial state probabilities <c>pi</c> to be used by this model.
        /// </param>
        /// <param name="symbols">The number of output symbols used for this model.</param>
        /// 
        public static HiddenMarkovModel<GeneralDiscreteDistribution, int> CreateDiscrete(ITopology topology, int symbols)
        {
            return CreateDiscrete(topology, symbols, false);
        }

        /// <summary>
        ///   Creates a discrete hidden Markov model using the generic interface.
        /// </summary>
        /// 
        /// <param name="topology">
        ///   A <see cref="Topology"/> object specifying the initial values of the matrix of transition 
        ///   probabilities <c>A</c> and initial state probabilities <c>pi</c> to be used by this model.
        /// </param>
        /// <param name="symbols">The number of output symbols used for this model.</param>
        /// <param name="random">Whether to initialize emissions with random probabilities
        ///   or uniformly with <c>1 / number of symbols</c>. Default is false (default is
        ///   to use <c>1/symbols</c>).</param>
        /// 
        public static HiddenMarkovModel<GeneralDiscreteDistribution, int> CreateDiscrete(ITopology topology, int symbols, bool random)
        {
            if (symbols <= 0)
            {
                throw new ArgumentOutOfRangeException("symbols",
                    "Number of symbols should be higher than zero.");
            }

            double[,] A;
            double[] pi;
            topology.Create(true, out A, out pi);

            // Initialize B with a uniform discrete distribution
            var B = new GeneralDiscreteDistribution[topology.States];

            if (random)
            {
                for (int i = 0; i < B.Length; i++)
                {
                    double[] probabilities = new double[symbols];

                    double sum = 0;
                    for (int j = 0; j < probabilities.Length; j++)
                        sum += probabilities[j] = global::Openize.Accord.Math.Random.Generator.Random.NextDouble();

                    for (int j = 0; j < probabilities.Length; j++)
                        probabilities[j] /= sum;

                    B[i] = new GeneralDiscreteDistribution(true, probabilities.Log());
                }
            }
            else
            {
                for (int i = 0; i < B.Length; i++)
                    B[i] = new GeneralDiscreteDistribution(logarithm: true, symbols: symbols);
            }


            return new HiddenMarkovModel<GeneralDiscreteDistribution, int>(A, B, pi, logarithm: true);
        }

        /// <summary>
        ///   Creates a discrete hidden Markov model using the generic interface.
        /// </summary>
        /// 
        /// <param name="states">The number of states for this model.</param>
        /// <param name="symbols">The number of output symbols used for this model.</param>
        /// 
        public static HiddenMarkovModel<GeneralDiscreteDistribution, int> CreateDiscrete(int states, int symbols)
        {
            return CreateDiscrete(new Ergodic(states), symbols);
        }

        /// <summary>
        ///   Creates a discrete hidden Markov model using the generic interface.
        /// </summary>
        /// 
        /// <param name="states">The number of states for this model.</param>
        /// <param name="symbols">The number of output symbols used for this model.</param>
        /// <param name="random">Whether to initialize emissions with random probabilities
        ///   or uniformly with <c>1 / number of symbols</c>. Default is false (default is
        ///   to use <c>1/symbols</c>).</param>
        /// 
        public static HiddenMarkovModel<GeneralDiscreteDistribution, int> CreateDiscrete(int states, int symbols, bool random)
        {
            return CreateDiscrete(new Ergodic(states, random), symbols, random);
        }



        #region IHiddenMarkovModel implementation
        int[] IHiddenMarkovModel.Decode(Array sequence, out double logLikelihood)
        {
            return this.Decode((int[])sequence, out logLikelihood);
        }

        double IHiddenMarkovModel.Evaluate(Array sequence)
        {
            return this.Evaluate((int[])sequence);
        }


        /// <summary>
        ///   Calculates the probability of each hidden state for each
        ///   observation in the observation vector.
        /// </summary>
        /// 
        /// <remarks>
        ///   If there are 3 states in the model, and the <paramref name="observations"/>
        ///   array contains 5 elements, the resulting vector will contain 5 vectors of
        ///   size 3 each. Each vector of size 3 will contain probability values that sum
        ///   up to one. By following those probabilities in order, we may decode those
        ///   probabilities into a sequence of most likely states. However, the sequence
        ///   of obtained states may not be valid in the model.
        /// </remarks>
        /// 
        /// <param name="observations">A sequence of observations.</param>
        /// 
        /// <returns>A vector of the same size as the observation vectors, containing
        ///  the probabilities for each state in the model for the current observation.
        ///  If there are 3 states in the model, and the <paramref name="observations"/>
        ///  array contains 5 elements, the resulting vector will contain 5 vectors of
        ///  size 3 each. Each vector of size 3 will contain probability values that sum
        ///  up to one.</returns>
        /// 
        double[][] IHiddenMarkovModel.Posterior(Array observations)
        {
            return this.Posterior((int[])observations);
        }

        /// <summary>
        ///   Calculates the probability of each hidden state for each observation 
        ///   in the observation vector, and uses those probabilities to decode the
        ///   most likely sequence of states for each observation in the sequence 
        ///   using the posterior decoding method. See remarks for details.
        /// </summary>
        /// 
        /// <remarks>
        ///   If there are 3 states in the model, and the <paramref name="observations"/>
        ///   array contains 5 elements, the resulting vector will contain 5 vectors of
        ///   size 3 each. Each vector of size 3 will contain probability values that sum
        ///   up to one. By following those probabilities in order, we may decode those
        ///   probabilities into a sequence of most likely states. However, the sequence
        ///   of obtained states may not be valid in the model.
        /// </remarks>
        /// 
        /// <param name="observations">A sequence of observations.</param>
        /// <param name="path">The sequence of states most likely associated with each
        ///   observation, estimated using the posterior decoding method.</param>
        /// 
        /// <returns>A vector of the same size as the observation vectors, containing
        ///  the probabilities for each state in the model for the current observation.
        ///  If there are 3 states in the model, and the <paramref name="observations"/>
        ///  array contains 5 elements, the resulting vector will contain 5 vectors of
        ///  size 3 each. Each vector of size 3 will contain probability values that sum
        ///  up to one.</returns>
        /// 
        double[][] IHiddenMarkovModel.Posterior(Array observations, out int[] path)
        {
            return this.Posterior((int[])observations, out path);
        }

        #endregion



        /// <summary>
        ///   Creates a new object that is a copy of the current instance.
        /// </summary>
        /// 
        /// <returns>
        ///   A new object that is a copy of this instance.
        /// </returns>
        /// 
        public override object Clone()
        {
            double[][] A = Vector.Copy(this.LogTransitions);
            double[][] B = Vector.Copy(this.LogEmissions);
            double[] pi = this.LogInitial.Copy();

            return new HiddenMarkovModel(A, B, pi, logarithm: true);
        }


        #region Save & Load methods

        /// <summary>
        ///   Saves the hidden Markov model to a stream.
        /// </summary>
        /// 
        /// <param name="stream">The stream to which the model is to be serialized.</param>
        /// 
        [Obsolete("Please use the Accord.Serializer class instead.")]
        public void Save(Stream stream)
        {
            BinaryFormatter b = new BinaryFormatter();
            b.Serialize(stream, this);
        }

        /// <summary>
        ///   Saves the hidden Markov model  to a stream.
        /// </summary>
        /// 
        /// <param name="path">The stream to which the model is to be serialized.</param>
        /// 
        [Obsolete("Please use the Accord.Serializer class instead.")]
        public void Save(string path)
        {
            using (FileStream fs = new FileStream(path, FileMode.Create))
            {
                this.Save(fs);
            }
        }

        /// <summary>
        ///   Loads a hidden Markov model  from a stream.
        /// </summary>
        /// 
        /// <param name="stream">The stream from which the model is to be deserialized.</param>
        /// 
        /// <returns>The deserialized classifier.</returns>
        /// 
        [Obsolete("Please use the Accord.Serializer class instead.")]
        public static HiddenMarkovModel Load(Stream stream)
        {
            BinaryFormatter b = new BinaryFormatter();
            return (HiddenMarkovModel)b.Deserialize(stream);
        }

        /// <summary>
        ///   Loads a hidden Markov model  from a file.
        /// </summary>
        /// 
        /// <param name="path">The path to the file from which the model is to be deserialized.</param>
        /// 
        /// <returns>The deserialized model.</returns>
        /// 
        [Obsolete("Please use the Accord.Serializer class instead.")]
        public static HiddenMarkovModel Load(string path)
        {
            using (FileStream fs = new FileStream(path, FileMode.Open))
            {
                return Load(fs);
            }
        }

        /// <summary>
        ///   Loads a hidden Markov model  from a stream.
        /// </summary>
        /// 
        /// <param name="stream">The stream from which the model is to be deserialized.</param>
        /// 
        /// <returns>The deserialized model.</returns>
        /// 
        [Obsolete("Please use the Accord.Serializer class instead.")]
        public static HiddenMarkovModel<TDistribution> Load<TDistribution>(Stream stream)
            where TDistribution : IDistribution
        {
            return HiddenMarkovModel<TDistribution>.Load(stream);
        }

        /// <summary>
        ///   Loads a hidden Markov model  from a file.
        /// </summary>
        /// 
        /// <param name="path">The path to the file from which the model is to be deserialized.</param>
        /// 
        /// <returns>The deserialized model.</returns>
        /// 
        [Obsolete("Please use the Accord.Serializer class instead.")]
        public static HiddenMarkovModel<TDistribution> Load<TDistribution>(string path)
            where TDistribution : IDistribution
        {
            return HiddenMarkovModel<TDistribution>.Load(path);
        }
#pragma warning restore 612, 618

        #endregion



        private void checkObservations(IndexOutOfRangeException ex, int[] observations)
        {
            for (int i = 0; i < observations.Length; i++)
            {
                if (observations[i] < 0 || observations[i] >= this.symbols)
                {
                    throw new ArgumentException("observations", "The observations vector must "
                    + "only contain values higher than or equal to 0 and less than " + this.symbols
                    + ". The value at the position " + i + " is " + observations[i] + ".", ex);
                }
            }
        }

    }
}
