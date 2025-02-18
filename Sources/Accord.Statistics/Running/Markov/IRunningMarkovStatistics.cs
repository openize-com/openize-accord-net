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

namespace Openize.Accord.Statistics.Running.Markov
{
    using System;
    using Openize.Accord.Math;
    using Openize.Accord.Math.Matrix;
    using Math;
    using Openize.Accord.Statistics.Models.Markov.Base;

    /// <summary>
    ///   Common interface for running Markov filters.
    /// </summary>
    public interface IRunningMarkovStatistics
    {
        /// <summary>
        ///   Gets whether the model has been initialized or not.
        /// </summary>
        /// 
        bool Started { get; set; }

         /// <summary>
        ///   Gets the current vector of probabilities of being in each state.
        /// </summary>
        /// 
        double[] Current { get; }

        /// <summary>
        ///   Gets the current most likely state (in the Viterbi path).
        /// </summary>
        /// 
        int CurrentState { get; }

        /// <summary>
        ///   Gets the current Viterbi probability
        ///   (along the most likely path).
        /// </summary>
        /// 
        double LogViterbi { get; }

        /// <summary>
        ///   Gets the current Forward probability
        ///   (along all possible paths).
        /// </summary>
        /// 
        double LogForward { get; }

        /// <summary>
        ///   Clears all measures previously computed
        ///   and indicate the sequence has ended.
        /// </summary>
        /// 
        void Clear();

    }

    /// <summary>
    ///   Base class for running hidden Markov filters.
    /// </summary>
    /// 
    [Serializable]
    internal class BaseRunningMarkovStatistics
    {

        private bool started;
        private double[] current;
        private int? currentState;

        private double? logViterbi;
        private double? logForward;


        /// <summary>
        ///   Initializes a new instance of the <see cref="BaseRunningMarkovStatistics"/> class.
        /// </summary>
        /// 
        /// <param name="model">The Markov model.</param>
        /// 
        public BaseRunningMarkovStatistics(IHiddenMarkovModel model)
        {
            if (model == null)
                throw new ArgumentNullException("model");

            this.current = new double[model.States];
        }

        /// <summary>
        ///   Gets whether the model has been initialized or not.
        /// </summary>
        /// 
        public bool Started
        {
            get { return this.started; }
            set { this.started = value; }
        }

        /// <summary>
        ///   Gets the current vector of probabilities of being in each state.
        /// </summary>
        /// 
        public double[] Current
        {
            get { return this.current; }
        }

        /// <summary>
        ///   Gets the current most likely state (in the Viterbi path).
        /// </summary>
        /// 
        public int CurrentState
        {
            get
            {
                if (this.currentState == null)
                {
                    int state; this.current.Max(out state);
                    this.currentState = state;
                }

                return this.currentState.Value;
            }
        }

        /// <summary>
        ///   Gets the current Viterbi probability
        ///   (along the most likely path).
        /// </summary>
        /// 
        public double LogViterbi
        {
            get
            {
                if (this.logViterbi == null)
                {
                    int state = this.CurrentState;
                    this.logViterbi = this.Current[state];
                }

                return this.logViterbi.Value;
            }
        }

        /// <summary>
        ///   Gets the current Forward probability
        ///   (along all possible paths).
        /// </summary>
        /// 
        public double LogForward
        {
            get
            {
                if (this.logForward == null)
                {
                    double sum = Double.NegativeInfinity;
                    for (int i = 0; i < this.Current.Length; i++)
                        sum = Special.LogSum(sum, this.Current[i]);
                    this.logForward = sum;
                }

                return this.logForward.Value;
            }
        }

        /// <summary>
        ///   Clears all measures previously computed
        ///   and indicate the sequence has ended.
        /// </summary>
        /// 
        public void Clear()
        {
            this.Started = false;
            for (int i = 0; i < this.current.Length; i++)
                this.current[i] = 0;

            this.currentState = null;
            this.logViterbi = null;
            this.logForward = null;
        }

        /// <summary>
        ///   Clears all measures previously computed.
        /// </summary>
        /// 
        public void Invalidate()
        {
            this.currentState = null;
            this.logViterbi = null;
            this.logForward = null;
        }

    }
}
