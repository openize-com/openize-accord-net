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

namespace Openize.Accord.Statistics.Models.Fields.Potential_Functions.Clique_Factor_Potentials.Specialized
{
    using System;
    using Base;
    using Features;
    using Features.Base;
    using Features.Multivariate;
    using Openize.Accord.Core;

    /// <summary>
    ///   Factor Potential function for a Markov model whose states are independent 
    ///   distributions composed of discrete and Normal distributed components.
    /// </summary>
    /// 
    [Serializable]
    public class MarkovIndependentFactor : FactorPotential<double[]>
    {

        private int[][] stateTable;


        /// <summary>
        ///   Creates a new factor (clique) potential function.
        /// </summary>
        /// 
        /// <param name="owner">The owner <see cref="IPotentialFunction{T}"/>.</param>
        /// <param name="states">The number of states in this clique potential.</param>
        /// <param name="factorIndex">The index of this factor potential in the <paramref name="owner"/>.</param>
        /// <param name="stateTable">The lookup table of states where the independent distributions begin.</param>
        /// <param name="classIndex">The index of the first class label feature in the <paramref name="owner"/>'s parameter vector.</param>
        /// <param name="classCount">The number of class label features in this factor.</param>
        /// <param name="edgeIndex">The index of the first edge feature in the <paramref name="owner"/>'s parameter vector.</param>
        /// <param name="edgeCount">The number of edge features in this factor.</param>
        /// <param name="stateIndex">The index of the first state feature in the <paramref name="owner"/>'s parameter vector.</param>
        /// <param name="stateCount">The number of state features in this factor.</param>
        /// 
        public MarkovIndependentFactor(IPotentialFunction<double[]> owner, int states, int factorIndex, int[][] stateTable,
            int edgeIndex, int edgeCount,
            int stateIndex, int stateCount,
            int classIndex = 0, int classCount = 0)
            : base(owner, states, factorIndex, edgeIndex, edgeCount, stateIndex, stateCount, classIndex, classCount)
        {
            this.stateTable = stateTable;
        }


        /// <summary>
        ///   Computes the factor potential function for the given parameters.
        /// </summary>
        /// 
        /// <param name="previousState">The previous state in a given sequence of states.</param>
        /// <param name="currentState">The current state in a given sequence of states.</param>
        /// <param name="observations">The observation vector.</param>
        /// <param name="index">The index of the observation in the current state of the sequence.</param>
        /// <param name="outputClass">The output class label for the sequence.</param>
        /// <returns>The value of the factor potential function evaluated for the given parameters.</returns>
        /// 
        public override double Compute(int previousState, int currentState, double[][] observations, int index, int outputClass = 0)
        {
            // PS: This code seems messy because it should be as fast as possible. Unfortunately, 
            // avoiding (virtual) function calls is the main objective for this method to exist;
            // thus, refactoring those sections would defeat the existence of this method in the
            // first place. An alternative approach might be reconsidered once the project fully
            // migrates to .NET 4.5 and we could then use the explicit in-line method attributes.


            if (outputClass != this.Index)
                return Double.NegativeInfinity;

            double[] parameters = this.Owner.Weights;
            IFeature<double[]>[] features = this.Owner.Features;


            double sum = 0;

            // Output (class) probability
            if (this.OutputParameters.Count != 0 && previousState == -1)
            {
                int cindex = this.OutputParameters.Offset;
                double w = parameters[cindex];

                if (Double.IsNegativeInfinity(w))
                    return Double.NegativeInfinity;

                if (Double.IsNaN(w))
                    return parameters[cindex] = Double.NegativeInfinity;

                sum += w;
            }

            if (previousState == -1)
            {
                // Initial state probability (pi)
                int pindex = this.EdgeParameters.Offset + currentState;
                double pi = parameters[pindex];

                if (Double.IsNegativeInfinity(pi))
                    return Double.NegativeInfinity;

                if (Double.IsNaN(pi))
                    return parameters[pindex] = Double.NegativeInfinity;

                sum += pi;
            }
            else
            {
                // State transition probabilities (A)
                int aindex = this.EdgeParameters.Offset + this.States + previousState * this.States + currentState;
                double a = parameters[aindex];

                if (Double.IsNegativeInfinity(a))
                    return Double.NegativeInfinity;

                if (Double.IsNaN(a))
                    return parameters[aindex] = Double.NegativeInfinity;

                sum += a;
            }


            // State emission probability (B)
            int[] lookup = this.stateTable[currentState];
            double[] observation = observations[index];

            double B = 0;

            // For each dimension of the observation
            for (int j = 0; j < observation.Length; j++)
            {
                int bindex = this.StateParameters.Offset + lookup[j];

                var discrete = features[bindex] as MultivariateEmissionFeature;
                if (discrete != null)
                {
                    int i = (int)observations[index][j];
                    double b = parameters[bindex + i];

                    if (Double.IsNegativeInfinity(b))
                        return Double.NegativeInfinity;

                    if (Double.IsNaN(b))
                        return parameters[bindex + i] = Double.NegativeInfinity;

                    B += b;
                    continue;
                }

                var occupancy = features[bindex] as OccupancyFeature<double[]>;
                if (occupancy != null)
                {
                    // State occupancy feature
                    double b = parameters[bindex];

                    if (Double.IsNegativeInfinity(b))
                        return Double.NegativeInfinity;

                    if (Double.IsNaN(b))
                        return parameters[bindex] = Double.NegativeInfinity;

                    double u = observation[j];

                    int m1index = ++bindex;
                    int m2index = ++bindex;

                    if (u != 0)
                    {
                        // First moment feature
                        double m1 = parameters[m1index];

                        if (Double.IsNegativeInfinity(m1))
                            return Double.NegativeInfinity;

                        if (Double.IsNaN(m1))
                            m1 = parameters[m1index] = Double.NegativeInfinity;

                        b += m1 * u;

                        // Second moment feature
                        double m2 = parameters[m2index];

                        if (Double.IsNegativeInfinity(m2))
                            return Double.NegativeInfinity;

                        if (Double.IsNaN(m2))
                            m2 = parameters[m2index] = Double.NegativeInfinity;

                        b += m2 * u * u;
                    }

                    B += b;
                    continue;
                }

                throw new InvalidOperationException("The feature vector contains an unsupported"
                + " feature type. The MarkovIndependentFactor (which was likely created by a "
                + " MarkovMultivariateFunction) supports only discrete MultivariateEmissionFeature "
                + " (for GeneralDiscreteDistributions) and Occupancy/FirstMoment/SecondMoment features "
                + " (for Normal distributions).");
            }

            sum += B;

            Debug.Assert(!Double.IsNaN(sum));

            return sum;
        }

    }
}
