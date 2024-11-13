// AForge Math Library
// AForge.NET framework
// http://www.aforgenet.com/framework/
//
// Copyright © AForge.NET, 2007-2010
// contacts@aforgenet.com
//

namespace FileFormat.Accord.Math.AForge.Math.Metrics
{
    using System;
    using Distances;

    /// <summary>
    /// Please use Accord.Math.Distances.Cosine instead.
    /// </summary>
    [Obsolete("Please use Accord.Math.Distances.Cosine instead.")]
    public sealed class CosineSimilarity : ISimilarity
    {
        /// <summary>
        /// Please use Accord.Math.Distances.Cosine instead.
        /// </summary>
        public double GetSimilarityScore(double[] p, double[] q)
        {
            return new Cosine().Similarity(p, q);
        }
    }
}
