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

namespace Openize.Accord.Statistics.Analysis
{
    using System;
    using System.Collections.ObjectModel;
    using System.ComponentModel;
    using Accord.MachineLearning;
    using Accord.MachineLearning.Classifiers.Multiple.Multiclass;
    using Accord.MachineLearning.Learning;
    using Base;
    using Openize.Accord.Math.Core;
    using Openize.Accord.Math.Matrix;
    using Openize.Accord.Statistics.Accord.MachineLearning;
    using Math.Decompositions;
    using Models.Regression.Linear;
    using Openize.Accord.Math.Accord.Statistics;

    /// <summary>
    ///   Linear Discriminant Analysis (LDA).
    /// </summary>
    /// 
    /// <remarks>
    /// <para>
    ///   Linear Discriminant Analysis (LDA) is a method of finding such a linear
    ///   combination of variables which best separates two or more classes.</para>
    /// <para>
    ///   In itself LDA is not a classification algorithm, although it makes use of class
    ///   labels. However, the LDA result is mostly used as part of a linear classifier.
    ///   The other alternative use is making a dimension reduction before using nonlinear
    ///   classification algorithms.</para>
    /// <para>
    ///   It should be noted that several similar techniques (differing in requirements to the sample)
    ///   go together under the general name of Linear Discriminant Analysis. Described below is one of
    ///   these techniques with only two requirements:</para>  
    ///   <list type="number">
    ///     <item><description>The sample size shall exceed the number of variables, and </description></item>
    ///     <item><description>Classes may overlap, but their centers shall be distant from each other. </description></item>
    ///   </list>
    ///   
    /// <para>
    ///   Moreover, LDA requires the following assumptions to be true:</para>
    ///   <list type="bullet">
    ///     <item><description>Independent subjects.</description></item>
    ///     <item><description>Normality: the variance-covariance matrix of the
    ///     predictors is the same in all groups.</description></item>
    ///   </list>
    ///   
    /// <para>
    ///   If the latter assumption is violated, it is common to use quadratic discriminant analysis in
    ///   the same manner as linear discriminant analysis instead.</para>
    ///   
    /// <para>
    ///   This class can also be bound to standard controls such as the 
    ///   <a href="http://msdn.microsoft.com/en-us/library/system.windows.forms.datagridview.aspx">DataGridView</a>
    ///   by setting their DataSource property to the analysis' <see cref="BaseDiscriminantAnalysis.Discriminants"/> property.</para>
    ///   
    /// <para>
    ///    References:
    ///    <list type="bullet">
    ///      <item><description>
    ///        R. Gutierrez-Osuna, Linear Discriminant Analysis. Available on:
    ///        http://research.cs.tamu.edu/prism/lectures/pr/pr_l10.pdf </description></item>
    ///     </list></para>     
    /// </remarks>
    /// 
    /// <example>
    /// <code source="Unit Tests\Accord.Tests.Statistics\Analysis\LinearDiscriminantAnalysisTest.cs" region="doc_learn" />
    /// </example>
    /// 
    [Serializable]
#pragma warning disable 612, 618
    public class LinearDiscriminantAnalysis : BaseDiscriminantAnalysis,
        IDiscriminantAnalysis, IProjectionAnalysis,
        ISupervisedLearning<LinearDiscriminantAnalysis.Pipeline, double[], int>
#pragma warning restore 612, 618
    {

        /// <summary>
        ///   Gets a classification pipeline that can be used to classify
        ///   new samples into one of the <see cref="BaseDiscriminantAnalysis.NumberOfClasses"/> 
        ///   learned in this discriminant analysis. This pipeline is
        ///   only available after a call to the <see cref="Learn"/> method.
        /// </summary>
        /// 
        public Pipeline Classifier { get; private set; }

#pragma warning disable 612, 618
        /// <summary>
        ///   Constructs a new Linear Discriminant Analysis object.
        /// </summary>
        /// 
        /// <param name="inputs">The source data to perform analysis. The matrix should contain
        /// variables as columns and observations of each variable as rows.</param>
        /// <param name="outputs">The labels for each observation row in the input matrix.</param>
        /// 
        [Obsolete("Please pass the 'inputs' and 'outputs' parameters to the Learn method instead.")]
        public LinearDiscriminantAnalysis(double[,] inputs, int[] outputs)
        {
            this.init(inputs, outputs);
            this.Threshold = 0;
        }

        /// <summary>
        ///   Constructs a new Linear Discriminant Analysis object.
        /// </summary>
        /// 
        /// <param name="inputs">The source data to perform analysis. The matrix should contain
        /// variables as columns and observations of each variable as rows.</param>
        /// <param name="outputs">The labels for each observation row in the input matrix.</param>
        /// 
        [Obsolete("Please pass the 'inputs' and 'outputs' parameters to the Learn method instead.")]
        public LinearDiscriminantAnalysis(double[][] inputs, int[] outputs)
        {
            this.init(inputs.ToMatrix(), outputs);
            this.Threshold = 0;
        }
#pragma warning disable 612, 618

        /// <summary>
        ///   Constructs a new Linear Discriminant Analysis object.
        /// </summary>
        /// 
        public LinearDiscriminantAnalysis()
        {
            this.Threshold = 0;
        }

        /// <summary>
        ///   Computes the Multi-Class Linear Discriminant Analysis algorithm.
        /// </summary>
        /// 
        [Obsolete("Please use the Learn method instead.")]
        public virtual void Compute()
        {
            this.Learn(this.Source.ToJagged(), this.Classifications);

            this.Result = this.Source.DotWithTransposed(this.DiscriminantVectors).ToMatrix();
        }

        /// <summary>
        /// Applies the transformation to an input, producing an associated output.
        /// </summary>
        /// <param name="input">The input data to which the transformation should be applied.</param>
        /// <param name="result">A location to store the output, avoiding unnecessary memory allocations.</param>
        /// <returns>
        /// The output generated by applying this transformation to the given input.
        /// </returns>
        public override double[][] Transform(double[][] input, double[][] result)
        {
            for (int i = 0; i < input.Length; i++)
                for (int j = 0; j < result[i].Length; j++)
                    for (int k = 0; k < input[i].Length; k++)
                        result[i][j] += input[i][k] * this.DiscriminantVectors[j][k];
            return result;
        }





        /// <summary>
        /// Learns a model that can map the given inputs to the given outputs.
        /// </summary>
        /// <param name="x">The model inputs.</param>
        /// <param name="y">The desired outputs associated with each <paramref name="x">inputs</paramref>.</param>
        /// <param name="weights">The weight of importance for each input-output pair (if supported by the learning algorithm).</param>
        /// <returns>
        /// A model that has learned how to produce <paramref name="y" /> given <paramref name="x" />.
        /// </returns>
        public Pipeline Learn(double[][] x, int[] y, double[] weights = null)
        {
            if (weights != null)
                throw new ArgumentException(global::Accord.Properties.Resources.NotSupportedWeights, "weights");

            this.Init(x, y);

            // Compute entire data set measures
            this.Means = Measures.Mean(x, dimension: 0);
            this.StandardDeviations = Measures.StandardDeviation(x, this.Means);

            // Initialize the scatter matrices
            var Sw = Jagged.Zeros(this.NumberOfInputs, this.NumberOfInputs);
            var Sb = Jagged.Zeros(this.NumberOfInputs, this.NumberOfInputs);

            // For each class
            for (int c = 0; c < this.Classes.Count; c++)
            {
                int[] idx = Matrix.Find(y, y_i => y_i == this.Classes[c].Number);

                // Get the class subset
                double[][] subset = x.Get(idx);
                int count = subset.GetLength(0);

                // Get the class mean
                double[] mean = Measures.Mean(subset, dimension: 0);

                // Continue constructing the Within-Class Scatter Matrix
                double[][] Swi = Measures.Scatter(subset, mean, (double)count);

                Sw.Add(Swi, result: Sw); // Sw = Sw + Swi

                // Continue constructing the Between-Class Scatter Matrix
                double[] d = mean.Subtract(this.Means);
                double[][] Sbi = Jagged.Outer(d, d);
                Sbi.Multiply((double)this.NumberOfInputs, result: Sbi);

                Sb.Add(Sbi, result: Sb); // Sb = Sb + Sbi

                // Store some additional information
                this.classScatter[c] = Swi;
                this.classCount[c] = count;
                this.classMeans[c] = mean;
                this.classStdDevs[c] = Measures.StandardDeviation(subset, mean);
            }


            // Compute the generalized eigenvalue decomposition
            var gevd = new JaggedGeneralizedEigenvalueDecomposition(Sb, Sw, sort: true);

            // Get the eigenvalues and corresponding eigenvectors
            double[] evals = gevd.RealEigenvalues;
            double[][] eigs = gevd.Eigenvectors;

            // Eliminate unwanted components
            int nonzero = x.Columns();
            if (this.Threshold > 0)
                nonzero = Math.Min(gevd.Rank, GetNonzeroEigenvalues(evals, this.Threshold));
            if (this.NumberOfInputs != 0)
                nonzero = Math.Min(nonzero, this.NumberOfInputs);
            if (this.NumberOfOutputs != 0)
                nonzero = Math.Min(nonzero, this.NumberOfOutputs);

            eigs = eigs.Get(null, 0, nonzero);
            evals = evals.Get(0, nonzero);

            // Store information
            this.Eigenvalues = evals;
            this.DiscriminantVectors = eigs.Transpose();
            base.ScatterBetweenClass = Sb;
            base.ScatterWithinClass = Sw;
            base.NumberOfInputs = x.Columns();
            base.NumberOfOutputs = evals.Length;

            // Compute feature space means for later classification
            for (int c = 0; c < this.projectedMeans.Length; c++)
                this.projectedMeans[c] = this.classMeans[c].Dot(eigs);

            // Computes additional information about the analysis and creates the
            //  object-oriented structure to hold the discriminants found.
            this.CreateDiscriminants();

            this.Classifier = this.CreateClassifier();

            return this.Classifier;
        }

        private Pipeline CreateClassifier()
        {
            double[][] eig = this.DiscriminantVectors.Transpose().Get(null, 0, this.NumberOfOutputs);

            return new Pipeline()
            {
                NumberOfInputs = this.NumberOfInputs,
                NumberOfOutputs = this.NumberOfClasses,
                NumberOfClasses = this.NumberOfClasses,
                First = new MultivariateLinearRegression()
                {
                    Weights = eig,
                    NumberOfInputs = this.NumberOfInputs,
                    NumberOfOutputs = this.NumberOfOutputs,
                },
                Second = new MinimumMeanDistanceClassifier()
                {
                    Means = this.projectedMeans.Get(null, 0, this.NumberOfOutputs),
                    NumberOfInputs = this.NumberOfOutputs,
                    NumberOfOutputs = this.NumberOfClasses,
                    NumberOfClasses = this.NumberOfClasses
                },
            };
        }

        /// <summary>Transform
        ///   Classifies a new instance into one of the available classes.
        /// </summary>
        /// 
        [Obsolete("Please use Classifier.Decide() instead.")]
        public override int Classify(double[] input)
        {
            return this.Classes[this.Classifier.Decide(input)].Number;
        }

        /// <summary>
        ///   Classifies a new instance into one of the available classes.
        /// </summary>
        /// 
        [Obsolete("Please use Classifier.Decide() or Classifier.Scores() instead.")]
        public override int Classify(double[] input, out double[] responses)
        {
            int decision;
            responses = this.Classifier.Scores(input, out decision);
            return this.Classes[decision].Number;
        }

        /// <summary>
        ///   Classifies new instances into one of the available classes.
        /// </summary>
        /// 
        [Obsolete("Please use Classifier.Decide() instead.")]
        public override int[] Classify(double[][] inputs)
        {
            int[] result = this.Classifier.Decide(inputs);
            for (int i = 0; i < result.Length; i++)
                result[i] = this.Classes[result[i]].Number;
            return result;
        }

        /// <summary>
        ///   Gets the output of the discriminant function for a given class.
        /// </summary>
        /// 
        public override double DiscriminantFunction(double[] input, int classIndex)
        {
            return this.Classifier.Score(input, classIndex);
        }

        /// <summary>
        ///   Standard regression and classification pipeline for <see cref="LinearDiscriminantAnalysis"/>.
        /// </summary>
        /// 
        [Serializable]
        public sealed class Pipeline : MulticlassScoreClassifierBase<double[]>
        {
            /// <summary>
            /// Gets or sets the first step in the pipeline.
            /// </summary>
            /// 
            public MultivariateLinearRegression First { get; set; }

            /// <summary>
            /// Gets or sets the second step in the pipeline.
            /// </summary>
            /// 
            public MinimumMeanDistanceClassifier Second { get; set; }

            /// <summary>
            /// Computes a numerical score measuring the association between
            /// the given <paramref name="input" /> vector and each class.
            /// </summary>
            /// <param name="input">The input vector.</param>
            /// <param name="result">An array where the result will be stored,
            /// avoiding unnecessary memory allocations.</param>
            /// <returns></returns>
            public override double[][] Scores(double[][] input, double[][] result)
            {
                return this.Second.Scores(this.First.Transform(input), result);
            }

            /// <summary>
            /// Computes a numerical score measuring the association between
            /// the given <paramref name="input" /> vector and a given
            /// <paramref name="classIndex" />.
            /// </summary>
            /// <param name="input">The input vector.</param>
            /// <param name="classIndex">The index of the class whose score will be computed.</param>
            /// <returns>System.Double.</returns>
            public override double Score(double[] input, int classIndex)
            {
                return this.Second.Score(this.First.Transform(input), classIndex);
            }
        }
    }

    #region Support Classes
    /// <summary>
    ///   Represents a class found during Discriminant Analysis, allowing it to
    ///   be bound to controls like the DataGridView.
    ///   
    ///   This class cannot be instantiated.
    /// </summary>
    /// 
    [Serializable]
    public class DiscriminantAnalysisClass
    {
        private BaseDiscriminantAnalysis analysis;
        private int classNumber;
        private int index;

        /// <summary>
        ///   Creates a new Class representation
        /// </summary>
        /// 
        internal DiscriminantAnalysisClass(BaseDiscriminantAnalysis analysis, int index, int classNumber)
        {
            this.analysis = analysis;
            this.index = index;
            this.classNumber = classNumber;
        }

        /// <summary>
        ///   Gets the Index of this class on the original analysis collection.
        /// </summary>
        /// 
        public int Index
        {
            get { return this.index; }
        }

        /// <summary>
        ///   Gets the number labeling this class.
        /// </summary>
        /// 
        public int Number
        {
            get { return this.classNumber; }
        }

        /// <summary>
        ///   Gets the prevalence of the class on the original data set.
        /// </summary>
        /// 
        public double Prevalence
        {
            get { return (double)this.Count / this.analysis.NumberOfSamples; }
        }

        /// <summary>
        ///   Gets the class' mean vector.
        /// </summary>
        /// 
        public double[] Mean
        {
            get { return this.analysis.classMeans[this.index]; }
        }

        /// <summary>
        ///   Gets the feature-space means of the projected data.
        /// </summary>
        /// 
        public double[] ProjectionMean
        {
            get { return this.analysis.projectedMeans[this.index]; }
        }

        /// <summary>
        ///   Gets the class' standard deviation vector.
        /// </summary>
        /// 
        public double[] StandardDeviation
        {
            get { return this.analysis.classStdDevs[this.index]; }
        }

        /// <summary>
        ///   Gets the Scatter matrix for this class.
        /// </summary>
        /// 
        public double[][] Scatter
        {
            get { return this.analysis.classScatter[this.index]; }
        }

#pragma warning disable 612, 618
        /// <summary>
        ///   Gets the indices of the rows in the original data which belong to this class.
        /// </summary>
        /// 
        [Obsolete("This property will be removed.")]
        public int[] Indices
        {
            get { return Matrix.Find(this.analysis.Classifications, y => y == this.classNumber); }
        }

        /// <summary>
        ///   Gets the subset of the original data spawned by this class.
        /// </summary>
        /// 
        [Obsolete("This property will be removed.")]
        public double[,] Subset
        {

            get { return this.analysis.Source.Submatrix(this.Indices); }
        }
#pragma warning restore 612, 618

        /// <summary>
        ///   Gets the number of observations inside this class.
        /// </summary>
        /// 
        public int Count
        {
            get { return this.analysis.classCount[this.index]; }
        }

        /// <summary>
        ///   Discriminant function for the class.
        /// </summary>
        /// 
        public double DiscriminantFunction(double[] projection)
        {
            return this.analysis.DiscriminantFunction(projection, classIndex: this.index);
        }
    }

    /// <summary>
    /// <para>
    ///   Represents a discriminant factor found during Discriminant Analysis,
    ///   allowing it to be bound to controls like the DataGridView.</para>
    /// <para>
    ///   This class cannot be instantiated.</para>  
    /// </summary>
    /// 
    [Serializable]
    public class Discriminant : IAnalysisComponent
    {
        private BaseDiscriminantAnalysis analysis;
        private int index;

        /// <summary>
        ///   Creates a new discriminant factor representation.
        /// </summary>
        /// 
        internal Discriminant(BaseDiscriminantAnalysis analysis, int index)
        {
            this.analysis = analysis;
            this.index = index;
        }

        /// <summary>
        ///   Gets the index of this discriminant factor
        ///   on the original analysis collection.
        /// </summary>
        /// 
        public int Index
        {
            get { return this.index; }
        }

        /// <summary>
        ///   Gets the Eigenvector for this discriminant factor.
        /// </summary>
        /// 
        public double[] Eigenvector
        {
            get { return this.analysis.DiscriminantVectors[this.index]; }
        }

        /// <summary>
        ///   Gets the Eigenvalue for this discriminant factor.
        /// </summary>
        /// 
        public double Eigenvalue
        {
            get { return this.analysis.Eigenvalues[this.index]; }
        }

        /// <summary>
        ///   Gets the proportion, or amount of information explained by this discriminant factor.
        /// </summary>
        /// 
        public double Proportion
        {
            get { return this.analysis.DiscriminantProportions[this.index]; }
        }

        /// <summary>
        ///   Gets the cumulative proportion of all discriminant factors until this factor.
        /// </summary>
        /// 
        [DisplayName("Cumulative")]
        public double CumulativeProportion
        {
            get { return this.analysis.CumulativeProportions[this.index]; }
        }

    }

    /// <summary>
    /// <para>
    ///   Represents a collection of Discriminants factors found in the Discriminant Analysis.</para>
    /// <para>
    ///   This class cannot be instantiated.</para>
    /// </summary>
    /// 
    [Serializable]
    public class DiscriminantCollection : ReadOnlyCollection<Discriminant>
    {
        internal DiscriminantCollection(Discriminant[] components)
            : base(components)
        {
        }
    }

    /// <summary>
    /// <para>
    ///   Represents a collection of classes found in the Discriminant Analysis.</para>
    /// <para>
    ///   This class cannot be instantiated.</para>  
    /// </summary>
    /// 
    [Serializable]
    public class DiscriminantAnalysisClassCollection : ReadOnlyCollection<DiscriminantAnalysisClass>
    {
        internal DiscriminantAnalysisClassCollection(DiscriminantAnalysisClass[] components)
            : base(components) { }
    }

    #endregion

}


