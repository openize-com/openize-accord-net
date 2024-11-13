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

namespace FileFormat.Accord.Statistics.Analysis
{
    using System;
    using Accord.MachineLearning;
    using Accord.MachineLearning.Classifiers.Multiple.Multiclass;
    using Accord.MachineLearning.Learning;
    using FileFormat.Accord.Core.Exceptions;
    using Kernels;
    using Kernels.Base;
    using Math.Accord.Statistics;
    using Math.Comparers;
    using Math.Core;
    using Math.Decompositions;
    using Models.Regression.Nonlinear;
    using Math.Matrix;

    /// <summary>
    ///   Kernel (Fisher) Discriminant Analysis.
    /// </summary>
    /// 
    /// <remarks>
    /// <para>
    ///   Kernel (Fisher) discriminant analysis (kernel FDA) is a non-linear generalization
    ///   of linear discriminant analysis (LDA) using techniques of kernel methods. Using a
    ///   kernel, the originally linear operations of LDA are done in a reproducing kernel
    ///   Hilbert space with a non-linear mapping.</para>
    /// <para>
    ///   The algorithm used is a multi-class generalization of the original algorithm by
    ///   Mika et al. in Fisher discriminant analysis with kernels (1999).</para>  
    ///   
    /// <para>
    ///   This class can also be bound to standard controls such as the 
    ///   <a href="http://msdn.microsoft.com/en-us/library/system.windows.forms.datagridview.aspx">DataGridView</a>
    ///   by setting their DataSource property to the analysis' <see cref="BaseDiscriminantAnalysis.Discriminants"/> property.</para>
    ///   
    /// <para>
    ///   References:
    ///   <list type="bullet">
    ///     <item><description>
    ///       Mika et al, Fisher discriminant analysis with kernels (1999). Available on
    ///       <a href="http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.35.9904">
    ///       http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.35.9904 </a></description></item>
    ///  </list></para>  
    /// </remarks>
    /// 
    /// <example>
    /// <code source="Unit Tests\Accord.Tests.Statistics\Analysis\KernelDiscriminantAnalysisTest.cs" region="doc_learn" />
    /// </example>
    /// 
    [Serializable]
    public class KernelDiscriminantAnalysis : BaseDiscriminantAnalysis,
        ISupervisedLearning<KernelDiscriminantAnalysis.Pipeline, double[], int>
    {
        private IKernel kernel;
        private double regularization = 1e-4;
        private double[][] input;

        /// <summary>
        ///   Gets a classification pipeline that can be used to classify
        ///   new samples into one of the <see cref="BaseDiscriminantAnalysis.NumberOfClasses"/> 
        ///   learned in this discriminant analysis. This pipeline is
        ///   only available after a call to the <see cref="Learn"/> method.
        /// </summary>
        /// 
        public Pipeline Classifier { get; private set; }

        /// <summary>
        ///   Gets or sets the matrix of original values used to create
        ///   this analysis. Those values are required to build kernel 
        ///  (Gram) matrices when classifying new samples.
        /// </summary>
        /// 
        public double[][] Input
        {
            get { return this.input; }
            set { this.input = value; }
        }

        /// <summary>
        ///   Constructs a new Kernel Discriminant Analysis object.
        /// </summary>
        /// 
        /// <param name="inputs">The source data to perform analysis. The matrix should contain
        /// variables as columns and observations of each variable as rows.</param>
        /// <param name="output">The labels for each observation row in the input matrix.</param>
        /// <param name="kernel">The kernel to be used in the analysis.</param>
        /// 
        [Obsolete("Please pass the 'inputs' and 'outputs' parameters to the Learn method instead.")]
        public KernelDiscriminantAnalysis(double[,] inputs, int[] output, IKernel kernel)
            : this(kernel)
        {
            this.init(inputs, output);
        }

        /// <summary>
        ///   Constructs a new Kernel Discriminant Analysis object.
        /// </summary>
        /// 
        /// <param name="inputs">The source data to perform analysis. The matrix should contain
        /// variables as columns and observations of each variable as rows.</param>
        /// <param name="output">The labels for each observation row in the input matrix.</param>
        /// <param name="kernel">The kernel to be used in the analysis.</param>
        /// 
        [Obsolete("Please pass the 'inputs' and 'outputs' parameters to the Learn method instead.")]
        public KernelDiscriminantAnalysis(double[][] inputs, int[] output, IKernel kernel)
            : this(kernel)
        {
            this.init(inputs.ToMatrix(), output);
        }

        /// <summary>
        ///   Constructs a new Kernel Discriminant Analysis object.
        /// </summary>
        /// 
        public KernelDiscriminantAnalysis(IKernel kernel)
        {
            this.kernel = kernel;
            this.Threshold = 0;
        }

        /// <summary>
        ///   Constructs a new Kernel Discriminant Analysis object.
        /// </summary>
        /// 
        public KernelDiscriminantAnalysis()
            : this(new Linear())
        {

        }

        /// <summary>
        ///   Gets or sets the Kernel used in the analysis.
        /// </summary>
        /// 
        public IKernel Kernel
        {
            get { return this.kernel; }
            set { this.kernel = value; }
        }

        /// <summary>
        ///   Gets or sets the regularization parameter to
        ///   avoid non-singularities at the solution.
        /// </summary>
        /// 
        public double Regularization
        {
            get { return this.regularization; }
            set
            {
                if (value < 0)
                    throw new ArgumentOutOfRangeException("value", "Value must be positive.");
                this.regularization = value;
            }
        }






        /// <summary>
        ///   Computes the Multi-Class Kernel Discriminant Analysis algorithm.
        /// </summary>
        /// 
        [Obsolete("Please use Learn(x, y) instead.")]
        public void Compute()
        {
            this.input = this.Source.ToJagged();

            // Get some initial information
            int dimension = this.Source.GetLength(0);
            double total = dimension;

            // Create the Gram (Kernel) Matrix
            double[,] K = new double[dimension, dimension];
            for (int i = 0; i < dimension; i++)
            {
                double[] row = this.Source.GetRow(i);
                for (int j = i; j < dimension; j++)
                {
                    double s = this.kernel.Function(row, this.Source.GetRow(j));
                    K[i, j] = s; // Assume K will be symmetric
                    K[j, i] = s;
                }
            }


            // Compute entire data set measures
            base.Means = Measures.Mean(K, dimension: 0);
            base.StandardDeviations = Measures.StandardDeviation(K, this.Means);


            // Initialize the kernel analogous scatter matrices
            double[,] Sb = new double[dimension, dimension];
            double[,] Sw = new double[dimension, dimension];


            // For each class
            for (int c = 0; c < this.Classes.Count; c++)
            {
                // Get the Kernel matrix class subset
                double[,] Kc = K.Submatrix(Matrix.Find(this.Classifications, y_i => y_i == this.Classes[c].Number));
                int count = Kc.GetLength(0);

                // Get the Kernel matrix class mean
                double[] mean = Measures.Mean(Kc, dimension: 0);

                // Construct the Kernel equivalent of the Within-Class Scatter matrix
                double[,] Swi = Measures.Scatter(Kc, mean, (double)count);

                // Sw = Sw + Swi
                for (int i = 0; i < dimension; i++)
                    for (int j = 0; j < dimension; j++)
                        Sw[i, j] += Swi[i, j];

                // Construct the Kernel equivalent of the Between-Class Scatter matrix
                double[] d = mean.Subtract(base.Means);
                double[,] Sbi = Matrix.Outer(d, d).Multiply(total);

                // Sb = Sb + Sbi
                for (int i = 0; i < dimension; i++)
                    for (int j = 0; j < dimension; j++)
                        Sb[i, j] += Sbi[i, j];

                // Store additional information
                base.ClassScatter[c] = Swi.ToJagged();
                base.ClassCount[c] = count;
                base.ClassMeans[c] = mean;
                base.ClassStandardDeviations[c] = Measures.StandardDeviation(Kc, mean);
            }


            // Add regularization to avoid singularity
            for (int i = 0; i < dimension; i++)
                Sw[i, i] += this.regularization;


            // Compute the generalized eigenvalue decomposition
            var gevd = new GeneralizedEigenvalueDecomposition(Sb, Sw);

            if (gevd.IsSingular) // check validity of the results
            {
                throw new SingularMatrixException("One of the matrices is singular. Please retry " +
                    "the method with a higher regularization constant.");
            }


            // Get the eigenvalues and corresponding eigenvectors
            double[] evals = gevd.RealEigenvalues;
            double[,] eigs = gevd.Eigenvectors;

            // Sort eigenvalues and vectors in descending order
            eigs = Matrix.Sort(evals, eigs, new GeneralComparer(ComparerDirection.Descending, true));


            if (this.Threshold > 0)
            {
                // We will be discarding less important
                // eigenvectors to conserve memory.

                // Calculate component proportions
                double sum = 0.0; // total variance
                for (int i = 0; i < dimension; i++)
                    sum += Math.Abs(evals[i]);

                if (sum > 0)
                {
                    int keep = 0;

                    // Now we will detect how many components we have
                    //  have to keep in order to achieve the level of
                    //  explained variance specified by the threshold.

                    while (keep < dimension)
                    {
                        // Get the variance explained by the component
                        double explainedVariance = Math.Abs(evals[keep]);

                        // Check its proportion
                        double proportion = explainedVariance / sum;

                        // Now, if the component explains an
                        // enough proportion of the variance,
                        if (proportion > this.Threshold)
                            keep++; // We can keep it.
                        else
                            break;  // Otherwise we can stop, since the
                        // components are ordered by variance.
                    }

                    if (keep > 0)
                    {
                        // Resize the vectors keeping only needed components
                        eigs = eigs.Submatrix(0, dimension - 1, 0, keep - 1);
                        evals = evals.Submatrix(0, keep - 1);
                    }
                    else
                    {
                        // No component will be kept.
                        eigs = new double[dimension, 0];
                        evals = new double[0];
                    }
                }
            }

            // Store information
            base.Eigenvalues = evals;
            base.DiscriminantVectors = eigs.ToJagged().Transpose();
            base.ScatterBetweenClass = Sb.ToJagged();
            base.ScatterWithinClass = Sw.ToJagged();
            this.NumberOfOutputs = eigs.Columns();
            this.NumberOfClasses = this.Classes.Count;

            // Project into the kernel discriminant space
            this.Result = Matrix.Dot(K, eigs);

            // Compute feature space means for later classification
            for (int c = 0; c < this.Classes.Count; c++)
                this.ProjectionMeans[c] = this.ClassMeans[c].Dot(eigs);

            // Computes additional information about the analysis and creates the
            //  object-oriented structure to hold the discriminants found.
            this.CreateDiscriminants();

            this.Classifier = this.CreateClassifier();
        }

        private Pipeline CreateClassifier()
        {
            if (this.NumberOfOutputs == 0)
                return null;

            double[][] eig = this.DiscriminantVectors;

            return new Pipeline()
            {
                NumberOfInputs = this.NumberOfInputs,
                NumberOfOutputs = this.NumberOfClasses,
                NumberOfClasses = this.NumberOfClasses,

                First = new MultivariateKernelRegression()
                {
                    Weights = eig,
                    BasisVectors = this.input,
                    Kernel = this.Kernel,
                    NumberOfInputs = this.NumberOfInputs,
                    NumberOfOutputs = this.NumberOfOutputs,
                },
                Second = new MinimumMeanDistanceClassifier()
                {
                    Means = this.projectedMeans,
                    NumberOfInputs = this.NumberOfOutputs,
                    NumberOfOutputs = this.NumberOfClasses,
                    NumberOfClasses = this.NumberOfClasses,
                }
            };
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
            // TODO: Do without forming the kernel matrix
            double[][] K = this.kernel.ToJagged2(x: input, y: this.input);
            // return K.DotWithTransposed(DiscriminantVectors);

            for (int i = 0; i < input.Length; i++)
                for (int j = 0; j < result[i].Length; j++)
                    for (int k = 0; k < K[i].Length; k++)
                        result[i][j] += K[i][k] * this.DiscriminantVectors[j][k];
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
        /// 
        public Pipeline Learn(double[][] x, int[] y, double[] weights = null)
        {
            if (weights != null)
                throw new ArgumentException(global::Accord.Properties.Resources.NotSupportedWeights, "weights");

            this.Init(x, y);

            // Create the Gram (Kernel) Matrix
            var K = this.kernel.ToJagged(x);

            // Compute entire data set measures
            base.Means = Measures.Mean(K, dimension: 0);
            base.StandardDeviations = Measures.StandardDeviation(K, this.Means);

            // Initialize the kernel analogous scatter matrices
            //int dimension = x.Columns();
            double[][] Sb = Jagged.Zeros(this.NumberOfSamples, this.NumberOfSamples);
            double[][] Sw = Jagged.Zeros(this.NumberOfSamples, this.NumberOfSamples);

            // For each class
            for (int c = 0; c < this.Classes.Count; c++)
            {
                var idx = Matrix.Find(y, y_i => y_i == c);

                // Get the Kernel matrix class subset
                double[][] Kc = K.Get(idx);
                int count = Kc.Rows();

                // Get the Kernel matrix class mean
                double[] mean = Measures.Mean(Kc, dimension: 0);

                // Construct the Kernel equivalent of the Within-Class Scatter matrix
                double[][] Swi = Measures.Scatter(Kc, dimension: 0, means: mean);
                Swi.Divide((double)count, result: Swi);
                Sw.Add(Swi, result: Sw); // Sw = Sw + Swi

                // Construct the Kernel equivalent of the Between-Class Scatter matrix
                double[] d = mean.Subtract(base.Means);
                double[][] Sbi = Jagged.Outer(d, d);
                Sbi.Multiply((double)this.NumberOfSamples, result: Sbi);

                Sb.Add(Sbi, result: Sb); // Sb = Sb + Sbi

                // Store additional information
                base.ClassScatter[c] = Swi;
                base.ClassCount[c] = count;
                base.ClassMeans[c] = mean;
                base.ClassStandardDeviations[c] = Measures.StandardDeviation(Kc, mean);
            }

            // Add regularization to avoid singularity
            Sw.AddToDiagonal(this.regularization, result: Sw);

            // Compute the generalized eigenvalue decomposition
            var gevd = new JaggedGeneralizedEigenvalueDecomposition(Sb, Sw, sort: true);

            if (gevd.IsSingular) // check validity of the results
            {
                throw new SingularMatrixException("One of the matrices is singular. Please retry " +
                    "the method with a higher regularization constant.");
            }

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
            this.input = x;
            base.Eigenvalues = evals;
            base.DiscriminantVectors = eigs.Transpose();
            base.ScatterBetweenClass = Sb;
            base.ScatterWithinClass = Sw;
            base.NumberOfOutputs = evals.Length;

            // Compute feature space means for later classification
            for (int c = 0; c < this.Classes.Count; c++)
                this.ProjectionMeans[c] = this.ClassMeans[c].Dot(eigs);

            // Computes additional information about the analysis and creates the
            //  object-oriented structure to hold the discriminants found.
            this.CreateDiscriminants();

            this.Classifier = this.CreateClassifier();

            return this.Classifier;
        }


        /// <summary>
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
            public MultivariateKernelRegression First { get; set; }

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
}
