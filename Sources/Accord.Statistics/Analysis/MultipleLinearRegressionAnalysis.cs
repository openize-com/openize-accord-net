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
    using System.Collections.Generic;
    using System.Collections.ObjectModel;
    using System.ComponentModel;
    using System.Threading;
    using Accord.MachineLearning.Classifiers;
    using Accord.MachineLearning.Learning;
    using Base;
    using Openize.Accord.Math.Matrix;
    using Models.Regression.Linear;
    using Models.Regression.Linear.Fitting;
    using Openize.Accord.Core.Ranges;
    using Openize.Accord.Math.Accord.Statistics;
    using Testing;
    using Testing.Multiple_Samples;

    /// <summary>
    ///   Multiple Linear Regression Analysis
    /// </summary>
    /// 
    /// <remarks>
    /// <para>
    ///   Linear regression is an approach to model the relationship between 
    ///   a single scalar dependent variable <c>y</c> and one or more explanatory
    ///   variables <c>x</c>. This class uses a <see cref="MultipleLinearRegression"/>
    ///   to extract information about a given problem, such as confidence intervals,
    ///   hypothesis tests and performance measures.</para>
    ///   
    /// <para>
    ///   This class can also be bound to standard controls such as the 
    ///   <a href="http://msdn.microsoft.com/en-us/library/system.windows.forms.datagridview.aspx">DataGridView</a>
    ///   by setting their DataSource property to the analysis' <see cref="Coefficients"/> property.</para>
    ///   
    /// <para>
    ///   References:
    ///   <list type="bullet">
    ///     <item><description>
    ///       Wikipedia contributors. "Linear regression." Wikipedia, The Free Encyclopedia, 4 Nov. 2012.
    ///       Available at: http://en.wikipedia.org/wiki/Linear_regression </description></item>
    ///  </list></para>  
    /// </remarks>
    /// 
    /// <example>
    /// <para>
    ///   The first example shows how to learn a multiple linear regression analysis 
    ///   from a dataset already given in matricial form (using jagged double[][] arrays).</para>
    /// <code source="Unit Tests\Accord.Tests.Statistics\Analysis\MultipleLinearRegressionAnalysisTest.cs" region="doc_learn_part2" />
    /// 
    /// <code>
    /// // Now we can show a summary of analysis
    /// // Accord.Controls.DataGridBox.Show(regression.Coefficients);
    /// </code>
    /// 
    ///   <img src="..\images\linear-regression.png" />
    /// 
    /// <code>
    /// // We can also show a summary ANOVA
    /// DataGridBox.Show(regression.Table);
    /// </code>
    /// 
    ///   <img src="..\images\linear-anova.png" />
    ///   
    /// <code source="Unit Tests\Accord.Tests.Statistics\Analysis\MultipleLinearRegressionAnalysisTest.cs" region="doc_learn_part2" />
    /// 
    /// <para>
    ///   The second example shows how to learn a multiple linear regression analysis using data 
    ///   given in the form of a System.Data.DataTable. This data is also heterogeneous, mixing 
    ///   both discrete (symbol) variables and continuous variables. This example is also available
    ///   for <see cref="LogisticRegressionAnalysis"/>.</para>
    ///   
    /// <code source="Unit Tests\Accord.Tests.Statistics\Analysis\MultipleLinearRegressionAnalysisTest.cs" region="doc_learn_database" />
    /// </example>
    /// 
    /// <seealso cref="LogisticRegressionAnalysis"/>
    /// 
    [Serializable]
    public class MultipleLinearRegressionAnalysis : TransformBase<double[], double>,
        IRegressionAnalysis, IAnova,
        ISupervisedLearning<MultipleLinearRegression, double[], double>
    {
        [NonSerialized]
        CancellationToken token = new CancellationToken();

        /// <summary>
        /// Gets or sets a cancellation token that can be used to
        /// stop the learning algorithm while it is running.
        /// </summary>
        public CancellationToken Token
        {
            get { return this.token; }
            set { this.token = value; }
        }


        internal MultipleLinearRegression regression;

        private string[] inputNames;
        private string outputName;

        [Obsolete]
        private double[][] inputData;
        [Obsolete]
        private double[] outputData;
        [Obsolete]
        private double[] results;


        private double[][] informationMatrix;

        private double SSe; // Error sum of squares
        private double SSr; // Regression sum of squares
        private double SSt; // Total sum of squares

        private double MSe; // Error mean sum of squares
        private double MSr; // Regression sum of squares
        private double MSt; // Total sum of squares

        private int DFe; // Error degrees of freedom
        private int DFr; // Regression degrees of freedom
        private int DFt; // Total degrees of freedom

        // Result related measures
        private double outputMean;
        private double stdError;
        private double rSquared;
        private double rAdjusted;
        private FTest ftest;
        private TTest ttest;
        private ZTest ztest;
        private ChiSquareTest chiSquareTest;

        // Coefficient measures
        internal double[] standardErrors;
        internal DoubleRange[] confidences;
        internal double confidencePercent = 0.95;
        internal TTest[] ttests;
        internal FTest[] ftests;

        private AnovaSourceCollection anovaTable;
        private LinearRegressionCoefficientCollection coefficientCollection;


        /// <summary>
        ///   Source data used in the analysis.
        /// </summary>
        /// 
        [Obsolete("This property will be removed.")]
        public double[,] Source { get; private set; }

        /// <summary>
        ///   Source data used in the analysis.
        /// </summary>
        /// 
        [Obsolete("This property will be removed.")]
        public double[][] Array { get { return this.inputData; } }

        /// <summary>
        ///   Gets the dependent variable value
        ///   for each of the source input points.
        /// </summary>
        /// 
        [Obsolete("This property will be removed.")]
        public double[] Outputs
        {
            get { return this.outputData; }
        }

        /// <summary>
        ///   Gets the resulting values obtained
        ///   by the linear regression model.
        /// </summary>
        /// 
        [Obsolete("This property will be removed.")]
        public double[] Results
        {
            get { return this.results; }
        }

        /// <summary>
        ///   Gets or sets the learning algorithm used to learn the <see cref="MultipleLinearRegression"/>.
        /// </summary>
        /// 
        public OrdinaryLeastSquares OrdinaryLeastSquares { get; set; }

        /// <summary>
        ///   Gets the standard deviation of the errors. 
        /// </summary>
        /// 
        public double StandardError
        {
            get { return this.stdError; }
        }

        /// <summary>
        /// Gets the information matrix obtained during learning.
        /// </summary>
        /// 
        public double[][] InformationMatrix
        {
            get { return this.informationMatrix; }
        }

        /// <summary>
        ///   Gets the coefficient of determination, as known as R²
        /// </summary>
        /// 
        public double RSquared
        {
            get { return this.rSquared; }
        }

        /// <summary>
        /// Gets the number of samples used to compute the analysis.
        /// </summary>
        /// 
        public int NumberOfSamples { get; private set; }

        /// <summary>
        ///   Gets the adjusted coefficient of determination, as known as R² adjusted
        /// </summary>
        /// 
        public double RSquareAdjusted
        {
            get { return this.rAdjusted; }
        }

        /// <summary>
        ///   Gets a F-Test between the expected outputs and results.
        /// </summary>
        /// 
        public FTest FTest
        {
            get { return this.ftest; }
        }

        /// <summary>
        ///   Gets a Z-Test between the expected outputs and the results.
        /// </summary>
        /// 
        public ZTest ZTest
        {
            get { return this.ztest; }
        }

        /// <summary>
        ///   Gets a Chi-Square Test between the expected outputs and the results.
        /// </summary>
        /// 
        public ChiSquareTest ChiSquareTest
        {
            get { return this.chiSquareTest; }
        }

        /// <summary>
        ///   Gets the Standard Error for each coefficient
        ///   found during the logistic regression.
        /// </summary>
        /// 
        public double[] StandardErrors
        {
            get { return this.standardErrors; }
        }

        /// <summary>
        ///   Gets the Regression model created
        ///   and evaluated by this analysis.
        /// </summary>
        /// 
        public MultipleLinearRegression Regression
        {
            get { return this.regression; }
        }

        /// <summary>
        ///   Gets the value of each coefficient.
        /// </summary>
        /// 
        public double[] CoefficientValues
        {
            get { return this.regression.Weights.Concatenate(this.regression.Intercept); }
        }

        /// <summary>
        ///   Gets or sets the name of the input variables for the model.
        /// </summary>
        /// 
        public String[] Inputs
        {
            get { return this.inputNames; }
            set { this.inputNames = value; }
        }

        /// <summary>
        ///   Gets or sets the name of the output variable for the model.
        /// </summary>
        /// 
        public String Output
        {
            get { return this.outputName; }
            set { this.outputName = value; }
        }

        /// <summary>
        ///   Gets the Confidence Intervals (C.I.)
        ///   for each coefficient found in the regression.
        /// </summary>
        /// 
        public DoubleRange[] Confidences
        {
            get { return this.confidences; }
        }

        /// <summary>
        ///   Gets the ANOVA table for the analysis.
        /// </summary>
        /// 
        public AnovaSourceCollection Table { get { return this.anovaTable; } }

        /// <summary>
        ///   Gets the collection of coefficients of the model.
        /// </summary>
        /// 
        public LinearRegressionCoefficientCollection Coefficients { get { return this.coefficientCollection; } }

        /// <summary>
        ///   Constructs a Multiple Linear Regression Analysis.
        /// </summary>
        /// 
        /// <param name="intercept">True to use an intercept term, false otherwise. Default is false.</param>
        /// 
        public MultipleLinearRegressionAnalysis(bool intercept = false)
        {
            this.OrdinaryLeastSquares = new OrdinaryLeastSquares() { UseIntercept = intercept };
        }

        /// <summary>
        ///   Constructs a Multiple Linear Regression Analysis.
        /// </summary>
        /// 
        /// <param name="inputs">The input data for the analysis.</param>
        /// <param name="outputs">The output data for the analysis.</param>
        /// <param name="intercept">True to use an intercept term, false otherwise. Default is false.</param>
        /// 
        [Obsolete("Please pass the 'inputs' and 'outputs' parameters to the Learn method instead.")]
        public MultipleLinearRegressionAnalysis(double[][] inputs, double[] outputs, bool intercept = false)
            : this(intercept)
        {
            // Initial argument checking
            this.init(inputs, outputs);

            // Store data sets
            this.Source = inputs.ToMatrix();
            this.inputData = inputs;
            this.outputData = outputs;
        }

        private void init(double[][] inputs, double[] outputs)
        {
            if (inputs == null)
                throw new ArgumentNullException("inputs");

            if (outputs == null)
                throw new ArgumentNullException("outputs");

            if (inputs.Length != outputs.Length)
                throw new ArgumentException("The number of rows in the input array must match the number of given outputs.");

            this.NumberOfInputs = inputs[0].Length;
            this.NumberOfOutputs = 1;

            for (int i = 0; i < inputs.Length; i++)
                if (inputs[i].Length != this.NumberOfInputs)
                    throw new ArgumentException("All input vectors must have the same length.");

            // Create the linear regression
            this.regression = new MultipleLinearRegression()
            {
                NumberOfInputs = this.NumberOfInputs
            };

            // Create additional structures
            int coefficientCount = this.NumberOfInputs + 1;
            this.standardErrors = new double[coefficientCount];
            this.confidences = new DoubleRange[coefficientCount];
            this.ftests = new FTest[coefficientCount];
            this.ttests = new TTest[coefficientCount];

            if (this.outputName == null)
                this.outputName = "Output";

            if (this.inputNames == null)
            {
                this.inputNames = new string[this.NumberOfInputs];
                for (int i = 0; i < this.inputNames.Length; i++)
                    this.inputNames[i] = "Input " + i;
            }

            // Create object-oriented structure to represent the analysis
            var coefs = new LinearRegressionCoefficient[coefficientCount];
            for (int i = 0; i < coefs.Length; i++)
                coefs[i] = new LinearRegressionCoefficient(this, i);
            this.coefficientCollection = new LinearRegressionCoefficientCollection(this, coefs);
        }

        /// <summary>
        ///   Constructs a Multiple Linear Regression Analysis.
        /// </summary>
        /// 
        /// <param name="inputs">The input data for the analysis.</param>
        /// <param name="outputs">The output data for the analysis.</param>
        /// <param name="intercept">True to use an intercept term, false otherwise. Default is false.</param>
        /// <param name="inputNames">The names of the input variables.</param>
        /// <param name="outputName">The name of the output variable.</param>
        /// 
        [Obsolete("Please pass the 'inputs' and 'outputs' parameters to the Learn method instead.")]
        public MultipleLinearRegressionAnalysis(double[][] inputs, double[] outputs,
            String[] inputNames, String outputName, bool intercept = false)
            : this(inputs, outputs, intercept)
        {
            if (inputNames.Length != this.inputNames.Length)
            {
                throw new ArgumentException("The input names vector should have the same length"
                  + " as the number of variables in the analysis. In this analysis, there are "
                  + this.inputNames.Length + " variables expected.");
            }

            this.inputNames = inputNames;
            this.outputName = outputName;
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
        public MultipleLinearRegression Learn(double[][] x, double[] y, double[] weights = null)
        {
            if (weights != null)
                throw new ArgumentException(global::Accord.Properties.Resources.NotSupportedWeights, "weights");

            this.init(x, y);
            this.compute(x, y);
            return this.regression;
        }

        /// <summary>
        ///   Computes the Multiple Linear Regression Analysis.
        /// </summary>
        /// 
        [Obsolete("Please use the Learn method instead.")]
        public void Compute()
        {
            this.compute(this.Source.ToJagged(), this.Outputs);
        }

        private void compute(double[][] x, double[] y)
        {
            int n = x.Length;
            int p = this.NumberOfInputs;

            this.SSt = 0;
            this.SSe = 0;
            this.outputMean = 0.0;
            this.NumberOfSamples = x.Length;

            // Compute the regression
            this.OrdinaryLeastSquares.Token = this.Token;
            this.regression = this.OrdinaryLeastSquares.Learn(x, y);
            this.informationMatrix = this.OrdinaryLeastSquares.GetInformationMatrix();

            // Calculate mean of the expected outputs
            this.outputMean = y.Mean();

            // Calculate actual outputs (results)
#pragma warning disable 612, 618
            this.results = this.regression.Transform(x);

            // Calculate SSe and SSt
            for (int i = 0; i < x.Length; i++)
            {
                double d;
                d = y[i] - this.results[i];
                this.SSe += d * d;

                d = y[i] - this.outputMean;
                this.SSt += d * d;
            }


            // Calculate SSr
            this.SSr = this.SSt - this.SSe;

            // Calculate R-Squared
            this.rSquared = (this.SSt != 0) ? 1.0 - (this.SSe / this.SSt) : 1.0;

            // Calculated Adjusted R-Squared
            if (this.rSquared == 1)
            {
                this.rAdjusted = 1;
            }
            else
            {
                if (n - p == 1)
                {
                    this.rAdjusted = double.NaN;
                }
                else
                {
                    this.rAdjusted = 1.0 - (1.0 - this.rSquared) * ((n - 1.0) / (n - p - 1.0));
                }
            }

            // Calculate Degrees of Freedom
            this.DFr = p;
            this.DFe = n - (p + 1);
            this.DFt = this.DFr + this.DFe;

            // Calculate Sum of Squares Mean
            this.MSe = this.SSe / this.DFe;
            this.MSr = this.SSr / this.DFr;
            this.MSt = this.SSt / this.DFt;

            // Calculate the F statistic
            this.ftest = new FTest(this.MSr / this.MSe, this.DFr, this.DFe);
            this.stdError = Math.Sqrt(this.MSe);

            // Create the ANOVA table
            List<AnovaVariationSource> table = new List<AnovaVariationSource>();
            table.Add(new AnovaVariationSource(this, "Regression", this.SSr, this.DFr, this.MSr, this.ftest));
            table.Add(new AnovaVariationSource(this, "Error", this.SSe, this.DFe, this.MSe, null));
            table.Add(new AnovaVariationSource(this, "Total", this.SSt, this.DFt, this.MSt, null));
            this.anovaTable = new AnovaSourceCollection(table);

            // Compute coefficient standard errors;
            this.standardErrors = new double[this.NumberOfInputs + 1];
            for (int i = 0; i < this.informationMatrix.Length; i++)
                this.standardErrors[i] = Math.Sqrt(this.MSe * this.informationMatrix[i][i]);

            // Compute coefficient tests
            for (int i = 0; i < this.CoefficientValues.Length; i++)
            {
                double tStatistic = this.CoefficientValues[i] / this.standardErrors[i];

                this.ttests[i] = new TTest(estimatedValue: this.CoefficientValues[i],
                    standardError: this.standardErrors[i], degreesOfFreedom: this.DFe);

                this.ftests[i] = new FTest(tStatistic * tStatistic, 1, this.DFe);

                this.confidences[i] = this.ttests[i].GetConfidenceInterval(this.confidencePercent);
            }


            // Compute model performance tests
            this.ttest = new TTest(this.results, this.outputMean);
            this.ztest = new ZTest(this.results, this.outputMean);
            this.chiSquareTest = new ChiSquareTest(y, this.results, n - p - 1);
#pragma warning restore 612, 618
        }

        internal void setConfidenceIntervals(double percent)
        {
            this.confidencePercent = percent;
            for (int i = 0; i < this.ttests.Length; i++)
                this.confidences[i] = this.ttest.GetConfidenceInterval(percent);
        }

        /// <summary>
        /// Applies the transformation to an input, producing an associated output.
        /// </summary>
        /// <param name="input">The input data to which the transformation should be applied.</param>
        /// <returns>
        /// The output generated by applying this transformation to the given input.
        /// </returns>
        public override double Transform(double[] input)
        {
            return this.regression.Transform(input);
        }

        /// <summary>
        /// Gets the prediction interval for a given input.
        /// </summary>
        /// 
        public DoubleRange GetPredictionInterval(double[] input)
        {
            return this.regression.GetPredictionInterval(input, Math.Sqrt(this.MSe), this.NumberOfSamples, this.InformationMatrix, this.confidencePercent);
        }

        /// <summary>
        /// Gets the confidence interval for a given input.
        /// </summary>
        /// 
        public DoubleRange GetConfidenceInterval(double[] input)
        {
            return this.regression.GetConfidenceInterval(input, Math.Sqrt(this.MSe), this.NumberOfSamples, this.InformationMatrix, this.confidencePercent);
        }
    }

    /// <summary>
    /// <para>
    ///   Represents a Linear Regression coefficient found in the Multiple
    ///   Linear Regression Analysis allowing it to be bound to controls like
    ///   the DataGridView. </para>
    /// <para>
    ///   This class cannot be instantiated.</para>   
    /// </summary>
    /// 
    [Serializable]
    public class LinearRegressionCoefficient
    {

        private int index;
        private MultipleLinearRegressionAnalysis analysis;


        /// <summary>
        ///   Creates a regression coefficient representation.
        /// </summary>
        /// 
        /// <param name="analysis">The analysis to which this coefficient belongs.</param>
        /// <param name="index">The coefficient index.</param>
        /// 
        internal LinearRegressionCoefficient(MultipleLinearRegressionAnalysis analysis, int index)
        {
            this.index = index;
            this.analysis = analysis;
        }


        /// <summary>
        ///   Gets the Index of this coefficient on the original analysis coefficient collection.
        /// </summary>
        /// 
        public int Index
        {
            get { return this.index; }
        }

        /// <summary>
        ///   Returns a reference to the parent analysis object.
        /// </summary>
        /// 
        [Browsable(false)]
        public MultipleLinearRegressionAnalysis Analysis
        {
            get { return this.analysis; }
        }

        /// <summary>
        ///   Gets the name for the current coefficient.
        /// </summary>
        /// 
        public string Name
        {
            get
            {
                if (this.IsIntercept) return "Intercept";
                else return this.analysis.Inputs[this.index];
            }
        }

        /// <summary>
        ///   Gets a value indicating whether this coefficient is an intercept term.
        /// </summary>
        /// 
        /// <value>
        ///   <c>true</c> if this coefficient is the intercept; otherwise, <c>false</c>.
        /// </value>
        /// 
        [DisplayName("Intercept?")]
        public bool IsIntercept { get { return this.index == this.Analysis.NumberOfInputs; } }

        /// <summary>
        ///   Gets the coefficient value.
        /// </summary>
        /// 
        [DisplayName("Value")]
        public double Value { get { return this.Analysis.CoefficientValues[this.index]; } }

        /// <summary>
        ///   Gets the Standard Error for the current coefficient.
        /// </summary>
        /// 
        [DisplayName("Std. Error")]
        public double StandardError { get { return this.Analysis.StandardErrors[this.index]; } }

        /// <summary>
        ///   Gets the T-test performed for this coefficient.
        /// </summary>
        /// 
        public TTest TTest { get { return this.Analysis.ttests[this.index]; } }

        /// <summary>
        ///   Gets the F-test performed for this coefficient.
        /// </summary>
        /// 
        public FTest FTest { get { return this.Analysis.ftests[this.index]; } }

        /// <summary>
        ///   Gets the confidence interval (C.I.) for the current coefficient.
        /// </summary>
        /// 
        [Browsable(false)]
        public DoubleRange Confidence
        {
            get { return this.analysis.Confidences[this.index]; }
        }

        /// <summary>
        ///   Gets the upper limit for the confidence interval.
        /// </summary>
        /// 
        [DisplayName("Upper confidence limit")]
        public double ConfidenceUpper { get { return this.Analysis.Confidences[this.index].Max; } }

        /// <summary>
        ///   Gets the lower limit for the confidence interval.
        /// </summary>
        /// 
        [DisplayName("Lower confidence limit")]
        public double ConfidenceLower { get { return this.Analysis.Confidences[this.index].Min; } }

    }

    /// <summary>
    ///   Represents a Collection of Linear Regression Coefficients found in the 
    ///   <see cref="MultipleLinearRegressionAnalysis"/>. This class cannot be instantiated.
    /// </summary>
    /// 
    [Serializable]
    public class LinearRegressionCoefficientCollection : ReadOnlyCollection<LinearRegressionCoefficient>
    {

        MultipleLinearRegressionAnalysis analysis;

        /// <summary>
        ///   Gets or sets the size of the confidence
        ///   intervals reported for the coefficients.
        ///   Default is 0.95.
        /// </summary>
        /// 
        public double ConfidencePercent
        {
            get { return this.analysis.confidencePercent; }
            set { this.analysis.setConfidenceIntervals(value); }
        }

        internal LinearRegressionCoefficientCollection(MultipleLinearRegressionAnalysis analysis,
            LinearRegressionCoefficient[] components)
            : base(components)
        {
            this.analysis = analysis;
        }
    }

}
