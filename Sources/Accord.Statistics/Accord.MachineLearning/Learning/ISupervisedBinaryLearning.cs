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

namespace Openize.Accord.Statistics.Accord.MachineLearning.Learning
{
    using Openize.Accord.Core.MachineLearning.Classifiers.Multiclass.Binary;

    /// <summary>
    ///   Common interface for supervised learning algorithms for
    ///   <see cref="IBinaryClassifier{TInput}">binary classifiers</see>.
    /// </summary>
    /// 
    /// <typeparam name="TModel">The type for the model being learned.</typeparam>
    /// <typeparam name="TInput">The type for the input data that enters the model.</typeparam>
    /// 
    public interface ISupervisedBinaryLearning<out TModel, in TInput> :
        ISupervisedMulticlassLearning<TModel, TInput>,
        ISupervisedLearning<TModel, TInput, bool>
        where TModel : IBinaryClassifier<TInput>
    {
    }

    /// <summary>
    ///   Common interface for supervised learning algorithms for
    ///   <see cref="IBinaryClassifier{TInput}">binary classifiers</see>.
    /// </summary>
    /// 
    /// <typeparam name="TModel">The type for the model being learned.</typeparam>
    /// 
    public interface ISupervisedBinaryLearning<out TModel> :
        ISupervisedBinaryLearning<TModel, int[]>,
        ISupervisedBinaryLearning<TModel, float[]>,
        ISupervisedBinaryLearning<TModel, double[]>
        where TModel : IBinaryClassifier
    {
    }

    /*
    public interface IWeightedSupervisedBinaryLearning<TModel, TInput> :
        IWeightedSupervisedMulticlassLearning<TModel, TInput>,
        IWeightedSupervisedLearning<TModel, TInput, bool>
        where TModel : IBinaryClassifier<TInput>
    {
    }

    public interface IWeightedSupervisedBinaryLearning<TModel> :
        IWeightedSupervisedBinaryLearning<TModel, int[]>,
        IWeightedSupervisedBinaryLearning<TModel, float[]>,
        IWeightedSupervisedBinaryLearning<TModel, double[]>
        where TModel : IBinaryClassifier
    {
    }*/
}
