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

namespace Openize.Accord.Statistics.Testing.Power.Base
{
    /// <summary>
    ///   Common interface for two-sample power analysis objects.
    /// </summary>
    /// 
    public interface ITwoSamplePowerAnalysis : IPowerAnalysis
    {

        /// <summary>
        ///   Gets the number of observations 
        ///   contained in the first sample.
        /// </summary>
        /// 
        double Samples1 { get; }

        /// <summary>
        ///   Gets the number of observations 
        ///   contained in the second sample.
        /// </summary>
        /// 
        double Samples2 { get; }

    }
}
