﻿// Accord Core Library
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

namespace Openize.Accord.Core
{
    using System;
    using System.Diagnostics;

    /// <summary>
    ///   Temporary internal framework class for handling debug assertions
    ///   while inside unit tests. Will be removed in a future release.
    /// </summary>
    /// 
    public static class Debug
    {
        /// <summary>
        ///   Throws an exception if a condition is false.
        /// </summary>
        /// 
        [Conditional("DEBUG")]
        public static void Assert(bool condition, string message = "Internal framework error.")
        {
            if (!condition)
                throw new Exception(message);
        }
    }
}
