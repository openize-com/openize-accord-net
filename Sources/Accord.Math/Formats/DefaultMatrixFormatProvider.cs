// Accord Math Library
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

namespace FileFormat.Accord.Math.Formats
{
    using System;
    using System.Globalization;
    using Base;
    using global::Accord.Math;

    /// <summary>
    ///   Gets the default matrix representation, where each row
    ///   is separated by a new line, and columns are separated by spaces.
    /// </summary>
    /// 
    /// <remarks>
    ///   This class can be used to convert to and from C#
    ///   matrices and their string representation. Please 
    ///   see the example for details.
    /// </remarks>
    /// 
    /// <example>
    /// <para>
    ///   Converting from a multidimensional matrix to a 
    ///   string representation:</para>
    ///   
    /// <code>
    ///   // Declare a number array
    ///   double[,] x = 
    ///   {
    ///      { 1, 2, 3, 4 },
    ///      { 5, 6, 7, 8 },
    ///   };
    ///   
    ///   // Convert the aforementioned array to a string representation:
    ///   string str = x.ToString(DefaultMatrixFormatProvider.CurrentCulture);
    ///   
    ///   // the final result will be equivalent to
    ///   @"1, 2, 3, 4
    ///     5, 6, 7, 8";
    /// </code>
    /// 
    /// <para>
    ///   Converting from strings to actual matrices:</para>
    /// 
    /// <code>
    ///   // Declare an input string
    ///   string str = @"1, 2, 3, 4
    ///                 "5, 6, 7, 8";
    ///   
    ///   // Convert the string representation to an actual number array:
    ///   double[,] matrix = Matrix.Parse(str, DefaultMatrixFormatProvider.InvariantCulture);
    ///   
    ///   // matrix will now contain the actual multidimensional 
    ///   // matrix representation of the given string.
    /// </code>
    /// </example>
    /// 
    /// <seealso cref="Matrix"/>
    /// <seealso cref="CSharpMatrixFormatProvider"/>
    /// 
    /// <seealso cref="CSharpJaggedMatrixFormatProvider"/>
    /// <seealso cref="CSharpArrayFormatProvider"/>
    /// 
    /// <seealso cref="OctaveMatrixFormatProvider"/>
    /// <seealso cref="OctaveArrayFormatProvider"/>
    /// 
    public sealed class DefaultMatrixFormatProvider : MatrixFormatProviderBase
    {

        /// <summary>
        ///   Initializes a new instance of the <see cref="DefaultMatrixFormatProvider"/> class.
        /// </summary>
        /// 
        public DefaultMatrixFormatProvider(IFormatProvider innerProvider)
            : base(innerProvider)
        {
            this.FormatMatrixStart = String.Empty;
            this.FormatMatrixEnd = String.Empty;
            this.FormatRowStart = String.Empty;
            this.FormatRowEnd = String.Empty;
            this.FormatColStart = String.Empty;
            this.FormatColEnd = String.Empty;
            this.FormatRowDelimiter = " \n";
            this.FormatColDelimiter = " ";

            this.ParseMatrixStart = String.Empty;
            this.ParseMatrixEnd = String.Empty;
            this.ParseRowStart = String.Empty;
            this.ParseRowEnd = String.Empty;
            this.ParseColStart = String.Empty;
            this.ParseColEnd = String.Empty;
            this.ParseRowDelimiter = "\n";
            this.ParseColDelimiter = " ";
        }

        /// <summary>
        ///   Gets the IMatrixFormatProvider which uses the CultureInfo used by the current thread.
        /// </summary>
        /// 
        public static DefaultMatrixFormatProvider CurrentCulture
        {
            get { return currentCulture; }
        }
        
        /// <summary>
        ///   Gets the IMatrixFormatProvider which uses the invariant system culture.
        /// </summary>
        /// 
        public static DefaultMatrixFormatProvider InvariantCulture
        {
            get { return invariantCulture; }
        }


        private static readonly DefaultMatrixFormatProvider currentCulture =
            new DefaultMatrixFormatProvider(CultureInfo.CurrentCulture);

        private static readonly DefaultMatrixFormatProvider invariantCulture = 
            new DefaultMatrixFormatProvider(CultureInfo.InvariantCulture);

    }
}
