// Accord Formats Library
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

namespace Openize.Accord.IO.Csv
{
    using System;
    using System.Data;
    using System.Globalization;
    using System.IO;
    using System.Text;
#if !NETSTANDARD1_4
#endif

    /// <summary>
    ///   Writer for CSV data.
    /// </summary>
    /// 
    /// <example>
    /// <para>
    ///   The following example shows how to use <see cref="CsvWriter"/> to write a matrix in .csv format.</para>
    ///   <code source="Unit Tests\Accord.Tests.IO\CsvWriterTest.cs" region="doc_matrix" />
    ///   
    /// <para>
    ///   The following example shows how to use <see cref="CsvWriter"/> to write a jagged array in .csv format.</para>
    ///   <code source="Unit Tests\Accord.Tests.IO\CsvWriterTest.cs" region="doc_jagged" />
    ///   
    /// <para>
    ///   The following example shows how to use <see cref="CsvWriter"/> to write a DataTable in .csv format.</para>
    ///   <code source="Unit Tests\Accord.Tests.IO\CsvWriterTest.cs" region="doc_table" />
    ///   
    /// <para>
    ///   It is also possible to use <see cref="CsvWriter"/> to write matrices (or jagged arrays) 
    ///   containing objects with mixed types:</para>
    ///   <code source="Unit Tests\Accord.Tests.IO\CsvWriterTest.cs" region="doc_objects" />
    /// </example>
    /// 
    public class CsvWriter : IDisposable
    {
        /// <summary>
        ///   Gets the writer.
        /// </summary>
        /// 
        /// <value>
        ///   The writer.
        /// </value>
        /// 
        public TextWriter Writer { get; private set; }

        /// <summary>
        ///   Gets or sets the comment character indicating that a line is commented out.
        /// </summary>
        /// 
        /// <value>The comment character indicating that a line is commented out.</value>
        /// 
        public char Comment { get; set; }

        /// <summary>
        ///   Gets or sets the escape character letting insert quotation characters inside a quoted field.
        /// </summary>
        /// 
        /// <value>The escape character letting insert quotation characters inside a quoted field.</value>
        /// 
        public char Escape { get; set; }

        /// <summary>
        ///   Gets or sets the delimiter character separating each field.
        /// </summary>
        /// 
        /// <value>The delimiter character separating each field.</value>
        /// 
        public char Delimiter { get; set; }

        /// <summary>
        ///   Gets or sets the quotation character wrapping every field.
        /// </summary>
        /// 
        /// <value>The quotation character wrapping every field.</value>
        /// 
        public char Quote { get; set; }

        /// <summary>
        ///   Gets or sets the format provider to use when converting 
        ///   data-types to text representations. Default is to use
        ///   CultureInfo.InvariantCulture.
        /// </summary>
        /// 
        /// <value>
        ///   The format provider.
        /// </value>
        /// 
        public IFormatProvider FormatProvider { get; set; }

        /// <summary>
        ///   Initializes a new instance of the <see cref="CsvWriter"/> class.
        /// </summary>
        /// 
        /// <param name="path">The path to the file to be written.</param>
        /// <param name="delimiter">The field delimiter character to separate values in the CSV file.
        ///   If set to zero, will use the system's default text separator. Default is '\0' (zero).</param>
        /// 
        public CsvWriter(String path, char delimiter = CsvReader.DefaultDelimiter)
        {
            this.init(new StreamWriter(new FileStream(path, FileMode.Create, FileAccess.Write)), delimiter);
        }

        /// <summary>
        ///   Initializes a new instance of the <see cref="CsvWriter"/> class.
        /// </summary>
        /// 
        /// <param name="writer">A <see cref="T:TextWriter"/> pointing to the CSV file.</param>
        /// <param name="delimiter">The field delimiter character to separate values in the CSV file.
        ///   If set to zero, will use the system's default text separator. Default is '\0' (zero).</param>
        /// 
        public CsvWriter(TextWriter writer, char delimiter = CsvReader.DefaultDelimiter)
        {
            this.init(writer, delimiter);
        }

        private void init(TextWriter writer, char delimiter)
        {
            this.Writer = writer;
            this.Quote = CsvReader.DefaultQuote;

            this.Comment = CsvReader.DefaultComment;
            this.Escape = CsvReader.DefaultEscape;
            this.FormatProvider = System.Globalization.CultureInfo.InvariantCulture;
            this.Delimiter = delimiter;

            if (delimiter == '\0')
                this.Delimiter = CultureInfo.CurrentCulture.TextInfo.ListSeparator[0];
        }

        /// <summary>
        ///   Initializes a new instance of the <see cref="CsvWriter"/> 
        ///   class to write the CSV fields to a in-memory string.
        /// </summary>
        /// 
        /// <param name="builder">A <see cref="T:StringBuilder"/> to write to.</param>
        /// <param name="delimiter">The field delimiter character to separate values in the CSV file.
        ///   If set to zero, will use the system's default text separator. Default is '\0' (zero).</param>
        /// 
        public static CsvWriter ToText(StringBuilder builder, char delimiter = CsvReader.DefaultDelimiter)
        {
            using (var writer = new StringWriter(builder))
                return new CsvWriter(writer, delimiter);
        }

        /// <summary>
        ///   Writes the column names of a data table as the headers of the CSV file.
        /// </summary>
        /// 
        /// <param name="columnNames">A list of column names to use.</param>
        /// 
        public void WriteHeaders(params string[] columnNames)
        {
            var headers = new string[columnNames.Length];
            for (int i = 0; i < headers.Length; i++)
                headers[i] = this.quote(columnNames[i]);

            this.write(headers, String.Empty);
        }

#if !NETSTANDARD1_4
        /// <summary>
        ///   Writes the column names of a data table as the headers of the CSV file.
        /// </summary>
        /// 
        /// <param name="table">A DataTable whose columns names will be written as headers.</param>
        /// 
        public void WriteHeaders(DataTable table)
        {
            var headers = new string[table.Columns.Count];
            for (int i = 0; i < headers.Length; i++)
                headers[i] = this.quote(table.Columns[i].ColumnName);

            this.write(headers, String.Empty);
        }
#endif

        /// <summary>
        ///   Writes the specified matrix in CSV format.
        /// </summary>
        /// 
        /// <typeparam name="T">The matrix data type.</typeparam>
        /// <param name="table">The table to be written.</param>
        /// 
        public void Write<T>(T[,] table)
        {
            int rows = table.GetLength(0);
            int cols = table.GetLength(1);

            string[] items = new string[cols];

            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < items.Length; j++)
                    items[j] = this.quote(table[i, j]);

                this.write(items, String.Empty);
            }

            this.Writer.Flush();
        }

        /// <summary>
        ///   Writes the specified matrix in CSV format.
        /// </summary>
        /// 
        /// <typeparam name="T">The matrix data type.</typeparam>
        /// <param name="table">The table to be written.</param>
        /// 
        public void Write<T>(T[][] table)
        {
            for (int i = 0; i < table.Length; i++)
            {
                string[] items = new string[table[i].Length];

                for (int j = 0; j < items.Length; j++)
                    items[j] = this.quote(table[i][j]);

                this.write(items, String.Empty);
            }

            this.Writer.Flush();
        }

#if !NETSTANDARD1_4
        /// <summary>
        ///   Writes the specified table in a CSV format.
        /// </summary>
        /// 
        /// <param name="table">The data table to be written.</param>
        /// 
        public void Write(DataTable table)
        {
            this.WriteHeaders(table);

            string[] items = new string[table.Columns.Count];

            foreach (DataRow row in table.Rows)
            {
                for (int i = 0; i < items.Length; i++)
                    items[i] = this.quote(row[i]);

                this.write(items, row.RowError);
            }

            this.Writer.Flush();
        }
#endif

        /// <summary>
        ///   Writes the specified fields in a CSV format.
        /// </summary>
        /// 
        /// <param name="fields">The fields to be written.</param>
        /// 
        public void WriteLine<T>(T[] fields)
        {
            this.WriteLine(fields, String.Empty);
        }

        /// <summary>
        ///   Writes the specified fields in a CSV format.
        /// </summary>
        /// 
        /// <param name="fields">The fields to be written.</param>
        /// <param name="comment">An optional comment for the line.</param>
        /// 
        public void WriteLine<T>(T[] fields, string comment)
        {
            string[] items = new string[fields.Length];

            for (int i = 0; i < items.Length; i++)
                items[i] = this.quote(fields[i]);

            this.write(items, comment);
        }





        private void write(string[] fields, string comment)
        {
            this.Writer.Write(String.Join(this.Delimiter.ToString(), fields));

            if (!String.IsNullOrEmpty(comment))
                this.Writer.Write(" {0} {1}", this.Comment, comment);

            this.Writer.WriteLine();
        }

        private string escape(object obj)
        {
            string text = String.Format(this.FormatProvider, "{0}", obj);

            text = text.Replace(this.Quote.ToString(), new String(new[] { this.Escape, this.Quote }));

            return text;
        }

        private string quote(object obj)
        {
            string text = this.escape(obj);
            return String.Format(this.FormatProvider, "{0}{1}{0}", this.Quote, text);
        }





        /// <summary>
        ///   Performs application-defined tasks associated with 
        ///   freeing,  releasing, or resetting unmanaged resources.
        /// </summary>
        /// 
        public void Dispose()
        {
            this.Dispose(true);
            GC.SuppressFinalize(this);
        }

        /// <summary>
        ///   Finalizes an instance of the <see cref="CsvWriter"/> class.
        /// </summary>
        /// 
        ~CsvWriter()
        {
            this.Dispose(false);
        }

        /// <summary>
        ///   Releases unmanaged and - optionally - managed resources.
        /// </summary>
        /// 
        /// <param name="disposing"><c>true</c> to release both managed and
        ///   unmanaged resources; <c>false</c> to release only unmanaged resources.
        /// </param>
        /// 
        protected virtual void Dispose(bool disposing)
        {
            if (disposing)
            {
                // free managed resources
                if (this.Writer != null)
                {
                    this.Writer.Dispose();
                    this.Writer = null;
                }
            }
        }

    }
}
