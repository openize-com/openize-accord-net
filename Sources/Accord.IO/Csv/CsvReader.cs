// Accord Formats Library
// The Accord.NET Framework
// http://accord-framework.net
//
// LumenWorks.Framework.IO.CSV.CsvReader
// Copyright (c) 2005 Sébastien Lorion
//
// Copyright © César Souza, 2009-2017
// cesarsouza at gmail.com
//
// This class has been based on the original work by Sébastien Lorion, originally
// published under the MIT license (and thus compatible with the LGPL). Original
// license text is reproduced below:
//
//    MIT license (http://en.wikipedia.org/wiki/MIT_License)
//
//    Permission is hereby granted, free of charge, to any person obtaining a copy
//    of this software and associated documentation files (the "Software"), to deal
//    in the Software without restriction, including without limitation the rights 
//    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
//    of the Software, and to permit persons to whom the Software is furnished to do so, 
//    subject to the following conditions:
//
//    The above copyright notice and this permission notice shall be included in all 
//    copies or substantial portions of the Software.
//
//    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
//    INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
//    PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
//    FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
//    ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//

namespace Openize.Accord.IO.Csv
{
#if !NETSTANDARD1_4
    using System.Data;
#endif
    using System;
    using System.Collections;
    using System.Collections.Generic;
    using System.Globalization;
    using System.IO;
    using System.Net;
    using Events;
    using Exceptions;
    using global::Accord.IO.Resources;
    using Openize.Accord.Core;
    using Debug = System.Diagnostics.Debug;

    /// <summary>
    ///   Represents a reader that provides fast, non-cached, forward-only access to CSV data.  
    /// </summary>
    /// 
    public partial class CsvReader :
#if !NETSTANDARD1_4
        IDataReader,
#endif
        IEnumerable<string[]>, IDisposable
    {
        /// <summary>
        ///   Defines the default buffer size.
        /// </summary>
        /// 
        public const int DefaultBufferSize = 0x1000;

        /// <summary>
        ///   Defines the default delimiter character separating each field.
        /// </summary>
        /// 
        public const char DefaultDelimiter = '\0';

        /// <summary>
        ///   Defines the default quote character wrapping every field.
        /// </summary>
        /// 
        public const char DefaultQuote = '"';

        /// <summary>
        ///   Defines the default escape character letting insert quotation characters inside a quoted field.
        /// </summary>
        /// 
        public const char DefaultEscape = '"';

        /// <summary>
        ///   Defines the default comment character indicating that a line is commented out.
        /// </summary>
        /// 
        public const char DefaultComment = '#';


        private char[] delimiters = { ',', ';', '\t', '|', ';', '^', ' ' };

        private static readonly StringComparer _fieldHeaderComparer =
            StringComparer.CurrentCultureIgnoreCase;


        private TextReader _reader; // pointing to the CSV file
        private int _bufferSize;

        private bool _initialized; // if the class is initialized.
        private string[] _fieldHeaders; // the field headers.

        //  Dictionary of field indexes by header. The key 
        //  is the field name and the value is its index.
        private Dictionary<string, int> _fieldHeaderIndexes;

        private long _currentRecordIndex; // the current record index in the CSV file.
        // A value of int.MinValue means that the reader has not been initialized yet.
        // Otherwise, a negative value means that no record has been read yet.

        private int _nextFieldStart; // starting position of the next unread field.
        private int _nextFieldIndex; // the index of the next unread field.

        private string[] _fields; // array of the field values for the current record. A
        // Null value indicates the field has not been parsed.

        private int _fieldCount; //  number of fields to retrieve for each record
        private char[] _buffer; // the read buffer

        private int _bufferLength; // current read buffer length
        private bool _eof; // if the end of the reader has been reached
        private bool _eol; // if the last read operation reached an EOL character

        private bool _firstRecordInCache; // if the first record is in cache.
        // This can happen when initializing a reader with no headers
        // because one record must be read to get the field count automatically

        private bool _missingFieldFlag; // if one or more field are missing 
        // for the current record. Resets after each successful record read.

        private bool _parseErrorFlag; // if a parse error occurred for the 
                                      // current record. Resets after each successful record read.


        /// <summary>
        ///   Initializes a new instance of the CsvReader class.
        /// </summary>
        /// 
        /// <param name="path">The path for the CSV file.</param>
        /// <param name="hasHeaders"><see langword="true"/> if field names are located on the first non commented line, otherwise, <see langword="false"/>.</param>
        /// 
        public CsvReader(string path, bool hasHeaders)
        {
            this.init(new StreamReader(new FileStream(path, FileMode.Open, FileAccess.Read)), hasHeaders);
        }

        /// <summary>
        ///   Initializes a new instance of the CsvReader class.
        /// </summary>
        /// 
        /// <param name="stream">A <see cref="T:Stream"/> pointing to the CSV file.</param>
        /// <param name="hasHeaders"><see langword="true"/> if field names are located on the first non commented line, otherwise, <see langword="false"/>.</param>
        /// 
        public CsvReader(Stream stream, bool hasHeaders)
        {
            this.init(new StreamReader(stream), hasHeaders);
        }

        /// <summary>
        ///   Initializes a new instance of the CsvReader class.
        /// </summary>
        /// 
        /// <param name="reader">A <see cref="T:TextReader"/> pointing to the CSV file.</param>
        /// <param name="hasHeaders"><see langword="true"/> if field names are located on the first non commented line, otherwise, <see langword="false"/>.</param>
        /// 
        public CsvReader(TextReader reader, bool hasHeaders)
        {
            this.init(reader, hasHeaders);
        }

        /// <summary>
        ///   Initializes a new instance of the CsvReader class.
        /// </summary>
        /// 
        /// <param name="reader">A <see cref="T:TextReader"/> pointing to the CSV file.</param>
        /// <param name="hasHeaders"><see langword="true"/> if field names are located on the first non commented line, otherwise, <see langword="false"/>.</param>
        /// <param name="bufferSize">The buffer size in bytes.</param>
        /// 
        public CsvReader(TextReader reader, bool hasHeaders, int bufferSize)
        {
            this.init(reader, hasHeaders, bufferSize: bufferSize);
        }

        /// <summary>
        ///   Initializes a new instance of the CsvReader class.
        /// </summary>
        /// 
        /// <param name="reader">A <see cref="T:TextReader"/> pointing to the CSV file.</param>
        /// <param name="hasHeaders"><see langword="true"/> if field names are located on the first non commented line, otherwise, <see langword="false"/>.</param>
        /// <param name="delimiter">The delimiter character separating each field. If set to zero, the
        ///   delimiter will be detected from the file automatically. Default is '\0' (zero).</param>
        /// 
        public CsvReader(TextReader reader, bool hasHeaders, char delimiter)
        {
            this.init(reader, hasHeaders, delimiter);
        }

        /// <summary>
        ///   Initializes a new instance of the CsvReader class.
        /// </summary>
        /// 
        /// <param name="reader">A <see cref="T:TextReader"/> pointing to the CSV file.</param>
        /// <param name="hasHeaders"><see langword="true"/> if field names are located on the first non commented line, otherwise, <see langword="false"/>.</param>
        /// <param name="delimiter">The delimiter character separating each field. If set to zero, the
        ///   delimiter will be detected from the file automatically. Default is '\0' (zero).</param>
        /// <param name="bufferSize">The buffer size in bytes.</param>
        /// 
        public CsvReader(TextReader reader, bool hasHeaders, char delimiter, int bufferSize)
        {
            this.init(reader, hasHeaders, delimiter, bufferSize);
        }

        private void init(TextReader reader, bool hasHeaders, char delimiter = DefaultDelimiter, int bufferSize = DefaultBufferSize)
        {
#if DEBUG && !NETSTANDARD1_4
            this._allocStack = new System.Diagnostics.StackTrace();
#endif

            if (reader == null)
                throw new ArgumentNullException("reader");

            if (bufferSize <= 0)
                throw new ArgumentOutOfRangeException("bufferSize", bufferSize, ExceptionMessage.BufferSizeTooSmall);

            this._bufferSize = bufferSize;

            StreamReader sr = reader as StreamReader;

            if (sr != null)
            {
                Stream stream = sr.BaseStream;

                if (stream.CanSeek)
                {
                    // Handle bad implementations returning 0 or less
                    if (stream.Length > 0)
                        this._bufferSize = (int)Math.Min(bufferSize, stream.Length);
                }
            }

            this._reader = reader;
            this.Delimiter = delimiter;
            this.Quote = DefaultQuote;
            this.Escape = DefaultEscape;
            this.Comment = DefaultComment;

            this.HasHeaders = hasHeaders;
            this.TrimmingOption = ValueTrimmingOptions.UnquotedOnly;
            this.SupportsMultiline = true;
            this.SkipEmptyLines = true;
            this.DefaultHeaderName = "Column";

            this._currentRecordIndex = -1;
            this.DefaultParseErrorAction = ParseErrorAction.RaiseEvent;
        }

        /// <summary>
        ///   Creates a new CsvReader to read from a Web URL.
        /// </summary>
        /// 
        /// <param name="url">The url pointing to the .csv file.</param>
        /// <param name="hasHeaders"><see langword="true"/> if field names are located on the first non commented line, otherwise, <see langword="false"/>.</param>
        /// 
        public static CsvReader FromUrl(string url, bool hasHeaders)
        {
            WebClient client = ExtensionMethods.NewWebClient();
            Console.WriteLine("Downloading {0}", url);
            byte[] bytes = client.DownloadData(url);
            MemoryStream stream = new MemoryStream(bytes);
            return new CsvReader(stream, hasHeaders);
        }

        /// <summary>
        ///   Creates a new CsvReader to read from a string.
        /// </summary>
        /// 
        /// <param name="text">The text containing the fields in the CSV format.</param>
        /// <param name="hasHeaders"><see langword="true"/> if field names are located on the first non commented line, otherwise, <see langword="false"/>.</param>
        /// 
        public static CsvReader FromText(string text, bool hasHeaders)
        {
            return new CsvReader(new StringReader(text), hasHeaders, text.Length);
        }


        /// <summary>
        ///   Occurs when there is an error while parsing the CSV stream.
        /// </summary>
        /// 
        public event EventHandler<ParseErrorEventArgs> ParseError;

        /// <summary>
        ///   Raises the <see cref="M:ParseError"/> event.
        /// </summary>
        /// 
        /// <param name="e">The <see cref="ParseErrorEventArgs"/> that contains the event data.</param>
        /// 
        protected virtual void OnParseError(ParseErrorEventArgs e)
        {
            EventHandler<ParseErrorEventArgs> handler = this.ParseError;

            if (handler != null)
                handler(this, e);
        }




        /// <summary>
        ///   Gets the comment character indicating that 
        ///   a line is commented out. Default is '#'.
        /// </summary>
        /// 
        /// <value>The comment character indicating that a line is commented out.</value>
        /// 
        public char Comment { get; set; }

        /// <summary>
        ///   Gets the escape character letting insert quotation 
        ///   characters inside a quoted field. Default is '"'.
        /// </summary>
        /// 
        /// <value>The escape character letting insert quotation characters inside a quoted field.</value>
        /// 
        public char Escape { get; set; }

        /// <summary>
        ///   Gets the delimiter character separating each field. If
        ///   set to zero ('\0') the reader will try to guess the
        ///   delimiter character automatically from the first line
        ///   of the file.
        /// </summary>
        /// 
        /// <value>The delimiter character separating each field.</value>
        /// 
        public char Delimiter { get; set; }

        /// <summary>
        ///   Gets the quotation character wrapping 
        ///   every field. Default is '"'.
        /// </summary>
        /// 
        /// <value>The quotation character wrapping every field.</value>
        /// 
        public char Quote { get; set; }

        /// <summary>
        ///   Indicates if field names are located on the first non commented line.
        /// </summary>
        /// 
        /// <value><see langword="true"/> if field names are located on the first non commented line, otherwise, <see langword="false"/>.</value>
        /// 
        public bool HasHeaders { get; private set; }

        /// <summary>
        ///   Indicates if spaces at the start and end of a field 
        ///   are trimmed. Default is to trim unquoted fields only.
        /// </summary>
        /// 
        /// <value><see langword="true"/> if spaces at the start and end of a field are trimmed, otherwise, <see langword="false"/>.</value>
        /// 
        public ValueTrimmingOptions TrimmingOption { get; set; }

        /// <summary>
        ///   Gets the buffer size.
        /// </summary>
        /// 
        public int BufferSize { get { return this._bufferSize; } }

        /// <summary>
        ///   Gets or sets the default action to take when a parsing error has occured.
        /// </summary>
        /// 
        /// <value>The default action to take when a parsing error has occured.</value>
        /// 
        public ParseErrorAction DefaultParseErrorAction { get; set; }

        /// <summary>
        ///   Gets or sets the action to take when a field is missing.
        /// </summary>
        /// 
        /// <value>The action to take when a field is missing.</value>
        /// 
        public MissingFieldAction MissingFieldAction { get; set; }

        /// <summary>
        ///   Gets or sets a value indicating if the reader supports multiline fields.
        /// </summary>
        /// 
        /// <value>A value indicating if the reader supports multiline field.</value>
        /// 
        public bool SupportsMultiline { get; set; }

        /// <summary>
        ///   Gets or sets a value indicating if the reader will skip empty lines.
        /// </summary>
        /// 
        /// <value>A value indicating if the reader will skip empty lines.</value>
        /// 
        public bool SkipEmptyLines { get; set; }

        /// <summary>
        ///   Gets or sets the default header name when it is an empty string or only whitespaces.
        ///   The header index will be appended to the specified name. Default is "Column".
        /// </summary>
        /// 
        /// <value>The default header name when it is an empty string or only whitespaces.</value>
        /// 
        public string DefaultHeaderName { get; set; }



        /// <summary>
        ///   Gets the maximum number of fields to retrieve for each record.
        /// </summary>
        /// 
        /// <value>The maximum number of fields to retrieve for each record.</value>
        /// 
        /// <exception cref="T:System.ComponentModel.ObjectDisposedException">
        ///	  The instance has been disposed of.
        /// </exception>
        /// 
        public int FieldCount
        {
            get
            {
                this.EnsureInitialize();
                return this._fieldCount;
            }
        }

        /// <summary>
        ///   Gets a value that indicates whether the current stream position is at the end of the stream.
        /// </summary>
        /// 
        /// <value><see langword="true"/> if the current stream position is at the end of the stream; otherwise <see langword="false"/>.</value>
        /// 
        public virtual bool EndOfStream
        {
            get { return this._eof; }
        }

        /// <summary>
        ///   Gets the field headers.
        /// </summary>
        /// 
        /// <returns>The field headers or an empty array if headers are not supported.</returns>
        /// 
        /// <exception cref="T:System.ComponentModel.ObjectDisposedException">
        ///	  The instance has been disposed of.
        /// </exception>
        /// 
        public string[] GetFieldHeaders()
        {
            this.EnsureInitialize();
            Debug.Assert(this._fieldHeaders != null, "Field headers must be non null.");

            string[] fieldHeaders = new string[this._fieldHeaders.Length];

            for (int i = 0; i < fieldHeaders.Length; i++)
                fieldHeaders[i] = this._fieldHeaders[i];

            return fieldHeaders;
        }

        /// <summary>
        ///   Gets the current record index in the CSV file.
        /// </summary>
        /// 
        /// <value>The current record index in the CSV file.</value>
        /// 
        public virtual long CurrentRecordIndex { get { return this._currentRecordIndex; } }

        /// <summary>
        ///   Indicates if one or more field are missing for the current record.
        ///   Resets after each successful record read.
        /// </summary>
        /// 
        public bool MissingFieldFlag { get { return this._missingFieldFlag; } }

        /// <summary>
        ///   Indicates if a parse error occurred for the current record.
        ///   Resets after each successful record read.
        /// </summary>
        /// 
        public bool ParseErrorFlag
        {
            get { return this._parseErrorFlag; }
        }


#if !NETSTANDARD1_4
        /// <summary>
        ///   Reads the entire stream into a DataTable.
        /// </summary>
        /// 
        /// <returns>A System.DataTable containing the read values.</returns>
        /// 
        public DataTable ToTable()
        {
            DataTable table = new DataTable();
            table.Locale = CultureInfo.InvariantCulture;
            table.Load(this);
            return table;
        }


        /// <summary>
        ///   Reads the entire stream into a DataTable.
        /// </summary>
        /// 
        /// <returns>A System.DataTable containing the read values.</returns>
        /// 
        public DataTable ToTable(params string[] columnNames)
        {
            return new DataView(this.ToTable()).ToTable(false, columnNames);
        }
#endif

        /// <summary>
        ///   Reads the entire stream into a list of records.
        /// </summary>
        /// 
        /// <returns>A list containing all records in the file.</returns>
        /// 
        public List<String[]> ReadToEnd()
        {
            this.EnsureInitialize();

            var lines = new List<string[]>();

            int fieldCount = this.FieldCount;

            while (this.ReadNextRecord())
            {
                var record = new string[fieldCount];

                for (int i = 0; i < this._fieldCount; i++)
                {
                    if (!this._parseErrorFlag)
                        record[i] = this[i];
                }

                lines.Add(record);
            }

            return lines;
        }

        /// <summary>
        ///   Reads the entire stream into a list of records.
        /// </summary>
        /// 
        /// <returns>A list containing all records in the file.</returns>
        /// 
        public String[] ReadLine()
        {
            this.EnsureInitialize();

            int fieldCount = this.FieldCount;

            if (!this.ReadNextRecord())
                throw new InvalidOperationException(ExceptionMessage.NoCurrentRecord);

            var record = new string[fieldCount];

            for (int i = 0; i < this._fieldCount; i++)
            {
                if (!this._parseErrorFlag)
                    record[i] = this[i];
            }

            return record;
        }


        /// <summary>
        ///   Reads the entire stream into a multi-dimensional matrix.
        /// </summary>
        /// 
        public double[,] ToMatrix()
        {
            return this.ToMatrix<double>();
        }

        /// <summary>
        ///   Reads the entire stream into a multi-dimensional matrix.
        /// </summary>
        /// 
        public T[,] ToMatrix<T>()
        {
            var lines = this.ReadToEnd();
            var m = new T[lines.Count, lines[0].Length];

            for (int i = 0; i < lines.Count; i++)
                for (int j = 0; j < lines[i].Length; j++)
                    m[i, j] = (T)System.Convert.ChangeType(lines[i][j], typeof(T));

            return m;
        }

        /// <summary>
        ///   Reads the entire stream into a jagged matrix.
        /// </summary>
        /// 
        public double[][] ToJagged()
        {
            return this.ToJagged<double>();
        }

        /// <summary>
        ///   Reads the entire stream into a jagged matrix.
        /// </summary>
        /// 
        public T[][] ToJagged<T>()
        {
            var lines = this.ReadToEnd();
            var m = new T[lines.Count][];

            for (int i = 0; i < lines.Count; i++)
            {
                m[i] = new T[lines[0].Length];
                for (int j = 0; j < lines[i].Length; j++)
                    m[i][j] = (T)System.Convert.ChangeType(lines[i][j], typeof(T));
            }

            return m;
        }

        /// <summary>
        ///   Gets the field with the specified name and record position. <see cref="M:hasHeaders"/> must be <see langword="true"/>.
        /// </summary>
        /// 
        /// <value>
        ///   The field with the specified name and record position.
        /// </value>
        /// 
        /// <exception cref="T:ArgumentNullException">
        ///   <paramref name="field"/> is <see langword="null"/> or an empty string.
        /// </exception>
        /// <exception cref="T:InvalidOperationException">
        ///   The CSV does not have headers (<see cref="M:HasHeaders"/> property is <see langword="false"/>).
        /// </exception>
        /// <exception cref="T:ArgumentException">
        ///   <paramref name="field"/> not found.
        /// </exception>
        /// <exception cref="T:ArgumentOutOfRangeException">
        ///   Record index must be > 0.
        /// </exception>
        /// <exception cref="T:InvalidOperationException">
        ///   Cannot move to a previous record in forward-only mode.
        /// </exception>
        /// <exception cref="T:EndOfStreamException">
        ///   Cannot read record at <paramref name="record"/>.
        /// </exception>
        /// <exception cref="T:MalformedCsvException">
        ///   The CSV appears to be corrupt at the current position.
        /// </exception>
        /// <exception cref="T:System.ComponentModel.ObjectDisposedException">
        ///   The instance has been disposed of.
        /// </exception>
        /// 
        public string this[int record, string field]
        {
            get
            {
                if (!this.MoveTo(record))
                    throw new InvalidOperationException(string.Format(CultureInfo.InvariantCulture, ExceptionMessage.CannotReadRecordAtIndex, record));

                return this[field];
            }
        }

        /// <summary>
        /// Gets the field at the specified index and record position.
        /// </summary>
        /// <value>
        /// The field at the specified index and record position.
        /// A <see langword="null"/> is returned if the field cannot be found for the record.
        /// </value>
        /// <exception cref="T:ArgumentOutOfRangeException">
        ///		<paramref name="field"/> must be included in [0, <see cref="M:FieldCount"/>[.
        /// </exception>
        /// <exception cref="T:ArgumentOutOfRangeException">
        ///		Record index must be > 0.
        /// </exception>
        /// <exception cref="T:InvalidOperationException">
        ///		Cannot move to a previous record in forward-only mode.
        /// </exception>
        /// <exception cref="T:EndOfStreamException">
        ///		Cannot read record at <paramref name="record"/>.
        /// </exception>
        /// <exception cref="T:MalformedCsvException">
        ///		The CSV appears to be corrupt at the current position.
        /// </exception>
        /// <exception cref="T:System.ComponentModel.ObjectDisposedException">
        ///	The instance has been disposed of.
        /// </exception>
        /// 
        public string this[int record, int field]
        {
            get
            {
                if (!this.MoveTo(record))
                    throw new InvalidOperationException(string.Format(CultureInfo.InvariantCulture, ExceptionMessage.CannotReadRecordAtIndex, record));

                return this[field];
            }
        }

        /// <summary>
        /// Gets the field with the specified name. <see cref="M:hasHeaders"/> must be <see langword="true"/>.
        /// </summary>
        /// <value>
        /// The field with the specified name.
        /// </value>
        /// <exception cref="T:ArgumentNullException">
        ///		<paramref name="field"/> is <see langword="null"/> or an empty string.
        /// </exception>
        /// <exception cref="T:InvalidOperationException">
        ///	The CSV does not have headers (<see cref="M:HasHeaders"/> property is <see langword="false"/>).
        /// </exception>
        /// <exception cref="T:ArgumentException">
        ///		<paramref name="field"/> not found.
        /// </exception>
        /// <exception cref="T:MalformedCsvException">
        ///		The CSV appears to be corrupt at the current position.
        /// </exception>
        /// <exception cref="T:System.ComponentModel.ObjectDisposedException">
        ///	The instance has been disposed of.
        /// </exception>
        public string this[string field]
        {
            get
            {
                if (string.IsNullOrEmpty(field))
                    throw new ArgumentNullException("field");

                if (!this.HasHeaders)
                    throw new InvalidOperationException(ExceptionMessage.NoHeaders);

                int index = this.GetFieldIndex(field);

                if (index < 0)
                    throw new ArgumentException(string.Format(CultureInfo.InvariantCulture, ExceptionMessage.FieldHeaderNotFound, field), "field");

                return this[index];
            }
        }

        /// <summary>
        /// Gets the field at the specified index.
        /// </summary>
        /// <value>The field at the specified index.</value>
        /// <exception cref="T:ArgumentOutOfRangeException">
        ///		<paramref name="field"/> must be included in [0, <see cref="M:FieldCount"/>[.
        /// </exception>
        /// <exception cref="T:InvalidOperationException">
        ///		No record read yet. Call ReadLine() first.
        /// </exception>
        /// <exception cref="T:MalformedCsvException">
        ///		The CSV appears to be corrupt at the current position.
        /// </exception>
        /// <exception cref="T:System.ComponentModel.ObjectDisposedException">
        ///	The instance has been disposed of.
        /// </exception>
        public virtual string this[int field]
        {
            get { return this.ReadField(field, false, false); }
        }



        /// <summary>
        /// Gets the field index for the provided header.
        /// </summary>
        /// <param name="header">The header to look for.</param>
        /// <returns>The field index for the provided header. -1 if not found.</returns>
        /// <exception cref="T:System.ComponentModel.ObjectDisposedException">
        ///	The instance has been disposed of.
        /// </exception>
        public int GetFieldIndex(string header)
        {
            this.EnsureInitialize();

            int index;

            if (this._fieldHeaderIndexes != null && this._fieldHeaderIndexes.TryGetValue(header, out index))
                return index;
            else
                return -1;
        }




        /// <summary>
        ///   Copies the field array of the current record to a one-dimensional array, starting at the beginning of the target array.
        /// </summary>
        /// 
        /// <param name="array"> The one-dimensional <see cref="T:Array"/> that is the destination of the fields of the current record.</param>
        /// 
        /// <exception cref="T:ArgumentNullException">
        ///		<paramref name="array"/> is <see langword="null"/>.
        /// </exception>
        ///
        public void CopyCurrentRecordTo(string[] array)
        {
            this.CopyCurrentRecordTo(array, 0);
        }

        /// <summary>
        ///   Copies the field array of the current record to a one-dimensional array, starting at the beginning of the target array.
        /// </summary>
        /// 
        /// <param name="array"> The one-dimensional <see cref="T:Array"/> that is the destination of the fields of the current record.</param>
        /// <param name="index">The zero-based index in <paramref name="array"/> at which copying begins.</param>
        /// 
        /// <exception cref="T:ArgumentNullException">
        ///		<paramref name="array"/> is <see langword="null"/>.
        /// </exception>
        /// <exception cref="T:ArgumentOutOfRangeException">
        ///		<paramref name="index"/> is les than zero or is equal to or greater than the length <paramref name="array"/>. 
        /// </exception>
        /// <exception cref="InvalidOperationException">
        ///	No current record.
        /// </exception>
        /// <exception cref="ArgumentException">
        ///		The number of fields in the record is greater than the available space from <paramref name="index"/> to the end of <paramref name="array"/>.
        /// </exception>
        /// 
        public void CopyCurrentRecordTo(string[] array, int index)
        {
            if (array == null)
                throw new ArgumentNullException("array");

            if (index < 0 || index >= array.Length)
                throw new ArgumentOutOfRangeException("index", index, string.Empty);

            if (this._currentRecordIndex < 0 || !this._initialized)
                throw new InvalidOperationException(ExceptionMessage.NoCurrentRecord);

            if (array.Length - index < this._fieldCount)
                throw new ArgumentException(ExceptionMessage.NotEnoughSpaceInArray, "array");

            for (int i = 0; i < this._fieldCount; i++)
            {
                if (this._parseErrorFlag)
                    array[index + i] = null;
                else
                    array[index + i] = this[i];
            }
        }



        /// <summary>
        ///   Gets the current raw CSV data.
        /// </summary>
        /// 
        /// <remarks>Used for exception handling purposes.</remarks>
        /// 
        /// <returns>The current raw CSV data.</returns>
        /// 
        private string GetCurrentRawData()
        {
            if (this._buffer != null && this._bufferLength > 0)
                return new string(this._buffer, 0, this._bufferLength);

            return string.Empty;
        }


        /// <summary>
        /// Ensures that the reader is initialized.
        /// </summary>
        private void EnsureInitialize()
        {
            if (!this._initialized)
                this.ReadNextRecord(true, false);

            Debug.Assert(this._fieldHeaders != null);
            Debug.Assert(this._fieldHeaders.Length > 0 || (this._fieldHeaders.Length == 0 && this._fieldHeaderIndexes == null));
        }

        /// <summary>
        /// Indicates whether the specified Unicode character is categorized as white space.
        /// </summary>
        /// <param name="c">A Unicode character.</param>
        /// <returns><see langword="true"/> if <paramref name="c"/> is white space; otherwise, <see langword="false"/>.</returns>
        private bool IsWhiteSpace(char c)
        {
            // Handle cases where the delimiter is a whitespace (e.g. tab)
            if (c == this.Delimiter)
                return false;

            // See char.IsLatin1(char c) in Reflector
            if (c <= '\x00ff')
                return (c == ' ' || c == '\t');

            return (System.Globalization.CharUnicodeInfo.GetUnicodeCategory(c) == System.Globalization.UnicodeCategory.SpaceSeparator);
        }



        /// <summary>
        ///   Moves to the specified record index.
        /// </summary>
        /// 
        /// <param name="record">The record index.</param>
        /// 
        /// <returns><c>true</c> if the operation was successful; otherwise, <c>false</c>.</returns>
        /// 
        /// <exception cref="T:System.ComponentModel.ObjectDisposedException">
        ///   The instance has been disposed of.
        /// </exception>
        /// 
        public virtual bool MoveTo(long record)
        {
            if (record < this._currentRecordIndex)
                return false;

            // Get number of record to read
            long offset = record - this._currentRecordIndex;

            while (offset > 0)
            {
                if (!this.ReadNextRecord())
                    return false;

                offset--;
            }

            return true;
        }

        /// <summary>
        ///   Reads the next record.
        /// </summary>
        /// 
        /// <returns><see langword="true"/> if a record has been successfully reads; otherwise, <see langword="false"/>.</returns>
        /// 
        /// <exception cref="T:System.ComponentModel.ObjectDisposedException">
        ///   The instance has been disposed of.
        /// </exception>
        /// 
        public bool ReadNextRecord()
        {
            return this.ReadNextRecord(false, false);
        }


        /// <summary>
        ///   Parses a new line delimiter.
        /// </summary>
        /// 
        /// <param name="pos">The starting position of the parsing. Will contain the resulting end position.</param>
        /// 
        /// <returns><see langword="true"/> if a new line delimiter was found; otherwise, <see langword="false"/>.</returns>
        /// 
        /// <exception cref="T:System.ComponentModel.ObjectDisposedException">
        ///   The instance has been disposed of.
        /// </exception>
        /// 
        private bool ParseNewLine(ref int pos)
        {
            Debug.Assert(pos <= this._bufferLength);

            // Check if already at the end of the buffer
            if (pos == this._bufferLength)
            {
                pos = 0;

                if (!this.ReadBuffer())
                    return false;
            }

            char c = this._buffer[pos];

            // Treat \r as new line only if it's not the delimiter

            if (c == '\r' && this.Delimiter != '\r')
            {
                pos++;

                // Skip following \n (if there is one)

                if (pos < this._bufferLength)
                {
                    if (this._buffer[pos] == '\n')
                        pos++;
                }
                else
                {
                    if (this.ReadBuffer())
                    {
                        if (this._buffer[0] == '\n')
                            pos = 1;
                        else
                            pos = 0;
                    }
                }

                if (pos >= this._bufferLength)
                {
                    this.ReadBuffer();
                    pos = 0;
                }

                return true;
            }
            else if (c == '\n')
            {
                pos++;

                if (pos >= this._bufferLength)
                {
                    this.ReadBuffer();
                    pos = 0;
                }

                return true;
            }

            return false;
        }

        /// <summary>
        ///   Determines whether the character at the specified position is a new line delimiter.
        /// </summary>
        /// 
        /// <param name="pos">The position of the character to verify.</param>
        /// 
        /// <returns>
        ///   <see langword="true"/> if the character at the specified position is a new line delimiter; otherwise, <see langword="false"/>.
        /// </returns>
        /// 
        private bool IsNewLine(int pos)
        {
            Debug.Assert(pos < this._bufferLength);

            char c = this._buffer[pos];

            if (c == '\n')
                return true;

            if (c == '\r' && this.Delimiter != '\r')
                return true;

            return false;
        }


        /// <summary>
        /// Fills the buffer with data from the reader.
        /// </summary>
        /// <returns><see langword="true"/> if data was successfully read; otherwise, <see langword="false"/>.</returns>
        /// <exception cref="T:System.ComponentModel.ObjectDisposedException">
        ///	The instance has been disposed of.
        /// </exception>
        private bool ReadBuffer()
        {
            if (this._eof)
                return false;

            this.CheckDisposed();

            this._bufferLength = this._reader.Read(this._buffer, 0, this._bufferSize);

            if (this._bufferLength > 0)
                return true;

            this._eof = true;
            this._buffer = null;

            return false;
        }

        private void GuessDelimiterFromBuffer()
        {
            if (this._buffer == null)
                return;

            // Find first unquoted non-empty end line
            int lineEndIndex = this._buffer.Length - 1;

            int start = 0;

            while (this._buffer[start] == '\r' || this._buffer[start] == '\n') start++;

            bool inquote = false;
            for (int i = start; i < this._buffer.Length; i++)
            {
                if (this._buffer[i] == this.Quote)
                    inquote = !inquote;

                if (!inquote && this._buffer[i] == '\n')
                {
                    lineEndIndex = i;
                    break;
                }
            }

            int imax = 0;
            int max = 0;
            inquote = false;

            for (int i = 0; i < this.delimiters.Length; i++)
            {
                char delimiter = this.delimiters[i];

                int count = 0;
                for (int j = start; j < lineEndIndex; j++)
                {
                    if (this._buffer[j] == this.Quote)
                        inquote = !inquote;

                    if (!inquote && this._buffer[j] == delimiter)
                        count++;
                }

                if (count > max)
                {
                    max = count;
                    imax = i;
                }
            }

            this.Delimiter = this.delimiters[imax];
        }

        /// <summary>
        /// Reads the field at the specified index.
        /// Any unread fields with an inferior index will also be read as part of the required parsing.
        /// </summary>
        /// <param name="field">The field index.</param>
        /// <param name="initializing">Indicates if the reader is currently initializing.</param>
        /// <param name="discardValue">Indicates if the value(s) are discarded.</param>
        /// <returns>
        /// The field at the specified index. 
        /// A <see langword="null"/> indicates that an error occured or that the last field has been reached during initialization.
        /// </returns>
        /// <exception cref="ArgumentOutOfRangeException">
        ///		<paramref name="field"/> is out of range.
        /// </exception>
        /// <exception cref="InvalidOperationException">
        ///		There is no current record.
        /// </exception>
        /// <exception cref="MissingFieldCsvException">
        ///		The CSV data appears to be missing a field.
        /// </exception>
        /// <exception cref="MalformedCsvException">
        ///		The CSV data appears to be malformed.
        /// </exception>
        /// <exception cref="T:System.ComponentModel.ObjectDisposedException">
        ///	The instance has been disposed of.
        /// </exception>
        private string ReadField(int field, bool initializing, bool discardValue)
        {
            if (!initializing)
            {
                if (field < 0 || field >= this._fieldCount)
                    throw new ArgumentOutOfRangeException("field", field, string.Format(CultureInfo.InvariantCulture, ExceptionMessage.FieldIndexOutOfRange, field));

                if (this._currentRecordIndex < 0)
                    throw new InvalidOperationException(ExceptionMessage.NoCurrentRecord);

                // Directly return field if cached
                if (this._fields[field] != null)
                    return this._fields[field];

                if (this._missingFieldFlag)
                    return this.HandleMissingField(null, field, ref this._nextFieldStart);
            }

            this.CheckDisposed();

            int index = this._nextFieldIndex;

            while (index < field + 1)
            {
                // Handle case where stated start of field is past buffer
                // This can occur because _nextFieldStart is simply 1 + last char position of previous field
                if (this._nextFieldStart == this._bufferLength)
                {
                    this._nextFieldStart = 0;

                    // Possible EOF will be handled later (see Handle_EOF1)
                    this.ReadBuffer();
                }

                if (initializing && (this.Delimiter == DefaultDelimiter)) // do not guess if delimiter was already set
                {
                    this.GuessDelimiterFromBuffer();
                }

                string value = null;

                if (this._missingFieldFlag)
                {
                    value = this.HandleMissingField(value, index, ref this._nextFieldStart);
                }
                else if (this._nextFieldStart == this._bufferLength)
                {
                    // Handle_EOF1: Handle EOF here

                    // If current field is the requested field, then the value of the field is "" as in "f1,f2,f3,(\s*)"
                    // otherwise, the CSV is malformed

                    if (index == field)
                    {
                        if (!discardValue)
                        {
                            value = string.Empty;
                            this._fields[index] = value;
                        }

                        this._missingFieldFlag = true;
                    }
                    else
                    {
                        value = this.HandleMissingField(value, index, ref this._nextFieldStart);
                    }
                }
                else
                {
                    // Trim spaces at start
                    if ((this.TrimmingOption & ValueTrimmingOptions.UnquotedOnly) != 0)
                        this.SkipWhiteSpaces(ref this._nextFieldStart);

                    if (this._eof)
                    {
                        value = string.Empty;
                        this._fields[field] = value;

                        if (field < this._fieldCount)
                            this._missingFieldFlag = true;
                    }
                    else if (this._buffer[this._nextFieldStart] != this.Quote)
                    {
                        // Non-quoted field

                        int start = this._nextFieldStart;
                        int pos = this._nextFieldStart;

                        for (; ; )
                        {
                            while (pos < this._bufferLength)
                            {
                                char c = this._buffer[pos];

                                if (c == this.Delimiter)
                                {
                                    this._nextFieldStart = pos + 1;

                                    break;
                                }
                                else if (c == '\r' || c == '\n')
                                {
                                    this._nextFieldStart = pos;
                                    this._eol = true;

                                    break;
                                }
                                else
                                {
                                    pos++;
                                }
                            }

                            if (pos < this._bufferLength)
                                break;

                            if (!discardValue)
                                value += new string(this._buffer, start, pos - start);

                            start = 0;
                            pos = 0;
                            this._nextFieldStart = 0;

                            if (!this.ReadBuffer())
                                break;
                        }

                        if (!discardValue)
                        {
                            if ((this.TrimmingOption & ValueTrimmingOptions.UnquotedOnly) == 0)
                            {
                                if (!this._eof && pos > start)
                                    value += new string(this._buffer, start, pos - start);
                            }
                            else
                            {
                                if (!this._eof && pos > start)
                                {
                                    // Do the trimming
                                    pos--;
                                    while (pos > -1 && this.IsWhiteSpace(this._buffer[pos]))
                                        pos--;
                                    pos++;

                                    if (pos > 0)
                                        value += new string(this._buffer, start, pos - start);
                                }
                                else
                                {
                                    pos = -1;
                                }

                                // If pos <= 0, that means the trimming went past buffer start,
                                // and the concatenated value needs to be trimmed too.
                                if (pos <= 0)
                                {
                                    pos = (value == null ? -1 : value.Length - 1);

                                    // Do the trimming
                                    while (pos > -1 && this.IsWhiteSpace(value[pos]))
                                        pos--;

                                    pos++;

                                    if (pos > 0 && pos != value.Length)
                                        value = value.Substring(0, pos);
                                }
                            }

                            if (value == null)
                                value = string.Empty;
                        }

                        if (this._eol || this._eof)
                        {
                            this._eol = this.ParseNewLine(ref this._nextFieldStart);

                            // Reaching a new line is ok as long as the parser is initializing or it is the last field
                            if (!initializing && index != this._fieldCount - 1)
                            {
                                if (value != null && value.Length == 0)
                                    value = null;

                                value = this.HandleMissingField(value, index, ref this._nextFieldStart);
                            }
                        }

                        if (!discardValue)
                            this._fields[index] = value;
                    }
                    else
                    {
                        // Quoted field

                        // Skip quote
                        int start = this._nextFieldStart + 1;
                        int pos = start;

                        bool quoted = true;
                        bool escaped = false;

                        if ((this.TrimmingOption & ValueTrimmingOptions.QuotedOnly) != 0)
                        {
                            this.SkipWhiteSpaces(ref start);
                            pos = start;
                        }

                        for (; ; )
                        {
                            while (pos < this._bufferLength)
                            {
                                char c = this._buffer[pos];

                                if (escaped)
                                {
                                    escaped = false;
                                    start = pos;
                                }
                                // IF current char is escape AND (escape and quote are different OR next char is a quote)
                                else if (c == this.Escape && (this.Escape != this.Quote || (pos + 1 < this._bufferLength && this._buffer[pos + 1] == this.Quote) || (pos + 1 == this._bufferLength && this._reader.Peek() == this.Quote)))
                                {
                                    if (!discardValue)
                                        value += new string(this._buffer, start, pos - start);

                                    escaped = true;
                                }
                                else if (c == this.Quote)
                                {
                                    quoted = false;
                                    break;
                                }

                                pos++;
                            }

                            if (!quoted)
                                break;

                            if (!discardValue && !escaped)
                                value += new string(this._buffer, start, pos - start);

                            start = 0;
                            pos = 0;
                            this._nextFieldStart = 0;

                            if (!this.ReadBuffer())
                            {
                                this.HandleParseError(new MalformedCsvException(this.GetCurrentRawData(), this._nextFieldStart, Math.Max(0, this._currentRecordIndex), index), ref this._nextFieldStart);
                                return null;
                            }
                        }

                        if (!this._eof)
                        {
                            // Append remaining parsed buffer content
                            if (!discardValue && pos > start)
                                value += new string(this._buffer, start, pos - start);

                            if (!discardValue && value != null && (this.TrimmingOption & ValueTrimmingOptions.QuotedOnly) != 0)
                            {
                                int newLength = value.Length;
                                while (newLength > 0 && this.IsWhiteSpace(value[newLength - 1]))
                                    newLength--;

                                if (newLength < value.Length)
                                    value = value.Substring(0, newLength);
                            }

                            // Skip quote
                            this._nextFieldStart = pos + 1;

                            // Skip whitespaces between the quote and the delimiter/eol
                            this.SkipWhiteSpaces(ref this._nextFieldStart);

                            // Skip delimiter
                            bool delimiterSkipped;
                            if (this._nextFieldStart < this._bufferLength && this._buffer[this._nextFieldStart] == this.Delimiter)
                            {
                                this._nextFieldStart++;
                                delimiterSkipped = true;
                            }
                            else
                            {
                                delimiterSkipped = false;
                            }

                            // Skip new line delimiter if initializing or last field
                            // (if the next field is missing, it will be caught when parsed)
                            if (!this._eof && !delimiterSkipped && (initializing || index == this._fieldCount - 1))
                                this._eol = this.ParseNewLine(ref this._nextFieldStart);

                            // If no delimiter is present after the quoted field and it is not the last field, then
                            // either we might be using the wrong the wrong field delimiter or it's a parsing error
                            if (!delimiterSkipped && !this._eof && !(this._eol || this.IsNewLine(this._nextFieldStart)))
                                this.HandleParseError(new MalformedCsvException(this.GetCurrentRawData(), this._nextFieldStart, Math.Max(0, this._currentRecordIndex), index), ref this._nextFieldStart);

                        }

                        if (!discardValue)
                        {
                            if (value == null)
                                value = string.Empty;

                            this._fields[index] = value;
                        }
                    }
                }

                this._nextFieldIndex = Math.Max(index + 1, this._nextFieldIndex);

                if (index == field)
                {
                    // If initializing, return null to signify the last field has been reached

                    if (initializing)
                    {
                        if (this._eol || this._eof)
                            return null;

                        return string.IsNullOrEmpty(value) ? string.Empty : value;
                    }

                    return value;
                }

                index++;
            }

            // Getting here is bad ...
            this.HandleParseError(new MalformedCsvException(this.GetCurrentRawData(), this._nextFieldStart, Math.Max(0, this._currentRecordIndex), index), ref this._nextFieldStart);
            return null;
        }



        /// <summary>
        /// Reads the next record.
        /// </summary>
        /// <param name="onlyReadHeaders">
        /// Indicates if the reader will proceed to the next record after having read headers.
        /// <see langword="true"/> if it stops after having read headers; otherwise, <see langword="false"/>.
        /// </param>
        /// <param name="skipToNextLine">
        /// Indicates if the reader will skip directly to the next line without parsing the current one. 
        /// To be used when an error occurs.
        /// </param>
        /// <returns><see langword="true"/> if a record has been successfully reads; otherwise, <see langword="false"/>.</returns>
        /// <exception cref="T:System.ComponentModel.ObjectDisposedException">
        ///	The instance has been disposed of.
        /// </exception>
        protected virtual bool ReadNextRecord(bool onlyReadHeaders, bool skipToNextLine)
        {
            if (this._eof)
            {
                if (this._firstRecordInCache)
                {
                    this._firstRecordInCache = false;
                    this._currentRecordIndex++;

                    return true;
                }

                return false;
            }

            this.CheckDisposed();

            if (!this._initialized)
            {
                this._buffer = new char[this._bufferSize];

                // will be replaced if and when headers are read
                this._fieldHeaders = new string[0];

                if (!this.ReadBuffer())
                    return false;

                if (!this.SkipEmptyAndCommentedLines(ref this._nextFieldStart))
                    return false;

                // Keep growing _fields array until the last field has been found
                // and then resize it to its final correct size

                this._fieldCount = 0;
                this._fields = new string[16];

                while (this.ReadField(this._fieldCount, true, false) != null)
                {
                    if (this._parseErrorFlag)
                    {
                        this._fieldCount = 0;
                        Array.Clear(this._fields, 0, this._fields.Length);
                        this._parseErrorFlag = false;
                        this._nextFieldIndex = 0;
                    }
                    else
                    {
                        this._fieldCount++;

                        if (this._fieldCount == this._fields.Length)
                            Array.Resize<string>(ref this._fields, (this._fieldCount + 1) * 2);
                    }
                }

                if (this.Delimiter == ' ' && this._fields[this._fieldCount] == String.Empty)
                {
                    // discard the last empty elements in this case
                    this._fieldCount--;
                }

                // _fieldCount contains the last field index, but it must contains the field count,
                // so increment by 1
                this._fieldCount++;

                if (this._fields.Length != this._fieldCount)
                    Array.Resize<string>(ref this._fields, this._fieldCount);

                this._initialized = true;

                // If headers are present, call ReadNextRecord again
                if (this.HasHeaders)
                {
                    // Don't count first record as it was the headers
                    this._currentRecordIndex = -1;

                    this._firstRecordInCache = false;

                    this._fieldHeaders = new string[this._fieldCount];
                    this._fieldHeaderIndexes = new Dictionary<string, int>(this._fieldCount, _fieldHeaderComparer);

                    for (int i = 0; i < this._fields.Length; i++)
                    {
                        string headerName = this._fields[i];
                        if (string.IsNullOrEmpty(headerName) || headerName.Trim().Length == 0)
                            headerName = this.DefaultHeaderName + i.ToString();

                        this._fieldHeaders[i] = headerName;
                        this._fieldHeaderIndexes.Add(headerName, i);
                    }

                    // Proceed to first record
                    if (!onlyReadHeaders)
                    {
                        // Calling again ReadNextRecord() seems to be simpler, 
                        // but in fact would probably cause many subtle bugs because a derived class does not expect a recursive behavior
                        // so simply do what is needed here and no more.

                        if (!this.SkipEmptyAndCommentedLines(ref this._nextFieldStart))
                            return false;

                        Array.Clear(this._fields, 0, this._fields.Length);
                        this._nextFieldIndex = 0;
                        this._eol = false;

                        this._currentRecordIndex++;
                        return true;
                    }
                }
                else
                {
                    if (onlyReadHeaders)
                    {
                        this._firstRecordInCache = true;
                        this._currentRecordIndex = -1;
                    }
                    else
                    {
                        this._firstRecordInCache = false;
                        this._currentRecordIndex = 0;
                    }
                }
            }
            else
            {
                if (skipToNextLine)
                {
                    this.SkipToNextLine(ref this._nextFieldStart);
                }
                else if (this._currentRecordIndex > -1 && !this._missingFieldFlag)
                {
                    // If not already at end of record, move there
                    if (!this._eol && !this._eof)
                    {
                        if (!this.SupportsMultiline)
                        {
                            this.SkipToNextLine(ref this._nextFieldStart);
                        }
                        else
                        {
                            // a dirty trick to handle the case where extra fields are present
                            while (this.ReadField(this._nextFieldIndex, true, true) != null)
                            {
                            }
                        }
                    }
                }

                if (!this._firstRecordInCache && !this.SkipEmptyAndCommentedLines(ref this._nextFieldStart))
                    return false;

                if (this.HasHeaders || !this._firstRecordInCache)
                    this._eol = false;

                // Check to see if the first record is in cache.
                // This can happen when initializing a reader with no headers
                // because one record must be read to get the field count automatically
                if (this._firstRecordInCache)
                {
                    this._firstRecordInCache = false;
                }
                else
                {
                    Array.Clear(this._fields, 0, this._fields.Length);
                    this._nextFieldIndex = 0;
                }

                this._missingFieldFlag = false;
                this._parseErrorFlag = false;
                this._currentRecordIndex++;
            }

            return true;
        }


        /// <summary>
        /// Skips empty and commented lines.
        /// If the end of the buffer is reached, its content be discarded and filled again from the reader.
        /// </summary>
        /// <param name="pos">
        /// The position in the buffer where to start parsing. 
        /// Will contains the resulting position after the operation.
        /// </param>
        /// <returns><see langword="true"/> if the end of the reader has not been reached; otherwise, <see langword="false"/>.</returns>
        /// <exception cref="T:System.ComponentModel.ObjectDisposedException">
        ///	The instance has been disposed of.
        /// </exception>
        private bool SkipEmptyAndCommentedLines(ref int pos)
        {
            if (pos < this._bufferLength)
                this.DoSkipEmptyAndCommentedLines(ref pos);

            while (pos >= this._bufferLength && !this._eof)
            {
                if (this.ReadBuffer())
                {
                    pos = 0;
                    this.DoSkipEmptyAndCommentedLines(ref pos);
                }
                else
                {
                    return false;
                }
            }

            return !this._eof;
        }

        /// <summary>
        /// <para>Worker method.</para>
        /// <para>Skips empty and commented lines.</para>
        /// </summary>
        /// <param name="pos">
        /// The position in the buffer where to start parsing. 
        /// Will contains the resulting position after the operation.
        /// </param>
        /// <exception cref="T:System.ComponentModel.ObjectDisposedException">
        ///	The instance has been disposed of.
        /// </exception>
        private void DoSkipEmptyAndCommentedLines(ref int pos)
        {
            while (pos < this._bufferLength)
            {
                if (this._buffer[pos] == this.Comment)
                {
                    pos++;
                    this.SkipToNextLine(ref pos);
                }
                else if (this.SkipEmptyLines && this.ParseNewLine(ref pos))
                {
                    continue;
                }

                break;
            }
        }


        /// <summary>
        /// Skips whitespace characters.
        /// </summary>
        /// <param name="pos">The starting position of the parsing. Will contain the resulting end position.</param>
        /// <returns><see langword="true"/> if the end of the reader has not been reached; otherwise, <see langword="false"/>.</returns>
        /// <exception cref="T:System.ComponentModel.ObjectDisposedException">
        ///	The instance has been disposed of.
        /// </exception>
        private bool SkipWhiteSpaces(ref int pos)
        {
            for (; ; )
            {
                while (pos < this._bufferLength && this.IsWhiteSpace(this._buffer[pos]))
                    pos++;

                if (pos < this._bufferLength)
                    break;
                else
                {
                    pos = 0;

                    if (!this.ReadBuffer())
                        return false;
                }
            }

            return true;
        }


        /// <summary>
        /// Skips ahead to the next NewLine character.
        /// If the end of the buffer is reached, its content be discarded and filled again from the reader.
        /// </summary>
        /// <param name="pos">
        /// The position in the buffer where to start parsing. 
        /// Will contain the resulting position after the operation.
        /// </param>
        /// <returns><see langword="true"/> if the end of the reader has not been reached; otherwise, <see langword="false"/>.</returns>
        /// <exception cref="T:System.ComponentModel.ObjectDisposedException">
        ///	The instance has been disposed of.
        /// </exception>
        private bool SkipToNextLine(ref int pos)
        {
            // ((pos = 0) == 0) is a little trick to reset position inline
            while ((pos < this._bufferLength || (this.ReadBuffer() && ((pos = 0) == 0))) && !this.ParseNewLine(ref pos))
                pos++;

            return !this._eof;
        }


        /// <summary>
        /// Handles a parsing error.
        /// </summary>
        /// <param name="error">The parsing error that occured.</param>
        /// <param name="pos">The current position in the buffer.</param>
        /// <exception cref="ArgumentNullException">
        ///	<paramref name="error"/> is <see langword="null"/>.
        /// </exception>
        private void HandleParseError(MalformedCsvException error, ref int pos)
        {
            if (error == null)
                throw new ArgumentNullException("error");

            this._parseErrorFlag = true;

            switch (this.DefaultParseErrorAction)
            {
                case ParseErrorAction.ThrowException:
                    throw error;

                case ParseErrorAction.RaiseEvent:
                    ParseErrorEventArgs e = new ParseErrorEventArgs(error, ParseErrorAction.ThrowException);
                    this.OnParseError(e);

                    switch (e.Action)
                    {
                        case ParseErrorAction.ThrowException:
                            throw e.Error;

                        case ParseErrorAction.RaiseEvent:
                            throw new InvalidOperationException(string.Format(CultureInfo.InvariantCulture, ExceptionMessage.ParseErrorActionInvalidInsideParseErrorEvent, e.Action), e.Error);

                        case ParseErrorAction.AdvanceToNextLine:
                            // already at EOL when fields are missing, so don't skip to next line in that case
                            if (!this._missingFieldFlag && pos >= 0)
                                this.SkipToNextLine(ref pos);
                            break;

                        default:
                            throw new NotSupportedException(string.Format(CultureInfo.InvariantCulture, ExceptionMessage.ParseErrorActionNotSupported, e.Action), e.Error);
                    }
                    break;

                case ParseErrorAction.AdvanceToNextLine:
                    // already at EOL when fields are missing, so don't skip to next line in that case
                    if (!this._missingFieldFlag && pos >= 0)
                        this.SkipToNextLine(ref pos);
                    break;

                default:
                    throw new NotSupportedException(string.Format(CultureInfo.InvariantCulture, ExceptionMessage.ParseErrorActionNotSupported, this.DefaultParseErrorAction), error);
            }
        }


        /// <summary>
        /// Handles a missing field error.
        /// </summary>
        /// <param name="value">The partially parsed value, if available.</param>
        /// <param name="fieldIndex">The missing field index.</param>
        /// <param name="currentPosition">The current position in the raw data.</param>
        /// <returns>
        /// The resulting value according to <see cref="M:MissingFieldAction"/>.
        /// If the action is set to <see cref="T:MissingFieldAction.TreatAsParseError"/>,
        /// then the parse error will be handled according to <see cref="DefaultParseErrorAction"/>.
        /// </returns>
        private string HandleMissingField(string value, int fieldIndex, ref int currentPosition)
        {
            if (fieldIndex < 0 || fieldIndex >= this._fieldCount)
                throw new ArgumentOutOfRangeException("fieldIndex", fieldIndex, string.Format(CultureInfo.InvariantCulture, ExceptionMessage.FieldIndexOutOfRange, fieldIndex));

            this._missingFieldFlag = true;

            for (int i = fieldIndex + 1; i < this._fieldCount; i++)
                this._fields[i] = null;

            if (value != null)
                return value;
            else
            {
                switch (this.MissingFieldAction)
                {
                    case MissingFieldAction.ParseError:
                        this.HandleParseError(new MissingFieldCsvException(this.GetCurrentRawData(), currentPosition, Math.Max(0, this._currentRecordIndex), fieldIndex), ref currentPosition);
                        return value;

                    case MissingFieldAction.ReplaceByEmpty:
                        return string.Empty;

                    case MissingFieldAction.ReplaceByNull:
                        return null;

                    default:
                        throw new NotSupportedException(String.Format(CultureInfo.InvariantCulture,
                            ExceptionMessage.MissingFieldActionNotSupported, this.MissingFieldAction));
                }
            }
        }

#if !NETSTANDARD1_4
        /// <summary>
        /// Validates the state of the data reader.
        /// </summary>
        /// <param name="validations">The validations to accomplish.</param>
        /// <exception cref="InvalidOperationException">
        ///	No current record.
        /// </exception>
        /// <exception cref="InvalidOperationException">
        ///	This operation is invalid when the reader is closed.
        /// </exception>
        private void ValidateDataReader(DataReaderValidations validations)
        {
            if ((validations & DataReaderValidations.IsInitialized) != 0 && !this._initialized)
                throw new InvalidOperationException(ExceptionMessage.NoCurrentRecord);

            if ((validations & DataReaderValidations.IsNotClosed) != 0 && this._isDisposed)
                throw new InvalidOperationException(ExceptionMessage.ReaderClosed);
        }
#endif

        /// <summary>
        /// Copy the value of the specified field to an array.
        /// </summary>
        /// <param name="field">The index of the field.</param>
        /// <param name="fieldOffset">The offset in the field value.</param>
        /// <param name="destinationArray">The destination array where the field value will be copied.</param>
        /// <param name="destinationOffset">The destination array offset.</param>
        /// <param name="length">The number of characters to copy from the field value.</param>
        /// <returns></returns>
        private long CopyFieldToArray(int field, long fieldOffset, Array destinationArray, int destinationOffset, int length)
        {
            this.EnsureInitialize();

            if (field < 0 || field >= this._fieldCount)
                throw new ArgumentOutOfRangeException("field", field, string.Format(CultureInfo.InvariantCulture, ExceptionMessage.FieldIndexOutOfRange, field));

            if (fieldOffset < 0 || fieldOffset >= int.MaxValue)
                throw new ArgumentOutOfRangeException("fieldOffset");

            // Array.Copy(...) will do the remaining argument checks

            if (length == 0)
                return 0;

            string value = this[field];

            if (value == null)
                value = string.Empty;

            Debug.Assert(fieldOffset < int.MaxValue);

            Debug.Assert(destinationArray.GetType() == typeof(char[]) || destinationArray.GetType() == typeof(byte[]));

            if (destinationArray.GetType() == typeof(char[]))
                Array.Copy(value.ToCharArray((int)fieldOffset, length), 0, destinationArray, destinationOffset, length);
            else
            {
                char[] chars = value.ToCharArray((int)fieldOffset, length);
                byte[] source = new byte[chars.Length];
                ;

                for (int i = 0; i < chars.Length; i++)
                    source[i] = Convert.ToByte(chars[i]);

                Array.Copy(source, 0, destinationArray, destinationOffset, length);
            }

            return length;
        }




        /// <summary>
        /// Returns an <see cref="T:RecordEnumerator"/>  that can iterate through CSV records.
        /// </summary>
        /// <returns>An <see cref="T:RecordEnumerator"/>  that can iterate through CSV records.</returns>
        /// <exception cref="T:System.ComponentModel.ObjectDisposedException">
        ///	The instance has been disposed of.
        /// </exception>
        public RecordEnumerator GetEnumerator()
        {
            return new RecordEnumerator(this);
        }

        IEnumerator<string[]> IEnumerable<string[]>.GetEnumerator()
        {
            return this.GetEnumerator();
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return this.GetEnumerator();
        }


#if DEBUG && !NETSTANDARD1_4
        /// <summary>
        /// Contains the stack when the object was allocated.
        /// </summary>
        private System.Diagnostics.StackTrace _allocStack;
#endif



        /// <summary>
        /// Contains the disposed status flag.
        /// </summary>
        private bool _isDisposed = false;

        /// <summary>
        /// Contains the locking object for multi-threading purpose.
        /// </summary>
        private readonly object _lock = new object();

        /// <summary>
        /// Occurs when the instance is disposed of.
        /// </summary>
        public event EventHandler Disposed;

        /// <summary>
        /// Gets a value indicating whether the instance has been disposed of.
        /// </summary>
        /// <value>
        /// 	<see langword="true"/> if the instance has been disposed of; otherwise, <see langword="false"/>.
        /// </value>
#if !NETSTANDARD1_4
        [System.ComponentModel.Browsable(false)]
#endif
        public bool IsDisposed
        {
            get { return this._isDisposed; }
        }

        /// <summary>
        /// Raises the <see cref="M:Disposed"/> event.
        /// </summary>
        /// <param name="e">A <see cref="T:System.EventArgs"/> that contains the event data.</param>
        protected virtual void OnDisposed(EventArgs e)
        {
            EventHandler handler = this.Disposed;

            if (handler != null)
                handler(this, e);
        }

        /// <summary>
        /// Checks if the instance has been disposed of, and if it has, throws an <see cref="T:System.ComponentModel.ObjectDisposedException"/>; otherwise, does nothing.
        /// </summary>
        /// <exception cref="T:System.ComponentModel.ObjectDisposedException">
        /// 	The instance has been disposed of.
        /// </exception>
        /// <remarks>
        /// 	Derived classes should call this method at the start of all methods and properties that should not be accessed after a call to <see cref="M:Dispose()"/>.
        /// </remarks>
        protected void CheckDisposed()
        {
            if (this._isDisposed)
                throw new ObjectDisposedException(this.GetType().FullName);
        }

        /// <summary>
        /// Releases all resources used by the instance.
        /// </summary>
        /// <remarks>
        /// 	Calls <see cref="M:Dispose(Boolean)"/> with the disposing parameter set to <see langword="true"/> to free unmanaged and managed resources.
        /// </remarks>
        public void Dispose()
        {
            if (!this._isDisposed)
            {
                this.Dispose(true);
                GC.SuppressFinalize(this);
            }
        }

        /// <summary>
        ///   Closes the <see cref="T:System.Data.IDataReader" /> Object.
        /// </summary>
        /// 
        public void Close()
        {
            this.Dispose();
        }

        /// <summary>
        /// Releases the unmanaged resources used by this instance and optionally releases the managed resources.
        /// </summary>
        /// <param name="disposing">
        /// 	<see langword="true"/> to release both managed and unmanaged resources; <see langword="false"/> to release only unmanaged resources.
        /// </param>
        protected virtual void Dispose(bool disposing)
        {
            // Refer to http://www.bluebytesoftware.com/blog/PermaLink,guid,88e62cdf-5919-4ac7-bc33-20c06ae539ae.aspx
            // Refer to http://www.gotdotnet.com/team/libraries/whitepapers/resourcemanagement/resourcemanagement.aspx

            // No exception should ever be thrown except in critical scenarios.
            // Unhandled exceptions during finalization will tear down the process.
            if (!this._isDisposed)
            {
                try
                {
                    // Dispose-time code should call Dispose() on all owned objects that implement the IDisposable interface. 
                    // "owned" means objects whose lifetime is solely controlled by the container. 
                    // In cases where ownership is not as straightforward, techniques such as HandleCollector can be used.  
                    // Large managed object fields should be nulled out.

                    // Dispose-time code should also set references of all owned objects to null, after disposing them. This will allow the referenced objects to be garbage collected even if not all references to the "parent" are released. It may be a significant memory consumption win if the referenced objects are large, such as big arrays, collections, etc. 
                    if (disposing)
                    {
                        // Acquire a lock on the object while disposing.

                        if (this._reader != null)
                        {
                            lock (this._lock)
                            {
                                if (this._reader != null)
                                {
                                    this._reader.Dispose();

                                    this._reader = null;
                                    this._buffer = null;
                                    this._eof = true;
                                }
                            }
                        }
                    }
                }
                finally
                {
                    // Ensure that the flag is set
                    this._isDisposed = true;

                    // Catch any issues about firing an event on an already disposed object.
                    try
                    {
                        this.OnDisposed(EventArgs.Empty);
                    }
                    catch
                    {
                    }
                }
            }
        }

#if !NETSTANDARD1_4
        /// <summary>
        /// Releases unmanaged resources and performs other cleanup operations before the instance is reclaimed by garbage collection.
        /// </summary>
        ~CsvReader()
        {
#if DEBUG
            Debug.WriteLine("FinalizableObject was not disposed" + this._allocStack.ToString());
#endif

            this.Dispose(false);
        }
#endif


    }
}