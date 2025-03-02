﻿// Accord Math Library
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
    using System.Collections;
    using System.Collections.Generic;
    using System.Linq;
    using System.Text;

    /// <summary>
    ///   Sparse vector representation (in LibSVM format).
    /// </summary>
    /// 
    /// <typeparam name="T">The type for the non-zero elements in this vector.</typeparam>
    /// 
    [Serializable]
    public sealed class Sparse<T> : IEnumerable<T>, ICloneable, IList<T>, IList, IFormattable
        where T : IEquatable<T>
    {
        private int[] indices;
        private T[] values;

        /// <summary>
        ///   Gets or sets the vector of indices indicating the location
        ///   of the non-zero elements contained in this sparse vector.
        /// </summary>
        /// 
        public int[] Indices
        {
            get { return this.indices; }
            set { this.indices = value; }
        }

        /// <summary>
        ///   Gets or sets the vector of values indicating which non-zero
        ///   value happens at each position indicated in <see cref="Indices"/>.
        /// </summary>
        /// 
        public T[] Values
        {
            get { return this.values; }
            set { this.values = value; }
        }

        /// <summary>
        ///   Creates a sparse vector with zero elements.
        /// </summary>
        /// 
        public Sparse()
            : this(0)
        {
        }

        /// <summary>
        ///   Creates a sparse vector with the maximum number of elements.
        /// </summary>
        /// 
        /// <param name="length">The maximum number of non-zero
        ///   elements that this vector can accomodate.</param>
        ///   
        public Sparse(int length)
        {
            this.indices = new int[length];
            this.values = new T[length];
        }

        /// <summary>
        ///   Creates a sparse vector from a vector of indices
        ///   and a vector of values occuring at those indices.
        /// </summary>
        /// 
        /// <param name="indices">The indices for non-zero entries.</param>
        /// <param name="values">The non-zero values happening at each index.</param>
        /// 
        public Sparse(int[] indices, T[] values)
        {
            this.indices = indices;
            this.values = values;
        }

        /// <summary>
        ///   Converts this sparse vector to a dense vector of the given length.
        /// </summary>
        /// 
        public T[] ToDense()
        {
            return this.ToDense(this.Indices.Max() + 1);
        }

        /// <summary>
        ///   Converts this sparse vector to a dense vector of the given length.
        /// </summary>
        /// 
        public T[] ToDense(int length)
        {
            T[] result = new T[length];
            for (int i = 0; i < this.Indices.Length; i++)
                result[this.Indices[i]] = this.Values[i];
            return result;
        }

        /// <summary>
        ///   Converts this sparse vector to a sparse representation where 
        ///   the indices are intertwined with their corresponding values.
        /// </summary>
        /// 
        public T[] ToSparse()
        {
            T[] result = new T[this.Indices.Length * 2];
            for (int i = 0; i < this.Indices.Length; i++)
            {
                result[2 * i + 0] = (T)System.Convert.ChangeType(this.Indices[i], typeof(T));
                result[2 * i + 1] = this.Values[i];
            }

            return result;
        }


        /// <summary>
        ///   Gets the the value stored at position <paramref name="i"/>.
        /// </summary>
        /// 
        public T this[int i]
        {
            get
            {
                int j = Array.IndexOf(this.Indices, i);
                if (j >= 0)
                    return this.Values[j];
                return default(T);
            }
            set
            {
                int j = Array.IndexOf(this.Indices, i);
                if (j >= 0)
                {
                    this.Values[j] = value;
                    return;
                }

                T[] newValues = new T[this.Values.Length + 1];
                int[] newIndices = new int[this.indices.Length + 1];

                int k;
                for (k = 0; k < this.indices.Length && this.indices[k] < i; k++)
                {
                    newIndices[k] = this.indices[k];
                    newValues[k] = this.values[k];
                }
                newIndices[k] = i;
                newValues[k] = value;
                k++;
                for (; k < newIndices.Length; k++)
                {
                    newIndices[k] = this.indices[k - 1];
                    newValues[k] = this.values[k - 1];
                }

                this.indices = newIndices;
                this.values = newValues;
            }
        }

        /// <summary>
        ///   Creates a new object that is a copy of the current instance.
        /// </summary>
        /// 
        /// <returns>
        ///   A new object that is a copy of this instance.
        /// </returns>
        /// 
        public object Clone()
        {
            return new Sparse<T>((int[])this.indices.Clone(), (T[])this.values.Clone());
        }

        /// <summary>
        ///   Performs an implicit conversion from <see cref="Sparse{T}"/> to <see cref="Array"/>.
        /// </summary>
        /// 
        public static implicit operator Array(Sparse<T> obj)
        {
            return obj.Values;
        }

        /// <summary>
        ///   Gets the maximum non-zero element index in the sparse vector.
        /// </summary>
        /// 
        public int Length
        {
            get
            {
                T zero = default(T);

                for (int i = this.Indices.Length - 1; i >= 0; i--)
                {
                    if (!this.Values[i].Equals(zero))
                        return this.Indices[i] + 1;
                }

                return 0;
            }
        }





        int IList<T>.IndexOf(T item)
        {
            throw new NotImplementedException();
        }

        void IList<T>.Insert(int index, T item)
        {
            throw new NotImplementedException();
        }

        void IList<T>.RemoveAt(int index)
        {
            throw new NotImplementedException();
        }

        T IList<T>.this[int index]
        {
            get { return this[index]; }
            set { throw new NotImplementedException(); }
        }

        void ICollection<T>.Add(T item)
        {
            throw new NotImplementedException();
        }

        void ICollection<T>.Clear()
        {
            Array.Clear(this.values, 0, this.values.Length);
        }

        bool ICollection<T>.Contains(T item)
        {
            return this.values.Contains(item);
        }

        /// <summary>
        /// Copies the elements of the <see cref="T:System.Collections.Generic.ICollection`1" /> to an <see cref="T:System.Array" />, starting at a particular <see cref="T:System.Array" /> index.
        /// </summary>
        /// <param name="array">The one-dimensional <see cref="T:System.Array" /> that is the destination of the elements copied from <see cref="T:System.Collections.Generic.ICollection`1" />. The <see cref="T:System.Array" /> must have zero-based indexing.</param>
        /// <param name="arrayIndex">The zero-based index in <paramref name="array" /> at which copying begins.</param>
        public void CopyTo(T[] array, int arrayIndex)
        {
            for (int i = 0; i < this.indices.Length; i++)
                array[arrayIndex + this.indices[i]] = this.values[i];
        }

        int ICollection<T>.Count
        {
            get { return 0; }
        }

        bool ICollection<T>.IsReadOnly
        {
            get { return false; }
        }

        bool ICollection<T>.Remove(T item)
        {
            throw new NotImplementedException();
        }

        IEnumerator<T> IEnumerable<T>.GetEnumerator()
        {
            foreach (T t in this.ToDense())
                yield return t;
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            foreach (T t in this.ToDense())
                yield return t;
        }

        int IList.Add(object value)
        {
            throw new NotImplementedException();
        }

        void IList.Clear()
        {
            Array.Clear(this.values, 0, this.values.Length);
        }

        bool IList.Contains(object value)
        {
            return ((IList)this.Values).Contains(value);
        }

        int IList.IndexOf(object value)
        {
            throw new NotImplementedException();
        }

        void IList.Insert(int index, object value)
        {
            throw new NotImplementedException();
        }

        bool IList.IsFixedSize
        {
            get { return true; }
        }

        bool IList.IsReadOnly
        {
            get { return false; }
        }

        void IList.Remove(object value)
        {
            throw new NotImplementedException();
        }

        void IList.RemoveAt(int index)
        {
            throw new NotImplementedException();
        }

        object IList.this[int index]
        {
            get { return this[index]; }
            set { throw new NotImplementedException(); }
        }

        void ICollection.CopyTo(Array array, int index)
        {
            for (int i = 0; i < this.indices.Length; i++)
                array.SetValue(this.values[i], index + this.indices[i]);
        }

        int ICollection.Count
        {
            get { return 0; }
        }

        bool ICollection.IsSynchronized
        {
            get { throw new NotImplementedException(); }
        }

        object ICollection.SyncRoot
        {
            get { throw new NotImplementedException(); }
        }

        /// <summary>
        /// Returns a <see cref="System.String" /> that represents this instance.
        /// </summary>
        /// <returns>
        /// A <see cref="System.String" /> that represents this instance.
        /// </returns>
        public override string ToString()
        {
            return this.ToString("g", System.Globalization.CultureInfo.CurrentCulture);
        }

        /// <summary>
        /// Returns a <see cref="System.String" /> that represents this instance.
        /// </summary>
        /// <param name="format">The format.</param>
        /// <param name="formatProvider">The format provider.</param>
        /// <returns>
        /// A <see cref="System.String" /> that represents this instance.
        /// </returns>
        public string ToString(string format, IFormatProvider formatProvider)
        {
            var sb = new StringBuilder();
            for (int i = 0; i < this.Indices.Length; i++)
            {
                sb.Append(this.Indices[i] + 1); // Note: LibSVM array format is one-based
                sb.Append(":");
                sb.AppendFormat(formatProvider, "{0:" + format + "}", this.Values[i]);
                if (i < this.Indices.Length - 1)
                    sb.Append(" ");
            }
            return sb.ToString();
        }

        /// <summary>
        /// Determines whether this Sparse vector has elements on all indices.
        /// </summary>
        /// 
        /// <returns><c>true</c> if this instance is full; otherwise, <c>false</c>.</returns>
        /// 
        public bool IsFull()
        {
            for (int i = 0; i < this.Indices.Length; i++)
            {
                if (this.Indices[i] != i)
                    return false;
            }

            return true;
        }
    }
}
