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

#if !NETSTANDARD1_4
namespace FileFormat.Accord.Statistics.Filters
{
    using System;
    using System.Collections.Generic;
    using System.Data;
    using System.Runtime.Serialization;
    using Base;
    using Math.Accord.Statistics;

    /// <summary>
    ///   Grouping filter.
    /// </summary>
    /// 
    [Serializable]
    public class Grouping : BaseFilter<Grouping.Options, Grouping>
    {

        [OptionalField]
        private bool lockGroups;

        [OptionalField]
        private int[] groupIndices;


        /// <summary>
        ///   Gets or sets a value indicating whether the group labels
        ///   are locked and should not be randomly re-selected.
        /// </summary>
        /// 
        /// <value><c>true</c> to lock groups; otherwise, <c>false</c>.</value>
        /// 
        public bool Lock
        {
            get { return this.lockGroups; }
            set { this.lockGroups = value; }
        }


        /// <summary>
        ///   Gets or sets the group index labels.
        /// </summary>
        /// 
        /// <value>The group indices.</value>
        /// 
        public int[] GroupIndices
        {
            get { return this.groupIndices; }
            set { this.groupIndices = value; }
        }

        /// <summary>
        ///   Gets or sets the two-group proportions.
        /// </summary>
        /// 
        public double Proportion { get; set; }

        /// <summary>
        ///   Gets or sets the name of the indicator
        ///   column which will be used to distinguish
        ///   samples from either group.
        /// </summary>
        /// 
        public string GroupIndicatorColumnName { get; set; }

        /// <summary>
        ///   Creates a new Grouping filter with equal group
        ///   proportions and default Group indicator column.
        /// </summary>
        /// 
        public Grouping()
        {
            this.Proportion = 0.5;
            this.GroupIndicatorColumnName = "Group";
        }

        /// <summary>
        ///   Creates a new Grouping filter.
        /// </summary>
        /// 
        public Grouping(string column)
        {
            this.Columns.Add(new Options(column));
        }


        /// <summary>
        ///   Processes the current filter.
        /// </summary>
        /// 
        protected override DataTable ProcessFilter(DataTable data)
        {

            if (!this.lockGroups)
            {
                // Check if we should balance label proportions
                if (this.Columns.Count == 0)
                {
                    // No. Just generate assign groups at random
                    this.groupIndices = Classes.Random(data.Rows.Count, this.Proportion);
                }

                else
                {
                    // Yes, we must balance the occurrences in a data column
                    this.groupIndices = this.balancedGroups(data);
                }
            }

            return this.apply(data);
        }

        private int[] balancedGroups(DataTable data)
        {
            // Works with only one column and for the binary case
            // TODO: Expand to multiple columns and multi-classes

            int[] classes = this.Columns[0].Classes;
            string column = this.Columns[0].ColumnName;
            int groupCount = 2;

            // Get subsets with 0 and 1
            List<DataRow>[] subsets = new List<DataRow>[classes.Length];
            for (int i = 0; i < subsets.Length; i++)
                subsets[i] = new List<DataRow>(data.Select("[" + column + "] = " + classes[i]));

            List<DataRow>[] groups = new List<DataRow>[groupCount];
            for (int i = 0; i < groups.Length; i++)
                groups[i] = new List<DataRow>();

            int totalPositives = subsets[0].Count;
            int totalNegatives = subsets[1].Count;


            int firstGroupPositives = (int)((subsets[0].Count / (double)groupCount) * 2 * this.Proportion);
            int firstGroupNegatives = (int)((subsets[1].Count / (double)groupCount) * 2 * this.Proportion);


            int[] groupIndices = new int[data.Rows.Count];

            // Put positives and negatives into first group
            for (int j = 0; j < firstGroupPositives; j++)
            {
                // Get first positive row
                DataRow row = subsets[0][0];
                subsets[0].Remove(row);
                groups[0].Add(row);
                groupIndices[row.Table.Rows.IndexOf(row)] = 0;
            }

            for (int j = 0; j < firstGroupNegatives; j++)
            {
                // Get first negative row
                DataRow row = subsets[1][0];
                subsets[1].Remove(row);
                groups[0].Add(row);
                groupIndices[row.Table.Rows.IndexOf(row)] = 0;
            }

            // Put positives and negatives into second group
            for (int j = 0; j < subsets[0].Count; j++)
            {
                // Get first positive row
                DataRow row = subsets[0][j];
                groups[1].Add(row);
                groupIndices[row.Table.Rows.IndexOf(row)] = 1;
            }

            for (int j = 0; j < subsets[1].Count; j++)
            {
                // Get first negative row
                DataRow row = subsets[1][j];
                groups[1].Add(row);
                groupIndices[row.Table.Rows.IndexOf(row)] = 1;
            }

            return groupIndices;
        }

        private DataTable apply(DataTable data)
        {
            DataTable result = data.Copy();
            if (!result.Columns.Contains(this.GroupIndicatorColumnName))
                result.Columns.Add(this.GroupIndicatorColumnName, typeof(int));

            for (int i = 0; i < result.Rows.Count; i++)
                result.Rows[i][this.GroupIndicatorColumnName] = this.groupIndices[i];

            return result;
        }

        /// <summary>
        ///   Options for the grouping filter.
        /// </summary>
        /// 
        [Serializable]
        public class Options : ColumnOptionsBase<Grouping>
        {
            /// <summary>
            ///   Gets or sets the labels used for each class contained in the column.
            /// </summary>
            /// 
            public int[] Classes { get; set; }

            /// <summary>
            ///   Constructs a new Options object for the given column.
            /// </summary>
            /// 
            /// <param name="name">
            ///   The name of the column to create this options for.
            /// </param>
            /// 
            public Options(String name)
                : base(name)
            {
                this.Classes = new int[] { 0, 1 };
            }

            /// <summary>
            ///   Constructs a new Options object.
            /// </summary>
            /// 
            public Options()
                : this("New column")
            {

            }
        }

    }
}
#endif