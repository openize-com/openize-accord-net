﻿namespace Openize.Accord.Tests.Math.AForge.Geometry
{
    using NUnit.Framework;
    using Openize.Accord.Math.AForge.Math.Geometry;

    [TestFixture]
    public class LineStraighteningOptimizerTest
    {
        private IShapeOptimizer optimizer = new LineStraighteningOptimizer( 3 );

        
        [TestCase( new int[] { 0, 0, 10, 0, 10, 10 }, new int[] { 0, 0, 10, 0, 10, 10 } )]
        [TestCase( new int[] { 0, 0, 10, 0, 5, 1 }, new int[] { 0, 0, 10, 0, 5, 1 } )]
        [TestCase( new int[] { 0, 0, 10, 0, 10, 10, 5, 5 }, new int[] { 0, 0, 10, 0, 10, 10 } )]
        [TestCase( new int[] { 0, 0, 10, 0, 10, 10, 4, 6 }, new int[] { 0, 0, 10, 0, 10, 10 } )]
        [TestCase( new int[] { 0, 0, 10, 0, 10, 10, 6, 4 }, new int[] { 0, 0, 10, 0, 10, 10 } )]
        [TestCase( new int[] { 0, 0, 10, 0, 10, 10, 7, 3 }, new int[] { 0, 0, 10, 0, 10, 10 } )]
        [TestCase( new int[] { 0, 0, 10, 0, 10, 10, 8, 2 }, new int[] { 0, 0, 10, 0, 10, 10, 8, 2 } )]
        [TestCase( new int[] { 4, 6, 0, 0, 5, 1, 10, 0, 10, 5, 10, 10 }, new int[] { 0, 0, 10, 0, 10, 10 } )]
        [TestCase( new int[] { 6, 4, 0, 0, 6, -1, 10, 0, 9, 4, 10, 10 }, new int[] { 0, 0, 10, 0, 10, 10 } )]
        public void OptimizationTest( int[] coordinates, int[] expectedCoordinates )
        {
            ShapeOptimizerTestBase.TestOptimizer( coordinates, expectedCoordinates, this.optimizer );
        }
    }
}
