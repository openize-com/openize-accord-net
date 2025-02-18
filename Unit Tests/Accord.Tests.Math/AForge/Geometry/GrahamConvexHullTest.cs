namespace Openize.Accord.Tests.Math.AForge.Geometry
{
    using System.Collections.Generic;
    using NUnit.Framework;
    using Openize.Accord.Core.AForge.Core;
    using Openize.Accord.Math.AForge.Math.Geometry;

    [TestFixture]
    public class GrahamConvexHullTest
    {
        private List<IntPoint> pointsList0 = new List<IntPoint>( );
        private List<IntPoint> pointsList1 = new List<IntPoint>( );
        private List<IntPoint> pointsList2 = new List<IntPoint>( );
        private List<IntPoint> pointsList3 = new List<IntPoint>( );
        private List<IntPoint> pointsList4 = new List<IntPoint>( );
        private List<IntPoint> pointsList5 = new List<IntPoint>( );
        private List<IntPoint> pointsList6 = new List<IntPoint>( );

        private List<IntPoint> pointsList7 = new List<IntPoint>( );
        private List<IntPoint> pointsList8 = new List<IntPoint>( );
        private List<IntPoint> pointsList9 = new List<IntPoint>( );

        private List<IntPoint> expectedHull8 = new List<IntPoint>( );

        private List<List<IntPoint>> pointsLists = new List<List<IntPoint>>( );
        private List<List<IntPoint>> expectedHulls = new List<List<IntPoint>>( );

        public GrahamConvexHullTest( )
        {
            // prepare 0st list
            this.pointsList0.Add( new IntPoint( 0, 0 ) );

            // prepare 1st list
            this.pointsList1.Add( new IntPoint( 0, 0 ) );
            this.pointsList1.Add( new IntPoint( 100, 0 ) );

            // prepare 2nd list
            this.pointsList2.AddRange( this.pointsList1 );
            this.pointsList2.Add( new IntPoint( 100, 100 ) );

            // prepare 3rd list
            this.pointsList3.AddRange( this.pointsList2 );
            this.pointsList3.Add( new IntPoint( 0, 100 ) );

            // prepare 4th list
            this.pointsList4.AddRange( this.pointsList2 );
            this.pointsList4.Add( new IntPoint( 60, 40 ) );

            // prepare 5th list
            this.pointsList5.AddRange( this.pointsList3 );
            this.pointsList5.Add( new IntPoint( 50, 50 ) );

            // prepare 6th list
            this.pointsList6.AddRange( this.pointsList3 );
            this.pointsList6.Add( new IntPoint( 0, 0 ) );

            // prepare 7th list
            this.pointsList7.AddRange( this.pointsList3 );
            this.pointsList7.AddRange( this.pointsList3 );

            // prepare 8th list
            this.pointsList8.AddRange( this.pointsList3 );
            this.pointsList8.Add( new IntPoint( 50, -10 ) );
            this.pointsList8.Add( new IntPoint( 110, 50 ) );
            this.pointsList8.Add( new IntPoint( 50, 110 ) );

            this.expectedHull8.AddRange( this.pointsList3 );
            this.expectedHull8.Insert( 1, new IntPoint( 50, -10 ) );
            this.expectedHull8.Insert( 3, new IntPoint( 110, 50 ) );
            this.expectedHull8.Insert( 5, new IntPoint( 50, 110 ) );

            // prepare 9th list
            this.pointsList9.AddRange( this.pointsList8 );
            this.pointsList9.Add( new IntPoint( 50, 10 ) );
            this.pointsList9.Add( new IntPoint( 90, 50 ) );
            this.pointsList9.Add( new IntPoint( 50, 90 ) );
            this.pointsList9.Add( new IntPoint( 10, 50 ) );

            // now prepare list of tests
            this.pointsLists.Add( this.pointsList0 );
            this.pointsLists.Add( this.pointsList1 );
            this.pointsLists.Add( this.pointsList2 );
            this.pointsLists.Add( this.pointsList3 );

            this.expectedHulls.AddRange( this.pointsLists );

            this.pointsLists.Add( this.pointsList4 );
            this.expectedHulls.Add( this.pointsList2 );

            this.pointsLists.Add( this.pointsList5 );
            this.expectedHulls.Add( this.pointsList3 );

            this.pointsLists.Add( this.pointsList6 );
            this.expectedHulls.Add( this.pointsList3 );

            this.pointsLists.Add( this.pointsList7 );
            this.expectedHulls.Add( this.pointsList3 );

            this.pointsLists.Add( this.pointsList8 );
            this.expectedHulls.Add( this.expectedHull8 );

            this.pointsLists.Add( this.pointsList9 );
            this.expectedHulls.Add( this.expectedHull8 );
        }

        [Test]        
        public void FindHullTest( )
        {
            GrahamConvexHull grahamHull = new GrahamConvexHull( );

            for ( int i = 0, n = this.pointsLists.Count; i < n; i++ )
            {
                this.ComparePointsLists( grahamHull.FindHull( this.pointsLists[i] ), this.expectedHulls[i] );
            }
        }

        private void ComparePointsLists( List<IntPoint> list1, List<IntPoint> list2 )
        {
            Assert.AreEqual( list1.Count, list2.Count );

            if ( list1.Count == list2.Count )
            {
                for ( int i = 0, n = list1.Count; i < n; i++ )
                {
                    Assert.AreEqual( list2[i], list1[i] );
                }
            }
        }
    }
}
