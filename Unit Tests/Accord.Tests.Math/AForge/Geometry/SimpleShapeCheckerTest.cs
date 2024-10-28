namespace FileFormat.Accord.Tests.Math.AForge.Geometry
{
    using System.Collections.Generic;
    using FileFormat.Accord.Core.AForge.Core;
    using FileFormat.Accord.Math.AForge.Math.Geometry;
    using NUnit.Framework;

    [TestFixture]
    public class SimpleShapeCheckerTest
    {
        private SimpleShapeChecker shapeChecker = new SimpleShapeChecker( );

        private List<IntPoint> idealCicle = new List<IntPoint>( );
        private List<IntPoint> distorredCircle = new List<IntPoint>( );

        private List<IntPoint> square1 = new List<IntPoint>( );
        private List<IntPoint> square1Test = new List<IntPoint>( );
        private List<IntPoint> square2 = new List<IntPoint>( );
        private List<IntPoint> square2Test = new List<IntPoint>( );
        private List<IntPoint> square3 = new List<IntPoint>( );
        private List<IntPoint> rectangle = new List<IntPoint>( );

        private List<IntPoint> triangle1 = new List<IntPoint>( );
        private List<IntPoint> isoscelesTriangle = new List<IntPoint>( );
        private List<IntPoint> equilateralTriangle = new List<IntPoint>( );
        private List<IntPoint> rectangledTriangle = new List<IntPoint>( );

        public SimpleShapeCheckerTest( )
        {
            System.Random rand = new System.Random( );

            // generate sample circles
            double radius = 100;

            for ( int i = 0; i < 360; i += 10 )
            {
                double angle = (double) i / 180 * System.Math.PI;

                // add point to ideal circle
                this.idealCicle.Add( new IntPoint(
                    (int) ( radius * System.Math.Cos( angle ) ),
                    (int) ( radius * System.Math.Sin( angle ) ) ) );

                // add a bit distortion for distorred cirlce
                double distorredRadius = radius + rand.Next( 7 ) - 3;

                this.distorredCircle.Add( new IntPoint(
                    (int) ( distorredRadius * System.Math.Cos( angle ) ),
                    (int) ( distorredRadius * System.Math.Sin( angle ) ) ) );
            }

            // generate sample squares
            this.square1.Add( new IntPoint( 0, 0 ) );
            this.square1.Add( new IntPoint( 50, 0 ) );
            this.square1.Add( new IntPoint( 100, 0 ) );
            this.square1.Add( new IntPoint( 100, 50 ) );
            this.square1.Add( new IntPoint( 100, 100 ) );
            this.square1.Add( new IntPoint( 50, 100 ) );
            this.square1.Add( new IntPoint( 0, 100 ) );
            this.square1.Add( new IntPoint( 0, 50 ) );

            this.square2.Add( new IntPoint( 50, 0 ) );
            this.square2.Add( new IntPoint( 75, 25 ) );
            this.square2.Add( new IntPoint( 100, 50 ) );
            this.square2.Add( new IntPoint( 75, 75 ) );
            this.square2.Add( new IntPoint( 50, 100 ) );
            this.square2.Add( new IntPoint( 25, 75 ) );
            this.square2.Add( new IntPoint( 0, 50 ) );
            this.square2.Add( new IntPoint( 25, 25 ) );

            // these should be obtained as corners
            this.square1Test.Add( new IntPoint( 0, 0 ) );
            this.square1Test.Add( new IntPoint( 100, 0 ) );
            this.square1Test.Add( new IntPoint( 100, 100 ) );
            this.square1Test.Add( new IntPoint( 0, 100 ) );

            this.square2Test.Add( new IntPoint( 50, 0 ) );
            this.square2Test.Add( new IntPoint( 100, 50 ) );
            this.square2Test.Add( new IntPoint( 50, 100 ) );
            this.square2Test.Add( new IntPoint( 0, 50 ) );

            // special square, which may look like circle, but should be recognized as circle
            this.square3.Add( new IntPoint( 50, 0 ) );
            this.square3.Add( new IntPoint( 100, 50 ) );
            this.square3.Add( new IntPoint( 50, 100 ) );
            this.square3.Add( new IntPoint( 0, 50 ) );

            // generate sample rectangle
            this.rectangle.Add( new IntPoint( 0, 0 ) );
            this.rectangle.Add( new IntPoint( 50, 0 ) );
            this.rectangle.Add( new IntPoint( 100, 0 ) );
            this.rectangle.Add( new IntPoint( 100, 20 ) );
            this.rectangle.Add( new IntPoint( 100, 40 ) );
            this.rectangle.Add( new IntPoint( 50, 40 ) );
            this.rectangle.Add( new IntPoint( 0, 40 ) );
            this.rectangle.Add( new IntPoint( 0, 20 ) );

            // generate some triangles
            this.triangle1.Add( new IntPoint( 0, 0 ) );
            this.triangle1.Add( new IntPoint( 50, 10 ) );
            this.triangle1.Add( new IntPoint( 100, 20 ) );
            this.triangle1.Add( new IntPoint( 90, 50 ) );
            this.triangle1.Add( new IntPoint( 80, 80 ) );
            this.triangle1.Add( new IntPoint( 40, 40 ) );

            this.isoscelesTriangle.Add( new IntPoint( 0, 0 ) );
            this.isoscelesTriangle.Add( new IntPoint( 50, 0 ) );
            this.isoscelesTriangle.Add( new IntPoint( 100, 0 ) );
            this.isoscelesTriangle.Add( new IntPoint( 75, 20 ) );
            this.isoscelesTriangle.Add( new IntPoint( 50, 40 ) );
            this.isoscelesTriangle.Add( new IntPoint( 25, 20 ) );

            this.equilateralTriangle.Add( new IntPoint( 0, 0 ) );
            this.equilateralTriangle.Add( new IntPoint( 50, 0 ) );
            this.equilateralTriangle.Add( new IntPoint( 100, 0 ) );
            this.equilateralTriangle.Add( new IntPoint( 75, 43 ) );
            this.equilateralTriangle.Add( new IntPoint( 50, 86 ) );
            this.equilateralTriangle.Add( new IntPoint( 25, 43 ) );

            this.rectangledTriangle.Add( new IntPoint( 0, 0 ) );
            this.rectangledTriangle.Add( new IntPoint( 20, 0 ) );
            this.rectangledTriangle.Add( new IntPoint( 40, 0 ) );
            this.rectangledTriangle.Add( new IntPoint( 20, 50 ) );
            this.rectangledTriangle.Add( new IntPoint( 0, 100 ) );
            this.rectangledTriangle.Add( new IntPoint( 0, 50 ) );
        }

        [Test]
        public void IsCircleTest( )
        {
            Assert.AreEqual( true, this.shapeChecker.IsCircle( this.idealCicle ) );
            Assert.AreEqual( true, this.shapeChecker.IsCircle( this.distorredCircle ) );

            Assert.AreEqual( false, this.shapeChecker.IsCircle( this.square1 ) );
            Assert.AreEqual( false, this.shapeChecker.IsCircle( this.square2 ) );
            Assert.AreEqual( false, this.shapeChecker.IsCircle( this.square3 ) );
            Assert.AreEqual( false, this.shapeChecker.IsCircle( this.rectangle ) );

            Assert.AreEqual( false, this.shapeChecker.IsCircle( this.triangle1 ) );
            Assert.AreEqual( false, this.shapeChecker.IsCircle( this.equilateralTriangle ) );
            Assert.AreEqual( false, this.shapeChecker.IsCircle( this.isoscelesTriangle ) );
            Assert.AreEqual( false, this.shapeChecker.IsCircle( this.rectangledTriangle ) );
        }

        [Test]
        public void IsQuadrilateralTest( )
        {
            Assert.AreEqual( true, this.shapeChecker.IsQuadrilateral( this.square1 ) );
            Assert.AreEqual( true, this.shapeChecker.IsQuadrilateral( this.square2 ) );
            Assert.AreEqual( true, this.shapeChecker.IsQuadrilateral( this.square3 ) );
            Assert.AreEqual( true, this.shapeChecker.IsQuadrilateral( this.rectangle ) );

            Assert.AreEqual( false, this.shapeChecker.IsQuadrilateral( this.idealCicle ) );
            Assert.AreEqual( false, this.shapeChecker.IsQuadrilateral( this.distorredCircle ) );

            Assert.AreEqual( false, this.shapeChecker.IsQuadrilateral( this.triangle1 ) );
            Assert.AreEqual( false, this.shapeChecker.IsQuadrilateral( this.equilateralTriangle ) );
            Assert.AreEqual( false, this.shapeChecker.IsQuadrilateral( this.isoscelesTriangle ) );
            Assert.AreEqual( false, this.shapeChecker.IsQuadrilateral( this.rectangledTriangle ) );
        }

        [Test]
        public void CheckQuadrilateralCornersTest( )
        {
            List<IntPoint> corners;

            Assert.AreEqual( true, this.shapeChecker.IsQuadrilateral( this.square1, out corners ) );
            Assert.AreEqual( 4, corners.Count );
            Assert.AreEqual( true, this.CompareShape( corners, this.square1Test ) );

            Assert.AreEqual( true, this.shapeChecker.IsQuadrilateral( this.square2, out corners ) );
            Assert.AreEqual( 4, corners.Count );
            Assert.AreEqual( true, this.CompareShape( corners, this.square2Test ) );
        }

        [Test]
        public void IsTriangleTest( )
        {
            Assert.AreEqual( true, this.shapeChecker.IsTriangle( this.triangle1 ) );
            Assert.AreEqual( true, this.shapeChecker.IsTriangle( this.equilateralTriangle ) );
            Assert.AreEqual( true, this.shapeChecker.IsTriangle( this.isoscelesTriangle ) );
            Assert.AreEqual( true, this.shapeChecker.IsTriangle( this.rectangledTriangle ) );

            Assert.AreEqual( false, this.shapeChecker.IsTriangle( this.idealCicle ) );
            Assert.AreEqual( false, this.shapeChecker.IsTriangle( this.distorredCircle ) );

            Assert.AreEqual( false, this.shapeChecker.IsTriangle( this.square1 ) );
            Assert.AreEqual( false, this.shapeChecker.IsTriangle( this.square2 ) );
            Assert.AreEqual( false, this.shapeChecker.IsTriangle( this.square3 ) );
            Assert.AreEqual( false, this.shapeChecker.IsTriangle( this.rectangle ) );
        }

        [Test]
        public void IsConvexPolygon( )
        {
            List<IntPoint> corners;

            Assert.AreEqual( true, this.shapeChecker.IsConvexPolygon( this.triangle1, out corners ) );
            Assert.AreEqual( 3, corners.Count );
            Assert.AreEqual( true, this.shapeChecker.IsConvexPolygon( this.equilateralTriangle, out corners ) );
            Assert.AreEqual( 3, corners.Count );
            Assert.AreEqual( true, this.shapeChecker.IsConvexPolygon( this.isoscelesTriangle, out corners ) );
            Assert.AreEqual( 3, corners.Count );
            Assert.AreEqual( true, this.shapeChecker.IsConvexPolygon( this.rectangledTriangle, out corners ) );
            Assert.AreEqual( 3, corners.Count );

            Assert.AreEqual( true, this.shapeChecker.IsConvexPolygon( this.square1, out corners ) );
            Assert.AreEqual( 4, corners.Count );
            Assert.AreEqual( true, this.shapeChecker.IsConvexPolygon( this.square2, out corners ) );
            Assert.AreEqual( 4, corners.Count );
            Assert.AreEqual( true, this.shapeChecker.IsConvexPolygon( this.square3, out corners ) );
            Assert.AreEqual( 4, corners.Count );
            Assert.AreEqual( true, this.shapeChecker.IsConvexPolygon( this.rectangle, out corners ) );
            Assert.AreEqual( 4, corners.Count );

            Assert.AreEqual( false, this.shapeChecker.IsConvexPolygon( this.idealCicle, out corners ) );
            Assert.AreEqual( false, this.shapeChecker.IsConvexPolygon( this.distorredCircle, out corners ) );
        }

        [Test]
        public void CheckShapeTypeTest( )
        {
            Assert.AreEqual( ShapeType.Circle, this.shapeChecker.CheckShapeType( this.idealCicle ) );
            Assert.AreEqual( ShapeType.Circle, this.shapeChecker.CheckShapeType( this.distorredCircle ) );

            Assert.AreEqual( ShapeType.Quadrilateral, this.shapeChecker.CheckShapeType( this.square1 ) );
            Assert.AreEqual( ShapeType.Quadrilateral, this.shapeChecker.CheckShapeType( this.square2 ) );
            Assert.AreEqual( ShapeType.Quadrilateral, this.shapeChecker.CheckShapeType( this.square3 ) );
            Assert.AreEqual( ShapeType.Quadrilateral, this.shapeChecker.CheckShapeType( this.rectangle ) );

            Assert.AreEqual( ShapeType.Triangle, this.shapeChecker.CheckShapeType( this.triangle1 ) );
            Assert.AreEqual( ShapeType.Triangle, this.shapeChecker.CheckShapeType( this.equilateralTriangle ) );
            Assert.AreEqual( ShapeType.Triangle, this.shapeChecker.CheckShapeType( this.isoscelesTriangle ) );
            Assert.AreEqual( ShapeType.Triangle, this.shapeChecker.CheckShapeType( this.rectangledTriangle ) );
        }

        private bool CompareShape( List<IntPoint> shape1, List<IntPoint> shape2 )
        {
            if ( shape1.Count != shape2.Count )
                return false;
            if ( shape1.Count == 0 )
                return true;

            int index = shape1.IndexOf( shape2[0] );

            if ( index == -1 )
                return false;

            index++;

            for ( int i = 1; i < shape2.Count; i++, index++ )
            {
                if ( index >= shape1.Count )
                    index = 0;

                if ( !shape1[index].Equals( shape2[i] ) )
                    return false;
            }

            return true;
        }

        
        [TestCase( PolygonSubType.Unknown, new int[] { 0, 0, 100, 0, 90, 10 } )]     // just a triangle
        [TestCase( PolygonSubType.IsoscelesTriangle, new int[] { 0, 0, 100, 0, 50, 10 } )]
        [TestCase( PolygonSubType.IsoscelesTriangle, new int[] { 0, 0, 100, 0, 50, 200 } )]
        [TestCase( PolygonSubType.EquilateralTriangle, new int[] { 0, 0, 100, 0, 50, 86 } )]
        [TestCase( PolygonSubType.RectangledIsoscelesTriangle, new int[] { 0, 0, 100, 0, 50, 50 } )]
        [TestCase( PolygonSubType.RectangledIsoscelesTriangle, new int[] { 0, 0, 100, 0, 0, 100 } )]
        [TestCase( PolygonSubType.RectangledTriangle, new int[] { 0, 0, 100, 0, 0, 50 } )]
        [TestCase( PolygonSubType.Unknown, new int[] { 0, 0, 100, 0, 90, 50, 10, 70 } )]     // just a quadrilateral
        [TestCase( PolygonSubType.Trapezoid, new int[] { 0, 0, 100, 0, 90, 50, 10, 50 } )]
        [TestCase( PolygonSubType.Trapezoid, new int[] { 0, 0, 100, 0, 90, 50, 0, 50 } )]
        [TestCase( PolygonSubType.Trapezoid, new int[] { 0, 0, 100, 0, 90, 50, 0, 53 } )]    // a bit disformed
        [TestCase( PolygonSubType.Parallelogram, new int[] { 0, 0, 100, 0, 120, 50, 20, 50 } )]
        [TestCase( PolygonSubType.Parallelogram, new int[] { 0, 0, 100, 0, 70, 50, -30, 50 } )]
        [TestCase( PolygonSubType.Rectangle, new int[] { 0, 0, 100, 0, 100, 50, 0, 50 } )]
        [TestCase( PolygonSubType.Rectangle, new int[] { 0, 0, 100, 0, 100, 52, -3, 50 } )]   // a bit disformed
        [TestCase( PolygonSubType.Square, new int[] { 0, 0, 100, 0, 100, 100, 0, 100 } )]
        [TestCase( PolygonSubType.Square, new int[] { 50, 0, 100, 50, 50, 100, 0, 50 } )]
        [TestCase( PolygonSubType.Square, new int[] { 51, 0, 100, 49, 50, 101, 1, 50 } )]    // a bit disformed
        [TestCase( PolygonSubType.Rhombus, new int[] { 30, 0, 60, 50, 30, 100, 0, 50 } )]
        [TestCase( PolygonSubType.Rhombus, new int[] { 0, 0, 100, 0, 130, 95, 30, 95 } )]
        [TestCase( PolygonSubType.Unknown, new int[] { 0, 0, 100, 0, 90, 50, 40, 70, 10, 40 } )]     // unknown if 5 corners or more
        public void CheckPolygonSubTypeTest( PolygonSubType expectedSubType, int[] corners )
        {
            Assert.AreEqual( expectedSubType, this.shapeChecker.CheckPolygonSubType( this.GetListOfPointFromArray( corners ) ) );
        }

        private List<IntPoint> GetListOfPointFromArray( int[] points )
        {
            List<IntPoint> list = new List<IntPoint>( );

            for ( int i = 0, n = points.Length; i < n; i += 2 )
            {
                list.Add( new IntPoint( points[i], points[i + 1] ) );
            }

            return list;
        }
    }
}
