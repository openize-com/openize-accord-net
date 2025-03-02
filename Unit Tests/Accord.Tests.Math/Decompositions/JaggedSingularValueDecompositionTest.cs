﻿// Accord Unit Tests
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

namespace Openize.Accord.Tests.Math.Decompositions
{
    using Openize.Accord.Math.Matrix;
    using NUnit.Framework;
    using Openize.Accord.Math.Decompositions;

    [TestFixture]
    public class JaggedSingularValueDecompositionTest
    {

        [Test]
        public void InverseTestNaN()
        {
            int n = 5;

            var I = Matrix.Identity(n).ToJagged();

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    double[][] value = Matrix.Magic(n).ToJagged();

                    value[i][j] = double.NaN;

                    var target = new JaggedSingularValueDecomposition(value);

                    double[][] solution = target.Solve(I);
                    double[][] inverse = target.Inverse();

                    Assert.IsTrue(Matrix.IsEqual(solution, inverse));
                }
            }
        }

        [Test]
        public void InverseTest()
        {
            var value = new double[][]
            { 
                  new double[] { 1.0, 1.0 },
                  new double[] { 2.0, 2.0 }
            };

            var target = new JaggedSingularValueDecomposition(value);

            var expected = new double[][]
            {
               new double[] { 0.1, 0.2 },
               new double[] { 0.1, 0.2 }
            };

            var actual = target.Solve(Matrix.JaggedIdentity(2));
            Assert.IsTrue(Matrix.IsEqual(expected, actual, 1e-3));
            Assert.IsTrue(Matrix.IsEqual(value, target.Reverse(), 1e-5));
            actual = target.Inverse();
            Assert.IsTrue(Matrix.IsEqual(expected, actual, 1e-3));
        }

        [Test]
        public void InverseTest2()
        {
            int n = 5;

            var I = Jagged.Identity(n);

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    double[][] value = Jagged.Magic(n);

                    var target = new JaggedSingularValueDecomposition(value);

                    double[][] solution = target.Solve(I);
                    double[][] inverse = target.Inverse();
                    double[][] reverse = target.Reverse();

                    Assert.IsTrue(Matrix.IsEqual(solution, inverse, 1e-4));
                    Assert.IsTrue(Matrix.IsEqual(value, reverse, 1e-4));
                }
            }
        }

        [Test]
        public void JaggedSingularValueDecompositionConstructorTest1()
        {
            // This test catches the bug in SingularValueDecomposition in the line
            //   for (int j = k + 1; j < nu; j++)
            // where it should be
            //   for (int j = k + 1; j < n; j++)


            // Test for m-x-n matrices where m < n. The available SVD
            // routine was not meant to be used in this case.

            var value = new double[][]
            { 
                 new double[] { 1, 2 },
                 new double[] { 3, 4 },
                 new double[] { 5, 6 },
                 new double[] { 7, 8 }
            }.Transpose(); // value is 2x4, having less rows than columns.

            var target = new JaggedSingularValueDecomposition(value, true, true, false);

            double[][] actual = Matrix.Multiply(Matrix.Multiply(
                target.LeftSingularVectors, target.DiagonalMatrix),
                target.RightSingularVectors.Transpose());

            // Checking the decomposition
            Assert.IsTrue(Matrix.IsEqual(actual, value, 1e-2));
            Assert.IsTrue(Matrix.IsEqual(value, target.Reverse(), 1e-2));

            // Checking values
            var U = new double[][]
            {
                new double[] { -0.641423027995072, -0.767187395072177 },
                new double[] { -0.767187395072177,  0.641423027995072 },
            };

            // U should be equal
            Assert.IsTrue(Matrix.IsEqual(target.LeftSingularVectors, U, 0.001));


            double[][] V = new double[][]// economy svd
            {
                new double[] { -0.152483233310201,  0.822647472225661,  },
                new double[] { -0.349918371807964,  0.421375287684580,  },
                new double[] { -0.547353510305727,  0.0201031031435023, },
                new double[] { -0.744788648803490, -0.381169081397574,  },
            };

            // V can be different, but for the economy SVD it is often equal
            Assert.IsTrue(Matrix.IsEqual(target.RightSingularVectors.Submatrix(0, 3, 0, 1), V, 0.0001));


            double[][] S = 
            {
                new double[] { 14.2690954992615, 0.000000000000000 },
                new double[] {  0.0000000000000,	0.626828232417543 },
            };

            // The diagonal values should be equal
            Assert.IsTrue(Matrix.IsEqual(target.Diagonal.First(2), Matrix.Diagonal(S), 0.001));
        }


        [Test]
        public void JaggedSingularValueDecompositionConstructorTest3()
        {
            // Test using SVD assumption auto-correction feature.

            // Test for m-x-n matrices where m < n. The available SVD
            // routine was not meant to be used in this case.

            double[][] value = new double[][]
            { 
                 new double[] { 1, 2 },
                 new double[] { 3, 4 },
                 new double[] { 5, 6 },
                 new double[] { 7, 8 }
            }.Transpose(); // value is 2x4, having less rows than columns.

            var target = new JaggedSingularValueDecomposition(value, true, true, true);

            double[][] actual = Matrix.Multiply(
                Matrix.Multiply(target.LeftSingularVectors, target.DiagonalMatrix), 
                target.RightSingularVectors.Transpose());

            // Checking the decomposition
            Assert.IsTrue(Matrix.IsEqual(actual, value, 1e-2));
            Assert.IsTrue(Matrix.IsEqual(value, target.Reverse(), 1e-2));

            // Checking values
            double[][] U =
            {
                new double[] { 0.641423027995072, -0.767187395072177 },
                new double[] { 0.767187395072177,  0.641423027995072 },
            };

            // U should be equal despite some sign changes
            Assert.IsTrue(Matrix.IsEqual(target.LeftSingularVectors, U, 0.001));


            double[][] V = // economy svd
            {
                new double[] {  0.152483233310201,  0.822647472225661,  },
                new double[] {  0.349918371807964,  0.421375287684580,  },
                new double[] {  0.547353510305727,  0.0201031031435023, },
                new double[] {  0.744788648803490, -0.381169081397574,  },
            };

            // V can be different, but for the economy SVD it is often equal
            Assert.IsTrue(Matrix.IsEqual(target.RightSingularVectors, V, 0.0001));


            double[][] S = 
            {
                new double[] { 14.2690954992615, 0.000000000000000 },
                new double[] {  0.0000000000000,	0.626828232417543 },
            };

            // The diagonal values should be equal
            Assert.IsTrue(Matrix.IsEqual(target.Diagonal, Matrix.Diagonal(S), 0.001));
        }


        [Test]
        public void JaggedSingularValueDecompositionConstructorTest2()
        {
            // test for m-x-n matrices where m > n (4 > 2)

            double[][] value = new double[][]
            { 
                 new double[] { 1, 2 },
                 new double[] { 3, 4 },
                 new double[] { 5, 6 },
                 new double[] { 7, 8 }
             }; // value is 4x2, thus having more rows than columns


            var target = new JaggedSingularValueDecomposition(value, true, true, false);

            double[][] actual = Matrix.Multiply(Matrix.Multiply(target.LeftSingularVectors,
                                target.DiagonalMatrix),
                                target.RightSingularVectors.Transpose());

            // Checking the decomposition
            Assert.IsTrue(Matrix.IsEqual(actual, value, 1e-2));
            Assert.IsTrue(Matrix.IsEqual(value, target.Reverse(), 1e-5));

            double[][] U = // economy svd
            {
                new double[] {  0.152483233310201,  0.822647472225661,  },
                new double[] {  0.349918371807964,  0.421375287684580,  },
                new double[] {  0.547353510305727,  0.0201031031435023, },
                new double[] {  0.744788648803490, -0.381169081397574,  },
            };

            // U should be equal except for some sign changes
            Assert.IsTrue(Matrix.IsEqual(target.LeftSingularVectors, U, 0.001));



            // Checking values
            double[][] V =
            {
                new double[] {  0.641423027995072, -0.767187395072177 },
                new double[] {  0.767187395072177,  0.641423027995072 },
            };

            // V should be equal except for some sign changes
            Assert.IsTrue(Matrix.IsEqual(target.RightSingularVectors, V, 0.0001));


            double[][] S = 
            {
                new double[] { 14.2690954992615, 0.000000000000000 },
                new double[] {  0.0000000000000,	0.626828232417543 },
            };

            // The diagonal values should be equal
            Assert.IsTrue(Matrix.IsEqual(target.Diagonal, Matrix.Diagonal(S), 0.001));
        }


        [Test]
        public void JaggedSingularValueDecompositionConstructorTest4()
        {
            // Test using SVD assumption auto-correction feature
            // without computing the right singular vectors.

            double[][] value = new double[][]
            { 
                 new double[] { 1, 2 },
                 new double[] { 3, 4 },
                 new double[] { 5, 6 },
                 new double[] { 7, 8 }
            }.Transpose(); // value is 2x4, having less rows than columns.

            var target = new JaggedSingularValueDecomposition(value, true, false, true);


            // Checking values
            double[][] U =
            {
                new double[] { 0.641423027995072, -0.767187395072177 },
                new double[] { 0.767187395072177,  0.641423027995072 },
            };

            // U should be equal despite some sign changes
            Assert.IsTrue(Matrix.IsEqual(target.LeftSingularVectors, U, 0.001));


            // Checking values
            double[][] V =
            {
                new double[] {  0.0, 0.0 },
                new double[] {  0.0, 0.0 },
                new double[] {  0.0, 0.0 },
                new double[] {  0.0, 0.0 },
            };

            // V should not have been computed.
            Assert.IsTrue(Matrix.IsEqual(target.RightSingularVectors, V));


            double[][] S = 
            {
                new double[] { 14.2690954992615, 0.000000000000000 },
                new double[] {  0.0000000000000,	0.626828232417543 },
            };

            // The diagonal values should be equal
            Assert.IsTrue(Matrix.IsEqual(target.Diagonal, Matrix.Diagonal(S), 0.001));
        }

        [Test]
        public void JaggedSingularValueDecompositionConstructorTest5()
        {
            // Test using SVD assumption auto-correction feature
            // without computing the left singular vectors.

            var value = new double[][]
            { 
                 new double[] { 1, 2 },
                 new double[] { 3, 4 },
                 new double[] { 5, 6 },
                 new double[] { 7, 8 }
            }.Transpose(); // value is 2x4, having less rows than columns.

            var target = new JaggedSingularValueDecomposition(value, false, true, true);


            // Checking values
            double[][] U =
            {
                new double[] { 0.0, 0.0 },
                new double[] { 0.0, 0.0 },
            };

            // U should not have been computed
            Assert.IsTrue(Matrix.IsEqual(target.LeftSingularVectors, U));


            double[][] V = // economy svd
            {
                new double[] {  0.152483233310201,  0.822647472225661,  },
                new double[] {  0.349918371807964,  0.421375287684580,  },
                new double[] {  0.547353510305727,  0.0201031031435023, },
                new double[] {  0.744788648803490, -0.381169081397574,  },
            };

            // V can be different, but for the economy SVD it is often equal
            Assert.IsTrue(Matrix.IsEqual(target.RightSingularVectors, V, 0.0001));



            double[][] S = 
            {
                new double[] { 14.2690954992615, 0.000000000000000 },
                new double[] {  0.0000000000000,	0.626828232417543 },
            };

            // The diagonal values should be equal
            Assert.IsTrue(Matrix.IsEqual(target.Diagonal, Matrix.Diagonal(S), 0.001));
        }


        [Test]
        public void JaggedSingularValueDecompositionConstructorTest6()
        {
            // Test using SVD assumption auto-correction feature in place

            double[][] value1 =
            {
                new double[] { 2.5,  2.4 },
                new double[] { 0.5,  0.7 },
                new double[] { 2.2,  2.9 },
                new double[] { 1.9,  2.2 },
                new double[] { 3.1,  3.0 },
                new double[] { 2.3,  2.7 },
                new double[] { 2.0,  1.6 },
                new double[] { 1.0,  1.1 },
                new double[] { 1.5,  1.6 },
                new double[] { 1.1,  0.9 }
            };

            var value2 = value1.Transpose();

            var cvalue1 = value1.Copy();
            var cvalue2 = value2.Copy();

            var target1 = new JaggedSingularValueDecomposition(cvalue1, true, true, true, true);
            var target2 = new JaggedSingularValueDecomposition(cvalue2, true, true, true, true);

            Assert.IsFalse(value1.IsEqual(cvalue1, 1e-3));
            Assert.IsTrue(value2.IsEqual(cvalue2, 1e-3)); // due to auto-transpose

            Assert.IsTrue(target1.LeftSingularVectors.IsEqual(target2.RightSingularVectors));
            Assert.IsTrue(target1.RightSingularVectors.IsEqual(target2.LeftSingularVectors));
            Assert.IsTrue(target1.DiagonalMatrix.IsEqual(target2.DiagonalMatrix));
            Assert.IsTrue(Matrix.IsEqual(value1, target1.Reverse(), 1e-2));
            Assert.IsTrue(Matrix.IsEqual(value2, target2.Reverse(), 1e-2));

            Assert.AreSame(target1.DiagonalMatrix, target1.DiagonalMatrix);
        }

        [Test]
        public void JaggedSingularValueDecompositionConstructorTest7()
        {
            int count = 100;
            double[][] value = new double[count][];
            double[] output = new double[count];

            for (int i = 0; i < count; i++)
            {
                value[i] = new double[3];

                double x = i + 1;
                double y = 2 * (i + 1) - 1;
                value[i][0] = x;
                value[i][1] = y;
                value[i][2] = 1;
                output[i] = 4 * x - y + 3;
            }



            var target = new JaggedSingularValueDecomposition(value,
                computeLeftSingularVectors: true,
                computeRightSingularVectors: true);

            {
                double[][] expected = value;
                double[][] actual = Matrix.Multiply(Matrix.Multiply(target.LeftSingularVectors,
                    target.DiagonalMatrix),
                    target.RightSingularVectors.Transpose());

                // Checking the decomposition
                Assert.IsTrue(Matrix.IsEqual(actual, expected, 1e-8));
            }

            {
                double[] solution = target.Solve(output);

                double[] expected = output;
                double[] actual = value.Multiply(solution);

                Assert.IsTrue(Matrix.IsEqual(actual, expected, 1e-8));
            }
        }


        [Test]
        public void solve_for_diagonal()
        {
            int count = 3;
            double[][] value = new double[count][];
            double[] output = new double[count];

            for (int i = 0; i < count; i++)
            {
                value[i] = new double[3];

                double x = i + 1;
                double y = 2 * (i + 1) - 1;
                value[i][0] = x;
                value[i][1] = y;
                value[i][2] = System.Math.Pow(x, 2);
                output[i] = 4 * x - y + 3;
            }



            var target = new JaggedSingularValueDecomposition(value,
                computeLeftSingularVectors: true,
                computeRightSingularVectors: true);

            {
                double[][] expected = value;
                double[][] actual = Matrix.Multiply(Matrix.Multiply(target.LeftSingularVectors,
                    target.DiagonalMatrix),
                    target.RightSingularVectors.Transpose());

                // Checking the decomposition
                Assert.IsTrue(Matrix.IsEqual(actual, expected, 1e-8));
            }

            {
                double[][] expected = value.Inverse();
                double[][] actual = Matrix.Multiply(Matrix.Multiply(
                    target.RightSingularVectors.Transpose().Inverse(),
                    target.DiagonalMatrix.Inverse()),
                    target.LeftSingularVectors.Inverse()
                    );

                // Checking the invers decomposition
                Assert.IsTrue(Matrix.IsEqual(actual, expected, 1e-8));
            }


            {
                double[][] solution = target.SolveForDiagonal(output);

                double[][] expected = Jagged.Diagonal(output);
                double[][] actual = value.Dot(solution);

                Assert.IsTrue(Matrix.IsEqual(actual, expected, 1e-8));
            }
        }

        [Test]
        public void issue_614()
        {
            // https://github.com/accord-net/framework/issues/614

            double[][] A =
            {
                new double[] { 1 },
                new double[] { 0 }
            };

            double[][] B =
            {
                new double[] { 1 },
                new double[] { 0 }
            };


            double[][] X = Matrix.Solve(A, B, true);

            double[][] expected =
            {
                new double[] { 1 }
            };

            Assert.IsTrue(expected.IsEqual(X));

            X = new JaggedSingularValueDecomposition(A).Solve(B);

            Assert.IsTrue(expected.IsEqual(X));
        }
    }
}
