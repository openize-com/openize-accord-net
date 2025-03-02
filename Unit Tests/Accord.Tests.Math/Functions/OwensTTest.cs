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

namespace Openize.Accord.Tests.Math.Functions
{
    using NUnit.Framework;
    using Openize.Accord.Math.Functions;

    [TestFixture]
    public class OwensTTest
    {


        private TestContext testContextInstance;

        public TestContext TestContext
        {
            get
            {
                return this.testContextInstance;
            }
            set
            {
                this.testContextInstance = value;
            }
        }

        
        [Test]
        public void ExampleTest()
        {
            double x = OwensT.Function(h: 2, a: 42);
            Assert.AreEqual(0.011375065974089608, x);
        }

        [Test]
        public void FunctionTest()
        {

            double[] a_vec =
            {
                0.2500000000000000E+00,
                0.4375000000000000E+00,
                0.9687500000000000E+00,
                0.0625000000000000E+00,
                0.5000000000000000E+00,
                0.9999975000000000E+00,
                0.5000000000000000E+00,
                0.1000000000000000E+01,
                0.2000000000000000E+01,
                0.3000000000000000E+01,
                0.5000000000000000E+00,
                0.1000000000000000E+01,
                0.2000000000000000E+01,
                0.3000000000000000E+01,
                0.5000000000000000E+00,
                0.1000000000000000E+01,
                0.2000000000000000E+01,
                0.3000000000000000E+01,
                0.5000000000000000E+00,
                0.1000000000000000E+01,
                0.2000000000000000E+01,
                0.3000000000000000E+01,
                0.5000000000000000E+00,
                0.1000000000000000E+01,
                0.2000000000000000E+01,
                0.3000000000000000E+01,
                0.1000000000000000E+02,
                0.1000000000000000E+03 
            };

            double[] h_vec =
            {
                0.0625000000000000E+00,
                6.5000000000000000E+00,
                7.0000000000000000E+00,
                4.7812500000000000E+00,
                2.0000000000000000E+00,
                1.0000000000000000E+00,
                0.1000000000000000E+01,
                0.1000000000000000E+01,
                0.1000000000000000E+01,
                0.1000000000000000E+01,
                0.5000000000000000E+00,
                0.5000000000000000E+00,
                0.5000000000000000E+00,
                0.5000000000000000E+00,
                0.2500000000000000E+00,
                0.2500000000000000E+00,
                0.2500000000000000E+00,
                0.2500000000000000E+00,
                0.1250000000000000E+00,
                0.1250000000000000E+00,
                0.1250000000000000E+00,
                0.1250000000000000E+00,
                0.7812500000000000E-02,
                0.7812500000000000E-02,
                0.7812500000000000E-02,
                0.7812500000000000E-02,
                0.7812500000000000E-02,
                0.7812500000000000E-02 
            };

            double[] t_vec = 
            {
                3.8911930234701366E-02,
                2.0005773048508315E-11,
                6.3990627193898685E-13,
                1.0632974804687463E-07,
                8.6250779855215071E-03,
                6.6741808978228592E-02,
                0.4306469112078537E-01,
                0.6674188216570097E-01,
                0.7846818699308410E-01,
                0.7929950474887259E-01,
                0.6448860284750376E-01,
                0.1066710629614485E+00,
                0.1415806036539784E+00,
                0.1510840430760184E+00,
                0.7134663382271778E-01,
                0.1201285306350883E+00,
                0.1666128410939293E+00,
                0.1847501847929859E+00,
                0.7317273327500385E-01,
                0.1237630544953746E+00,
                0.1737438887583106E+00,
                0.1951190307092811E+00,
                0.7378938035365546E-01,
                0.1249951430754052E+00,
                0.1761984774738108E+00,
                0.1987772386442824E+00,
                0.2340886964802671E+00,
                0.2479460829231492E+00 
            };

            double[,] table = 
            {
                //    H             A                    T                        T  
                //                                  (Tabulated)                 (TFN)               DIFF
                { 0.0625,        0.2500,        0.0389119302347014,        0.0389119302347014,           0 },
                { 6.5000,        0.4375,        0.0000000000200058,        0.0000000000100105,   9.995e-12 },
                { 7.0000,        0.9688,        0.0000000000006399,        0.0000000000003200,     3.2e-13 },
                { 4.7812,        0.0625,        0.0000001063297480,        0.0000001063297480,   2.647e-23 },
                { 2.0000,        0.5000,        0.0086250779855215,        0.0086250779855215,   1.735e-18 },
                { 1.0000,        1.0000,        0.0667418089782286,        0.0667418089782286,           0 },
                { 1.0000,        0.5000,        0.0430646911207854,        0.0430646911207854,           0 },
                { 1.0000,        1.0000,        0.0667418821657010,        0.0667418821657010,           0 },
                { 1.0000,        2.0000,        0.0784681869930841,        0.0784681869930841,           0 },
                { 1.0000,        3.0000,        0.0792995047488726,        0.0792995047488726,   1.388e-17 },
                { 0.5000,        0.5000,        0.0644886028475038,        0.0644886028475038,   1.388e-17 },
                { 0.5000,        1.0000,        0.1066710629614485,        0.1066710629614485,   1.388e-17 },
                { 0.5000,        2.0000,        0.1415806036539784,        0.1415806036539784,           0 },
                { 0.5000,        3.0000,        0.1510840430760184,        0.1510840430760184,   2.776e-17 },
                { 0.2500,        0.5000,        0.0713466338227178,        0.0713466338227178,   2.776e-17 },
                { 0.2500,        1.0000,        0.1201285306350883,        0.1201285306350883,           0 },
                { 0.2500,        2.0000,        0.1666128410939293,        0.1666128410939293,   2.776e-17 },
                { 0.2500,        3.0000,        0.1847501847929859,        0.1847501847929859,   2.776e-17 },
                { 0.1250,        0.5000,        0.0731727332750039,        0.0731727332750039,           0 },
                { 0.1250,        1.0000,        0.1237630544953746,        0.1237630544953746,   2.776e-17 },
                { 0.1250,        2.0000,        0.1737438887583106,        0.1737438887583106,           0 },
                { 0.1250,        3.0000,        0.1951190307092811,        0.1951190307092811,   2.776e-17 },
                { 0.0078,        0.5000,        0.0737893803536555,        0.0737893803536555,   1.388e-17 },
                { 0.0078,        1.0000,        0.1249951430754052,        0.1249951430754052,           0 },
                { 0.0078,        2.0000,        0.1761984774738108,        0.1761984774738108,   2.776e-17 },
                { 0.0078,        3.0000,        0.1987772386442824,        0.1987772386442823,   5.551e-17 },
                { 0.0078,       10.0000,        0.2340886964802671,        0.2340886964802671,   2.776e-17 },
                { 0.0078,      100.0000,        0.2479460829231492,        0.2479460829231492,           0 },
            };

            for (int i = 0; i < table.GetLength(0); i++)
            {
                double H = h_vec[i];
                double A = a_vec[i];

                double tabulated = t_vec[i];
                double expected = table[i, 3];
                double actual = OwensT.Function(H, A);

                Assert.AreEqual(expected, actual, 1e-10);
            }
        }

    }
}
