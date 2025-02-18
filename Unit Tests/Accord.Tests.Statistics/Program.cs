namespace Openize.Accord.Tests.Statistics
{
    using System;

    public static class Program
    {
        public static void Main2(string[] args)
        {
            Console.ReadLine();
#if NETCORE2_0
            // Dynamically enumerate all univariate distributions and make
            // sure the distribution contains properties for every exposed
            // constructors' parameters

            var baseType = typeof(IUnivariateDistribution);

            // Get all univariate distributions in Accord.NET:
            var distributions = GetDerivedConcreteTypes(baseType);

            var check = new HashSet<Type>();

            foreach (Type type in distributions)
            {
                if (!baseType.IsAssignableFrom(type))
                    continue;

                foreach (var ctor in type.GetConstructors())
                {
                    foreach (var param in ctor.GetParameters())
                    {
                        string distName = type.Name;
                        string paramName = param.Name;

                        // Search for a property with a similar name
                        var p = type.GetProperties().FirstOrDefault(
                            prop => cmp(distName, paramName, prop.Name));

                        bool found = p != null;

                        Assert.IsTrue(found, $"{distName}.{paramName}");
                        check.Add(type);
                    }
                }
            }
#endif
        }
    }
}
