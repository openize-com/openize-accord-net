// AForge Image Processing Library
// AForge.NET framework
//
// Copyright � Andrew Kirillov, 2005-2008
// andrew.kirillov@aforgenet.com
//

namespace Openize.Accord.Imaging.AForge.Imaging.Filters.Convolution
{
	/// <summary>
	/// Simple edge detector.
	/// </summary>
	/// 
    /// <remarks><para>The filter performs <see cref="Convolution">convolution filter</see> using
    /// the edges kernel:</para>
    /// 
    /// <code lang="none">
    ///  0  -1   0
    /// -1   4  -1
    ///  0  -1   0
    /// </code>
    /// 
    /// <para>For the list of supported pixel formats, see the documentation to <see cref="Convolution"/>
    /// filter.</para>
    /// 
    /// <para>Sample usage:</para>
    /// <code>
    /// // create filter
    /// Edges filter = new Edges( );
    /// // apply the filter
    /// filter.ApplyInPlace( image );
    /// </code>
    ///
    /// <para><b>Initial image:</b></para>
    /// <img src="..\images\imaging\sample1.jpg" width="480" height="361" />
    /// <para><b>Result image:</b></para>
    /// <img src="..\images\imaging\edges.png" width="480" height="361" />
    /// </remarks>
    /// 
    /// <seealso cref="Convolution"/>
    ///
    public sealed class Edges : Convolution
	{
		/// <summary>
		/// Initializes a new instance of the <see cref="Edges"/> class.
		/// </summary>
		public Edges( ) : base( new int[,] {
										{  0, -1,  0 },
										{ -1,  4, -1 },
										{  0, -1,  0 } } )
		{
		}
	}
}
