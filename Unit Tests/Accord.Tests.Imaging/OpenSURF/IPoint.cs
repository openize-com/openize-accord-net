﻿namespace Openize.Accord.Tests.Imaging.OpenSURF
{
    public class IPoint
  {
    /// <summary>
    /// Default ctor
    /// </summary>
    public IPoint()
    {
      this.orientation = 0;
    }

    /// <summary>
    /// Coordinates of the detected interest point
    /// </summary>
    public float x, y;

    /// <summary>
    /// Detected scale
    /// </summary>
    public float scale;

    /// <summary>
    /// Response of the detected feature (strength)
    /// </summary>
    public float response;

    /// <summary>
    /// Orientation measured anti-clockwise from the x-axis
    /// </summary>
    public float orientation;

    /// <summary>
    /// Sign of Laplacian for fast matching purposes
    /// </summary>
    public int laplacian;

    /// <summary>
    /// Descriptor vector
    /// </summary>
    public int descriptorLength;
    public float [] descriptor = null;
    public void SetDescriptorLength(int Size)
    {
      this.descriptorLength = Size;
      this.descriptor = new float[Size];
    }
  }
}
