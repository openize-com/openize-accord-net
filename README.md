# Openize.Accord

**[Openize.Aсcord](https://github.com/Openize-Accord/Openize.Accord-for-.NET)**  This is a fork of the [Accord.Net](https://github.com/accord-net) project, which includes the Imaging library. This version only  NetStandard 2.0 framework, and uses [Aspose.Drawing](https://docs.aspose.com/drawing/net/) as a graphics engine, which allows you to create cross-platform applications using lastest .Net platforms.

## Platform dependence

Openize.Accord can be used to develop applications on Windows Desktop (x86, x64), Windows Server (x86, x64), Windows Azure, Windows Embedded (CE 6.0 R2), as well as Linux x64. The supported platforms include Net Core 3.1, Net6.0, Net7.0, Net8.0.

## New Features & Enhancements in Version 25.2
 - Change FileFormat to Openize
 - Update Aspose.Drawing Engine

## Getting Started with Openize.Accord for .NET

Are you ready to give Openize.Accord a try? Simply execute 

```
Install-Package Openize.Accord
```

from Package Manager Console in Visual Studio to fetch the NuGet package. If you already have Accord.Imaging.Net and want to upgrade the version, please execute 

```
Update-Package Openize.Accord
```

 to get the latest version.
 
## Product License
 - [Openize.Accord is distributed under LGPL license](http://www.gnu.org/licenses/lgpl.html)
 - [Aspose.Drawing .NET is distributed under Aspose EULA license](https://www.conholdate.app/viewer/view/4Y8UNm7laVFjMAd0r/aspose_end-user-license-agreement_2024-05-16.pdf);


## Usage example

### Apply GrayWorld filter

```csharp
using System.Drawing.AsposeDrawing;
using Openize.Accord.Imaging.Filters;
using Image = Openize.Accord.Imaging.AForge.Imaging.Image;

//Set license Aspose.Drawing
System.Drawing.AsposeDrawing.License lic = new License();
lic.SetLicense("license.lic");

//Load image
using (var image = Image.FromFile("lena_color.jpg"))
{   
	//Create filter
    var grayWorld = new GrayWorld();
	
	//apply filter
    grayWorld.ApplyInPlace(image);
	
	//save
    image.Save("lena_gray.jpg");
}
```

You can see other examples directly on the [Accord website](http://accord-framework.net/samples.html)