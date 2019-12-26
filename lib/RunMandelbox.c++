#include "Mandelbox.h"
#include <tclap/CmdLine.h>

using namespace TCLAP;

int main( int argc, char **argv )
{
  try
  {
    CmdLine cmd("Renders a fractal", ' ', "2.0");

    ValueArg<string> phongArg( "p", "phong", "Phong data", false, " ", "string", cmd );

    SwitchArg dontPathTraceArg( "", "raycast", "Use raycasting instead of pathtracing", cmd, false );

    SwitchArg lightingArg( "l", "lighting", "Use direct lighting", cmd, false );

    ValueArg<string> backgroundArg( "b", "background", "Background data", false, "backgrounds/desert.dat", "string", cmd);

    ValueArg<string> gradientArg( "g", "gradient", "Gradient data", false, "gradients/default.dat", "string", cmd);

    ValueArg<string> rectangleArg( "r", "rectangle", "Rectangle to render", false, " ", "string", cmd );

    ValueArg<float> zoomArg( "z", "zoom", "Zoom multiplier", false, 1.0, "float", cmd );

    ValueArg<float> materialArg( "m", "material", "Material reflectance", false, 0.0, "float", cmd );

    ValueArg<float> colorArg( "c", "color", "Color multiplier", false, 0.0, "float", cmd );

    ValueArg<float> valueArg( "v", "value", "Value to use in fractal equation", false, 1, "float", cmd );

    ValueArg<int> maxStepsArg( "", "maxSteps", "Max steps of ray marching", false, 500, "int", cmd );

    ValueArg<int> imageSizeArg( "s", "size", "Size of image", false, 256, "int", cmd );

    ValueArg<int> depthArg( "d", "depth", "Max depth of ray", false, 2, "int", cmd );

    ValueArg<float> iterArg( "i", "iter", "Minimum iteration size to use", false, 0.001, "float", cmd );

    ValueArg<int> sampleArg( "n", "numSamples", "Number of samples to render", false, 100, "int", cmd );

    ValueArg<int> fractalTypeArg( "f", "fractal", "0=Mandelbulb,1=Mandelbox,2=Menger", false, 0, "int", cmd );

    ValueArg<string> outputArg( "o", "output", "Output image name", false, "out.png", "string", cmd);

    cmd.parse( argc, argv );

    string outputName                 = outputArg.getValue();
    string backgroundName             = backgroundArg.getValue();
    string gradientName               = gradientArg.getValue();
    string rectangleData              = rectangleArg.getValue();
    string phongData                  = phongArg.getValue();
    int fractalType                   = fractalTypeArg.getValue();
    int numSamples                    = sampleArg.getValue();
    int maxDepth                      = depthArg.getValue();
    int maxSteps                      = maxStepsArg.getValue();
    int imageSize                     = imageSizeArg.getValue();
    float zoom                        = zoomArg.getValue();
    float minIter                     = iterArg.getValue();
    float value                       = valueArg.getValue();
    float color                       = colorArg.getValue();
    float reflectance                 = materialArg.getValue();
    bool directLighting               = lightingArg.getValue();
    bool dontPathTrace                = dontPathTraceArg.getValue();

    int bx, by, xSize, ySize;

    if( rectangleData == " " )
    {
      bx = 0;
      by = 0;
      xSize = imageSize;
      ySize = imageSize;
    }
    else
    {
      stringstream ss(rectangleData);

      ss >> bx;
      ss.ignore();
      ss >> by;
      ss.ignore();
      ss >> xSize;
      ss.ignore();
      ss >> ySize;
    }

    imageSize *= zoom;
    bx *= zoom;
    by *= zoom;
    xSize *= zoom;
    ySize *= zoom;

    float kd, ks, ka, alpha;

    if( phongData == " " )
    {
      kd = 0.33333333;
      ks = 0.33333333;
      ka = 0.33333333;
      alpha = 10.0;
    }
    else
    {
      stringstream ss(phongData);

      ss >> kd;
      ss.ignore();
      ss >> ks;
      ss.ignore();
      ss >> ka;
      ss.ignore();
      ss >> alpha;

      float kSum = kd + ks + ka;
      kd /= kSum;
      ks /= kSum;
      ka /= kSum;
    }

    if( VIPS_INIT( argv[0] ) ) return( -1 );

    RunMandelbox( outputName, backgroundName, gradientName, directLighting, fractalType, numSamples, maxDepth, minIter, value, color, reflectance, imageSize, bx, by, xSize, ySize, kd, ks, ka, alpha, maxSteps, dontPathTrace );
  }
  catch (ArgException &e)  // catch any exceptions
  {
    cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
  }

  return 0;
}