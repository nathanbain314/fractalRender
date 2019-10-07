#include "Gradient.h"
#include "progress_bar.hpp"

void Gradient( string imageName, string outputName )
{
  // Load input image
  VImage image = VImage::vipsload( (char *)imageName.c_str() ).autorot();

  // Convert to a three band image
  if( image.bands() == 1 )
  {
    image = image.bandjoin(image).bandjoin(image);
  }
  if( image.bands() == 4 )
  {
    image = image.flatten();
  }

  image = image.resize( 1024.0/image.width() );

  unsigned char * gradientData = ( unsigned char * )image.data();

  ofstream data( outputName, ios::binary );
  data.write( (char *)gradientData, sizeof(unsigned char)*3*1024 );

  data.close();
}

void RunGradient( string imageName, string outputName )
{
  Gradient( imageName, outputName );
}