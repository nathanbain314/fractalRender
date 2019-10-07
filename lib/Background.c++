#include "Background.h"
#include "progress_bar.hpp"

int backgroundWidth, backgroundHeight;

int equiProject( float* d )
{
  float x = d[0];
  float y = d[1];
  float z = d[2];

  float lat = acos( y )/3.141592658979;
  float lng = atan2( z, x )/(2.0*3.141592658979)+0.5;

  float j1 = lng*backgroundWidth;
  float i1 = lat*backgroundHeight;
  
  int j = j1;
  int i = i1;

  return i * backgroundWidth + j;
}

float interpolate(float r1, float r2, float a)
{
  return (1 - a) * r1 + a * r2;
}

float length( float* z )
{
  return sqrt(z[0]*z[0] + z[1]*z[1] + z[2]*z[2]);
}

void crossProduct( float* a, float* b, float* c )
{
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}

float dotProduct( float* a, float* b )
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

void normalize( float* n )
{
  float l = length( n );

  for( int k1 = 0; k1 < 3; ++k1 )
  {
    n[k1] /= l;
  }
}

void randomFromNormal( float* n, float r1, float r2 )
{
  float t[3];
  if( abs( n[0] ) > abs( n[1] ) )
  {
    float nl = sqrt( n[0]*n[0] + n[2]*n[2] );
    t[0] = n[2]/nl;
    t[1] = 0;
    t[2] = -n[0]/nl;
  }
  else
  {
    float nl = sqrt( n[1]*n[1] + n[2]*n[2] );
    t[0] = 0;
    t[1] = -n[2]/nl;
    t[2] = n[1]/nl;
  }

  float b[3];

  crossProduct( n, t, b );

  float sinTheta = sqrt( 1 - r1*r1 );
  float phi = 2 * 3.141592658979 * r2;

  float x = sinTheta * cos( phi );
  float y = r1;
  float z = sinTheta * sin( phi );

  float nx = x*b[0] + y*n[0] + z*t[0];
  float ny = x*b[1] + y*n[1] + z*t[1];
  float nz = x*b[2] + y*n[2] + z*t[2];

  n[0] = nx;
  n[1] = ny;
  n[2] = nz;

  normalize(n);
}

void Background( string imageName, string outputName )
{
  VImage image = VImage::vipsload((char *)imageName.c_str());

  unsigned char * backgroundData = ( unsigned char * )image.data();

  backgroundWidth = image.width();
  backgroundHeight = image.height();

  float maxColor = 0.0;
  float maxDirection[3];
  int maxColor2[3];

  for( float r1 = 0.0; r1 <= 1.0; r1 += 0.0001 )
  {
    for( float r2 = 0.0; r2 <= 1.0; r2 += 0.0001 )
    {
      float n[3] = {0,1,0};

      randomFromNormal( n, r1, r2 );

      int index = equiProject( n );

      int m = backgroundData[4*index+3];
      float mantissa = m - 128;

      float totalColor = 0;

      float color2[3];

      for( int k1 = 0; k1 < 3; ++k1 )
      {
        int f = backgroundData[4*index+k1];
        float f1 = f * pow(2.0,mantissa);
        totalColor += f1 * f1;
        color2[k1] = f;
      }

      if( totalColor > maxColor )
      {
        maxDirection[0] = n[0];
        maxDirection[1] = n[1];
        maxDirection[2] = n[2];
        maxColor2[0] = color2[0];
        maxColor2[1] = color2[1];
        maxColor2[2] = color2[2];
        maxColor = totalColor;
      }
    }
  }

  for( float r1 = 0.0; r1 <= 1.0; r1 += 0.0001 )
  {
    for( float r2 = 0.0; r2 <= 1.0; r2 += 0.0001 )
    {
      float n[3] = {0,-1,0};

      randomFromNormal( n, r1, r2 );

      int index = equiProject( n );

      int m = backgroundData[4*index+3];
      float mantissa = m - 128;

      float totalColor = 0;

      float color2[3];

      for( int k1 = 0; k1 < 3; ++k1 )
      {
        int f = backgroundData[4*index+k1];
        float f1 = f * pow(2.0,mantissa);
        totalColor += f1 * f1;
        color2[k1] = f;
      }

      if( totalColor > maxColor )
      {
        maxDirection[0] = n[0];
        maxDirection[1] = n[1];
        maxDirection[2] = n[2];
        maxColor2[0] = color2[0];
        maxColor2[1] = color2[1];
        maxColor2[2] = color2[2];
        maxColor = totalColor;
      }
    }
  }

  ofstream data( outputName, ios::binary );
  data.write( (char *)&backgroundWidth, sizeof(int) );
  data.write( (char *)&backgroundHeight, sizeof(int) );
  data.write( (char *)backgroundData, sizeof(unsigned char)*4*backgroundWidth*backgroundHeight );
  data.write( (char *)maxDirection, sizeof(float)*3 );
  data.write( (char *)maxColor2, sizeof(int)*3 );

  data.close();
}

void RunBackground( string imageName, string outputName )
{
  Background( imageName, outputName );
}