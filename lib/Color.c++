#include "Color.h"
#include "progress_bar.hpp"

int backgroundWidth, backgroundHeight;

/*
int equiProject( double* d )
{
  double x = d[0];
  double y = d[1];
  double z = d[2];

  double r = sqrt( x * x + y * y + z * z );
  double lat = acos( z / r )-3.141592658979/2;
  double lng = atan2( y, x );

  double j = (lng + 3.141592658979) * 2000 / (2*3.141592658979);
  double i = (3.141592658979/2 - lat) * 1000 / 3.141592658979;

  return i * 2000 + j;
}
*/

int minJ=1000000, maxJ= -1000000, minI=1000000, maxI=-1000000;

int equiProject( double* d )
{
  double x = d[0];
  double y = d[1];
  double z = d[2];

  double lat = acos( y )/3.141592658979;
  double lng = atan2( z, x )/(2.0*3.141592658979)+0.5;

//  lng += 0.25;
//  if( lng > 1.0 ) lng -= 1.0;

  double j1 = lng*backgroundWidth;
  double i1 = lat*backgroundHeight;
  
  int j = j1;
  int i = i1;

//  minJ = min(j,minJ);
//  maxJ = max(j,maxJ);

//  minI = min(i,minI);
//  maxI = max(i,maxI);

//  cout << minJ << " " << maxJ << " " << minI << " " << maxI << endl;

//  j *= 2000;
//  i *= 1000;

//  cout << d[0] << " " << d[1] << " " << d[2] << endl;
//  cout << j << " " << i << " " << i * 2000 + j << endl;

//  cout << lng << " " << lat << endl;

  return i * backgroundWidth + j;
}

double interpolate(double r1, double r2, double a)
{
  return (1 - a) * r1 + a * r2;
}

double length( double* z )
{
  return sqrt(z[0]*z[0] + z[1]*z[1] + z[2]*z[2]);
}

void crossProduct( double* a, double* b, double* c )
{
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}

double dotProduct( double* a, double* b )
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

void normalize( double* n )
{
  double l = length( n );

  for( int k1 = 0; k1 < 3; ++k1 )
  {
    n[k1] /= l;
  }
}

void randomFromNormal( double* n, double r1, double r2 )
{
  double t[3];
  if( abs( n[0] ) > abs( n[1] ) )
  {
    double nl = sqrt( n[0]*n[0] + n[2]*n[2] );
    t[0] = n[2]/nl;
    t[1] = 0;
    t[2] = -n[0]/nl;
  }
  else
  {
    double nl = sqrt( n[1]*n[1] + n[2]*n[2] );
    t[0] = 0;
    t[1] = -n[2]/nl;
    t[2] = n[1]/nl;
  }

  double b[3];

  crossProduct( n, t, b );

//  double r1 = ((double) rand() / (RAND_MAX));
//  double r2 = ((double) rand() / (RAND_MAX));


//  double seedX = seed - 1;
//  double seedY = seed + 1;

/*
  seed[0] -= 1;
  seed[1] += 1;

  double r1 = sin((seed[0] * 12.9898 + seed[1] * 78.233))*43758.5453;
  r1 = r1 - floor(r1);

  double r2 = cos((seed[0] * 4.898 + seed[1] * 7.23))*23421.631;
  r2 = r2 - floor(r2);

  r1 = sqrt(r1);
*/
//  cout << r1 << " " << r2 << endl;

//  r1 = r1*2.0 - 1.0;
//  r2 = r2*2.0 - 1.0;

  double sinTheta = sqrt( 1 - r1*r1 );
  double phi = 2 * 3.141592658979 * r2;

  double x = sinTheta * cos( phi );
  double y = r1;
  double z = sinTheta * sin( phi );

  double nx = x*b[0] + y*n[0] + z*t[0];
  double ny = x*b[1] + y*n[1] + z*t[1];
  double nz = x*b[2] + y*n[2] + z*t[2];

  n[0] = nx;
  n[1] = ny;
  n[2] = nz;

  normalize(n);
}

void Color( string imageName, string outputName )
{
  VImage image = VImage::vipsload((char *)imageName.c_str());

  unsigned char * backgroundData = ( unsigned char * )image.data();

  backgroundWidth = image.width();
  backgroundHeight = image.height();

  double maxColor = 0.0;
  double maxDirection[3];
  int maxColor2[3];

  for( double r1 = 0.0; r1 <= 1.0; r1 += 0.0001 )
  {
    for( double r2 = 0.0; r2 <= 1.0; r2 += 0.0001 )
    {
      double n[3] = {0,1,0};

      randomFromNormal( n, r1, r2 );

      int index = equiProject( n );

      int m = backgroundData[4*index+3];
      float mantissa = m - 128;

      double totalColor = 0;

      double color2[3];

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

  for( double r1 = 0.0; r1 <= 1.0; r1 += 0.0001 )
  {
    for( double r2 = 0.0; r2 <= 1.0; r2 += 0.0001 )
    {
      double n[3] = {0,-1,0};

      randomFromNormal( n, r1, r2 );

      int index = equiProject( n );

      int m = backgroundData[4*index+3];
      float mantissa = m - 128;

      double totalColor = 0;

      double color2[3];

      for( int k1 = 0; k1 < 3; ++k1 )
      {
        int f = backgroundData[4*index+k1];
        float f1 = f * pow(2.0,mantissa);
//        color[k1] *= f1;
//        cout << 255.0*direct[k1] << endl;
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
  data.write( (char *)maxDirection, sizeof(double)*3 );
  data.write( (char *)maxColor2, sizeof(int)*3 );

  data.close();

//  cout << maxColor << " " << maxDirection[0] << " " << maxDirection[1] << " " << maxDirection[2] << endl;
//  cout << maxColor << " " << maxColor2[0] << " " << maxColor2[1] << " " << maxColor2[2] << endl;
}

void RunColor( string imageName, string outputName )
{
  Color( imageName, outputName );
}