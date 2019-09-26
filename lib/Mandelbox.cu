#include "Mandelbox.h"
#include "progress_bar.hpp"

__device__
int equiProject( double* d, int backgroundWidth, int backgroundHeight )
{
  double x = d[0];
  double y = d[1];
  double z = d[2];

  double lat = acos( y )/3.141592658979;
  double lng = atan2( z, x )/(2.0*3.141592658979)+0.5;

  double j1 = lng*backgroundWidth;
  double i1 = lat*backgroundHeight;
  
  int j = j1;
  int i = i1;

  return i * backgroundWidth + j;
}

__device__
double interpolate(double r1, double r2, double a)
{
  return (1 - a) * r1 + a * r2;
}

__device__
double length( double* z )
{
  return sqrt(z[0]*z[0] + z[1]*z[1] + z[2]*z[2]);
}

__device__
void crossProduct( double* a, double* b, double* c )
{
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}

__device__
double dotProduct( double* a, double* b )
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}





// Mandelbulb
__device__
double MandelbulbDE( double* z, double n )
{
  z[2] *= -1;

  double z2[3];

  z2[0] = z[0];
  z2[1] = z[1];
  z2[2] = z[2];

  int maxIter = 100;

  double dr = 1.0;

  double r;

  int k;

  for( k = 0; k < maxIter; ++k )
  {
    r = length(z2);

    if( r > 4 ) break;

    double t = acos(z2[2]/r);
    double p = atan2(z2[1],z2[0]);

    dr = pow( r, n-1.0 ) * n * dr + 1.0;

    double r1 = pow(r,n);

    z2[0] = r1*sin(t*n)*cos(p*n) + z[0];
    z2[1] = r1*sin(t*n)*sin(p*n) + z[1];
    z2[2] = r1*cos(t*n) + z[2];
  }

  z[2] *= -1;

  return 0.25 * log(r) * r / dr;
}

// Mandelbox
__device__
double MandelboxDE( double* z1, double scale )
{
  double sideLength = (scale < -1 ? 4 : 4*(scale+1)/(scale-1))*1.25;

  int iter = 0;

  double DEfactor = scale;

  double fR2 = 1.0;
  double mR2 = 0.25;

  double x = z1[0]*sideLength*0.5;
  double y = z1[1]*sideLength*0.5;
  double z = z1[2]*sideLength*0.5;

  double cx = x;
  double cy = y;
  double cz = z;

  for( ; iter < 100; ++iter )
  {    
    if( sqrt(x*x+y*y+z*z) > 1000 ) break;

    if (x > 1.0)
    x = 2.0 - x;
    else if (x < -1.0) x = -2.0 - x;
    if (y > 1.0)
    y = 2.0 - y;
    else if (y < -1.0) y = -2.0 - y;
    if (z > 1.0)
    z = 2.0 - z;
    else if (z < -1.0) z = -2.0 - z;

    double r2 = x*x + y*y + z*z;

    if (r2 < mR2)
    {
       x = x * fR2 / mR2;
       y = y * fR2 / mR2;
       z = z * fR2 / mR2;
       DEfactor = DEfactor * fR2 / mR2;
    }
    else if (r2 < fR2)
    {
       x = x * fR2 / r2;
       y = y * fR2 / r2;
       z = z * fR2 / r2;
       DEfactor *= fR2 / r2;
    }

    x = x * scale + cx;
    y = y * scale + cy;
    z = z * scale + cz;
    DEfactor *= scale;
  }

  return sqrt(x*x+y*y+z*z)/abs(sideLength*DEfactor);
}


// Menger
__device__
double MengerDE( double* point, double n10 )
{
  int n = 6;

  double x=point[0]*1.2, y=point[1]*1.2, z=point[2]*1.2;
  x=x*0.5+0.5;y=y*0.5+0.5;z=z*0.5+0.5;

  double xx=abs(x-0.5)-0.5, yy=abs(y-0.5)-0.5, zz=abs(z-0.5)-0.5;
  double d1=max(xx,max(yy,zz));
  double d=d1;
  double p=1.0;
  for (int i=1; i<=n; ++i) {
    double xa = fmod(3.0*x*p,3.0);
    double ya = fmod(3.0*y*p,3.0);
    double za = fmod(3.0*z*p,3.0);
    p*=3.0;

    double xx=0.5-abs(xa-1.5), yy=0.5-abs(ya-1.5), zz=0.5-abs(za-1.5);
    d1=min(max(xx,zz),min(max(xx,yy),max(yy,zz))) / p;

    d=max(d,d1);
  }

  return d;
}


__device__
double de( double* z, double value, int fractalType )
{
  switch( fractalType )
  {
    case 0:
      return MandelbulbDE( z, value );
    case 1:
      return MandelboxDE( z, value );
    case 2:
      return MengerDE( z, value );
  }

  return MandelbulbDE( z, value );
}



__device__
void normalize( double* n )
{
  double l = length( n );

  for( int k1 = 0; k1 < 3; ++k1 )
  {
    n[k1] /= l;
  }
}

__device__
void computeNormal( double* p, double* d, double value, int fractalType, double *n )
{
  double t[3];

  for( int k1 = 0; k1 < 3; ++k1 )
  {
    t[0] = p[0];
    t[1] = p[1];
    t[2] = p[2];

    t[k1] = p[k1];// + d[k1];

    double nt = de(t,value,fractalType);

    t[k1] = p[k1] - d[k1];
    
    nt -= de(t,value,fractalType);

    n[k1] = nt;
  }

  normalize( n );
}

__device__
void coneSample( double* n, double* seed, double extent )
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

  seed[0] -= 1;
  seed[1] += 1;

  double r1 = sin((seed[0] * 12.9898 + seed[1] * 78.233))*43758.5453;
  r1 = r1 - floor(r1);

  double r2 = cos((seed[0] * 4.898 + seed[1] * 7.23))*23421.631;
  r2 = r2 - floor(r2);

  r1 = 1.0 - r1*extent;

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

__device__
void randomFromNormal( double* n, double* seed )
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

  seed[0] -= 1;
  seed[1] += 1;

  double r1 = sin((seed[0] * 12.9898 + seed[1] * 78.233))*43758.5453;
  r1 = r1 - floor(r1);

  double r2 = cos((seed[0] * 4.898 + seed[1] * 7.23))*23421.631;
  r2 = r2 - floor(r2);

  r1 = sqrt(r1);

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

__device__
void reflect( double* d, double* n )
{
  double v = 2.0 * dotProduct( d, n );
  for( int k1 = 0; k1 < 3; ++k1 )
  {
    d[k1] = d[k1] - v * n[k1];
  }
  normalize(d);
}

__device__
bool findHit( double *p, double* d, double value, int fractalType, double minSize, double* output )
{
  while( abs(p[0]) < 15 && abs(p[1]) < 15 && abs(p[2]) < 15 )
  {
    double iter = de( p, value, fractalType );

    if( iter < minSize )
    {
      normalize(d);

      for( int k1 = 0; k1 < 3; ++k1 )
      {
        d[k1] = minSize/1000;
      }

      computeNormal( p, d, value, fractalType, output );

      return true;
    }

    for( int k1 = 0; k1 < 3; ++k1 )
    {
      p[k1] += d[k1] * iter;
    }
  }

  return false;
}

__global__
void generateMandelboxPoint( int start, int stride, int aliasIndex, int numAlias, int bx, int by, int xSize, int ySize, int size, int sample, double *startData, int *imageData, unsigned char * backgroundData, bool directLighting, double value, double reflectance, int fractalType, double minSize, double ia, double ja, int maxDepth, int backgroundWidth, int backgroundHeight, double* sunDirect2, int* sunColor )
{
  int idx1 = blockIdx.x * blockDim.x + threadIdx.x;
  int idx2 = 7 * ( idx1 * numAlias * numAlias + aliasIndex );
  int idx = idx1*stride + start;
  int yPixel = by + idx / xSize;
  int xPixel = bx + idx % xSize;

  double num1 = 1.732050807569;
  double num2 = 2.121320343560;

  if( idx < xSize*ySize || ( sample > 1 && startData[idx2] > 1.0 ) )
  {
    double i = (double)yPixel + ia - (double) size / 2.0;
    double j = (double)xPixel + ja - (double) size / 2.0;

    double x = 10.0 + num1*(double)i/(double)size + num2*(double)j/(double)size;
    double y = 10.0 - num1*(double)i/(double)size;
    double z = 10.0 + num1*(double)i/(double)size - num2*(double)j/(double)size;

    double minOffset = min( x - 1.5, min( y - 1.5, z - 1.5 ) );

    x -= minOffset;
    y -= minOffset;
    z -= minOffset;

    double p[3] = {x,y,z};
    double d[3] = { -0.57735026919, -0.57735026919, -0.57735026919 };

    double n[3], n2[3];

    if( sample < 2 )
    {
      if( findHit( p, d, value, fractalType, minSize, n ) )
      {
        startData[idx2] = 2.0;

        startData[idx2+1] = p[0];
        startData[idx2+2] = p[1];
        startData[idx2+3] = p[2];

        startData[idx2+4] = n[0];
        startData[idx2+5] = n[1];
        startData[idx2+6] = n[2];

        d[0] = -0.57735026919;
        d[1] = -0.57735026919;
        d[2] = -0.57735026919;
      }
      else
      {
        startData[idx2] = 0.0;
      }
    }
    else
    {
      p[0] = startData[idx2+1];
      p[1] = startData[idx2+2];
      p[2] = startData[idx2+3];

      n[0] = startData[idx2+4];
      n[1] = startData[idx2+5];
      n[2] = startData[idx2+6];
    }



    if( startData[idx2] > 1.0 )
    {
      double seed[2] = { (double)(xPixel*sample)+ja, (double)(yPixel*sample)+ia };

      double color[3] = {1.0,1.0,1.0};

      double direct = 0;

      for( int depth = 0; depth <= maxDepth; ++depth )
      {
        if( !depth || findHit( p, d, value, fractalType, minSize, n ) )
        {
          if( depth == maxDepth ) break;

          seed[0] -= 1;
          seed[1] += 1;

          double r1 = sin((seed[0] * 12.9898 + seed[1] * 78.233))*43758.5453;
          r1 = r1 - floor(r1);

          if( r1 < reflectance )
          {
            reflect( d, n );
          }
          else
          {
            d[0] = n[0];
            d[1] = n[1];
            d[2] = n[2];

            randomFromNormal( d, seed );
          }

          for( int k1 = 0; k1 < 3; ++k1 )
          {
            p[k1] += d[k1] * 2.0 * minSize;
          }

          if( directLighting )
          {
            double sunDirect[3] = { sunDirect2[0], sunDirect2[1], sunDirect2[2] };

            coneSample( sunDirect, seed, 0.00001 );

            double sunLight = dotProduct( sunDirect, n );

            double p2[3] = {p[0],p[1],p[2]};

            if( sunLight > 0 && !findHit( p2, sunDirect, value, fractalType, minSize, n2 ) )
            {
              direct += sunLight;
            }
          }
        }
        else
        {
          int index = equiProject( d, backgroundWidth, backgroundHeight );

          int m = backgroundData[4*index+3];
          int mantissa = m - 128;

          for( int k1 = 0; k1 < 3; ++k1 )
          {
            color[k1] += min(255.0,sunColor[k1]*direct + color[k1] * ldexp((double)backgroundData[4*index+k1],mantissa));
          }

          break;
        }
      }

      for( int k1 = 0; k1 < 3; ++k1 )
      {
        imageData[3*idx1+k1] += color[k1];
      }
    }
  }
}

void MandelboxThread( int start, int stride, int bx, int by, int xSize, int ySize, int size, unsigned char * imageData, unsigned char * backgroundData, bool directLighting, double value, double minSize, double reflectance, int numSamples, int numAlias, int maxDepth, int fractalType, int backgroundWidth, int backgroundHeight, double* sunDirect, int* sunColor, ProgressBar *generatingMandelbrot )
{
  int numToCompute = xSize*ySize/stride;

  int blockSize = 256;
  int numBlocks = (numToCompute+blockSize-1)/blockSize;

  cudaSetDevice(start);
  
  int *imageData2 = (int *)malloc(numToCompute*3*sizeof(int));

  int *imageDataCuda;
  cudaMalloc((void**)&imageDataCuda, numToCompute*3*sizeof(int));

  double *startDataCuda;
  cudaMalloc((void**)&startDataCuda, numToCompute*7*numAlias*numAlias*sizeof(double));

  unsigned char *backgroundDataCuda;
  cudaMalloc((void**)&backgroundDataCuda, backgroundWidth*backgroundHeight*4*sizeof(unsigned char));
  cudaMemcpy(backgroundDataCuda,backgroundData,backgroundWidth*backgroundHeight*4*sizeof(unsigned char),cudaMemcpyHostToDevice);

  double *sunDirectCuda;
  cudaMalloc((void**)&sunDirectCuda, 3*sizeof(double));
  cudaMemcpy(sunDirectCuda,sunDirect,3*sizeof(double),cudaMemcpyHostToDevice);

  int *sunColorCuda;
  cudaMalloc((void**)&sunColorCuda, 3*sizeof(int));
  cudaMemcpy(sunColorCuda,sunColor,3*sizeof(int),cudaMemcpyHostToDevice);

  for( int sample = 0; sample++ < numSamples; )
  {
    for( int ia = 0, aliasIndex = 0; ia < numAlias; ++ia )
    {
      for( int ja = 0; ja < numAlias; ++ja, ++aliasIndex )
      {
        double ia1 = (double)ia / (double)numAlias;
        double ja1 = (double)ja / (double)numAlias;

        generateMandelboxPoint<<<numBlocks, blockSize>>>( start, stride, aliasIndex, numAlias, bx, by, xSize, ySize, size, sample, startDataCuda, imageDataCuda, backgroundDataCuda, directLighting, value, reflectance, fractalType, minSize, ia1, ja1, maxDepth, backgroundWidth, backgroundHeight, sunDirectCuda, sunColorCuda );

        cudaDeviceSynchronize();

        if( start == 0 ) generatingMandelbrot->Increment();
      }
    }
  }

  cudaMemcpy(imageData2,imageDataCuda,numToCompute*3*sizeof(int),cudaMemcpyDeviceToHost);

  for( int index = 0; index < numToCompute; ++index )
  {
    int imageIndex = index*stride + start;
    if( imageIndex < xSize*ySize )
    {
      for( int k1 = 0; k1 < 3; ++k1 )
      {
        imageData2[3*index+k1] /= (numSamples*numAlias*numAlias);
        imageData[3*imageIndex+k1] = imageData2[3*index+k1];
      }
    }
  }
}

void Mandelbox( string outputName, string backgroundName, bool directLighting, int fractalType, int numSamples, int numAlias, int maxDepth, double minIter, double value, double reflectance, int imageSize, int bx, int by, int xSize, int ySize )
{
  ProgressBar *generatingMandelbrot = new ProgressBar( numSamples*numAlias*numAlias, "Generating mandelbox" );

  generatingMandelbrot->Progressed(0);

  unsigned char *imageData = (unsigned char *)malloc(xSize*ySize*3*sizeof(unsigned char));

  int numGpus = 1;

  int threads = numGpus;

  int backgroundWidth, backgroundHeight;

  ifstream data( backgroundName, ios::binary );
  data.read( (char *)&backgroundWidth, sizeof(int) );
  data.read( (char *)&backgroundHeight, sizeof(int) );

  unsigned char *backgroundData = (unsigned char *)malloc(4*backgroundWidth*backgroundHeight*sizeof(unsigned char));

  data.read( (char *)backgroundData, sizeof(unsigned char)*4*backgroundWidth*backgroundHeight );

  double *sunDirect = (double *)malloc(3*sizeof(double));
  int *sunColor = (int *)malloc(3*sizeof(int));

  data.read( (char *)sunDirect, sizeof(double)*3 );
  data.read( (char *)sunColor, sizeof(int)*3 );

  data.close();

  future< void > ret[threads];

  for( int k = 0; k < threads; ++k )
  {
    ret[k] = async( launch::async, &MandelboxThread, k, threads, bx, by, xSize, ySize, imageSize, imageData, backgroundData, directLighting, value, minIter, reflectance, numSamples, numAlias, maxDepth, fractalType, backgroundWidth, backgroundHeight, sunDirect, sunColor, generatingMandelbrot ); 
  }

  // Wait for threads to finish
  for( int k = 0; k < threads; ++k )
  {
    ret[k].get();
  }

  generatingMandelbrot->Finish();

  VImage::new_from_memory( imageData, xSize*ySize*3, xSize, ySize, 3, VIPS_FORMAT_UCHAR ).vipssave((char *)outputName.c_str());
}

void RunMandelbox( string outputName, string backgroundName, bool directLighting, int fractalType, int numSamples, int numAlias, int maxDepth, double minIter, double value, double reflectance, int imageSize, int bx, int by, int xSize, int ySize )
{
  Mandelbox( outputName, backgroundName, directLighting, fractalType, numSamples, numAlias, maxDepth, minIter, value, reflectance, imageSize, bx, by, xSize, ySize );
}