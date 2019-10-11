#include "Mandelbox.h"
#include "Zoomable.h"
#include "progress_bar.hpp"

int equiProject( float* d, int backgroundWidth, int backgroundHeight )
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

  j = min( j, backgroundWidth - 1 );
  i = min( i, backgroundHeight - 1 );

  j = max( j, 0 );
  i = max( i, 0 );
  
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





// Mandelbulb
float MandelbulbDE( float* z, float n, float &colorIteration )
{
  z[2] *= -1;

  float z2[3];

  z2[0] = z[0];
  z2[1] = z[1];
  z2[2] = z[2];

  int maxIter = 100;

  float dr = 1.0;

  float r;

  int k;

  for( k = 0; k < maxIter; ++k )
  {
    r = length(z2);

    if( r > 4 ) break;

    float t = acos(z2[2]/r);
    float p = atan2(z2[1],z2[0]);

    dr = pow( r, n-1.0 ) * n * dr + 1.0;

    float r1 = pow(r,n);

    z2[0] = r1*sin(t*n)*cos(p*n) + z[0];
    z2[1] = r1*sin(t*n)*sin(p*n) + z[1];
    z2[2] = r1*cos(t*n) + z[2];
  }

  colorIteration = (float)k - log2f(log2f(r));

  z[2] *= -1;

  return 0.25 * log(r) * r / dr;
}

// Mandelbox
float MandelboxDE( float* z1, float scale, float &colorIteration )
{
  float sideLength = (scale < -1 ? 4 : 4*(scale+1)/(scale-1))*1.25;

  int iter = 0;

  float DEfactor = scale;

  float fR2 = 1.0;
  float mR2 = 0.25;

  float x = z1[0]*sideLength*0.5;
  float y = z1[1]*sideLength*0.5;
  float z = z1[2]*sideLength*0.5;

  float cx = x;
  float cy = y;
  float cz = z;

  float r = sqrt(x*x+y*y+z*z);

  for( ; iter < 100; ++iter )
  {    
    if( r > 1000 ) break;

    if (x > 1.0)
    x = 2.0 - x;
    else if (x < -1.0) x = -2.0 - x;
    if (y > 1.0)
    y = 2.0 - y;
    else if (y < -1.0) y = -2.0 - y;
    if (z > 1.0)
    z = 2.0 - z;
    else if (z < -1.0) z = -2.0 - z;

    float r2 = x*x + y*y + z*z;

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

    r = sqrt(x*x+y*y+z*z);
  }

  colorIteration = (float)iter - log2f(log2f(r));

  return r/abs(sideLength*DEfactor);
}


// Menger
float MengerDE( float* point, float n10, float &colorIteration )
{
  int n = 6;

  float x=point[0]*1.2, y=point[1]*1.2, z=point[2]*1.2;
  x=x*0.5+0.5;y=y*0.5+0.5;z=z*0.5+0.5;

  float xx=abs(x-0.5)-0.5, yy=abs(y-0.5)-0.5, zz=abs(z-0.5)-0.5;
  float d1=max(xx,max(yy,zz));
  float d=d1;
  float p=1.0;
  for (int i=1; i<=n; ++i) {
    float xa = fmod(3.0*x*p,3.0);
    float ya = fmod(3.0*y*p,3.0);
    float za = fmod(3.0*z*p,3.0);
    p*=3.0;

    float xx=0.5-abs(xa-1.5), yy=0.5-abs(ya-1.5), zz=0.5-abs(za-1.5);
    d1=min(max(xx,zz),min(max(xx,yy),max(yy,zz))) / p;

    d=max(d,d1);
  }

  colorIteration = 0;

  return d;
}


// Koch
void fold( float &x, float &y, float angle )
{
  float a1 = cos( -angle );
  float a2 = sin( -angle );

  float d = 2.0 * min( 0.0f, a1 * x + a2 * y );

  x -= d * a1;
  y -= d * a2;
}

void trifold( float *p, float angle )
{
  float pi = 3.14159265358979;
  fold( p[0], p[1], pi/3.0 - cos(angle)/10.0 );
  fold( p[0], p[1], -pi/3.0 );
  fold( p[1], p[2], -pi/6.0 + sin(angle)/2.0 );
  fold( p[1], p[2], pi/6.0 );
}

void tricurve( float *p, float angle )
{
  for( int i = 0; i < 7; ++i )
  {
    p[0] = 2 * p[0] - 2.6;
    p[1] *= 2;
    p[2] *= 2;

    trifold( p, angle );
  }
}

float kochDE( float* z, float angle, float &colorIteration )
{
  float p[3] = {z[0],z[1],z[2]};

  p[0] += 1.5;

  tricurve( p, angle );

  colorIteration = 0;

  return 0.004 * length( p ) - 0.01;
}

float de( float* z, float value, int fractalType, float &colorIteration )
{
  switch( fractalType )
  {
    case 0:
      return MandelbulbDE( z, value, colorIteration );
    case 1:
      return MandelboxDE( z, value, colorIteration );
    case 2:
      return MengerDE( z, value, colorIteration );
    case 3:
      return kochDE( z, value, colorIteration );
  }

  return MandelbulbDE( z, value, colorIteration );
}



void normalize( float* n )
{
  float l = length( n );

  for( int k1 = 0; k1 < 3; ++k1 )
  {
    n[k1] /= l;
  }
}

void computeNormal( float* p, float* d, float value, int fractalType, float *n )
{
  float colorIteration;
  float t[3];

  for( int k1 = 0; k1 < 3; ++k1 )
  {
    t[0] = p[0];
    t[1] = p[1];
    t[2] = p[2];

    t[k1] = p[k1];// + d[k1];

    float nt = de(t,value,fractalType,colorIteration);

    t[k1] = p[k1] - d[k1];
    
    nt -= de(t,value,fractalType,colorIteration);

    n[k1] = nt;
  }

  normalize( n );
}

void coneSample( float* n, float* seed, float extent )
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

  seed[0] -= 1;
  seed[1] += 1;

  float r1 = sin((seed[0] * 12.9898 + seed[1] * 78.233))*43758.5453;
  r1 = r1 - floor(r1);

  float r2 = cos((seed[0] * 4.898 + seed[1] * 7.23))*23421.631;
  r2 = r2 - floor(r2);

  r1 = 1.0 - r1*extent;

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

void randomFromNormal( float* n, float* seed )
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

  seed[0] -= 1;
  seed[1] += 1;

  float r1 = sin((seed[0] * 12.9898 + seed[1] * 78.233))*43758.5453;
  r1 = r1 - floor(r1);

  float r2 = cos((seed[0] * 4.898 + seed[1] * 7.23))*23421.631;
  r2 = r2 - floor(r2);

  r1 = sqrt(r1);

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

void reflect( float* d, float* n )
{
  float v = 2.0 * dotProduct( d, n );
  for( int k1 = 0; k1 < 3; ++k1 )
  {
    d[k1] = d[k1] - v * n[k1];
  }
  normalize(d);
}

bool findHit( float *p, float* d, float value, int fractalType, float minSize, float* output, float &colorIteration, int &numSteps, int maxSteps )
{
  numSteps = 0;

  while( abs(p[0]) < 15 && abs(p[1]) < 15 && abs(p[2]) < 15 && ++numSteps < maxSteps )
  {
    float iter = de( p, value, fractalType, colorIteration );

    if( iter < minSize )
    {
      normalize(d);

      for( int k1 = 0; k1 < 3; ++k1 )
      {
        d[k1] = minSize;
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

void generateMandelboxPoint( int start, int stride, size_t idx1, int aliasIndex, int numAlias, int bx, int by, int xSize, int ySize, int size, int sample, float *startData, int *imageData, unsigned char * backgroundData, unsigned char * gradientData, bool directLighting, float value, float colorMultiplier, float reflectance, int fractalType, float minSize, float ia, float ja, int maxDepth, int backgroundWidth, int backgroundHeight, float* sunDirect2, int* sunColor, int maxSteps, bool dontPathTrace )
{
  size_t idx2 = 8 * ( idx1 * numAlias * numAlias + aliasIndex );
  size_t idx = idx1*stride + start;
  int yPixel = by + idx / xSize;
  int xPixel = bx + idx % xSize;

  float num1 = 1.732050807569;
  float num2 = 2.121320343560;

  if( idx < xSize*ySize || ( sample > 1 && startData[idx2] > 1.0 ) )
  {
    float i = (float)yPixel + ia - (float) size / 2.0;
    float j = (float)xPixel + ja - (float) size / 2.0;

    float x = 10.0 + num1*(float)i/(float)size + num2*(float)j/(float)size;
    float y = 10.0 - num1*(float)i/(float)size;
    float z = 10.0 + num1*(float)i/(float)size - num2*(float)j/(float)size;

    float minOffset = min( x - 1.5, min( y - 1.5, z - 1.5 ) );

    x -= minOffset;
    y -= minOffset;
    z -= minOffset;

    float p[3] = {x,y,z};
    float d[3] = { -0.57735026919, -0.57735026919, -0.57735026919 };

    float n[3], n2[3];

    float colorIteration;

    int numSteps;

    if( sample < 2 )
    {
      if( !aliasIndex )
      {
        imageData[3*idx1+0] = 0;
        imageData[3*idx1+1] = 0;
        imageData[3*idx1+2] = 0;
      }

      if( findHit( p, d, value, fractalType, minSize, n, colorIteration, numSteps, maxSteps ) )
      {
        startData[idx2] = 2.0;

        startData[idx2+1] = p[0];
        startData[idx2+2] = p[1];
        startData[idx2+3] = p[2];

        startData[idx2+4] = n[0];
        startData[idx2+5] = n[1];
        startData[idx2+6] = n[2];

        startData[idx2+7] = colorIteration;

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

      colorIteration = startData[idx2+7];
    }



    if( startData[idx2] > 1.0 )
    {
      float seed[2] = { (float)(xPixel*sample)+ja, (float)(yPixel*sample)+ia };

      float color[3] = {1.0,1.0,1.0};

      float direct[3] = {0.0,0.0,0.0};

      for( int depth = 0; depth <= maxDepth; ++depth )
      {
        if( !depth || findHit( p, d, value, fractalType, minSize, n, colorIteration, numSteps, maxSteps ) )
        {
          if( depth == maxDepth ) break;

          if( dontPathTrace )
          {
            float ashadow = 1.0 - double( numSteps ) / double( maxSteps );

            if( colorMultiplier > 0.0 )
            {
              while( colorIteration < 1024 ) colorIteration += 1024;

              colorIteration *= colorMultiplier;

              int colorIndex = colorIteration;

              colorIndex %= 1024;

              colorIteration -= floor( colorIteration );

              for( int k1 = 0; k1 < 3; ++k1 )
              {
                imageData[3*idx1+k1] += min( 255.0f, ashadow * interpolate(gradientData[3*colorIndex+k1],gradientData[3*((colorIndex+1)%1024)+k1],colorIteration) );
              }
            }
            else
            {
              for( int k1 = 0; k1 < 3; ++k1 )
              {
                imageData[3*idx1+k1] += ashadow * 255.0f;
              }
            }

            return;
          }

          seed[0] -= 1;
          seed[1] += 1;

          float r1 = sin((seed[0] * 12.9898 + seed[1] * 78.233))*43758.5453;
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


          if( colorMultiplier > 0.0 )
          {
            while( colorIteration < 1024 ) colorIteration += 1024;

            colorIteration *= colorMultiplier;

            int colorIndex = colorIteration;

            colorIndex %= 1024;

            colorIteration -= floor( colorIteration );

            for( int k1 = 0; k1 < 3; ++k1 )
            {
              color[k1] *= interpolate(gradientData[3*colorIndex+k1],gradientData[3*((colorIndex+1)%1024)+k1],colorIteration) / 255.0;
            }
          }


          if( directLighting )
          {
            float sunDirect[3] = { sunDirect2[0], sunDirect2[1], sunDirect2[2] };

            coneSample( sunDirect, seed, 0.00001 );

            float sunLight = dotProduct( sunDirect, n );

            float p2[3] = {p[0],p[1],p[2]};

            if( sunLight > 0 && !findHit( p2, sunDirect, value, fractalType, minSize, n2, colorIteration, numSteps, maxSteps ) )
            {
              for( int k1 = 0; k1 < 3; ++k1 )
              {
                direct[k1] += color[k1] * sunLight;
              }
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
            imageData[3*idx1+k1] += min(255.0f,sunColor[k1]*direct[k1] + color[k1] * ldexp((float)backgroundData[4*index+k1],mantissa));
          }

          break;
        }
      }
    }
  }
}

void MandelboxBlock( int start, int stride, int bx, int by, int xSize, int ySize, int size, int * imageData2, float * startData, unsigned char * backgroundData, unsigned char * gradientData, bool directLighting, float value, float minSize, float color, float reflectance, int numSamples, int numAlias, int maxDepth, int fractalType, int backgroundWidth, int backgroundHeight, float* sunDirect, int* sunColor, int maxSteps, bool dontPathTrace, bool isDeepZoom, ProgressBar *generatingMandelbrot )
{
  int numToCompute = xSize*ySize/stride;
  
  for( int sample = 0; sample++ < numSamples; )
  {
    for( int ia = 0, aliasIndex = 0; ia < numAlias; ++ia )
    {
      for( int ja = 0; ja < numAlias; ++ja, ++aliasIndex )
      {
        float ia1 = (float)ia / (float)numAlias;
        float ja1 = (float)ja / (float)numAlias;

        for( int index = 0; index < numToCompute; ++index )
        {
          generateMandelboxPoint( start, stride, index, aliasIndex, numAlias, bx, by, xSize, ySize, size, sample, startData, imageData2, backgroundData, gradientData, directLighting, value, color, reflectance, fractalType, minSize, ia1, ja1, maxDepth, backgroundWidth, backgroundHeight, sunDirect, sunColor, maxSteps, dontPathTrace );
        }

        if( !isDeepZoom && start == 0 ) generatingMandelbrot->Increment();
      }
    }
  }
}

void MandelboxThread( int start, int stride, int bx, int by, int xSize, int ySize, int size, unsigned char * imageData, unsigned char * backgroundData, unsigned char * gradientData, bool directLighting, float value, float minSize, float color, float reflectance, int numSamples, int numAlias, int maxDepth, int fractalType, int backgroundWidth, int backgroundHeight, float* sunDirect, int* sunColor, int maxSteps, bool dontPathTrace, bool isDeepZoom, string outputName, ProgressBar *generatingMandelbrot )
{
  if( isDeepZoom )
  {
    int numAcross = size/256;
    int numImagesToCompute = numAcross * numAcross;
    numImagesToCompute /= stride;

    unsigned char *imageData3 = (unsigned char *)malloc(256*256*3*sizeof(unsigned char));

    int *imageData2 = (int*)malloc(256*256*3*sizeof(int));

    float *startData = (float*)malloc(256*256*8*numAlias*numAlias*sizeof(float));

    if( !start ) generatingMandelbrot->Progressed(0);

    for( int index1 = 0; index1 < numImagesToCompute; ++index1 )
    {
      int idx = index1*stride + start;
      int i = idx/numAcross;
      int j = idx%numAcross;
    
      MandelboxBlock( 0, 1, j*256, i*256, 256, 256, size, imageData2, startData, backgroundData, gradientData, directLighting, value, minSize, color, reflectance, numSamples, numAlias, maxDepth, fractalType, backgroundWidth, backgroundHeight, sunDirect, sunColor, maxSteps, dontPathTrace, isDeepZoom, generatingMandelbrot ); 

      for( int index = 0; index < 256*256; ++index )
      {
        for( int k1 = 0; k1 < 3; ++k1 )
        {
          imageData2[3*index+k1] /= (numSamples*numAlias*numAlias);
          imageData3[3*index+k1] = imageData2[3*index+k1];
        }
      }

      if( start == 0 ) generatingMandelbrot->Increment();

      VImage::new_from_memory( imageData3, 256*256*3, 256, 256, 3, VIPS_FORMAT_UCHAR ).vipssave((char *)string(outputName).append(to_string(j)+"_"+to_string(i)+".jpeg").c_str());
    }
  }
  else
  {
    int numToCompute = xSize*ySize/stride;

    int *imageData2 = (int*)malloc(numToCompute*3*sizeof(int));

    float *startData = (float*)malloc(numToCompute*8*numAlias*numAlias*sizeof(float));

    if( !start ) generatingMandelbrot->Progressed(0);

    MandelboxBlock( start, stride, bx, by, xSize, ySize, size, imageData2, startData, backgroundData, gradientData, directLighting, value, minSize, color, reflectance, numSamples, numAlias, maxDepth, fractalType, backgroundWidth, backgroundHeight, sunDirect, sunColor, maxSteps, dontPathTrace, isDeepZoom, generatingMandelbrot ); 

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
}

void Mandelbox( string outputName, string backgroundName, string gradientName, bool directLighting, int fractalType, int numSamples, int numAlias, int maxDepth, float minIter, float value, float color, float reflectance, int imageSize, int bx, int by, int xSize, int ySize, int maxSteps, bool dontPathTrace )
{
  bool isDeepZoom = (vips_foreign_find_save( outputName.c_str() ) == NULL);


  int threads = sysconf(_SC_NPROCESSORS_ONLN);

  ifstream gradientFile( gradientName, ios::binary );

  unsigned char *gradientData = (unsigned char *)malloc(3*1024*sizeof(unsigned char));

  gradientFile.read( (char *)gradientData, 3*1024*sizeof(unsigned char) );

  gradientFile.close();

  int backgroundWidth, backgroundHeight;

  ifstream data( backgroundName, ios::binary );
  data.read( (char *)&backgroundWidth, sizeof(int) );
  data.read( (char *)&backgroundHeight, sizeof(int) );

  unsigned char *backgroundData = (unsigned char *)malloc(4*backgroundWidth*backgroundHeight*sizeof(unsigned char));

  data.read( (char *)backgroundData, sizeof(unsigned char)*4*backgroundWidth*backgroundHeight );

  float *sunDirect = (float *)malloc(3*sizeof(float));
  int *sunColor = (int *)malloc(3*sizeof(int));

  data.read( (char *)sunDirect, sizeof(float)*3 );
  data.read( (char *)sunColor, sizeof(int)*3 );

  data.close();



  if( isDeepZoom )
  {
    string outputName2 = outputName;
    outputName.append("_files/");
    g_mkdir(outputName.c_str(), 0777);
    outputName.append(to_string(int(ceil(log2(imageSize))))+"/");
    g_mkdir(outputName.c_str(), 0777);

    int total = imageSize/256;

    ProgressBar *generatingMandelbrot = new ProgressBar( total*total/threads, "Generating mandelbox" );

    unsigned char *imageData;

    future< void > ret[threads];

    for( int k = 0; k < threads; ++k )
    {
      ret[k] = async( launch::async, &MandelboxThread, k, threads, bx, by, xSize, ySize, imageSize, imageData, backgroundData, gradientData, directLighting, value, minIter, color, reflectance, numSamples, numAlias, maxDepth, fractalType, backgroundWidth, backgroundHeight, sunDirect, sunColor, maxSteps, dontPathTrace, isDeepZoom, outputName, generatingMandelbrot ); 
    }

    // Wait for threads to finish
    for( int k = 0; k < threads; ++k )
    {
      ret[k].get();
    }

    generatingMandelbrot->Finish();

    RunZoomable( total, outputName2 );
  }
  else
  {
    ProgressBar *generatingMandelbrot = new ProgressBar( numSamples*numAlias*numAlias, "Generating mandelbox" );

    generatingMandelbrot->Progressed(0);

    unsigned char *imageData = (unsigned char *)malloc(xSize*ySize*3*sizeof(unsigned char));

    future< void > ret[threads];

    for( int k = 0; k < threads; ++k )
    {
      ret[k] = async( launch::async, &MandelboxThread, k, threads, bx, by, xSize, ySize, imageSize, imageData, backgroundData, gradientData, directLighting, value, minIter, color, reflectance, numSamples, numAlias, maxDepth, fractalType, backgroundWidth, backgroundHeight, sunDirect, sunColor, maxSteps, dontPathTrace, isDeepZoom, outputName, generatingMandelbrot ); 
    }

    // Wait for threads to finish
    for( int k = 0; k < threads; ++k )
    {
      ret[k].get();
    }

    generatingMandelbrot->Finish();

    VImage::new_from_memory( imageData, xSize*ySize*3, xSize, ySize, 3, VIPS_FORMAT_UCHAR ).vipssave((char *)outputName.c_str());
  }
}

void RunMandelbox( string outputName, string backgroundName, string gradientName, bool directLighting, int fractalType, int numSamples, int numAlias, int maxDepth, float minIter, float value, float color, float reflectance, int imageSize, int bx, int by, int xSize, int ySize, int maxSteps, bool dontPathTrace )
{
  Mandelbox( outputName, backgroundName, gradientName, directLighting, fractalType, numSamples, numAlias, maxDepth, minIter, value, color, reflectance, imageSize, bx, by, xSize, ySize, maxSteps, dontPathTrace );
}