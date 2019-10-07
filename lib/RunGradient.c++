#include "Gradient.h"

int main( int argc, char **argv )
{
  if( VIPS_INIT( argv[0] ) ) return( -1 );

  RunGradient( string(argv[1]), string(argv[2]) );

  return 0;
}