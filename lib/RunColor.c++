#include "Color.h"

int main( int argc, char **argv )
{
  if( VIPS_INIT( argv[0] ) ) return( -1 );

  RunColor( string(argv[1]), string(argv[2]) );

  return 0;
}