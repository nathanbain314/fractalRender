#include "Zoomable.h"
#include "progress_bar.hpp"

void RunZoomable( int outputWidth, int outputHeight, string outputImage )
{
  int level = (int)ceil(log2( max(outputWidth,outputHeight) ) );

  ofstream dzi_file;
  dzi_file.open(string(outputImage).append(".dzi").c_str());
  dzi_file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
  dzi_file << "<Image xmlns=\"http://schemas.microsoft.com/deepzoom/2008\" Format=\"jpeg\" Overlap=\"0\" TileSize=\"256\">" << endl;
  dzi_file << "    <Size Height=\"" << outputHeight << "\" Width=\"" << outputWidth << "\"/>" << endl;
  dzi_file << "</Image>" << endl;
  dzi_file.close();

  int numLevels = 0;
  for( int o = ceil((double)outputWidth/256.0); o > 1; o = ceil((double)o/2.0) ) numLevels += o;

  ProgressBar *lowerLevels = new ProgressBar(numLevels, "Building lower levels");

  outputWidth = (int)ceil((double)outputWidth/ 128.0 );
  outputHeight = (int)ceil((double)outputHeight/ 128.0 );

  for( ; level > 0; --level )
  {
    outputWidth = (int)ceil((double)outputWidth/ 2.0 );
    outputHeight = (int)ceil((double)outputHeight/ 2.0 );

    string current = string(outputImage).append("_files/"+to_string(level-1)+"/");
    string upper = string(outputImage).append("_files/"+to_string(level)+"/");

    g_mkdir((char *)current.c_str(), 0777);

    for(int i = 1; i < outputWidth; i+=2)
    {
      lowerLevels->Increment();

      for(int j = 1; j < outputHeight; j+=2)
      {
        (VImage::jpegload((char *)string(upper).append(to_string(i-1)+"_"+to_string(j-1)+".jpeg").c_str(), VImage::option()->set( "shrink", 2)).
        join(VImage::jpegload((char *)string(upper).append(to_string(i)+"_"+to_string(j-1)+".jpeg").c_str(), VImage::option()->set( "shrink", 2)),VIPS_DIRECTION_HORIZONTAL)).
        join(VImage::jpegload((char *)string(upper).append(to_string(i-1)+"_"+to_string(j)+".jpeg").c_str(), VImage::option()->set( "shrink", 2)).
        join(VImage::jpegload((char *)string(upper).append(to_string(i)+"_"+to_string(j)+".jpeg").c_str(), VImage::option()->set( "shrink", 2)),VIPS_DIRECTION_HORIZONTAL),VIPS_DIRECTION_VERTICAL).
        jpegsave((char *)string(current).append(to_string(i>>1)+"_"+to_string(j>>1)+".jpeg").c_str() );
      }
    }
    if(outputWidth%2 == 1)
    {
      for(int j = 1; j < outputHeight; j+=2)
      {
        (VImage::jpegload((char *)string(upper).append(to_string(outputWidth-1)+"_"+to_string(j-1)+".jpeg").c_str(), VImage::option()->set( "shrink", 2)).
        join(VImage::jpegload((char *)string(upper).append(to_string(outputWidth-1)+"_"+to_string(j)+".jpeg").c_str(), VImage::option()->set( "shrink", 2)),VIPS_DIRECTION_VERTICAL)).
        jpegsave((char *)string(current).append(to_string(outputWidth>>1)+"_"+to_string(j>>1)+".jpeg").c_str(), VImage::option()->set( "optimize_coding", true )->set( "strip", true ) );
      }
    }
    if(outputHeight%2 == 1)
    {
      for(int j = 1; j < outputWidth; j+=2)
      {
        (VImage::jpegload((char *)string(upper).append(to_string(j-1)+"_"+to_string(outputHeight-1)+".jpeg").c_str(), VImage::option()->set( "shrink", 2)).
        join(VImage::jpegload((char *)string(upper).append(to_string(j)+"_"+to_string(outputHeight-1)+".jpeg").c_str(), VImage::option()->set( "shrink", 2)),VIPS_DIRECTION_HORIZONTAL)).
        jpegsave((char *)string(current).append(to_string(j>>1)+"_"+to_string(outputHeight>>1)+".jpeg").c_str(), VImage::option()->set( "optimize_coding", true )->set( "strip", true ) );
      }
    }
    if(outputWidth%2 == 1 && outputHeight%2 == 1)
    {
      VImage::jpegload((char *)string(upper).append(to_string(outputWidth-1)+"_"+to_string(outputHeight-1)+".jpeg").c_str()).resize(0.5).
      jpegsave((char *)string(current).append(to_string(outputWidth>>1)+"_"+to_string(outputHeight>>1)+".jpeg").c_str(), VImage::option()->set( "optimize_coding", true )->set( "strip", true ) );  
    }
  }

  lowerLevels->Finish();

  // Generate html file to view deep zoom image
  ofstream htmlFile(string(outputImage).append(".html").c_str());
  htmlFile << "<!DOCTYPE html>\n<html>\n<head><script src=\"js/openseadragon.min.js\"></script></head>\n<body>\n<style>\nhtml,\nbody,\n#collage\n{\nposition: fixed;\nleft: 0;\ntop: 0;\nwidth: 100%;\nheight: 100%;\n}\n</style>\n\n<div id=\"collage\"></div>\n\n<script>\nvar viewer = OpenSeadragon({\nid: 'collage',\nprefixUrl: 'icons/',\ntileSources:   \"" + outputImage + ".dzi\",\nminZoomImageRatio: 0,\nmaxZoomImageRatio: 1\n});\n</script>\n</body>\n</html>";
  htmlFile.close();
}