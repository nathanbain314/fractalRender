#include <iostream>
#include <vector>
#include <cmath>
#include <climits>
#include <dirent.h>
#include <unistd.h>
#include <future>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <random>
#include <set>
#include <tuple>
#include <time.h>

#include <vips/vips8>

#include "progress_bar.hpp"

using namespace vips;
using namespace std;

void RunMandelbox( string outputName, string backgroundName, bool directLighting, int fractalType, int numSamples, int numAlias, int maxDepth, double minIter, double value, double reflectance, int imageSize, int bx, int by, int xSize, int ySize );