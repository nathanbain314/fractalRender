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

void RunMandelbox( string outputName, string backgroundName, string gradientName, bool directLighting, int fractalType, int numSamples, int numAlias, int maxDepth, float minIter, float value, float color, float reflectance, int imageSize, int bx, int by, int xSize, int ySize, float kd, float ks, float ka, float alpha, int maxSteps, bool dontPathTrace );