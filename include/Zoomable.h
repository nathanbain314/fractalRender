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

void RunZoomable( int outputWidth, int outputHeight, string outputImage );