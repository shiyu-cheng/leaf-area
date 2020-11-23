#pragma once

#include <iostream>
#include <fstream>
#include <cstdio>
using namespace std;
#include <omp.h>
#include <ctime>
#include "Eigen/Dense"



using namespace Eigen;



struct angle_extent {
	float azimuth_min;
	float azimuth_max;
	float zenith_min;
	float zenith_max;
};