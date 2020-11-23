#pragma once

#include "PTXreader.h"
#include "OBJreader.h"
#include "LAIPath.h"


using namespace Eigen;

double CalcFAVD(double gap_probability_within_canopy, VectorXd path_lengths_whinin_canopy, double _G = 0.5);
double CalcMaxPathLength(double gap_probability_within_canopy, VectorXd path_lengths_whinin_canopy, double _G = 0.5);
double CalcFAVD2(double gap_probability_within_canopy, VectorXd ray_lengths_whinin_canopy, double _G = 0.5);
double CalcFAVDmltLmax(double gap_probability_within_canopy, VectorXd path_lengths_whinin_canopy, double _G = 0.5);

ArrayXf CalcPathLengths(OBJreader envelope, PTXreader ray, ArrayXXb & points_outside_canopy, ArrayXXb zenith_mask);
ArrayXf CalcRayLengths(OBJreader envelope, PTXreader ray, ArrayXXb& points_outside_canopy, ArrayXXb zenith_mask);

float CalculatePathLength(OBJreader & envelope, const Vector3f& dir, bool &point_outside_canopy);
float CalculateRayLength(OBJreader& envelope, const Vector3f& dir, bool& point_outside_canopy);

double CalcRayArea(OBJreader& envelope, PTXreader& ray, Index row, Index next_row, double& height);
double CalcRayArea2(ArrayXXb within_canopy_image, OBJreader& envelope, PTXreader& ray, Index row, Index next_row, int number_of_divisions, bool use_gap_p,
	                Map<ArrayXXf> path_lengths_image, ArrayXXf gap_p, ArrayXXb gap_p_valid, double& height, double& volumn);

bool IntersectTriangle(const Vector3f& orig, const Vector3f& dir,
	Vector3f & v0, Vector3f& v1, Vector3f& v2,
	float* t, float* u, float* v);

bool IntersectTriangle(const Vector3d& orig, const Vector3d& dir, 
	Vector3d & v0, Vector3d& v1, Vector3d& v2, 
	double* t, double* u, double* v);
