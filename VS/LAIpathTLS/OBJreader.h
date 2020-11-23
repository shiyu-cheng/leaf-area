#pragma once
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif // !1

#ifndef OBJREADER_H
#define OBJREADER_H

using namespace std;

#include "commdef.h"

class OBJreader
{
public:
	bool OBJreader::read();
	OBJreader() {};
	OBJreader(const wchar_t* _filename) :filename(_filename) {};

	Vector3d get_triangle_vective(unsigned index_triangle, unsigned index_vectice);
	int num_triangles();

	void set_scanner_position (Vector3d _scanner_position);
	void shift_to_observational_frame();
	bool is_scanner_too_near();
	angle_extent get_angle_extent();
	

private:
	const wchar_t* filename;
	angle_extent angle_extent_envelope_;

	MatrixXd vertices;
	MatrixXi triangles;

	unsigned n_vertices;
	unsigned n_triangles;

	Vector3d scanner_position_;
	bool is_shifted = false;

};
#endif
