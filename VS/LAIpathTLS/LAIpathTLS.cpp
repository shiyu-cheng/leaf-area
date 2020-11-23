#include "LAIpathTLS.h"


double CalcFAVD(double gap_probability_within_canopy, VectorXd path_lengths_whinin_canopy, double _G)
{
	printf("Calculate FAVD...");
	clock_t launch = clock();
	gsl_histogram * gsl_hist_path = gsl_histogram_alloc(NUM_BINS);			// path length distribution

	double max_pathlength = Stat_hist(path_lengths_whinin_canopy.data(), (unsigned long)path_lengths_whinin_canopy.size(), gsl_hist_path);			// running statistics to get path length distribution
	/*START-----------------------11/24/2019-----------------------------*/
	//max_pathlength unknown    
	//max_pathlength = 1;
	// maybe is the max ray length?
	// double max_raylength = Stat_hist(path_lengths_whinin_canopy.data(), (unsigned long)path_lengths_whinin_canopy.size(), gsl_hist_path);
	/*END-----------------------11/24/2019-----------------------------*/

	FILE* file_debug = _wfopen(L"debug.txt", L"a");
	fprintf(file_debug, "Gap probability: %lf\n", gap_probability_within_canopy);
	fprintf(file_debug, "Max path length: %lf\n", max_pathlength);
	//cout << "max_pathlength:" << max_pathlength << endl;
	gsl_histogram_fprintf(file_debug, gsl_hist_path, "%.2f", "%.3f");
	fprintf(file_debug, "\n\n");
	fclose(file_debug);
	file_debug = 0;

	double FAVD = LAI_PATH_FAVD(gsl_hist_path, gap_probability_within_canopy, max_pathlength, _G);
	clock_t done = clock();
	printf(" %lf s\n", double(done - launch) / CLOCKS_PER_SEC);

	return FAVD;
}

double CalcFAVDmltLmax(double gap_probability_within_canopy, VectorXd path_lengths_whinin_canopy, double _G)
{
	printf("Calculate FAVD...");
	clock_t launch = clock();
	gsl_histogram* gsl_hist_path = gsl_histogram_alloc(NUM_BINS);			// path length distribution

	double max_pathlength = Stat_hist(path_lengths_whinin_canopy.data(), (unsigned long)path_lengths_whinin_canopy.size(), gsl_hist_path);			// running statistics to get path length distribution
	/*START-----------------------11/24/2019-----------------------------*/
	//max_pathlength unknown    
	max_pathlength = 1;
	/*END-----------------------11/24/2019-----------------------------*/

	FILE* file_debug = _wfopen(L"debug.txt", L"a");
	fprintf(file_debug, "Gap probability: %lf\n", gap_probability_within_canopy);
	fprintf(file_debug, "Max path length: %lf\n", max_pathlength);
	//cout << "max_pathlength:" << max_pathlength << endl;
	gsl_histogram_fprintf(file_debug, gsl_hist_path, "%.2f", "%.3f");
	fprintf(file_debug, "\n\n");
	fclose(file_debug);
	file_debug = 0;

	double FAVD = LAI_PATH_FAVD(gsl_hist_path, gap_probability_within_canopy, max_pathlength, _G);
	clock_t done = clock();
	printf(" %lf s\n", double(done - launch) / CLOCKS_PER_SEC);

	return FAVD;
}

double CalcFAVD2(double gap_probability_within_canopy, VectorXd ray_lengths_whinin_canopy, double _G)
{
	printf("Calculate FAVD2...");
	clock_t launch = clock();
	gsl_histogram* gsl_hist_path = gsl_histogram_alloc(NUM_BINS);			// path length distribution

	/*START-----------------------11/24/2019-----------------------------*/
	// maybe is the max ray length?
	double max_pathlength = Stat_hist(ray_lengths_whinin_canopy.data(), (unsigned long)ray_lengths_whinin_canopy.size(), gsl_hist_path);			// running statistics to get path length distribution
	//NONONONONONO!!!!
	/*END-----------------------11/24/2019-----------------------------*/

	FILE* file_debug = _wfopen(L"debug.txt", L"a");
	fprintf(file_debug, "Gap probability: %lf\n", gap_probability_within_canopy);
	fprintf(file_debug, "Max path length: %lf\n", max_pathlength);
	//cout << "max_pathlength:" << max_pathlength << endl;
	gsl_histogram_fprintf(file_debug, gsl_hist_path, "%.2f", "%.3f");
	fprintf(file_debug, "\n\n");
	fclose(file_debug);
	file_debug = 0;

	//double FAVD = LAI_PATH_FAVD(gsl_hist_path, gap_probability_within_canopy, max_pathlength, _G);
	clock_t done = clock();
	cout << "max_pathlength = " << max_pathlength  << endl;
	printf(" %lf s\n", double(done - launch) / CLOCKS_PER_SEC);

	return max_pathlength;
}

double CalcMaxPathLength(double gap_probability_within_canopy, VectorXd path_lengths_whinin_canopy, double _G)
{
	printf("Calculate MaxPathLength...");
	clock_t launch = clock();
	gsl_histogram* gsl_hist_path = gsl_histogram_alloc(NUM_BINS);			// path length distribution

	/*START-----------------------11/24/2019-----------------------------*/
	// maybe is the max ray length?
	double max_pathlength = Stat_hist(path_lengths_whinin_canopy.data(), (unsigned long)path_lengths_whinin_canopy.size(), gsl_hist_path);			// running statistics to get path length distribution
	/*END-----------------------11/24/2019-----------------------------*/

	FILE* file_debug = _wfopen(L"debug.txt", L"a");
	fprintf(file_debug, "Gap probability: %lf\n", gap_probability_within_canopy);
	fprintf(file_debug, "Max path length: %lf\n", max_pathlength);
	//cout << "max_pathlength:" << max_pathlength << endl;
	gsl_histogram_fprintf(file_debug, gsl_hist_path, "%.2f", "%.3f");
	fprintf(file_debug, "\n\n");
	fclose(file_debug);
	file_debug = 0;

	//double FAVD = LAI_PATH_FAVD(gsl_hist_path, gap_probability_within_canopy, max_pathlength, _G);
	clock_t done = clock();
	cout << "max_pathlength = " << max_pathlength << endl;
	printf(" %lf s\n", double(done - launch) / CLOCKS_PER_SEC);

	return max_pathlength;
}


ArrayXf CalcPathLengths(OBJreader envelope, PTXreader ray, ArrayXXb & points_outside_canopy, ArrayXXb zenith_mask)
{
	printf("Calculate Path Length...");
	clock_t launch = clock();

	ArrayXf path_distribution = VectorXf::Constant(ray.n_points(), -99);
	MatrixXf dir = ray.get_directions();
	points_outside_canopy = ray.get_gaps_image(); //input: isgaps, output: points_outside_canopy

#pragma omp parallel for 
	for (int i = 0; i < dir.rows(); ++i)
	{
		//bool point_outside_canopy = false;
		if (zenith_mask(i))
		{
			path_distribution(i) = CalculatePathLength(envelope, dir.row(i), points_outside_canopy(i));
			//printf("%d / %d\t%lf\n", i, dir.rows(), path_distribution(i));
		}
	}

	clock_t done = clock();
	printf(" %lf s\n", double(done - launch) / CLOCKS_PER_SEC);

	return path_distribution;

}

ArrayXf CalcRayLengths(OBJreader envelope, PTXreader ray, ArrayXXb& points_outside_canopy, ArrayXXb zenith_mask)
{
	printf("Calculate Path Length...");
	clock_t launch = clock();

	ArrayXf path_distribution = VectorXf::Constant(ray.n_points(), -99);
	MatrixXf dir = ray.get_directions();
	points_outside_canopy = ray.get_gaps_image(); //input: isgaps, output: points_outside_canopy

#pragma omp parallel for 
	for (int i = 0; i < dir.rows(); ++i)
	{
		//bool point_outside_canopy = false;
		if (zenith_mask(i))
		{
			/*START-----------------------11/24/2019-----------------------------*/
			//alter the path length to calculate the max path length of method  2 
			//path_distribution(i) = CalculatePathLength(envelope, dir.row(i), points_outside_canopy(i));
			path_distribution(i) = CalculateRayLength(envelope, dir.row(i), points_outside_canopy(i));
			/*END-----------------------11/24/2019-----------------------------*/
			//printf("%d / %d\t%lf\n", i, dir.rows(), path_distribution(i));
		}
	}

	clock_t done = clock();
	printf(" %lf s\n", double(done - launch) / CLOCKS_PER_SEC);

	return path_distribution;

}



float CalculatePathLength(OBJreader & envelope, const Vector3f & dir, bool & point_outside_canopy)
{
	int n_intersects = 0;
	float tmp_path = 0;
	//Vector3d intersect_coord;
	ArrayXf tt = ArrayXf::Zero(20);
	Array2f sign(-1, 1);
	ArrayXf tt_sign = sign.replicate(10, 1);
	//tt_sign << -1, 1, -1, 1, -1, 1, -1, 1, -1, 1;
	float d = dir.norm();
	bool isgap = point_outside_canopy;
	point_outside_canopy = false;

	for (int j = 0; j < envelope.num_triangles(); ++j)
	{
		//j = 336;
		/*float t, u, v;
		Vector3f v0 = envelope.get_triangle_vective(j, 0).cast<float>();
		Vector3f v1 = envelope.get_triangle_vective(j, 1).cast<float>();
		Vector3f v2 = envelope.get_triangle_vective(j, 2).cast<float>();*/
	//	if (IntersectTriangle(ray.get_scanner_position().cast<float>(), dir, v0, v1, v2, &t, &u, &v))

		float t, u, v;
		Vector3f v0 = envelope.get_triangle_vective(j, 0).cast<float>();
		Vector3f v1 = envelope.get_triangle_vective(j, 1).cast<float>();
		Vector3f v2 = envelope.get_triangle_vective(j, 2).cast<float>();
		if (IntersectTriangle(Vector3f::Zero(3), dir, v0, v1, v2, &t, &u, &v))  // Vectives shifted to observational frame of reference
		{
			// IntersectTriangle：Determine whether a ray intersect with a triangle
			//intersect_coord = (1 - u - v)*v0 + u * v1 + v * v2;
			tt(n_intersects++) = t * d;    //向量叉乘，四边形面积，原点到交点的距离；d：表达式向量的长度
			//https://www.cnblogs.com/graphics/archive/2010/08/09/1795348.html

			for (int k = 0; k < n_intersects - 1; ++k) 
			{
				if (abs(tt(n_intersects) - tt(n_intersects-1)) < 0.0001)
				{
					--n_intersects;	//remove duplicate records
				}
			}
		}

	}
	
	
	if (n_intersects == 0)
	{
		tmp_path = -1;
	} 
	else if (n_intersects == 1)
	{
		tmp_path = 0;
	}
	else 
	{ 
		std::sort(tt.data(), tt.data() + n_intersects);
		if (n_intersects % 2 == 0)
		{
			tmp_path = (tt * tt_sign).sum();   //tt_sign 1，-1，1，-1，相当于两两做差
			//tmp_path = -n_intersects;
		}
		else
		{
			tmp_path = tt(n_intersects - 1 ) - tt(0);
		}

		if (d > tt(n_intersects - 1) + 0.2)  // point pass through the canopy
		{
			point_outside_canopy = true;
		}
		else if (d < tt(0) - 0.2)	// point do not reach the canopy
		{
			tmp_path = -2;
		}
	}
		


	//switch (n_intersects)
	//{
	//case 0:
	//	tmp_path = -1;
	//	break;
	//case 1:
	//	tmp_path = 0;
	//	break;
	//case 2:
	//	tmp_path = tt(1) - tt(0);

	//	//if (!isgap && ((tt(0)-d)*(tt(1) - d) > 0))
	//	if (!isgap && d < 9998 && d > tt(1) + 0.1)
	//	{
	//		point_outside_canopy = true;
	//	}

	//	break;
	//default:
	//	if (n_intersects % 2 == 0)
	//	{
	//		tmp_path = (tt * tt_sign).sum();
	//
	//	} 
	//	else
	//	{
	//		tmp_path = tt(n_intersects) - tt(0);

	//	}

	//	if (!isgap && d < 9998 && d > tt(n_intersects - 1) + 0.1)
	//	{
	//		point_outside_canopy = true;
	//	}


	//	//tmp_path = -n_intersects;
	//	//printf("Warning: %d points of intersections\n", n_intersects);
	//	break;
	//}
	return tmp_path;
}



/*START-----------------------11/24/2019-----------------------------*/
/*Calculate ray length*/
float CalculateRayLength(OBJreader& envelope, const Vector3f& dir, bool& point_outside_canopy)
{
	int n_intersects = 0;
	float tmp_path = 0;
	//Vector3d intersect_coord;
	ArrayXf tt = ArrayXf::Zero(20);
	Array2f sign(-1, 1);
	ArrayXf tt_sign = sign.replicate(10, 1);
	//tt_sign << -1, 1, -1, 1, -1, 1, -1, 1, -1, 1;
	float d = dir.norm();
	bool isgap = point_outside_canopy;
	point_outside_canopy = false;

	for (int j = 0; j < envelope.num_triangles(); ++j)
	{
		//j = 336;
		/*float t, u, v;
		Vector3f v0 = envelope.get_triangle_vective(j, 0).cast<float>();
		Vector3f v1 = envelope.get_triangle_vective(j, 1).cast<float>();
		Vector3f v2 = envelope.get_triangle_vective(j, 2).cast<float>();*/
		//	if (IntersectTriangle(ray.get_scanner_position().cast<float>(), dir, v0, v1, v2, &t, &u, &v))

		float t, u, v;
		Vector3f v0 = envelope.get_triangle_vective(j, 0).cast<float>();
		Vector3f v1 = envelope.get_triangle_vective(j, 1).cast<float>();
		Vector3f v2 = envelope.get_triangle_vective(j, 2).cast<float>();
		if (IntersectTriangle(Vector3f::Zero(3), dir, v0, v1, v2, &t, &u, &v))  // Vectives shifted to observational frame of reference
		{
			// IntersectTriangle：Determine whether a ray intersect with a triangle
			//intersect_coord = (1 - u - v)*v0 + u * v1 + v * v2;
			tt(n_intersects++) = t * d;    //向量叉乘，四边形面积，原点到交点的距离；d：表达式向量的长度
			//https://www.cnblogs.com/graphics/archive/2010/08/09/1795348.html

			for (int k = 0; k < n_intersects - 1; ++k)
			{
				if (abs(tt(n_intersects) - tt(n_intersects - 1)) < 0.0001)
				{
					--n_intersects;	//remove duplicate records
				}
			}
		}

	}


	if (n_intersects == 0)
	{
		tmp_path = -1;
	}
	else if (n_intersects == 1)
	{
		tmp_path = 0;
	}
	else
	{
		std::sort(tt.data(), tt.data() + n_intersects);
		
		tmp_path = 0.5 * (tt(n_intersects - 1) + tt(0));
		//tmp_path = tt(n_intersects - 1);


		if (d > tt(n_intersects - 1) + 0.2)  // point pass through the canopy
		{
			point_outside_canopy = true;
		}
		else if (d < tt(0) - 0.2)	// point do not reach the canopy
		{
			tmp_path = -2;
		}
	}
	return tmp_path;
}



// Determine whether a ray intersect with a triangle
// Parameters
// orig: origin of the ray
// dir: direction of the ray
// v0, v1, v2: vertices of triangle
// t(out): weight of the intersection for the ray
// u(out), v(out): barycentric coordinate of intersection

bool IntersectTriangle(const Vector3f& orig, const Vector3f& dir,
	Vector3f & v0, Vector3f& v1, Vector3f& v2,
	float* t, float* u, float* v)
{
	// E1
	Vector3f E1 = v1 - v0;

	// E2
	Vector3f E2 = v2 - v0;

	// P
	Vector3f P = dir.cross(E2);

	// determinant
	float det = E1.dot(P);

	// keep det > 0, modify T accordingly
	Vector3f T;
	if (det > 0)
	{
		T = orig - v0;
	}
	else
	{
		T = v0 - orig;
		det = -det;
	}

	// If determinant is near zero, ray lies in plane of triangle
	if (det < 0.0001f)
		return false;

	// Calculate u and make sure u <= 1
	*u = T.dot(P);
	if (*u < 0.0f || *u > det)
		return false;

	// Q
	Vector3f Q = T.cross(E1);

	// Calculate v and make sure u + v <= 1
	*v = dir.dot(Q);
	if (*v < 0.0f || *u + *v > det)
		return false;

	// Calculate t, scale parameters, ray intersects triangle
	*t = E2.dot(Q);
	if (*t < 0.0f)
		return false;

	float fInvDet = 1.0f / det;
	*t *= fInvDet;
	*u *= fInvDet;
	*v *= fInvDet;

	return true;
}


bool IntersectTriangle(const Vector3d& orig, const Vector3d& dir,
	Vector3d & v0, Vector3d& v1, Vector3d& v2,
	double* t, double* u, double* v)
{
	// E1
	Vector3d E1 = v1 - v0;

	// E2
	Vector3d E2 = v2 - v0;

	// P
	Vector3d P = dir.cross(E2);

	// determinant
	double det = E1.dot(P);

	// keep det > 0, modify T accordingly
	Vector3d T;
	if (det > 0)
	{
		T = orig - v0;
	}
	else
	{
		T = v0 - orig;
		det = -det;
	}

	// If determinant is near zero, ray lies in plane of triangle
	if (det < 0.0001f)
		return false;

	// Calculate u and make sure u <= 1
	*u = T.dot(P);
	if (*u < 0.0f || *u > det)
		return false;

	// Q
	Vector3d Q = T.cross(E1);

	// Calculate v and make sure u + v <= 1
	*v = dir.dot(Q);
	if (*v < 0.0f || *u + *v > det)
		return false;

	// Calculate t, scale parameters, ray intersects triangle
	*t = E2.dot(Q);
	if (*t < 0.0f)
		return false;

	double fInvDet = 1.0f / det;
	*t *= fInvDet;
	*u *= fInvDet;
	*v *= fInvDet;

	
	return true;
}

double CalcRayArea(OBJreader& envelope, PTXreader& ray, Index row, Index next_row, double & height) {
	if (next_row >= ray.n_rows())
		return 0;
	int mid_col = 1 + ray.n_cols() / 2;
	MatrixXf dir = ray.get_directions();
	ArrayXf azimuths = ray.get_azimuths();
	ArrayXf zeniths = ray.get_zeniths();
	Index dir_index = 0;
	double start_loc = 0, end_loc = 0, start_azimuth = 0, end_azimuth = 0, start_zenith = 0, end_zenith = 0;
	double top_width, bottom_width;
	double length, l = 0;
	for (dir_index = mid_col; dir_index < ray.n_cols(); dir_index++) {
		Vector3f direction = dir.row(dir_index * ray.n_rows() + row);
		bool point_outside = false;
		length = CalculateRayLength(envelope, direction, point_outside);
		if (length > 0 && !point_outside) {
			start_zenith = zeniths(dir_index * ray.n_rows() + row);
			l = length;
			break;
		}
	}
	for (dir_index = 0; dir_index < ray.n_cols(); dir_index ++) {
		Vector3f direction = dir.row(dir_index * ray.n_rows() + row);
		bool point_outside = false;
		length = CalculateRayLength(envelope, direction, point_outside);
		if (length > 0 && !point_outside) {
			if (start_loc == 0) {
				start_azimuth = azimuths(dir_index * ray.n_rows() + row);
				start_loc = start_azimuth * length * sin(start_zenith);
			}
			end_azimuth = azimuths(dir_index * ray.n_rows() + row);
			end_loc = end_azimuth * length * sin(start_zenith);
		}
	}
	bottom_width = 2 * tan(abs(-end_azimuth + start_azimuth) / 2) * l * sin(start_zenith);

	start_loc = 0;
	end_loc = 0;
	for (dir_index = mid_col; dir_index < ray.n_cols(); dir_index++) {
		Vector3f direction = dir.row(dir_index * ray.n_rows() + next_row);
		bool point_outside = false;
		length = CalculateRayLength(envelope, direction, point_outside);
		if (length > 0 && !point_outside) {
			end_zenith = zeniths(dir_index * ray.n_rows() + next_row);
			l = length;
			break;
		}
	}
	for (dir_index = 0; dir_index < ray.n_cols(); dir_index++) {
		Vector3f direction = dir.row(dir_index * ray.n_rows() + next_row);
		bool point_outside = false;
		length = CalculateRayLength(envelope, direction, point_outside);
		if (length > 0 && !point_outside) {
			if (start_loc == 0) {
				start_azimuth = azimuths(dir_index * ray.n_rows() + next_row);
				start_loc = start_azimuth * length * sin(end_zenith);
			}
			end_azimuth = azimuths(dir_index * ray.n_rows() + next_row);
			end_loc = end_azimuth * length * sin(end_zenith);
		}
	}
	top_width = 2 * tan(abs(-end_azimuth + start_azimuth)/2) * l * sin(end_zenith);

	height = (start_zenith - end_zenith) * l / sin(end_zenith);
	// cout << start_azimuth << ' ' << end_azimuth << ' ' << l << ' ' << sin(end_zenith) << endl;
	// cout << top_width << ' ' << bottom_width << ' ' << l << ' ' << height << endl;
	if (end_zenith == 0) {	// 光打在树外，找不到了
		height = 0;
		return 0;
	}
	if (height < 0) {
		height = 0;
		return 0;
	}
	return height * (bottom_width + top_width) / 2;
}


double CalcRayArea2(ArrayXXb within_canopy_image, OBJreader& envelope, PTXreader& ray, Index row, Index next_row, int number_of_divisions, bool use_gap_p,
	               Map<ArrayXXf> path_lengths_image, ArrayXXf gap_p, ArrayXXb gap_p_valid, double & height, double & volumn) {
	if (next_row >= ray.n_rows())
		return 0;
	MatrixXf dir = ray.get_directions();
	Map<ArrayXXf> azimuths = ray.get_azimuths_image();
	Map<ArrayXXf> zeniths = ray.get_zeniths_image();
	double area = 0;
	bool has_canopy[100], has_canopy_next[100];
	Index col_step = ray.n_cols() / number_of_divisions;

	for (Index i = 0, start_col = 0; i < number_of_divisions; i++, start_col += col_step) {
		has_canopy[i] = false;
		Index dir_index = 0;
		for (dir_index = start_col; dir_index < start_col + col_step; dir_index++) {
			Vector3f direction = dir.row(dir_index * ray.n_rows() + row);
			bool point_outside = false;
			double length = CalculateRayLength(envelope, direction, point_outside);
			if (length > 0 && !point_outside) {
				has_canopy[i] = true;
				break;
			}
		}
	}
	for (Index i = 0, start_col = 0; i < number_of_divisions; i++, start_col += col_step) {
		has_canopy_next[i] = false;
		Index dir_index = 0;
		for (dir_index = start_col; dir_index < start_col + col_step; dir_index++) {
			Vector3f direction = dir.row(dir_index * ray.n_rows() + next_row);
			bool point_outside = false;
			double length = CalculateRayLength(envelope, direction, point_outside);
			if (length > 0 && !point_outside) {
				has_canopy_next[i] = true;
				break;
			}
		}
	}

	for (Index i = 0, start_col = 0; i < number_of_divisions; i ++, start_col += col_step) {
		int mid_col = start_col + col_step / 2;
		Index dir_index = 0;
		double start_azimuth = 0, end_azimuth = 0, start_zenith = 0, end_zenith = 0;
		double top_width, bottom_width;
		double length, l = 0;
		for (int delta = 0; delta < col_step / 2; delta++) {
			dir_index = mid_col + delta;
			Vector3f direction = dir.row(dir_index * ray.n_rows() + row);
			bool point_outside = false;
			length = CalculateRayLength(envelope, direction, point_outside);
			if (length > 0 && !point_outside) {
				start_zenith = zeniths(row, dir_index);
				l = length;
				break;
			}
			dir_index = mid_col - delta;
			direction = dir.row(dir_index * ray.n_rows() + row);
			length = CalculateRayLength(envelope, direction, point_outside);
			if (length > 0 && !point_outside) {
				start_zenith = zeniths(row, dir_index);
				l = length;
				break;
			}
		}
		if (i > 0 && has_canopy[i - 1]) {
			start_azimuth = azimuths(row, start_col);
		}
		if (i < number_of_divisions - 1 && has_canopy[i + 1]) {
			end_azimuth = azimuths(row, start_col + col_step - 1);
		}
		if (start_azimuth == 0 || end_azimuth == 0)
			for (dir_index = start_col; dir_index < start_col + col_step; dir_index++) {
				Vector3f direction = dir.row(dir_index * ray.n_rows() + row);
				bool point_outside = false;
				length = CalculateRayLength(envelope, direction, point_outside);
				if (length > 0 && !point_outside) {
					if (start_azimuth == 0)
						start_azimuth = azimuths(row, dir_index);
					if (end_azimuth != azimuths(row, start_col + col_step - 1))
						end_azimuth = azimuths(row, dir_index);

				}
			}
		bottom_width = abs(sin(end_azimuth - start_azimuth)) * l * sin(start_zenith);

		start_azimuth = 0;
		end_azimuth = 0;
		for (int delta = 0; delta < col_step / 2; delta++) {
			dir_index = mid_col + delta;
			Vector3f direction = dir.row(dir_index * ray.n_rows() + next_row);
			bool point_outside = false;
			length = CalculateRayLength(envelope, direction, point_outside);
			if (length > 0 && !point_outside) {
				end_zenith = zeniths(next_row, dir_index);
				l = length;
				break;
			}
			dir_index = mid_col - delta;
			direction = dir.row(dir_index * ray.n_rows() + next_row);
			length = CalculateRayLength(envelope, direction, point_outside);
			if (length > 0 && !point_outside) {
				end_zenith = zeniths(next_row, dir_index);
				l = length;
				break;
			}
		}
		double path_length = 0, counter = 0;
		if (use_gap_p) {
			for (int row_index = row; row_index < next_row; row_index++)
				for (int col_index = start_col; col_index < start_col + col_step; col_index++)
					if (gap_p_valid(row_index, col_index)) {
						path_length += gap_p(row_index, col_index);
						counter++;
					}
		} else {
			for (int row_index = row; row_index < next_row; row_index++)
				for (int col_index = start_col; col_index < start_col + col_step; col_index++) {
					if (path_lengths_image(row_index, col_index) > 0) {
						path_length += path_lengths_image(row_index, col_index);
					}
					counter++;
				}
		}
		path_length /= counter;


		if (i > 0 && has_canopy_next[i - 1]) {
			start_azimuth = azimuths(next_row, start_col);
		}
		if (i < number_of_divisions - 1 && has_canopy_next[i + 1]) {
			end_azimuth = azimuths(next_row, start_col + col_step - 1);
		}
		if (start_azimuth == 0 || end_azimuth == 0) {
			for (dir_index = start_col; dir_index < start_col + col_step; dir_index++) {
				Vector3f direction = dir.row(dir_index * ray.n_rows() + next_row);
				bool point_outside = false;
				length = CalculateRayLength(envelope, direction, point_outside);
				if (length > 0 && !point_outside) {
					if (start_azimuth == 0)
						start_azimuth = azimuths(next_row, dir_index);
					if (end_azimuth != azimuths(next_row, start_col + col_step - 1))
						end_azimuth = azimuths(next_row, dir_index);
				}
			}
		}
		top_width = abs(sin(end_azimuth - start_azimuth)) * l * sin(end_zenith);
		double sub_height = (start_zenith - end_zenith) * l / sin(end_zenith);
		double sub_area = sub_height * (bottom_width + top_width) / 2;
		//cout << top_width << ' ' << bottom_width << ' ' << sub_area << ' ' <<sub_height << ' ' << l << ' ' << path_length << endl;
		// cout << start_azimuth << ' ' << end_azimuth << ' ' << start_azimuth / 3.14 * 180 << ' ' << end_azimuth / 3.14 * 180 << ' ' << l << endl;
		if (end_zenith == 0 || l == 0) {	// 光打在树外，找不到了
			sub_area = 0;
			path_length = 0;
		}
		else if (sub_height < 0 || sub_area < 0) {
			sub_area = 0;
			path_length = 0;
		}
		else if (path_length < 0) {
			path_length = 0;
		}
		else {
			height = sub_height;
		}
		area += sub_area;
		// cout << sub_area * path_length << ' ' << sub_area << ' ' << path_length << endl;
		if (sub_area * path_length > 0)
			volumn += sub_area * path_length;
	}
	// cout << endl;
	return area;
}