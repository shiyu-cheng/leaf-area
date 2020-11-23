#include "OBJreader.h"

bool OBJreader::read()
{
	printf("Envelope data (.obj): Reading data...");
	clock_t launch = clock();

	if (filename == NULL)
	{
		printf("ERROR: File name not set\n");
		return false;
	}

	FILE *fp;
	if (0 != _wfopen_s(&fp, filename, L"rb"))
	{
		printf("ERROR: Open file failed\n");
		return false;
	}

	n_vertices = 0;
	n_triangles = 0;

	char tmpline[100];
	while (fgets(tmpline, 100, fp))
	{
		
		if (tmpline[0] == 'v')
		{	
			++n_vertices;
		}
		else if (tmpline[0] == 'f')
		{
			++n_triangles;
		}
		else if (tmpline[0] == 'g' )
		{
			break;
		}
	} 
	rewind(fp);

	this->vertices.resize(n_vertices, 3);
	this->triangles.resize(n_triangles, 3);
	
	n_vertices = 0;
	n_triangles = 0;
	while (fgets(tmpline, 100, fp))
	{
		if (tmpline[0] == 'v')
		{
			sscanf_s(tmpline, "v %lf %lf %lf",  &this->vertices(n_vertices,0), &this->vertices(n_vertices,1), &this->vertices(n_vertices,2));
			++n_vertices;
		}
		else if (tmpline[0] == 'f')
		{
			sscanf_s(tmpline, "f %d %d %d", &this->triangles(n_triangles, 0), &this->triangles(n_triangles, 1), &this->triangles(n_triangles, 2));
			++n_triangles;
		} 
		else if (tmpline[0] == 'g')
		{
			clock_t done = clock();
			printf(" %.3lf s\n", double(done - launch) / CLOCKS_PER_SEC);
			return 1;
			
		} 
		else
		{
			printf("Warning: Unexpected data\n");
		}
	
	} 

	

	return true;
}

Vector3d OBJreader::get_triangle_vective(unsigned index_triangle, unsigned index_vectice)
{
	Vector3d v = this->vertices.row(this->triangles(index_triangle, index_vectice)-1);
	//Vector3d v2 = this->vertices.row(this->triangles(index_triangle, index_vectice) - 1).transpose();
	return v;
}

int OBJreader::num_triangles()
{
	return this->n_triangles;
}

void OBJreader::set_scanner_position(Vector3d _scanner_position)
{
	this->scanner_position_ = _scanner_position;
}



void OBJreader::shift_to_observational_frame()
{
	if (!this->is_shifted)
	{
		this->vertices.rowwise() -= this->scanner_position_.adjoint();
		this->scanner_position_ = Vector3d(0.0, 0.0, 0.0);
		this->is_shifted = true;
	}
}

bool OBJreader::is_scanner_too_near()
{
	double _x_max = vertices.col(0).maxCoeff();
	double _y_max = vertices.col(1).maxCoeff();
	double _z_max = vertices.col(2).maxCoeff();

	double _x_min = vertices.col(0).minCoeff();
	double _y_min = vertices.col(1).minCoeff();
	double _z_min = vertices.col(2).minCoeff();


	Vector3d _size_envelope = Vector3d((_x_max - _x_min) / 2, (_y_max - _y_min) / 2, (_z_max - _z_min) / 2);
	Vector3d _center_envelope = Vector3d((_x_max + _x_min) / 2, (_y_max + _y_min) / 2, (_z_max + _z_min) / 2);
	
	double _r_envelope = sqrt(_size_envelope.squaredNorm()/3);
	
	double _distance = (_center_envelope - this->scanner_position_).norm();
	cout << "d  "<<_distance<<endl;
	cout << "r  "<<_r_envelope <<endl;
	

	if (_distance < _r_envelope + 0.6)
	{
		return true;
	}

	return false;
}

// Azimuth min, max, zenith min, max
angle_extent OBJreader::get_angle_extent()
{
	this->shift_to_observational_frame();
	ArrayXf zeniths = this->vertices.col(2).cwiseQuotient(this->vertices.rowwise().norm()).array().cast<float>().acos();   //cos(zenith) = z
	ArrayXf azimuths = this->vertices.col(1).cwiseQuotient(this->vertices.col(0)).array().cast<float>().atan();   //tan(azimuth) = y/x
	azimuths = (this->vertices.col(0).array() < 0).select(azimuths + M_PI, azimuths);
	azimuths = (this->vertices.col(0).array() > 0 && this->vertices.col(1).array() < 0).select(azimuths + 2*M_PI, azimuths);


	float azimuth_min = azimuths.minCoeff(); 
	float azimuth_max = azimuths.maxCoeff();

	//ArrayXf min_azimuth_col = azimuth_image.isNaN().select(99, azimuth_image).colwise().minCoeff();

	if ((azimuth_max - azimuth_min) > M_PI)
	{
		azimuth_min = (azimuths > (float)M_PI).select(azimuths, 10).minCoeff();
		azimuth_max = (azimuths < (float)M_PI).select(azimuths, 0).maxCoeff();
	} 


	//ofstream file("env_azimuth.txt");
	//if (file.is_open())
	//{
	//	file << azimuths;
	//	file.close();
	//}


	this->angle_extent_envelope_ = { azimuth_min, azimuth_max, zeniths.minCoeff(), zeniths.maxCoeff() };

	//printf("TLS data (.ptx): Extent mask set: Azimuth = [%.2f, %.2f], Zenith = [%.2f, %.2f] (rad)\n", \
	//	this->angle_extent_mask_.azimuth_min, \
	//	this->angle_extent_mask_.azimuth_max, \
	//	this->angle_extent_mask_.zenith_min, \
	//	this->angle_extent_mask_.zenith_max);
	printf("                      Extent: Azimuth = [%.2f, %.2f], Zenith = [%.2f, %.2f] (rad)\n", \
		this->angle_extent_envelope_.azimuth_min, \
		this->angle_extent_envelope_.azimuth_max, \
		this->angle_extent_envelope_.zenith_min, \
		this->angle_extent_envelope_.zenith_max);

	return this->angle_extent_envelope_;
}


