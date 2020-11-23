#include "LAIpathTLS.h"
#include <fstream>
#include <ctime>
#include <list>
#include <windows.h>
#include <cstdio>
#include <vector>
#include "Freeimage.h"
#include <direct.h>  
#include "GFun.h"

using namespace std;



ArrayXXb RemoveEdgeGap(ArrayXXb within_canopy_image, ArrayXXb all_gaps_image);
double GapWithinCrownProbabilityByZenith(ArrayXXb _within_canopy, ArrayXXb _is_gaps, ArrayXXb _zenith_layer);
ArrayXd  DistributionByMask(ArrayXXd _path_lengths, ArrayXXb _mask);
int CalcFromData(const wchar_t* ptx_path, const wchar_t* envelope_path, wchar_t* append, bool refine_envelope = false, float intensity_threshold = 0.2f, float zenith_low = 0.0f, float zenith_high = 180.0f,
	float zenith_increment = 10.0f, float zenith_layers = 0.0f, float window_width = 0.04, float window_height = 0.04, int n_divisions = 1, bool use_finite_avg = false, bool use_gap_p = false, bool cache_pathlen = false);

void get_debug_img_path(char* img_path, const wchar_t* ptx_path, const wchar_t* envelope_path, const char* debug_dir, const char* append = "");

void write_debug_img(const char* img_path, Eigen::ArrayXXf& _image);
void write_debug_img(const char* img_path, Map<ArrayXXf>& _image);
void write_debug_img(const char* img_path, ArrayXXb& _image);

list<wstring> GetFileNames(wchar_t* lpPath);
list<wstring> GetFileNames(wchar_t* lpPath, wchar_t* ext);
wchar_t* output_path(wchar_t* new_path, const wchar_t* original_path, const wchar_t* fname_append, const wchar_t* new_ext);

bool img_debug = false;
//bool set_resolution = false;
float g_resolution_TLS = 8.f / 100.0f;
//float g_resolution_TLS = 360.f / 3600.0f;
unsigned g_n_skip = 1;
bool b_zenith_layer = false;


double* g_Gfun_values = 0;
int g_n_bins_gMes = -1;		//实测g函数行数
double g_angle_interval_G = 0.1;

int wmain(int argc, wchar_t* argv[])
{
	float zenith_low = 0.0; //degrees for input, rad for calculation 
	float zenith_high = 180.0;
	float zenith_increment = 10.0;
	float zenith_layers = 0.0;
	float window_width = 0.04;
	float window_height = 0.04;
	int n_divisions = 1;

	float intensity_threshold = 0.2f;

	//strcpy_s(ptx_path, "d:\\Data\\Jardin\\Meta\\LiDAR\\Ronghai_Elena\\20150525\\input_arbre\\arbre004_20m_YXZ000.ptx");
	//strcpy_s(envelope_path, "d:\\Data\\Jardin\\Meta\\LiDAR\\Ronghai_Elena\\20150525\\TA4_envelope_convexe.obj");

#pragma region Command line

	bool verbose = false;
	bool output_to_file = false;
	bool refine_envelope = false;
	bool cache_pathlen = false;
	bool use_gap_p = false, use_finite_avg = false;
	wchar_t ptx_path[_MAX_PATH] = L"";

	wchar_t envelope_path[_MAX_PATH] = L"";
	wchar_t gMes_path[_MAX_PATH] = L"";

	double* gMes = 0;			//实测g函数


	wchar_t append[_MAX_PATH] = L"_FAVD";



	if (argc == 1)
	{
		fprintf(stderr, "%ls is better run in the command line\n", argv[0]);

		fprintf(stderr, "enter input ptx file: "); fgetws(ptx_path, _MAX_PATH, stdin);

		fprintf(stderr, "enter input envelope_path file: "); fgetws(envelope_path, _MAX_PATH, stdin);

		//fprintf(stderr, "enter output file: "); fgetws(out_path, _MAX_PATH, stdin);
	}

	for (int i = 1; i < argc; i++)
	{
		if (argv[i][0] == '\0')
		{
			continue;
		}
		else if (wcscmp(argv[i], L"-v") == 0 || wcscmp(argv[i], L"-verbose") == 0)
		{
			verbose = true;
		}
		else if (wcscmp(argv[i], L"-refine") == 0)
		{
			refine_envelope = true;
			wcscat(append, L"_REFINED");
		}
		else if (wcscmp(argv[i], L"-img_debug") == 0)
		{
			img_debug = true;
		}
		else if (wcscmp(argv[i], L"-cache_pathlen") == 0)
		{
			cache_pathlen = true;
		}
		else if (wcscmp(argv[i], L"-use_gap_p") == 0)
		{
			use_gap_p = true;
		}
		else if (wcscmp(argv[i], L"-use_finite_avg") == 0)
		{
			use_finite_avg = true;
		}
		else if (wcscmp(argv[i], L"-ptx") == 0)
		{

			if ((i + 1) >= argc)
			{
				fprintf(stderr, "ERROR: '%ls' needs 1 argument: stop\n", argv[i]);
			}
			wcscpy_s(ptx_path, argv[i + 1]);

			++i;
		}
		else if (wcscmp(argv[i], L"-g") == 0)  //new: 2018-01-22
		{

			if ((i + 2) >= argc)
			{
				fprintf(stderr, "ERROR: '%ls' needs 1 argument: stop\n", argv[i]);
			}
			wcscpy_s(gMes_path, argv[i + 1]);
			g_n_bins_gMes = _wtoi(argv[i + 2]);
			i += 2;

			gMes = new double[g_n_bins_gMes];
			errno_t err;
			if (err = ReadMes(gMes_path, gMes, g_n_bins_gMes))
			{
				fprintf(stderr, "ERROR: Read %ls failed, error code '%d'\n", argv[i + 1], err);
			}

		}
		else if (wcscmp(argv[i], L"-envelope") == 0)
		{

			if ((i + 1) >= argc)
			{
				fprintf(stderr, "ERROR: '%ls' needs 1 argument: stop\n", argv[i]);
			}
			wcscpy_s(envelope_path, argv[i + 1]);
			++i;
		}
		else if (wcscmp(argv[i], L"-intensity_threshold") == 0)
		{

			if ((i + 1) >= argc)
			{
				fprintf(stderr, "ERROR: '%ls' needs 1 argument: stop\n", argv[i]);
			}
			intensity_threshold = (float)_wtof(argv[i + 1]);
			++i;
		}
		else if (wcscmp(argv[i], L"-resolution_degree") == 0)
		{

			if ((i + 1) >= argc)
			{
				fprintf(stderr, "ERROR: '%ls' needs 1 argument: stop\n", argv[i]);
			}
			g_resolution_TLS = (float)_wtof(argv[i + 1]);
			++i;
		}
		else if (wcscmp(argv[i], L"-width") == 0)
		{

			if ((i + 1) >= argc)
			{
				fprintf(stderr, "ERROR: '%ls' needs 1 argument: stop\n", argv[i]);
			}
			window_width = (float)_wtof(argv[i + 1]);
			++i;
		}
		else if (wcscmp(argv[i], L"-height") == 0)
		{

			if ((i + 1) >= argc)
			{
				fprintf(stderr, "ERROR: '%ls' needs 1 argument: stop\n", argv[i]);
			}
			window_height = (float)_wtof(argv[i + 1]);
			++i;
		}
		else if (wcscmp(argv[i], L"-n_divisions") == 0)
		{

			if ((i + 1) >= argc)
			{
				fprintf(stderr, "ERROR: '%ls' needs 1 argument: stop\n", argv[i]);
			}
			n_divisions = (int)_wtof(argv[i + 1]);
			++i;
		}
		else if (wcscmp(argv[i], L"-resolution_rad") == 0)
		{

			if ((i + 1) >= argc)
			{
				fprintf(stderr, "ERROR: '%ls' needs 1 argument: stop\n", argv[i]);
			}
			g_resolution_TLS = (float)_wtof(argv[i + 1]) / M_PI * 180.0f;
			++i;
		}
		else if (wcscmp(argv[i], L"-skip") == 0)
		{

			if ((i + 1) >= argc)
			{
				fprintf(stderr, "ERROR: '%ls' needs 1 argument: stop\n", argv[i]);
			}
			g_n_skip = (unsigned)_wtoi(argv[i + 1]);
			++i;
		}
		else if (wcscmp(argv[i], L"-o") == 0)
		{

			if ((i + 1) >= argc)
			{
				fprintf(stderr, "ERROR: '%ls' needs 1 argument: stop\n", argv[i]);
			}
			wcscpy_s(append, argv[i + 1]);
			wcscat(append, L"\\out.txt");
			++i;
		}
		else if (wcscmp(argv[i], L"-otxt") == 0)
		{
			output_to_file = true;
		}
		else if (wcscmp(argv[i], L"-zenith_ranges") == 0)
		{
			if ((i + 3) >= argc)
			{
				fprintf(stderr, "ERROR: '%ls' needs 3 argument: stop\n", argv[i]);
			}
			b_zenith_layer = true;
			zenith_low = (float)_wtof(argv[i + 1]);
			zenith_high = (float)_wtof(argv[i + 2]);
			zenith_increment = (float)_wtof(argv[i + 3]);

			zenith_low = (float)(zenith_low / 180.0f * M_PI);
			zenith_high = (float)(zenith_high / 180.0f * M_PI);
			zenith_increment = (float)(zenith_increment / 180.0f * M_PI);

			i = i + 3;
		}
		else if (wcscmp(argv[i], L"-zenith_layers") == 0)
		{
			if ((i + 1) >= argc)
			{
				fprintf(stderr, "ERROR: '%ls' needs 1 argument: stop\n", argv[i]);
			}
			b_zenith_layer = true;
			zenith_layers = (float)_wtof(argv[i + 1]);

			i = i + 1;
		}
		else
		{
			fprintf(stderr, "ERROR: cannot understand argument '%ls'\n", argv[i]);
		}
	}

	if (!output_to_file)
	{
		wcscpy_s(append, L"-stdout");
	}



	//GMes

	int panels_G = 1800;
	double max_angle_G = M_PI;
	double min_angle_G = 0.0;
	g_angle_interval_G = (max_angle_G - min_angle_G) / panels_G;

	if (g_n_bins_gMes > 0)
	{
		g_Gfun_values = new double[panels_G + 1];

		#pragma omp parallel for 
		for (int i = 0; i <= panels_G; ++i)
		{
			double u1 = cos(max_angle_G + i * g_angle_interval_G);
			g_Gfun_values[i] = GFun(abs(u1), gMes, g_n_bins_gMes);
		}

		//G_value = GMean(g_Gfun_values, envelope_extent.zenith_min, envelope_extent.zenith_max, angleInterval);
	}

#pragma endregion Command line


	CalcFromData(ptx_path, envelope_path, append, refine_envelope, intensity_threshold, zenith_low, zenith_high, zenith_increment,
		         zenith_layers, window_width, window_height, n_divisions, use_finite_avg, use_gap_p, cache_pathlen);
	
	delete[] g_Gfun_values;
	g_Gfun_values = 0;
	delete[] gMes;
	gMes = 0;
	return 0;
	
#pragma region batch processing
	list <wstring> fpaths_ptx, fpaths_envelope;

	//fpaths_ptx = GetFileNames(ptx_path);
	fpaths_ptx = GetFileNames(ptx_path, L"ptx");

	////fpaths_envelope = GetFileNames(envelope_path);
	fpaths_envelope = GetFileNames(envelope_path, L"obj");

	for (list<wstring>::iterator it_path_env = fpaths_envelope.begin(); it_path_env != fpaths_envelope.end(); it_path_env++)
	{
		for (list<wstring>::iterator it_path_ptx = fpaths_ptx.begin(); it_path_ptx != fpaths_ptx.end(); it_path_ptx++)
		{
			printf("Input (.ptx): %ls\n", (*it_path_ptx).c_str());
			printf("      (envelope): %ls\n", (*it_path_env).c_str());

			CalcFromData((*it_path_ptx).c_str(), (*it_path_env).c_str(), append, refine_envelope, intensity_threshold, zenith_low, zenith_high, zenith_increment);
		}
	}

	delete[] g_Gfun_values;
	g_Gfun_values = 0;
	delete[] gMes;
	gMes = 0;

	//
	//VectorXd zeniths = envelope
	//cout << path_lengths_image << within_canopy_image << all_gaps_image << points_outside_canopy_image << zenith_mask_image << points_not_reach_canopy_image << endl;
	//cout << intensity_image(0) <<  endl;
	//getc(stdin);
	return 0;
#pragma endregion batch processing
}


int CalcFromData(const wchar_t* ptx_path, const wchar_t* envelope_path, wchar_t* append, bool refine_envelope, float intensity_threshold, float zenith_low, float zenith_high,
	             float zenith_increment, float zenith_layers, float window_width, float window_height, int n_divisions, bool use_finite_avg, bool use_gap_p, bool cache_pathlen)
{
	ArrayXXb within_canopy_image;	// pulses that contact canopy 
	ArrayXXb all_gaps_image;	// returns without xyz (x,y,z) = (0,0,0)
	ArrayXXb points_not_reach_canopy_image; // // returns with xyz but do not reach canopy ( too near )

	FILE* file_out;
	wchar_t out_path[_MAX_PATH] = L"";
	if (wcscmp(append, L"-stdout") == 0)
	{
		file_out = stdout;
	}
	else
	{
		output_path(out_path, envelope_path, append, L"txt");

		file_out = _wfopen(out_path, L"a");
		if (file_out == 0)
		{
			fprintf(stderr, "ERROR: could not open '%ls' for write\n", out_path);
			file_out = stdout;
		}
	}


	time_t timeNow = time(NULL);
	struct tm p;
	localtime_s(&p, &timeNow);
	fprintf_s(file_out, "Process at %4d.%02d.%02d %02d:%02d:%02d Copyright by Ronghai HU (sea@mail.bnu.edu.cn)\n", p.tm_year + 1900, p.tm_mon + 1, p.tm_mday, p.tm_hour, p.tm_min, p.tm_sec);

	fprintf(file_out, "Input (.ptx): %ls\n", ptx_path);
	fprintf(file_out, "      (envelope): %ls\n", envelope_path);
	fprintf(file_out, "Output : %ls\n", out_path);
	fprintf(file_out, "Intensity threshold : %.2f\n", intensity_threshold);
	fprintf(file_out, "Skip : %d\n", g_n_skip);
	if (refine_envelope) fprintf(file_out, "-refine\n");

	//fclose(file_out);

	//cout << "zenith_low: " << zenith_low << endl;
	//cout << "zenith_high: " << zenith_high << endl;
	//cout << "zenith_increment: " << zenith_increment << endl;

	size_t n_bins = (size_t)round((zenith_high - zenith_low) / zenith_increment);		//zenith bin number
	VectorXf zenith_ranges = VectorXf::LinSpaced(n_bins + 1, zenith_low, zenith_high);

#pragma region Read point cloud and envelope data 

	// Read point cloud data
	//PTXreader ray("c:\\Work\\ICube\\TLS\\arbre004_20m_YXZ000.ptx");

#ifdef _DEBUG
	//PTXreader ray(ptx_path, 10, 10);
	PTXreader ray(ptx_path);
#else
	PTXreader ray(ptx_path, g_n_skip, g_n_skip);
#endif // _DEBUG
	ray.set_resoultion_degree(g_resolution_TLS);
	if (!ray.open()) return -1;
	if (!ray.read_header()) return -1;
	Vector3d scanner_position = ray.get_scanner_position();


	// Read envelope data
	//OBJreader envelope("c:\\Work\\ICube\\TLS\\TA4_envelope_convexe.obj");
	OBJreader envelope(envelope_path);
	if (!envelope.read()) return -1;
	envelope.set_scanner_position(scanner_position);

	if (envelope.is_scanner_too_near())
	{
		fprintf(file_out, "Laser scanner is too near\n\n");
		ray.close();
		fclose(file_out);
		return 1;
	}


	envelope.shift_to_observational_frame();

	angle_extent envelope_extent = envelope.get_angle_extent();


	ray.set_angle_extent_mask(envelope_extent);
	ray.set_intensity_threshold(intensity_threshold);

	if (ray.read_data() == -1)
	{
		fprintf(file_out, "No intersection\n\n");
		ray.close();
		fclose(file_out);
		return 2;
	}
	ray.close();

	angle_extent data_extent = ray.get_angle_extent();

	ArrayXXf zeniths = ray.get_zeniths_image();
	ArrayXXf azimuths = ray.get_azimuths_image();

	char debug_dir[_MAX_DIR] = "debug";

	//if (img_debug)
	//{
	//	//ArrayXXf norms_image = ray.get_norms_image();
	//	//char debug_dir[_MAX_DIR] = "debug";
	//	_mkdir(debug_dir);

	//	char img_path_zenith[_MAX_PATH];
	//	get_debug_img_path(img_path_zenith, ptx_path, envelope_path, debug_dir, "_zenith_0");
	//	write_debug_img(img_path_zenith, zeniths);

	//	char img_path_azimuth[_MAX_PATH];
	//	get_debug_img_path(img_path_azimuth, ptx_path, envelope_path, debug_dir, "_azimuth_0");
	//	write_debug_img(img_path_azimuth, azimuths);
	//	return 0;
	//}


#ifdef _DEBUG

	//ArrayXXf norms_image = ray.get_norms_image();
	ArrayXXf intensity_image = ray.get_intensity_image();
#else
	ArrayXXf intensity_image = ray.get_intensity_image();
#endif


	//New: 2017-11-16 add image debug output
	if (img_debug)
	{
		ArrayXXf norms_image = ray.get_norms_image();
		//char debug_dir[_MAX_DIR] = "debug";
		_mkdir(debug_dir);

		char img_path_zenith[_MAX_PATH];
		get_debug_img_path(img_path_zenith, ptx_path, envelope_path, debug_dir, "_zenith");
		write_debug_img(img_path_zenith, zeniths);

		char img_path_azimuth[_MAX_PATH];
		get_debug_img_path(img_path_azimuth, ptx_path, envelope_path, debug_dir, "_azimuth");
		write_debug_img(img_path_azimuth, azimuths);

		char img_path_intensity[_MAX_PATH];
		get_debug_img_path(img_path_intensity, ptx_path, envelope_path, debug_dir, "_intensity");
		write_debug_img(img_path_intensity, intensity_image);

		char img_path_norms[_MAX_PATH];
		get_debug_img_path(img_path_norms, ptx_path, envelope_path, debug_dir, "_norms");
		write_debug_img(img_path_norms, norms_image);
	}


	// Zeniths interpolation
	Map<ArrayXXf> zeniths_interpolated_image = ray.get_interpolated_zeniths_image();
	Map<ArrayXXf> azimuth_interpolated_image = ray.get_azimuths_image();

	ArrayXXf norms_image2 = ray.get_norms_image();


	//wchar_t fname[_MAX_FNAME];

		/*FreeImage_ConvertFromRawBits(BYTE *bits, int width, int
			height, int pitch, unsigned bpp, unsigned red_mask, unsigned green_mask, unsigned
			blue_mask, BOOL topdown FI_DEFAULT(FALSE));*/

			//unsigned bpp = 32; 
			//int pitch = ((((bpp * intensity_image.rows()) + 31) / 32) * 4);

			//FIBITMAP* dib = FreeImage_ConvertFromRawBits((BYTE*)intensity_image.data(), intensity_image.rows(), intensity_image.cols(), pitch, bpp , 0, 0, 0, FALSE);
			//FreeImage_Save(FIF_TIFF, dib, "blob.tif", TIFF_DEFAULT);
			//FreeImage_Unload(dib);

			// attach the binary data to a memory stream
			//FIMEMORY *hmem = FreeImage_OpenMemory((BYTE*)intensity_image.data(), intensity_image.size()*sizeof(float)/sizeof(BYTE));

			// get the file type
			//FREE_IMAGE_FORMAT fif = FreeImage_GetFileTypeFromMemory(hmem, 0);
			// load an image from the memory stream
			//FIBITMAP *check = FreeImage_LoadFromMemory(fif, hmem, 0);
			// save as a regular file
			//FreeImage_Save(FIF_TIFF, check, "blob.tif", TIFF_DEFAULT);
			//FreeImage_Unload(check);
			//FreeImage_CloseMemory(hmem);

		//}


		//fprintf(file_out, "Data\t%16d\t%11.1f\t%11.1f\t%10.1f\t%10.1f\n", 1000000000, 10.0f, 350.0f, 20.0f, 70.0f);
		//fprintf(file_out, "Data\t%16d\t(%4.2f) %4.1f\t(%4.2f) %4.1f\t(%4.2f) %4.1f\t(%4.2f) %4.1f\n", 1000000000, 10.0f / 180 * M_PI, 10.0f, 350.0f / 180 * M_PI, 350.0f, 20.0f / 180 * M_PI, 20.0f, 70.0f / 180 * M_PI, 70.0f);

		//fprintf(file_out, "%.1lf-%.1lf\t%6.4lf\t%15.4lf\t%16lld%15.2lf\t%15.2lf\n",
	fprintf(file_out, "\n\tNumber_of_pulses\tAzimuth_min\tAzimuth_max\tzenith_min\tzenith_max\n");
	fprintf(file_out, "Data\t%16lld\t%11.1f\t%11.1f\t%10.1f\t%10.1f\n", \
		ray.n_points_data(), \
		data_extent.azimuth_min * 180 / M_PI, \
		data_extent.azimuth_max * 180 / M_PI, \
		data_extent.zenith_min * 180 / M_PI, \
		data_extent.zenith_max * 180 / M_PI);


#pragma endregion Read point cloud and envelope data (.ptx)

	// 1. zenith mask
	ArrayXXb zenith_mask_image = (zeniths_interpolated_image > zenith_low) && (zeniths_interpolated_image < zenith_high);

	//size_t tmp_s = sizeof(zenith_mask_image(0, 0));
	//size_t tmp_s2 = sizeof(bool);
	//size_t tmp_s3 = sizeof(UINT16);

	string PathLengthsFile = "cached_pathlen.txt";
	ArrayXf path_lengths = VectorXf::Constant(ray.n_points(), -99);
	ArrayXXb points_outside_canopy_image = ray.get_gaps_image();
	std::ifstream file(PathLengthsFile);
	if (file.is_open() && cache_pathlen) {
		for (int i = 0; i < ray.n_points(); ++i)
			if (zenith_mask_image(i))
			{
				file >> path_lengths(i) >> points_outside_canopy_image(i);
			}
		file.close();
	}
	else {
		path_lengths = CalcPathLengths(envelope, ray, points_outside_canopy_image, zenith_mask_image);

		std::ofstream file(PathLengthsFile);
		if (file.is_open())
		{
			for (int i = 0; i < ray.n_points(); ++i)
				if (zenith_mask_image(i))
				{
					file << path_lengths(i) << ' ';
					file << points_outside_canopy_image(i) << ' ';
				}
			file.close();
		}
	}
	Map<ArrayXXf> path_lengths_image(path_lengths.data(), ray.n_rows(), ray.n_cols());


	// 2. within_canopy mask
	within_canopy_image = (path_lengths_image > 0);

	// remove gaps outside envelope
	ArrayXXb within_canopy_image_new;
	if (refine_envelope)
	{
		within_canopy_image_new = RemoveEdgeGap(within_canopy_image, all_gaps_image);
	}
	else
	{
		within_canopy_image_new = within_canopy_image;
	}

	points_not_reach_canopy_image = ((path_lengths_image).cast<int>() == -2);



	fprintf(file_out, "Canopy\t%16lld\t%11.1f\t%11.1f\t%10.1f\t%10.1f\n", \
		within_canopy_image.count(), \
		envelope_extent.azimuth_min * 180 / M_PI, \
		envelope_extent.azimuth_max * 180 / M_PI, \
		envelope_extent.zenith_min * 180 / M_PI, \
		envelope_extent.zenith_max * 180 / M_PI);


	//float zenith_max_within_canopy = within_canopy_image.select(zeniths_interpolated_image, -99).maxCoeff();
	//float zenith_min_within_canopy = within_canopy_image.select(zeniths_interpolated_image, 99).minCoeff();
	//fprintf(file_out, "zenith_min_within_canopy: %f (%f degrees)\n", zenith_min_within_canopy, zenith_min_within_canopy * 180 / M_PI);
	//fprintf(file_out, "zenith_max_within_canopy: %f (%f degrees)\n", zenith_max_within_canopy, zenith_max_within_canopy * 180 / M_PI);
	if (zenith_layers > 0) {
		n_bins = zenith_layers;
		zenith_ranges = VectorXf::LinSpaced(n_bins + 1, envelope_extent.zenith_min, envelope_extent.zenith_max);
	}


	// 3. all gap mask
	all_gaps_image = ray.get_gaps_image() || points_outside_canopy_image;
	//cout << "n_gaps: " << is_gaps.count() << endl;
	//cout << "n_canopy_gap1 (is_gaps): " << ( within_canopy * is_gaps).count() << endl;
	//cout << "n_canopy_gap2 (points_outside_canopy): " << (points_outside_canopy * within_canopy).count() << endl;

	double G_value = 0.5;
	/*START-----------------------12/18/2019-----------------------------*/


	/* print envelope shape
	for (Index i = 0; i < all_gaps_image.rows(); ++i) {
		for (Index j = 0; j < all_gaps_image.cols(); ++j) {
			if (rand() % 10000 < 5)
				cout << "i:" << i << " j:" << j << " zeniths:" << zeniths_interpolated_image(i, j) << "cos(zenith): " << cos(zeniths_interpolated_image(i, j)) << endl;
		}
	}*/


	ArrayXXf gap_p(all_gaps_image.rows(), all_gaps_image.cols());
	ArrayXXb gap_p_valid_region(all_gaps_image.rows(), all_gaps_image.cols());
	cout << "Image size: " << gap_p.rows() << ", " << gap_p.cols() << endl;

	double WIDTH = window_width;
	double HEIGHT = window_height;
	//double total_path_len = 0, total_LAI = 0;
	double total_costheta = 0, total_sintheta = 0, total_LAI = 0, total_LA = 0, total_area = 0, total_height = 0, total_volumn = 0;
	Index v = 0, u = 0;
	int mid_col = 1 + ray.n_cols() / 2;
	int mid_row = 1 + ray.n_rows() / 2;
	FILE* file_debug = _wfopen(L"distribution.txt", L"a");
	fprintf(file_debug, "Start\n");
	//fprintf(file_debug, "window width:%lf window height::%lf\n", window_width, window_height);
	if (use_gap_p) {
		fprintf(file_debug, "Path length: -ln(slide window gap probability)\n");
		fprintf(file_debug, "Window width: %lfrad, height: %lfrad\n", window_width, window_height);
	}
	else if (use_finite_avg)
	{
		fprintf(file_debug, "Finite Average\n");
	}
	else
		fprintf(file_debug, "Path length: beam length through crown\n");
	fprintf(file_debug, "data = %ls, number of divisions = %d\n", ptx_path, n_divisions);
	fclose(file_debug);
	
	// for each row slice
	for (size_t slice_id = 0; slice_id < n_bins; slice_id++)
	{
		Index row_upper = 0, row_lower = 0;
		// The two bounds of the slice. row_upper has smaller zenith than row_lower
		cout << "zenith range for slice: " << zenith_ranges(slice_id) * 180 / 3.14159 << " to " << zenith_ranges(slice_id + 1) * 180 / 3.14159 << endl;
		while (!(zeniths_interpolated_image(row_upper, mid_col) > zenith_ranges(slice_id) - 0.01
			&& zeniths_interpolated_image(row_upper, mid_col) < zenith_ranges(slice_id) + 0.01) && row_upper < gap_p.rows()) row_upper++;
		while (!(zeniths_interpolated_image(row_lower, mid_col) > zenith_ranges(slice_id+1) - 0.01
			&& zeniths_interpolated_image(row_lower, mid_col) < zenith_ranges(slice_id+1) + 0.01) && row_lower < gap_p.rows()) row_lower++;
		cout << "row_lower: " << row_lower << " row_upper: " << row_upper << endl;

		// zenith mask
		ArrayXXb zenith_layer_slice = (zeniths_interpolated_image > zenith_ranges(slice_id)) && (zeniths_interpolated_image < zenith_ranges(slice_id + 1));

		double all_samples = 0, valid_samples = 0;
		// reset valid region mask for the whole image
		for (Index i = 0; i < all_gaps_image.rows(); i += 1)
			for (Index j = 0; j < all_gaps_image.cols(); j += 1) {
				gap_p_valid_region(i, j) = 0;
			}

		float neg_log_p_sum = 0;
		cout << WIDTH << ' ' << HEIGHT << endl;
		vector<int> n_col_each_row;
		// sliding window in the row slice
		for (Index i = row_lower, j = 0; i < all_gaps_image.rows() && zeniths_interpolated_image(i, mid_col) - HEIGHT >= zenith_ranges(slice_id);) {
			j = 14;//7：28 14：38.8568
			n_col_each_row.push_back(0);
			//while (!within_canopy_image_new(i, j)) j+=5;
			for (; j < all_gaps_image.cols();) {
				all_samples++;
				double sum_path_len = 0;
				double p;
				if (use_finite_avg)
				{
					cout << "i & j: " << i << ' ' << j;
					// cout << "zenith & azimuth: " << zeniths_interpolated_image(i, mid_col) << ' ' << azimuth_interpolated_image(mid_row, j) << endl;
					p = P_of_the_window_finite_avg(i, j, all_gaps_image, within_canopy_image_new, path_lengths_image, WIDTH,
						HEIGHT, azimuth_interpolated_image, zeniths_interpolated_image, sum_path_len);
					cout << " P: " << p << endl;
				}
				else
				{
					p = P_of_the_window(i, j, all_gaps_image, within_canopy_image_new, path_lengths_image, WIDTH,
						HEIGHT, azimuth_interpolated_image, zeniths_interpolated_image, sum_path_len);
				}
				gap_p(i, j) = -log(p);
				if (p>0)
					neg_log_p_sum += -log(p);
				gap_p_valid_region(i, j) = p >= 0;
				valid_samples += p >= 0;
				//valid_samples ++;
				int v = 0;
				for (v = 0; j + v < all_gaps_image.cols() && azimuth_interpolated_image(mid_row, j + v) > azimuth_interpolated_image(mid_row, j) - WIDTH; v++);
				if (use_finite_avg) j += v;
				else j += 1;
				n_col_each_row[n_col_each_row.size() - 1] ++;
			}
			int u = 0;
			for (u = 0; i + u < all_gaps_image.rows() && zeniths_interpolated_image(i + u, mid_col) > zeniths_interpolated_image(i, mid_col) - HEIGHT; u++);
			if (use_finite_avg) i += u;
			else i += 1;
			// cout << i << ' ' << zeniths_interpolated_image(i, mid_col) - HEIGHT << ' ' << zeniths_interpolated_image(row, mid_col) - STEP_SIZE << endl;
		}
		
		cout << endl << "window numbers each row: " << endl;
		for (int i = 0; i < n_col_each_row.size(); i++)
			cout << n_col_each_row[i] << endl;
		//cout << "valid windows: " << valid_samples << endl;
		//cout << "all windows: " << all_samples << endl;

		//cout << "row_step: " << row_step << endl;
		if (valid_samples / all_samples < 0.01)
			continue;

		// The previous way to calculate gap fraction
		double gapFraction_ = GapWithinCrownProbabilityByZenith(within_canopy_image_new, all_gaps_image, zenith_layer_slice);
		// calculate the gap probability and the sum_path_len of the whole row slice
		double sum_path_len = 0;
		double gapFraction = P_of_the_window((row_lower+row_upper)/2, mid_col, all_gaps_image, within_canopy_image_new, path_lengths_image, 100,
			zenith_ranges(slice_id + 1) - zenith_ranges(slice_id), azimuth_interpolated_image, zeniths_interpolated_image, sum_path_len);
		//total_path_len += sum_path_len;
		//cout << "sum_path_len of the row slice: " << sum_path_len << endl;


		VectorXd neg_log_gap = DistributionByMask(gap_p.cast<double>(), gap_p_valid_region);
		VectorXd path_lengths_whinin_row = DistributionByMask(path_lengths_image.cast<double>(), zenith_layer_slice && within_canopy_image_new);
		//cout << "gap probability of the row slice: " << p << endl;

		if (neg_log_gap.size() <= 0) {
			cout << "no valid data" << endl << endl;
			continue;
		}
		gsl_histogram * gsl_hist = gsl_histogram_alloc(NUM_BINS);
		double max_val;
		if (use_gap_p)
			max_val = Stat_hist(neg_log_gap.data(), (unsigned long)neg_log_gap.size(), gsl_hist);
		else
			max_val = Stat_hist(path_lengths_whinin_row.data(), (unsigned long)path_lengths_whinin_row.size(), gsl_hist);

		gap_p = gap_p / max_val;
		Index delta = 0;
		double area = 0, height_sum = 0, volumn_sum = 0;
		for (delta = 0; row_lower + delta + 10 < row_upper; delta += 10)
		{
			double volumn = 0, height = 0;
			area += CalcRayArea2(within_canopy_image_new, envelope, ray, row_lower + delta, row_lower + delta + 10, n_divisions, use_gap_p, path_lengths_image, gap_p, gap_p_valid_region, height, volumn);
			height_sum += height;
			volumn_sum += volumn;
		}
		double volumn = 0, height = 0;
		area += CalcRayArea2(within_canopy_image_new, envelope, ray, row_lower + delta, row_upper, n_divisions, use_gap_p, path_lengths_image, gap_p, gap_p_valid_region, height, volumn);
		volumn_sum += volumn;
		height_sum += height;
		total_height += height_sum;
		total_volumn += volumn_sum;
		cout << "area: " << area << " height: " << height_sum << " volumn: " << volumn_sum << endl;
		total_area += area;
		// recently added
		if (use_finite_avg)
		{
			file_debug = _wfopen(L"distribution.txt", L"a");
			fprintf(file_debug, "zenith range for slice: %lf to %lf ", zenith_ranges(slice_id) * 180 / 3.1415926, zenith_ranges(slice_id + 1) * 180 / 3.1415926);  // * 180 / 3.14159

			cout << "valid samples: " << valid_samples << " all samples: " << all_samples << "G: " << G_value << endl;
			float finite_average_leaf_area = neg_log_p_sum / valid_samples / G_value * area;
			cout << "finite_average_LAI: " << neg_log_p_sum / valid_samples / G_value << "finite_average_LA: " << finite_average_leaf_area << endl;
			total_LA += finite_average_leaf_area;
			fprintf(file_debug, "LA: %lf\n", finite_average_leaf_area);
			fclose(file_debug);
		}
		else {
			file_debug = _wfopen(L"distribution.txt", L"a");
			fprintf(file_debug, "zenith range for slice: %lf to %lf ", zenith_ranges(slice_id) * 180 / 3.1415926, zenith_ranges(slice_id + 1) * 180 / 3.1415926);  // * 180 / 3.14159
			//fprintf(file_debug, "zenith range for slice: %lf to %lf\n", zenith_ranges(slice_id) * 180 / 3.1415926, zenith_ranges(slice_id + 1) * 180 / 3.1415926);  // * 180 / 3.14159
			//fprintf(file_debug, "Gap probability: new %lf, previous %lf\n", gapFraction, gapFraction_);
			//fprintf(file_debug, "Max value in gsl_hist: %lf\n", max_val);
			//gsl_histogram_fprintf(file_debug, gsl_hist, "%.2f", "%.3f");

			// calculate LAI for this row slice
			double zenith = zeniths_interpolated_image((row_lower + row_upper) / 2, all_gaps_image.cols() / 2);
			//cout << "zenith: " << zenith << endl;
			double costheta = 0;
			costheta = cos(zenith);
			double sintheta = 0;
			sintheta = sin(zenith);
			double LAImax = solve_LAImax(gsl_hist, gapFraction, G_value);

			double FAVD;
			if (use_gap_p)
				FAVD = LAImax;
			else
				FAVD = LAImax / max_val;

			//fprintf(file_debug, "LAImax: %lf\n", LAImax);
			//fprintf(file_debug, "FAVD: %lf\n", FAVD);
			fprintf(file_debug, "LA: %lf\n", FAVD* volumn_sum);
			double LAItrue = calc_LAItrue(gsl_hist, zenith, LAImax) * sin(zenith); // * cos(zenith);
			//cout << "LAI true (with cos(zenith)): " << LAItrue << endl;
			//cout << "sum_path_len: " << sum_path_len << endl;
			//total_LAI += sintheta * LAItrue;
			total_LAI += LAItrue;
			//total_LA += LAItrue * area;
			total_LA += FAVD * volumn_sum;
			//cout << "area & LAI: " << area << ' ' << LAItrue << endl;
			cout << "volumn & FAVD: " << volumn_sum << ' ' << FAVD << endl;
			//total_sintheta += sintheta;
			//fprintf(file_debug, "\n\n");
			fclose(file_debug);
		}

		
	}
	//total_LAI /= total_sintheta;
	//cout << "total height: " << total_height << endl;
	//cout << endl << "LAI: " << total_LAI << endl;
	cout << endl << "LA: " << total_LA << " total area: " << total_area << " total volumn: " << total_volumn << endl;

	FILE* file_result = _wfopen(L"E:\\result.txt", L"a");
	fprintf_s(file_result, "Process at %4d.%02d.%02d %02d:%02d:%02d \n", p.tm_year + 1900, p.tm_mon + 1, p.tm_mday, p.tm_hour, p.tm_min, p.tm_sec);

	fprintf(file_result, "Input (.ptx): %ls\n", ptx_path);
	fprintf(file_result, "      (envelope): %ls\n", envelope_path);
	//fprintf(file_result, "Output : %ls\n", out_path);
	fprintf(file_result, "Intensity threshold : %.2f\n", intensity_threshold);
	fprintf(file_result, "Skip : %d\n", g_n_skip);
	if (refine_envelope) fprintf(file_result, "-refine\n");
	if (cache_pathlen) fprintf(file_result, "-cache_pathlen\n");
	fprintf(file_result, "Number of divisions: %d\n",  n_divisions);
	if (use_gap_p){
		fprintf(file_result, "Path length: -ln(slide window gap probability)\n");
		fprintf(file_result, "Window width: %llfrad, height: %llfrad\n", window_width, window_height);
		}
	else
		fprintf(file_result, "Path length: beam length through crown\n");
	fprintf(file_result, "Zenith range: %llf - %llf, increment: %llf\n", zenith_low, zenith_high, zenith_increment);
	fprintf(file_result, "LA: %llf, total area: %llf, total volumn: %llf\n\n", total_LA, total_area, total_volumn);
	fclose(file_result);

	int pause;
	cin >> pause;
	file_debug = _wfopen(L"distribution.txt", L"a");
	fprintf(file_debug, "LA: %llf, total area: %llf, total volumn: %llf\n", total_LA, total_area, total_volumn);
	fprintf(file_debug, "Done\n\n");
	fclose(file_debug);
	// return 0;
	/*END-------------------------12/18/2019-----------------------------*/


	// 5. Calcuate within-canopy gap probability in zenith mask
	double gap_probability_within_canopy = GapWithinCrownProbabilityByZenith(within_canopy_image_new, all_gaps_image, zenith_mask_image);

	if (gap_probability_within_canopy < 0)
	{
		fprintf(file_out, "\nFAVD: %d\tGap probability: %d\n\n", -1, -1);
		fclose(file_out);
		return 3;
	}
	//cout << "gap_probability_within_canopy: " << gap_probability_within_canopy << endl;

	// 6. Get path length distribution within canopy
	/*Modified-------------------11/24/2019-----------------------------*/
	bool USE_PATH_LENGTH = true;
	ArrayXd data_distribution, ray_lengths_whinin_canopy, path_length_distribution;
	if (USE_PATH_LENGTH)
		data_distribution = DistributionByMask(path_lengths_image.cast<double>(), zenith_mask_image && within_canopy_image_new);
	else {
		//间隙率负对数分布
		data_distribution = DistributionByMask(gap_p.cast<double>(), zenith_mask_image && gap_p_valid_region);
	}
	//没用的东西，激光雷达到点云位置的分布
	// ray_lengths_whinin_canopy = DistributionByMask(ray_lengths_image.cast<double>(), zenith_mask_image && within_canopy_image_new);
	//第一种形式的路径长度分布
	path_length_distribution = DistributionByMask(path_lengths_image.cast<double>(), zenith_mask_image && within_canopy_image_new);


	//New: 2017-11-16 add image debug output
	if (img_debug)
	{

		//_mkdir(debug_dir);

		char img_path_norms2[_MAX_PATH];
		get_debug_img_path(img_path_norms2, ptx_path, envelope_path, debug_dir, "_norms2");
		write_debug_img(img_path_norms2, norms_image2);

		char img_path_gap[_MAX_PATH];
		get_debug_img_path(img_path_gap, ptx_path, envelope_path, debug_dir, "_gap");
		write_debug_img(img_path_gap, all_gaps_image);


		char img_path_canopy[_MAX_PATH];
		if (refine_envelope)
			get_debug_img_path(img_path_canopy, ptx_path, envelope_path, debug_dir, "_canopy_refined");
		else
			get_debug_img_path(img_path_canopy, ptx_path, envelope_path, debug_dir, "_canopy");
		write_debug_img(img_path_canopy, within_canopy_image_new);

		char img_path_pl[_MAX_PATH];
		get_debug_img_path(img_path_pl, ptx_path, envelope_path, debug_dir, "_path_length");
		write_debug_img(img_path_pl, path_lengths_image);

		char img_path_gap_1[_MAX_PATH];
		get_debug_img_path(img_path_gap_1, ptx_path, envelope_path, debug_dir, "_not_reach");
		write_debug_img(img_path_gap_1, points_not_reach_canopy_image);

		char img_path_gap_2[_MAX_PATH];
		get_debug_img_path(img_path_gap_2, ptx_path, envelope_path, debug_dir, "_pass_through");
		write_debug_img(img_path_gap_2, points_outside_canopy_image);
	}





	//cout << "path_lengths_min:" << data_distribution.minCoeff()  << endl;
	//cout << "path_lengths_max:" << data_distribution.maxCoeff()  << endl;

	if (g_n_bins_gMes > 0)
		G_value = GMean(g_Gfun_values, envelope_extent.zenith_min, envelope_extent.zenith_max, g_angle_interval_G);

	// 7. Calculate FAVD using path length distribution method
	double FAVD = CalcFAVD(gap_probability_within_canopy, data_distribution, G_value);
	/*START-----------------------11/24/2019-----------------------------*/
	//double max_raylength = CalcFAVD2(gap_probability_within_canopy, ray_lengths_whinin_canopy, G_value);
	double max_pathlength = CalcMaxPathLength(gap_probability_within_canopy, path_length_distribution, G_value);
	double FAVD_envelope = CalcFAVD(gap_probability_within_canopy, path_length_distribution, G_value);
	FAVD_envelope /= max_pathlength;
	//double FAVD2 = FAVD / max_raylength;
	double FAVD2 = FAVD / max_pathlength;
	/*END-----------------------11/24/2019-----------------------------*/


	//cout << "FAVD: " << FAVD << endl;

	fprintf(file_out, "\nZenith_Range\tFAVD\tGap_probability\tG_value\tNumber_of_pulses\tN_pulses_total\tSum_of_Path_Length\n");
	fprintf(file_out, "%.1lf-%.1lf\t%6.4lf\t%15.4lf\t%7.4lf\t%16lld\t%14lld\t%12.2lf\t%12.2lf\n", max(envelope_extent.zenith_min, zenith_low) / M_PI * 180, \
		min(envelope_extent.zenith_max, zenith_high) / M_PI * 180, \
		FAVD, gap_probability_within_canopy, G_value, \
		data_distribution.count(), \
		data_distribution.count() + points_not_reach_canopy_image.count(), \
		within_canopy_image_new.select(path_lengths_image, 0).sum(), \
		within_canopy_image.select(path_lengths_image, 0).sum());
	fprintf(file_out, "\nFAVD*maxpathlength\tFAVD\tFAVDenvelope\n");
	fprintf(file_out, "%6.4lf\t%6.4lf\t%6.4lf\n", FAVD, FAVD2, FAVD_envelope);
	fprintf(file_out, "\nFAVD calculated by slice\n");
	//fprintf(file_out, "%6.4lf\n", favd_ttl_cul);

	FILE* file_debug2 = _wfopen(L"debug.txt", L"w");
	fclose(file_debug2);
	if (b_zenith_layer)
	{
		VectorXd FAVDs(n_bins);
		for (size_t i = 0; i < n_bins; i++)
		{
			// zenith mask
			ArrayXXb zenith_layer_i = (zeniths_interpolated_image > zenith_ranges(i)) && (zeniths_interpolated_image < zenith_ranges(i + 1));
			file_debug2 = _wfopen(L"debug.txt", L"a");
			fprintf(file_debug2, "zenith range for slice: %lf to %lf\n", zenith_ranges(i) * 180 / 3.14159, zenith_ranges(i + 1) * 180 / 3.14159);
			fclose(file_debug2);
			//if ((zenith_layer_i * within_canopy).count() == 0) 
			//{
			//	fprintf(file_out, "%.1lf-%.1lf\t%d\t%d\n", zenith_ranges(i) / M_PI * 180, zenith_ranges(i + 1) / M_PI * 180, -1, -1);
			//	continue;
			//}
			double gap_probability_within_canopy_i = GapWithinCrownProbabilityByZenith(within_canopy_image_new, all_gaps_image, zenith_layer_i);

			ArrayXXb all_mask = within_canopy_image_new && zenith_layer_i;

			Index n_pulses_total = (points_not_reach_canopy_image && zenith_layer_i).count() + all_mask.count();

			if (gap_probability_within_canopy_i < -1e-5)	//no pulse
			{
				fprintf(file_out, "%.1lf-%.1lf\t%6d\t%15d\t%7d\t%16d\t%14lld\t%12.2lf\n", zenith_ranges(i) / M_PI * 180, zenith_ranges(i + 1) / M_PI * 180, -1, -1, -1, 0, n_pulses_total, 0.0f);
				continue;
			}
			else if (abs(gap_probability_within_canopy_i) < 1e-5) //Gap probability == 0
			{
				fprintf(file_out, "%.1lf-%.1lf\t%6d\t%15d\t%7d\t%16lld\t%14lld\t%12.2lf\n", zenith_ranges(i) / M_PI * 180, zenith_ranges(i + 1) / M_PI * 180, -1, 0, -1, \
					all_mask.count(), n_pulses_total, all_mask.select(path_lengths_image, 0).sum());
				continue;
			}
			else if (abs(gap_probability_within_canopy_i - 1) < 1e-5) //Gap probability == 1
			{
				fprintf(file_out, "%.1lf-%.1lf\t%6d\t%15d\t%7d\t%16lld\t%14lld\t%12.2lf\n", zenith_ranges(i) / M_PI * 180, zenith_ranges(i + 1) / M_PI * 180, 0, 1, -1, \
					all_mask.count(), n_pulses_total, all_mask.select(path_lengths_image, 0).sum());
				continue;
			}
			else
			{
				VectorXd data_distribution = DistributionByMask(path_lengths.cast<double>(), all_mask);
				if (g_n_bins_gMes > 0)
					G_value = GMean(g_Gfun_values, zenith_ranges(i), zenith_ranges(i + 1), g_angle_interval_G);
				FAVDs(i) = CalcFAVD(gap_probability_within_canopy_i, data_distribution, G_value);
				fprintf(file_out, "%.1lf-%.1lf\t%6.4lf\t%15.4lf\t%7.4lf\t%16lld\t%14lld\t%12.2lf\n", zenith_ranges(i) / M_PI * 180, zenith_ranges(i + 1) / M_PI * 180, FAVDs(i), \
					gap_probability_within_canopy_i, G_value, \
					data_distribution.count(), \
					n_pulses_total, \
					all_mask.select(path_lengths_image, 0).sum());
			}

		}
	}


	printf("\n");
	fprintf(file_out, "\n\n");
	fclose(file_out);
	return 0;
}

void get_debug_img_path(char* img_path, const wchar_t* ptx_path, const wchar_t* envelope_path, const char* debug_dir, const char* append)
{
	wchar_t fname[_MAX_FNAME];
	wchar_t fname_append[_MAX_FNAME];

	_wsplitpath_s(ptx_path, NULL, 0, NULL, 0, fname, _MAX_FNAME, NULL, 0);
	_wsplitpath_s(envelope_path, NULL, 0, NULL, 0, fname_append, _MAX_FNAME, NULL, 0);
	sprintf(img_path, "%s\\%ws_%ws%s.%s", debug_dir, fname, fname_append, append, "tif");
}

void write_debug_img(const char* img_path, Eigen::ArrayXXf& _image)
{
	FIBITMAP* dib = FreeImage_AllocateT(FIT_FLOAT, static_cast<int>(_image.rows()), static_cast<int>(_image.cols()), 32);
	BYTE* bits = FreeImage_GetBits(dib);
	memcpy(bits, (BYTE*)_image.data(), sizeof(_image(0, 0)) * _image.size());
	bool state = FreeImage_Save(FIF_TIFF, dib, img_path, TIFF_DEFAULT);
	FreeImage_Unload(dib);
}

void write_debug_img(const char* img_path, Map<ArrayXXf>& _image)
{
	FIBITMAP* dib = FreeImage_AllocateT(FIT_FLOAT, static_cast<int>(_image.rows()), static_cast<int>(_image.cols()), 32);
	BYTE* bits = FreeImage_GetBits(dib);
	memcpy(bits, (BYTE*)_image.data(), sizeof(_image(0, 0)) * _image.size());
	bool state = FreeImage_Save(FIF_TIFF, dib, img_path, TIFF_DEFAULT);
	FreeImage_Unload(dib);
}
/*
void write_debug_img(const char * img_path, ArrayXXb &_image)
{
	FIBITMAP* dib = FreeImage_AllocateT(FIT_BITMAP, static_cast<int>(_image.rows()), static_cast<int>(_image.cols()), 8);
	BYTE *bits = FreeImage_GetBits(dib);
	memcpy(bits, (BYTE*)_image.data(), sizeof(_image(0, 0)) * _image.size());
	dib = FreeImage_Threshold(dib, 1);
	bool state = FreeImage_Save(FIF_TIFF, dib, img_path, TIFF_DEFAULT);
	FreeImage_Unload(dib);
}
*/
void write_debug_img(const char* img_path, ArrayXXb& _image)
{
	FIBITMAP* dib = FreeImage_AllocateT(FIT_BITMAP, static_cast<int>(_image.rows()), static_cast<int>(_image.cols()), 8);
	BYTE* bits = FreeImage_GetBits(dib);
	if (_image.cols() % 4 == 0)
	{
		memcpy(bits, (BYTE*)_image.data(), sizeof(_image(0, 0)) * _image.size());
	}
	else
	{
		for (int y = 0; y < static_cast<int>(_image.cols()); y++) {
			BYTE* bits_line = (BYTE*)FreeImage_GetScanLine(dib, y);
			memcpy(bits_line, (BYTE*)_image.data() + static_cast<int>(_image.rows()) * y, sizeof(_image(0, 0)) * static_cast<int>(_image.rows()));
		}
	}

	dib = FreeImage_Threshold(dib, 1);
	bool state = FreeImage_Save(FIF_TIFF, dib, img_path, TIFF_DEFAULT);
	FreeImage_Unload(dib);
}



ArrayXXb RemoveEdgeGap(ArrayXXb within_canopy_image, ArrayXXb all_gaps_image)
{
	ArrayXXb _within_canopy_image_new = within_canopy_image;
	for (Index i = 0; i < all_gaps_image.rows(); i++)
	{
		for (Index j = 1; j < all_gaps_image.cols(); j++)
		{

			if (_within_canopy_image_new(i, j) && all_gaps_image(i, j) && (j == 1 || all_gaps_image(i, j - 1) || !_within_canopy_image_new(i, j - 1)))
			{
				_within_canopy_image_new(i, j) = false;
			}
			else if (_within_canopy_image_new(i, j) && !all_gaps_image(i, j))
			{
				break;
			}

		}

		for (size_t j = all_gaps_image.cols() - 2; j > 0; j--)
		{

			if (_within_canopy_image_new(i, j) && all_gaps_image(i, j) && (j == all_gaps_image.cols() - 2 || all_gaps_image(i, j + 1) || !_within_canopy_image_new(i, j + 1)))
			{
				_within_canopy_image_new(i, j) = false;
			}
			else if (_within_canopy_image_new(i, j) && !all_gaps_image(i, j))
			{
				break;
			}

		}

	}
	return _within_canopy_image_new;
}


list<wstring> GetFileNames(wchar_t* lpPath)
{
	list<wstring> fPaths;
	wchar_t file_path[_MAX_PATH];

	wchar_t drive[_MAX_DRIVE];
	wchar_t dir[_MAX_DIR];
	wchar_t fname[_MAX_FNAME];
	//wchar_t ext[_MAX_EXT];
	_wsplitpath(lpPath, drive, dir, fname, NULL);
	//_wmakepath_s(file_path, _MAX_PATH, drive, dir, fname, ext);

	wchar_t szFind[_MAX_PATH];
	WIN32_FIND_DATA FindFileData;
	wcscpy_s(szFind, _MAX_PATH, lpPath);
	HANDLE hFind = FindFirstFile((LPCTSTR)szFind, &FindFileData);
	if (INVALID_HANDLE_VALUE == hFind)    return fPaths;


	if (FindFileData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)//如果是目录
	{
		if (FindFileData.cFileName[0] != '.')//排除.和..文件夹
		{
			wcscat_s(szFind, _MAX_PATH, L"\\*.*");
			hFind = FindFirstFile((LPCTSTR)szFind, &FindFileData);
			if (INVALID_HANDLE_VALUE == hFind)    return fPaths;
		}
	}

	while (TRUE)
	{
		if (FindFileData.cFileName[0] != '.')
		{
			_wmakepath_s(file_path, _MAX_PATH, drive, dir, FindFileData.cFileName, NULL);
			fPaths.push_back(file_path);
		}

		if (!FindNextFile(hFind, &FindFileData))    break;

	}
	FindClose(hFind);
	return fPaths;
}


list<wstring> GetFileNames(wchar_t* lpPath, wchar_t* ext)
{
	list<wstring> fPaths;
	wchar_t file_path[_MAX_PATH];

	wchar_t drive[_MAX_DRIVE];
	wchar_t dir[_MAX_DIR];
	wchar_t fname[_MAX_FNAME];
	//wchar_t ext[_MAX_EXT];
	_wsplitpath(lpPath, drive, dir, fname, NULL);


	wchar_t szFind[_MAX_PATH];
	WIN32_FIND_DATA FindFileData;
	wcscpy_s(szFind, _MAX_PATH, lpPath);
	HANDLE hFind = FindFirstFile((LPCTSTR)szFind, &FindFileData);
	if (INVALID_HANDLE_VALUE == hFind)    return fPaths;


	if (FindFileData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)//如果是目录
	{
		if (FindFileData.cFileName[0] != '.')//排除.和..文件夹
		{
			wcscat_s(szFind, _MAX_PATH, L"\\*.");
			wcscat_s(szFind, _MAX_PATH, ext);
			hFind = FindFirstFile((LPCTSTR)szFind, &FindFileData);
			if (INVALID_HANDLE_VALUE == hFind)    return fPaths;
		}
	}
	else
	{
		_wmakepath_s(lpPath, _MAX_PATH, drive, dir, fname, ext);
	}


	while (TRUE)
	{
		if (FindFileData.cFileName[0] != '.')
		{
			_wmakepath_s(file_path, _MAX_PATH, drive, dir, FindFileData.cFileName, NULL);
			fPaths.push_back(file_path);
		}

		if (!FindNextFile(hFind, &FindFileData))    break;

	}
	FindClose(hFind);
	return fPaths;
}


wchar_t* output_path(wchar_t* new_path, const wchar_t* original_path, const wchar_t* fname_append, const wchar_t* new_ext)
{

	wchar_t drive[_MAX_DRIVE];
	wchar_t dir[_MAX_DIR];
	wchar_t fname[_MAX_FNAME];
	wchar_t ext[_MAX_EXT];


	//_wsplitpath_s(original_path, drive, dir, fname, ext);
	//wcscat_s(fname, _MAX_FNAME, fname_append);
	_wsplitpath_s(fname_append, drive, dir, fname, ext);

	if (new_ext != 0)
	{
		wcscpy_s(ext, _MAX_EXT, new_ext);
	}

	_wmakepath_s(new_path, _MAX_PATH, drive, dir, fname, ext);
	return new_path;
}


double GapWithinCrownProbabilityByZenith(ArrayXXb _within_canopy, ArrayXXb _is_gaps, ArrayXXb _zenith_layer)
{
	size_t n_canopy = (_zenith_layer * _within_canopy).count();
	if (n_canopy == 0)
	{
		return -1;
	}
	//cout << "(is_gaps * points_outside_canopy).count(): " << (_is_gaps * _points_outside_canopy).count() << endl;
	size_t n_canopy_gap = (_zenith_layer && _within_canopy && _is_gaps).count();
	double gap_probability_within_canopy_by_zenith = double(n_canopy_gap) / n_canopy;
	cout << "previous: n_gaps / n_rays: " << n_canopy_gap << ' ' << n_canopy << endl;

	return gap_probability_within_canopy_by_zenith;
}


ArrayXd  DistributionByMask(ArrayXXd _path_lengths, ArrayXXb _mask)
{
	ArrayXd _data_distribution(_path_lengths.size());
	//cout << "path_lengths.size(): " << _data_distribution.size() << endl;
	size_t ii = 0;
	for (Index i = 0; i < _path_lengths.size(); i++)
	{
		if (_path_lengths(i) > 0 && _mask(i))
		{
			_data_distribution(ii++) = _path_lengths(i);
		}
	}
	_data_distribution.conservativeResize(ii);
	return _data_distribution;
	//cout << "data_distribution.size(): " << ii << " " << data_distribution.size() << endl;
}