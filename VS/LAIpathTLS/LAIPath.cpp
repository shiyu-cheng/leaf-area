/*!
 * \file LAIPath.cpp
 * \date 
 *			2012/02/23 (standalone version)	Initial: Path length distribution method based on ellipse section assumption 
 *			2013/06/21 (standalone version)	New : Path length distribution inversed from measured gap data 
 *			2013/06/23 (GSL version)		Rewrite: Migrate to GSL library to improve speed   
 *			2016/06/19 (GSL version)		Distribute
 *
 * \author Ronghai Hu
 * Contact: rhu@unistra.fr
 *
 * \brief 
 *		"Path length distribution method for indirect leaf area index measurement"
 *
 * \note
 *		Theory: 
 *		Hu, R., Yan, G., Mu, X., & Luo, J. (2014). Indirect measurement of leaf area index 
 *		on the basis of path length distribution. REMOTE SENSING OF ENVIRONMENT, 155, 239-247. 
 *		doi: http://dx.doi.org/10.1016/j.rse.2014.08.032
 *
*/
#include "LAIPath.h"


//************************************
// Method:    LAI_PATH	Calculate Leaf Area Index (LAI) on the basis of Measured Path length distribution (GSL version)
// FullName:  LAI_PATH
// Access:    public 
// Returns:   double					True LAI of Path length distribution model
// Qualifier:
// Parameter: gsl_histogram * pathHist	Path length distribution (Histrgram format)
// Parameter: double gapFraction		Total gap fraction
// Parameter: double zenith				Zenith angle (degree) of data 
// Parameter: double G					Leaf projection function G
//************************************
double LAI_PATH(gsl_histogram *pathHist, double gapFraction, double zenith, double G)
{
	double effLAI = -log(gapFraction);
	struct params_GapBiasFromLAImax params = { pathHist, gapFraction };

	//1.Resolve LAImax
	int status;
	int iter = 0, max_iter = 100;

	double r_LAImax = 0, r_expected = effLAI;
	double x_lo = effLAI, x_hi = effLAI * 40;

	double f_lo = Func_GapBiasFromLAImax(x_lo, &params);
	double f_hi = Func_GapBiasFromLAImax(x_hi, &params);
	if ((f_lo < 0.0 && f_hi < 0.0) || (f_lo > 0.0 && f_hi > 0.0))
	{
		return LAI_MAX;
		printf("endpoints do not straddle y=0\n");
	}

	gsl_function F;
	F.function = &Func_GapBiasFromLAImax;

	F.params = &params;

	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc(T);
	status = gsl_root_fsolver_set(s, &F, x_lo, x_hi);

	do
	{
		iter++;
		status = gsl_root_fsolver_iterate(s);
		r_LAImax = gsl_root_fsolver_root(s);
		x_lo = gsl_root_fsolver_x_lower(s);
		x_hi = gsl_root_fsolver_x_upper(s);
		status = gsl_root_test_interval(x_lo, x_hi,
			0, 0.0001);
	} while (status == GSL_CONTINUE && iter < max_iter);

	gsl_root_fsolver_free(s);

	//2.Interation: lr*P(lr)
	F.function = &Func_WeightedPath;
	F.params = pathHist;

	gsl_integration_workspace * w
		= gsl_integration_workspace_alloc(1000);
	double integralWeightedPath, error;
	double *pts = pathHist->range;

	gsl_integration_qagp(&F, pts, pathHist->n + 1, 0, 1e-7, 1000,
		w, &integralWeightedPath, &error);

	gsl_integration_workspace_free(w);

	//Return true LAI
	return r_LAImax / G * integralWeightedPath * cos(zenith*M_PI / 180);
}


double solve_LAImax(gsl_histogram *pathHist, double gapFraction, double G)
{
	double effLAI = -log(gapFraction);
	struct params_GapBiasFromLAImax params = { pathHist, gapFraction };

	//1.Resolve LAImax
	int status;
	int iter = 0, max_iter = 100;

	double r_LAImax = 0, r_expected = effLAI;
	double x_lo = effLAI, x_hi = effLAI * 40;

	double f_lo = Func_GapBiasFromLAImax(x_lo, &params);
	double f_hi = Func_GapBiasFromLAImax(x_hi, &params);
	if ((f_lo < 0.0 && f_hi < 0.0) || (f_lo > 0.0 && f_hi > 0.0))
	{
		return LAI_MAX;
		printf("endpoints do not straddle y=0\n");
	}

	gsl_function F;
	F.function = &Func_GapBiasFromLAImax;

	F.params = &params;

	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc(T);
	status = gsl_root_fsolver_set(s, &F, x_lo, x_hi);

	do
	{
		iter++;
		status = gsl_root_fsolver_iterate(s);
		r_LAImax = gsl_root_fsolver_root(s);
		x_lo = gsl_root_fsolver_x_lower(s);
		x_hi = gsl_root_fsolver_x_upper(s);
		status = gsl_root_test_interval(x_lo, x_hi,
			0, 0.0001);
	} while (status == GSL_CONTINUE && iter < max_iter);

	gsl_root_fsolver_free(s);
	return r_LAImax / G;
}


double calc_LAItrue(gsl_histogram *pathHist, double zenith, double LAImax) {
	gsl_function F;
	F.function = &Func_WeightedPath;
	F.params = pathHist;

	gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);
	double integralWeightedPath, error;
	double *pts = pathHist->range;

	gsl_integration_qagp(&F, pts, pathHist->n + 1, 0, 1e-7, 1000,
		w, &integralWeightedPath, &error);

	gsl_integration_workspace_free(w);

	//Return true LAI
	//cout << "sin(zenith):" << sin(zenith) << "     cos(zenith):" << cos(zenith) << endl;
	//cout<<"LAImax * integralWeightedPath * sin(zenith):"<< LAImax * integralWeightedPath * sin(zenith)<<"     cos(zenith):" << LAImax * integralWeightedPath * cos(zenith) << endl;
	return LAImax * integralWeightedPath;
}


double LAI_PATH_FAVD(gsl_histogram *pathHist, double gapFraction, double max_pathlength, double G)
{
	double effLAI = -log(gapFraction);
	struct params_GapBiasFromLAImax params = { pathHist, gapFraction };

	//1.Resolve LAImax * 
	int status;
	int iter = 0, max_iter = 100;

	double r_LAImax = 0, r_expected = effLAI;
	double x_lo = effLAI, x_hi = effLAI * 40;

	double f_lo = Func_GapBiasFromLAImax(x_lo, &params);
	double f_hi = Func_GapBiasFromLAImax(x_hi, &params);
	if ((f_lo < 0.0 && f_hi < 0.0) || (f_lo > 0.0 && f_hi > 0.0))
	{
		return LAI_MAX;
		printf("endpoints do not straddle y=0\n");
	}

	gsl_function F;
	F.function = &Func_GapBiasFromLAImax;

	F.params = &params;

	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc(T);
	status = gsl_root_fsolver_set(s, &F, x_lo, x_hi);

	do
	{
		iter++;
		status = gsl_root_fsolver_iterate(s);
		r_LAImax = gsl_root_fsolver_root(s);
		x_lo = gsl_root_fsolver_x_lower(s);
		x_hi = gsl_root_fsolver_x_upper(s);
		status = gsl_root_test_interval(x_lo, x_hi,
			0, 0.0001);
	} while (status == GSL_CONTINUE && iter < max_iter);

	gsl_root_fsolver_free(s);


	//Return true LAI
	return r_LAImax / G / max_pathlength;
}

//************************************
// Method:    Calculate LAI on the basis of Path length distribution of ellipse section assumption (GSL version)
// Attention: Only used when path length distribution is unavailable (e.g. LAI-2000)
// FullName:  LAI_PATH_Circle
// Access:    public 
// Returns:   double				True LAI
// Qualifier:
// Parameter: double gapFraction	Total gap fraction
// Parameter: double zenith			Zenith angle (degree) of data 
// Parameter: double G				Leaf projection function G
//************************************
double LAI_PATH_Circle(double gapFraction, double zenith, double G)
{
	double effLAI = -log(gapFraction) / G * cos(zenith*M_PI / 180);

	//calculate the normalization coefficient
	gsl_function F;
	F.function = &Func_PathProb_Circle;

	gsl_integration_workspace * w
		= gsl_integration_workspace_alloc(1000);
	double normalizedScale, error;

	gsl_integration_qags(&F, 0, 1, 0, 1e-7, 1000,
		w, &normalizedScale, &error);

	gsl_integration_workspace_free(w);

	//1.Resolve LAImax
	int status;
	int iter = 0, max_iter = 100;
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
	double r_LAImax = 0, r_expected = effLAI;
	double x_lo = effLAI * G / cos(zenith*M_PI / 180), x_hi = x_lo * 20;

	F.function = &Func_GapBiasFromLAImax_Circle;

	struct dualParams params = { gapFraction,normalizedScale };
	F.params = &params;

	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc(T);
	gsl_root_fsolver_set(s, &F, x_lo, x_hi);

	do
	{
		iter++;
		status = gsl_root_fsolver_iterate(s);
		r_LAImax = gsl_root_fsolver_root(s);
		x_lo = gsl_root_fsolver_x_lower(s);
		x_hi = gsl_root_fsolver_x_upper(s);
		status = gsl_root_test_interval(x_lo, x_hi,
			0, 0.0001);

	} while (status == GSL_CONTINUE && iter < max_iter);

	gsl_root_fsolver_free(s);


	//2.Interation: lr*P(lr)
	F.function = &Func_WeightedPath_Circle;

	gsl_integration_workspace * w2
		= gsl_integration_workspace_alloc(1000);
	double integralWeightedPath;

	gsl_integration_qags(&F, 0, 1, 0, 1e-7, 1000,
		w2, &integralWeightedPath, &error);

	gsl_integration_workspace_free(w2);

	integralWeightedPath /= normalizedScale;

	//Return true LAI
	return r_LAImax * integralWeightedPath / G * cos(zenith*M_PI / 180);
}


//************************************
// Method:    Convert Effective LAI to true LAI on the basis of Path length distribution of ellipse section assumption (GSL version)
// Attention: Only used when path length distribution is unavailable (e.g. LAI-2000)
// FullName:  CalcTrueLAI_Circle
// Access:    public 
// Returns:   double				True LAI
// Qualifier:
// Parameter: double effLAI			Effective LAI (LAIe)
// Parameter: double zenith			Zenith angle (degree) of data 
// Parameter: double G				Leaf projection function G
//************************************
double LAIe2LAI_PATH_Circle(double effLAI, double zenith, double G)
{
	double gapFraction = exp(-effLAI * G / cos(zenith*M_PI / 180));
	
	return LAI_PATH_Circle(gapFraction, zenith, G);
}

//************************************
// Method:    Running statistics to get path length distribution
// FullName:  Stat_hist	
// Access:    public 
// Qualifier:
// Parameter: double* data			path length distribution 
// Parameter: DWORD ndata,			Number of path length distribution (size of path_length_data)
// Parameter: gsl_histogram *		path length distribution histogram (for return)
//************************************
double Stat_hist(double * path_length_data, unsigned long ndata, gsl_histogram *out_hist)
{
	// normalize path length to [0,1]
	double maxPathLen = gsl_stats_max(path_length_data, 1, ndata);
	for (unsigned long i = 0; i < ndata; i++)  path_length_data[i] /= maxPathLen;

	// obtain path length distribution
	gsl_histogram_set_ranges_uniform (out_hist, 0, 1);
	out_hist->range[out_hist->n] = 1 + 1e-6;
	for(unsigned long i = 0; i < ndata; i++)
		gsl_histogram_increment(out_hist, path_length_data[i]);

	// normalize total probability to 1
	gsl_histogram_scale(out_hist, static_cast<double>(out_hist->n) / ndata);

	return maxPathLen;
	//Modified: 2016-03-24 by rhhu  improve speed
}


//************************************
// Method:    Probability of a specific path length 
// FullName:  Func_PathProb			
// Access:    public 
// Returns:   double				Probability of path length = _pathLen
// Qualifier:
// Parameter: double _pathLen		path length
// Parameter: void * params			path length distribution histogram
//************************************
double Func_PathProb(double _pathLen, void *params)  
{
	gsl_histogram h = *(gsl_histogram*) params;
	size_t i;
	gsl_histogram_find(&h, _pathLen, &i);	//find the bin of _pathLen in histogram
	return gsl_histogram_get(&h, i);		//return its probability
}


//************************************
// Method:    Calculate gap fraction of a specific path length 
// FullName:  Func_PathLen2GapF
// Access:    public 
// Returns:   double				Gap fraction 
// Qualifier:
// Parameter: double _pathLen		Path length
// Parameter: void * params			Path length distribution histogram, LAImax
//************************************
double Func_PathLen2GapF(double _pathLen, void *params)
{
	struct params_PathLen2GapF *_params = (struct params_PathLen2GapF*)params;
	gsl_histogram * _pathHist = _params->pathHist;
	double _maxLAI = _params->maxLAI;
	return exp(-_maxLAI * _pathLen) * Func_PathProb(_pathLen, _pathHist);
}


//************************************
// Method:    lr*P(lr)
// FullName:  Func_WeightedPath
// Access:    public 
// Returns:   double				lr*P(lr)
// Qualifier:
// Parameter: double _pathLen		path length
// Parameter: void * params			Path length distribution histogram
//************************************
double Func_WeightedPath(double _pathLen, void *params)
{
	gsl_histogram *_pathHist = (gsl_histogram*)params;

	return _pathLen * Func_PathProb(_pathLen,_pathHist);
}


//************************************
// Method:    Func_GapBiasFromLAImax		cost function of a specific LAImax 
// FullName:  Func_GapBiasFromLAImax
// Access:    public 
// Returns:   double						Simulated gap fraction of LAImax and path length distribution - measured gap fraction
//											Eq(13) in (Hu et al., 2014)  - measured gap fraction
// Qualifier:
// Parameter: double LAImax					LAImax
// Parameter: void * params					Path length distribution histogram
//************************************
double Func_GapBiasFromLAImax(double LAImax, void* params)
{
	struct params_GapBiasFromLAImax _params = *(struct params_GapBiasFromLAImax*) params;
	
	double gapF = _params.gapF;
	struct params_PathLen2GapF par_F = {_params.pathHist, LAImax};

	gsl_function F;
	F.function = &Func_PathLen2GapF;	
	F.params = &par_F;

	//Interation: Equation(13) (Hu et al., 2014, RSE)
	gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (1000);
	double resGap, error;
	double *pts = _params.pathHist->range;

	gsl_integration_qagp(&F,pts,_params.pathHist->n + 1,0, 1e-7, 1000,
		w, &resGap, &error);

	gsl_integration_workspace_free (w);

	//���ؼ�϶�ʲ�ֵ
	return resGap - gapF;
}

//************************************
// Method:    cost function of a specific LAImax (ellipse section assumption)
// FullName:  Func_GapBiasFromLAImax_Circle
// Access:    public 
// Returns:   double						Simulated gap fraction of LAImax and path length distribution - measured gap fraction
// Qualifier:
// Parameter: double LAImax					LAImax
// Parameter: void * params					Path length distribution histogram
//************************************
double Func_GapBiasFromLAImax_Circle(double LAImax, void* params)
{
	//��ȡ�����̶�����
	struct dualParams _params = *(struct dualParams*) params;
	double gapF = _params.par;
	double normalizedScale = _params.normalizedScale;

	//���û��ַ���
	gsl_function F;
	F.function = &Func_PathLen2GapF_Circle;	
	F.params = &LAImax;

	//�������
	gsl_integration_workspace * w 
		= gsl_integration_workspace_alloc (1000);
	double resGap, error;

	gsl_integration_qags (&F, 0, 1, 0, 1e-7, 1000,
		w, &resGap, &error); 

	gsl_integration_workspace_free (w);

	resGap /= normalizedScale;
	//���ؼ�϶�ʲ�ֵ
	return resGap - gapF;
}


double P_of_the_window(Index i, Index j, ArrayXXb& all_gaps_image, ArrayXXb& canopy_region, Map<ArrayXXf>& path_lengths_image, double width, double height,
	Map<ArrayXXf>& azimuth_interpolated_image, Map<ArrayXXf>& zeniths_interpolated_image, double& sum_path_len)
{
	float MIN_RAYS = 0.1;
	if (width == 100) MIN_RAYS = 0;
	Index v = 0, u = 0;
	float n_gaps = 0;
	float n_rays = 0;
	float n_space = 0;
	//�м�ĵ�ĽǶȱȽ�׼
	double mid_col = all_gaps_image.cols() / 2;
	double mid_row = all_gaps_image.rows() / 2;

	double v_start = 0, u_start = 0;
	while (j - v_start > 0 && azimuth_interpolated_image(mid_row, j - v_start) < azimuth_interpolated_image(mid_row, j) + width / 2) v_start++;
	while (i - u_start > 0 && zeniths_interpolated_image(i - u_start, mid_col) < zeniths_interpolated_image(i, mid_col) + height / 2) u_start++;
	for (v = -v_start; j + v < all_gaps_image.cols() && azimuth_interpolated_image(mid_row, j + v) > azimuth_interpolated_image(mid_row, j) - width / 2; v++)
		for (u = -u_start; i + u < all_gaps_image.rows() && zeniths_interpolated_image(i + u, mid_col) > zeniths_interpolated_image(i, mid_col) - height / 2; u++) {
			// cout << u << ' ' << v << ' ' << zeniths_interpolated_image(i + u, j + v) << ' ' << azimuth_interpolated_image(i + u, j + v) << endl;
			if (canopy_region(i + u, j + v)) {
				if (all_gaps_image(i + u, j + v))
					n_gaps++;
				n_rays++;
				sum_path_len += path_lengths_image(i + u, j + v);
			}
			n_space++;
		}
	// cout << u << ' ' << v << ' ' << n_rays << ' ' << n_gaps << endl;
	//boundary remove
	if (n_rays < n_space * MIN_RAYS || n_rays < 10) {
		return -1;
	}
	else if (n_gaps > 0) {
		//��϶�ʷֲ�
		return n_gaps / n_rays;
	}
	else {
		// expand the window until there is a gap
		int edge = 1;
		while (1) {
			Index y_start = i - edge >= 0 ? i - edge : 0;
			Index y_end = i + u + edge <= all_gaps_image.rows() ? i + u + edge : all_gaps_image.rows();
			Index x_start = j - edge >= 0 ? j - edge : 0;
			Index x_end = j + v + edge <= all_gaps_image.cols() ? j + v + edge : all_gaps_image.cols();
			n_gaps = 0;
			n_rays = 0;
			for (Index y = y_start; y < y_end; y++)
				for (Index x = x_start; x < x_end; x++) {
					if (canopy_region(y, x)) {
						if (all_gaps_image(y, x))
							n_gaps++;
						n_rays++;
					}
					n_space++;
				}
			if (n_gaps == 0)
				edge++;
			else
				break;
		}
		if (n_rays < n_space * MIN_RAYS || n_rays < 10) {
			return -1;
		}
		return n_gaps / n_rays;
	}
}
	double P_of_the_window_finite_avg(Index i, Index j, ArrayXXb & all_gaps_image, ArrayXXb & canopy_region, Map<ArrayXXf> & path_lengths_image, double width, double height,
		Map<ArrayXXf> & azimuth_interpolated_image, Map<ArrayXXf> & zeniths_interpolated_image, double& sum_path_len)
	{
		float MIN_RAYS = 0.1;
		if (width == 100) MIN_RAYS = 0;
		Index v = 0, u = 0;
		float n_gaps = 0;
		float n_rays = 0;
		float n_space = 0;
		//�м�ĵ�ĽǶȱȽ�׼
		double mid_col = all_gaps_image.cols() / 2;
		double mid_row = all_gaps_image.rows() / 2;

		double v_start = 0, u_start = 0;
		while (j - v_start > 0 && azimuth_interpolated_image(mid_row, j - v_start) < azimuth_interpolated_image(mid_row, j) + width / 2) v_start++;
		while (i - u_start > 0 && zeniths_interpolated_image(i - u_start, mid_col) < zeniths_interpolated_image(i, mid_col) + height / 2) u_start++;
		for (v = -v_start; j + v < all_gaps_image.cols() && azimuth_interpolated_image(mid_row, j + v) > azimuth_interpolated_image(mid_row, j) - width / 2; v++)
			for (u = -u_start; i + u < all_gaps_image.rows() && zeniths_interpolated_image(i + u, mid_col) > zeniths_interpolated_image(i, mid_col) - height / 2; u++) {
				// cout << u << ' ' << v << ' ' << zeniths_interpolated_image(i + u, j + v) << ' ' << azimuth_interpolated_image(i + u, j + v) << endl;
				if (canopy_region(i + u, j + v)) {
					if (all_gaps_image(i + u, j + v))
						n_gaps++;
					n_rays++;
					sum_path_len += path_lengths_image(i + u, j + v);
				}
				n_space++;
			}
		// cout << u << ' ' << v << ' ' << n_rays << ' ' << n_gaps << endl;
		//boundary remove
		if (n_rays < n_space * MIN_RAYS || n_rays < 10) {
			return -1;
		}
		else if (n_gaps > 0) {
			//��϶�ʷֲ�
			return n_gaps / n_rays;
		}
		else {
			n_gaps = 1;
			if (n_rays < n_space * MIN_RAYS || n_rays < 10) {
				return -1;
			}
			return n_gaps / n_rays;
		}
}



