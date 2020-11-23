/*!
* \file LAIPath.h
* \date
*			2012/02/23 (standalone version)	Initial: Path length distribution method based on ellipse section assumption
*			2013/06/21 (standalone version)	New : Path length distribution inversed from measured gap data
*			2013/06/23 (GSL version)		Rewrite: Migrate to GSL library to improve speed
*			2016/06/19 (GSL version)		Modify for distribute
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

#include <gsl/gsl_histogram.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_statistics.h>
//#include "GFun.h"
#include "commdef.h"
#include "PTXreader.h"

#define LAI_MAX		10
#define NUM_BINS	25			//number of bins in histogram (not sensitive)

struct params_PathLen2GapF
{
	gsl_histogram *pathHist;
	double maxLAI;

};

struct params_GapBiasFromLAImax
{
	gsl_histogram *pathHist;
	double gapF;
};

struct dualParams 
{
	double par;
	double normalizedScale;
};


inline double neglog(double x) { return -log(x);};

double Stat_hist( double * data , unsigned long ndata, gsl_histogram *out_hist);

//路径长度为实测直方图
double LAI_PATH(gsl_histogram *pathHist, double gapFraction, double zenith = 0, double G = 0.5);
double solve_LAImax(gsl_histogram *pathHist, double gapFraction, double G);
double calc_LAItrue(gsl_histogram *pathHist,double zenith, double LAImax);
double LAI_PATH_FAVD(gsl_histogram *pathHist, double gapFraction, double max_pathlength, double G = 0.5);
double Func_PathProb(double _pathLen, void *params);
double Func_PathLen2GapF(double _pathLen, void *params);
double Func_WeightedPath(double _pathLen, void *params);
double Func_GapBiasFromLAImax(double x, void* params);

//路径长度为默认圆函数
double LAI_PATH_Circle(double gapFraction, double zenith = 0, double G = 0.5);
double LAIe2LAI_PATH_Circle(double effLAI, double zenith = 0, double G = 0.5);

inline double Func_PathProb_Circle(double _pathLen, void *params = 0) {return _pathLen / sqrt( 1- _pathLen * _pathLen );};
inline double Func_PathLen2GapF_Circle(double _pathLen, void *params){ return exp(-*(double *)params * _pathLen) * Func_PathProb_Circle(_pathLen);};
inline double Func_WeightedPath_Circle(double _pathLen, void *params = 0) {return _pathLen * Func_PathProb_Circle(_pathLen);};
double Func_GapBiasFromLAImax_Circle(double LAImax, void* params);

double P_of_the_window(Index i, Index j, ArrayXXb & all_gaps_image, ArrayXXb & canopy_region, Map<ArrayXXf> & path_lengths_image, double width, double height,
	Map<ArrayXXf> & azimuth_interpolated_image, Map<ArrayXXf> & zeniths_interpolated_image, double & sum_path_len);
double P_of_the_window_finite_avg(Index i, Index j, ArrayXXb& all_gaps_image, ArrayXXb& canopy_region, Map<ArrayXXf>& path_lengths_image, double width, double height,
	Map<ArrayXXf>& azimuth_interpolated_image, Map<ArrayXXf>& zeniths_interpolated_image, double& sum_path_len);
