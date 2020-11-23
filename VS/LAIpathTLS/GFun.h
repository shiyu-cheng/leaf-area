#pragma once
// G����������� by ���ݺ� 2011.11.30
// ��ȡ��������	new: 2017-12-14 rhhu n_bins


#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif // !1



#include <cstdio>
#include <cmath>
using namespace std;




double gfun(double theta, double a, double b, double c);	//�����ηֲ���ʾg����
double gfun(double u1, double sigma, double v);			//Beta�ֲ���ʾg����
double gfun(double u1, double gd[], const int _n_bins);	  //���ݲ���ֵ����   new: 2017-12-14 rhhu n_bins

double AFun(double u1, double up);

double GFun(double up);									//���Ǻ������ηֲ�
double GFun(double up, double _a, double _b, double _c);	//���������κ����ֲ���������������G
double GFun(double up, double sigma, double v);			//����Beta�ֲ���������������G
double GFun(double up, double gd[], const int _n_bins);	//���ݲ���ֵ����	new: 2017-12-14 rhhu n_bins
double GMean(double * G, double _zenith_min, double _zenith_max, double _angle_interval);	//���ݲ���ֵ����	new: 2017-12-14 rhhu n_bins
int ReadMes(wchar_t *fMes, double *tmpg, const int _n_bins);	//��ȡ��������	new: 2017-12-14 rhhu n_bins

double gama(double x);

