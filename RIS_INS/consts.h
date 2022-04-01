#pragma once
#include<cmath>
#include<vector>
#include<iostream>

#define PI 3.1416

enum typeOfCorrection
{
	WITHOUTCORR = 0,
	SPEEDCORR = 1,
	POSITIONCORR = 2
};

enum typeOfInfluence
{
	COMBINE = 2,
	ONLYSISTEM = 1,
	ONLYNOISE = 0
};

struct orientationangles
{
	double pitch, roll, heading;
};
struct matrix
{
	double c00, c01, c02, c10, c11, c12, c20, c21, c22;
};


struct databody
{
	double X, Y, Z;
};

struct dataenup
{
	double E, N, Up;
};

struct SINS_SOLUTION
{
	orientationangles angles;
	double Ve, Vn, Vup;
	double E, N, Up, fi, lmd;
	dataenup omega;
	matrix C;
};

typedef std::vector<SINS_SOLUTION> SOLUTIONvector;
typedef std::vector<databody> datavector;
typedef std::vector<double> doublevector;

struct conststr
{
	double f, T, R, u, g, a, b, e;
};


class consts
{
public:
	
	static double deg2rad(double angle)
	{
		return angle * PI / 180;
	}

	static double degPerHours2radPerSeconds(double angularVelocity)
	{
		return angularVelocity * PI / 180 / 3600;
	}


	static double rad2deg(double angle)
	{
		return angle * 180 / PI;
	}

	static conststr getconsts()
	{
		conststr str;
		str.f = 10;
		str.T = 3.0 * 3600;
		str.R = 6371 * pow(10, 3);
		str.u = degPerHours2radPerSeconds(15.0);
		str.g = 9.81;
		str.a = 6378245;
		str.b = 6356863;
		str.e = sqrt(1 - pow((str.b / str.a), 2));
		return str;
	}

	//функция нормировки кватерниона
	static void norm(double* q0, double* q1, double* q2, double* q3)
	{
		double n;
		n = sqrt(pow(*q0, 2) + pow(*q1, 2) + pow(*q2, 2) + pow(*q3, 2));
		*q0 = *q0 / n;
		*q1 = *q1 / n;
		*q2 = *q2 / n;
		*q3 = *q3 / n;

	}

};
