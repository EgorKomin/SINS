#pragma once
#include"consts.h"
#include "CalculateNavParams.h"
#include "CalculateNavParamsSpeedCorr.h"
#include"CalculateNavParamsWithoutCorr.h"
#include"CalculateNavParamsPositionCorr.h"





class SINS_algorithm
{

private:
	double q0, q1, q2, q3;
	bool apply_second_platform;

	CalculateNavParams* model_calculate;
	CalculateNavParams* model_base_calculate;
	CalculateNavParams* model_corr_calculate;
	
	double tcorr;
	double Ts;
	SINS_SOLUTION SOLinit;

	void multiply(double* q0, double* q1, double* q2, double* q3, double matrix[4][4], double a, double b, double c, double d);
	dataenup multiply(matrix matr, databody acc_body_takt);
	
	matrix getmatrix(double q0, double q1, double q2, double q3);
	orientationangles getangles(double q0, double q1, double q2, double q3);
	
	void CalculateRadius(double* Rfi, double *Rlmd, double fi_prev);

	void CalculateQuat1(databody gyro_takt, double* q0, double* q1, double* q2, double* q3);
	void CalculateQuat2(dataenup omega_takt, double* q0, double* q1, double* q2, double* q3);

	SINS_SOLUTION Navigation(CalculateNavParams* model_calculate, databody acc, databody gyro, SINS_SOLUTION SOL_prev, double* q0, double* q1, double* q2, double* q3, double dataCorr_E_lmd, double dataCorr_N_fi);
	
public:
	SINS_algorithm(double q0, double q1, double q2, double q3, typeOfCorrection typeCorr, double tcorr, double Ts, bool apply_second_platform);
	SOLUTIONvector start(const datavector& acc, const datavector& gyro, const doublevector& dataCorr_E_lmd, const doublevector& dataCorr_N_fi);
	void setInit (double Ve, double Vn, double Vup, double fi, double lmd, double Up);

};
