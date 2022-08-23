#pragma once
#include"consts.h"
#include "CalculateNavParams.h"
#include "CalculateNavParamsSpeedCorr.h"
#include"CalculateNavParamsWithoutCorr.h"
#include"CalculateNavParamsPositionCorr.h"





class SINS_algorithm
{

private:
	Quat quat_init;
	bool layout_base_platform, layout_second;

	CalculateNavParams* model_base_calculate;
	CalculateNavParams* model_corr_calculate;
	
	double tcorr;
	double Ts;

	double initEcoord, initNcoord;

	void UpdateQuat(Quat* quat, matrix<double> matr, double a, double b, double c, double d);
	dataenup acc2ENUp(matrix<double> matr, databody acc_body_takt);
	
	matrix<double> getmatrix(Quat quat);
	orientationangles getangles(Quat quat);
	
	void CalculateRadius(double* Rfi, double *Rlmd, double fi_prev);

	void CalculateQuat1(databody gyro_takt, Quat* quat);
	void CalculateQuat2(dataenup omega_takt, Quat* quat);

	void Navigation(CalculateNavParams* model_calculate, databody acc, databody gyro, double dataCorr_E_lmd, double dataCorr_N_fi);
	
public:
	SINS_algorithm(Quat quat_init, typeOfCorrection typeCorr, double tcorr, double Ts, bool apply_second_platform);
	void navigationPass(const databody acc, const databody gyro, const double dataCorr_E_lmd, const double dataCorr_N_fi);
	void setInit (double Ve, double Vn, double Vup, double fi, double lmd, double Up);
	SINS_SOLUTION getSolution();
};
