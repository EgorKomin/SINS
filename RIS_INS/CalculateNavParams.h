#pragma once
#include "consts.h"



class CalculateNavParams
{
protected:
	double k1, k2, k3;
	Quat quat; // образ гироплатформы
	SINS_SOLUTION sol;

public:
	void setQuat(Quat quat)
	{
		this->quat = quat;
	}

	Quat getQuat()
	{
		return quat;
	}

	void setSolution(SINS_SOLUTION sol)
	{
		this->sol = sol;
	}

	SINS_SOLUTION getSolution()
	{
		return sol;
	}
	
	CalculateNavParams()
	{
		quat = { 0 };
		k1 = 0;
		k2 = 0;
		k3 = 0;
	}

	CalculateNavParams(double Ts): CalculateNavParams()
	{
		double omegas = 2 * 3.14 / Ts;

		auto str = consts::getconsts();
		double nu2 = str.g / str.R;

		k1 = 1.75 * omegas;
		k2 = 2.15 * pow(omegas, 2) / nu2 - 1;
		k3 = pow(omegas, 3) / nu2 - 1.75 * omegas;
	}

	virtual ~CalculateNavParams()
	{

	}

	// модель вычислений по гироплатформе, реализуется в наследниках этого абстрактного класса
	virtual dataenup getAngularVelocity(double Ve, double Vn, double Rfi, double Rlmd, double fi, double lmd, double dataCorr_E_lmd, double dataCorr_N_fi) = 0;

	virtual dataenup getLinearVelocity(dataenup omega_takt, dataenup acc_takt, dataenup velocity_prev, double fi_prev, double lmd_prev) = 0;
	
	virtual SINS_SOLUTION getCoord(double Up_prev, double fi_prev, double lmd_prev, double Rfi, double Rlmd, dataenup velocity_takt) = 0;
};

