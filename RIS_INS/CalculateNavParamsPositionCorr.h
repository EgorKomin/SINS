#pragma once
#include "CalculateNavParams.h"

class CalculateNavParamsPositionCorr :public CalculateNavParams
{
public:

	CalculateNavParamsPositionCorr(double Ts);

private:
	double FI_corr, LMD_corr;

	dataenup getAngularVelocity(double Ve, double Vn, double Rfi, double Rlmd, double fi, double lmd, double dataCorr_E_lmd, double dataCorr_N_fi) override;
	dataenup getLinearVelocity(dataenup omega_takt, dataenup acc_takt, dataenup velocity_prev, double fi_prev, double lmd_prev) override;
	SINS_SOLUTION getCoord(double Up_prev, double fi_prev, double lmd_prev, double Rfi, double Rlmd, dataenup velocity_takt) override;
};
