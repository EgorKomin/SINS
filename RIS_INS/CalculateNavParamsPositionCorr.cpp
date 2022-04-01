#include "CalculateNavParamsPositionCorr.h"



CalculateNavParamsPositionCorr::CalculateNavParamsPositionCorr(double Ts):CalculateNavParams(Ts)
{
	FI_corr = 0;
	LMD_corr = 0;	
}

dataenup CalculateNavParamsPositionCorr::getAngularVelocity(double Ve, double Vn, double Rfi, double Rlmd, double fi, double lmd, double dataCorr_E_lmd, double dataCorr_N_fi)
{
	FI_corr = dataCorr_N_fi;
	LMD_corr = dataCorr_E_lmd;

	auto str = consts::getconsts();
	dataenup omega_takt;
	
	omega_takt.E = -Vn / Rfi - k3 * (fi- FI_corr);
	omega_takt.N = Ve / Rlmd + str.u * cos(fi) + k3 * (lmd - LMD_corr)*cos(fi);
	omega_takt.Up = Ve / Rlmd * tan(fi) + str.u * sin(fi);
	return omega_takt;
}


dataenup CalculateNavParamsPositionCorr::getLinearVelocity(dataenup omega_takt, dataenup acc_takt, dataenup velocity_prev, double fi_prev, double lmd_prev)
{
	auto str = consts::getconsts();
	auto dt = 1 / str.f;

	auto omegaxn = omega_takt.E;
	auto omegayn = omega_takt.N;
	auto omegazn = omega_takt.Up;

	auto fe = acc_takt.E;
	auto fn = acc_takt.N;
	auto fup = acc_takt.Up;

	dataenup velocity_takt;

	velocity_takt.E = velocity_prev.E + (fe + (str.u * sin(fi_prev) + omegazn) * velocity_prev.N) * dt - k2 * str.g * (lmd_prev - LMD_corr) * cos(fi_prev) * dt;
	velocity_takt.N = velocity_prev.N + (fn - (str.u * sin(fi_prev) + omegazn) * velocity_prev.E) * dt - k2 * str.g * (fi_prev - FI_corr) * dt;
	velocity_takt.Up = velocity_prev.Up + (fup + (omegayn + str.u * cos(fi_prev)) * velocity_prev.E - omegaxn * velocity_prev.N - str.g) * dt;
	return velocity_takt;
}


SINS_SOLUTION CalculateNavParamsPositionCorr::getCoord(double Up_prev, double fi_prev, double lmd_prev, double Rfi, double Rlmd, dataenup velocity_takt)
{
	auto str = consts::getconsts();
	auto dt = 1 / str.f;

	SINS_SOLUTION coord;

	coord.Up = Up_prev + velocity_takt.Up * dt;
	coord.fi = fi_prev + velocity_takt.N / Rfi * dt - k1 * (fi_prev - FI_corr) * dt;
	coord.lmd = lmd_prev + velocity_takt.E / Rlmd / cos(coord.fi) * dt - k1 * (lmd_prev - LMD_corr) * dt;
	return coord;
}
