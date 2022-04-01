#include "CalculateNavParamsWithoutCorr.h"


dataenup CalculateNavParamsWithoutCorr::getAngularVelocity(double Ve, double Vn, double Rfi, double Rlmd, double fi, double lmd, double dataCorr_E_lmd, double dataCorr_N_fi)
{
	auto str = consts::getconsts();
	dataenup omegatakt;
	omegatakt.E = -Vn / Rfi;
	omegatakt.N = Ve / Rlmd + str.u * cos(fi);
	omegatakt.Up = Ve / Rlmd * tan(fi) + str.u * sin(fi);
	return omegatakt;
}


dataenup CalculateNavParamsWithoutCorr::getLinearVelocity(dataenup omega_takt, dataenup acc_takt, dataenup velocity_prev, double fi_prev, double lmd_prev)
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

	velocity_takt.E = velocity_prev.E + (fe + (str.u * sin(fi_prev) + omegazn) * velocity_prev.N) * dt;
	velocity_takt.N = velocity_prev.N + (fn - (str.u * sin(fi_prev) + omegazn) * velocity_prev.E) * dt;
	velocity_takt.Up = velocity_prev.Up + (fup + (omegayn + str.u * cos(fi_prev)) * velocity_prev.E - omegaxn * velocity_prev.N - str.g) * dt;
	return velocity_takt;
}


SINS_SOLUTION CalculateNavParamsWithoutCorr::getCoord(double Up_prev, double fi_prev, double lmd_prev, double Rfi, double Rlmd, dataenup velocity_takt)
{
	auto str = consts::getconsts();
	auto dt = 1 / str.f;

	SINS_SOLUTION coord;

	coord.Up = Up_prev + velocity_takt.Up * dt;
	coord.fi = fi_prev + velocity_takt.N / Rfi * dt;
	coord.lmd = lmd_prev + velocity_takt.E / Rlmd / cos(coord.fi) * dt;
	return coord;
}