#include "SINS_algorithm.h"

//конструктор класса
SINS_algorithm::SINS_algorithm(Quat quat_init, typeOfCorrection typeCorr, double tcorr, double Ts, bool layout_base_platform)
{
	this->quat_init = quat_init;

	model_base_calculate = new CalculateNavParamsWithoutCorr();
	model_base_calculate->setQuat(quat_init);
	model_corr_calculate = nullptr;
	
	
	this->layout_base_platform = layout_base_platform;
	this->tcorr = tcorr;
	this->Ts = Ts;
	layout_second = false;

	initEcoord = 0;
	initNcoord = 0;

	switch (typeCorr)
	{
	case WITHOUTCORR:
	{
		break;
	}
	case SPEEDCORR:
	{
		model_corr_calculate = new CalculateNavParamsSpeedCorr(Ts);
		break;
	}
	case POSITIONCORR:
	{
		model_corr_calculate = new CalculateNavParamsPositionCorr(Ts);
		break;
	}	
	default:
		break;
	}
	
}


void SINS_algorithm::navigationPass(const databody acc, const databody gyro, const double dataCorr_E_lmd, const double dataCorr_N_fi)
{
	
	auto str = consts::getconsts();
	double dt = 1 / str.f;
	

	static bool switched = false;

	static int step = 0;

	// НАВИГАЦИОННЫЙ АЛГОРИТМ

	// ВЫЧИСЛЕНИЯ ПО БАЗОВОЙ ПЛАТФОРМЕ
	Navigation(model_base_calculate, acc, gyro, NULL, NULL);

	// ПЕРЕКЛЮЧЕНИЕ НА КОРРЕКЦИЮ
	if ((step * dt > tcorr) && (model_corr_calculate != nullptr))
	{

		if (!switched)
		{
			// начинаем вычисления по новой платформе с последних вычислений по старой платформе
			switched = true;
			model_corr_calculate->setQuat(model_base_calculate->getQuat());
			model_corr_calculate->setSolution(model_base_calculate->getSolution());
		}


		// ВЫЧИСЛЕНИЯ ПО НОВОЙ ПЛАТФОРМЕ
		Navigation(model_corr_calculate, acc, gyro, dataCorr_E_lmd, dataCorr_N_fi);
		layout_second = true;

	}

	// во время переходного процесса продолжаем выводить показания базовой платформы
	if ((layout_base_platform) && (step * dt < tcorr + 4 * Ts))
	{
		layout_second = false;
	}

}

//функция установки начальных условий
void SINS_algorithm::setInit(double Ve, double Vn, double Vup, double fi, double lmd, double Up)
{
	auto str = consts::getconsts();

	SINS_SOLUTION SOLinit;
	SOLinit.angles = getangles(quat_init);
	SOLinit.C = getmatrix(quat_init);

	SOLinit.Ve = Ve;
	SOLinit.Vn = Vn;
	SOLinit.Vup = Vup;

	SOLinit.fi = fi;
	SOLinit.lmd = lmd;
	SOLinit.E = 0;
	SOLinit.N = 0;
	SOLinit.Up = Up;

	SOLinit.omega = model_base_calculate->getAngularVelocity(Ve, Vn, str.R, str.R, fi, lmd, NULL, NULL);

	double  Rfi, Rlmd;
	CalculateRadius(&Rfi, &Rlmd, SOLinit.fi);

	initEcoord = SOLinit.lmd * Rlmd * cos(SOLinit.fi);
	initNcoord = SOLinit.fi * Rfi;

	model_base_calculate->setSolution(SOLinit);
}

SINS_SOLUTION SINS_algorithm::getSolution()
{
	// ВЫВОД
	if (layout_second)
	{
		return model_corr_calculate->getSolution();
	}
	else
	{
		return model_base_calculate->getSolution();
	}
}

//функция обновления кватерниона
void SINS_algorithm::UpdateQuat(Quat* quat, matrix<double> matr, double a, double b, double c, double d)
{
	vector<double> col(4);
	col(0) = a;
	col(1) = b;
	col(2) = c;
	col(3) = d;

	vector<double> result = prod(matr, col);
	quat->q0 = result(0);
	quat->q1 = result(1);
	quat->q2 = result(2);
	quat->q3 = result(3);
}

//функция пересчета ускорений в ENUp
dataenup SINS_algorithm::acc2ENUp(matrix<double> matr, databody acc_body_takt)
{
	vector<double> acc_body(3); 
	acc_body(0) = acc_body_takt.X;
	acc_body(1) = acc_body_takt.Y;
	acc_body(2) = acc_body_takt.Z;

	vector<double> acc = prod(matr, acc_body);

	dataenup acc_takt;
	acc_takt.E = acc(0);
	acc_takt.N = acc(1);
	acc_takt.Up = acc(2);
	return acc_takt;
}

//функция получения матрицы направляющих косинусов
matrix<double> SINS_algorithm::getmatrix(Quat quat)
{
	double q0 = quat.q0;
	double q1 = quat.q1;
	double q2 = quat.q2;
	double q3 = quat.q3;

	matrix<double> MNK(3,3);
	MNK(0, 0) = pow(q0, 2) + pow(q1, 2) - pow(q2, 2) - pow(q3, 2);
	MNK(0, 1) = 2 * (q1 * q2 - q0 * q3);
	MNK(0, 2) = 2 * (q1 * q3 + q0 * q2);
	MNK(1, 0) = 2 * (q1 * q2 + q0 * q3);
	MNK(1, 1) = pow(q0, 2) - pow(q1, 2) + pow(q2, 2) - pow(q3, 2);
	MNK(1, 2) = 2 * (q2 * q3 - q0 * q1);
	MNK(2, 0) = 2 * (q1 * q3 - q0 * q2);
	MNK(2, 1) = 2 * (q2 * q3 + q0 * q1);
	MNK(2, 2) = pow(q0, 2) - pow(q1, 2) - pow(q2, 2) + pow(q3, 2);

	return MNK;
}

//функция получения углов ориентации
orientationangles SINS_algorithm::getangles(Quat quat)
{
	orientationangles angles;
	
	double q0 = quat.q0;
	double q1 = quat.q1;
	double q2 = quat.q2;
	double q3 = quat.q3;

	angles.pitch = atan(2 * (q2 * q3 + q0 * q1) / sqrt(4 * pow((q1 * q3 - q0 * q2), 2) + pow((pow(q0, 2) - pow(q1, 2) - pow(q2, 2) + pow(q3, 2)), 2)));
	angles.roll = -atan(2 * (q1 * q3 - q0 * q2) / (pow(q0, 2) - pow(q1, 2) - pow(q2, 2) + pow(q3, 2)));
	angles.heading = atan2(2 * (q1 * q2 - q0 * q3), (pow(q0, 2) - pow(q1, 2) + pow(q2, 2) - pow(q3, 2)));
	if (angles.heading < 0)
		angles.heading = angles.heading + 2 * PI;

	return angles;
}

void SINS_algorithm::CalculateRadius(double *Rfi, double *Rlmd, double fi_prev)
{
	auto str = consts::getconsts();

	*Rfi = str.R * (1 - pow(str.e, 2)) / pow(1 - pow(str.e, 2) * pow(sin(fi_prev), 2), (3 / 2));
	*Rlmd = str.R / pow(1 - pow(str.e, 2) * pow(sin(fi_prev), 2), (1 / 2));
}

void SINS_algorithm::CalculateQuat1(databody gyro_takt, Quat* quat)
{
	auto str = consts::getconsts();
	double dt = 1 / str.f;

	// ВЫЧИСЛЕНИЕ ПРИРАЩЕНИЙ УГЛОВ ЗА ТАКТ
	double Fx = gyro_takt.X * dt;
	double Fy = gyro_takt.Y * dt;
	double Fz = gyro_takt.Z * dt;


	double deltaF = sqrt(pow(Fx, 2) + pow(Fy, 2) + pow(Fz, 2));

	matrix<double> A(4, 4);

	// ВЫЧИСЛЕНИЕ КВАТЕРНИОНА - 1
	if (deltaF != 0)
	{

		double deltalmd0 = cos(deltaF / 2);
		double deltalmd1 = Fx / deltaF * sin(deltaF / 2);
		double deltalmd2 = Fy / deltaF * sin(deltaF / 2);
		double deltalmd3 = Fz / deltaF * sin(deltaF / 2);
		A(0, 0) = deltalmd0;
		A(0, 1) = -deltalmd1;
		A(0, 2) = -deltalmd2;
		A(0, 3) = -deltalmd3;
		A(1, 0) = deltalmd1;
		A(1, 1) = deltalmd0;
		A(1, 2) = deltalmd3;
		A(1, 3) = -deltalmd2;
		A(2, 0) = deltalmd2;
		A(2, 1) = -deltalmd3;
		A(2, 2) = deltalmd0;
		A(2, 3) = deltalmd1;
		A(3, 0) = deltalmd3;
		A(3, 1) = deltalmd2;
		A(3, 2) = -deltalmd1;
		A(3, 3) = deltalmd0;



		UpdateQuat(quat, A, quat->q0, quat->q1, quat->q2, quat->q3);

		// НОРМИРОВКА КВАТЕРНИОНА
		consts::norm(quat);
	}

}

void SINS_algorithm::CalculateQuat2(dataenup omega_takt, Quat* quat)
{
	auto str = consts::getconsts();
	double dt = 1 / str.f;


	double omegaxn = omega_takt.E;
	double omegayn = omega_takt.N;
	double omegazn = omega_takt.Up;

	double omegan = sqrt(pow(omegaxn, 2) + pow(omegayn, 2) + pow(omegazn, 2));

	matrix<double> D(4, 4);
	// ВЫЧИСЛЕНИЕ КВАТЕРНИОНА - 2 КОРРЕКЦИИ
	if (omegan != 0)
	{

		double deltam0 = cos(omegan * dt / 2);
		double deltam1 = -omegaxn / omegan * sin(omegan * dt / 2);
		double deltam2 = -omegayn / omegan * sin(omegan * dt / 2);
		double deltam3 = -omegazn / omegan * sin(omegan * dt / 2);

		D(0, 0) = quat->q0;
		D(0, 1) = -quat->q1;
		D(0, 2) = -quat->q2;
		D(0, 3) = -quat->q3;
		D(1, 0) = quat->q1;
		D(1, 1) = quat->q0;
		D(1, 2) = quat->q3;
		D(1, 3) = -quat->q2;
		D(2, 0) = quat->q2;
		D(2, 1) = -quat->q3;
		D(2, 2) = quat->q0;
		D(2, 3) = quat->q1;
		D(3, 0) = quat->q3;
		D(3, 1) = quat->q2;
		D(3, 2) = -quat->q1;
		D(3, 3) = quat->q0;


		UpdateQuat(quat, D, deltam0, deltam1, deltam2, deltam3);

		// НОРМИРОВКА КВАТЕРНИОНА
		consts::norm(quat);
	}
	
}

void SINS_algorithm::Navigation(CalculateNavParams* model_calculate, databody acc, databody gyro, double dataCorr_E_lmd, double dataCorr_N_fi)
{
	double  Rfi, Rlmd;
	dataenup acc_takt, velocity_prev, velocity_takt;
	SINS_SOLUTION coord;

	// получаем текущий кватернион и нав решение
	Quat quat = model_calculate->getQuat();
	SINS_SOLUTION SOL = model_calculate->getSolution();

	// ВЫЧИСЛЕНИЕ КВАТЕРНИОНА - 1
	CalculateQuat1(gyro, &quat);

	// ВЫЧИСЛЕНИЕ РАДИУСОВ ГЛАВНЫХ НОРМАЛЬНЫХ СЕЧЕНИЙ
	CalculateRadius(&Rfi, &Rlmd, SOL.fi);

	// ВЫЧИСЛЕНИЕ СКОРОСТЕЙ ENUp
	SOL.omega = model_calculate->getAngularVelocity(SOL.Ve, SOL.Vn, Rfi, Rlmd, SOL.fi, SOL.lmd, dataCorr_E_lmd, dataCorr_N_fi);

	// ВЫЧИСЛЕНИЕ КВАТЕРНИОНА - 2 КОРРЕКЦИИ
	CalculateQuat2(SOL.omega, &quat);

	// пушим обновленный кватернион
	model_calculate->setQuat(quat);

	// ВЫЧИСЛЕНИЕ МАТРИЦЫ ОРИЕНТАЦИИ
	SOL.C = getmatrix(quat);

	// ПЕРЕСЧЕТ УСКОРЕНИЙ В ENUP
	acc_takt = acc2ENUp(SOL.C, acc);

	// ВЫЧИСЛЕНИЕ УГЛОВ ОРИЕНТАЦИИ
	SOL.angles = getangles(quat);

	// ВЫЧИСЛЕНИЕ СКОРОСТЕЙ
	velocity_prev.E = SOL.Ve;
	velocity_prev.N = SOL.Vn;
	velocity_prev.Up = SOL.Vup;

	velocity_takt = model_calculate->getLinearVelocity(SOL.omega, acc_takt, velocity_prev, SOL.fi, SOL.lmd);
	SOL.Ve = velocity_takt.E;
	SOL.Vn = velocity_takt.N;
	SOL.Vup = velocity_takt.Up;

	// ВЫЧИСЛЕНИЕ КООРДИНАТ
	coord = model_calculate->getCoord(SOL.Up, SOL.fi, SOL.lmd, Rfi, Rlmd, velocity_takt);
	SOL.Up = coord.Up;
	SOL.fi = coord.fi;
	SOL.lmd = coord.lmd;

	SOL.E = SOL.lmd * Rlmd * cos(SOL.fi) - initEcoord;
	SOL.N = SOL.fi * Rfi - initNcoord;

	// пушим обновленное решение
	model_calculate->setSolution(SOL);
}

