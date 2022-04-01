#include "SINS_algorithm.h"

//конструктор класса
SINS_algorithm::SINS_algorithm(double q0, double q1, double q2, double q3, typeOfCorrection typeCorr, double tcorr, double Ts, bool apply_second_platform)
{
	this->q0 = q0;
	this->q1 = q1;
	this->q2 = q2;
	this->q3 = q3;

	model_base_calculate = new CalculateNavParamsWithoutCorr();
	model_corr_calculate = nullptr;
	model_calculate = nullptr;
	
	SOLinit.angles = { 0 };
	SOLinit.C = { 0 };

	SOLinit.Ve = 0;
	SOLinit.Vn = 0;
	SOLinit.Vup = 0;

	SOLinit.omega = { 0 };

	SOLinit.fi = 0;
	SOLinit.lmd = 0;

	SOLinit.E = 0;
	SOLinit.N = 0;
	SOLinit.Up = 0;
	
	this->apply_second_platform = apply_second_platform;
	this->tcorr = tcorr;
	this->Ts = Ts;


	switch (typeCorr)
	{
	case WITHOUTCORR:
	{
		model_corr_calculate = model_base_calculate;
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


SOLUTIONvector SINS_algorithm::start(const datavector &acc, const datavector  &gyro, const doublevector &dataCorr_E_lmd, const doublevector &dataCorr_N_fi)
{
	
	
	double q0 = this->q0;
	double q1 = this->q1;
	double q2 = this->q2;
	double q3 = this->q3;

	double q0_second;
	double q1_second;
	double q2_second;
	double q3_second;
	
	auto str = consts::getconsts();
	auto dt = 1 / str.f;
	
	double  Rfi, Rlmd;

	SOLUTIONvector SOLVEC, SOLVEC_second, SOLVEC_layout;
	SOLVEC.push_back(SOLinit);
	SOLVEC_layout.push_back(SOLinit);

	SINS_SOLUTION SOL, SOL_second;

	CalculateRadius(&Rfi, &Rlmd, SOLinit.fi);

	SOLinit.E = SOLinit.lmd * Rlmd * cos(SOLinit.fi);
	SOLinit.N = SOLinit.fi * Rfi;

	bool switched = false;
	int iter_second = 0;

	// НАВИГАЦИОННЫЙ АЛГОРИТМ
	for (size_t i = 1; i < acc.size(); ++i)
	{

		// ПЕРЕКЛЮЧЕНИЕ НА КОРРЕКЦИЮ (ВЫБОР ПЛАТФОРМЫ)
		if (i * dt > tcorr)
		{
			model_calculate = model_corr_calculate;
		}
		else
		{
			model_calculate = model_base_calculate;	
		}

		
		// ВЫЧИСЛЕНИЯ ПО ВЫБРАННОЙ ПЛАТФОРМЕ
		SOL = Navigation(model_calculate, acc[i], gyro[i], SOLVEC[i - 1], &q0, &q1, &q2, &q3, dataCorr_E_lmd[i], dataCorr_N_fi[i]);
		SOLVEC.push_back(SOL);


		// ЗАПУСК РАСЧЕТА ПО ВТОРОЙ ПЛАТФОРМЕ (В ДАННОМ СЛУЧАЕ ВО ВРЕМЯ ПЕРЕХОДНОГО ПРОЦЕССА)
		if ((apply_second_platform) && (i * dt >= tcorr) && (i * dt <= tcorr + 4 * Ts))
		{
			iter_second++;
			if (!switched)
			{
				switched = true;
				q0_second = q0;
				q1_second = q1;
				q2_second = q2;
				q3_second = q3;
				SOLVEC_second.push_back(SOL);
			}

			// ВЫЧИСЛЕНИЯ ПО БАЗОВОЙ ПЛАТФОРМЕ
			SOL_second = Navigation(model_base_calculate, acc[i], gyro[i], SOLVEC_second[iter_second - 1], &q0_second, &q1_second, &q2_second, &q3_second, NULL, NULL);
			SOLVEC_second.push_back(SOL_second);

			SOLVEC_layout.push_back(SOL_second);
		}
		else
		{
			SOLVEC_layout.push_back(SOL);
		}
	}


	return SOLVEC_layout;
}

//функция обновления кватерниона
void SINS_algorithm::multiply(double *q0, double *q1, double *q2, double *q3, double matrix[4][4], double a, double b, double c, double d)
{
	*q0 = matrix[0][0] * a + matrix[0][1] * b + matrix[0][2] * c + matrix[0][3] * d;
	*q1 = matrix[1][0] * a + matrix[1][1] * b + matrix[1][2] * c + matrix[1][3] * d;
	*q2 = matrix[2][0] * a + matrix[2][1] * b + matrix[2][2] * c + matrix[2][3] * d;
	*q3 = matrix[3][0] * a + matrix[3][1] * b + matrix[3][2] * c + matrix[3][3] * d;
}

//функция пересчета ускорений в ENUp
dataenup SINS_algorithm::multiply(matrix matr, databody acc_body_takt)
{
	dataenup acc_takt;
	acc_takt.E = matr.c00 * acc_body_takt.X + matr.c01 * acc_body_takt.Y + matr.c02 * acc_body_takt.Z;
	acc_takt.N = matr.c10 * acc_body_takt.X + matr.c11 * acc_body_takt.Y + matr.c12 * acc_body_takt.Z;
	acc_takt.Up = matr.c20 * acc_body_takt.X + matr.c21 * acc_body_takt.Y + matr.c22 * acc_body_takt.Z;
	return acc_takt;
}

//функция установки начальных условий
void SINS_algorithm::setInit(double Ve, double Vn, double Vup, double fi, double lmd, double Up)
{
	auto str = consts::getconsts();

	SOLinit.angles = getangles(q0, q1, q2, q3);
	SOLinit.C = getmatrix(q0, q1, q2, q3);

	SOLinit.Ve = Ve;
	SOLinit.Vn = Vn;
	SOLinit.Vup = Vup;

	SOLinit.fi = fi;
	SOLinit.lmd = lmd;
	SOLinit.Up = Up;

	SOLinit.omega = model_base_calculate->getAngularVelocity(Ve, Vn, str.R, str.R, fi, lmd, NULL, NULL);
}

//функция получения матрицы направляющих косинусов
matrix SINS_algorithm::getmatrix(double q0, double q1, double q2, double q3)
{
	matrix MNK;
	MNK.c00 = pow(q0, 2) + pow(q1, 2) - pow(q2, 2) - pow(q3, 2);
	MNK.c01 = 2 * (q1 * q2 - q0 * q3);
	MNK.c02 = 2 * (q1 * q3 + q0 * q2);
	MNK.c10 = 2 * (q1 * q2 + q0 * q3);
	MNK.c11 = pow(q0, 2) - pow(q1, 2) + pow(q2, 2) - pow(q3, 2);
	MNK.c12 = 2 * (q2 * q3 - q0 * q1);
	MNK.c20 = 2 * (q1 * q3 - q0 * q2);
	MNK.c21 = 2 * (q2 * q3 + q0 * q1);
	MNK.c22 = pow(q0, 2) - pow(q1, 2) - pow(q2, 2) + pow(q3, 2);

	return MNK;
}

//функция получения углов ориентации
orientationangles SINS_algorithm::getangles(double q0, double q1, double q2, double q3)
{
	orientationangles angles;

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

void SINS_algorithm::CalculateQuat1(databody gyro_takt, double* q0, double* q1, double* q2, double* q3)
{
	auto str = consts::getconsts();
	double dt = 1 / str.f;

	// ВЫЧИСЛЕНИЕ ПРИРАЩЕНИЙ УГЛОВ ЗА ТАКТ
	double Fx = gyro_takt.X * dt;
	double Fy = gyro_takt.Y * dt;
	double Fz = gyro_takt.Z * dt;


	double deltaF = sqrt(pow(Fx, 2) + pow(Fy, 2) + pow(Fz, 2));

	double A[4][4];

	// ВЫЧИСЛЕНИЕ КВАТЕРНИОНА - 1
	if (deltaF != 0)
	{

		double deltalmd0 = cos(deltaF / 2);
		double deltalmd1 = Fx / deltaF * sin(deltaF / 2);
		double deltalmd2 = Fy / deltaF * sin(deltaF / 2);
		double deltalmd3 = Fz / deltaF * sin(deltaF / 2);
		A[0][0] = deltalmd0;
		A[0][1] = -deltalmd1;
		A[0][2] = -deltalmd2;
		A[0][3] = -deltalmd3;
		A[1][0] = deltalmd1;
		A[1][1] = deltalmd0;
		A[1][2] = deltalmd3;
		A[1][3] = -deltalmd2;
		A[2][0] = deltalmd2;
		A[2][1] = -deltalmd3;
		A[2][2] = deltalmd0;
		A[2][3] = deltalmd1;
		A[3][0] = deltalmd3;
		A[3][1] = deltalmd2;
		A[3][2] = -deltalmd1;
		A[3][3] = deltalmd0;



		multiply(q0, q1, q2, q3, A, *q0, *q1, *q2, *q3);

		// НОРМИРОВКА КВАТЕРНИОНА
		consts::norm(q0, q1, q2, q3);
	}

}

void SINS_algorithm::CalculateQuat2(dataenup omega_takt, double* q0, double* q1, double* q2, double* q3)
{
	auto str = consts::getconsts();
	double dt = 1 / str.f;


	double omegaxn = omega_takt.E;
	double omegayn = omega_takt.N;
	double omegazn = omega_takt.Up;

	double omegan = sqrt(pow(omegaxn, 2) + pow(omegayn, 2) + pow(omegazn, 2));

	double  D[4][4];
	// ВЫЧИСЛЕНИЕ КВАТЕРНИОНА - 2 КОРРЕКЦИИ
	if (omegan != 0)
	{

		double deltam0 = cos(omegan * dt / 2);
		double deltam1 = -omegaxn / omegan * sin(omegan * dt / 2);
		double deltam2 = -omegayn / omegan * sin(omegan * dt / 2);
		double deltam3 = -omegazn / omegan * sin(omegan * dt / 2);

		D[0][0] = *q0;
		D[0][1] = -*q1;
		D[0][2] = -*q2;
		D[0][3] = -*q3;
		D[1][0] = *q1;
		D[1][1] = *q0;
		D[1][2] = *q3;
		D[1][3] = -*q2;
		D[2][0] = *q2;
		D[2][1] = -*q3;
		D[2][2] = *q0;
		D[2][3] = *q1;
		D[3][0] = *q3;
		D[3][1] = *q2;
		D[3][2] = -*q1;
		D[3][3] = *q0;


		multiply(q0, q1, q2, q3, D, deltam0, deltam1, deltam2, deltam3);

		// НОРМИРОВКА КВАТЕРНИОНА
		consts::norm(q0, q1, q2, q3);
	}
	
}

SINS_SOLUTION SINS_algorithm::Navigation(CalculateNavParams* model_calculate, databody acc, databody gyro, SINS_SOLUTION SOL_prev, double *q0, double* q1, double* q2, double* q3, double dataCorr_E_lmd, double dataCorr_N_fi)
{
	double  Rfi, Rlmd;
	dataenup omega_takt, acc_takt, velocity_prev, velocity_takt;
	SINS_SOLUTION SOL, coord;

	// ВЫЧИСЛЕНИЕ КВАТЕРНИОНА - 1
	CalculateQuat1(gyro, q0, q1, q2, q3);

	// ВЫЧИСЛЕНИЕ РАДИУСОВ ГЛАВНЫХ НОРМАЛЬНЫХ СЕЧЕНИЙ
	CalculateRadius(&Rfi, &Rlmd, SOL_prev.fi);

	// ВЫЧИСЛЕНИЕ СКОРОСТЕЙ ENUp
	omega_takt = model_calculate->getAngularVelocity(SOL_prev.Ve, SOL_prev.Vn, Rfi, Rlmd, SOL_prev.fi, SOL_prev.lmd, dataCorr_E_lmd, dataCorr_N_fi);
	SOL.omega = omega_takt;

	// ВЫЧИСЛЕНИЕ КВАТЕРНИОНА - 2 КОРРЕКЦИИ
	CalculateQuat2(omega_takt, q0, q1, q2, q3);

	// ВЫЧИСЛЕНИЕ МАТРИЦЫ ОРИЕНТАЦИИ
	SOL.C = getmatrix(*q0, *q1, *q2, *q3);

	// ПЕРЕСЧЕТ УСКОРЕНИЙ В ENUP
	acc_takt = multiply(SOL.C, acc);

	// ВЫЧИСЛЕНИЕ УГЛОВ ОРИЕНТАЦИИ
	SOL.angles = getangles(*q0, *q1, *q2, *q3);

	// ВЫЧИСЛЕНИЕ СКОРОСТЕЙ
	velocity_prev.E = SOL_prev.Ve;
	velocity_prev.N = SOL_prev.Vn;
	velocity_prev.Up = SOL_prev.Vup;

	velocity_takt = model_calculate->getLinearVelocity(omega_takt, acc_takt, velocity_prev, SOL_prev.fi, SOL_prev.lmd);
	SOL.Ve = velocity_takt.E;
	SOL.Vn = velocity_takt.N;
	SOL.Vup = velocity_takt.Up;

	// ВЫЧИСЛЕНИЕ КООРДИНАТ
	coord = model_calculate->getCoord(SOL_prev.Up, SOL_prev.fi, SOL_prev.lmd, Rfi, Rlmd, velocity_takt);
	SOL.Up = coord.Up;
	SOL.fi = coord.fi;
	SOL.lmd = coord.lmd;

	SOL.E = SOL.lmd * Rlmd * cos(SOL.fi) - SOLinit.E;
	SOL.N = SOL.fi * Rfi - SOLinit.N;

	return SOL;
}

