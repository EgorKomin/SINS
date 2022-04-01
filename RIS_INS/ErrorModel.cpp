#include "ErrorModel.h"
#define DEBUG
//Модель ошибок реализована только для разомкнутой выставки!!!
ErrorModel::ErrorModel()
{
}


errorvector ErrorModel::geterror(const datavector& acc, const SOLUTIONvector& sol, double deltapitch, double deltaroll, double deltaheading, databody deltaa, databody deltaomega)
{
	double Fe0, Fn0, Fup0;
	double A[9][9];
	double E[9][9];
	double B[9][6];
	double matr[9][9];
	double matr1[9];
	double matr2[9];
	double U[6];
	double omegaE, omegaN, omegaUp, fi, Ve, Vn, Vup, heading, pitch;
	matrix C;
	auto str = consts::getconsts();
	auto dt = 1 / str.f;

	heading = sol[0].angles.heading;
	pitch = sol[0].angles.pitch;
	Fe0 = -(deltapitch * cos(heading) + deltaroll * cos(pitch) * sin(heading));
	Fn0 = deltapitch * sin(heading) - deltaroll * cos(pitch) * cos(heading);
	Fup0 = deltaheading + (Fe0 * sin(heading) + Fn0 * cos(heading)) * tan(pitch);

	Error err;
	err.Fe = Fe0;
	err.Fn = Fn0;
	err.Fup = Fup0;
	err.deltaVe = 0;
	err.deltaVn = 0;
	err.deltaVup = 0;
	err.deltafi = 0;
	err.deltalmd = 0;
	err.deltaE = 0;
	err.deltaN = 0;
	err.deltaUp = 0;
	err.deltapitch = deltapitch;
	err.deltaroll = deltaroll;
	err.deltaheading = deltaheading;

	U[0] = deltaa.X;
	U[1] = deltaa.Y;
	U[2] = deltaa.Z;
	U[3] = deltaomega.X;
	U[4] = deltaomega.Y;
	U[5] = deltaomega.Z;

	errorvector X;
	X.push_back(err);

	//нулевая  матрица А 9х9
	for (int i = 0; i < 9; i++)
	{
		for (int j = 0; j < 9; j++)
		{
			A[i][j] = 0;
		}

	}

	//нулевая  матрица B 9х6
	for (int i = 0; i < 9; i++)
	{
		for (int j = 0; j < 6; j++)
		{
			B[i][j] = 0;
		}

	}

	//единичная матрица Е 9х9
	for (int i = 0; i < 9; i++)
	{
		for (int j = 0; j < 9; j++)
		{
			if (i == j)
				E[i][j] = 1;
			else
				E[i][j] = 0;
		}
	}

	//либо берем только начальные значения(для случая неподвижного объекта)
#ifdef DEBUG
	omegaE = sol[0].omega.E;
	omegaN = sol[0].omega.N;
	omegaUp = sol[0].omega.Up;
	fi = sol[0].fi;
	Ve = sol[0].Ve;
	Vn = sol[0].Vn;
	Vup = sol[0].Vup;

	C.c00 = sol[0].C.c00;
	C.c01 = sol[0].C.c01;
	C.c02 = sol[0].C.c02;
	C.c10 = sol[0].C.c10;
	C.c11 = sol[0].C.c11;
	C.c12 = sol[0].C.c12;
	C.c20 = sol[0].C.c20;
	C.c21 = sol[0].C.c21;
	C.c22 = sol[0].C.c22;
#endif

	//цикл расчета ошибок
	for (size_t i = 1; i < acc.size(); i++)
	{
		//либо берем значения из алгоритма навигации
#ifndef DEBUG
		omegaE = sol[i].omega.E;
		omegaN = sol[i].omega.N;
		omegaUp = sol[i].omega.Up;
		fi = sol[i].fi;
		Ve = sol[i].Ve;
		Vn = sol[i].Vn;
		Vup = sol[i].Vup;

		C.c00 = sol[i].C.c00;
		C.c01 = sol[i].C.c01;
		C.c02 = sol[i].C.c02;
		C.c10 = sol[i].C.c10;
		C.c11 = sol[i].C.c11;
		C.c12 = sol[i].C.c12;
		C.c20 = sol[i].C.c20;
		C.c21 = sol[i].C.c21;
		C.c22 = sol[i].C.c22;
		heading = sol[i].angles.heading;
		pitch = sol[i].angles.pitch;
#endif

		//матрица А, остальные элементы нулевые
		A[0][1] = omegaUp;
		A[0][2] = -omegaN;
		A[0][4] = -1 / str.R;
		A[1][0] = -omegaUp;
		A[1][2] = omegaE;
		A[1][3] = 1 / str.R;
		A[1][6] = -str.u * sin(fi);
		A[2][0] = omegaN;
		A[2][1] = -omegaE;
		A[2][3] = 1 / str.R * tan(fi);
		A[2][6] = str.u * cos(fi) + Ve / str.R / pow(cos(fi), 2);
		A[3][1] = -(C.c20 * acc[i].X + C.c21 * acc[i].Y + C.c22 * acc[i].Z);
		A[3][2] = C.c10 * acc[i].X + C.c11 * acc[i].Y + C.c12 * acc[i].Z;
		A[3][3] = Vn / str.R * tan(fi);
		A[3][4] = Ve / str.R * tan(fi) + 2 * str.u * sin(fi);

		A[3][6] = (Ve / str.R / pow(cos(fi), 2) + 2 * str.u * cos(fi)) * Vn;
		A[4][0] = C.c20 * acc[i].X + C.c21 * acc[i].Y + C.c22 * acc[i].Z;
		A[4][2] = -(C.c00 * acc[i].X + C.c01 * acc[i].Y + C.c02 * acc[i].Z);
		A[4][3] = -2 * (Ve / str.R * tan(fi) + str.u * sin(fi));


		A[4][6] = -(Ve / str.R / pow(cos(fi), 2) + 2 * str.u * cos(fi)) * Ve;
		A[5][0] = -(C.c10 * acc[i].X + C.c11 * acc[i].Y + C.c12 * acc[i].Z);
		A[5][1] = C.c00 * acc[i].X + C.c01 * acc[i].Y + C.c02 * acc[i].Z;
		A[5][3] = 2 * (Ve / str.R + str.u * cos(fi));
		A[5][4] = 2 * Vn / str.R;
		A[5][6] = -2 * str.u * sin(fi) * Ve;
		A[5][8] = 2 * str.g / str.R;
		A[6][4] = 1 / str.R;
		A[7][3] = 1 / str.R / cos(fi);
		A[7][6] = Ve / str.R / cos(fi) * tan(fi);
		A[8][5] = 1;



		//матрица B, остальные элементы нулевые
		B[0][3] = -C.c00;
		B[0][4] = -C.c01;
		B[0][5] = -C.c02;

		B[1][3] = -C.c10;
		B[1][4] = -C.c11;
		B[1][5] = -C.c12;

		B[2][3] = -C.c20;
		B[2][4] = -C.c21;
		B[2][5] = -C.c22;

		B[3][0] = C.c00;
		B[3][1] = C.c01;
		B[3][2] = C.c02;

		B[4][0] = C.c10;
		B[4][1] = C.c11;
		B[4][2] = C.c12;

		B[5][0] = C.c20;
		B[5][1] = C.c21;
		B[5][2] = C.c22;

		//матрица (E+A*dt) 9х9
		for (int j = 0; j < 9; j++)
		{
			for (int k = 0; k < 9; k++)
			{
				matr[j][k] = E[j][k] + A[j][k] * dt;
			}
		}
		//вектор-столбец (E+A*dt)*X(n-1) 9x1
		for (int j = 0; j < 9; j++)
		{
			matr1[j] = matr[j][0] * X[i - 1].Fe + matr[j][1] * X[i - 1].Fn + matr[j][2] * X[i - 1].Fup + matr[j][3] * X[i - 1].deltaVe + matr[j][4] * X[i - 1].deltaVn + matr[j][5] * X[i - 1].deltaVup + matr[j][6] * X[i - 1].deltafi + matr[j][7] * X[i - 1].deltalmd + matr[j][8] * X[i - 1].deltaUp;
		}

		//вектор-столбец B*U 9x1
		for (int j = 0; j < 9; j++)
		{
			matr2[j] = B[j][0] * U[0] + B[j][1] * U[1] + B[j][2] * U[2] + B[j][3] * U[3] + B[j][4] * U[4] + B[j][5] * U[5];
		}

		//вектор-столбец B*U*dt 9x1
		for (int j = 0; j < 9; j++)
		{
			matr2[j] = matr2[j] * dt;
		}


		//результирующий вектор-столбец (E+A*dt)*X(n-1)+B*U*dt 9x1
		for (int j = 0; j < 9; j++)
		{
			matr2[j] = matr1[j] + matr2[j];
		}

		err.Fe = matr2[0];
		err.Fn = matr2[1];
		err.Fup = matr2[2];
		err.deltaVe = matr2[3];
		err.deltaVn = matr2[4];
		err.deltaVup = matr2[5];
		err.deltafi = matr2[6];
		err.deltalmd = matr2[7];
		err.deltaUp = matr2[8];
		err.deltaE = err.deltalmd * str.R * cos(fi);
		err.deltaN = err.deltafi * str.R;
		err.deltapitch = -(err.Fe * cos(heading) - err.Fn * sin(heading));
		err.deltaroll = -(err.Fn * cos(heading) + err.Fe * sin(heading)) / cos(pitch);
		err.deltaheading = err.Fup - (err.Fe * sin(heading) + err.Fn * cos(heading)) * tan(pitch);
		X.push_back(err);
	}

	return X;
}