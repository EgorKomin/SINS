#include<iostream>
#include"consts.h"
#include"SINS_algorithm.h"
#include<fstream>
#include<ctime>
#include"ErrorModel.h"
#include"Alignment.h"
#include"Noise.h"

#define DEBUG



void inputError(databody &deltaa, databody &deltaomega)
{
	auto str = consts::getconsts();

	std::cout << "введите дрейф акселерометра х в единицах mg" << std::endl;
	std::cin >> deltaa.X;
	deltaa.X = deltaa.X * str.g * pow(10, -3);

	std::cout << "введите дрейф акселерометра y в единицах mg" << std::endl;
	std::cin >> deltaa.Y;
	deltaa.Y = deltaa.Y * str.g * pow(10, -3);

	std::cout << "введите дрейф акселерометра z в единицах mg" << std::endl;
	std::cin >> deltaa.Z;
	deltaa.Z = deltaa.Z * str.g * pow(10, -3);

	std::cout << "введите дрейф гироскопа х в градусах в час" << std::endl;
	std::cin >> deltaomega.X;
	deltaomega.X = consts::degPerHours2radPerSeconds(deltaomega.X);


	std::cout << "введите дрейф гироскопа y в градусах в час" << std::endl;
	std::cin >> deltaomega.Y;
	deltaomega.Y = consts::degPerHours2radPerSeconds(deltaomega.Y);

	std::cout << "введите дрейф гироскопа z в градусах в час" << std::endl;
	std::cin >> deltaomega.Z;
	deltaomega.Z = consts::degPerHours2radPerSeconds(deltaomega.Z);
	std::cout << "=========================================" << std::endl;

}

void inputNoiseParamsIMU(double & bettaacc, double& sigmaacc, double& bettagyro, double& sigmagyro)
{
	auto str = consts::getconsts();

	std::cout << "введите бетта акселерометра" << std::endl;
	std::cin >> bettaacc;

	std::cout << "введите СКО акселерометра в mg" << std::endl;
	std::cin >> sigmaacc;
	sigmaacc = pow(10, -3) * str.g;

	std::cout << "введите бетта гироскопа" << std::endl;
	std::cin >> bettagyro;

	std::cout << "введите СКО гироскопа в градусах в час" << std::endl;
	std::cin >> sigmagyro;
	sigmagyro = consts::degPerHours2radPerSeconds(sigmagyro);

	std::cout << "=========================================" << std::endl;
}

//умножение матрицы 3х3 на вектор столбец
databody multiply(double C[3][3], dataenup obj)
{
	databody databody;
	databody.X = C[0][0] * obj.E + C[0][1] * obj.N + C[0][2] * obj.Up;
	databody.Y = C[1][0] * obj.E + C[1][1] * obj.N + C[1][2] * obj.Up;
	databody.Z = C[2][0] * obj.E + C[2][1] * obj.N + C[2][2] * obj.Up;
	return databody;
}

//поиск среднего значения массива за 10 мин
databody avg(const datavector& obj,int *size)
{ 
	auto str = consts::getconsts();
	double Tsr = 10.0 * 60;
	auto dt = 1 / str.f;
	double sum, sum1, sum2;
	
	sum = 0;
	sum1 = 0;
	sum2 = 0;
	databody databody;
	for (int i = 0; i < Tsr/dt; ++i)
	{
		sum += obj[i].X;
		sum1 += obj[i].Y;
		sum2 += obj[i].Z;
		*size = i;
	}
	databody.X = sum / *size;
	databody.Y = sum1 / *size;
	databody.Z = sum2 / *size;
	return databody;
}



int main()
{
#ifndef DEBUG
	srand(int(time(NULL)));
#endif // !DEBUG

	setlocale(LC_ALL, "Ru");

	double pitchetalon, rolletalon, headingetalon;

	double C[3][3];
	double fi, lmd, Ve, Vn, Vup, Up;

	double q0 = 1;
	double q1 = 0;
	double q2 = 0;
	double q3 = 0;
	
	auto str = consts::getconsts();
	auto dt = 1 / str.f;
	lmd = consts::deg2rad(37);
	Ve = 0;
	Vn = 0;
	Vup = 0;
	Up = 0;

	double bettagyro = 0;
	double sigmagyro = 0;

	double bettaacc = 0;
	double sigmaacc = 0;

	dataenup Fenup, omegaenup;
	databody Fxyz, omegaxyz;
	databody deltaa = { 0 };
	databody deltaomega = { 0 };
	
	Fenup.E = 0;
	Fenup.N = 0;
	Fenup.Up = str.g;
	

	

	std::cout << "введите угол тангажа в градусах" << std::endl;
	std::cin >> pitchetalon;
	pitchetalon = consts::deg2rad(pitchetalon);

	std::cout << "введите угол крена в градусах" << std::endl;
	std::cin >> rolletalon;
	rolletalon = consts::deg2rad(rolletalon);

	std::cout << "введите угол курса в градусах" << std::endl;
	std::cin >> headingetalon;
	headingetalon = consts::deg2rad(headingetalon);

	std::cout << "введите широту" << std::endl;
	std::cin >> fi;
	fi = consts::deg2rad(fi);


	omegaenup.E = 0;
	omegaenup.N = str.u * cos(fi);
	omegaenup.Up = str.u * sin(fi);

	
	int selectInfluence;
	std::cout << "Выберите тип ошибок ИИБ: Комбинированные-2/Только систем.ошибки-1/Только шум-0" << std::endl;
	std::cin >> selectInfluence;
	typeOfInfluence typeInfluence = (typeOfInfluence)selectInfluence;
	std::cout << "=========================================" << std::endl;



	doublevector sum_noise_gyro(int(str.T / dt), NULL), sum_noise_acc(int(str.T / dt), NULL);

	doublevector white_noise;

	Noise noise;

	switch (typeInfluence)
	{

	case COMBINE:
	{
		inputError(deltaa, deltaomega);

		white_noise = noise.getWhiteNoise(1);
		inputNoiseParamsIMU(bettaacc, sigmaacc, bettagyro, sigmagyro);
		sum_noise_gyro = noise.getColorNoise(white_noise, bettagyro, sigmagyro);
		sum_noise_acc = noise.getColorNoise(white_noise, bettaacc, sigmaacc);
		break;
	}
	case ONLYSISTEM:
	{
		inputError(deltaa, deltaomega);
		break;
	}
	case ONLYNOISE:
	{
		white_noise = noise.getWhiteNoise(1);
		inputNoiseParamsIMU(bettaacc, sigmaacc, bettagyro, sigmagyro);
		sum_noise_gyro = noise.getColorNoise(white_noise, bettagyro, sigmagyro);
		sum_noise_acc = noise.getColorNoise(white_noise, bettaacc, sigmaacc);
		break;
	}
	default:
		break;
	}


	int corr;
	std::cout << "Тип коррекции: Позиционная-2/Скоростная-1/Нет-0" << std::endl;
	std::cin >> corr;
	std::cout << "=========================================" << std::endl;

	typeOfCorrection typeCorr = (typeOfCorrection)corr;

	double tcorr = 0;
	double deltaCorr = 0;
	double Ts = 0;
	bool err = false;
	bool apply_second_platform = false;
	double sigmaSNS = 0;

	if (typeCorr != WITHOUTCORR)
	{
		std::cout << "введите время начала коррекции в минутах" << std::endl;
		std::cin >> tcorr;
		tcorr = tcorr * 60;

		std::cout << "введите постоянную времени скорректированной системы в минутах" << std::endl;
		std::cin >> Ts;
		Ts = Ts * 60;

		
		std::cout << "введите систематичекую ошибку корректора в м(м/с)" << std::endl;
		std::cin >> deltaCorr;

		std::cout << "введите СКО корректора в м(м/с)" << std::endl;
		std::cin >> sigmaSNS;

		std::cout << "подключать ли вторую платформу? Да-1/Нет-0" << std::endl;
		std::cin >> apply_second_platform;
		std::cout << "=========================================" << std::endl;
		
	}
	else
	{
		std::cout << "подключать ли модель ошибок(только без коррекции)? Да-1/Нет-0" << std::endl;
		std::cin >> err;
		std::cout << "=========================================" << std::endl;
	}
	
	// ФОРМИРОВАНИЕ ЭТАЛОНА И ОШИБОК СНС В ЗАВИСИМОСТИ ОТ ТИПА КОРРЕКЦИИ
	double etalonCorr_E_lmd = 0;
	double etalonCorr_N_fi = 0;

	double deltaCorr_E_lmd = 0;
	double deltaCorr_N_fi = 0;

	doublevector sum_noise_SNS_E_lmd(int(str.T / dt), NULL), sum_noise_SNS_N_fi(int(str.T / dt), NULL);

	switch (typeCorr)
	{
	case WITHOUTCORR:
		break;

	case SPEEDCORR:
	{
		etalonCorr_E_lmd = Ve;
		etalonCorr_N_fi = Vn;

		deltaCorr_E_lmd = deltaCorr;
		deltaCorr_N_fi = deltaCorr;

		sum_noise_SNS_E_lmd = noise.getWhiteNoise(sigmaSNS);
		sum_noise_SNS_N_fi = noise.getWhiteNoise(sigmaSNS);
		break;
	}
	case POSITIONCORR:
	{
		etalonCorr_E_lmd = lmd;
		etalonCorr_N_fi = fi;

		deltaCorr_N_fi = deltaCorr / str.R;
		deltaCorr_E_lmd = (deltaCorr + deltaCorr_N_fi * lmd * str.R * sin(fi)) / str.R / cos(fi);

		double sigmaSNS_N_fi = sigmaSNS / str.R;
		double sigmaSNS_E_lmd = (sigmaSNS + sigmaSNS_N_fi * lmd * str.R * sin(fi)) / str.R / cos(fi);

		sum_noise_SNS_E_lmd = noise.getWhiteNoise(sigmaSNS_E_lmd);
		sum_noise_SNS_N_fi = noise.getWhiteNoise(sigmaSNS_N_fi);
		break;
	}
	default:
		break;
	}

	// ФОРМИРОВАНИЕ ПОКАЗАНИЙ СНС
	doublevector dataCorr_E_lmd, dataCorr_N_fi;
	double dataCorr_E_lmd_takt, dataCorr_N_fi_takt;

	for (int i = 0; i < str.T / dt; ++i)
	{
		dataCorr_E_lmd_takt = etalonCorr_E_lmd + deltaCorr_E_lmd + sum_noise_SNS_E_lmd[i];
		dataCorr_N_fi_takt = etalonCorr_N_fi + deltaCorr_N_fi + sum_noise_SNS_N_fi[i];

		dataCorr_E_lmd.push_back(dataCorr_E_lmd_takt);
		dataCorr_N_fi.push_back(dataCorr_N_fi_takt);
	}


	bool calculateSKO = false;
	if ((typeCorr != WITHOUTCORR) && (typeInfluence != ONLYSISTEM))
	{
		calculateSKO = true;
	}

	//матрица перехода из опорной ск в body
	C[0][0] = cos(rolletalon)*cos(headingetalon) + sin(pitchetalon)*sin(rolletalon)*sin(headingetalon);
	C[0][1] = -cos(rolletalon)*sin(headingetalon) + sin(pitchetalon)*sin(rolletalon)*cos(headingetalon);
	C[0][2] = -cos(pitchetalon)*sin(rolletalon);
	C[1][0] = cos(pitchetalon)*sin(headingetalon);
	C[1][1] = cos(pitchetalon)*cos(headingetalon);
	C[1][2] = sin(pitchetalon);
	C[2][0] = sin(rolletalon)*cos(headingetalon) - sin(pitchetalon)*cos(rolletalon)*sin(headingetalon);
	C[2][1] = -sin(rolletalon)*sin(headingetalon) - sin(pitchetalon)*cos(rolletalon)*cos(headingetalon);
	C[2][2] = cos(pitchetalon)*cos(rolletalon);

	
	//показания акселерометров без ошибок
	Fxyz = multiply(C, Fenup);

	//показания гироскопов без ошибок
	omegaxyz = multiply(C, omegaenup);
	
	
	//формирование массива входных данных с акселерометров и гироскопов
	datavector acc, gyro;
	databody acc_takt, gyro_takt;

	for (int i = 0; i < str.T / dt; ++i)
	{
		acc_takt.X = Fxyz.X + deltaa.X + sum_noise_acc[i];
		acc_takt.Y = Fxyz.Y + deltaa.Y + sum_noise_acc[i];
		acc_takt.Z = Fxyz.Z + deltaa.Z + sum_noise_acc[i];

		gyro_takt.X = omegaxyz.X + deltaomega.X + sum_noise_gyro[i];
		gyro_takt.Y = omegaxyz.Y + deltaomega.Y + sum_noise_gyro[i];
		gyro_takt.Z = omegaxyz.Z + deltaomega.Z + sum_noise_gyro[i];
		
	
		acc.push_back(acc_takt);
		gyro.push_back(gyro_takt);
	}

	databody avg_acc, avg_gyro;

	if (typeInfluence != typeOfInfluence::ONLYSISTEM)
	{
		// осреднение за 10 мин
		int size = 0;

		avg_acc = avg(acc, &size);
		avg_gyro = avg(gyro, &size);

		// обрезаем показания, соответствующие первым 10 мин
		acc.erase(acc.cbegin(), acc.cbegin() + size);
		gyro.erase(gyro.cbegin(), gyro.cbegin() + size);

		dataCorr_E_lmd.erase(dataCorr_E_lmd.cbegin(), dataCorr_E_lmd.cbegin() + size);
		dataCorr_N_fi.erase(dataCorr_N_fi.cbegin(), dataCorr_N_fi.cbegin() + size);
	}
	else
	{
		avg_acc = acc[0];
		avg_gyro = gyro[0];
	}


	// ВЫСТАВКА
	Alignment align(avg_acc, avg_gyro, fi);
	// определение начального кватерниона
	align.startAlignment(&q0, &q1, &q2, &q3);

	// получение начальных углов ориентации
	auto angles = align.getAngles();

	std::cout << "углы ориентации после выставки, градусы" << std::endl;
	std::cout << consts::rad2deg(angles.pitch) << std::endl;
	std::cout << consts::rad2deg(angles.roll) << std::endl;
	std::cout << consts::rad2deg(angles.heading) << std::endl;
	std::cout << std::endl;

	std::cout << "ошибка в углах ориентации, градусы" << std::endl;
	std::cout << consts::rad2deg(angles.pitch - pitchetalon) << std::endl;
	std::cout << consts::rad2deg(angles.roll - rolletalon) << std::endl;
	std::cout << consts::rad2deg(angles.heading - headingetalon) << std::endl;
	std::cout << std::endl;

	// получение начальных ошибок ориентации
	auto deltaAngles = align.getAlignmentError(deltaa, deltaomega);

	std::cout << "ошибки в углах ориентации после выставки по формулам, градусы" << std::endl;
	std::cout << consts::rad2deg(deltaAngles.pitch) << std::endl;
	std::cout << consts::rad2deg(deltaAngles.roll) << std::endl;
	std::cout << consts::rad2deg(deltaAngles.heading) << std::endl;
	std::cout << std::endl;

	

	//формирование массива времени
	std::vector<double> t;
	double tcur = 0;

	for (size_t i = 0; i < acc.size(); ++i)
	{
		t.push_back(tcur);
		tcur = tcur + dt;
	}
		
	// ЗАПИСЬ КОНФИГУРАЦИИ
	std::cout << "creating output file <config>" << std::endl;
	std::ofstream file_cfg;
	file_cfg.open("config.txt");
	file_cfg <<
		(tcorr + 4 * Ts) / dt << "\n" <<	// Такт конца переходного процесса
		err << "\n" <<						// Флаг применения модели ошибок
		calculateSKO << std::endl;			// Флаг вычисления СКО

	file_cfg.close();


	//инициализация навигационного алгоритма
	SINS_algorithm alg(q0, q1, q2, q3, typeCorr, tcorr, Ts, apply_second_platform);

	//ввод начальных условий
	alg.setInit(Ve, Vn, Vup, fi, lmd, Up);

	std::cout << "working navigation algrorithm" << std::endl;
	auto SOLVEC = alg.start(acc, gyro, dataCorr_E_lmd, dataCorr_N_fi);
	

	// ЗАПИСЬ РЕЗУЛЬТАТОВ МОДЕЛИРОВАНИЯ В ФАЙЛ
	std::cout << "creating output file <solution>" << std::endl;
	std::ofstream file_sol;
	file_sol.open("solution.txt");

	for (size_t i = 0; i < SOLVEC.size(); ++i)
	{
		file_sol <<
			SOLVEC[i].angles.pitch - pitchetalon << "\t" <<
			SOLVEC[i].angles.roll - rolletalon << "\t" <<
			SOLVEC[i].angles.heading - headingetalon << "\t" <<
			SOLVEC[i].Ve << "\t" <<
			SOLVEC[i].Vn << "\t" <<
			SOLVEC[i].Vup << "\t" <<
			SOLVEC[i].E << "\t" <<
			SOLVEC[i].N << "\t" <<
			SOLVEC[i].Up << "\t" <<
			t[i] << std::endl;
	}
	file_sol.close();

	if (err)
	{
		//модель ошибок
		ErrorModel erralg;
		std::cout << "working error model" << std::endl;
		auto ERRVEC = erralg.geterror(acc, SOLVEC, deltaAngles.pitch, deltaAngles.roll, deltaAngles.heading, deltaa, deltaomega);

		std::cout << "creating output file <error>" << std::endl;
		std::ofstream file_err;
		file_err.open("error.txt");
		for (size_t i = 0; i < ERRVEC.size(); ++i)
		{
			file_err <<
				ERRVEC[i].deltapitch << "\t" << 
				ERRVEC[i].deltaroll << "\t" << 
				ERRVEC[i].deltaheading << "\t" << 
				ERRVEC[i].deltaVe << "\t" << 
				ERRVEC[i].deltaVn << "\t" << 
				ERRVEC[i].deltaVup << "\t" << 
				ERRVEC[i].deltaE << "\t" << 
				ERRVEC[i].deltaN << "\t" << 
				ERRVEC[i].deltaUp << "\t" << std::endl;
		}
		file_err.close();

	}
	return 0;
}