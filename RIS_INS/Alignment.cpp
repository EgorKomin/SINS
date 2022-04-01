#include "Alignment.h"



Alignment::Alignment(databody avg_acc, databody avg_gyro, double fi)
{
	this->avg_acc = avg_acc;
	this->avg_gyro = avg_gyro;
	this->fi = fi;
	C = { 0 };
	orientation = { 0 };
	align_completed = false;
}

void Alignment::startAlignment(double *q0, double *q1, double *q2, double *q3)
{
	// разомкнута€ выставка
	auto str = consts::getconsts();

	dataenup omegaenup;
	omegaenup.Up = str.u * sin(fi);

	
	double pitch = atan(avg_acc.Y / sqrt(pow(str.g, 2) - pow(avg_acc.Y, 2)));
	double roll = -atan(avg_acc.X / avg_acc.Z);
	double heading = atan2(avg_acc.X * avg_gyro.Z - avg_acc.Z * avg_gyro.X, str.g * avg_gyro.Y - avg_acc.Y * omegaenup.Up);

	if (heading < 0)
	{
		heading = heading + 2.0 * PI;
	}
	//test
	/*double q0_test, q1_test, q2_test, q3_test;

	q0_test = cos(heading / 2) * cos(pitch / 2) * cos(roll / 2) + sin(heading / 2) * sin(pitch / 2) * sin(roll / 2);

	q1_test = cos(heading / 2) * sin(pitch / 2) * cos(roll / 2) + sin(heading / 2) * cos(pitch / 2) * sin(roll / 2);

	q2_test = cos(heading / 2) * cos(pitch / 2) * sin(roll / 2) - sin(heading / 2) * sin(pitch / 2) * cos(roll / 2);

	q3_test = cos(heading / 2) * sin(pitch / 2) * sin(roll / 2) - sin(heading / 2) * cos(pitch / 2) * cos(roll / 2);*/

	//матрица перехода из body в опорную ск
	C.c00 = cos(roll) * cos(heading) + sin(pitch) * sin(roll) * sin(heading);
	C.c01 = cos(pitch) * sin(heading);
	C.c02 = sin(roll) * cos(heading) - sin(pitch) * cos(roll) * sin(heading);
	C.c10 = -cos(roll) * sin(heading) + sin(pitch) * sin(roll) * cos(heading);
	C.c11 = cos(pitch) * cos(heading);
	C.c12 = -sin(roll) * sin(heading) - sin(pitch) * cos(roll) * cos(heading);
	C.c20 = -cos(pitch) * sin(roll);
	C.c21 = sin(pitch);
	C.c22 = cos(pitch) * cos(roll);

	*q0 = sqrt(C.c00 + C.c11 + C.c22 + 1) / 2;

	if (*q0 == 0)
	{
		*q1 = sqrt((C.c00 + 1) / 2);
		*q2 = sqrt((C.c11 + 1) / 2);
		*q3 = sqrt((C.c22 + 1) / 2);
	}
	else
	{
		*q1 = -(C.c12 - C.c21) / 4 / *q0;
		*q2 = -(C.c20 - C.c02) / 4 / *q0;
		*q3 = -(C.c01 - C.c10) / 4 / *q0;
	}

	//test
	//std::cout << *q0 - q0_test << std::endl;
	//std::cout << *q1 - q1_test << std::endl;
	//std::cout << *q2 - q2_test << std::endl;
	//std::cout << *q3 - q3_test << std::endl;

	consts::norm(q0, q1, q2, q3);

	orientation.heading = heading;
	orientation.pitch = pitch;
	orientation.roll = roll;

	align_completed = true;
}



orientationangles Alignment::getAlignmentError(databody deltaa, databody deltaomega)
{
	align_check();

	auto str = consts::getconsts();

	dataenup deltaaenup, deltaomegaenup, omegaenup;
	double heading = orientation.heading;

	omegaenup.N = str.u * cos(fi);

	deltaaenup.E = deltaa.X * cos(heading) + deltaa.Y * sin(heading);
	deltaomegaenup.E = deltaomega.X * cos(heading) + deltaomega.Y * sin(heading);

	orientationangles orientationError;

	orientationError.pitch = deltaa.Y / sqrt(pow(str.g, 2) - pow(avg_acc.Y, 2));
	orientationError.roll = (deltaa.Z * avg_acc.X - deltaa.X * avg_acc.Z) / (pow(avg_acc.X, 2) + pow(avg_acc.Z, 2));
	orientationError.heading = -deltaomegaenup.E / omegaenup.N + (deltaaenup.E / str.g) * tan(fi) - (deltaa.Z / str.g) * (sin(2 * heading) / 2);

	return orientationError;
}

orientationangles Alignment::getAngles()
{
	align_check();
	return orientation;
}

matrix Alignment::getMatrix()
{
	align_check();
	return C;
}

void Alignment::align_check()
{
	if (!align_completed)
	{
		std::cout << "¬ыставка не произведена" << std::endl;
		exit(1);
	}
}