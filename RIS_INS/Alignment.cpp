#include "Alignment.h"



Alignment::Alignment(databody avg_acc, databody avg_gyro, double fi)
{
	orientation = { 0 };
	this->avg_acc = avg_acc;
	this->avg_gyro = avg_gyro;
	this->fi = fi;
	C.resize(3, 3);
	align_completed = false;
}

void Alignment::startAlignment(Quat *quat)
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
	C(0, 0) = cos(roll) * cos(heading) + sin(pitch) * sin(roll) * sin(heading);
	C(0, 1) = cos(pitch) * sin(heading);
	C(0, 2) = sin(roll) * cos(heading) - sin(pitch) * cos(roll) * sin(heading);
	C(1, 0) = -cos(roll) * sin(heading) + sin(pitch) * sin(roll) * cos(heading);
	C(1, 1) = cos(pitch) * cos(heading);
	C(1, 2) = -sin(roll) * sin(heading) - sin(pitch) * cos(roll) * cos(heading);
	C(2, 0) = -cos(pitch) * sin(roll);
	C(2, 1) = sin(pitch);
	C(2, 2) = cos(pitch) * cos(roll);

	quat->q0 = sqrt(C(0, 0) + C(1, 1) + C(2, 2) + 1) / 2;

	if (quat->q0 == 0)
	{
		quat->q1 = sqrt((C(0, 0) + 1) / 2);
		quat->q2 = sqrt((C(1, 1) + 1) / 2);
		quat->q3 = sqrt((C(2, 2) + 1) / 2);
	}
	else
	{
		quat->q1 = -(C(1, 2) - C(2, 1)) / 4 / quat->q0;
		quat->q2 = -(C(2, 0) - C(0, 2)) / 4 / quat->q0;
		quat->q3 = -(C(0, 1) - C(1, 0)) / 4 / quat->q0;
	}

	//test
	//std::cout << *q0 - q0_test << std::endl;
	//std::cout << *q1 - q1_test << std::endl;
	//std::cout << *q2 - q2_test << std::endl;
	//std::cout << *q3 - q3_test << std::endl;

	consts::norm(quat);

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

matrix<double> Alignment::getMatrix()
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