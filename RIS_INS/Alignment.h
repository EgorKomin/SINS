#pragma once
#include"consts.h"
#include<iostream>

class Alignment
{
	matrix<double> C;
	orientationangles orientation;
	databody avg_acc, avg_gyro;
	double fi;

	bool align_completed;

	void align_check();

public:
	Alignment(databody avg_acc, databody avg_gyro, double fi);
	
	void startAlignment(Quat* quat);

	orientationangles getAlignmentError(databody deltaa, databody deltaomega);

	orientationangles getAngles();
	matrix<double> getMatrix();
};

