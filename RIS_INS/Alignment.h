#pragma once
#include"consts.h"
#include<iostream>

class Alignment
{
	matrix C;
	orientationangles orientation;
	databody avg_acc, avg_gyro;
	double fi;

	bool align_completed;

	void align_check();

public:
	Alignment(databody avg_acc, databody avg_gyro, double fi);
	
	void startAlignment(double* q0, double* q1, double* q2, double* q3);

	orientationangles getAlignmentError(databody deltaa, databody deltaomega);

	orientationangles getAngles();
	matrix getMatrix();
};

