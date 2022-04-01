#pragma once
#include"consts.h"

struct Error
{
	double Fe, Fn, Fup, deltaVe, deltaVn, deltaVup, deltafi, deltalmd, deltaE, deltaN, deltaUp, deltapitch, deltaroll, deltaheading;
};
typedef std::vector<Error> errorvector;
class ErrorModel
{
public:
	ErrorModel();
	errorvector geterror(const datavector& acc, const SOLUTIONvector& sol, double deltapitch, double deltaroll, double deltaheading, databody deltaa, databody deltaomega);
};

