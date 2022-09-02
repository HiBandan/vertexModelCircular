#ifndef QUARTIC_H_INCLUDED
#define QUARTIC_H_INCLUDED
#include <Eigen/Dense>

#include <complex>

const double eps=1e-12;
const double M_2PI = 2*M_PI;

using namespace Eigen;

typedef std::complex<double> DComplex;

 inline DComplex polinom_2(DComplex x, double a, double b)
 {
	 return x * (x + a) + b;
 }

 inline DComplex polinom_3(DComplex x, double a, double b, double c)
 {
	 return x * (x * (x + a) + b) + c;
 }

 inline DComplex polinom_4(DComplex x, double a, double b, double c, double d)
 {
	 return x * (x * (x * (x + a) + b) + c) + d;
 }

unsigned int solveP3(double* x, double a, double b, double c);
DComplex* solve_quartic(double a, double b, double c, double d);

int HF(double);
int isPointBelowHorizontalLine(double,double);
int isPointInsideEllipse(double,double,double,double);
double DistancePoints(double,double,double,double);
double DistancePointEllipse (double,double,Array2d,Array2d&,int);
bool intersectionLineSegments(Vector2d,Vector2d,Vector2d,Vector2d,Vector2d&);

template <typename T> int sgn(T val)
{
    int sign = (T(0) < val) - (val < T(0));
    return sign == 0 ? 1 : sign;
}

#endif
