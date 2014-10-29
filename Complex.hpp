#ifndef COMPLEX_H
#define COMPLEX_H
#include <iostream>
#include <cmath>
#include <cstring>
class Complex
{
    private:
	double real;
	double imag;
    public:
	Complex();
	Complex(double real, double imag);
	Complex(double real);
	Complex(int parameter);
	Complex(const Complex &parameter);
	double norm() const;
	double argument() const;
	Complex conjugate() const;
	double getReal() const;
	double getImag() const;
	double doubleValue() const;
	const Complex& operator= (const Complex &parameter);
	const Complex& operator= (const double parameter);
	const Complex& operator= (const int parameter);
	bool operator== (const Complex &parameter) const;
	bool operator== (const double parameter) const;
	bool operator== (const int parameter) const;
	static Complex toComplex(double real, double imag);
	static Complex toComplex(double real);
	Complex operator+ (const Complex &parameter) const;
	Complex operator- (const Complex &parameter) const;
	Complex operator* (const Complex &parameter) const;
	Complex operator* (double factor) const;
	Complex operator/ (const Complex &parameter) const;
	Complex operator/ (double factor) const;
	std::string getTypeName() const;
	friend std::ostream & operator<< (std::ostream &os, const Complex &parameter);
	friend Complex exp(const Complex &parameter);
	friend double dble(const Complex &parameter);
	friend double dble(double a);
	friend double aimag(const Complex &parameter);
	friend Complex operator* (double factor, const Complex &parameter);
	friend Complex operator+ (double a, const Complex &parameter);
	friend Complex operator- (double a, const Complex &parameter);
	friend Complex operator- (const Complex &parameter);
	friend Complex operator/ (double a, const Complex &parameter);
};

Complex::Complex()
{
    this->real = 0;
    this->imag = 0;
}

Complex::Complex(double real, double imag)
{
    this->real = real;
    this->imag = imag;
}

Complex::Complex(double real)
{
    this->real = real;
    this->imag = 0;
}

Complex::Complex(int parameter)
{
    this->real = parameter;
    this->imag = 0;
}

Complex::Complex(const Complex &parameter)
{
    this->real = parameter.real;
    this->imag = parameter.imag;
}

double Complex::norm() const
{
    return sqrt(real*real + imag*imag);
}

double Complex::argument() const
{
    return atan2(imag, real);
}

Complex Complex::conjugate() const
{
    return Complex(real, -imag);
}

double Complex::getReal() const
{
    return this->real;
}

double Complex::getImag() const
{
    return this->imag;
}

double Complex::doubleValue() const
{
    return this->real;
}

const Complex& Complex::operator= (const Complex &parameter)
{
    this->real = parameter.real;
    this->imag = parameter.imag;
    return *this;
}

const Complex& Complex::operator= (const double parameter)
{
    this->real = parameter;
    this->imag = 0;
    return *this;
}

const Complex& Complex::operator= (const int parameter)
{
    this->real = parameter;
    this->imag = 0;
    return *this;
}

bool Complex::operator== (const Complex &parameter) const
{
    return (real == parameter.real) && (imag == parameter.imag);
}

bool Complex::operator== (const double parameter) const
{
    return (imag == 0)&&(real == parameter);
}

bool Complex::operator== (const int parameter) const
{
    return (imag == 0)&&(real == parameter);
}

Complex Complex::toComplex(double real, double imag)
{
    return Complex(real, imag);
}

Complex Complex::toComplex(double real)
{
    return Complex(real, 0.0);
}

Complex Complex::operator+ (const Complex &parameter) const
{
    return Complex(real + parameter.real, imag + parameter.imag);
}

Complex Complex::operator- (const Complex &parameter) const
{
    return Complex(real - parameter.real, imag - parameter.imag);
}

Complex Complex::operator* (const Complex &parameter) const
{
    double x, y;
    x = real*parameter.real - imag*parameter.imag;
    y = real*parameter.imag + imag*parameter.real;
    return Complex(x, y);
}

Complex Complex::operator* (double factor) const
{
    return Complex(real*factor, imag*factor);
}

Complex Complex::operator/ (const Complex &parameter) const
{
    Complex temp;
    temp = (*this) * (parameter.conjugate());
    temp = temp/(parameter.norm() * parameter.norm());
    return temp;
}

Complex Complex::operator/ (double factor) const
{
    return Complex(real/factor, imag/factor);
}

std::string Complex::getTypeName() const
{
    return std::string("complex");
}

std::ostream & operator<< (std::ostream &os, const Complex &parameter)
{
    os << "(" << parameter.real << ", " << parameter.imag << ")";
    return os;
}

Complex exp(const Complex &parameter)
{
    double x, y;
    x = parameter.real;
    y = parameter.imag;
    Complex temp(exp(x)*cos(y), exp(x)*sin(y));
    return temp;
}

double dble(const Complex &parameter)
{
    return parameter.real;
}

double dble(double a)
{
    return a;
}

double aimag(const Complex &parameter)
{
    return parameter.imag;
}

Complex operator* (double factor, const Complex &parameter)
{
    return parameter*factor;
}

Complex operator+ (double a, const Complex &parameter)
{
    return parameter + Complex(a, 0);
}

Complex operator- (double a, const Complex &parameter)
{
    return Complex(a - parameter.real, -parameter.imag);
}

Complex operator- (const Complex &parameter)
{
    return Complex(-parameter.real, -parameter.imag);
}

Complex operator/ (double a, const Complex &parameter)
{
    return Complex(a, 0)/parameter;
}
#endif
