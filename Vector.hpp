#ifndef VECTOR_H
#define VECTOR_H
#include <iostream>
#include <cmath>
#include <cstdlib>
#include "Complex.hpp"
template <typename type>
class Vector
{
    private:
	int length;
	type *array;
    public:
	Vector();
	explicit Vector(int length);
	Vector(type *array, int length);
	Vector(const Vector<type> &parameter);
	~Vector();
	double norm() const;
	void setLength(int length);
	int getLength() const;
	const Vector<type> & operator= (const Vector<type> &parameter);
	Vector<type> operator+ (const Vector<type> &parameter) const;
	Vector<type> operator- (const Vector<type> &parameter) const;
	void checkBounds(int i) const;
	type operator[] (int i) const;
	type & operator[] (int i);
	type operator* (const Vector<type> &parameter) const;
	Vector<type> operator* (type factor) const;
	Vector<type> scale(double factor) const;
	Vector<type> normalize() const;
	template<typename T>
	friend std::ostream & operator<< (std::ostream &os, const Vector<T> & parameter);
	template<typename T>
	friend Vector<T> operator* (T factor, const Vector<T> &parameter);
};

template<typename type>
Vector<type>::Vector()
{
    length = 1;
    array = new type[length];
    for (int i = 0; i < length; ++i)
	array[i] = 0;
}

template<typename type>
Vector<type>::Vector(int length)
{
    this->length = length;
    array = new type[length];
    for (int i = 0; i < length; ++i)
    {
	array[i] = 0;
    }
}

template<typename type>
Vector<type>::Vector(type *array, int length)
{
    this->length = length;
    this->array = new type[length];
    for (int i = 0; i < length; i++)
    {
	this->array[i] = array[i];
    }
}

template<typename type>
Vector<type>::Vector(const Vector<type> &parameter)
{
    this->length = parameter.length;
    array = new type[length];
    for (int i = 0; i < length; ++i)
    {
	array[i] = parameter.array[i];
    }
}

template<typename type>
Vector<type>::~Vector()
{
    delete []array;
}

template<typename type>
double Vector<type>::norm() const
{
    double sum = 0;
    for (int i = 0; i < length; ++i)
    {
	sum = sum + (array[i].conjugate() * array[i]).getReal();
    }
    return sqrt(sum);
}

template<>
double Vector<double>::norm() const
{
    double sum = 0;
    for (int i = 0;  i < length; ++i)
    {
	sum = sum + array[i]*array[i];
    }
    return sqrt(sum);
}

template <typename type>
void Vector<type>::setLength(int length)
{
    delete []array;
    this->length = length;
    array = new type[length];
    for (int i = 0; i < length; ++i)
    {
	array[i] = 0;
    }
}

template<typename type>
int Vector<type>::getLength() const
{
    return this->length;
}

template<typename type>
const Vector<type> & Vector<type>::operator= (const Vector<type> &parameter)
{
    delete []array;
    this->length = parameter.length;
    array = new type[length];
    for (int i = 0; i < length; ++i)
    {
	array[i] = parameter.array[i];
    }
    return *this;
}

template <typename type>
Vector<type> Vector<type>::operator+ (const Vector<type> &parameter) const
{
    if (this->length != parameter.length)
    {
	std::cout << "Vector lengths are incompatible. " << std::endl;
	std::exit(-1);
    }
    Vector<type> temp(length);
    for (int i = 0; i < length; ++i)
    {
	temp.array[i] = array[i] + parameter.array[i];
    }
    return temp;
}

template <typename type>
Vector<type> Vector<type>::operator- (const Vector<type> &parameter) const
{
    if (this->length != parameter.length)
    {
	std::cout << "Vector lengths are incompatible. " << std::endl;
	std::exit(-1);
    }
    Vector<type> temp(length);
    for (int i = 0; i < length; ++i)
    {
	temp.array[i] = array[i] - parameter.array[i];
    }
    return temp;
}

template<typename type>
void Vector<type>::checkBounds(int i) const
{
    if (i < 0 || i >= length)
    {
	std::cout << "Array index out of bounds. " << std::endl;
	std::exit(-1);
    }
}

template<typename type>
type Vector<type>::operator[] (int i) const
{
    this->checkBounds(i);
    return this->array[i];
}

template<typename type>
type & Vector<type>::operator[] (int i)
{
    this->checkBounds(i);
    return this->array[i];
}

template<typename type>
type Vector<type>::operator* (const Vector<type> &parameter) const
{
    type sum = 0;
    for (int i = 0; i < length; ++i)
    {
	sum = sum + array[i].conjugate() * parameter.array[i];
    }
    return sum;
}

template<>
double Vector<double>::operator* (const Vector<double> &parameter) const
{
    double sum = 0;
    for (int i = 0; i < length; i++)
    {
	sum = sum + array[i]*parameter.array[i];
    }
    return sum;
}

template<typename type>
Vector<type> Vector<type>::operator* (type factor) const
{
    Vector<type> temp(this->length);
    for (int i = 0; i < length; ++i)
    {
	temp.array[i] = factor*array[i];
    }
    return temp;
}

template<typename type>
Vector<type> Vector<type>::scale(double factor) const
{
    Vector<type> temp(this->length);
    for (int i = 0; i < length; ++i)
    {
	temp.array[i] = array[i]*factor;
    }
    return temp;
}

template<typename type>
Vector<type> Vector<type>::normalize() const
{
    return (*this).scale(1.0/this->norm());
}

template<typename type>
std::ostream & operator<< (std::ostream &os, const Vector<type> &parameter)
{
    for (int i = 0; i < parameter.length; ++i)
    {
	os << parameter.array[i] << "  ";
    }
    return os;
}

template <typename type>
Vector<type> operator* (type factor, const Vector<type> &parameter)
{
    int length = parameter.length;
    Vector<type> temp(length);
    for (int i = 0; i < length; ++i)
	temp.array[i] = factor*parameter.array[i];
    return temp;
}
#endif
