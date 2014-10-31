#ifndef MATRIX_H
#define MATRIX_H
#include <iostream>
#include <cstdlib>
#include "Vector.hpp"
template <typename type>
class Matrix
{
    private:
	int dimension;
	type **matrix;
	type determinant(Matrix<type> &upperTriangle) const;	
	void checkBounds(int i) const;
    public:
	Matrix();
	explicit Matrix(int dimension);
	Matrix(type **array, int dimension);
	Matrix(const Matrix<type> &parameter);
	~Matrix();
	int getDimension() const;
	void setDimension(int dimension);
	const Matrix<type> & operator= (const Matrix<type> &parameter);
	Matrix<type> operator+ (const Matrix<type> &parameter) const;
	Matrix<type> operator- (const Matrix<type> &parameter) const;
	Matrix<type> operator* (const Matrix<type> &parameter) const;
	Matrix<type> operator* (type factor) const;
	type operator() (int i, int j) const;
	type& operator() (int i, int j);
	void QRDecomposition(Matrix<type> &Q, Matrix<type> &R) const;
	type determinant() const;
	Matrix<type> inverse() const;
	Matrix<type> Hermite() const;
	Matrix<type> QRIteration() const;
	type cofactor(int i, int j) const;
	Matrix<type> adjugate() const;
	template <typename T>
	friend std::ostream & operator<< (std::ostream &os, const Matrix<T> &parameter);
};

template <typename type>
Matrix<type>::Matrix()
{
    dimension = 2;
    matrix = new type *[dimension];
    for (int i = 0; i < dimension; ++i)
    {
	matrix[i] = new type[dimension];
    }
    for (int i = 0; i < dimension; ++i)
    {
	for (int j = 0; j < dimension; ++j)
	{
	    matrix[i][j] = 0;
	}
    }
}

template <typename type>
Matrix<type>::Matrix(int dimension)
{
    this->dimension = dimension;
    matrix = new type*[dimension];
    for (int i = 0; i < dimension; ++i)
    {
	matrix[i] = new type[dimension];
    }
    for (int i = 0; i < dimension; ++i)
    {
	for (int j = 0; j < dimension; ++j)
	{
	    matrix[i][j] = 0;
	}
    }
}

template <typename type>
Matrix<type>::Matrix(type **array, int dimension)
{
    this->dimension = dimension;
    matrix = new type *[this->dimension];
    for (int i = 0; i < dimension; ++i)
    {
	matrix[i] = new type[dimension];
    }

    for (int i = 0; i < dimension; ++i)
    {
	for (int j = 0; j < dimension; ++j)
	{
	    matrix[i][j] = array[i][j];
	}
    }
}

template <typename type>
Matrix<type>::Matrix(const Matrix<type> &parameter)
{
    for (int i = 0; i < this->dimension; ++i)
    {
	delete []matrix[i];
    }
    delete []matrix;

    this->dimension = parameter.dimension;
    matrix = new type*[dimension];
    for (int i = 0 ; i < dimension; ++i)
    {
	matrix[i] = new type[dimension];
    }

    for (int i = 0; i < dimension; ++i)
    {
	for (int j = 0; j < dimension; ++j)
	{
	    this->matrix[i][j] = parameter.matrix[i][j];
	}
    }
}

template <typename type>
Matrix<type>::~Matrix()
{
    for (int i = 0; i < dimension; ++i)
    {
	delete []matrix[i];
    }
    delete []matrix;
}

template <typename type>
int Matrix<type>::getDimension() const
{
    return this -> dimension;
}

template <typename type>
void Matrix<type>::setDimension(int dimension)
{
    for (int i = 0; i < this->dimension; ++i)
    {
	delete []matrix[i];
    }
    delete []matrix;

    this->dimension = dimension;
    matrix = new type*[dimension];
    for (int i = 0; i < dimension; ++i)
    {
	matrix[i] = new type[dimension];
    }

    for (int i = 0; i < dimension; ++i)
    {
	for (int j = 0; j < dimension; ++j)
	{
	    matrix[i][j] = 0;
	}
    }
}

template <typename type>
const Matrix<type>& Matrix<type>::operator= (const Matrix<type> &parameter)
{
    for (int i = 0; i < this->dimension; ++i)
    {
	delete []matrix[i];
    }
    delete []matrix;

    this->dimension = parameter.dimension;
    matrix = new type*[dimension];
    for (int i = 0; i < dimension; ++i)
    {
	matrix[i] = new type[dimension];
    }

    for (int i = 0; i < dimension; ++i)
    {
	for (int j = 0; j < dimension; ++j)
	{
	    matrix[i][j] = parameter.matrix[i][j];
	}
    }
    return *this;
}

template <typename type>
Matrix<type> Matrix<type>::operator+ (const Matrix<type> &parameter) const
{
    if (this->dimension != parameter.dimension)
    {
	std::cout << "Matrix dimensions are imcompatible. " << std::endl;
	std::exit(-1);
    }

    Matrix<type> temp(dimension);
    for (int i = 0; i < dimension; ++i)
    {
	for (int j = 0; j < dimension; ++j)
	{
	    temp.matrix[i][j] = this->matrix[i][j] + parameter.matrix[i][j];
	}
    }
    return temp;
}

template <typename type>
Matrix<type> Matrix<type>::operator- (const Matrix<type> &parameter) const
{
    if (this->dimension != parameter.dimension)
    {
	std::cout << "Matrix dimensions are incompatible. " << std::endl;
	std::exit(-1);
    }

    Matrix<type> temp(dimension);
    for (int i = 0; i < dimension; ++i)
    {
	for (int j = 0; j < dimension; ++j)
	{
	    temp.matrix[i][j] = this->matrix[i][j] - parameter.matrix[i][j];
	}
    }
    return temp;
}

template <typename type>
Matrix<type> Matrix<type>::operator* (const Matrix<type> &parameter) const
{
    if (this->dimension != parameter.dimension)
    {
	std::cout << "Matrix dimensions are incompatible. " << std::endl;
	std::exit(-1);
    }

    Matrix<type> temp(dimension);
    for (int i = 0; i < dimension; ++i)
    {
	for (int j = 0; j < dimension; ++j)
	{
	    temp.matrix[i][j] = 0;
	    for (int k = 0; k < dimension; ++k)
	    {
		temp.matrix[i][j] = temp.matrix[i][j] + this->matrix[i][k]*parameter.matrix[k][j];
	    }
	}
    }
    return temp;
}

template <typename type>
Matrix<type> Matrix<type>::operator* (type factor) const
{
    Matrix<type> temp(dimension);
    for (int i = 0; i < dimension; ++i)
    {
	for (int j = 0; j < dimension; ++j)
	{
	    temp.matrix[i][j] = this->matrix[i][j]*factor;
	}
    }
    return temp;
}

template<typename type>
void Matrix<type>::checkBounds(int i) const
{
    if (i < 0 || i >= dimension)
    {
	std::cout << "Index out of bounds. " << std::endl;
	std::exit(-1);
    }
}

template <typename type>
type Matrix<type>::operator() (int i, int j)const
{
    this->checkBounds(i);
    this->checkBounds(j);

    return this->matrix[i][j];
}

template <typename type>
type & Matrix<type>::operator() (int i, int j)
{
    this->checkBounds(i);
    this->checkBounds(j);

    return this->matrix[i][j];
}

template <typename type>
void Matrix<type>::QRDecomposition(Matrix<type> &Q, Matrix<type> &R) const
{
    Vector<type> *e;
    Vector<type> *xi;
    e = new Vector<type>[dimension];
    xi = new Vector<type>[dimension];
    for (int i = 0; i < dimension; ++i)
    {
	e[i].setLength(dimension);
	xi[i].setLength(dimension);
    }
    for (int col = 0; col < dimension; ++col)
    {
	for (int row = 0; row < dimension; ++row)
	{
	    xi[col][row] = matrix[row][col];
	}
    }

    e[0] = xi[0].normalize();
    for (int col = 1; col < dimension; ++col)
    {
    	Vector<type> projection(dimension);
	for (int i = 0; i < dimension; ++i)
	    projection[i] = 0;
	for (int i = 0; i < col; ++i)
	{
	    projection = projection + e[i]*(e[i]*xi[col]);
	}
	e[col] = (xi[col] - projection).normalize();
    }

    for (int row = 0; row < dimension; ++row)
    {
	for (int col = 0; col < dimension; ++col)
	{
	    Q(row, col) = e[col][row];
	}
    }

    for (int row = 0; row < dimension; ++row)
    {
	for (int col = 0; col < dimension; ++col)
	{
	    R(row, col) = 0;
	}
    }
    for (int row = 0; row < dimension; ++row)
    {
	for (int col = row; col < dimension; ++col)
	{
	    R(row, col) = e[row]*xi[col];
	}
    }

    delete []e;
    delete []xi;
}

template <typename type>
type Matrix<type>::determinant(Matrix<type> &upperTriangle) const
{
    if (dimension != upperTriangle.dimension)
    {
	std::cout << "Matrix dimension is wrong. " << std::endl;
	exit(-1);
    }
    int count = 0;
    type **array;
    array = new type*[dimension];
    for (int i = 0; i < dimension; ++i)
    {
	array[i] = new type[dimension];
    }
    for (int i = 0; i < dimension; ++i)
	for (int j = 0; j < dimension; ++j)
	    array[i][j] = matrix[i][j];

    int col = 0;
    while (col < dimension)
    {
	int nonZero = col;
	while(array[nonZero][col] == 0)
	{
	    nonZero++;
	    if (nonZero == dimension)
	    {
		return 0;
	    }
	}

	if (nonZero != col)
	{
	    count++;
	    type temp;
	    for (int j = 0; j < dimension; ++j)
	    {
		temp = array[col][j];
		array[col][j] = array[nonZero][j];
		array[nonZero][j] = temp;
	    }
	}

	type factor = 1;
	for (int i = col+1; i < dimension; ++i)
	{
	    factor = -array[i][col]/array[col][col];
	    for (int j = 0; j < dimension; ++j)
	    {
		array[i][j] = array[i][j] + factor*array[col][j];
	    }
	}
	col++;
    }	
    type product = 1;
    for (int i = 0; i < dimension; ++i)
	product = product*array[i][i];

    upperTriangle = Matrix(array, dimension);

    for (int i = 0; i < dimension; ++i)
	delete []array[i];
    delete []array;

    if (count%2 == 0)
	return product;
    else
	return -product;
}

template <typename type>
type Matrix<type>::determinant() const
{
    Matrix<type> upperTriangle(dimension);
    type det = this->determinant(upperTriangle);
    return det;
}

template <typename type>
Matrix<type> Matrix<type>::inverse() const
{
    Matrix<type> upperTriangle(dimension);
    type det = this->determinant(upperTriangle);
    if (det == 0)
    {
	std::cout << "This matrix is not invertible. " << std::endl;
	std::exit(-1);
    }
    Matrix<type> Q(dimension);
    this->QRDecomposition(Q, upperTriangle);
    Vector<type> *auxMatrix;
    auxMatrix = new Vector<type>[dimension];
    for (int i = 0; i < dimension; ++i)
	auxMatrix[i].setLength(2*dimension);
    for (int i = 0; i < dimension; i++)
    {
	for (int j = 0; j < 2*dimension; j++)
	{
	    if (j < dimension)
	    {
		auxMatrix[i][j] = upperTriangle(i, j);
	    }
	    else
	    {
		if (j == i + dimension)
		    auxMatrix[i][j] = 1;
		else
		    auxMatrix[i][j] = 0;
	    }
	}
    }

    for (int row = dimension - 1; row >= 0; row--)
    {
	type factor = 1/auxMatrix[row][row];
	auxMatrix[row] = factor*auxMatrix[row];
	for (int i = row - 1; i >= 0; --i)
	{
	    factor = -auxMatrix[i][row]/auxMatrix[row][row];
	    auxMatrix[i] = auxMatrix[i] + factor*auxMatrix[row];
	}
    }

    for (int i = 0; i < dimension; ++i)
    {
	for (int j = 0; j < dimension; ++j)
	{
	    upperTriangle(i, j) = auxMatrix[i][j + dimension];
	}
    }
    delete []auxMatrix;
    Matrix<type> inverseMatrix(dimension);
    inverseMatrix = upperTriangle*Q.Hermite();
    return inverseMatrix;
}

template<typename type>
Matrix<type> Matrix<type>::Hermite() const
{
    Matrix<type> temp(dimension);
    for (int i = 0; i < dimension; ++i)
    {
	for (int j = 0; j < dimension; ++j)
	{
	    temp(j, i) = matrix[i][j].conjugate();
	}
    }
    return temp;
}

template<>
Matrix<double> Matrix<double>::Hermite() const
{
    Matrix<double> temp(dimension);
    for (int i = 0; i < dimension; ++i)
    {
	for (int j = 0; j < dimension; ++j)
	{
	    temp(j, i) = matrix[i][j];
	}
    }
    return temp;
}

template<typename type>
Matrix<type> Matrix<type>::QRIteration() const
{
    const int iterationMax = 100;
    Matrix<type> Q(dimension);
    Matrix<type> R(dimension);
    Matrix<type> A(dimension);
    int iteration = 0;
    this->QRDecomposition(Q, R);
    while(iteration < iterationMax)
    {
	A = R*Q;
	A.QRDecomposition(Q, R);
	iteration++;
    }
    return A;
}

template<typename type>
type Matrix<type>::cofactor(int i, int j) const
{
    this->checkBounds(i);
    this->checkBounds(j);
    Matrix<type> cofactorMatrix(dimension - 1);
    for (int row = 0; row < i; ++row)
    {
	for (int col = 0; col < j; ++col)
	{
	    cofactorMatrix(row, col) = matrix[row][col];
	}
    }
    for (int row = i + 1; row < dimension; ++row)
    {
	for (int col = j + 1; col < dimension; ++col)
	{
	    cofactorMatrix(row-1, col-1) = matrix[row][col];
	}
    }
    for (int row = 0; row < i; ++row)
    {
	for (int col = j + 1; col < dimension; ++col)
	{
	    cofactorMatrix(row, col - 1) = matrix[row][col];
	}
    }
    for (int row = i + 1; row < dimension; ++row)
	for (int col = 0; col < j; ++col)
	    cofactorMatrix(row-1, col) = matrix[row][col];
    if ((i + j)%2 == 0)
	return cofactorMatrix.determinant();
    else
	return -cofactorMatrix.determinant();
}

template<typename type>
Matrix<type> Matrix<type>::adjugate() const
{
    Matrix<type> adj(dimension);
    for (int i = 0; i < dimension; ++i)
    {
	for (int j = 0; j < dimension; ++j)
	{
	    adj(i, j) = this->cofactor(j, i);
	}
    }
    return adj;
}

template <typename T>
std::ostream & operator<< (std::ostream &os, const Matrix<T> &parameter)
{
    int dimension = parameter.dimension;
    for (int i = 0; i < dimension; ++i)
    {
	for (int j = 0; j < dimension; ++j)
	{
	    os << parameter.matrix[i][j] << "  ";
	}
	os << std::endl;
    }
    return os;
}
#endif
