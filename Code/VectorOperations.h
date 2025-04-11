#ifndef __VOper_h__
#define __VOper_h__

#include <cassert>
#include <iostream>
#include <numeric>
#include <tuple>
#include <vector>

#include "VectorStat.h"

using namespace std;

// ===============================================================================
// Sign of a scalar
// ===============================================================================

template <typename T> inline int sgn(T val) { return (T(0) < val) - (val < T(0)); }

// ===============================================================================
// Sum of two vectors: sum component by component
// ===============================================================================

template <typename T> inline vector<T> operator+(const vector<T>& a, const vector<T>& b) {
	assert(a.size() == b.size());
	vector<T> result(a.size());

	for (int i = 0; i < static_cast<int>(a.size()); i++)
		result[i] = a[i] + b[i];

	return result;
}

// ===============================================================================
// Difference of two vectors component by component
// ===============================================================================

template <typename T> inline vector<T> operator-(const vector<T>& a, const vector<T>& b) {
	assert(a.size() == b.size());
	vector<T> result(a.size());

	for (int i = 0; i < static_cast<int>(a.size()); i++)
		result[i] = a[i] - b[i];

	return result;
}

// ===============================================================================
// Negation of a vector by component
// ===============================================================================

template <typename T> inline vector<T> operator-(const vector<T>& a) {
	vector<T> result(a.size());

	for (int i = 0; i < static_cast<int>(a.size()); i++)
		result[i] = -a[i];

	return result;
}

// ===============================================================================
// Scalar product between two vectors
// ===============================================================================

template <typename T> inline T operator*(const vector<T>& a, const vector<T>& b) {
	assert(a.size() == b.size());
	return inner_product(a.begin(), a.end(), b.begin(), T(0));
}

// ===============================================================================
// Cross product between two vectors
// ===============================================================================
template <typename T> inline vector<T> CrossProduct(const vector<T>& a, const vector<T>& b) {
	assert(a.size() == b.size() && a.size() == 3);

	vector<T> result(3);
	result[0] = a[1] * b[2] - a[2] * b[1];
	result[1] = a[2] * b[0] - a[0] * b[2];
	result[2] = a[0] * b[1] - a[1] * b[0];

	return result;
}

// ===============================================================================
// Matrix size
// ===============================================================================
template <typename T> tuple<size_t, size_t> Size(const vector<vector<T>>& matrix) {
	size_t n_rows = matrix.size();
	size_t n_cols = (n_rows == 0) ? 0 : matrix[0].size();

	// Checks that rows are all the same length
	for (size_t i = 1; i < n_rows; i++) {
		if (matrix[i].size() != n_cols) {
			assert(matrix[i].size() == n_cols && "Rows must be the same size!");
		}
	}

	return make_tuple(n_rows, n_cols);
}

// ===============================================================================
// Matrix-vector multiplication
// ===============================================================================
template <typename T>
inline vector<T> operator*(const vector<vector<T>>& matrix, const vector<T>& vec) {
	size_t n_rows, n_cols;
	tie(n_rows, n_cols) = Size(matrix);

	if (n_rows <= 0 or n_cols != vec.size()) {
		cerr << "Error: tried applying matrix: " << endl;
		Print(matrix);
		cerr << "To vector : " << endl;
		Print(vec);
		assert(false);
	}

	vector<T> result(n_rows, T(0)); // Initialize result vector with zeros

	for (size_t i = 0; i < n_rows; ++i) {
		for (size_t j = 0; j < vec.size(); ++j) {
			result[i] += matrix[i][j] * vec[j];
		}
	}

	return result;
}

// ===============================================================================
// Norm of a vector
// ===============================================================================
template <typename T> inline T Norm(const vector<T>& a) { return sqrt(a * a); }

// ===============================================================================
// Product between a scalar and a vector
// ===============================================================================

template <typename T> inline vector<T> operator*(T c, const vector<T>& a) {
	vector<T> result(a.size());

	for (int i = 0; i < static_cast<int>(a.size()); i++)
		result[i] = c * a[i];

	return result;
}

// ===============================================================================
// Product between a vector and a scalar
// ===============================================================================

template <typename T> inline vector<T> operator*(const vector<T>& a, T c) {
	vector<T> result(a.size());
	for (int i = 0; i < static_cast<int>(a.size()); i++)
		result[i] = c * a[i];

	return result;
}

// ===============================================================================
// Division between a vector and a scalar
// ===============================================================================

template <typename T> inline vector<T> operator/(const vector<T>& a, T c) {
	vector<T> result(a.size());
	for (int i = 0; i < static_cast<int>(a.size()); i++)
		result[i] = a[i] / c;
	return result;
}

// ===============================================================================
// Add to a vector a, vector b, and the result is stored in a
// ===============================================================================

template <typename T> inline vector<T>& operator+=(vector<T>& a, const vector<T>& b) {
	assert(a.size() == b.size());

	for (int i = 0; i < static_cast<int>(a.size()); i++)
		a[i] += b[i];

	return a;
}

// ===============================================================================
// Subtract from a vector a, vector b, and the result is stored in a
// ===============================================================================

template <typename T> vector<T>& operator-=(vector<T>& a, const vector<T>& b) {
	assert(a.size() == b.size());

	for (int i = 0; i < static_cast<int>(a.size()); i++)
		a[i] -= b[i];

	return a;
}

// ===============================================================================
// Matrix-matrix multiplication
// ===============================================================================
template <typename T>
inline vector<vector<T>> operator*(const vector<vector<T>>& A, const vector<vector<T>>& B) {
	// Check if the matrices are compatible for multiplication
	size_t A_rows, A_cols, B_rows, B_cols;
	tie(A_rows, A_cols) = Size(A);
	tie(B_rows, B_cols) = Size(B);
	if (A_rows == 0 || B_rows == 0 || A_cols != B_rows) {
		cerr << "Error: Incompatible matrix dimensions for multiplication." << endl;
		cerr << "Matrix A: " << A_rows << "x" << (A_rows > 0 ? A_cols : 0) << endl;
		cerr << "Matrix B: " << B_rows << "x" << (B_rows > 0 ? B_cols : 0) << endl;
		assert(false); // Ensure proper dimensions
	}

	size_t rows = A_rows;	   // Number of rows in the result matrix
	size_t cols = B_cols;	   // Number of columns in the result matrix
	size_t inner_dim = B_rows; // Common dimension for multiplication

	// Initialize the result matrix with zeros
	vector<vector<T>> result(rows, vector<T>(cols, T(0)));

	// Perform matrix-matrix multiplication
	for (size_t i = 0; i < rows; ++i) {				 // Iterate over rows of A
		for (size_t j = 0; j < cols; ++j) {			 // Iterate over columns of B
			for (size_t k = 0; k < inner_dim; ++k) { // Summing over the shared dimension
				result[i][j] += A[i][k] * B[k][j];
			}
		}
	}

	return result;
}

// ===============================================================================
// Transpose of a matrix
// ===============================================================================

template <typename T> vector<vector<T>> Transpose(const vector<vector<T>>& matrix) {
	if (matrix.empty())
		return {};

	size_t n_rows, n_cols;
	tie(n_rows, n_cols) = Size(matrix);
	vector<vector<T>> transposed(n_cols, vector<T>(n_rows));

	for (size_t i = 0; i < n_rows; ++i) {
		for (size_t j = 0; j < n_cols; ++j) {
			transposed[j][i] = (j <= matrix[i].size()) ? matrix[i][j] : 0;
		}
	}

	return transposed;
}

// ===============================================================================
// Minor of a matrix (remover row "row" and column "col")
// ===============================================================================

template <typename T>
vector<vector<T>> Minor(const vector<vector<T>>& matrix, size_t row, size_t col) {
	size_t n_rows, n_cols;
	tie(n_rows, n_cols) = Size(matrix);
	assert(n_rows == n_cols);
	vector<vector<T>> minor_matrix(n_rows - 1, vector<T>(n_rows - 1));

	for (size_t i = 0, minor_i = 0; i < n_rows; ++i) {
		if (i == row)
			continue; // Skip the given row
		for (size_t j = 0, minor_j = 0; j < n_cols; ++j) {
			if (j == col)
				continue; // Skip the given column
			minor_matrix[minor_i][minor_j] = matrix[i][j];
			++minor_j;
		}
		++minor_i;
	}
	return minor_matrix;
}

// ===============================================================================
// Determinant of a matrix
// ===============================================================================
template <typename T> T Det(const vector<vector<T>>& matrix) {
	size_t n_rows, n_cols;
	tie(n_rows, n_cols) = Size(matrix);

	// Ensure the matrix is square
	assert(n_rows > 0 && n_rows == n_cols);

	// Base case: 1x1 matrix
	if (n_rows == 1) {
		return matrix[0][0];
	}

	// Base case: 2x2 matrix
	if (n_rows == 2) {
		return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
	}

	// Recursive case: Laplace expansion along the first row
	T determinant = T(0);
	for (size_t j = 0; j < n_rows; ++j) {
		// Calculate the minor matrix by removing row 0 and column j
		auto minor_matrix = Minor(matrix, 0, j);

		// Add the current term to the determinant using alternating signs
		T sign = (j % 2 == 0) ? 1 : -1;
		determinant += sign * matrix[0][j] * Det(minor_matrix);
	}

	return determinant;
}

// ===============================================================================
// Cofactor matrix
// ===============================================================================
template <typename T> vector<vector<T>> Cofactor(const vector<vector<T>>& matrix) {
	size_t n_rows, n_cols;
	tie(n_rows, n_cols) = Size(matrix);

	// Ensure the matrix is square
	assert(n_rows > 0 && n_rows == n_cols);

	// Initialize the cofactor matrix
	vector<vector<T>> cofactor_matrix(n_rows, vector<T>(n_cols));

	// Compute the cofactor for each element
	for (size_t i = 0; i < n_rows; ++i) {
		for (size_t j = 0; j < n_rows; ++j) {
			// Calculate the minor matrix
			auto minor_matrix = Minor(matrix, i, j);

			// Compute the cofactor with alternating sign
			T sign = ((i + j) % 2 == 0) ? 1 : -1;
			cofactor_matrix[i][j] = sign * Det(minor_matrix);
		}
	}

	return cofactor_matrix;
}

// ===============================================================================
// Inverse of a square matrix
// ===============================================================================
template <typename T> vector<vector<T>> Inverse(const vector<vector<T>>& matrix) {
	size_t n_rows, n_cols;
	tie(n_rows, n_cols) = Size(matrix);

	// Ensure the matrix is square
	assert(n_rows == n_cols && "Matrix must be square to compute the inverse!");

	// Calculate the determinant of the matrix
	T determinant = Det(matrix);
	if (determinant == T(0)) {
		cerr << "Error: Matrix is singular and cannot be inverted." << endl;
		assert(false);
	}

	// Compute the cofactor matrix
	vector<vector<T>> cofactor_matrix = Cofactor(matrix);

	// Compute the adjugate matrix (transpose of the cofactor matrix)
	vector<vector<T>> adjugate_matrix = Transpose(cofactor_matrix);

	// Compute the inverse matrix by dividing each element by the determinant
	vector<vector<T>> inverse_matrix(n_rows, vector<T>(n_cols));
	for (size_t i = 0; i < n_rows; ++i) {
		for (size_t j = 0; j < n_cols; ++j) {
			inverse_matrix[i][j] = adjugate_matrix[i][j] / determinant;
		}
	}

	return inverse_matrix;
}

#endif