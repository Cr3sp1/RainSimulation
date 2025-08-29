#ifndef __VStats_h__
#define __VStats_h__

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>

using namespace std;

template <typename T> double Avg(const vector<T>& v) {
	double sum = accumulate(v.begin(), v.end(), 0.0);
	return sum / v.size();
}

template <typename T> double Var(const vector<T>& v) {
	if (v.size() == 1)
		return 0;
	double va = 0;
	for (int i = 0; i < v.size(); i++) {
		T avg = Avg<T>(v);
		va = va + pow(avg - v[i], 2);
	}
	return va / (v.size() - 1);
}

template <typename T> double StDev(const vector<T>& v) { return sqrt(Var(v)); }

template <typename T> double Med(vector<T> v) {
	double med;
	sort(v.begin(), v.end());

	if (v.size() % 2 == 0) {
		med = (v[v.size() / 2 - 1] + v[v.size() / 2]) / 2;
	} else {
		med = v[v.size() / 2];
	}

	return med;
}

template <typename T> double Covariance(const vector<T>& x, const vector<T>& y) {
	double mediax = Avg(x);
	double mediay = Avg(y);
	double result{0.};

	if (x.size() != y.size()) {
		cout << "ERROR Covariance, vectors of different size!" << endl;
		return -1;
	}

	if (x.size() == 0) {
		return result;
	}

	for (int i = 0; i < x.size(); i++) {
		result += (x[i] - mediax) * (y[i] - mediay);
	}

	return result / (StDev(x) * StDev(y));
}

template <typename T> vector<T> Read(unsigned int size, char* filename) {
	ifstream fin(filename);
	if (!fin) {
		cout << "File not found" << endl;
		exit(-1);
	}
	vector<T> v;
	T appoggio;
	for (int k = 0; k < size; k++) {
		fin >> appoggio;
		v.push_back(appoggio);
		if (fin.eof()) {
			cout << "End of file reached before: " << size << endl;
			break;
		}
	}
	fin.close();
	return v;
}

template <typename T> void Print(ostream& fout, const vector<T>& v) {
	for (const auto& element : v) {
		fout << element << "\t";
	}
	fout << endl;
}

template <typename T> void Print(const vector<T>& v) {
	for (const auto& element : v) {
		cout << element << "\t";
	}
	cout << endl;
}

template <typename T> void Print(ostream& fout, const vector<T>& v, unsigned int precision) {
	for (const auto& element : v) {
		fout << setprecision(precision) << element << "\t";
	}
	fout << endl;
}

template <typename T> void Print(const vector<T>& v, unsigned int precision) {
	for (const auto& element : v) {
		cout << setprecision(precision) << element << "\t";
	}
	cout << endl;
}

template <typename T>
void Print(ostream& fout, const vector<vector<T>>& matrix, unsigned int precision) {
	for (const auto& row : matrix) {
		for (const auto& element : row) {
			fout << setprecision(precision) << element << "\t";
		}
		fout << endl;
	}
}

template <typename T>
void Print(string outfile, const vector<vector<T>>& matrix, unsigned int precision) {
	ofstream fout(outfile);
	Print(fout, matrix, precision);
	fout.close();
}

template <typename T> void Print(ostream& fout, const vector<vector<T>>& matrix) {
	for (const auto& row : matrix) {
		for (const auto& element : row) {
			fout << element << "\t";
		}
		fout << endl;
	}
}

template <typename T> void Print(string outfile, const vector<vector<T>>& matrix) {
	ofstream fout(outfile);
	Print(fout, matrix);
	fout.close();
}

template <typename T> void Print(const vector<vector<T>>& matrix) {
	for (const auto& row : matrix) {
		for (const auto& element : row) {
			cout << element << "\t";
		}
		cout << endl;
	}
}

template <typename T> void Print(const vector<vector<T>>& matrix, unsigned int precision) {
	for (const auto& row : matrix) {
		for (const auto& element : row) {
			cout << setprecision(precision) << element << "\t";
		}
		cout << endl;
	}
}

#endif