#ifndef __mins_h__
#define __mins_h__

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace std;

struct Bracketmethod {
  public:
	// Member variables
	double ax, bx, cx, fa, fb, fc;
	const double tol;

	// Constructor with default tolerance
	Bracketmethod(const double toll = 3.0e-8) : tol(toll) {}

	// Bracket method, return true if it finds b so that f(b)<f(a) and f(b)<f(c)
	template <class T>
	bool bracket(const double a, const double c, const unsigned int n_tries, T& func) {
		ax = a + tol;
		cx = c;
		fa = func(ax);
		fc = func(cx);

		for (unsigned int i = 1; i <= n_tries; i++) {
			bx = c - i * tol;
			fb = func(bx);
			if (fb < fa and fb < fc) {
				return true;
			}
		}
		return false;
	}
};

struct Brent : public Bracketmethod {
  public:
	// Member variables for Brent's method
	double xmin, fmin;

	// Constructor with default tolerance
	Brent(const double toll = 3.0e-8) : Bracketmethod(toll) {}

	// Brent's minimization method, return all the points evaluated in order of evaluation
	template <class T> vector<vector<double>> minimize(T& func) {
		vector<double> xvals = {ax, cx, bx};
		vector<double> fvals = {fa, fc, fb};

		// cout << "a = " << ax << ", b = " << bx << ", c = " << cx << endl;

		const int ITMAX = 100;			// Maximum number of iterations
		const double CGOLD = 0.3819660; // Golden ratio
		const double ZEPS =
			numeric_limits<double>::epsilon() * 1.0e-3; // Small number to prevent division by zero

		// Initial setup
		double a, b, d = 0.0, etemp, fu, fv, fw, fx;
		double p, q, r, tol1, tol2, u, v, w, x, xm;
		double e = 0.0; // Distance moved on the previous step

		// Ensure a and b are in ascending order based on ax and cx
		a = (ax < cx) ? ax : cx;
		b = (ax > cx) ? ax : cx;

		// Initial point setup
		x = w = v = bx;
		fw = fv = fx = func(x);

		// Main iteration loop
		for (int iter = 0; iter < ITMAX; ++iter) {
			xm = 0.5 * (a + b);
			tol1 = tol * abs(x) + ZEPS;
			tol2 = 2.0 * tol1;

			// Check convergence
			if (abs(x - xm) <= (tol2 - 0.5 * (b - a))) {
				fmin = fx;
				xmin = x;
				// cout << "Number of Brent iterations: " << iter << endl;
				return {xvals, fvals};
			}

			// Construct a parabolic fit
			if (abs(e) > tol1) {
				r = (x - w) * (fx - fv);
				q = (x - v) * (fx - fw);
				p = (x - v) * q - (x - w) * r;
				q = 2.0 * (q - r);
				if (q > 0.0)
					p = -p;
				q = abs(q);
				etemp = e;
				e = d;

				// Conditions for acceptability of parabolic step
				if (abs(p) >= abs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x)) {
					// Take the golden section step
					d = CGOLD * (e = (x >= xm ? a - x : b - x));
				} else {
					// Take the parabolic step
					d = p / q;
					u = x + d;
					if ((u - a) < tol2 || (b - u) < tol2)
						d = SIGN(tol1, xm - x);
				}
			} else {
				// Golden section step
				d = CGOLD * (e = (x >= xm ? a - x : b - x));
			}

			u = (abs(d) >= tol1) ? x + d : x + SIGN(tol1, d);
			fu = func(u);
			xvals.push_back(u);
			fvals.push_back(fu);

			// Evaluation and updating of bounds
			if (fu <= fx) {
				if (u >= x)
					a = x;
				else
					b = x;
				shft3(v, w, x, u);
				shft3(fv, fw, fx, fu);
			} else {
				if (u < x)
					a = u;
				else
					b = u;
				if (fu <= fw || w == x) {
					v = w;
					w = u;
					fv = fw;
					fw = fu;
				} else if (fu <= fv || v == x || v == w) {
					v = u;
					fv = fu;
				}
			}
		}

		// If the loop completes without convergence
		throw runtime_error("Too many iterations in Brent's method");
	}

  private:
	// Utility functions
	inline void shft3(double& a, double& b, double& c, const double d) {
		a = b;
		b = c;
		c = d;
	}

	// Sign function, returning `a` with the sign of `b`
	inline double SIGN(const double a, const double b) { return (b >= 0.0) ? abs(a) : -abs(a); }
};

#endif // __mins_h__