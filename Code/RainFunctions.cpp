#include "RainFunctions.h"

#include "body.h"
#include "ray.h"
#include "tinyxml2.h"

using namespace std;
using namespace tinyxml2;

// Projects the Point on a plane perpendicular to v and passing through p
vector<double> Project(vector<double> Point, vector<double> p, vector<double> v) {
	vector<double> diff = p - Point;
	return (Point + v * (diff * v) / (v * v));
}

// Finds the vertex in the middle of the three seen faces of a parallelepiped defined by a point p
// and three sides
vector<double> FindMiddle(vector<double> p, vector<vector<double>> sides, vector<double> v) {
	for (size_t i = 0; i < sides.size(); i++) {
		if (sides[i] * v < 0)
			p += sides[i];
	}
	return p;
}

// Finds the hexagonal projection H of a parallelepiped defined by a vertex p and sides on a plane
// perpendicular to v and passing through P0
vector<vector<double>> FindHexProj(vector<double> p, vector<vector<double>> Side, vector<double> v,
								   vector<double> P0) {
	vector<vector<double>> H = {FindMiddle(p, Side, v)};
	vector<vector<double>> delta(
		3, vector<double>(3,
						  0.0)); // Used to calculate the position of the vertices to project
	for (int i = 0; i < 3; i++) {
		delta[i] = Side[i] * v < 0 ? ((double)-1) * Side[i] : Side[i];
	}
	H.push_back(H[0] + delta[0]);
	H.push_back(H[0] + delta[0] + delta[1]);
	H.push_back(H[0] + delta[1]);
	H.push_back(H[0] + delta[1] + delta[2]);
	H.push_back(H[0] + delta[2]);
	H.push_back(H[0] + delta[2] + delta[0]);
	for (int i = 0; i < 7; i++) {
		H[i] = Project(H[i], P0, v); // We project them
	}

	return H;
}

// Return the highest absolute value of the projections of the vertices of H on a line in direction
// u1 passing through H[0]
double MaxU(vector<vector<double>> H, vector<double> u) {
	double result = 0;
	for (int i = 1; i < 7; i++) {
		H[i] -= H[0];
		double proj = abs(H[i] * u / Norm(u));
		if (proj > result)
			result = proj;
	}
	// cout << result << endl;
	return result;
}

// Return wether the Point is inside the hexagon H using triangles and barycentric coordinates
bool PointIsInsideT(vector<double> Point, vector<vector<double>> H) {
	// Centers all other points on H[1]
	Point -= H[1];
	for (int i = 2; i < 7; i++) {
		H[i] -= H[1];
	}

	// Checks if Point is inside the triangle with vertices H[1], H[i], H[i+1]
	for (int i = 2; i < 6; i++) {
		int i_next = i + 1;
		double epsilon = 1e-10;
		double A = Norm(CrossProduct(H[i], H[i_next]));
		double alpha = Norm(CrossProduct(Point, H[i_next])) / A;
		double beta = Norm(CrossProduct(Point, H[i])) / A;
		double gamma = Norm(CrossProduct(Point - H[i], Point - H[i_next])) / A;
		if (0 <= alpha and alpha <= 1 and 0 <= beta and beta <= 1 and 0 <= gamma and gamma <= 1 and
			abs(alpha + beta + gamma - 1) < epsilon) {
			return true;
		}
	}
	return false;
}

// Periodic Boundary conditions for the index of H, keeps it between 1 and 6
int PBCH(int i) {
	while (i > 6)
		i -= 6;
	while (i < 1)
		i += 6;
	return i;
}

// Return the minimum distance between the point p and the segment line with extremes l1 and l2
double PointSegDist(vector<double> p, vector<double> l1, vector<double> l2) {
	// Changes frame of reference to l1 = 0
	p -= l1;
	l2 -= l1;
	double l2_norm2 = l2 * l2;
	// Makes sure not to divide by zero
	if (l2_norm2 == 0)
		return Norm(p);
	// Calculates projection of p on line passing through l1 (0) and l2, value of projection is
	// multiplied by l2_norm
	double proj = p * l2;
	// Return distance between p and the closest point belonging to the segment
	if (proj <= 0)
		// p closest to l1
		return Norm(p);
	if (proj >= l2_norm2)
		// p closest to l2
		return Norm(p - l2);

	// p closest to its projection on the segment, calculates distance with pythagoras
	return sqrt(max(0., p * p - proj * proj / l2_norm2));
}

// Return NxN identity matrix
vector<vector<double>> IdMat(unsigned int N) {
	vector<vector<double>> idmat(N, vector<double>(N));
	for (size_t i = 0; i < N; i++) {
		for (size_t j = 0; j < N; j++) {
			idmat[i][j] = (i == j) ? 1 : 0;
		}
	}
	return idmat;
}

// Return matrix representing the rotation around axis of theta radians
vector<vector<double>> RotMat(vector<double> axis, double theta) {
	if (Norm(axis) == 0 or axis.size() != 3)
		return IdMat(axis.size()); // Handle zero-length vector to avoid division by zero
	axis = axis / Norm(axis);
	double s = sin(theta);
	double c = cos(theta);
	double G = 1 - c;

	return {{axis[0] * axis[0] * G + c, axis[0] * axis[1] * G - axis[2] * s,
			 axis[0] * axis[2] * G + axis[1] * s},
			{axis[1] * axis[0] * G + axis[2] * s, axis[1] * axis[1] * G + c,
			 axis[1] * axis[2] * G - axis[0] * s},
			{axis[2] * axis[0] * G - axis[1] * s, axis[2] * axis[1] * G + axis[0] * s,
			 axis[2] * axis[2] * G + c}};
}

// Rotates a Point relative to the point rot0
void RotatePoint(vector<double>& Point, const vector<double>& rot0,
				 const vector<vector<double>>& rotmat) {
	if (rotmat == IdMat(3))
		return;
	Point -= rot0;
	Point = rotmat * Point;
	Point += rot0;
}

// Prints the shadow of a body at nstep different time steps in [tmin, tmax)
void PrintDynShadow(vector<double> box, Body& body, vector<double> relvel, double dx, double tmin,
					double tmax, unsigned int nstep, string outfile) {
	ProjSurface canvas(box, relvel, dx);
	double dt = nstep < 2 ? 0 : (tmax - tmin) / nstep;
	double t = tmin;

	for (unsigned int i = 0; i < nstep; i++) {
		body.Move(t);
		canvas.reset();
		canvas.BodyProj(body);
		string out = outfile + to_string(t) + ".dat";
		cout << "Printing to " << out << endl;
		canvas.PrintRaysFlat(out);
		t += dt;
	}
}

// Prints the state of a body at nstep different time steps in [tmin, tmax)
void PrintDynState(Body& body, double tmin, double tmax, unsigned int nstep, string outfile) {
	double dt = nstep < 2 ? 0 : (tmax - tmin) / nstep;
	double t = tmin;

	for (unsigned int i = 0; i < nstep; i++) {
		body.Move(t);
		string out = outfile + to_string(t) + ".dat";
		cout << "Printing to " << out << endl;
		body.PrintState(out);
		t += dt;
	}
}

// Transforms distance into a value in [0, 1]
double d_to_w(double delta_r, double dx) {
	delta_r /= dx;
	if (delta_r >= 1)
		return 0.;
	if (delta_r <= -1)
		return 1.;
	return (1 - sin(delta_r * M_PI / 2)) / 2;
}

// Estimates wetness
double Wetness(vector<double> box, Body& body, vector<double> rain_v, double vb, double dx) {
	vector<double> relvel = rain_v;
	relvel[0] -= vb;
	return Norm(relvel) * ProjSurface(box, relvel, dx).BodyProj(body) / vb;
}

// Estimates wetness of the dynamic body
double Wetness(vector<double> box, Body& body, vector<double> rain_v, double vb, double dx,
			   double tmin, double tmax, unsigned int nstep) {
	vector<double> relvel = rain_v;
	relvel[0] -= vb;
	return Norm(relvel) * ProjSurface(box, relvel, dx).BodyProj(body, tmin, tmax, nstep) / vb;
}

// Estimates wetness for N velocities of the dynamic body between vmin and vmax, and return
// a matrix with the velocities as the first colunmn and the respective wetness as the second column
vector<vector<double>> Simulate(vector<double> box, Body& body, vector<double> rain_v, double vmin,
								double vmax, unsigned int N, double dx, double tmin, double tmax,
								unsigned int nstep) {
	if (vmin > vmax or vmin < 0)
		cout << "Error: Vmin and Vmax have to be positive and Vmax > Vmin!" << endl;
	vector<double> body_v(N);
	vector<double> wetness(N);
	for (size_t i = 0; i < N; i++) {
		body_v[i] = (N == 1 ? vmin : vmin + (vmax - vmin) * (double)i / ((double)N - 1));
		wetness[i] = Wetness(box, body, rain_v, body_v[i], dx, tmin, tmax, nstep);
	}
	return Transpose(vector<vector<double>>{body_v, wetness});
}

// Estimates wetness for N velocities of the body between vmin and vmax, and return a matrix
// with the velocities as the first colunmn and the respective wetness as the second column
vector<vector<double>> Simulate(vector<double> box, Body& body, vector<double> rain_v, double vmin,
								double vmax, unsigned int N, double dx) {
	if (vmin > vmax or vmin < 0)
		cout << "Error: Vmin and Vmax have to be positive and Vmax > Vmin!" << endl;
	vector<double> body_v(N);
	vector<double> wetness(N);
	for (size_t i = 0; i < N; i++) {
		body_v[i] = (N == 1 ? vmin : vmin + (vmax - vmin) * (double)i / ((double)N - 1));
		wetness[i] = Wetness(box, body, rain_v, body_v[i], dx);
	}
	return Transpose(vector<vector<double>>{body_v, wetness});
}

// Estimates wetness for N values of dx between dxmin and dxmax, and return a matrix with dx
// as the first colunmn and the respective wetness as the second column
vector<vector<double>> SimErr(vector<double> box, Body& body, vector<double> rain_vel,
							  double body_vel, unsigned int N, double dxmin, double dxmax) {
	vector<double> d = {dxmax};
	vector<double> S = {Wetness(box, body, rain_vel, body_vel, dxmax)};
	double k = N == 0 ? 0 : pow(dxmin / dxmax, (double)1 / (N - 1));

	for (size_t i = 1; i < N; i++) {
		d.push_back(d[i - 1] * k);
		cout << "dx = " << d[i] << endl;
		S.push_back(Wetness(box, body, rain_vel, body_vel, d[i]));
	}

	return Transpose(vector<vector<double>>{d, S});
}

// Estimates wetness for N values of dx between dxmin and dxmax for dynamic body, and return
// a matrix with dx as the first colunmn and the respective wetness as the second column
vector<vector<double>> SimErr(vector<double> box, Body& body, vector<double> rain_vel,
							  double body_vel, unsigned int N, double dxmin, double dxmax,
							  double tmin, double tmax, unsigned int nstep) {
	vector<double> dx = {dxmax};
	vector<double> S = {Wetness(box, body, rain_vel, body_vel, dxmax, tmin, tmax, nstep)};
	double k = N == 0 ? 0 : pow(dxmin / dxmax, (double)1 / (N - 1));

	for (size_t i = 1; i < N; i++) {
		dx.push_back(dx[i - 1] * k);
		cout << "dx = " << dx[i] << endl;
		S.push_back(Wetness(box, body, rain_vel, body_vel, dx[i], tmin, tmax, nstep));
	}

	return Transpose(vector<vector<double>>{dx, S});
}

// Estimates wetness for N values of nstep between nstepmin and nstepmax, and return a
// matrix with nstep as the first colunmn and the respective wetness as the second column
vector<vector<double>> SimErrT(vector<double> box, Body& body, vector<double> rain_vel,
							   double body_vel, double dx, double tmin, double tmax, unsigned int N,
							   unsigned int nstepmin, unsigned int nstepmax) {
	double dtmin = 1. / nstepmax;
	double dtmax = 1. / nstepmin;
	vector<double> dt = {dtmax};
	vector<double> nstep = {(double)nstepmin};
	vector<double> S = {Wetness(box, body, rain_vel, body_vel, dx, tmin, tmax, nstepmin)};
	double k = N == 0 ? 0 : pow(dtmin / dtmax, (double)1 / (N - 1));

	for (size_t i = 1; i < N; i++) {
		dt.push_back(dt[i - 1] * k);
		double nstepnew = round(1. / dt[i]);
		if (nstepnew > nstep.back()) {
			cout << "Nstep = " << nstepnew << endl;
			nstep.push_back(nstepnew);
			S.push_back(Wetness(box, body, rain_vel, body_vel, dx, tmin, tmax, nstepnew));
		}
	}

	return Transpose(vector<vector<double>>{nstep, S});
}

// Estimates wetness for N_v velocities of the dynamic body between vmin and vmax, and
// N_nstep values of nstep between nstep_min and nstep_max and return a matrix with the nsteps
// velocities as the first colunmn, the velocities as the second column and the respective wetness
// as the third column
vector<vector<double>> SimulateNstep(vector<double> box, Body& body, vector<double> rain_v,
									 double vmin, double vmax, unsigned int N_v, double dx,
									 unsigned int nstep_min, unsigned int nstep_max,
									 unsigned int N_nstep) {
	double dtmin = 1. / nstep_max;
	double dtmax = 1. / nstep_min;
	double k = N_nstep < 2 ? 1 : pow(dtmin / dtmax, (double)1 / (N_nstep - 1));
	vector<double> nstep, vb, wetness;

	for (size_t i = 0; i < N_nstep; i++) {
		double dt_i = dtmax * pow(k, i);
		unsigned int nstep_i = round(1. / dt_i);
		for (size_t j = 0; j < N_v; j++) {
			double vb_j = (N_v < 2 ? vmin : vmin + (vmax - vmin) * (double)j / ((double)N_v - 1));
			nstep.push_back(nstep_i);
			vb.push_back(vb_j);
			wetness.push_back(Wetness(box, body, rain_v, vb_j, dx, 0., 1., nstep_i));
		}
	}

	return Transpose(vector<vector<double>>{nstep, vb, wetness});
}

// Fits points with a parabola y = k(x - x0)^2 + y0, using least squares minimization, and return a
// tuple containing (k, k_std, x0, x0_std, y0, y0_std)
tuple<double, double, double, double, double, double> ParabolicFit(vector<double> x_vals,
																   vector<double> y_vals) {
	size_t n = x_vals.size();
	assert(y_vals.size() == n && "Number of x and y values must be equal!");

	// Convert to long double
	vector<long double> x_long(n), y_long(n);
	transform(x_vals.begin(), x_vals.end(), x_long.begin(),
			  [](double val) { return static_cast<long double>(val); });
	transform(y_vals.begin(), y_vals.end(), y_long.begin(),
			  [](double val) { return static_cast<long double>(val); });

	// Build design matrix
	vector<vector<long double>> X(n, vector<long double>(3));
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < 3; j++) {
			X[i][j] = pow(x_long[i], j);
		}
	}

	// Transpose X
	vector<vector<long double>> X_t = Transpose(X);

	// Compute (X^T * X)^(-1)
	vector<vector<long double>> XtX_inv = Inverse(X_t * X);

	// Evaluate coefficients of y =  beta[0] + beta[1] * x + beta[2] * x^2
	vector<long double> beta = XtX_inv * X_t * y_long;

	// Ensures that beta[2] != 0
	if (beta[2] == 0)
		beta[2] = __DBL_EPSILON__;
	// cout << "a=" << beta[2] << " b=" << beta[1] << " c=" << beta[0] << endl;

	// Compute residuals
	vector<long double> residuals(n);
	for (size_t i = 0; i < n; ++i) {
		long double y_pred = beta[0] + beta[1] * x_long[i] + beta[2] * x_long[i] * x_long[i];
		residuals[i] = y_vals[i] - y_pred;
	}

	// Estimate variance of residuals (sigma^2)
	long double rss = 0.0;
	for (long double r : residuals)
		rss += r * r;
	long double sigma2 = rss / (n - 3); // degrees of freedom: n - number of parameters

	// Compute standard errors: sqrt(diag(sigma^2 * (X^T * X)^-1))
	// cout << " +-" << sqrt(sigma2 * XtX_inv[2][2]) << " +-" << sqrt(sigma2 * XtX_inv[1][1]) << " +-" << sqrt(sigma2 * XtX_inv[0][0]) << endl;

	// Evaluate parameters and errors
	long double k = beta[2];
	long double x0 = -beta[1] / (2 * beta[2]);
	long double y0 = beta[0] - beta[1] * beta[1] / (4 * beta[2]);

	// Evaluate partial derivatives
	long double dx0_db2 = beta[1] / (2 * beta[2] * beta[2]);
	long double dx0_db1 = -1 / (2 * beta[2]);
	long double dy0_db2 = beta[1] * beta[1] / (4 * beta[2] * beta[2]);
	long double dy0_db1 = -beta[1] / (2 * beta[2]);
	long double dy0_db0 = 1;

	// Evaluate erroros via error propagation
	long double k_std = sqrt(sigma2 * XtX_inv[2][2]);
	long double x0_std =
		sqrt(sigma2 * ((dx0_db2 * dx0_db2 * XtX_inv[2][2]) + dx0_db1 * dx0_db1 * XtX_inv[1][1] +
					   2 * dx0_db2 * dx0_db1 * XtX_inv[2][1]));
	long double y0_std = sqrt(
		sigma2 * (dy0_db2 * dy0_db2 * XtX_inv[2][2] + dy0_db1 * dy0_db1 * XtX_inv[1][1] +
				  dy0_db0 * dy0_db0 * XtX_inv[0][0] + 2 * dy0_db2 * dy0_db1 * XtX_inv[2][1] +
				  2 * dy0_db1 * dy0_db0 * XtX_inv[1][0] + 2 * dy0_db0 * dy0_db2 * XtX_inv[0][2]));

	// Return both parameters and errors
	return make_tuple(k, k_std, x0, x0_std, y0, y0_std);
}

// Finds minimums of wetness using Brent algorithm, calculates wetness for nfit values spaced
// dv around it, and return a tuple containing the optimal velocity, its error, the the minimum rain,
// its error, and a matrix containing the fit points, each row is a point, in the first colum are velocities
// and in the second the wetnesses
tuple<double, double, double, double, vector<vector<double>>>
MinFit(vector<double> box, Body& body, double vmax, double dx, unsigned int nstep, double vcross,
	   double vtail, int nfit, double dv) {
	vector<double> vb, wetness;
	Brent mins(dv);

	vector<double> rain_vel = {vtail, vcross, -1};
	auto wetfunc = [&box, &body, &rain_vel, dx, nstep](double x) {
		return Wetness(box, body, rain_vel, x, dx, 0., 1., nstep);
	};

	double k, k_std, vopt, vopt_std, Rmin, Rmin_std;
	if (mins.bracket(0, vmax, 3, wetfunc)) {
		mins.minimize(wetfunc);

		// Evaluate nfit points around minimum
		for (int j = -(nfit - 1) / 2; j <= nfit / 2; j++) {
			double vb_j = mins.xmin + j * dv;
			double wetness_j = (j == 0) ? mins.fmin : wetfunc(vb_j);
			vb.push_back(vb_j);
			wetness.push_back(wetness_j);
		}

		// Try quadratic fit, if resulting minimum isn't in the middle of evaluated points move
		// towards it and repeat until it is
		int max_iter = 100;
		int n_half = nfit / 2;
		bool moving_forward = false, moving_backwards = false;

		for (int iter = 0; iter < max_iter; iter++) {

			tie(k, k_std, vopt, vopt_std, Rmin, Rmin_std) = ParabolicFit(vb, wetness);

			// Avoid following minimum outside of range
			if (k > 0 && (vopt <= 0 || vopt >= vmax)) {
				break;
			}

			// Find number of values of vb that are lower (and higher) than vopt
			int n_lower = 0;
			while (n_lower < nfit && vb[n_lower] < vopt)
				n_lower++;
			int n_higher = nfit - n_lower;

			// Account for failed fit, if k < 0 move in  opposite direction
			if (k < 0)
				swap(n_lower, n_higher);

			// cout << "Iteration " << iter + 1 << ", n_lower = " << n_lower
			// 	 << ", n_higher = " << n_higher << endl;

			// Move towards minimum
			if (n_lower < n_half) {
				moving_backwards = true;
				double new_vb = vb[0] - dv;
				vb.insert(vb.begin(), new_vb);
				wetness.insert(wetness.begin(), wetfunc(new_vb));
				vb.pop_back();
				wetness.pop_back();
			} else if (n_higher < n_half) {
				moving_forward = true;
				double new_vb = vb.back() + dv;
				vb.push_back(new_vb);
				wetness.push_back(wetfunc(new_vb));
				vb.erase(vb.begin());
				wetness.erase(wetness.begin());
			} else
				break;

			// Avoid moving back and forth
			if (moving_backwards && moving_forward)
				break;
		}

		if (k <= 0 || vopt > vmax) {
			vopt = mins.cx;
			vopt_std = 0;
			Rmin = mins.fc;
			Rmin_std = 0;
		}

	} else { // If bracketing fails the bound minimum is in vb = vmax
		vb.push_back(mins.cx);
		wetness.push_back(mins.fc);
		vopt = mins.cx;
		vopt_std = 0;
		Rmin = mins.fc;
		Rmin_std = 0;
	}

	vector<vector<double>> fitPoints = Transpose(vector<vector<double>>{vb, wetness});

	return make_tuple(vopt, vopt_std, Rmin, Rmin_std, fitPoints);
}

// Finds minimums of wetness for a fixed vcross and [vtail_min, vtail_max] using Brent
// algorithm, and calculates wetness for nfit values spaced dv around it, return all these values
vector<vector<double>> FindMinFit(vector<double> box, Body& body, double vmax, double dx,
								  unsigned int nstep, double vcross, double vtail_min,
								  double vtail_max, unsigned int n_tail, int nfit, double dv) {
	vector<vector<double>> res;

	for (size_t i = 0; i < n_tail; i++) {
		double vtail_i =
			n_tail > 1 ? vtail_min + i * (vtail_max - vtail_min) / (n_tail - 1) : vtail_min;

		// Calculate fit points
		vector<vector<double>> res_i;
		tie(ignore, ignore, ignore, ignore, res_i) =
			MinFit(box, body, vmax, dx, nstep, vcross, vtail_i, nfit, dv);

		// Add vtail_i to each row
		for (size_t i = 0; i < res_i.size(); ++i) {
			res_i[i].insert(res_i[i].begin(), vtail_i);
		}
		// Add resulting points to output
		res.insert(res.end(), res_i.begin(), res_i.end());

		cout << "Step " << i + 1 << "/" << n_tail << " completed" << endl;
	}

	return res;
}

// Finds minimums of wetness for a fixed vcross and vtail_min using Brent algorithm with
// nstep in [nstep_min, nstep_max], and calculates wetness for nfit values spaced dv around it,
// return all these values
vector<vector<double>> FindMinFitNstep(vector<double> box, Body& body, double vmax, double dx,
									   unsigned int nstep_min, unsigned int nstep_max,
									   unsigned int N_nstep, double vcross, double vtail, int nfit,
									   double dv) {
	double dtmin = 1. / nstep_max;
	double dtmax = 1. / nstep_min;
	double k = N_nstep < 2 ? 1 : pow(dtmin / dtmax, (double)1 / (N_nstep - 1));

	vector<vector<double>> res;

	for (size_t i = 0; i < N_nstep; i++) {
		double dt_i = dtmax * pow(k, i);
		unsigned int nstep_i = round(1. / dt_i);

		// Calculate fit points
		vector<vector<double>> res_i;
		tie(ignore, ignore, ignore, ignore, res_i) =
			MinFit(box, body, vmax, dx, nstep_i, vcross, vtail, nfit, dv);

		// Add vtail_i to each row
		for (size_t i = 0; i < res_i.size(); ++i) {
			res_i[i].insert(res_i[i].begin(), nstep_i);
		}

		// Add resulting points to output
		res.insert(res.end(), res_i.begin(), res_i.end());

		cout << "Step " << i + 1 << "/" << N_nstep << " completed" << endl;
	}

	return res;
}

// Finds minimums of wetness for a fixed vcross and [vtail_min, vtail_max]x[vcross_min,
// vcross_max] with brent, calculates wetness for nfit values around it, return all these values
vector<vector<double>> OptMapFit(vector<double> box, Body& body, double vmax, double dx,
								 unsigned int nstep, unsigned int nfit, double dv, double vtail_min,
								 double vtail_max, unsigned int n_tail, double vcross_min,
								 double vcross_max, unsigned int n_cross) {
	vector<vector<double>> res;

	for (size_t i = 0; i < n_cross; i++) {
		double vcross_i =
			n_cross > 1 ? vcross_min + i * (vcross_max - vcross_min) / (n_cross - 1) : vcross_min;

		// Calculate fit points
		vector<vector<double>> res_i = FindMinFit(box, body, vmax, dx, nstep, vcross_i, vtail_min,
												  vtail_max, n_tail, nfit, dv);

		// Add vcross_i to each row
		for (size_t i = 0; i < res_i.size(); ++i) {
			res_i[i].insert(res_i[i].begin(), vcross_i);
		}

		// Add resulting points to output
		res.insert(res.end(), res_i.begin(), res_i.end());

		cout << "Super step " << i + 1 << "/" << n_cross << " completed" << endl;
	}

	return res;
}

// Write header file for results of minimization with varying vtail and vcross
void WriteHeadRes(ostream& out, string bodyName, double vmax, double dx, int nstep, int nfit,
				  double dv) {
	out << "#######################################################################################"
		   "#####################\n";
	out << "# Minimums of wetness found for the following parameters:\n";
	out << "# Body = " + bodyName + ", vmax = " << vmax << " vfall\n";
	out << "# dx = " << dx << " m, nstep = " << nstep << "\n";
	out << "# nfit = " << nfit << ", dv = " << dv << " vfall\n";
	out << "# Columns:\n";
	out << "# vcross (vfall)   vtail (vfall)   vopt (vfall)   vopt_std (vfall)   Rmin (m^2)   Rmin_std "
		   "(m^2)\n";
	out << "#######################################################################################"
		   "#####################"
		<< endl;
}

// Write header file for fit points of minimization with varying vtail and vcross
void WriteHeadFit(ostream& out, string bodyName, double vmax, double dx, int nstep, int nfit,
				  double dv) {
	out << "#######################################################################################"
		   "#####################\n";
	out << "# Fit points used to estimate minimums of wetness found for the following "
		   "parameters:\n";
	out << "# Body = " + bodyName + ", vmax = " << vmax << " vfall\n";
	out << "# dx = " << dx << " m, nstep = " << nstep << "\n";
	out << "# nfit = " << nfit << ", dv = " << dv << " vfall\n";
	out << "# Columns:\n";
	out << "# vcross (vfall)   vtail (vfall)   vb (vfall)   Rb (m^2)\n";
	out << "#######################################################################################"
		   "#####################"
		<< endl;
}

// Perform a test calcultation and check that results are right
bool AllGood() {
	ManyBody testBody;

	vector<double> centS = {0., 0., 0.};
	double radS = 0.3;
	vector<double> rotcentS(3), axisS(3);
	vector<vector<double>> transS = {{0.2, 0.1, 0}, {0., 0.2, 0.1}};
	vector<double> transPhiS = {M_PI / 2, 0.3};
	vector<double> delta0S = {0., -0.2, -0.1};
	testBody.AddBody(
		Sphere(centS, radS, "testSphere", rotcentS, axisS, {}, {}, transS, transPhiS, 0, delta0S));

	vector<double> centP = {0., 0.1, 0.};
	vector<vector<double>> sidesP = {{0.1, 0, 0}, {0, 0.4, 0}, {0, 0, 0.2}};
	vector<double> rotcentP = {0., 0.45, 0.};
	vector<double> axisP = {1., 0., 0.};
	double w0P = -M_PI / 2;
	vector<double> wP = {60., 90.};
	vector<double> wPhiP = {0.4, 2};
	vector<double> delta0P = {0., -0.5, -0.1};
	vector<vector<double>> transP = {{0., 0.1, 0.05}, {0., 0.5, 0.1}};
	vector<double> transPhiP = {1, M_PI / 2};
	testBody.AddBody(Parallelepiped(centP, sidesP, "testPippo", rotcentP, axisP, wP, wPhiP, transP,
									transPhiP, w0P, delta0P));

	vector<double> l1 = {0., 0.1, 0.1};
	vector<double> l2 = {0., 0.0, -0.2};
	double radC = 0.15;
	vector<double> rotcentC = l1;
	vector<double> axisC = {0., 0., 1};
	vector<double> wC = {15., 20.};
	vector<double> wPhiC = {4, 3 * M_PI / 2};
	double w0C = -20;
	vector<vector<double>> transC;
	testBody.AddBody(
		Capsule(l1, l2, radC, "testCapsule", rotcentC, axisC, wC, wPhiC, transC, {}, w0C));
	testBody.Attach("testCapsule", "testPippo");

	int tmin = 0, tmax = 1, nstep = 10;
	double dx = 0.01;
	vector<double> testBox = testBody.GetBox(tmin, tmax, nstep, 2 * dx);

	vector<double> rainVel = {0.1, 0.2, -1};
	double vb = 0.5;

	double expected = 0.904234010057522;
	double evaluated = Wetness(testBox, testBody, rainVel, vb, dx, tmin, tmax, nstep);

	double epsilon = expected * 1e-10;
	if (abs(expected - evaluated) < epsilon)
		return true;
	cout << setprecision(15) << "WARNING: SOMETHING IS NOT WORKING!\n"
		 << "Expected test value:  " << expected << "\n"
		 << "Evaluated test value: " << evaluated << endl;
	return false;
}

// Parse a vector of doubles from a string
vector<double> parseDoubles(const char* text) {
	vector<double> result;
	if (text == nullptr)
		return result;
	istringstream iss(text);
	double val;
	while (iss >> val)
		result.push_back(val);
	return result;
}

// Try to get a vector of three doubles from an XMLElement
vector<double> SafeGet3Vec(XMLElement* parent, const char* tag, string bodyName) {
	XMLElement* child = parent->FirstChildElement(tag);
	if (!child || !child->GetText())
		throw std::logic_error(string("Missing or empty <" + string(tag) + "> in " + bodyName));
	vector<double> res = parseDoubles(child->GetText());
	if (res.size() != 3)
		throw logic_error(
			string("Invalid <" + string(tag) + "> value in " + bodyName + ", must be 3 numbers"));
	return res;
}

// Try to get text from an XMLElement, throw error if it is null
string SafeGetText(XMLElement* parent, const char* tag) {
	XMLElement* child = parent->FirstChildElement(tag);
	if (!child || !child->GetText())
		throw std::logic_error(string("Missing or empty <") + tag + ">");
	return string(child->GetText());
}