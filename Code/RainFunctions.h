#ifndef __RainFunctions_h__
#define __RainFunctions_h__

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <tuple>
#include <vector>

#include "VectorOperations.h"
#include "mins.h"

using namespace std;

// Forward declaration
class Body;
class Ray;
class ProjSurface;

// Projects the Point on a plane perpendicular to v and passing through p
vector<double> Project(vector<double> Point, vector<double> p, vector<double> v);

// Finds the vertex in the middle of the three seen faces of a parallelepiped defined by a point p
// and three sides
vector<double> FindMiddle(vector<double> p, vector<vector<double>> sides, vector<double> v);

// Finds the hexagonal projection H of a parallelepiped on a plane perpendicular to v and passing
// through P
vector<vector<double>> FindHexProj(vector<double> p, vector<vector<double>> Side, vector<double> v,
								   vector<double> P);

// Returns the highest absolute value of the projections of the vertices of H on a line in direction
// u1 passing through H[0]
double MaxU(vector<vector<double>> H, vector<double> u);

// Checks wether a point is inside the hexagon using triangles and baycentric
// coordinates
bool PointIsInsideT(vector<double> Point, vector<vector<double>> H);

// Periodic Boundary conditions for the index of H
int PBCH(int i);

// Checks rays generation
void RayGenCheck(string outfile, vector<double> box, vector<double> rel_vel);

// Estimates wetness
double Wetness(vector<double> box, Body& body, vector<double> rain_v, double vb, double dx);

// Estimates wetness of the dynamic body
double Wetness(vector<double> box, Body& body, vector<double> rain_v, double vb, double dx,
			   double tmin, double tmax, unsigned int nstep);

// Estimates wetness for N velocities of the body between vmin and vmax, and returns a matrix with
// the velocities as the first colunmn and the respective wetness as the second column
vector<vector<double>> Simulate(vector<double> box, Body& body, vector<double> rain_v, double vmin,
								double vmax, unsigned int N, double dx);

// Estimates wetness for N velocities of the dynamic body between vmin and vmax, and returns a
// matrix with the velocities as the first colunmn and the respective wetness as the second column
vector<vector<double>> Simulate(vector<double> box, Body& body, vector<double> rain_v, double vmin,
								double vmax, unsigned int N, double dx, double tmin, double tmax,
								unsigned int nstep);

// Estimates wetness for N values of dx between dxmin and dxmax, and returns a matrix with dx as the
// first colunmn and the respective wetness as the second column
vector<vector<double>> SimErr(vector<double> box, Body& body, vector<double> rain_vel,
							  double body_vel, unsigned int N, double dxmin, double dxmax);

// Estimates wetness for N values of dx between dxmin and dxmax for dynamic body, and returns a
// matrix with dx as the first colunmn and the respective wetness as the second column
vector<vector<double>> SimErr(vector<double> box, Body& body, vector<double> rain_vel,
							  double body_vel, unsigned int N, double dxmin, double dxmax,
							  double tmin, double tmax, unsigned int nstep);

// Estimates wetness for N values of nstep between nstepmin and nstepmax, and returns a matrix with
// nstep as the first colunmn and the respective wetness as the second column
vector<vector<double>> SimErrT(vector<double> box, Body& body, vector<double> rain_vel,
							   double body_vel, double dx, double tmin, double tmax, unsigned int N,
							   unsigned int nstepmin, unsigned int nstepmax);

// Returns the minimum distance between the point p and the segment line with extremes l1 and l2
double PointSegDist(vector<double> p, vector<double> l1, vector<double> l2);

// Returns NxN identity matrix
vector<vector<double>> IdMat(unsigned int N);

// Returns the rotation matrix
vector<vector<double>> RotMat(vector<double> axis, double theta);

// Rotates a Point relative to the point Rot0
void Rotate(vector<double>& Point, const vector<double>& Rot0,
			const vector<vector<double>>& Rotmat);

// Prints the shadow of a body at nstep different time steps in [tmin, tmax)
void PrintDynShadow(vector<double> box, Body& body, vector<double> relvel, double dx, double tmin,
					double tmax, unsigned int nstep, string outfile);

// Prints the smooth shadow of a body at nstep different time steps in [tmin, tmax)
void PrintDynShadowSmooth(vector<double> box, Body& body, vector<double> relvel, double dx,
						  double tmin, double tmax, unsigned int nstep, string outfile);

// Prints the state of a body at nstep different time steps in [tmin, tmax)
void PrintDynState(Body& body, double tmin, double tmax, unsigned int nstep, string outfile);

// Looks for the minimum of wetness between vmin and vmax with grid search with N steps and returns
// its value if it finds it, else returns -1.
double FindMin(vector<double> box, Body& body, vector<double> rain_v, double vmin, double vmax,
			   unsigned int N, double dx, unsigned int nstep);

// Finds minimums of wetness for a fixed vcross and [vtail_min, vtail_max] with grid search, and
// calculates wetness for n_fit values around it, returns all these values
vector<vector<double>> FindMinFit(vector<double> box, Body& body, double vmin, double vmax,
								  unsigned int N, double dx, unsigned int nstep, unsigned int n_fit,
								  double vcross, double vtail_min, double vtail_max,
								  unsigned int n_tail);

// Finds minimums of wetness for a fixed vcross and [vtail_min, vtail_max]x[vcross_min, vcross_max]
// with grid search, and calculates wetness for n_fit values around it, returns all these values
vector<vector<double>> OptMapFit(vector<double> box, Body& body, double vmin, double vmax,
								 unsigned int N, double dx, unsigned int nstep, unsigned int n_fit,
								 double vtail_min, double vtail_max, unsigned int n_tail,
								 double vcross_min, double vcross_max, unsigned int n_cross);

// Transforms distance into a value in [0, 1]
double smooth_w(double delta_r, double dx);

// Estimates smooth wetness
double WetnessSmooth(vector<double> box, Body& body, vector<double> rain_v, double vb, double dx);

// Estimates smooth wetness of the dynamic body
double WetnessSmooth(vector<double> box, Body& body, vector<double> rain_v, double vb, double dx,
					 double tmin, double tmax, unsigned int nstep);

// Estimates smooth wetness for N velocities of the body between vmin and vmax, and returns a matrix
// with the velocities as the first colunmn and the respective wetness as the second column
vector<vector<double>> SimulateSmooth(vector<double> box, Body& body, vector<double> rain_v,
									  double vmin, double vmax, unsigned int N, double dx);

// Estimates smooth wetness for N velocities of the dynamic body between vmin and vmax, and returns
// a matrix with the velocities as the first colunmn and the respective wetness as the second column
vector<vector<double>> SimulateSmooth(vector<double> box, Body& body, vector<double> rain_v,
									  double vmin, double vmax, unsigned int N, double dx,
									  double tmin, double tmax, unsigned int nstep);

// Estimates smooth wetness for N values of dx between dxmin and dxmax, and returns a matrix with dx
// as the first colunmn and the respective wetness as the second column
vector<vector<double>> SimSmoothErr(vector<double> box, Body& body, vector<double> rain_vel,
									double body_vel, unsigned int N, double dxmin, double dxmax);

// Estimates smooth wetness for N values of dx between dxmin and dxmax for dynamic body, and returns
// a matrix with dx as the first colunmn and the respective wetness as the second column
vector<vector<double>> SimSmoothErr(vector<double> box, Body& body, vector<double> rain_vel,
									double body_vel, unsigned int N, double dxmin, double dxmax,
									double tmin, double tmax, unsigned int nstep);

// Estimates smooth wetness for N values of nstep between nstepmin and nstepmax, and returns a
// matrix with nstep as the first colunmn and the respective wetness as the second column
vector<vector<double>> SimSmoothErrT(vector<double> box, Body& body, vector<double> rain_vel,
									 double body_vel, double dx, double tmin, double tmax,
									 unsigned int N, unsigned int nstepmin, unsigned int nstepmax);

// Estimates smooth wetness for N_v velocities of the dynamic body between vmin and vmax, and
// nN_nstep values of nstep between nstep_min and nstep_max and returns a matrix with the nsteps
// velocities as the first colunmn, the velocities as the second column and the respective wetness
// as the third column
vector<vector<double>> SimulateNstepSmooth(vector<double> box, Body& body, vector<double> rain_v,
										   double vmin, double vmax, unsigned int N_v, double dx,
										   unsigned int nstep_min, unsigned int nstep_max,
										   unsigned int N_nstep);

// Fits points with a parabola y = ax^2 + bx + c, using least squares minimization, and returns a
// tuple containing (a, b, c)
tuple<double, double, double> QuadraticRegression(vector<double> x_vals, vector<double> y_vals);

// Finds minimums of smooth wetness using Brent algorithm,calculates wetness for n_fit values spaced
// dv around it, and returns two vectors, the first containing the values of vb and the second the
// respective wetness
tuple<vector<double>, vector<double>> MinFitSmooth(vector<double> box, Body& body, double vmin,
												   double vmax, double dx, unsigned int nstep,
												   double vcross, double vtail, int n_fit,
												   double dv);

// Finds minimums of smooth wetness for a fixed vcross and [vtail_min,
// vtail_max] using Brent algorithm, and calculates wetness for n_fit values
// spaced dv around it, returns all these values
vector<vector<double>> FindMinFitSmooth(vector<double> box, Body& body, double vmin, double vmax,
										double dx, unsigned int nstep, double vcross,
										double vtail_min, double vtail_max, unsigned int n_tail,
										int n_fit, double dv);

// Finds minimums of smooth wetness for a fixed vcross and vtail_min using Brent algorithm with
// nstep in [nstep_min, nstep_max], and calculates wetness for n_fit values spaced dv around it,
// returns all these values
vector<vector<double>> FindMinFitSmoothNstep(vector<double> box, Body& body, double vmin,
											 double vmax, double dx, unsigned int nstep_min,
											 unsigned int nstep_max, unsigned int N_nstep,
											 double vcross, double vtail, int n_fit, double dv);

// Finds minimums of smooth wetness for a fixed vcross and [vtail_min, vtail_max]x[vcross_min,
// vcross_max] with brent, calculates wetness for n_fit values around it, returns all these values
vector<vector<double>> OptMapFitSmooth(vector<double> box, Body& body, double vmin, double vmax,
									   double dx, unsigned int nstep, unsigned int n_fit, double dv,
									   double vtail_min, double vtail_max, unsigned int n_tail,
									   double vcross_min, double vcross_max, unsigned int n_cross);

#endif