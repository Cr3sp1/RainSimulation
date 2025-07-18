#ifndef __RainFunctions_h__
#define __RainFunctions_h__

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <ostream>
#include <tuple>
#include <vector>

#include "VectorOperations.h"
#include "mins.h"
#include "tinyxml2.h"

using namespace std;
using namespace tinyxml2;

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

// Return the highest absolute value of the projections of the vertices of H on a line in direction
// u1 passing through H[0]
double MaxU(vector<vector<double>> H, vector<double> u);

// Checks wether a point is inside the hexagon using triangles and baycentric
// coordinates
bool PointIsInsideT(vector<double> Point, vector<vector<double>> H);

// Periodic Boundary conditions for the index of H
int PBCH(int i);

// Checks rays generation
void RayGenCheck(string outfile, vector<double> box, vector<double> rel_vel);

// Return the minimum distance between the point p and the segment line with extremes l1 and l2
double PointSegDist(vector<double> p, vector<double> l1, vector<double> l2);

// Return NxN identity matrix
vector<vector<double>> IdMat(unsigned int N);

// Return the rotation matrix
vector<vector<double>> RotMat(vector<double> axis, double theta);

// Rotates a Point relative to the point rot0
void RotatePoint(vector<double>& Point, const vector<double>& rot0,
				 const vector<vector<double>>& rotmat);

// Prints the shadow of a body at nstep different time steps in [tmin, tmax)
void PrintDynShadow(vector<double> box, Body& body, vector<double> relvel, double dx, double tmin,
					double tmax, unsigned int nstep, string outfile);

// Prints the state of a body at nstep different time steps in [tmin, tmax)
void PrintDynState(Body& body, double tmin, double tmax, unsigned int nstep, string outfile);

// Transforms distance into a value in [0, 1]
double d_to_w(double delta_r, double dx);

// Estimates wetness
double Wetness(vector<double> box, Body& body, vector<double> rain_v, double vb, double dx);

// Estimates wetness of the dynamic body
double Wetness(vector<double> box, Body& body, vector<double> rain_v, double vb, double dx,
			   double tmin, double tmax, unsigned int nstep);

// Estimates wetness for N velocities of the body between vmin and vmax, and return a matrix
// with the velocities as the first colunmn and the respective wetness as the second column
vector<vector<double>> Simulate(vector<double> box, Body& body, vector<double> rain_v, double vmin,
								double vmax, unsigned int N, double dx);

// Estimates wetness for N velocities of the dynamic body between vmin and vmax, and return
// a matrix with the velocities as the first colunmn and the respective wetness as the second column
vector<vector<double>> Simulate(vector<double> box, Body& body, vector<double> rain_v, double vmin,
								double vmax, unsigned int N, double dx, double tmin, double tmax,
								unsigned int nstep);

// Estimates wetness for N values of dx between dxmin and dxmax, and return a matrix with dx
// as the first colunmn and the respective wetness as the second column
vector<vector<double>> SimErr(vector<double> box, Body& body, vector<double> rain_vel,
							  double body_vel, unsigned int N, double dxmin, double dxmax);

// Estimates wetness for N values of dx between dxmin and dxmax for dynamic body, and return
// a matrix with dx as the first colunmn and the respective wetness as the second column
vector<vector<double>> SimErr(vector<double> box, Body& body, vector<double> rain_vel,
							  double body_vel, unsigned int N, double dxmin, double dxmax,
							  double tmin, double tmax, unsigned int nstep);

// Estimates wetness for N values of nstep between nstepmin and nstepmax, and return a
// matrix with nstep as the first colunmn and the respective wetness as the second column
vector<vector<double>> SimErrT(vector<double> box, Body& body, vector<double> rain_vel,
							   double body_vel, double dx, double tmin, double tmax, unsigned int N,
							   unsigned int nstepmin, unsigned int nstepmax);

// Estimates wetness for N_v velocities of the dynamic body between vmin and vmax, and
// nN_nstep values of nstep between nstep_min and nstep_max and return a matrix with the nsteps
// velocities as the first colunmn, the velocities as the second column and the respective wetness
// as the third column
vector<vector<double>> SimulateNstep(vector<double> box, Body& body, vector<double> rain_v,
									 double vmin, double vmax, unsigned int N_v, double dx,
									 unsigned int nstep_min, unsigned int nstep_max,
									 unsigned int N_nstep);

// Fits points with a parabola y = k(x - x0)^2 + y0, using least squares minimization, and return a
// tuple containing (k, k_std, x0, x0_std, y0, y0_std)
tuple<double, double, double, double, double, double> ParabolicFit(vector<double> x_vals,
																   vector<double> y_vals);

// Finds minimums of wetness using Brent algorithm, calculates wetness for nfit values spaced
// dv around it, and return a tuple containing the optimal velocity, its error, the the minimum rain,
// its error, and a matrix containing the fit points, each row is a point, in the first colum are velocities
// and in the second the wetnesses
tuple<double, double, double, double, vector<vector<double>>>
MinFit(vector<double> box, Body& body, double vmax, double dx, unsigned int nstep, double vcross,
	   double vtail, int nfit, double dv);

// Finds minimums of wetness for a fixed vcross and [vtail_min,
// vtail_max] using Brent algorithm, and calculates wetness for nfit values
// spaced dv around it, return all these values
vector<vector<double>> FindMinFit(vector<double> box, Body& body, double vmax, double dx,
								  unsigned int nstep, double vcross, double vtail_min,
								  double vtail_max, unsigned int n_tail, int nfit, double dv);

// Finds minimums of wetness for a fixed vcross and vtail_min using Brent algorithm with
// nstep in [nstep_min, nstep_max], and calculates wetness for nfit values spaced dv around it,
// return all these values
vector<vector<double>> FindMinFitNstep(vector<double> box, Body& body, double vmax, double dx,
									   unsigned int nstep_min, unsigned int nstep_max,
									   unsigned int N_nstep, double vcross, double vtail, int nfit,
									   double dv);

// Finds minimums of wetness for a fixed vcross and [vtail_min, vtail_max]x[vcross_min,
// vcross_max] with brent, calculates wetness for nfit values around it, return all these values
vector<vector<double>> OptMapFit(vector<double> box, Body& body, double vmax, double dx,
								 unsigned int nstep, unsigned int nfit, double dv, double vtail_min,
								 double vtail_max, unsigned int n_tail, double vcross_min,
								 double vcross_max, unsigned int n_cross);

// Write header file for results of minimization with varying vtail and vcross
void WriteHeadRes(ostream& out, string bodyName, double vmax, double dx, int nstep, int nfit,
				  double dv);

// Write header file for fit points of minimization with varying vtail and vcross
void WriteHeadFit(ostream& out, string bodyName, double vmax, double dx, int nstep, int nfit,
				  double dv);

// Perform a test calcultation and check that results are right
bool AllGood();

// Parse a vector of doubles from a string
vector<double> parseDoubles(const char* text);

// Try to get a vector of three doubles from an XMLElement
vector<double> SafeGet3Vec(XMLElement* parent, const char* tag, string bodyName);

// Try to get text from an XMLElement, throw error if it is null
string SafeGetText(XMLElement* parent, const char* tag);

#endif