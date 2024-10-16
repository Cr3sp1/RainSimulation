#ifndef __RainFunctions_h__
#define __RainFunctions_h__


#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm> 
#include "VectorOperations.h"


using namespace std;

// Forward declaration
class Body;
class Ray;
class ProjSurface;


// Projects the Point on a plane perpendicular to v and passing through p
vector<double> Project( vector<double> Point, vector<double> p, vector<double> v );

// Finds the vertex in the middle of the three seen faces of a parallelepiped defined by a point p and three sides 
vector<double> FindMiddle( vector<double> p, vector<vector<double>> sides, vector<double> v );

// Finds the hexagonal projection H of a parallelepiped on a plane perpendicular to v and passing through P
vector<vector<double>> FindHexProj( vector<double> p, vector<vector<double>> Side, vector<double> v, vector<double> P);

/// Returns the highest absolute value of the projections of the vertices of H on a line in direction u1 passing through H[0]
double MaxU(vector<vector<double>> H, vector<double> u );

// Checks wether a point is inside the hexagon using triangles and baycentric coordinates
bool PointIsInsideT( vector<double> Point, vector<vector<double>> H );

// Periodic Boundary conditions for the index of H
int PBCH( int i );

// Checks rays generation 
void RayGenCheck( string outfile, vector<double> box, vector<double> rel_vel );

// Estimates wetness for N velocities of the body between vmin and vmax, and returns a matrix with the velocities as the first colunmn and the respective wetness as the second column
vector<vector<double>> Simulate( vector<double> box, Body& body, vector<double> rain_v, double vmin, double vmax, unsigned int N, double dx );

// Estimates wetness for N velocities of the dynamic body between vmin and vmax, and returns a matrix with the velocities as the first colunmn and the respective wetness as the second column
vector<vector<double>> Simulate( vector<double> box, Body& body, vector<double> rain_v, double vmin, double vmax, unsigned int N, double dx, double tmin, double tmax, unsigned int nstep );

// Estimates wetness for N values of dx between dxmin and dxmax, and returns a matrix with dx as the first colunmn and the respective wetness as the second column
vector<vector<double>> SimErr( vector<double> box, Body& body, vector<double> relvel, double bodyvel, unsigned int N, double dxmin, double dxmax);

// Estimates wetness for N values of dx between dxmin and dxmax for dynamic body, and returns a matrix with dx as the first colunmn and the respective wetness as the second column
vector<vector<double>> SimErr( vector<double> box, Body& body, vector<double> relvel, double bodyvel, unsigned int N, double dxmin, double dxmax, double tmin, double tmax, unsigned int nstep);

// Estimates wetness for N values of nstep between nstepmin and nstepmax, and returns a matrix with nstep as the first colunmn and the respective wetness as the second column
vector<vector<double>> SimErrT( vector<double> box, Body& body, vector<double> relvel, double bodyvel, double dx, double tmin, double tmax, unsigned int N, unsigned int nstepmin, unsigned int nstepmax);

// Estimates wetness for N_dx values of dx between dxmin and dxmax, and N_T values of nstep between nstepmin and nstepmax, and returns a matrix where first column is the nstep and first row is the dx
vector<vector<double>> SimErrTdx( vector<double> box, Body& body, vector<double> relvel, double bodyvel, unsigned int N_dx, double dxmin, double dxmax, unsigned int N_T, unsigned int nstepmin, unsigned int nstepmax);

// Estimates wetness for N velocities of the body between vmin and vmax, and returns a matrix with the velocities as the first colunmn and the respective theorical wetness as the second column and the estimated wetness as the third
vector<vector<double>> CompareAN( vector<double> box, Body& body, vector<double> rain_v, double vmin, double vmax, unsigned int N, double dx);

// Estimates wetness for N velocities of two body between vmin and vmax, and returns a matrix with the velocities as the first colunmn and the wetness of the first body as the second column and of the second body as the third column
vector<vector<double>> CompareBB( vector<double> box, Body& body1, Body& body2, vector<double> rain_v, double vmin, double vmax, unsigned int N, double dx);

// Returns the minimum distance between the point p and the segment line with extremes l1 and l2
double PointSegDist( vector<double> p, vector<double> l1, vector<double> l2 );

// Returns NxN identity matrix
vector<vector<double>> IdMat( unsigned int N );

// Returns the rotation matrix
vector<vector<double>> RotMat( vector<double> axis, double theta );

// Rotates a Point relative to the point Rot0
void Rotate( vector<double>& Point, const vector<double>& Rot0, const vector<vector<double>>& Rotmat );

// Prints the shadow of a body at nstep different time steps in [tmin, tmax)
void PrintDynShadow( vector<double> box, Body& body, vector<double> relvel, double dx, double tmin, double tmax, unsigned int nstep, string outfile );

// Prints the state of a body at nstep different time steps in [tmin, tmax)
void PrintDynState( Body& body, double tmin, double tmax, unsigned int nstep, string outfile );

// Looks for the minimum of wetness between vmin and vmax with N steps and returns its value if it finds it, else returns -1.
double FindMin( vector<double> box, Body& body, vector<double> rain_v, double vmin, double vmax, unsigned int N, double dx, unsigned int nstep );

// Finds minimums of wetness on a square lattice in the space of coordinates [vtail_min, vtail_max]x[vcross_min, vcross_max], returns a matrix with vtail as first column, vcross as second and in the third column best vb, or -1 if it doesn't exist
vector<vector<double>> OptMap( vector<double> box, Body& body, double vmin, double vmax, unsigned int N, double dx, unsigned int nstep, double vtail_min, double vtail_max, unsigned int n_tail, double vcross_min, double vcross_max, unsigned int n_cross );

// Looks for the minimum of wetness for two bodies between vmin and vmax with N steps and returns in result[0] its value, and in result[1] returns number corresponding to the body with less wetness
vector<double> FindMinCompare( vector<double> box1, Body& body1, vector<double> box2, Body& body2, vector<double> rain_v, double vmin, double vmax1, double vmax2, unsigned int N, double dx, unsigned int nstep );

// Finds minimums of wetness between two bodies on a square lattice in the space of coordinates [vtail_min, vtail_max]x[vcross_min, vcross_max], returns a matrix with vtail as first column, vcross as second and in the third column best vb, on the fourth column the number corresponding to the body with less wetness
vector<vector<double>> OptMapCompare( vector<double> box1, Body& body1, vector<double> box2, Body& body2, double vmin, double vmax1, double vmax2, unsigned int N, double dx, unsigned int nstep, double vtail_min, double vtail_max, unsigned int n_tail, double vcross_min, double vcross_max, unsigned int n_cross );

// Finds minimums of wetness for a fixed vcross and [vtail_min, vtail_max], and calculates wetness for n_fit values around it, returns all these values
vector<vector<double>> FindMinFit(vector<double> box, Body& body, double vmin, double vmax, unsigned int N, double dx, unsigned int nstep, unsigned int n_fit, double vcross, double vtail_min, double vtail_max, unsigned int n_tail );

// Finds minimums of wetness for a fixed vcross and [vtail_min, vtail_max]x[vcross_min, vcross_max], and calculates wetness for n_fit values around it, returns all these values
vector<vector<double>> OptMapFit(vector<double> box, Body& body, double vmin, double vmax, unsigned int N, double dx, unsigned int nstep, unsigned int n_fit, double vtail_min, double vtail_max, unsigned int n_tail, double vcross_min, double vcross_max, unsigned int n_cross );

// Evaluates the wetness for N vb in the range [vmin, vmax] for a set vtail and vcross
vector<vector<double>> WetFit(vector<double> box, Body& body, double vmin, double vmax, unsigned int N, double dx, unsigned int nstep, double vtail, double vcross );


#endif