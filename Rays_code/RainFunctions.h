#ifndef __RainFunctions_h__
#define __RainFunctions_h__


#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include "VectorOperations.h"


using namespace std;

// Forward declaration
class Body;
class Ray;
class ProjSurface;


// Projects the Point on a plane perpendicular to v and passing through p
vector<long double> Project( vector<long double> Point, vector<long double> p, vector<long double> v );

// Finds the vertex in the middle of the three seen faces of a parallelepiped defined by a pont and three sides 
vector<long double> FindMiddle( vector<long double> p, vector<vector<long double>> sides, vector<long double> v );

// Finds the hexagonal projection H of a parallelepiped on a plane perpendicular to v and passing through P
vector<vector<long double>> FindHexProj(  vector<long double> p, vector<vector<long double>> Side, vector<long double> v, vector<long double> P);

/// Returns the highest absolute value of the projections of the vertices of H on a line in direction u1 passing through H[0]
long double MaxU(vector<vector<long double>> H, vector<long double> u );

// (NOTE: DOESN'T ALWAYS WORK)  Auxiliary function used only in the surface constructor that checks wether a point is inside the hexagon checking parallelograms
// bool PointIsInsideP( vector<long double> Point, vector<vector<long double>> H );

// Auxiliary function used only in the surface constructor that checks wether a point is inside the hexagon checking triangles
bool PointIsInsideT( vector<long double> Point, vector<vector<long double>> H );

// Periodic Boundary conditions for the index of H
int PBCH( int i );

// Checks rays generation 
void RayGenCheck( string outfile, vector<long double> box, vector<long double> rel_vel );

// Estimates wetness for N velocities of the body between vmin and vmax, and returns a matrix with the velocities as the first colunmn and the respective wetness as the second column
vector<vector<long double>> Simulate( vector<long double> box, Body& body, vector<long double> rain_v, long double vmin, long double vmax, unsigned int N, long double dx );

// Estimates wetness for N velocities of the body between vmin and vmax, and returns a matrix with the velocities as the first colunmn and the respective theorical wetness as the second column and the estimated wetness as the third
vector<vector<long double>> CompareAN( vector<long double> box, Body& body, vector<long double> rain_v, long double vmin, long double vmax, unsigned int N, long double dx);

// Returns the minimum distance between the point p and the segment line with extremes l1 and l2
long double PointSegDist( vector<long double> p, vector<long double> l1, vector<long double> l2 );

#endif