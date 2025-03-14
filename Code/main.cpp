#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "VectorStat.h"
#include "body.h"
#include "mins.h"
#include "ray.h"

using namespace std;

int main(int argc, char* argv[]) {
	// Declare stuff
	vector<double> rain_vel(3);
	vector<double> rel_vel(3);
	vector<double> box(3);
	double body_vel, dx;
	unsigned int nstep_v, nstep_t;
	int N_fit;
	double dv;
	vector<vector<double>> results;
	double dxmin;
	int nstepmax;
	double voptW;
	double voptR;
	int nstep_min;
	int nstep_max;
	int N_nstep;

	// Reads the input file
	ifstream ReadInput;
	if (!ReadInput) {
		std::cerr << "Error opening input file." << std::endl;
		return 1;
	}
	ReadInput.open("input.in");
	ReadInput >> box[0] >> box[1] >> box[2];
	ReadInput >> dx;
	ReadInput >> nstep_t;
	ReadInput >> body_vel;
	ReadInput >> nstep_v;
	ReadInput >> rain_vel[0] >> rain_vel[1] >> rain_vel[2];
	ReadInput.close();

	rel_vel = rain_vel;
	rel_vel[0] -= body_vel;

	// Outputs settings and stuff
	cout << "Rain velocity = [ " << rain_vel[0] << ", " << rain_vel[1] << ", " << rain_vel[2]
		 << " ], " << " dx = " << dx << ", nstep_t = " << nstep_t << endl;
	cout << "Object velocity = " << body_vel << endl;
	cout << "Relative velocity = [ " << rel_vel[0] << ", " << rel_vel[1] << ", " << rel_vel[2]
		 << " ]" << endl;

	// Base bodies
	ManyBody TrialS("../Bodies/Sphere.in");
	ManyBody TrialP("../Bodies/Parallelepiped.in");
	ManyBody TrialC("../Bodies/Capsule.in");
	box = {1.1, 1.1, 1.1};

	// Walking man
	ManyBody Walk("../Bodies/WalkingMan.in");
	vector<double> boxW = {0.90, 0.56, 1.74};

	// Running man
	ManyBody Run("../Bodies/RunningMan.in");
	vector<double> boxR = {1.13, 0.58, 1.82};

	// Dynamic Sphere
	ManyBody DynSphere("../Bodies/DynamicSphere.in");
	vector<double> boxDS = {1, 1, 2};

	// Check that body is moving
	cout << "Still: " << WetnessSmooth(boxW, Walk, {0.5, 0.15, -1}, 0.60, 0.01) << endl;
	cout << "Moving: " << WetnessSmooth(boxW, Walk, {0.5, 0.15, -1}, 0.60, 0.01, 0, 1, 50) << endl;

	// // Check boxes
	// PrintDynShadow(boxW, Walk, {0, 0, -1}, dx, 0, 1, 60, "../data/Walk/Proj/Walk_xy");
	// PrintDynShadow(boxW, Walk, {0, -100, -1}, dx, 0, 1, 60, "../data/Walk/Proj/Walk_xz");
	// PrintDynShadow(boxW, Walk, {-100, 0, -1}, dx, 0, 1, 60, "../data/Walk/Proj/Walk_yz");
	// PrintDynShadow(boxR, Run, {0, 0, -1}, dx, 0, 1, 60, "../data/Run/Proj/Run_xy");
	// PrintDynShadow(boxR, Run, {0, -100, -1}, dx, 0, 1, 60, "../data/Run/Proj/Run_xz");
	// PrintDynShadow(boxR, Run, {-100, 0, -1}, dx, 0, 1, 60, "../data/Run/Proj/Run_yz");

	// // Print dynamic state
	// PrintDynState(Walk, 0, 1, 60, "../data/Walk/Status/Walk");
	// PrintDynState(Run, 0, 1, 60, "../data/Run/Status/Run");
	// PrintDynState(DynSphere, 0, 1, 60, "../data/Sphere/Status/DynSphere");

	// // Error analysis with discrete points
	// dx = 0.001;
	// dxmin = 0.0001;
	// nstep_t = 50;
	// nstepmax = 1000;
	// rain_vel = {0.5, 0.25, -1};
	// body_vel = 2.;

	// results = SimErr(box, TrialS, rain_vel, body_vel, 200, dxmin, 0.1);
	// Print("../data/Sphere/ErrorS.dat", results, 15);
	// results = Simulate(box, TrialS, rain_vel, 1, 10, nstep_v, dx);
	// Print("../data/Sphere/CompareS.dat", results, 15);

	// results = SimErr(box, TrialP, rain_vel, body_vel, 200, dxmin, 0.1);
	// Print("../data/Parallelepiped/ErrorP.dat", results, 15);
	// results = Simulate(box, TrialP, rain_vel, 1, 10, nstep_v, dx);
	// Print("../data/Parallelepiped/CompareP.dat", results, 15);

	// results = SimErr(box, TrialC, rain_vel, body_vel, 200, dxmin, 0.1);
	// Print("../data/Capsule/ErrorC.dat", results, 15);
	// results = Simulate(box, TrialC, rain_vel, 1, 10, nstep_v, dx);
	// Print("../data/Capsule/CompareC.dat", results, 15);

	// results = SimErr(boxW, Walk, rain_vel, body_vel, 100, dxmin, 0.1, 0, 1, nstep_t);
	// Print("../data/Walk/ErrDx.dat", results, 12);
	// results = SimErrT(boxW, Walk, rain_vel, body_vel, dx, 0, 1, 100, 1, nstepmax);
	// Print("../data/Walk/ErrT.dat", results, 12);

	// results = SimErr(boxR, Run, rain_vel, body_vel, 100, dxmin, 0.1, 0, 1, nstep_t);
	// Print("../data/Run/ErrDx.dat", results, 12);
	// results = SimErrT(boxR, Run, rain_vel, body_vel, dx, 0, 1, 100, 1, nstepmax);
	// Print("../data/Run/ErrT.dat", results, 12);

	// // Error analysis smooth
	// dx = 0.001;
	// dxmin = 0.0001;
	// nstep_t = 50;
	// nstepmax = 1000;
	// rain_vel = {0.5, 0.25, -1};
	// body_vel = 2.;

	// results = SimSmoothErr(box, TrialS, rain_vel, body_vel, 200, dxmin, 0.1);
	// Print("../data/Sphere/ErrorS_Smooth.dat", results, 15);
	// results = SimulateSmooth(box, TrialS, rain_vel, 1, 10, nstep_v, dx);
	// Print("../data/Sphere/CompareS_Smooth.dat", results, 15);

	// results = SimSmoothErr(box, TrialP, rain_vel, body_vel, 200, dxmin, 0.1);
	// Print("../data/Parallelepiped/ErrorP_Smooth.dat", results, 15);
	// results = SimulateSmooth(box, TrialP, rain_vel, 1, 10, nstep_v, dx);
	// Print("../data/Parallelepiped/CompareP_Smooth.dat", results, 15);

	// results = SimSmoothErr(box, TrialC, rain_vel, body_vel, 200, dxmin, 0.1);
	// Print("../data/Capsule/ErrorC_Smooth.dat", results, 15);
	// results = SimulateSmooth(box, TrialC, rain_vel, 1, 10, nstep_v, dx);
	// Print("../data/Capsule/CompareC_Smooth.dat", results, 15);

	// results = SimSmoothErr(boxW, Walk, rain_vel, body_vel, 100, dxmin, 0.1, 0, 1, nstep_t);
	// Print("../data/Walk/ErrDx_Smooth.dat", results, 12);
	// results = SimSmoothErrT(boxW, Walk, rain_vel, body_vel, dx, 0, 1, 100, 1, nstepmax);
	// Print("../data/Walk/ErrT_Smooth.dat", results, 12);

	// results = SimSmoothErr(boxR, Run, rain_vel, body_vel, 100, dxmin, 0.1, 0, 1, nstep_t);
	// Print("../data/Run/ErrDx_Smooth.dat", results, 12);
	// results = SimSmoothErrT(boxR, Run, rain_vel, body_vel, dx, 0, 1, 100, 1, nstepmax);
	// Print("../data/Run/ErrT_Smooth.dat", results, 12);

	// // Draw shadow of capsules
	// ManyBody OneCap("../Bodies/Capsule.in");
	// ManyBody TwoCap("../Bodies/DoubleCapsule.in");

	// ProjSurface Canv(box, {0, 0, -1}, dx);
	// Canv.BodyProj(OneCap);
	// Canv.PrintRaysFlat("../data/Capsule/OneCapProj.dat");
	// Canv.BodyProj(TwoCap);
	// Canv.PrintRaysFlat("../data/Capsule/TwoCapProj.dat");

	// // Simulation of two Parallelepipeds compenetrating
	// ManyBody Trial2P("../Bodies/DoubleParallelepiped.in");
	// rain_vel = {0, 0.03, -1};
	// body_vel = 0.52;
	// vector<double> dist;
	// vector<double> wet2P;
	// for (size_t i = 0; i < nstep_t; i++) {
	// 	Trial2P.Move(asin((double)i / (nstep_t - 1)) / (2 * M_PI)); // it just works ;)
	// 	vector<double> cent1 = dynamic_cast<Parallelepiped*>(Trial2P.Find("Still"))->GetCent();
	// 	vector<double> cent2 = dynamic_cast<Parallelepiped*>(Trial2P.Find("Moving"))->GetCent();
	// 	dist.push_back(Norm(cent1 - cent2));
	// 	wet2P.push_back(WetnessSmooth(box, Trial2P, rain_vel, body_vel, 0.001));
	// }
	// results = {dist, wet2P};
	// results = Transpose(results);
	// Print("../data/Parallelepiped/DoubleP.dat", results, 12);

	// // For fitting graph
	// vector<vector<double>> RunFitGraph = Simulate(boxR, Run, rain_vel, 0, 2, 9, dx, 0, 1, nstep_t);
	// Print("../data/Run/GraphFit.dat", RunFitGraph, 12);

	// // Wetness around minimum
	// dx = 0.001;
	// nstep_t = 10;
	// voptW = 0.59191;
	// results = Simulate(boxW, Walk, vector<double>{0.5, 0.15, -1}, voptW - 0.03, voptW + 0.03, 600,
	// 				   dx, 0, 1, nstep_t);
	// Print("../data/Walk/Wetfun.dat", results, 12);

	// voptR = 0.884849;
	// results = Simulate(boxR, Run, vector<double>{0.45, 0.3, -1}, voptR - 0.015, voptR + 0.015, 300,
	// 				   dx, 0, 1, nstep_t);
	// Print("../data/Run/Wetfun.dat", results, 12);

	// dx = 0.0005;
	// nstep_t = 20;
	// voptW = 0.59;
	// results = Simulate(boxW, Walk, vector<double>{0.5, 0.15, -1}, voptW - 0.03, voptW + 0.03, 600,
	// 				   dx, 0, 1, nstep_t);
	// Print("../data/Walk/Wetfun_dx_dt.dat", results, 12);
	// voptR = 0.885;
	// results = Simulate(boxR, Run, vector<double>{0.45, 0.3, -1}, voptR - 0.015, voptR + 0.015, 300,
	// 				   dx, 0, 1, nstep_t);
	// Print("../data/Run/Wetfun_dx_dt.dat", results, 12);

	// dx = 0.0005;
	// nstep_t = 10;
	// voptW = 0.59;
	// results = Simulate(boxW, Walk, vector<double>{0.5, 0.15, -1}, voptW - 0.03, voptW + 0.03, 600,
	// 				   dx, 0, 1, nstep_t);
	// Print("../data/Walk/Wetfun_dx.dat", results, 12);
	// voptR = 0.885;
	// results = Simulate(boxR, Run, vector<double>{0.45, 0.3, -1}, voptR - 0.015, voptR + 0.015, 300,
	// 				   dx, 0, 1, nstep_t);
	// Print("../data/Run/Wetfun_dx.dat", results, 12);

	// dx = 0.001;
	// nstep_t = 20;
	// voptW = 0.59;
	// results = Simulate(boxW, Walk, vector<double>{0.5, 0.15, -1}, voptW - 0.03, voptW + 0.03, 600,
	// 				   dx, 0, 1, nstep_t);
	// Print("../data/Walk/Wetfun_dt.dat", results, 12);
	// voptR = 0.885;
	// results = Simulate(boxR, Run, vector<double>{0.45, 0.3, -1}, voptR - 0.015, voptR + 0.015, 300,
	// 				   dx, 0, 1, nstep_t);
	// Print("../data/Run/Wetfun_dt.dat", results, 12);

	// // Trials
	// double wet = Wetness(boxR, Run, rain_vel, body_vel, dx, 0, 1, nstep_t);
	// cout << "Wentness: " << setprecision(15) << wet << endl;

	// double wetSmooth = WetnessSmooth(boxR, Run, rain_vel, body_vel, dx, 0, 1, nstep_t);
	// cout << "Smooth wentness: " << setprecision(15) << wetSmooth << endl;

	// // Smooth wetness around minimum
	// dx = 0.001;
	// nstep_t = 10;
	// voptW = 0.59;
	// results = SimulateSmooth(boxW, Walk, vector<double>{0.5, 0.15, -1}, voptW - 0.03, voptW + 0.03,
	// 						 600, dx, 0, 1, nstep_t);
	// Print("../data/Walk/WetfunS.dat", results, 12);
	// voptR = 0.885;
	// results = SimulateSmooth(boxR, Run, vector<double>{0.45, 0.3, -1}, voptR - 0.015, voptR + 0.015,
	// 						 300, dx, 0, 1, nstep_t);
	// Print("../data/Run/WetfunS.dat", results, 12);

	// dx = 0.0005;
	// nstep_t = 20;
	// voptW = 0.59;
	// results = SimulateSmooth(boxW, Walk, vector<double>{0.5, 0.15, -1}, voptW - 0.03, voptW + 0.03,
	// 						 600, dx, 0, 1, nstep_t);
	// Print("../data/Walk/WetfunS_dx_dt.dat", results, 12);
	// voptR = 0.885;
	// results = SimulateSmooth(boxR, Run, vector<double>{0.45, 0.3, -1}, voptR - 0.015, voptR + 0.015,
	// 						 300, dx, 0, 1, nstep_t);
	// Print("../data/Run/WetfunS_dx_dt.dat", results, 12);

	// dx = 0.0005;
	// nstep_t = 10;
	// voptW = 0.59;
	// results = SimulateSmooth(boxW, Walk, vector<double>{0.5, 0.15, -1}, voptW - 0.03, voptW + 0.03,
	// 						 600, dx, 0, 1, nstep_t);
	// Print("../data/Walk/WetfunS_dx.dat", results, 12);
	// voptR = 0.885;
	// results = SimulateSmooth(boxR, Run, vector<double>{0.45, 0.3, -1}, voptR - 0.015, voptR + 0.015,
	// 						 300, dx, 0, 1, nstep_t);
	// Print("../data/Run/WetfunS_dx.dat", results, 12);

	// dx = 0.001;
	// nstep_t = 20;
	// voptW = 0.59;
	// results = SimulateSmooth(boxW, Walk, vector<double>{0.5, 0.15, -1}, voptW - 0.03, voptW + 0.03,
	// 						 600, dx, 0, 1, nstep_t);
	// Print("../data/Walk/WetfunS_dt.dat", results, 12);
	// voptR = 0.885;
	// results = SimulateSmooth(boxR, Run, vector<double>{0.45, 0.3, -1}, voptR - 0.015, voptR + 0.015,
	// 						 300, dx, 0, 1, nstep_t);
	// Print("../data/Run/WetfunS_dt.dat", results, 12);

	// // Brent path
	// Brent mins(0.001);
	// auto wetfuncW = [&boxW, &Walk, dx, nstep_t](double x) {
	// 	return WetnessSmooth(boxW, Walk, vector<double>{0.5, 0.15, -1}, x, dx, 0., 1., nstep_t);
	// };
	// mins.bracket(0., 0.7, 3, wetfuncW);
	// results = mins.minimize(wetfuncW);
	// Print("../data/Walk/MinPath.dat", Transpose(results), 12);
	// auto wetfuncR = [&boxR, &Run, dx, nstep_t](double x) {
	// 	return WetnessSmooth(boxR, Run, vector<double>{0.45, 0.3, -1}, x, dx, 0., 1., nstep_t);
	// };
	// mins.bracket(0., 2., 3, wetfuncR);
	// results = mins.minimize(wetfuncR);
	// Print("../data/Run/MinPath.dat", Transpose(results), 12);

	// // Simulate different nstep
	// dx = 0.001;
	// nstep_min = 2;
	// nstep_max = 100;
	// N_nstep = 10;

	// voptW = 0.59;
	// results = SimulateNstepSmooth(boxW, Walk, vector<double>{0.5, 0.15, -1}, voptW - 0.05,
	// 							  voptW + 0.05, 40, dx, nstep_min, nstep_max, N_nstep);
	// Print("../data/Walk/WetfunS_nstep.dat", results, 12);

	// voptR = 0.89;
	// results = SimulateNstepSmooth(boxR, Run, vector<double>{0.45, 0.3, -1}, voptR - 0.015,
	// 							  voptR + 0.015, 40, dx, nstep_min, nstep_max, N_nstep);
	// Print("../data/Run/WetfunS_nstep.dat", results, 12);

	// double voptDS = 2.872;
	// results = SimulateNstepSmooth(boxDS, DynSphere, vector<double>{0.45, 0.3, -1}, voptDS - 0.03,
	// 							  voptDS + 0.03, 40, dx, nstep_min, nstep_max, N_nstep);
	// Print("../data/Sphere/WetfunS_nstep.dat", results, 12);

	// dx = 0.0005;
	// results = SimulateNstepSmooth(boxDS, DynSphere, vector<double>{0.45, 0.3, -1}, voptDS - 0.03,
	// 							  voptDS + 0.03, 40, dx, nstep_min, nstep_max, N_nstep);
	// Print("../data/Sphere/WetfunS_nstep0005.dat", results, 12);

	// // Check shadows of dynamic sphere
	// nstep_t = 20;
	// dx = 0.001;
	// PrintDynShadowSmooth(boxDS, DynSphere, vector<double>{0.45, 0.3, -1}, dx, 0, 1, nstep_t,
	// 					 "../data/Sphere/Shadows/Smooth001_");
	// dx = 0.002;
	// PrintDynShadowSmooth(boxDS, DynSphere, vector<double>{0.45, 0.3, -1}, dx, 0, 1, nstep_t,
	// 					 "../data/Sphere/Shadows/Smooth002_");

	// // Smooth Brent minimization fit
	// int N_vtail = 50;
	// N_fit = 9;
	// dv = 0.006;
	// dx = 0.001;
	// nstep_t = 50;

	// results = FindMinFitSmooth(boxW, Walk, 0., 0.7, dx, nstep_t, 0., 0., 0.7, N_vtail, N_fit, dv);
	// Print("../data/Walk/OptFitSmoothW0.dat", results, 12);

	// results = FindMinFitSmooth(boxW, Walk, 0., 0.7, dx, nstep_t, 0., 0.236, 0.25, 2, N_fit, dv);
	// Print("../data/Walk/OptFitSmoothW0_Extra.dat", results, 12);

	// results = FindMinFitSmooth(boxW, Walk, 0., 0.7, dx, nstep_t, 0.15, 0., 0.7, N_vtail, N_fit, dv);
	// Print("../data/Walk/OptFitSmoothW015.dat", results, 12);

	// results = FindMinFitSmooth(boxW, Walk, 0., 0.7, dx, nstep_t, 0.15, 0.250, 0.250, 1, N_fit, dv);
	// Print("../data/Walk/OptFitSmoothW015_Extra.dat", results, 12);

	// results = FindMinFitSmooth(boxW, Walk, 0., 0.7, dx, nstep_t, 0.25, 0., 0.7, N_vtail, N_fit, dv);
	// Print("../data/Walk/OptFitSmoothW025.dat", results, 12);

	// results = FindMinFitSmooth(boxW, Walk, 0., 0.7, dx, nstep_t, 0.30, 0., 0.7, N_vtail, N_fit, dv);
	// Print("../data/Walk/OptFitSmoothW030.dat", results, 12);

	// results = FindMinFitSmooth(boxW, Walk, 0., 0.7, dx, nstep_t, 0.30, 0.3785714285, 0.3785714285,
	// 						   1, N_fit, dv);
	// Print("../data/Walk/OptFitSmoothW030_Extra.dat", results, 12);

	// results = FindMinFitSmooth(boxR, Run, 0., 2., dx, nstep_t, 0., 0., 2, N_vtail, N_fit, dv);
	// Print("../data/Run/OptFitSmoothR0.dat", results, 12);

	// results = FindMinFitSmooth(boxR, Run, 0., 2., dx, nstep_t, 0., 0.143, 0.143, 1, N_fit, dv);
	// Print("../data/Run/OptFitSmoothR0_Extra.dat", results, 12);

	// results = FindMinFitSmooth(boxR, Run, 0., 2., dx, nstep_t, 0.3, 0., 2, N_vtail, N_fit, dv);
	// Print("../data/Run/OptFitSmoothR030.dat", results, 12);

	// results = FindMinFitSmooth(boxR, Run, 0., 2., dx, nstep_t, 0.3, 0.183, 0.183, 1, N_fit, dv);
	// Print("../data/Run/OptFitSmoothR030_Extra.dat", results, 12);

	// results = FindMinFitSmooth(boxR, Run, 0., 2., dx, nstep_t, 0.6, 0., 2, N_vtail, N_fit, dv);
	// Print("../data/Run/OptFitSmoothR060.dat", results, 12);

	// results = FindMinFitSmooth(boxR, Run, 0., 2., dx, nstep_t, 0.8, 0., 2, N_vtail, N_fit, dv);
	// Print("../data/Run/OptFitSmoothR080.dat", results, 12);

	// results = FindMinFitSmooth(boxR, Run, 0., 2., dx, nstep_t, 1., 0., 2, N_vtail, N_fit, dv);
	// Print("../data/Run/OptFitSmoothR1.dat", results, 12);

	// results = FindMinFitSmooth(boxW, Walk, 0., 0.7, dx, nstep_t, 0.15, 0.5, 0.7, 1, N_fit, dv);
	// Print("../data/Walk/OptFitSmoothW_Ex.dat", results, 12);

	// results = FindMinFitSmooth(boxR, Run, 0., 2., dx, nstep_t, 0.3, 0.45, 2, 1, N_fit, dv);
	// Print("../data/Run/OptFitSmoothR_Ex.dat", results, 12);

	// // Smooth Brent minimization fit with variable nstep
	// nstep_min = 5;
	// nstep_max = 100;
	// N_nstep = 30;
	// N_fit = 9;
	// dv = 0.006;
	// dx = 0.001;

	// results = FindMinFitSmoothNstep(boxW, Walk, 0., 0.7, dx, nstep_min, nstep_max, N_nstep, 0.15,
	// 								0.5, N_fit, dv);
	// Print("../data/Walk/OptFitSmoothW_nstep.dat", results, 12);

	// results = FindMinFitSmoothNstep(boxR, Run, 0., 2., dx, nstep_min, nstep_max, N_nstep, 0.30,
	// 								0.45, N_fit, dv);
	// Print("../data/Run/OptFitSmoothR_nstep.dat", results, 12);

	// // Smooth Brent minimization fit map
	// nstep_v = 20;
	// N_fit = 9;
	// dv = 0.006;
	// dx = 0.001;
	// nstep_t = 50;

	// results = OptMapFitSmooth(boxW, Walk, 0, 0.7, dx, nstep_t, N_fit, dv, 0, 0.7, nstep_v, 0, 0.40,
	// 						  nstep_v);
	// Print("../data/Walk/OptMapFitSmoothW.dat", results, 12);

	// results =
	// 	OptMapFitSmooth(boxR, Run, 0, 2, dx, nstep_t, N_fit, dv, 0, 2, nstep_v, 0, 1.2, nstep_v);
	// Print("../data/Run/OptMapFitSmoothR.dat", results, 12);

	// results = OptMapFitSmooth(boxW, Walk, 0, 0.7, dx, nstep_t, N_fit, dv, 0., 1.2, nstep_v, 0, 1.2,
	// 						  nstep_v);
	// Print("../data/Walk/OptMapCompareSmoothW.dat", results, 12);

	// results =
	// 	OptMapFitSmooth(boxR, Run, 0, 2, dx, nstep_t, N_fit, dv, 0., 1.2, nstep_v, 0, 1.2, nstep_v);
	// Print("../data/Run/OptMapCompareSmoothR.dat", results, 12);

	// // Points for example fit
	// dx = 0.001;
	// nstep_t = 50;
	// nstep_v = 300;
	// rain_vel = {0.5, 0.15, -1};
	// double vmin = 0.55, vmax = 0.64;

	// results = SimulateSmooth(boxW, Walk, rain_vel, vmin, vmax, nstep_v, dx, 0, 1, nstep_t);
	// Print("../data/Walk/EgWetFun.dat", results, 12);


	// Running with straight torso
	int N_vtail = 50;
	N_fit = 9;
	dv = 0.006;
	dx = 0.001;
	nstep_t = 50;
	ManyBody RunStr("../Bodies/RunningStraight.in");
	// PrintDynState(RunStr, 0, 1, 60, "../data/Run/Status/RunStr");

	// results = FindMinFitSmooth(boxR, RunStr, 0., 2., dx, nstep_t, 0., 0., 2, N_vtail, N_fit, dv);
	// Print("../data/Run/OptFitSmoothRSTR0.dat", results, 12);

	results = FindMinFitSmooth(boxR, RunStr, 0., 2., dx, nstep_t, 0.3, 0., 2, N_vtail, N_fit, dv);
	Print("../data/Run/OptFitSmoothRSTR030.dat", results, 12);
}