#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <random> 
#include "body.h"
#include "ray.h"
#include "VectorStat.h"
#include "mins.h"


using namespace std;


int main (int argc, char *argv[]){

    // Declare stuff
    vector<double> rain_vel(3);
    vector<double> rel_vel(3);
    vector<double> box(3); 
    double body_vel, dx;
    unsigned int nstep_v, nstep_t;
    // int N_vb, N_fit;
    

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
    cout << "Rain velocity = [ " << rain_vel[0] << ", " << rain_vel[1] << ", " << rain_vel[2] << " ], " << " dx = " << dx << ", nstep_t = " << nstep_t << endl;
    cout << "Object velocity = " << body_vel << endl;
    cout << "Relative velocity = [ " << rel_vel[0] << ", " << rel_vel[1] << ", " << rel_vel[2] << " ]" << endl;


    // Base bodies
    ManyBody TrialS("../Bodies/Sphere.in");
    ManyBody TrialP("../Bodies/Parallelepiped.in");
    ManyBody TrialC("../Bodies/Capsule.in");
    box = {1.1, 1.1, 1.1};

    // Walking man
    ManyBody Walk("../Bodies/WalkingMan.in");
    vector<double> boxW = {0.90, 0.56, 1.74 };

    // Running man
    ManyBody Run("../Bodies/RunningMan.in");
    vector<double> boxR = {1.13, 0.58, 1.82};

    // Dynamic Sphere
    ManyBody DynSphere("../Bodies/DynamicSphere.in");
    vector<double> boxDS = {1, 1, 2};

    // Check that body is moving
    cout << "Still: " << WetnessSmooth(boxW, Walk, {0.5, 0.15, -1}, 0.60, 0.01 ) << endl;
    cout << "Moving: " << WetnessSmooth(boxW, Walk, {0.5, 0.15, -1}, 0.60, 0.01, 0, 1, 50 ) << endl;


    // // Check boxes
    // PrintDynShadow(boxW, Walk, {0, 0, -1}, dx, 0, 1, 60, "../data/Walk/Proj/Walk_xy" );
    // PrintDynShadow(boxW, Walk, {0, -100, -1}, dx, 0, 1, 60, "../data/Walk/Proj/Walk_xz" );
    // PrintDynShadow(boxW, Walk, {-100, 0, -1}, dx, 0, 1, 60, "../data/Walk/Proj/Walk_yz" );
    // PrintDynShadow(boxR, Run, {0, 0, -1}, dx, 0, 1, 60, "../data/Run/Proj/Run_xy" );
    // PrintDynShadow(boxR, Run, {0, -100, -1}, dx, 0, 1, 60, "../data/Run/Proj/Run_xz" );
    // PrintDynShadow(boxR, Run, {-100, 0, -1}, dx, 0, 1, 60, "../data/Run/Proj/Run_yz" );

    // // Print dynamic state
    // PrintDynState( Walk, 0, 1, 60, "../data/Walk/Status/Walk");
    // PrintDynState( Run, 0, 1, 60, "../data/Run/Status/Run");
    // PrintDynState( DynSphere, 0, 1, 60, "../data/Sphere/Status/DynSphere");



    // // Error analysis with discrete points
    // dx = 0.001;
    // double dxmin = 0.0001;
    // nstep_t = 50;
    // int nstepmax = 1000;
    // rain_vel = {0.5, 0.25, -1};
    // body_vel = 2.;

    // vector<vector<double>> resultsS;
    // resultsS = SimErr( box, TrialS, rain_vel, body_vel, 200, dxmin, 0.1 );
    // Print( "../data/Sphere/ErrorS.dat", resultsS, 15);
    // resultsS = Simulate( box, TrialS, rain_vel, 1, 10, nstep_v, dx );
    // Print("../data/Sphere/CompareS.dat", resultsS, 15);
        
    // vector<vector<double>> resultsP;
    // resultsP = SimErr( box, TrialP, rain_vel, body_vel, 200, dxmin, 0.1 );
    // Print( "../data/Parallelepiped/ErrorP.dat", resultsP, 15);
    // resultsP = Simulate( box, TrialP, rain_vel, 1, 10, nstep_v, dx );
    // Print("../data/Parallelepiped/CompareP.dat", resultsP, 15);

    // vector<vector<double>> resultsC;
    // resultsC = SimErr( box, TrialC, rain_vel, body_vel, 200, dxmin, 0.1 );
    // Print( "../data/Capsule/ErrorC.dat", resultsC, 15);
    // resultsC = Simulate( box, TrialC, rain_vel, 1, 10, nstep_v, dx );
    // Print("../data/Capsule/CompareC.dat", resultsC, 15);

    // vector<vector<double>> WalkResDx = SimErr(boxW, Walk, rain_vel, body_vel, 100, dxmin, 0.1, 0, 1, nstep_t);
    // Print( "../data/Walk/ErrDx.dat", WalkResDx, 12 );
    // vector<vector<double>> WalkResT = SimErrT(boxW, Walk, rain_vel, body_vel, dx, 0, 1, 100, 1, nstepmax);
    // Print( "../data/Walk/ErrT.dat", WalkResT, 12 );

    // vector<vector<double>> RunResDx = SimErr(boxR, Run, rain_vel, body_vel, 100, dxmin, 0.1, 0, 1, nstep_t);
    // Print( "../data/Run/ErrDx.dat", RunResDx, 12 );
    // vector<vector<double>> RunResT = SimErrT(boxR, Run, rain_vel, body_vel, dx, 0, 1, 100, 1, nstepmax);
    // Print( "../data/Run/ErrT.dat", RunResT, 12 );


    // Error analysis smooth
    // dx = 0.001;
    // double dxmin = 0.0001;
    // nstep_t = 50;
    // int nstepmax = 1000;
    // rain_vel = {0.5, 0.25, -1};
    // body_vel = 2.;

    // vector<vector<double>> resultsS;
    // resultsS = SimSmoothErr( box, TrialS, rain_vel, body_vel, 200, dxmin, 0.1 );
    // Print( "../data/Sphere/ErrorS_Smooth.dat", resultsS, 15);
    // resultsS = SimulateSmooth( box, TrialS, rain_vel, 1, 10, nstep_v, dx );
    // Print("../data/Sphere/CompareS_Smooth.dat", resultsS, 15);
        
    // vector<vector<double>> resultsP;
    // resultsP = SimSmoothErr( box, TrialP, rain_vel, body_vel, 200, dxmin, 0.1 );
    // Print( "../data/Parallelepiped/ErrorP_Smooth.dat", resultsP, 15);
    // resultsP = SimulateSmooth( box, TrialP, rain_vel, 1, 10, nstep_v, dx );
    // Print("../data/Parallelepiped/CompareP_Smooth.dat", resultsP, 15);

    // vector<vector<double>> resultsC;
    // resultsC = SimSmoothErr( box, TrialC, rain_vel, body_vel, 200, dxmin, 0.1 );
    // Print( "../data/Capsule/ErrorC_Smooth.dat", resultsC, 15);
    // resultsC = SimulateSmooth( box, TrialC, rain_vel, 1, 10, nstep_v, dx );
    // Print("../data/Capsule/CompareC_Smooth.dat", resultsC, 15);

    // vector<vector<double>> WalkResDx = SimSmoothErr(boxW, Walk, rain_vel, body_vel, 100, dxmin, 0.1, 0, 1, nstep_t);
    // Print( "../data/Walk/ErrDx_Smooth.dat", WalkResDx, 12 );
    // vector<vector<double>> WalkResT = SimSmoothErrT(boxW, Walk, rain_vel, body_vel, dx, 0, 1, 100, 1, nstepmax);
    // Print( "../data/Walk/ErrT_Smooth.dat", WalkResT, 12 );

    // vector<vector<double>> RunResDx = SimSmoothErr(boxR, Run, rain_vel, body_vel, 100, dxmin, 0.1, 0, 1, nstep_t);
    // Print( "../data/Run/ErrDx_Smooth.dat", RunResDx, 12 );
    // vector<vector<double>> RunResT = SimSmoothErrT(boxR, Run, rain_vel, body_vel, dx, 0, 1, 100, 1, nstepmax);
    // Print( "../data/Run/ErrT_Smooth.dat", RunResT, 12 );



    // // Draw shadow of capsules
    // ManyBody OneCap("../Bodies/Capsule.in");
    // ManyBody TwoCap("../Bodies/DoubleCapsule.in");

    // ProjSurface Canv( box, {0, 0, -1}, dx );
    // Canv.BodyProj(OneCap);
    // Canv.PrintRaysFlat("../data/Capsule/OneCapProj.dat");
    // Canv.BodyProj(TwoCap);
    // Canv.PrintRaysFlat("../data/Capsule/TwoCapProj.dat");
    


    // // Simulation of two Parallelepipeds compenetrating
    // ManyBody Trial2P("../Bodies/DoubleParallelepiped.in");
    // rain_vel = {0, 0.03, -1 };
    // body_vel = 0.52;
    // vector<double> dist;
    // vector<double> wet2P;
    // for( size_t i = 0; i < nstep_t; i++ ) {
    //     Trial2P.Move( asin( (double)i/(nstep_t-1))/(2*M_PI));   // it just works ;)
    //     vector<double> cent1 = dynamic_cast<Parallelepiped*>(Trial2P.Find("Still"))->GetCent();
    //     vector<double> cent2 = dynamic_cast<Parallelepiped*>(Trial2P.Find("Moving"))->GetCent();
    //     dist.push_back( Norm( cent1 - cent2 ));
    //     wet2P.push_back( WetnessSmooth( box, Trial2P, rain_vel, body_vel, 0.001 ));
    // }
    // vector<vector<double>> results2P = { dist, wet2P };
    // results2P = Transpose(results2P);
    // Print("../data/Parallelepiped/DoubleP.dat", results2P, 12 );
    


    // // Find v_opt fit map
    // vector<vector<double>> WalkOptMapFit = OptMapFit( boxW, Walk, 0, 0.7, N_vb, dx, nstep_t, N_fit,  0, 0.7, nstep_v, 0, 0.35, nstep_v );
    // Print( "../data/Walk/OptMapFitW.dat", WalkOptMapFit, 12 );

    // vector<vector<double>> RunOptMapFit = OptMapFit( boxR, Run, 0, 2, N_vb, dx, nstep_t, N_fit, 0, 2, nstep_v, 0, 1.2, nstep_v );
    // Print( "../data/Run/OptMapFitR.dat", RunOptMapFit, 12 ); 



    // Compare walking and running
    // vector<vector<double>> WalkOptMapComp = OptMapFit( boxW, Walk, 0, 0.7, N_vb, dx, nstep_t, N_fit,  -0.5, 2, nstep_v, 0, 1.2, nstep_v );
    // Print( "../data/Walk/OptMapCompareW.dat", WalkOptMapComp, 12 );

    // vector<vector<double>> RunOptMapComp = OptMapFit( boxR, Run, 0, 2, N_vb, dx, nstep_t, N_fit,  -0.5, 2, nstep_v, 0, 1.2, nstep_v );
    // Print( "../data/Run/OptMapCompareR.dat", RunOptMapComp, 12 );



    // // V_opt fit graphs
    // N_vb = 50;
    // N_fit = 5;
    // vector<vector<double>> WalkMins1 = FindMinFit( boxW, Walk, 0, 0.7, N_vb, dx, nstep_t, N_fit, 0., 0., 0.7, nstep_v );
    // Print( "../data/Walk/OptFitW0.dat", WalkMins1, 12 );

    // vector<vector<double>> WalkMins2 = FindMinFit( boxW, Walk, 0, 0.7, N_vb, dx, nstep_t, N_fit,  0.15, 0., 0.7, nstep_v );
    // Print( "../data/Walk/OptFitW015.dat", WalkMins2, 12 );

    // vector<vector<double>> WalkMins3 = FindMinFit( boxW, Walk, 0, 0.7, N_vb, dx, nstep_t, N_fit, 0.25, 0., 0.7, nstep_v );
    // Print( "../data/Walk/OptFitW025.dat", WalkMins3, 12 );

    // vector<vector<double>> WalkMins4 = FindMinFit( boxW, Walk, 0, 0.7, N_vb, dx, nstep_t, N_fit, 0.30, 0., 0.7, nstep_v );
    // Print( "../data/Walk/OptFitW030.dat", WalkMins4, 12 );

    // vector<vector<double>> RunMins1 = FindMinFit( boxR, Run, 0, 2, N_vb, dx, nstep_t, N_fit,  0, 0, 2, nstep_v );
    // Print( "../data/Run/OptFitR0.dat", RunMins1, 12 );

    // vector<vector<double>> RunMins2 = FindMinFit( boxR, Run, 0, 2, N_vb, dx, nstep_t, N_fit,  0.5, 0, 2, nstep_v );
    // Print( "../data/Run/OptFitR05.dat", RunMins2, 12 );

    // vector<vector<double>> RunMins3 = FindMinFit( boxR, Run, 0, 2, N_vb, dx, nstep_t, N_fit,  0.75, 0, 2, nstep_v );
    // Print( "../data/Run/OptFitR075.dat", RunMins3, 12 );

    // vector<vector<double>> RunMins4 = FindMinFit( boxR, Run, 0, 2, N_vb, dx, nstep_t, N_fit,  1, 0, 2, nstep_v );
    // Print( "../data/Run/OptFitR1.dat", RunMins4, 12 );



    // v_opt mins brent
    // vector<vector<double>> WalkMins1 = FindMinBrent( boxW, Walk, 0, 0.7, dx, nstep_t, 0.001, 0., 0., 0.7, nstep_v );
    // Print( "../data/Walk/OptW0.dat", WalkMins1, 12 );

    // vector<vector<double>> WalkMins2 = FindMinBrent( boxW, Walk, 0, 0.7, dx, nstep_t, 0.001, 0.15, 0., 0.7, nstep_v );
    // Print( "../data/Walk/OptW015.dat", WalkMins2, 12 );

    // vector<vector<double>> WalkMins3 = FindMinBrent( boxW, Walk, 0, 0.7, dx, nstep_t, 0.001, 0.25, 0., 0.7, nstep_v );
    // Print( "../data/Walk/OptW025.dat", WalkMins3, 12 );

    // vector<vector<double>> WalkMins4 = FindMinBrent( boxW, Walk, 0, 0.7, dx, nstep_t, 0.001, 0.30, 0., 0.7, nstep_v );
    // Print( "../data/Walk/OptW030.dat", WalkMins4, 12 );

    // vector<vector<double>> RunMins1 = FindMinBrent( boxR, Run, 0, 2, dx, nstep_t, 0.001,  0., 0, 2, nstep_v );
    // Print( "../data/Run/OptR0.dat", RunMins1, 12 );

    // vector<vector<double>> RunMins2 = FindMinBrent( boxR, Run, 0, 2, dx, nstep_t, 0.001,  0.3, 0, 2, nstep_v );
    // Print( "../data/Run/OptR03.dat", RunMins2, 12 );

    // vector<vector<double>> RunMins3 = FindMinBrent( boxR, Run, 0, 2, dx, nstep_t, 0.001,  0.6, 0, 2, nstep_v );
    // Print( "../data/Run/OptR06.dat", RunMins3, 12 );

    // vector<vector<double>> RunMins4 = FindMinBrent( boxR, Run, 0, 2, dx, nstep_t, 0.001,  0.8, 0, 2, nstep_v );
    // Print( "../data/Run/OptR08.dat", RunMins4, 12 );

    // vector<vector<double>> RunMins5 = FindMinBrent( boxR, Run, 0, 2, dx, nstep_t, 0.001,  1, 0, 2, nstep_v );
    // Print( "../data/Run/OptR1.dat", RunMins5, 12 );



    // For fitting graph 
    // vector<vector<double>> RunFitGraph = Simulate( boxR,Run, rain_vel, 0, 2, N_vb, dx, 0, 1, nstep_t );
    // Print( "../data/Run/GraphFit.dat", RunFitGraph, 12 );


    // // Trial
    // auto wetfunc = [&boxR, &Run, &rain_vel, dx, nstep_t] (double x) {return Wetness( boxR, Run, rain_vel, x, dx, 0., 1., nstep_t );};
    // Brent trial(0.001);
    // if(trial.bracket( 0.01, 1.99, 2, wetfunc )) {
    //     double min = trial.minimize(wetfunc);
    //     cout << "Minumum found in vb = " << min << endl;
    // } else {
    //     cout << "Minimum not found." << endl;
    // }



    // Wetness around minimum
    // dx = 0.001;
    // nstep_t = 10;
    // double voptW = 0.59191;
    // vector<vector<double>> wetfunW = Simulate( boxW, Walk, vector<double>{0.5, 0.15, -1}, voptW - 0.03, voptW + 0.03, 600, dx, 0, 1, nstep_t );
    // Print( "../data/Walk/Wetfun.dat", wetfunW, 12 );
    
    // double voptR = 0.884849;
    // vector<vector<double>> wetfunR = Simulate( boxR, Run, vector<double>{0.45, 0.3, -1}, voptR - 0.015, voptR + 0.015, 300, dx, 0, 1, nstep_t );
    // Print( "../data/Run/Wetfun.dat", wetfunR, 12 );


    // dx = 0.0005;
    // nstep_t = 20;
    // double voptW = 0.59;
    // vector<vector<double>> wetfunW = Simulate( boxW, Walk, vector<double>{0.5, 0.15, -1}, voptW - 0.03, voptW + 0.03, 600, dx, 0, 1, nstep_t );
    // Print( "../data/Walk/Wetfun_dx_dt.dat", wetfunW, 12 );
    // double voptR = 0.885;
    // vector<vector<double>> wetfunR = Simulate( boxR, Run, vector<double>{0.45, 0.3, -1}, voptR - 0.015, voptR + 0.015, 300, dx, 0, 1, nstep_t );
    // Print( "../data/Run/Wetfun_dx_dt.dat", wetfunR, 12 );

    // dx = 0.0005;
    // nstep_t = 10;
    // double voptW = 0.59;
    // vector<vector<double>> wetfunW = Simulate( boxW, Walk, vector<double>{0.5, 0.15, -1}, voptW - 0.03, voptW + 0.03, 600, dx, 0, 1, nstep_t );
    // Print( "../data/Walk/Wetfun_dx.dat", wetfunW, 12 );
    // double voptR = 0.885;
    // vector<vector<double>> wetfunR = Simulate( boxR, Run, vector<double>{0.45, 0.3, -1}, voptR - 0.015, voptR + 0.015, 300, dx, 0, 1, nstep_t );
    // Print( "../data/Run/Wetfun_dx.dat", wetfunR, 12 );

    // dx = 0.001;
    // nstep_t = 20;
    // double voptW = 0.59;
    // vector<vector<double>> wetfunW = Simulate( boxW, Walk, vector<double>{0.5, 0.15, -1}, voptW - 0.03, voptW + 0.03, 600, dx, 0, 1, nstep_t );
    // Print( "../data/Walk/Wetfun_dt.dat", wetfunW, 12 );
    // double voptR = 0.885;
    // vector<vector<double>> wetfunR = Simulate( boxR, Run, vector<double>{0.45, 0.3, -1}, voptR - 0.015, voptR + 0.015, 300, dx, 0, 1, nstep_t );
    // Print( "../data/Run/Wetfun_dt.dat", wetfunR, 12 );


    // Trials
    // double wet = Wetness( boxR, Run, rain_vel, body_vel, dx, 0, 1, nstep_t );
    // cout << "Wentness: "<< setprecision(15) << wet << endl;

    // double wetSmooth = WetnessSmooth( boxR, Run, rain_vel, body_vel, dx, 0, 1, nstep_t );
    // cout << "Smooth wentness: "<< setprecision(15) << wetSmooth << endl;



    // Smooth wetness around minimum
    // dx = 0.001;
    // nstep_t = 10;
    // double voptW = 0.59;
    // vector<vector<double>> wetfunSW = SimulateSmooth( boxW, Walk, vector<double>{0.5, 0.15, -1}, voptW - 0.03, voptW + 0.03, 600, dx, 0, 1, nstep_t );
    // Print( "../data/Walk/WetfunS.dat", wetfunSW, 12 );
    // double voptR = 0.885;
    // vector<vector<double>> wetfunSR = SimulateSmooth( boxR, Run, vector<double>{0.45, 0.3, -1}, voptR - 0.015, voptR + 0.015, 300, dx, 0, 1, nstep_t );
    // Print( "../data/Run/WetfunS.dat", wetfunSR, 12 );

    // dx = 0.0005;
    // nstep_t = 20;
    // double voptW = 0.59;
    // vector<vector<double>> wetfunSW = SimulateSmooth( boxW, Walk, vector<double>{0.5, 0.15, -1}, voptW - 0.03, voptW + 0.03, 600, dx, 0, 1, nstep_t );
    // Print( "../data/Walk/WetfunS_dx_dt.dat", wetfunSW, 12 );
    // double voptR = 0.885;
    // vector<vector<double>> wetfunSR = SimulateSmooth( boxR, Run, vector<double>{0.45, 0.3, -1}, voptR - 0.015, voptR + 0.015, 300, dx, 0, 1, nstep_t );
    // Print( "../data/Run/WetfunS_dx_dt.dat", wetfunSR, 12 );

    // dx = 0.0005;
    // nstep_t = 10;
    // double voptW = 0.59;
    // vector<vector<double>> wetfunSW = SimulateSmooth( boxW, Walk, vector<double>{0.5, 0.15, -1}, voptW - 0.03, voptW + 0.03, 600, dx, 0, 1, nstep_t );
    // Print( "../data/Walk/WetfunS_dx.dat", wetfunSW, 12 );
    // double voptR = 0.885;
    // vector<vector<double>> wetfunSR = SimulateSmooth( boxR, Run, vector<double>{0.45, 0.3, -1}, voptR - 0.015, voptR + 0.015, 300, dx, 0, 1, nstep_t );
    // Print( "../data/Run/WetfunS_dx.dat", wetfunSR, 12 );

    // dx = 0.001;
    // nstep_t = 20;
    // double voptW = 0.59;
    // vector<vector<double>> wetfunSW = SimulateSmooth( boxW, Walk, vector<double>{0.5, 0.15, -1}, voptW - 0.03, voptW + 0.03, 600, dx, 0, 1, nstep_t );
    // Print( "../data/Walk/WetfunS_dt.dat", wetfunSW, 12 );
    // double voptR = 0.885;
    // vector<vector<double>> wetfunSR = SimulateSmooth( boxR, Run, vector<double>{0.45, 0.3, -1}, voptR - 0.015, voptR + 0.015, 300, dx, 0, 1, nstep_t );
    // Print( "../data/Run/WetfunS_dt.dat", wetfunSR, 12 );

    // Brent path 
    // Brent mins(0.001);
    // auto wetfuncW = [&boxW, &Walk, dx, nstep_t] (double x) {return WetnessSmooth( boxW, Walk, vector<double>{0.5, 0.15, -1}, x, dx, 0., 1., nstep_t );};
    // mins.bracket( 0., 0.7, 3, wetfuncW );
    // vector<vector<double>> pathW = mins.minimize( wetfuncW );
    // Print("../data/Walk/MinPath.dat", Transpose(pathW), 12 );
    // auto wetfuncR = [&boxR, &Run, dx, nstep_t] (double x) {return WetnessSmooth( boxR, Run, vector<double>{0.45, 0.3, -1}, x, dx, 0., 1., nstep_t );};
    // mins.bracket( 0., 2., 3, wetfuncR );
    // vector<vector<double>> pathR = mins.minimize( wetfuncR );
    // Print("../data/Run/MinPath.dat", Transpose(pathR), 12 );

    



    // Smooth v_opt mins brent
    // vector<vector<double>> WalkMins1 = FindMinBrentSmooth( boxW, Walk, 0., 0.7, dx, nstep_t, 0.001, 3, 0., 0., 0.7, nstep_v );
    // Print( "../data/Walk/OptSmoothW0_dx_dt.dat", WalkMins1, 12 );

    // vector<vector<double>> WalkMins2 = FindMinBrentSmooth( boxW, Walk, 0, 0.7, dx, nstep_t, 0.001, 3, 0.15, 0., 0.7, nstep_v );
    // Print( "../data/Walk/OptSmoothW015_dx_dt.dat", WalkMins2, 12 );

    // vector<vector<double>> WalkMins3 = FindMinBrentSmooth( boxW, Walk, 0, 0.7, dx, nstep_t, 0.001, 3, 0.25, 0., 0.7, nstep_v );
    // Print( "../data/Walk/OptSmoothW025_dx_dt.dat", WalkMins3, 12 );

    // vector<vector<double>> WalkMins4 = FindMinBrentSmooth( boxW, Walk, 0, 0.7, dx, nstep_t, 0.001, 3, 0.30, 0., 0.7, nstep_v );
    // Print( "../data/Walk/OptSmoothW030_dx_dt.dat", WalkMins4, 12 );

    // vector<vector<double>> RunMins1 = FindMinBrentSmooth( boxR, Run, 0, 2, dx, nstep_t, 0.001, 3, 0., 0, 2, nstep_v );
    // Print( "../data/Run/OptSmoothR0_dx_dt.dat", RunMins1, 12 );

    // vector<vector<double>> RunMins2 = FindMinBrentSmooth( boxR, Run, 0, 2, dx, nstep_t, 0.001, 3, 0.3, 0, 2, nstep_v );
    // Print( "../data/Run/OptSmoothR03_dx_dt.dat", RunMins2, 12 );

    // vector<vector<double>> RunMins3 = FindMinBrentSmooth( boxR, Run, 0, 2, dx, nstep_t, 0.001, 3, 0.6, 0, 2, nstep_v );
    // Print( "../data/Run/OptSmoothR06_dx_dt.dat", RunMins3, 12 );

    // vector<vector<double>> RunMins4 = FindMinBrentSmooth( boxR, Run, 0, 2, dx, nstep_t, 0.001, 3, 0.8, 0, 2, nstep_v );
    // Print( "../data/Run/OptSmoothR08_dx_dt.dat", RunMins4, 12 );

    // vector<vector<double>> RunMins5 = FindMinBrentSmooth( boxR, Run, 0, 2, dx, nstep_t, 0.001, 3, 1, 0, 2, nstep_v );
    // Print( "../data/Run/OptSmoothR1_dx_dt.dat", RunMins5, 12 );



    // Simulate different nstep
    // dx = 0.001;
    // int nstep_min = 2;
    // int nstep_max = 100;
    // int N_nstep = 10;

    // double voptW = 0.59;
    // vector<vector<double>> wetfunSW = SimulateNstepSmooth( boxW, Walk, vector<double>{0.5, 0.15, -1}, voptW - 0.05, voptW + 0.05, 40, dx, nstep_min, nstep_max, N_nstep );
    // Print( "../data/Walk/WetfunS_nstep.dat", wetfunSW, 12 );

    // double voptR = 0.89;
    // vector<vector<double>> wetfunSR = SimulateNstepSmooth( boxR, Run, vector<double>{0.45, 0.3, -1}, voptR - 0.015, voptR + 0.015, 40, dx, nstep_min, nstep_max, N_nstep );
    // Print( "../data/Run/WetfunS_nstep.dat", wetfunSR, 12 );

    // double voptDS = 2.872;
    // vector<vector<double>> wetfunDS = SimulateNstepSmooth( boxDS, DynSphere, vector<double>{0.45, 0.3, -1}, voptDS - 0.03, voptDS + 0.03, 40, dx, nstep_min, nstep_max, N_nstep );
    // Print( "../data/Sphere/WetfunS_nstep.dat", wetfunDS, 12 );

    // dx = 0.0005;
    // vector<vector<double>> wetfunDS = SimulateNstepSmooth( boxDS, DynSphere, vector<double>{0.45, 0.3, -1}, voptDS - 0.03, voptDS + 0.03, 40, dx, nstep_min, nstep_max, N_nstep );
    // Print( "../data/Sphere/WetfunS_nstep0005.dat", wetfunDS, 12 );



    // Check shadows of dynamic sphere
    // nstep_t = 20;
    // dx =  0.001;
    // PrintDynShadowSmooth( boxDS, DynSphere, vector<double>{0.45, 0.3, -1}, dx, 0, 1, nstep_t, "../data/Sphere/Shadows/Smooth001_" );
    // dx =  0.002;
    // PrintDynShadowSmooth( boxDS, DynSphere, vector<double>{0.45, 0.3, -1}, dx, 0, 1, nstep_t, "../data/Sphere/Shadows/Smooth002_" );
    


    // Smooth Brent minimization fit
    int N_vtail = 50;
    int N_fit = 9;
    double dv = 0.006;
    dx = 0.001;
    nstep_t = 50;

    // vector<vector<double>> WalkMinsS0 = FindMinFitSmooth( boxW, Walk, 0., 0.7, dx, nstep_t, 0., 0., 0.7, N_vtail, N_fit, dv );
    // Print( "../data/Walk/OptFitSmoothW0.dat", WalkMinsS0, 12 );

    // vector<vector<double>> WalkMinsS0_Extra = FindMinFitSmooth( boxW, Walk, 0., 0.7, dx, nstep_t, 0., 0.236, 0.25, 2, N_fit, dv );
    // Print( "../data/Walk/OptFitSmoothW0_Extra.dat", WalkMinsS0_Extra, 12 );

    // vector<vector<double>> WalkMinsS015 = FindMinFitSmooth( boxW, Walk, 0., 0.7, dx, nstep_t, 0.15, 0., 0.7, N_vtail, N_fit, dv );
    // Print( "../data/Walk/OptFitSmoothW015.dat", WalkMinsS015, 12 );

    // vector<vector<double>> WalkMinsS015_Extra = FindMinFitSmooth( boxW, Walk, 0., 0.7, dx, nstep_t, 0.15, 0.250, 0.250, 1, N_fit, dv );
    // Print( "../data/Walk/OptFitSmoothW015_Extra.dat", WalkMinsS015_Extra, 12 );

    // vector<vector<double>> WalkMinsS025 = FindMinFitSmooth( boxW, Walk, 0., 0.7, dx, nstep_t, 0.25, 0., 0.7, N_vtail, N_fit, dv );
    // Print( "../data/Walk/OptFitSmoothW025.dat", WalkMinsS025, 12 );

    // vector<vector<double>> WalkMinsS030 = FindMinFitSmooth( boxW, Walk, 0., 0.7, dx, nstep_t, 0.30, 0., 0.7, N_vtail, N_fit, dv );
    // Print( "../data/Walk/OptFitSmoothW030.dat", WalkMinsS030, 12 );

    vector<vector<double>> WalkMinsS030_Extra = FindMinFitSmooth( boxW, Walk, 0., 0.7, dx, nstep_t, 0.30, 0.3785714285, 0.3785714285, 1, N_fit, dv );
    Print( "../data/Walk/OptFitSmoothW030_Extra.dat", WalkMinsS030_Extra, 12 );

    // vector<vector<double>> RunMinsS0 = FindMinFitSmooth( boxR, Run, 0., 2., dx, nstep_t, 0., 0., 2, N_vtail, N_fit, dv );
    // Print( "../data/Run/OptFitSmoothR0.dat", RunMinsS0, 12 );

    // vector<vector<double>> RunMinsS0_Extra = FindMinFitSmooth( boxR, Run, 0., 2., dx, nstep_t, 0., 0.143, 0.143, 1, N_fit, dv );
    // Print( "../data/Run/OptFitSmoothR0_Extra.dat", RunMinsS0_Extra, 12 );

    // vector<vector<double>> RunMinsS030 = FindMinFitSmooth( boxR, Run, 0., 2., dx, nstep_t, 0.3, 0., 2, N_vtail, N_fit, dv );
    // Print( "../data/Run/OptFitSmoothR030.dat", RunMinsS030, 12 );

    // vector<vector<double>> RunMinsS030_Extra = FindMinFitSmooth( boxR, Run, 0., 2., dx, nstep_t, 0.3, 0.183, 0.183, 1, N_fit, dv );
    // Print( "../data/Run/OptFitSmoothR030_Extra.dat", RunMinsS030_Extra, 12 );

    // vector<vector<double>> RunMinsS060 = FindMinFitSmooth( boxR, Run, 0., 2., dx, nstep_t, 0.6, 0., 2, N_vtail, N_fit, dv );
    // Print( "../data/Run/OptFitSmoothR060.dat", RunMinsS060, 12 );

    // vector<vector<double>> RunMinsS080 = FindMinFitSmooth( boxR, Run, 0., 2., dx, nstep_t, 0.8, 0., 2, N_vtail, N_fit, dv );
    // Print( "../data/Run/OptFitSmoothR080.dat", RunMinsS080, 12 );

    // vector<vector<double>> RunMinsS1 = FindMinFitSmooth( boxR, Run, 0., 2., dx, nstep_t, 1., 0., 2, N_vtail, N_fit, dv );
    // Print( "../data/Run/OptFitSmoothR1.dat", RunMinsS1, 12 );
    

    // vector<vector<double>> WalkMinsSEx= FindMinFitSmooth( boxW, Walk, 0., 0.7, dx, nstep_t, 0.15, 0.5, 0.7, 1, N_fit, dv );
    // Print( "../data/Walk/OptFitSmoothW_Ex.dat", WalkMinsSEx, 12 );

    // vector<vector<double>> RunMinsSEx = FindMinFitSmooth( boxR, Run, 0., 2., dx, nstep_t, 0.3, 0.45, 2, 1, N_fit, dv );
    // Print( "../data/Run/OptFitSmoothR_Ex.dat", RunMinsSEx, 12 );


    // Smooth Brent minimization fit with variable nstep
    // int nstep_min = 5;
    // int nstep_max = 100;
    // int N_nstep = 30;
    // int N_fit = 9;
    // double dv = 0.006;
    // dx = 0.001;


    // vector<vector<double>> WalkMinsSnstep= FindMinFitSmoothNstep( boxW, Walk, 0., 0.7, dx, nstep_min, nstep_max, N_nstep, 0.15, 0.5, N_fit, dv );
    // Print( "../data/Walk/OptFitSmoothW_nstep.dat", WalkMinsSnstep, 12 );

    // vector<vector<double>> RunMinsSnstep = FindMinFitSmoothNstep( boxR, Run, 0., 2., dx, nstep_min, nstep_max, N_nstep, 0.30, 0.45, N_fit, dv );
    // Print( "../data/Run/OptFitSmoothR_nstep.dat", RunMinsSnstep, 12 );


    // Smooth Brent minimization fit map
    // nstep_v = 20;
    // int N_fit = 9;
    // double dv = 0.006;
    // dx = 0.001;
    // nstep_t = 50;

    // vector<vector<double>> WalkOptMapFit = OptMapFitSmooth( boxW, Walk, 0, 0.7, dx, nstep_t, N_fit, dv, 0, 0.7, nstep_v, 0, 0.40, nstep_v);
    // Print( "../data/Walk/OptMapFitSmoothW.dat", WalkOptMapFit, 12 );

    // vector<vector<double>> RunOptMapFit = OptMapFitSmooth( boxR, Run, 0, 2, dx, nstep_t, N_fit, dv, 0, 2, nstep_v, 0, 1.2, nstep_v );
    // Print( "../data/Run/OptMapFitSmoothR.dat", RunOptMapFit, 12 ); 

    // vector<vector<double>> WalkOptMapComp = OptMapFitSmooth( boxW, Walk, 0, 0.7, dx, nstep_t, N_fit, dv, 0., 1.2, nstep_v, 0, 1.2, nstep_v );
    // Print( "../data/Walk/OptMapCompareSmoothW.dat", WalkOptMapComp, 12 );

    // vector<vector<double>> RunOptMapComp = OptMapFitSmooth( boxR, Run, 0, 2, dx, nstep_t, N_fit, dv, 0., 1.2, nstep_v, 0, 1.2, nstep_v );
    // Print( "../data/Run/OptMapCompareSmoothR.dat", RunOptMapComp, 12 );


    // vector<vector<double>> WalkOptMapComp = OptMapFitSmooth( boxW, Walk, 0, 0.7, dx, nstep_t, N_fit, dv, 0., 1.2, nstep_v, 0, 1.2*(4./19), 5 );
    // Print( "../data/Walk/OptMapCompareSmoothW1.dat", WalkOptMapComp, 12 );
    // vector<vector<double>> WalkOptMapComp = OptMapFitSmooth( boxW, Walk, 0, 0.7, dx, nstep_t, N_fit, dv, 0., 1.2, nstep_v, 1.2*(5./19), 1.2*(9./19), 5 );
    // Print( "../data/Walk/OptMapCompareSmoothW2.dat", WalkOptMapComp, 12 );
    // vector<vector<double>> WalkOptMapComp = OptMapFitSmooth( boxW, Walk, 0, 0.7, dx, nstep_t, N_fit, dv, 0., 1.2, nstep_v, 1.2*(10./19), 1.2*(14./19), 5 );
    // Print( "../data/Walk/OptMapCompareSmoothW3.dat", WalkOptMapComp, 12 );
    // vector<vector<double>> WalkOptMapComp = OptMapFitSmooth( boxW, Walk, 0, 0.7, dx, nstep_t, N_fit, dv, 0., 1.2, nstep_v, 1.2*(15./19), 1.2*(19./19), 5 );
    // Print( "../data/Walk/OptMapCompareSmoothW4.dat", WalkOptMapComp, 12 );

    // vector<vector<double>> RunOptMapComp = OptMapFitSmooth( boxR, Run, 0, 2, dx, nstep_t, N_fit, dv, 0., 1.2, nstep_v, 0, 1.2*(1./19), 2 );
    // Print( "../data/Run/OptMapCompareSmoothR1.dat", RunOptMapComp, 12 ); 
    // vector<vector<double>> RunOptMapComp = OptMapFitSmooth( boxR, Run, 0, 2, dx, nstep_t, N_fit, dv, 0., 1.2, nstep_v, 1.2*(2./19), 1.2*(3./19), 2 );
    // Print( "../data/Run/OptMapCompareSmoothR2.dat", RunOptMapComp, 12 ); 
    // vector<vector<double>> RunOptMapComp = OptMapFitSmooth( boxR, Run, 0, 2, dx, nstep_t, N_fit, dv, 0., 1.2, nstep_v, 1.2*(4./19), 1.2*(5./19), 2 );
    // Print( "../data/Run/OptMapCompareSmoothR3.dat", RunOptMapComp, 12 ); 
    // vector<vector<double>> RunOptMapComp = OptMapFitSmooth( boxR, Run, 0, 2, dx, nstep_t, N_fit, dv, 0., 1.2, nstep_v, 1.2*(6./19), 1.2*(7./19), 2 );
    // Print( "../data/Run/OptMapCompareSmoothR4.dat", RunOptMapComp, 12 ); 
    // vector<vector<double>> RunOptMapComp = OptMapFitSmooth( boxR, Run, 0, 2, dx, nstep_t, N_fit, dv, 0., 1.2, nstep_v, 1.2*(8./19), 1.2*(9./19), 2 );
    // Print( "../data/Run/OptMapCompareSmoothR5.dat", RunOptMapComp, 12 );
    // vector<vector<double>> RunOptMapComp = OptMapFitSmooth( boxR, Run, 0, 2, dx, nstep_t, N_fit, dv, 0., 1.2, nstep_v, 1.2*(10./19), 1.2*(11./19), 2 );
    // Print( "../data/Run/OptMapCompareSmoothR6.dat", RunOptMapComp, 12 );
    // vector<vector<double>> RunOptMapComp = OptMapFitSmooth( boxR, Run, 0, 2, dx, nstep_t, N_fit, dv, 0., 1.2, nstep_v, 1.2*(12./19), 1.2*(13./19), 2 );
    // Print( "../data/Run/OptMapCompareSmoothR7.dat", RunOptMapComp, 12 );
    // vector<vector<double>> RunOptMapComp = OptMapFitSmooth( boxR, Run, 0, 2, dx, nstep_t, N_fit, dv, 0., 1.2, nstep_v, 1.2*(14./19), 1.2*(15./19), 2 );
    // Print( "../data/Run/OptMapCompareSmoothR8.dat", RunOptMapComp, 12 );
    // vector<vector<double>> RunOptMapComp = OptMapFitSmooth( boxR, Run, 0, 2, dx, nstep_t, N_fit, dv, 0., 1.2, nstep_v, 1.2*(16./19), 1.2*(17./19), 2 );
    // Print( "../data/Run/OptMapCompareSmoothR9.dat", RunOptMapComp, 12 );
    // vector<vector<double>> RunOptMapComp = OptMapFitSmooth( boxR, Run, 0, 2, dx, nstep_t, N_fit, dv, 0., 1.2, nstep_v, 1.2*(18./19), 1.2*(19./19), 2 );
    // Print( "../data/Run/OptMapCompareSmoothR10.dat", RunOptMapComp, 12 );
    


    // Points for example fit
    // dx = 0.001;
    // nstep_t = 50;
    // nstep_v = 300;
    // rain_vel = { 0.5, 0.15, -1};
    // double vmin = 0.55, vmax = 0.64;

    // vector<vector<double>> EgWetfun = SimulateSmooth( boxW, Walk, rain_vel, vmin, vmax, nstep_v, dx, 0, 1, nstep_t );
    // Print( "../data/Walk/EgWetFun.dat", EgWetfun , 12 ); 
}