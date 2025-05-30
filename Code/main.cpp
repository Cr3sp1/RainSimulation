#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <ostream>
#include <random>
#include <string>
#include <vector>

#include "VectorStat.h"
#include "body.h"
#include "mins.h"
#include "ray.h"

using namespace std;

int main(int argc, char* argv[]) {

	// Check that everything is working properly
	if (AllGood() == false)
		return 1;

	// Declare stuff
	string bodyPath, bodyName;
	double vmax;
	double dx;
	int nstep;
	int nfit;
	double dv;
	double vtail_min, vtail_max;
	int n_vtail;
	double vcross_min, vcross_max;
	int n_vcross;
	string resPath;
	string fitPath;

	// Read the input file
	ifstream ReadInput("input.in");
	if (!ReadInput) {
		std::cerr << "Error opening input file." << std::endl;
		return 1;
	}
	ReadInput >> bodyPath;
	bodyName = filesystem::path(bodyPath).filename().string();
	ReadInput >> vmax;
	ReadInput >> dx;
	ReadInput >> nstep;
	ReadInput >> nfit;
	ReadInput >> dv;
	ReadInput >> vtail_min >> vtail_max >> n_vtail;
	ReadInput >> vcross_min >> vcross_max >> n_vcross;
	ReadInput >> resPath;
	ReadInput >> fitPath;
	ReadInput.close();

	// Output parameters
	cout << "######################################################################################"
			"############\n";
	cout << "Body = " + bodyName + ", vmax = " << vmax << " vfall\n";
	cout << "dx = " << dx << " m, nstep = " << nstep << "\n";
	cout << "nfit = " << nfit << ", dv = " << dv << " vfall\n";
	cout << "vtail_min = " << vtail_min << " vfall, vtail_max = " << vtail_max
		 << " vfall, n_vtail = " << n_vtail << "\n";
	cout << "vcross_min = " << vcross_min << " vfall, vcross_max = " << vcross_max
		 << " vfall, n_vcross = " << n_vcross << "\n";
	cout << "Printing fit results to " + resPath + " and fit points to " + fitPath + "\n";
	cout << "######################################################################################"
			"############"
		 << endl;

	// Get Body and box
	ManyBody body(bodyPath);
	vector<double> box = body.GetBox(0, 1, nstep, dx);

	// Prepare output streams and files
	ofstream outRes(resPath);
	if (!outRes) {
		cerr << "Error opening results output file." << endl;
		return 1;
	}
	WriteHeadRes(outRes, bodyName, vmax, dx, nstep, nfit, dv);
	outRes << fixed << setprecision(12);

	ostream* outFit;
	ofstream outFitFile;
	ostringstream nullStream;
	if (fitPath == "None") {
		outFit = &nullStream;
		cout << "Not writing fit points to file." << endl;
	} else {
		outFitFile.open(fitPath);
		if (!outFitFile) {
			cerr << "Error opening results output file." << endl;
			return 1;
		}
		outFit = &outFitFile;
	}
	WriteHeadFit(*outFit, bodyName, vmax, dx, nstep, nfit, dv);
	*outFit << fixed << setprecision(12);

	// Evaluate and write to file results
	for (size_t i = 0; i < n_vcross; i++) {
		double vcross_i =
			n_vcross > 1 ? vcross_min + i * (vcross_max - vcross_min) / (n_vcross - 1) : vcross_min;

		for (size_t j = 0; j < n_vtail; j++) {
			double vtail_j =
				n_vtail > 1 ? vtail_min + j * (vtail_max - vtail_min) / (n_vtail - 1) : vtail_min;

			double vopt, vopt_std, Rmin, Rmin_std;
			vector<vector<double>> fitPoints;

			tie(vopt, vopt_std, Rmin, Rmin_std, fitPoints) =
				MinFit(box, body, vmax, dx, nstep, vcross_i, vtail_j, nfit, dv);

			outRes << vcross_i << "\t" << vtail_j << "\t" << vopt << "\t" << vopt_std << "\t"
				   << Rmin << "\t" << Rmin_std << endl;

			// Add vcross_i and vtail_j and to each row and write to file
			for (size_t i = 0; i < fitPoints.size(); ++i) {
				fitPoints[i].insert(fitPoints[i].begin(), vtail_j);
				fitPoints[i].insert(fitPoints[i].begin(), vcross_i);
			}
			Print(*outFit, fitPoints, 12);

			cout << "Step " << i * n_vtail + j + 1 << "/" << n_vtail * n_vcross << " completed!"
				 << endl;
		}
	}



	cout << "All done!" << endl;

	return 0;
}