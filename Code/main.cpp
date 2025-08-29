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
using namespace tinyxml2;

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
	XMLDocument doc;
	if (doc.LoadFile("Input.xml") != XML_SUCCESS) {
		cerr << "Error loading Input.xml file." << endl;
		return 0;
	}
	XMLElement* root = doc.FirstChildElement("Settings");
	if (!root) {
		cerr << "Error: No <Settings> root found in Input.xml." << endl;
		return 0;
	}
	XMLElement* elem = nullptr;
	elem = root->FirstChildElement("BodyPath");
	if (elem && elem->GetText())
		bodyPath = elem->GetText();
	else {
		cout << "Error: Missing or empty <BodyPath>." << endl;
		return 0;
	}
	bodyName = filesystem::path(bodyPath).filename().string();
	elem = root->FirstChildElement("Vmax");
	if (!elem || elem->QueryDoubleText(&vmax) != XML_SUCCESS) {
		cout << "Error: Missing or invalid <Vmax>." << endl;
		return 0;
	}
	elem = root->FirstChildElement("Dx");
	if (!elem || elem->QueryDoubleText(&dx) != XML_SUCCESS) {
		cout << "Error: Missing or invalid <Dx>." << endl;
		return 0;
	}
	elem = root->FirstChildElement("Nstep");
	if (!elem || elem->QueryIntText(&nstep) != XML_SUCCESS) {
		cout << "Error: Missing or invalid <Nstep>." << endl;
		return 0;
	}
	elem = root->FirstChildElement("Nfit");
	if (!elem || elem->QueryIntText(&nfit) != XML_SUCCESS) {
		cout << "Error: Missing or invalid <Nfit>." << endl;
		return 0;
	}
	elem = root->FirstChildElement("Dv");
	if (!elem || elem->QueryDoubleText(&dv) != XML_SUCCESS) {
		cout << "Error: Missing or invalid <Dv>." << endl;
		return 0;
	}
	elem = root->FirstChildElement("Vcross");
	if (!elem || elem->QueryDoubleAttribute("min", &vcross_min) != XML_SUCCESS ||
		elem->QueryDoubleAttribute("max", &vcross_max) != XML_SUCCESS ||
		elem->QueryIntAttribute("n", &n_vcross) != XML_SUCCESS) {
		cerr << "Error: Missing or invalid <Vcross> (with attributes min, max, n)." << endl;
		return 1;
	}
	elem = root->FirstChildElement("Vtail");
	if (!elem || elem->QueryDoubleAttribute("min", &vtail_min) != XML_SUCCESS ||
		elem->QueryDoubleAttribute("max", &vtail_max) != XML_SUCCESS ||
		elem->QueryIntAttribute("n", &n_vtail) != XML_SUCCESS) {
		cerr << "Error: Missing or invalid <Vtail> (with attributes min, max, n)." << endl;
		return 1;
	}
	elem = root->FirstChildElement("ResPath");
	if (elem && elem->GetText())
		resPath = elem->GetText();
	else
		resPath = "";
	elem = root->FirstChildElement("FitPath");
	if (elem && elem->GetText())
		fitPath = elem->GetText();
	else
		fitPath = "";

	// Output parameters
	cout << "######################################################################################"
			"############\n";
	cout << "Body = " + bodyName + ", vmax = " << vmax << " vfall\n";
	cout << "dx = " << dx << " m, nstep = " << nstep << "\n";
	cout << "nfit = " << nfit << ", dv = " << dv << " vfall\n";
	cout << "vcross_min = " << vcross_min << " vfall, vcross_max = " << vcross_max
		 << " vfall, n_vcross = " << n_vcross << "\n";
	cout << "vtail_min = " << vtail_min << " vfall, vtail_max = " << vtail_max
		 << " vfall, n_vtail = " << n_vtail << "\n";
	if (resPath != "") {
		cout << "Printing fit results to " << resPath << "\n";
	} else {
		cout << "Not printing fit results\n";
	}
	if (fitPath != "") {
		cout << "Printing fit points to " << fitPath << "\n";
	} else {
		cout << "Not printing fit points\n";
	}
	cout << "######################################################################################"
			"############"
		 << endl;

	// Get Body and box
	ManyBody* body;
	try {
		body = new ManyBody(bodyPath);
	} catch (logic_error e) {
		cout << "Error reading the body input file:" << endl;
		cout << e.what() << endl;
		return 0;
	}

	vector<double> box = body->GetBox(0, 1, nstep, dx);

	// Prepare output streams and files
	ofstream outRes;
	if (resPath != "") {
		filesystem::path resDir = filesystem::path(resPath).parent_path();
        if (!resDir.empty() && !filesystem::exists(resDir)) {
            filesystem::create_directories(resDir);
        }
		outRes.open(resPath);
		if (!outRes) {
			cerr << "Error opening results output file." << endl;
			return 0;
		}
		WriteHeadRes(outRes, bodyName, vmax, dx, nstep, nfit, dv);
		outRes << fixed << setprecision(12);
	}

	ofstream outFit;
	if (fitPath != "") {
		filesystem::path fitDir = filesystem::path(fitPath).parent_path();
        if (!fitDir.empty() && !filesystem::exists(fitDir)) {
            filesystem::create_directories(fitDir);
        }
		outFit.open(fitPath);
		if (!outFit) {
			cerr << "Error opening fit output file." << endl;
			return 0;
		}
		WriteHeadFit(outFit, bodyName, vmax, dx, nstep, nfit, dv);
		outFit << fixed << setprecision(12);
	}

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
				MinFit(box, *body, vmax, dx, nstep, vcross_i, vtail_j, nfit, dv);

			if (resPath != "")
				outRes << vcross_i << "\t" << vtail_j << "\t" << vopt << "\t" << vopt_std << "\t"
					   << Rmin << "\t" << Rmin_std << endl;

			// Add vcross_i and vtail_j and to each row and write to file
			for (size_t i = 0; i < fitPoints.size(); ++i) {
				fitPoints[i].insert(fitPoints[i].begin(), vtail_j);
				fitPoints[i].insert(fitPoints[i].begin(), vcross_i);
			}
			if (fitPath != "")
				Print(outFit, fitPoints, 12);

			cout << "Step " << i * n_vtail + j + 1 << "/" << n_vtail * n_vcross << " completed!"
				 << endl;
		}
	}

	cout << "All done!" << endl;

	return 0;
}