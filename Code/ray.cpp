#include "ray.h"

#include "body.h"

using namespace std;

// Declares V
vector<double> Ray::V(3);

// Resets the surface between steps (turns on all rays)
void ProjSurface::reset() {
	for (size_t i = 0; i < rays.size(); i++) {
		rays[i].On();
		rays[i].SetWeight(0.);
	}
}

// Complete constructor
ProjSurface::ProjSurface(vector<double> box, vector<double> vel, double Dx) : dx(Dx) {
	// Checks validity of arguments
	if (box.size() != 3 or box[0] <= 0 or box[1] <= 0 or box[2] <= 0) {
		throw invalid_argument("Box must be of size 3 with all components positive!");
	}
	if (vel.size() != 3) {
		throw invalid_argument("Velocity must be of size 3!");
	}

	// Finds vertex between three "seen" faces and projects the other vertices on
	// the plane passing by it perpendicular to vel
	vector<vector<double>> sides = {{box[0], 0, 0}, {0, box[1], 0}, {0, 0, box[2]}};
	H = FindHexProj({0, 0, 0}, sides, vel, FindMiddle({0, 0, 0}, sides, vel));

	// Generates the rays on a square grid along directions u1 and u2 (u1 and ud perpendicular and
	// belonging to surface) cout << "Generating rays" << endl;
	rays = {};
	// cout << "dx = " << dx << endl;
	vector<double> u1 = (Norm(H[5] - H[0]) > Norm(H[3] - H[0]))
							? H[5] - H[0]
							: H[3] - H[0]; // makes sure that u1 isn't infinitesimal
	u1 = (u1 / Norm(u1)) * dx;
	double maxu1 = MaxU(H, u1);
	// cout << "|u1| = " << Norm(u1) << endl;
	// cout << "steps1 = " << Norm(u1) << endl;
	vector<double> u2 = CrossProduct(u1, vel);
	u2 = (u2 / Norm(u2)) * dx;
	double maxu2 = MaxU(H, u2);
	// cout << "|u2| = " << Norm(u2) << endl;

	for (double sign1 = -1; sign1 < 2; sign1 += 2) {
		vector<double> point1 = (sign1 == -1) ? H[0] : H[0] + u1;

		while (Norm(point1 - H[0]) < maxu1) {
			for (double sign2 = -1; sign2 < 2; sign2 += 2) {
				vector<double> point2 = (sign2 == -1) ? point1 : point1 + u2;

				while (Norm(point2 - point1) < maxu2) {
					if (PointIsInsideT(point2, H)) {
						Ray temp(point2);
						rays.push_back(temp);
					}
					point2 += (sign2 * u2);
				}
			}
			point1 += (sign1 * u1);
		}
	}

	// Sets the rain speed
	Ray::V = vel;

	cout << "Number of rays generated: " << rays.size() << endl;
}

// Prints all the origins of the active rays to file
void ProjSurface::PrintR(ofstream& fout) {
	for (Ray& ray : rays) {
		vector<double> r = ray.GetR0();
		fout << r[0] << ", " << r[1] << ", " << r[2] << ", " << ray.IsOn() << endl;
	}
}

void ProjSurface::PrintR(string outfile) {
	ofstream fout(outfile);
	PrintR(fout);
	fout.close();
}

// Prints all the origins of the rays projected on the x-y plane to file, last column is 1 if
// active, 0 if not
void ProjSurface::PrintRaysFlat(ofstream& fout) {
	// Finds rotation to make points parallel to x-y
	vector<double> uz = {0, 0, 1};
	vector<double> v = Ray::V / Norm(Ray::V);
	vector<double> axis = CrossProduct(v, uz);
	double theta = acos(v * uz);
	vector<vector<double>> rotmat = RotMat(axis, theta);

	// Finds translation that brings H[0] to 0 after rotation
	vector<double> trans = rotmat * H[0];

	for (Ray& ray : rays) {
		vector<double> r = ray.GetR0();
		r = rotmat * r;
		r -= trans;
		fout << r[0] << "," << r[1] << "," << ray.IsOn() << endl;
	}
}

void ProjSurface::PrintRaysFlat(string outfile) {
	ofstream fout(outfile);
	PrintRaysFlat(fout);
	fout.close();
}

// Prints all the origins of the active rays projected on the x-y plane and their weight to file
void ProjSurface::PrintRaysFlatSmooth(ofstream& fout) {
	// Finds rotation to make points parallel to x-y
	vector<double> uz = {0, 0, 1};
	vector<double> v = Ray::V / Norm(Ray::V);
	vector<double> axis = CrossProduct(v, uz);
	double theta = acos(v * uz);
	vector<vector<double>> rotmat = RotMat(axis, theta);

	// Finds translation that brings H[0] to 0 after rotation
	vector<double> trans = rotmat * H[0];

	for (Ray& ray : rays) {
		vector<double> r = ray.GetR0();
		r = rotmat * r;
		r -= trans;
		fout << r[0] << "," << r[1] << "," << ray.GetWeight() << endl;
	}
}

void ProjSurface::PrintRaysFlatSmooth(string outfile) {
	ofstream fout(outfile);
	PrintRaysFlatSmooth(fout);
	fout.close();
}

// Prints H to file
void ProjSurface::PrintH(ofstream& fout) {
	for (size_t i = 0; i < H.size(); i++) {
		fout << H[i][0] << "," << H[i][1] << "," << H[i][2] << endl;
	}
}

void ProjSurface::PrintH(string outfile) {
	ofstream fout(outfile);
	PrintH(fout);
	fout.close();
}

// Returns an estimate of the projection of the body on the plane
double ProjSurface::BodyProj(Body& body) {
	unsigned int nhit = 0;
	// cout << "Projecting on " << rays.size() << " rays" << endl;
	body.Prime(H[0], Ray::V);
	for (Ray& ray : rays) {
		if (body.Check(ray)) {
			ray.Off();
			nhit++;
		}
	}
	// cout << "nhit = " << nhit << endl;
	return nhit * dx * dx;
}

// Returns an estimate of the projection of the dynamic body on the plane
double ProjSurface::BodyProj(Body& body, double tmin, double tmax, unsigned int nstep) {
	if (nstep == 0)
		return 0;
	double dt = (tmax - tmin) / nstep;
	unsigned int nhit = 0;

	// cout << "Projecting on " << rays.size() << " rays" << endl;
	for (unsigned int i = 0; i < nstep; i++) {
		double t = tmin + i * dt;
		body.Move(t);
		body.Prime(H[0], Ray::V);
		for (Ray& ray : rays) {
			if (body.Check(ray))
				nhit++;
		}
	}
	// cout << "nhit = " << nhit << endl;
	return nhit * dx * dx / nstep;
}

// Returns a smooth estimate of the projection of the body on the plane
double ProjSurface::BodyProjSmooth(Body& body) {
	double tot_w = 0;
	// cout << "Projecting on " << rays.size() << " rays" << endl;
	body.Prime(H[0], Ray::V);
	for (Ray& ray : rays) {
		double w = body.CheckSmooth(ray, dx);
		tot_w += w;
		ray.SetWeight(w);
	}
	// cout << "nhit = " << nhit << endl;
	return tot_w * dx * dx;
}

// Returns a smooth estimate of the projection of the dynamic body on the plane
double ProjSurface::BodyProjSmooth(Body& body, double tmin, double tmax, unsigned int nstep) {
	if (nstep == 0)
		return 0;
	double dt = (tmax - tmin) / nstep;
	double tot_w = 0;

	// cout << "Projecting on " << rays.size() << " rays" << endl;
	for (unsigned int i = 0; i < nstep; i++) {
		double t = tmin + i * dt;
		body.Move(t);
		body.Prime(H[0], Ray::V);
		for (Ray& ray : rays) {
			tot_w += body.CheckSmooth(ray, dx);
		}
	}
	// cout << "nhit = " << nhit << endl;
	return tot_w * dx * dx / nstep;
}