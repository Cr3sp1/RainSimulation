#include "body.h"

#include "RainFunctions.h"
#include "ray.h"

using namespace std;

// Virtual destructor
Body::~Body() {
	// cout << "Destroying " << name << endl;
}

// Evaluate theta at a time t
double Body::getTheta(double t) {
	double res = 0;
	// Add sin terms
	for (size_t i = 0; 2 * i < w.size(); i++) {
		res += w[2 * i] * sin(t * (i + 1) * 2 * M_PI);
	}
	// Add cos terms
	for (size_t i = 0; 2 * i + 1 < w.size(); i++) {
		res += w[2 * i + 1] * cos(t * (i + 1) * 2 * M_PI);
	}
	return res;
}

// Evaluate delta at a time t
vector<double> Body::getDelta(double t) {
	vector<double> res = {0, 0, 0};
	// Add sin terms
	for (size_t i = 0; 2 * i < trans.size(); i++) {
		res += trans[2 * i] * sin(t * (i + 1) * 2 * M_PI);
	}
	// Add cos terms
	for (size_t i = 0; 2 * i + 1 < trans.size(); i++) {
		res += trans[2 * i + 1] * cos(t * (i + 1) * 2 * M_PI);
	}
	return res;
}

// Time evolution of the body in its own frame of reference, also propagates to the sub-bodies
void Body::Move(double tNew) {
	// Calculate the total translation and rotation for the step
	vector<double> deltaNew = getDelta(tNew);
	vector<double> deltaShift = deltaNew - delta;
	delta = deltaNew;

	double thetaNew = getTheta(tNew);
	double thetaShift = thetaNew - theta;
	theta = thetaNew;
	vector<vector<double>> rotmat = RotMat(rotax, thetaShift);

	// Move
	Translate(deltaShift);
	Rotate(rotcent, rotmat);

	// Propagate motion to sub-bodies
	for (Body* body : SubBodies)
		body->BeMoved(deltaShift, rotcent, rotmat);
}

// Time evolution caused by the super-body, affects the whole frame of reference, also propagates to
// the sub-bodie
void Body::BeMoved(vector<double> shift, vector<double> rot0, vector<vector<double>> rotmat) {
	// Translate
	Translate(shift);

	// Rotate
	Rotate(rot0, rotmat);
	// Rotate frame of reference
	if (w.size() != 0) {
		RotatePoint(rotcent, rot0, rotmat);
		rotax = rotmat * rotax;
	}
	for (vector<double>& vec : trans)
		vec = rotmat * vec;

	// Move sub-bodies
	for (Body* body : SubBodies)
		body->BeMoved(shift, rot0, rotmat);
}

// Find smallest bounds containing the body throughout movement, return a vector containing lower
// values of x, y, z and a vector containing higher values of x, y, z.
tuple<vector<double>, vector<double>> Body::GetBounds(double tmin, double tmax,
													  unsigned int nstep) {
	vector<double> low = vector<double>(3, numeric_limits<double>::infinity());
	vector<double> high = vector<double>(3, -numeric_limits<double>::infinity());
	return make_tuple(low, high);
}

// Return box starting on origin that contains the body with allowance epsilon and move body in it
vector<double> Body::GetBox(double tmin, double tmax, unsigned int nstep, double epsilon) {
	// Get bounds
	vector<double> low, high;
	tie(low, high) = GetBounds(tmin, tmax, nstep);
	low -= vector<double>(3, epsilon);
	high += vector<double>(3, epsilon);

	// Translate Body
	Translate(-low);

	return high - low;
}

// Prints to file the state of the body
void Body::PrintState(ofstream& fout) {}

void Body::PrintState(string outfile) {
	ofstream fout(outfile);
	PrintState(fout);
	fout.close();
}

// Primes the body to be checked (projects the center of the sphere onto the surface)
void Sphere::Prime(vector<double> p, vector<double> v) { Hcent = Project(cent, p, v); }

// Return a value in [0, 1] describing how close the ray is to the body, 0 if the ray is at least a
// distance dx from the body, 1 if the ray is at least dx inside the body
double Sphere::Check(Ray& ray, double dx) {
	vector<double> Xrel = ray.GetR0() - Hcent;
	double delta_r = Norm(Xrel) - rad;
	return d_to_w(delta_r, dx);
}

// Translates the sphere by shift
void Sphere::Translate(vector<double> shift) {
	if (w.size() > 0)
		rotcent += shift;
	cent += shift;
}

// Rotates the sphere around point rot0 according to rotation matrix rotmat
void Sphere::Rotate(vector<double> rot0, vector<vector<double>> rotmat) {
	RotatePoint(cent, rot0, rotmat);
}

// Find smallest bounds containing the sphere throughout movement, return a vector containing lower
// values of x, y, z and a vector containing higher values of x, y, z.
tuple<vector<double>, vector<double>> Sphere::GetBounds(double tmin, double tmax,
														unsigned int nstep) {
	vector<double> low = vector<double>(3, numeric_limits<double>::infinity());
	vector<double> high = vector<double>(3, -numeric_limits<double>::infinity());
	if (nstep == 0)
		return make_tuple(low, high);
	double dt = (tmax - tmin) / nstep;

	for (unsigned int i = 0; i < nstep; i++) {
		double t = tmin + i * dt;
		Move(t);
		// Update bounds
		for (size_t j = 0; j < 3; j++) {
			low[j] = min(low[j], cent[j] - rad);
			high[j] = max(high[j], cent[j] - rad);
		}
	}

	return make_tuple(low, high);
}

// Prints to file the state of the body
void Sphere::PrintState(ofstream& fout) {
	fout << setprecision(4);
	fout << "S," << cent[0] << "," << cent[1] << "," << cent[2] << "," << rad << endl;
}

void Sphere::PrintState(string outfile) {
	ofstream fout(outfile);
	PrintState(fout);
	fout.close();
}

// Primes the body to be checked (finds hexagonal projection on the same plane as the origins of the
// rays)
void Parallelepiped::Prime(vector<double> P, vector<double> V) {
	vector<double> p = cent - (double)0.5 * (side[0] + side[1] + side[2]);
	H = FindHexProj(p, side, V, P);
}

// Return a value in [0, 1] describing how close the ray is to the body, 0 if the ray is at least a
// distance dx from the body, 1 if the ray is at least dx inside the body
double Parallelepiped::Check(Ray& ray, double dx) {
	// Finds smallest distance from all side of the hexagon
	double delta_r = dx; // Sentinel value
	vector<double> point = ray.GetR0();
	for (int i = 1; i < 7; i++) {
		delta_r = min(delta_r, PointSegDist(point, H[i], H[PBCH(i + 1)]));
	}
	// Changes sign of the distance if the ray intersects the hexagon
	if (PointIsInsideT(point, H))
		delta_r = -delta_r;
	return d_to_w(delta_r, dx);
}

// Translates the parallelepiped by shift
void Parallelepiped::Translate(vector<double> shift) {
	if (w.size() > 0)
		rotcent += shift;
	cent += shift;
}

// Rotates the parallelepiped around point rot0 according to rotation matrix rotmat
void Parallelepiped::Rotate(vector<double> rot0, vector<vector<double>> rotmat) {
	RotatePoint(cent, rot0, rotmat);
	for (vector<double>& point : side)
		point = rotmat * point;
}

// Return all 8 vertices of the parallelepiped
vector<vector<double>> Parallelepiped::GetVertices() {
	vector<vector<double>> vertices{};
	for (double i = -0.5; i <= 0.5; i++) {
		for (double j = -0.5; j <= 0.5; j++) {
			for (double k = -0.5; k <= 0.5; k++) {
				vertices.push_back(cent + i * side[0] / 2. + j * side[1] / 2. + k * side[2] / 2.);
			}
		}
	}
	return vertices;
}

// Find smallest bounds containing the parallelepiped throughout movement, return a vector containing lower
// values of x, y, z and a vector containing higher values of x, y, z.
tuple<vector<double>, vector<double>> Parallelepiped::GetBounds(double tmin, double tmax,
																unsigned int nstep) {
	vector<double> low = vector<double>(3, numeric_limits<double>::infinity());
	vector<double> high = vector<double>(3, -numeric_limits<double>::infinity());
	if (nstep == 0)
		return make_tuple(low, high);
	double dt = (tmax - tmin) / nstep;

	for (unsigned int i = 0; i < nstep; i++) {
		double t = tmin + i * dt;
		Move(t);
		// Update bounds
		vector<vector<double>> vertices = GetVertices();
		for (vector<double> vertex : vertices) {
			for (size_t j = 0; j < 3; j++) {
				low[j] = min(low[j], vertex[i]);
				high[j] = max(high[j], vertex[i]);
			}
		}
	}

	return make_tuple(low, high);
}

// Prints to file the state of the body
void Parallelepiped::PrintState(ofstream& fout) {
	fout << setprecision(4);
	fout << "P," << cent[0] << "," << cent[1] << "," << cent[2] << ",";
	for (size_t i = 0; i < side.size(); i++) {
		fout << side[i][0] << "," << side[i][1] << "," << side[i][2];
		if (i + 1 != side.size())
			fout << ",";
	}
	fout << endl;
}

void Parallelepiped::PrintState(string outfile) {
	ofstream fout(outfile);
	PrintState(fout);
	fout.close();
}

// Primes the body to be checked (projects l1 and l2 onto the surface)
void Capsule::Prime(vector<double> p, vector<double> v) {
	H1 = Project(l1, p, v);
	H2 = Project(l2, p, v);
}

// Return a value in [0, 1] describing how close the ray is to the body, 0 if the ray is at least a
// distance dx from the body, 1 if the ray is at least dx inside the body
double Capsule::Check(Ray& ray, double dx) {
	double delta_r = PointSegDist(ray.GetR0(), H1, H2) - rad;
	return d_to_w(delta_r, dx);
}

// Translates the sphere by shift
void Capsule::Translate(vector<double> shift) {
	if (w.size() > 0)
		rotcent += shift;
	l1 += shift;
	l2 += shift;
}

// Rotates the Capsule around point rot0 according to rotation matrix rotmat
void Capsule::Rotate(vector<double> rot0, vector<vector<double>> rotmat) {
	RotatePoint(l1, rot0, rotmat);
	RotatePoint(l2, rot0, rotmat);
}

// Find smallest bounds containing the sphere throughout movement, return a vector containing lower
// values of x, y, z and a vector containing higher values of x, y, z.
tuple<vector<double>, vector<double>> Capsule::GetBounds(double tmin, double tmax,
														 unsigned int nstep) {
	vector<double> low = vector<double>(3, numeric_limits<double>::infinity());
	vector<double> high = vector<double>(3, -numeric_limits<double>::infinity());
	if (nstep == 0)
		return make_tuple(low, high);
	double dt = (tmax - tmin) / nstep;

	for (unsigned int i = 0; i < nstep; i++) {
		double t = tmin + i * dt;
		Move(t);
		// Update bounds
		for (size_t j = 0; j < 3; j++) {
			low[j] = min(min(low[j], l1[j] - rad), l2[j] - rad);
			high[j] = max(max(low[j], l1[j] + rad), l2[j] + rad);
		}
	}

	return make_tuple(low, high);
}

// Prints to file the state of the body
void Capsule::PrintState(ofstream& fout) {
	fout << setprecision(4);
	fout << "C," << l1[0] << "," << l1[1] << "," << l1[2] << "," << l2[0] << "," << l2[1] << ","
		 << l2[2] << "," << rad << endl;
}

void Capsule::PrintState(string outfile) {
	ofstream fout(outfile);
	PrintState(fout);
	fout.close();
}

// Complete ManyBody constructor
ManyBody::ManyBody(const vector<Sphere>& Spheres, const vector<Parallelepiped>& Parallelepipeds,
				   const vector<Capsule>& Capsules)
	: Body() {
	for (Sphere sphere : Spheres)
		AddBody(sphere);
	for (Parallelepiped Parallelepiped : Parallelepipeds)
		AddBody(Parallelepiped);
	for (Capsule capsule : Capsules)
		AddBody(capsule);
}

// Constructor from file
ManyBody::ManyBody(string filename) : Body() {
	ifstream file(filename);
	if (!file.is_open()) {
		std::cerr << "Error opening file: " << filename << std::endl;
		return;
	}

	string line;

	// Get bodies
	while (getline(file, line)) {
		istringstream iss(line);
		string type, name, superbody;
		iss >> type;

		if (type == "Sphere") {
			vector<vector<double>> rot, trans;
			vector<double> cent(3);
			vector<double> rotcent(3), axis(3);
			double rad;
			size_t trans_size, w_size;
			for (size_t i = 0; i < 3; i++)
				iss >> cent[i];
			iss >> rad;

			iss >> name;
			iss >> superbody;
			iss >> w_size;
			vector<double> w(w_size);
			for (size_t i = 0; i < w_size; i++)
				iss >> w[i];
			w = w * (M_PI / 180);
			if (w_size != 0) {
				for (size_t i = 0; i < 3; i++)
					iss >> rotcent[i];
				for (size_t i = 0; i < 3; i++)
					iss >> axis[i];
			}
			iss >> trans_size;
			for (size_t i = 0; i < trans_size; i++) {
				vector<double> temp(3);
				for (size_t j = 0; j < 3; j++) {
					iss >> temp[j];
				}
				trans.push_back(temp);
			}

			AddBody(Sphere(cent, rad, name, rotcent, axis, w, trans));
			if (superbody != "None")
				Attach(name, superbody);
		}

		if (type == "Parallelepiped") {
			vector<vector<double>> rot, trans, sides;
			vector<double> cent(3);
			vector<double> rotcent(3), axis(3);
			size_t trans_size, w_size;
			for (size_t i = 0; i < 3; i++)
				iss >> cent[i];
			for (size_t i = 0; i < 3; i++) {
				vector<double> temp(3);
				for (size_t j = 0; j < 3; j++) {
					iss >> temp[j];
				}
				sides.push_back(temp);
			}

			iss >> name;
			iss >> superbody;
			iss >> w_size;
			vector<double> w(w_size);
			for (size_t i = 0; i < w_size; i++)
				iss >> w[i];
			w = w * (M_PI / 180);
			if (w_size != 0) {
				for (size_t i = 0; i < 3; i++)
					iss >> rotcent[i];
				for (size_t i = 0; i < 3; i++)
					iss >> axis[i];
			}
			iss >> trans_size;
			for (size_t i = 0; i < trans_size; i++) {
				vector<double> temp(3);
				for (size_t j = 0; j < 3; j++) {
					iss >> temp[j];
				}
				trans.push_back(temp);
			}

			AddBody(Parallelepiped(cent, sides, name, rotcent, axis, w, trans));
			if (superbody != "None")
				Attach(name, superbody);
		}

		if (type == "Capsule") {
			vector<vector<double>> rot, trans;
			vector<double> cent(3), l1(3), l2(3);
			vector<double> rotcent(3), axis(3);
			double rad;
			size_t trans_size, w_size;
			for (size_t i = 0; i < 3; i++)
				iss >> l1[i];
			for (size_t i = 0; i < 3; i++)
				iss >> l2[i];
			iss >> rad;

			iss >> name;
			iss >> superbody;
			iss >> w_size;
			vector<double> w(w_size);
			for (size_t i = 0; i < w_size; i++)
				iss >> w[i];
			w = w * (M_PI / 180);
			if (w_size != 0) {
				for (size_t i = 0; i < 3; i++)
					iss >> rotcent[i];
				for (size_t i = 0; i < 3; i++)
					iss >> axis[i];
			}
			iss >> trans_size;
			for (size_t i = 0; i < trans_size; i++) {
				vector<double> temp(3);
				for (size_t j = 0; j < 3; j++) {
					iss >> temp[j];
				}
				trans.push_back(temp);
			}

			AddBody(Capsule(l1, l2, rad, name, rotcent, axis, w, trans));
			if (superbody != "None")
				Attach(name, superbody);
		}
	}

	file.close();
}

// Destructor
ManyBody::~ManyBody() {
	// cout << "Destroying ManyBody" << endl;
	// Delete all dynamically allocated Body objects
	for (Body* body : bodies) {
		delete body;
	}
}

// Primes the body to be checked. Primes each body
void ManyBody::Prime(vector<double> p, vector<double> v) {
	for (Body* body : bodies)
		body->Prime(p, v);
}

// Return a value in [0, 1] describing how close the ray is to the body, 0 if the ray is at least a
// distance dx from the body, 1 if the ray is at least dx inside the body
double ManyBody::Check(Ray& ray, double dx) {
	double w = 0;
	// Iterates over all bodies or until a body makes full contact
	for (Body* body : bodies) {
		w = max(w, body->Check(ray, dx));
		if (w == 1)
			return w;
	}
	return w;
}

// Translates the sphere by shift
void ManyBody::Translate(vector<double> shift) {
	if (w.size() > 0)
		rotcent += shift;
	for (Body* body : bodies)
		body->Translate(shift);
}

// Rotates the Capsule around point rot0 according to rotation matrix rotmat
void ManyBody::Rotate(vector<double> rot0, vector<vector<double>> rotmat) {
	for (Body* body : bodies)
		body->Rotate(rot0, rotmat);
}

// Time evolution of the body
void ManyBody::Move(double tNew) {
	// Move parts
	for (Body* body : bodies)
		body->Move(tNew);
}

// Return a pointer to the body with that name in the ManyBody
Body* ManyBody::Find(string name) {
	for (Body* body : bodies)
		if (body->GetName() == name)
			return body;

	cout << name << " not found!" << endl;
	return nullptr;
}

// Attaches the sub-body to the super-body
void ManyBody::Attach(Body SubBody, string SuperName) {
	Body* Super = Find(SuperName);
	if (Super)
		Super->AddSubBody(SubBody);
}

void ManyBody::Attach(string SubName, string SuperName) {
	Body* Super = Find(SuperName);
	Body* Sub = Find(SubName);
	if (Super and Sub)
		Super->AddSubBody(*Sub);
}

// Find smallest bounds containing the object throughout movement, return a vector containing lower
// values of x, y, z and a vector containing higher values of x, y, z.
tuple<vector<double>, vector<double>> ManyBody::GetBounds(double tmin, double tmax,
														  unsigned int nstep) {
	vector<double> low = vector<double>(3, numeric_limits<double>::infinity());
	vector<double> high = vector<double>(3, -numeric_limits<double>::infinity());
	if (nstep == 0)
		return make_tuple(low, high);
	double dt = (tmax - tmin) / nstep;

	for (unsigned int i = 0; i < nstep; i++) {
		double t = tmin + i * dt;
		Move(t);
		// Update bounds
		for (Body* body : bodies) {
			vector<double> bodyLow;
			vector<double> bodyHigh;
			tie(bodyLow, bodyHigh) = body->GetBounds(t, t, 1); // Get bounds without moving body
			for (size_t j = 0; j < 3; j++) {
				low[j] = min(low[j], bodyLow[j]);
				high[j] = max(high[j], bodyHigh[j]);
			}
		}
	}

	return make_tuple(low, high);
}

// Prints to file the state (all the bodies and their parameters)
void ManyBody::PrintState(ofstream& fout) {
	for (Body* body : bodies) {
		body->PrintState(fout);
	}
}

void ManyBody::PrintState(string outfile) {
	ofstream fout(outfile);
	PrintState(fout);
	fout.close();
}