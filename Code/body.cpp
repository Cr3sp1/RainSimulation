#include "body.h"

#include "RainFunctions.h"
#include "ray.h"

using namespace std;

// Virtual destructor
Body::~Body() {
	// cout << "Destroying " << name << endl;
}

// Time evolution of the body in its own frame of reference, also propagates to the sub-bodies
void Body::Move(double T) {
	if (T == t)
		return;

	// Calculate the total translation for the step
	vector<double> delta({0, 0, 0});
	// Add sin terms
	for (size_t i = 0; 2 * i < trans.size(); i++) {
		delta += trans[2 * i] * (sin(T * 2 * M_PI / (i + 1)) - sin(t * 2 * M_PI / (i + 1)));
	}
	// Add cos terms
	for (size_t i = 0; 2 * i + 1 < trans.size(); i++) {
		delta += trans[2 * i + 1] * (cos(T * 2 * M_PI / (i + 1)) - cos(t * 2 * M_PI / (i + 1)));
	}
	// Translate
	if (w.size() > 0)
		rotcent += delta;

	// Generate rotation matrix
	vector<vector<double>> rotmat;
	if (w.size() > 0) {
		double theta = 0;
		// Add sin terms
		for (size_t i = 0; 2 * i < w.size(); i++) {
			theta += w[2 * i] * (sin(T * 2 * M_PI / (i + 1)) - sin(t * 2 * M_PI / (i + 1)));
		}
		// Add cos terms
		for (size_t i = 0; 2 * i + 1 < w.size(); i++) {
			theta += w[2 * i + 1] * (cos(T * 2 * M_PI / (i + 1)) - cos(t * 2 * M_PI / (i + 1)));
		}
		rotmat = RotMat(rotax, theta);
	} else {
		rotmat = IdMat(3);
	}

	// Move sub-bodies
	for (Body* body : SubBodies)
		body->BeMoved(delta, rotcent, rotmat);

	t = T;
}

// Time evolution caused by the super-body, affects the whole frame of reference, also propagates to
// the sub-bodie
void Body::BeMoved(vector<double> Delta, vector<double> RotCent, vector<vector<double>> Rotmat) {
	// Translate
	if (w.size() != 0)
		rotcent += Delta;

	// Rotate
	if (w.size() != 0) {
		Rotate(rotcent, RotCent, Rotmat);
		rotax = Rotmat * rotax;
	}
	for (vector<double>& vec : trans)
		vec = Rotmat * vec;

	// Move sub-bodies
	for (Body* body : SubBodies)
		body->BeMoved(Delta, RotCent, Rotmat);
}

// Finds smallest box around body
void Body::FindBox(vector<double>& min, vector<double>& max) {
	if (min.size() != 3 or max.size() != 3) {
		cout << "Error in FindBox: min and max must be of size 3" << endl;
	}
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

// Checks if the Sphere is making contact with a ray
bool Sphere::Check(Ray& ray) {
	vector<double> Xrel = ray.GetR0() - Hcent;
	if (Xrel * Xrel <= rad2) { // Use std::Norm instead, store rad2 only once
		return true;
	}
	return false;
}

// Returns a value in [0, 1] describing how close the ray is to the body, 0 if the ray is at least a
// distance dx from the body, 1 if the ray is at least dx inside the body
double Sphere::CheckSmooth(Ray& ray, double dx) {
	vector<double> Xrel = ray.GetR0() - Hcent;
	double delta_r = Norm(Xrel) - rad;
	return smooth_w(delta_r, dx);
}

// Time evolution of the body in its own frame of reference, also propagates to the sub-bodies
void Sphere::Move(double T) {
	if (T == t)
		return;

	// Calculate the total translation for the step
	vector<double> delta({0, 0, 0});
	// Add sin terms
	for (size_t i = 0; 2 * i < trans.size(); i++) {
		delta += trans[2 * i] * (sin(T * 2 * M_PI / (i + 1)) - sin(t * 2 * M_PI / (i + 1)));
	}
	// Add cos terms
	for (size_t i = 0; 2 * i + 1 < trans.size(); i++) {
		delta += trans[2 * i + 1] * (cos(T * 2 * M_PI / (i + 1)) - cos(t * 2 * M_PI / (i + 1)));
	}
	// Translate
	if (w.size() > 0)
		rotcent += delta;
	;
	cent += delta;

	// Generate rotation matrix
	vector<vector<double>> rotmat;
	if (w.size() > 0) {
		double theta = 0;
		// Add sin terms
		for (size_t i = 0; 2 * i < w.size(); i++) {
			theta += w[2 * i] * (sin(T * 2 * M_PI / (i + 1)) - sin(t * 2 * M_PI / (i + 1)));
		}
		// Add cos terms
		for (size_t i = 0; 2 * i + 1 < w.size(); i++) {
			theta += w[2 * i + 1] * (cos(T * 2 * M_PI / (i + 1)) - cos(t * 2 * M_PI / (i + 1)));
		}
		rotmat = RotMat(rotax, theta);
	} else {
		rotmat = IdMat(3);
	}

	// Rotate
	Rotate(cent, rotcent, rotmat);

	// Move sub-bodies
	for (Body* body : SubBodies)
		body->BeMoved(delta, rotcent, rotmat);

	t = T;
}

// Time evolution caused by the super-body, affects the whole frame of reference, also propagates to
// the sub-bodie
void Sphere::BeMoved(vector<double> Delta, vector<double> RotCent, vector<vector<double>> Rotmat) {
	// Translate
	if (w.size() != 0)
		rotcent += Delta;
	cent += Delta;

	// Rotate
	// Rotate
	if (w.size() != 0) {
		Rotate(rotcent, RotCent, Rotmat);
		rotax = Rotmat * rotax;
	}
	for (vector<double>& vec : trans)
		vec = Rotmat * vec;
	Rotate(cent, RotCent, Rotmat);

	// Move sub-bodies
	for (Body* body : SubBodies)
		body->BeMoved(Delta, RotCent, Rotmat);
}

// Finds smallest box around body
void Sphere::FindBox(vector<double>& min, vector<double>& max) {
	if (min.size() != 3 or max.size() != 3) {
		cout << "Error in FindBox: min and max must be of size 3" << endl;
	}
	for (size_t i = 0; i < 3; i++) {
		if (cent[i] - rad < min[i])
			min[i] = cent[i] - rad;
		if (cent[i] + rad > max[i])
			max[i] = cent[i] + rad;
	}
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

// Checks if the body is making contact with a ray
bool Parallelepiped::Check(Ray& ray) {
	if (PointIsInsideT(ray.GetR0(), H)) { // Use less triangles (4 not 6)
		return true;
	}
	return false;
}

// Returns a value in [0, 1] describing how close the ray is to the body, 0 if the ray is at least a
// distance dx from the body, 1 if the ray is at least dx inside the body
double Parallelepiped::CheckSmooth(Ray& ray, double dx) {
	// Finds smallest distance from all side of the hexagon
	double delta_r = dx; // Sentinel value
	vector<double> point = ray.GetR0();
	for (int i = 1; i < 7; i++) {
		delta_r = min(delta_r, PointSegDist(point, H[i], H[PBCH(i + 1)]));
	}
	// Changes sign of the distance if the ray intersects the hexagon
	if (PointIsInsideT(point, H))
		delta_r = -delta_r;
	return smooth_w(delta_r, dx);
}

// Time evolution of the body in its own frame of reference, also propagates to the sub-bodies
void Parallelepiped::Move(double T) {
	if (T == t)
		return;

	// Calculate the total translation for the step
	vector<double> delta({0, 0, 0});
	// Add sin terms
	for (size_t i = 0; 2 * i < trans.size(); i++) {
		delta += trans[2 * i] * (sin(T * 2 * M_PI / (i + 1)) - sin(t * 2 * M_PI / (i + 1)));
	}
	// Add cos terms
	for (size_t i = 0; 2 * i + 1 < trans.size(); i++) {
		delta += trans[2 * i + 1] * (cos(T * 2 * M_PI / (i + 1)) - cos(t * 2 * M_PI / (i + 1)));
	}
	// Translate
	if (w.size() > 0)
		rotcent += delta;
	cent += delta;

	// Generate rotation matrix
	vector<vector<double>> rotmat;
	if (w.size() > 0) {
		double theta = 0;
		// Add sin terms
		for (size_t i = 0; 2 * i < w.size(); i++) {
			theta += w[2 * i] * (sin(T * 2 * M_PI / (i + 1)) - sin(t * 2 * M_PI / (i + 1)));
		}
		// Add cos terms
		for (size_t i = 0; 2 * i + 1 < w.size(); i++) {
			theta += w[2 * i + 1] * (cos(T * 2 * M_PI / (i + 1)) - cos(t * 2 * M_PI / (i + 1)));
		}
		rotmat = RotMat(rotax, theta);
	} else {
		rotmat = IdMat(3);
	}

	// Rotate
	for (vector<double>& point : side)
		point = rotmat * point;
	Rotate(cent, rotcent, rotmat);

	// Move sub-bodies
	for (Body* body : SubBodies)
		body->BeMoved(delta, rotcent, rotmat);

	t = T;
}

// Time evolution caused by the super-body, affects the whole frame of reference, also propagates to
// the sub-bodie
void Parallelepiped::BeMoved(vector<double> Delta, vector<double> RotCent,
							 vector<vector<double>> Rotmat) {
	// Translate
	if (w.size() != 0)
		rotcent += Delta;
	cent += Delta;

	// Rotate
	if (w.size() != 0) {
		Rotate(rotcent, RotCent, Rotmat);
		rotax = Rotmat * rotax;
	}
	for (vector<double>& vec : trans)
		vec = Rotmat * vec;
	for (vector<double>& point : side)
		point = Rotmat * point;
	Rotate(cent, RotCent, Rotmat);

	// Move sub-bodies
	for (Body* body : SubBodies)
		body->BeMoved(Delta, RotCent, Rotmat);
}

// Returns all 8 vertices of the parallelepiped
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

// Finds smallest box around body
void Parallelepiped::FindBox(vector<double>& min, vector<double>& max) {
	if (min.size() != 3 or max.size() != 3) {
		cout << "Error in FindBox: min and max must be of size 3" << endl;
	}
	vector<vector<double>> vertices = GetVertices();
	for (vector<double> vertex : vertices) {
		for (size_t i = 0; i < 3; i++) {
			if (vertex[i] < min[i])
				min[i] = vertex[i];
			if (vertex[i] > max[i])
				max[i] = vertex[i];
		}
	}
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

// Checks if the Capsule is making contact with a ray
bool Capsule::Check(Ray& ray) {
	if (PointSegDist(ray.GetR0(), H1, H2) <= rad) {
		return true;
	}
	return false;
}

// Returns a value in [0, 1] describing how close the ray is to the body, 0 if the ray is at least a
// distance dx from the body, 1 if the ray is at least dx inside the body
double Capsule::CheckSmooth(Ray& ray, double dx) {
	double delta_r = PointSegDist(ray.GetR0(), H1, H2) - rad;
	return smooth_w(delta_r, dx);
}

// Time evolution of the body in its own frame of reference, also propagates to the sub-bodies
void Capsule::Move(double T) {
	if (T == t)
		return;

	// Calculate the total translation for the step
	vector<double> delta({0, 0, 0});
	// Add sin terms
	for (size_t i = 0; 2 * i < trans.size(); i++) {
		delta += trans[2 * i] * (sin(T * 2 * M_PI / (i + 1)) - sin(t * 2 * M_PI / (i + 1)));
	}
	// Add cos terms
	for (size_t i = 0; 2 * i + 1 < trans.size(); i++) {
		delta += trans[2 * i + 1] * (cos(T * 2 * M_PI / (i + 1)) - cos(t * 2 * M_PI / (i + 1)));
	}
	// Translate
	if (w.size() > 0)
		rotcent += delta;
	l1 += delta;
	l2 += delta;

	// Generate rotation matrix
	vector<vector<double>> rotmat;
	if (w.size() > 0) {
		double theta = 0;
		// Add sin terms
		for (size_t i = 0; 2 * i < w.size(); i++) {
			theta += w[2 * i] * (sin(T * 2 * M_PI / (i + 1)) - sin(t * 2 * M_PI / (i + 1)));
		}
		// Add cos terms
		for (size_t i = 0; 2 * i + 1 < w.size(); i++) {
			theta += w[2 * i + 1] * (cos(T * 2 * M_PI / (i + 1)) - cos(t * 2 * M_PI / (i + 1)));
		}
		rotmat = RotMat(rotax, theta);
	} else {
		rotmat = IdMat(3);
	}

	// Rotate
	Rotate(l1, rotcent, rotmat);
	Rotate(l2, rotcent, rotmat);

	// Move sub-bodies
	for (Body* body : SubBodies)
		body->BeMoved(delta, rotcent, rotmat);

	t = T;
}

// Time evolution caused by the super-body, affects the whole frame of reference, also propagates to
// the sub-bodie
void Capsule::BeMoved(vector<double> Delta, vector<double> RotCent, vector<vector<double>> Rotmat) {
	// Translate
	if (w.size() != 0)
		rotcent += Delta;
	l1 += Delta;
	l2 += Delta;

	// Rotate
	if (w.size() != 0) {
		Rotate(rotcent, RotCent, Rotmat);
		rotax = Rotmat * rotax;
	}
	for (vector<double>& vec : trans)
		vec = Rotmat * vec;
	Rotate(l1, RotCent, Rotmat);
	Rotate(l2, RotCent, Rotmat);
	// Move sub-bodies
	for (Body* body : SubBodies)
		body->BeMoved(Delta, RotCent, Rotmat);
}

// Finds smallest box around body
void Capsule::FindBox(vector<double>& min, vector<double>& max) {
	if (min.size() != 3 or max.size() != 3) {
		cout << "Error in FindBox: min and max must be of size 3" << endl;
	}
	for (size_t i = 0; i < 3; i++) {
		if (l1[i] - rad < min[i])
			min[i] = l1[i] - rad;
		if (l1[i] + rad > max[i])
			max[i] = l1[i] + rad;
	}
	for (size_t i = 0; i < 3; i++) {
		if (l2[i] - rad < min[i])
			min[i] = l2[i] - rad;
		if (l2[i] + rad > max[i])
			max[i] = l2[i] + rad;
	}
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

// Checks if the ManyBody is making contact with a ray
bool ManyBody::Check(Ray& ray) {
	// Iterates over all bodies or until a body makes contact
	for (Body* body : bodies)
		if (body->Check(ray))
			return true;
	return false;
}

// Returns a value in [0, 1] describing how close the ray is to the body, 0 if the ray is at least a
// distance dx from the body, 1 if the ray is at least dx inside the body
double ManyBody::CheckSmooth(Ray& ray, double dx) {
	double w = 0;
	// Iterates over all bodies or until a body makes full contact
	for (Body* body : bodies) {
		w = max(w, body->CheckSmooth(ray, dx));
		if (w == 1)
			return w;
	}
	return w;
}

// Time evolution of the body
void ManyBody::Move(double T) {
	if (T == t)
		return;
	// Move parts
	for (Body* body : bodies)
		body->Move(T);

	t = T;
}

// Returns a pointer to the body with that name in the ManyBody
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

// Finds smallest box around body
void ManyBody::FindBox(vector<double>& min, vector<double>& max) {
	if (min.size() != 3 or max.size() != 3) {
		cout << "Error in FindBox: min and max must be of size 3" << endl;
	}
	for (Body* body : bodies) {
		body->FindBox(min, max);
	}
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