#include "body.h"

#include "RainFunctions.h"
#include "ray.h"
#include "tinyxml2.h"

using namespace std;
using namespace tinyxml2;

// Virtual destructor
Body::~Body() {
	// cout << "Destroying " << name << endl;
}

// Evaluate theta at a time t
double Body::getTheta(double t) {
	double res = theta0;
	// Add terms
	int imax = min(w.size(), wPhi.size()) - 1;
	for (int i = 0; i <= imax; i++) {
		res += w[i] * sin(2 * M_PI * (i + 1) * t + wPhi[i]);
	}
	return res;
}

// Evaluate delta at a time t
vector<double> Body::getDelta(double t) {
	vector<double> res = delta0;
	// Add terms
	int imax = min(trans.size(), transPhi.size()) - 1;
	for (int i = 0; i <= imax; i++) {
		res += trans[i] * sin(2 * M_PI * (i + 1) + transPhi[i]);
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
			high[j] = max(high[j], cent[j] + rad);
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
				vertices.push_back(cent + i * side[0] + j * side[1] + k * side[2]);
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
				low[j] = min(low[j], vertex[j]);
				high[j] = max(high[j], vertex[j]);
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
			low[j] = min(min(l1[j] - rad, l2[j] - rad), low[j]);
			high[j] = max(max(l1[j] + rad, l2[j] + rad), high[j]);
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

// Constuctor
ManyBody::ManyBody(string filename) : Body() {
	XMLDocument doc;
	if (doc.LoadFile(filename.c_str()) != XML_SUCCESS) {
		throw logic_error(string("Error loading XML file: " + filename));
	}

	XMLElement* root = doc.FirstChildElement("ManyBody");
	if (!root) {
		throw logic_error(string("No <ManyBody> root found."));
	}

	for (XMLElement* bodyElem = root->FirstChildElement("Body"); bodyElem;
		 bodyElem = bodyElem->NextSiblingElement("Body")) {
		string type = bodyElem->Attribute("type");

		XMLElement* nameElem = bodyElem->FirstChildElement("Name");
		if (!nameElem || !nameElem->GetText())
			throw logic_error(string("Missing or empty <Name> in <Body> starting at line " +
									 to_string(bodyElem->GetLineNum())));
		name = nameElem->GetText();

		double theta0 = 0;
		vector<double> w, wPhi, rotcent(3), axis(3), delta0(3), transPhi;
		vector<vector<double>> trans;

		XMLElement* rotElem = bodyElem->FirstChildElement("Rotation");
		if (rotElem) {
			XMLElement* theta0Elem = rotElem->FirstChildElement("Theta0");
			if (theta0Elem && theta0Elem->GetText()) {
				vector<double> theta0vals = parseDoubles(theta0Elem->GetText());
				if (theta0vals.size() != 1)
					throw logic_error(string("Invalid <Theta0> value in <Rotation> in " + name +
											 ", must be a single number"));
				theta0 = theta0vals[0] * M_PI / 180.0;

				XMLElement* wElem = rotElem->FirstChildElement("W");
				if (!wElem)
					throw logic_error(string("Missing <W> in <Rotation> in " + name));
				w = parseDoubles(wElem->GetText());
				for (double& val : w)
					val *= M_PI / 180.0;

				XMLElement* wPhiElem = rotElem->FirstChildElement("WPhi");
				if (!wPhiElem)
					throw logic_error(string("Missing <WPhi> in <Rotation> in " + name));
				wPhi = parseDoubles(wPhiElem->GetText());
				for (double& val : wPhi)
					val *= M_PI / 180.0;
				if (w.size() != wPhi.size())
					throw logic_error(
						string("<W> and <WPhi> in " + name +
							   " have a different number of elements, must be equal"));

				rotcent = SafeGet3Vec(rotElem, "RotCenter", name);

				axis = SafeGet3Vec(rotElem, "Axis", name);
			}
		}

		XMLElement* transElem = bodyElem->FirstChildElement("Translations");
		if (transElem) {
			XMLElement* delta0Elem = transElem->FirstChildElement("Delta0");
			if (delta0Elem && delta0Elem->GetText()) {
				vector<double> delta0vals = parseDoubles(delta0Elem->GetText());
				if (delta0vals.size() != 3)
					throw logic_error(string("Invalid <Delta0> value in <Translations> in " + name +
											 ", must be 3 number"));
				delta0 = delta0vals;
			}

			for (XMLElement* t = transElem->FirstChildElement("Translation"); t;
				 t = t->NextSiblingElement("Translation")) {
				vector<double> tran = parseDoubles(t->GetText());
				if (tran.size() != 3)
					throw logic_error(string("Invalid <\"Translation\"> value in " + name +
											 ", must be three numbers"));
				trans.push_back(tran);
			}

			XMLElement* transPhiElem = transElem->FirstChildElement("TransPhi");
			if (!transPhiElem)
				throw logic_error(string("Missing <TransPhi> in <Translation> in " + name));
			transPhi = parseDoubles(transPhiElem->GetText());
			for (double& val : transPhi)
				val *= M_PI / 180.0;
			if (trans.size() != transPhi.size())
				throw logic_error(string(
					"Number of <Rotation> elements and number of elements in <TransPhi> in " +
					name + " is different, must be equal"));
		}

		if (type == "Sphere") {
			vector<double> cent = SafeGet3Vec(bodyElem, "Center", name);
			if (cent.size() != 3)
				throw logic_error(
					string("Invalid <Center> value in " + name + ", must be three numbers"));

			double rad = atof(bodyElem->FirstChildElement("Radius")->GetText());
			if (rad <= 0)
				throw logic_error(
					string("Invalid value of <Radius> in " + name + ", must be a positive number"));

			AddBody(
				Sphere(cent, rad, name, rotcent, axis, w, wPhi, trans, transPhi, theta0, delta0));
		} else if (type == "Parallelepiped") {
			vector<double> cent = SafeGet3Vec(bodyElem, "Center", name);

			vector<vector<double>> sides;
			if (!bodyElem->FirstChildElement("Sides"))
				throw logic_error(string("Missing or empty <Sides> in " + name));
			for (XMLElement* s = bodyElem->FirstChildElement("Sides")->FirstChildElement("Side"); s;
				 s = s->NextSiblingElement("Side")) {
				sides.push_back(parseDoubles(s->GetText()));
			}
			if (sides.size() != 3)
				throw logic_error(string("Invalid number (" + to_string(sides.size()) +
										 ") of <Side> elements in " + name + ", must be 3"));

			AddBody(Parallelepiped(cent, sides, name, rotcent, axis, w, wPhi, trans, transPhi,
								   theta0, delta0));
		} else if (type == "Capsule") {
			vector<double> l1 = SafeGet3Vec(bodyElem, "L1", name);
			vector<double> l2 = SafeGet3Vec(bodyElem, "L2", name);

			if (bodyElem->FirstChildElement("Radius") == nullptr)
				throw logic_error(string("Missing value of <Radius> in " + name));

			double rad = atof(bodyElem->FirstChildElement("Radius")->GetText());
			if (rad <= 0)
				throw logic_error(
					string("Invalid value of <Radius> in " + name + ", must be a positive number"));

			AddBody(Capsule(l1, l2, rad, name, rotcent, axis, w, wPhi, trans, transPhi, theta0,
							delta0));
		} else
			throw logic_error(string("Missing or invalid Body type in " + name +
									 ", must be one of "
									 "\"Sphere\", \"Parallelepiped\" or \"Capsule\"!"));

		XMLElement* superElement = bodyElem->FirstChildElement("Superbody");
		if (superElement) {
			Attach(name, superElement->GetText());
		}
	}
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

	// cout << name << " not found!" << endl;
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
	if (Super == nullptr)
		cout << "Warning: Couldn't find SuperBody \"" + SuperName + "\"!" << endl;
	if (Sub == nullptr)
		cout << "Warning: Couldn't find SubBody \"" + SubName + "\"!" << endl;
	if (Super and Sub) {
		Super->AddSubBody(*Sub);
	}
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