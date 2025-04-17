#ifndef __Body_h__
#define __Body_h__

#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include "RainFunctions.h"
#include "VectorOperations.h"

using namespace std;

// Forward declaration
class Ray;
class ProjSurface;
class ManyBody;

// Body class
class Body {
  protected:
	// Name (optional)
	string name;
	// Time (should go from 0 to 1)
	double t;
	// Center of rotation
	vector<double> rotcent;
	// Axis of rotation
	vector<double> rotax;
	// Amplitude of oscillation in radiants (coefficients in fourier series
	// expansion of periodic motion, structured as [0]sin [1]cos [2]sin ...)
	vector<double> w;
	// Vectors used for periodic translation in time evolution (coefficients in
	// sin expasion of periodic motion)
	vector<vector<double>> trans;
	// Pointers to all the bodies that move relative to this one
	vector<Body*> SubBodies;

  public:
	// Default constructor
	Body() : t(0), w(0){};
	// Dynamic constructor
	Body(string Name, vector<double> Rotcent, vector<double> Rotax, vector<double> W,
		 vector<vector<double>> Trans)
		: name(Name), t(0), rotcent(Rotcent), rotax(Rotax), w(W), trans(Trans) {}
	// Virtual destructor
	virtual ~Body();
	// Primes the body to be checked. p is a point on the surface containing the ray origins and v
	// is the relative velocity
	virtual void Prime(vector<double> p, vector<double> v) {}
	// Return a value in [0, 1] describing how close the ray is to the body, 0 if the ray is at
	// least a distance dx from the body, 1 if the ray is at least dx inside the body
	virtual double Check(Ray& ray, double dx) { return 0.; }
	// Translates the body by Delta
	virtual void Translate(vector<double> Delta) {
		if (w.size() > 0)
			rotcent += Delta;
	}
	// Rotates the body around point Rot0 according to rotation matrix Rotmat
	virtual void Rotate(vector<double> Rot0, vector<vector<double>> Rotmat) {}
	// Time evolution of the body in its own frame of reference, also propagates to the sub-bodies
	virtual void Move(double T);
	// Time evolution caused by the super-body, affects the whole frame of reference, also
	// propagates to the sub-bodies
	virtual void BeMoved(vector<double> Delta, vector<double> Rot0, vector<vector<double>> Rotmat);
	// Adds a Body to SubBodies
	virtual void AddSubBody(Body& SubBody) { SubBodies.push_back(&SubBody); }
	// Attaches the Body to a SuperBody
	virtual void AttachTo(Body& SupBody) { SupBody.AddSubBody(*this); }
	// Get stuff
	virtual string GetName() { return name; }
	// Find smallest bounds containing the body throughout movement, return a vector containing lower
	// values of x, y, z and a vector containing higher values of x, y, z.
	virtual tuple<vector<double>, vector<double>> GetBounds(double tmin, double tmax,
															unsigned int nstep);
	// Return box starting on origin that contains the body with allowance epsilon and move body in it
	vector<double> GetBox(double tmin, double tmax, unsigned int nstep, double epsilon);
	// Prints to file the state of the body
	virtual void PrintState(ofstream& fout);
	virtual void PrintState(string outfile);
};

// Sphere class
class Sphere : public Body {
  protected:
	vector<double> cent;  // Position of the center of the sphere
	double rad;			  // Radius of the sphere
	double rad2;		  // Square radius of the sphere
	vector<double> Hcent; // Projection of the center of the sphere

  public:
	// Complete static constructor
	Sphere(vector<double> center, double radius)
		: Body(), cent(center), rad(radius), rad2(radius * radius) {}
	// Complete dynamic constructor
	Sphere(vector<double> center, double radius, string Name, vector<double> Rotcent,
		   vector<double> Rotax, vector<double> W, vector<vector<double>> Trans)
		: Body(Name, Rotcent, Rotax, W, Trans), cent(center), rad(radius), rad2(radius * radius) {}
	// Primes the body to be checked. p is a point on the surface containing the ray origins and v
	// is the relative velocity
	void Prime(vector<double> p, vector<double> v) override;
	// Return a value in [0, 1] describing how close the ray is to the body, 0 if the ray is at
	// least a distance dx from the body, 1 if the ray is at least dx inside the body
	double Check(Ray& ray, double dx) override;
	// Translates the sphere by Delta
	void Translate(vector<double> Delta) override;
	// Rotates the sphere around point Rot0 according to rotation matrix Rotmat
	void Rotate(vector<double> Rot0, vector<vector<double>> Rotmat) override;
	// Gets stuff
	vector<double> GetCent() { return cent; }
	double GetRad() { return rad; }
	// Find smallest bounds containing the sphere throughout movement, return a vector containing lower
	// values of x, y, z and a vector containing higher values of x, y, z.
	tuple<vector<double>, vector<double>> GetBounds(double tmin, double tmax,
													unsigned int nstep) override;
	// Prints to file the state of the body
	void PrintState(ofstream& fout) override;
	void PrintState(string outfile) override;
};

// Parallelepiped class
class Parallelepiped : public Body {
  protected:
	vector<double> cent;		 // Position of the center
	vector<vector<double>> side; // Sides of the parallelepiped
	vector<vector<double>> H;	 // Hexagonal projection on surface

  public:
	// Complete static constructor
	Parallelepiped(vector<double> Center, vector<vector<double>> Side)
		: Body(), cent(Center), side(Side) {}
	// Complete dynamic constructor
	Parallelepiped(vector<double> Center, vector<vector<double>> Side, string Name,
				   vector<double> Rotcent, vector<double> Rotax, vector<double> W,
				   vector<vector<double>> Trans)
		: Body(Name, Rotcent, Rotax, W, Trans), cent(Center), side(Side) {}
	// Primes the body to be checked. p is a point on the surface containing the ray origins and v
	// is the relative velocity
	void Prime(vector<double> p, vector<double> v) override;
	// Return a value in [0, 1] describing how close the ray is to the body, 0 if the ray is at
	// least a distance dx from the body, 1 if the ray is at least dx inside the body
	double Check(Ray& ray, double dx) override;
	// Translates the parallelepiped by Delta
	void Translate(vector<double> Delta) override;
	// Rotates the parallelepiped around point Rot0 according to rotation matrix Rotmat
	void Rotate(vector<double> Rot0, vector<vector<double>> Rotmat) override;
	// Gets stuff
	vector<double> GetCent() { return cent; }
	vector<vector<double>> GetSide() { return side; }
	vector<vector<double>> GetVertices();
	// Find smallest bounds containing the parallelepiped throughout movement, return a vector containing
	// lower values of x, y, z and a vector containing higher values of x, y, z.
	tuple<vector<double>, vector<double>> GetBounds(double tmin, double tmax,
													unsigned int nstep) override;
	// Prints to file the state of the body
	void PrintState(ofstream& fout) override;
	void PrintState(string outfile) override;
};

// Capsule class
class Capsule : public Body {
  protected:
	vector<double> l1, l2; // Position of the two extremes of the axis
	double rad;			   // Radius of the sphere
	vector<double> H1, H2; // Projections of the two extremes of the axise

  public:
	// Complete static constructor
	Capsule(vector<double> L1, vector<double> L2, double Radius)
		: Body(), l1(L1), l2(L2), rad(Radius) {}
	// Complete dynamic constructor
	Capsule(vector<double> L1, vector<double> L2, double Radius, string Name,
			vector<double> Rotcent, vector<double> Rotax, vector<double> W,
			vector<vector<double>> Trans)
		: Body(Name, Rotcent, Rotax, W, Trans), l1(L1), l2(L2), rad(Radius) {}
	// Primes the body to be checked. p is a point on the surface containing theray origins and v is
	// the relative velocity
	void Prime(vector<double> p, vector<double> v) override;
	// Return a value in [0, 1] describing how close the ray is to the body, 0 ifthe ray is at
	// least a distance dx from the body, 1 if the ray is at leastdx inside the body
	double Check(Ray& ray, double dx) override;
	// Translates capsule body by Delta
	void Translate(vector<double> Delta) override;
	// Rotates the capsule around point Rot0 according to rotation matrix Rotmat
	void Rotate(vector<double> Rot0, vector<vector<double>> Rotmat) override;
	// Gets stuff
	vector<double> GetL1() { return l1; }
	vector<double> GetL2() { return l2; }
	double GetRad() { return rad; }
	// Find smallest bounds containing the capsule throughout movement, return a vector containing
	// lower values of x, y, z and a vector containing higher values of x, y, z.
	tuple<vector<double>, vector<double>> GetBounds(double tmin, double tmax,
													unsigned int nstep) override;
	// Prints to file the state of the body
	void PrintState(ofstream& fout) override;
	void PrintState(string outfile) override;
};

// ManyBody class
class ManyBody : public Body {
  protected:
	// Bodies contained in ManyBody
	vector<Body*> bodies;

  public:
	// Empty constructor
	ManyBody() : Body() {}
	// Complete constructor
	ManyBody(const vector<Sphere>& Spheres, const vector<Parallelepiped>& Parallelepipeds,
			 const vector<Capsule>& Capsules);
	// Constructor from file
	ManyBody(string filename);
	// Destructor
	~ManyBody() override;
	// Primes the body to be checked. p is a point on the surface containing theray origins and v is
	// the relative velocity
	void Prime(vector<double> p, vector<double> v) override;
	// Return a value in [0, 1] describing how close the ray is to the body, 0 ifthe ray is at
	// least a distance dx from the body, 1 if the ray is at leastdx inside the body
	double Check(Ray& ray, double dx) override;
	// Translates all the bodies by Delta
	void Translate(vector<double> Delta) override;
	// Rotates all the bodies around point Rot0 according to rotation matrix Rotmat
	virtual void Rotate(vector<double> Rot0, vector<vector<double>> Rotmat) override;
	// Time evolution of all the bodies
	void Move(double T) override;
	// Time evolution caused by the super-body, affects the whole frame ofreference, also propagates
	// to the sub-bodies
	void BeMoved(vector<double> Delta, vector<double> Rot0,
				 vector<vector<double>> Rotmat) override {};
	// Add bodies
	void AddBody(Sphere sphere) { bodies.push_back(new Sphere(sphere)); }
	void AddBody(Parallelepiped parallelepiped) {
		bodies.push_back(new Parallelepiped(parallelepiped));
	}
	void AddBody(Capsule capsule) { bodies.push_back(new Capsule(capsule)); }
	// Pointer to the body with that name
	Body* Find(string name);
	// Attaches the sub-body to the super-body
	void Attach(Body SubBody, string SuperName);
	void Attach(string SubName, string SuperName);
	// Find smallest box containing the object throughout movement, return higher values of x, y, z. First
	// translate the object so that lower values of x, y, z are all 0. Increase box dimensions by epsilon
	tuple<vector<double>, vector<double>> GetBounds(double tmin, double tmax,
													unsigned int nstep) override;
	// Prints to file the state (all the bodies and their parameters)
	void PrintState(ofstream& fout) override;
	void PrintState(string outfile) override;
};

#endif // __Body_h__
