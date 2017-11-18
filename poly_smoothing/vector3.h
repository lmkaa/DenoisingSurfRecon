#ifndef __VECTOR3_H__
#define __VECTOR3_H__

#include <cmath>
#include <iostream>
using namespace std;

template <class C> class Vector3;
typedef Vector3<float> Vector3f;
typedef Vector3<double> Vector3d;

template <class C>
bool operator== (const Vector3<C> & v, const Vector3<C> & r) { return (v.X()==r.X() && v.Y()==r.Y() && v.Z()==r.Z()); }

template <class C>
bool operator!= (const Vector3<C> & v, const Vector3<C> & r) { return (v.X()!=r.X() || v.Y()!=r.Y() || v.Z()==r.Z()); }

template <class C>
C __min(C a, C b) {
	if (a < b)
		return a;
	else
		return b;
};

template <class C>
C __max(C a, C b) {
	if (a > b)
		return a;
	else
		return b;
};

template <class C>
class Vector3 {
private:
	C v[3];

public:
	// constructors
	Vector3() { v[0] = v[1]= v[2] = 0; }
	Vector3(const C & value) { v[0] = v[1] = v[2] = value; }
	Vector3(const C & a, const C & b, const C & c) {
		v[0] = a;
		v[1] = b;
		v[2] = c;
	}
	Vector3(C arr[]) {
		v[0] = arr[0];
		v[1] = arr[1];
		v[2] = arr[2];
	}
	Vector3(const C arr[]) {
		v[0] = arr[0];
		v[1] = arr[1];
		v[2] = arr[2];
	}
	Vector3(const Vector3<C> & right) {
		v[0] = right.v[0];
		v[1] = right.v[1];
		v[2] = right.v[2];
	}

	const C & X() const { return v[0]; }
	const C & Y() const { return v[1]; }
	const C & Z() const { return v[2]; }
	C & X() { return v[0]; }
	C & Y() { return v[1]; }
	C & Z() { return v[2]; }

	C & operator[] (int index) { return v[index]; }
	const C & operator[] (int index) const { return v[index]; }

	Vector3<C> operator+ (const Vector3<C> & r) const { return Vector3<C>(v[0]+r[0], v[1]+r[1], v[2]+r[2]); }
	Vector3<C> operator- (const Vector3<C> & r) const { return Vector3<C>(v[0]-r[0], v[1]-r[1], v[2]-r[2]); }
	Vector3<C> operator* (const C & r) const { return Vector3<C>(v[0]*r, v[1]*r, v[2]*r); }
	Vector3<C> operator/ (const C & r) const { return Vector3<C>(v[0]/r, v[1]/r, v[2]/r); }
	Vector3<C> operator-() const { return Vector3<C>(-v[0], -v[1], -v[2]); }
	Vector3<C> & operator= (const Vector3<C> & r) { v[0]=r[0]; v[1]=r[1]; v[2]=r[2]; return *this; }
	Vector3<C> & operator+= (const Vector3<C> & r) { return (*this) = (*this) + r; }
	Vector3<C> & operator-= (const Vector3<C> & r) { return (*this) = (*this) - r; }
	Vector3<C> & operator*= (const C & r) { return (*this) = (*this) * r; }
	Vector3<C> & operator/= (const C & r) { return (*this) = (*this) / r; }

	bool operator== (const Vector3<C> & r) { return (v[0]==r[0] && v[1]==r[1] && v[2]==r[2]); }

	C Dot(const Vector3<C> & r) const { return v[0]*r[0] + v[1]*r[1] + v[2]*r[2]; }
	Vector3<C> Cross(const Vector3<C> & r) const {
		return Vector3<C>(	v[1] * r[2] - v[2] * r[1],
							v[2] * r[0] - v[0] * r[2],
							v[0] * r[1] - v[1] * r[0] );
	}

	Vector3d unitcross(const Vector3d& r) const {Vector3d vr(v[1] * r.v[2] - v[2] * r.v[1],
											 v[2] * r.v[0] - v[0] * r.v[2],
											 v[0] * r.v[1] - v[1] * r.v[0]); vr.normalize(); return vr;};

	void normalize() {double a = double(sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])); if (a!=0.0) {v[0]/=a; v[1]/=a; v[2]/=a;};};
	C L1Norm() const { return fabs(v[0]) + fabs(v[1]) + fabs(v[2]); }
	C L2Norm() const { return sqrt(Dot(*this)); }
	C length() {return C(sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]));}
	C length() const {return C(sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]));}
	C SquareLength() const { return Dot(*this); }
	C Distance(const Vector3<C> & r) const { return (*this-r).L2Norm(); }
	double angle(const Vector3<C> & r) const {
		double ang = acos((*this).Dot(r)/((*this).length()*r.length()));
		if (isnan(ang))
			return 0;
		return acos((*this).Dot(r)/((*this).length()*r.length()));
	}

	Vector3<C> Min(const Vector3<C> & r) {
		return Vector3<C>(__min(v[0], r[0]), __min(v[1], r[1]), __min(v[2], r[2]));
	}
	Vector3<C> Max(const Vector3<C> & r) {
		return Vector3<C>(__max(v[0], r[0]), __max(v[1], r[1]), __max(v[2], r[2]));
	}

	const C * ToArray() const { return v; }
	C * ToArray() { return v; }

	friend Vector3<C> operator*(const C & l, const Vector3<C> & r) {
		return Vector3<C>(l*r[0], l*r[1], l*r[2]);
	}
	friend ostream & operator<< (ostream & out, const Vector3<C> & r) {
		return out << r[0] << " " << r[1] << " " << r[2];
	}
};

#endif // __VECTOR3_H__
