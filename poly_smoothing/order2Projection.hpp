#include <sys/time.h>
#include <fstream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <cmath>

void convertToOrder2Poly(const double* orgP, double* newP) {
	newP[0] = orgP[0]*orgP[0] + orgP[1]*orgP[1] + orgP[2]*orgP[2];
	newP[1] = orgP[0];
	newP[2] = orgP[1];
	newP[3] = orgP[2];
	newP[4] = 1;
}

void getJacM(const double* orgP, double jacM[][3]) {
	jacM[0][0] = 2*orgP[1];
	jacM[0][1] = 2*orgP[2];
	jacM[0][2] = 2*orgP[3];
	jacM[1][0] = 1;
	jacM[2][1] = 1;
	jacM[3][2] = 1;
}

bool isLocallyFlat(const vector<double*>& points, int pointSize, const Vector3d& cp, const Vector3d& normal) {
	double myPi = 3.14159265358979323846;
	for (int i = 0; i < pointSize; i++) {
		Vector3d p1(points[i]);
		Vector3d dv = p1 - cp;
		double angle = normal.angle(dv);
		if (angle < (myPi/2 - myPi/8) || angle > (myPi/2 + myPi/8))
			return false;
	}
	return true;
}

void solvePCA(const vector<double*>& tmpPoint, int nvsize, vector<double>& weight, double totalWeight, gsl_vector* eval, gsl_matrix* evec, int dim) {
	double centroid[dim];
	for (int d = 0; d < dim; d++)	{
		centroid[d] = 0;
		for (int i = 0; i < nvsize; i++) {
			centroid[d] = centroid[d] + weight[i]*tmpPoint[i][d];
		}
		centroid[d] = centroid[d]/totalWeight;
	}
	double cov[dim][dim];
	for (int d1 = 0; d1 < dim; d1++) {
		for (int d2 = 0; d2 < dim; d2++) {
			cov[d1][d2] = 0;
			for (int i = 0; i < nvsize; i++) {
				cov[d1][d2] = cov[d1][d2] + weight[i]*(tmpPoint[i][d1] - centroid[d1])*(tmpPoint[i][d2] - centroid[d2]);
			}
		}
	}
	for (int d1 = 0; d1 < dim; d1++) {
		for (int d2 = 0; d2 < dim; d2++) {
			cov[d1][d2] /= totalWeight;
		}
	}

	gsl_matrix * m = gsl_matrix_alloc(dim, dim);
	for (int d1 = 0; d1 < dim; d1++) {
		for (int d2 = 0; d2 < dim; d2++) {
			gsl_matrix_set (m, d1, d2, cov[d1][d2]);
		}
	}


	gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (dim);
	gsl_eigen_symmv (m, eval, evec, w);
	gsl_eigen_symmv_free (w);
	gsl_matrix_free(m);
}

void solveAD(const vector<double*>& tmpPoint, int nvsize, vector<double>& weight, double totalWeight, gsl_vector* eval, gsl_matrix* evec, int dim) {
	double mm[dim][dim];
	double nn[dim][dim];
	for (int d1 = 0; d1 < dim; d1++) {
		for (int d2 = 0; d2 < dim; d2++) {
			mm[d1][d2] = 0;
			nn[d1][d2] = 0;
		}
	}
	for (int i = 0; i < nvsize; i++) {
		double cJac[dim][3];
		for (int d1 = 0; d1 < dim; d1++) {
			for (int d2 = 0; d2 < 3; d2++) {
				cJac[d1][d2] = 0;
			}
		}
		getJacM(tmpPoint[i],cJac);
		for (int d1 = 0; d1 < dim; d1++) {
			for (int d2 = 0; d2 < dim; d2++) {
				mm[d1][d2] += (weight[i]*(tmpPoint[i][d1]*tmpPoint[i][d2]));
				if (d1 == dim-1 || d2 == dim-1)
					nn[d1][d2] += weight[i]*0.0000001;
				else {
					for (int d = 0; d < 3; d++) {
						nn[d1][d2] += (weight[i]*(cJac[d1][d]*cJac[d2][d]));
					}
				}
			}
		}
	}
	for (int d1 = 0; d1 < dim; d1++) {
		for (int d2 = 0; d2 < dim; d2++) {
			mm[d1][d2] /= (totalWeight);
			nn[d1][d2] /= (totalWeight);
		}
	}

	gsl_matrix * mp = gsl_matrix_alloc(dim, dim);
	for (int d1 = 0; d1 < dim; d1++) {
		for (int d2 = 0; d2 < dim; d2++) {
			gsl_matrix_set (mp, d1, d2, mm[d1][d2]);
		}
	}
	gsl_matrix * np = gsl_matrix_alloc(dim, dim);
	for (int d1 = 0; d1 < dim; d1++) {
		for (int d2 = 0; d2 < dim; d2++) {
			gsl_matrix_set (np, d1, d2, nn[d1][d2]);
		}
	}
	gsl_eigen_gensymmv_workspace * w = gsl_eigen_gensymmv_alloc (dim);
	gsl_eigen_gensymmv (mp, np, eval, evec, w);
	gsl_eigen_gensymmv_free (w);
	gsl_matrix_free(mp);
	gsl_matrix_free(np);
}

void computePolyFunction(const vector<const double*>& neighborPoints, const vector<double*>& neighborPolyPoints, const double* refPoint, double sideLength, int dim, double* polyFunc, double& polyb) {
	vector<double> weight;
	double totalWeight = 0;
	int neighborSize = neighborPoints.size();
	for (int i = 0; i < neighborSize; i++) {
		double cd = 0;
		for (int d = 0; d < 3; d++)
			cd += ((neighborPoints[i][d] - refPoint[d])*(neighborPoints[i][d] - refPoint[d]));
		double cw = exp(-1*(cd*cd)/(2*sideLength*sideLength));
		weight.push_back(cw);
		totalWeight += cw;
	}

	gsl_vector *eval = gsl_vector_alloc (dim);
	gsl_matrix *evec = gsl_matrix_alloc (dim, dim);
	solveAD(neighborPolyPoints, neighborSize, weight, totalWeight, eval, evec, dim);

	double minEigenValue = 100;
	int minIndex = 0;
	for (int d = 0; d < dim; d++) {
		if (d == 0 || fabs(gsl_vector_get (eval, d)) < minEigenValue) {
			minEigenValue = fabs(gsl_vector_get (eval, d));
			minIndex = d;
		}
	}

	gsl_vector *smallestEigenVector = gsl_vector_alloc(dim);
	gsl_matrix_get_col (smallestEigenVector, evec, minIndex);
	for (int d = 0; d < dim-1; d++) {
		polyFunc[d] = gsl_vector_get(smallestEigenVector, d);
	}
	polyb = gsl_vector_get(smallestEigenVector, dim-1);

	gsl_vector_free(eval);
	gsl_matrix_free(evec);
}

void normalOfOrder2Poly(double* polyFunc, double* p, double* normal) {
	normal[0] = 2*polyFunc[0]*p[0] + polyFunc[1];
	normal[1] = 2*polyFunc[0]*p[1] + polyFunc[2];
	normal[2] = 2*polyFunc[0]*p[2] + polyFunc[3];
}

bool projectToOrder2Poly(double* polyFunc, double b, double* normal, double* orgP, double* projectedP) {
	double A = polyFunc[0]*(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
	double B = 2*polyFunc[0]*(orgP[0]*normal[0] + orgP[1]*normal[1] + orgP[2]*normal[2]) + polyFunc[1]*normal[0] + polyFunc[2]*normal[1] + polyFunc[3]*normal[2];
	double C = polyFunc[0]*(orgP[0]*orgP[0] + orgP[1]*orgP[1] + orgP[2]*orgP[2]) + polyFunc[1]*orgP[0] + polyFunc[2]*orgP[1] + polyFunc[3]*orgP[2] + b;

	double discrim = B*B - 4*A*C;

	double ip1[3];
	double ip2[3];

	if (discrim >= 0) {
		double sol1 = (-B + sqrt(discrim))/(2*A);
		double sol2 = (-B - sqrt(discrim))/(2*A);

		for (int d = 0; d < 3; d++) {
			ip1[d] = orgP[d] + sol1*normal[d];
			ip2[d] = orgP[d] + sol2*normal[d];
		}
		double dist1 = 0;
		double dist2 = 0;
		for (int d = 0; d < 3; d++) {
			dist1 += ((orgP[d] - ip1[d])*(orgP[d] - ip1[d]));
			dist2 += ((orgP[d] - ip2[d])*(orgP[d] - ip2[d]));
		}

		if (fabs(dist2) < fabs(dist1)) {
			for (int d = 0; d < 3; d++)
				projectedP[d] = ip2[d];
		}
		else {
			for (int d = 0; d < 3; d++)
				projectedP[d] = ip1[d];
		}
		return true;
	}
	return false;
}
