#include "grid.hpp"
#include <cmath>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>

using namespace std;

grid::grid(octreeNode* _node, const double* _center, double _sideLength) {
	workingNode = _node;
	sideLength = _sideLength;
	for (int d = 0; d < 3; d++)
	{
		center[d] = _center[d];
		minCoord[d] = center[d] - 1.5*sideLength;
		plane[d] = 0;
	}
	b = 0;
	octreeNode* const * workingNeighbor = workingNode->getSameLevelNeighbors();
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				double tmpCenter[3];
				tmpCenter[0] = center[0] + (i - 1)*sideLength;
				tmpCenter[1] = center[1] + (j - 1)*sideLength;
				tmpCenter[2] = center[2] + (k - 1)*sideLength;
				gc[i][j][k].setCenter(tmpCenter);
				gc[i][j][k].setSideLength(sideLength);
			}
		}
	}
	gc[2][2][2].setWorkingNode(workingNeighbor[0]);
	gc[1][2][2].setWorkingNode(workingNeighbor[1]);
	gc[0][2][2].setWorkingNode(workingNeighbor[2]);
	gc[2][1][2].setWorkingNode(workingNeighbor[3]);
	gc[1][1][2].setWorkingNode(workingNeighbor[4]);
	gc[0][1][2].setWorkingNode(workingNeighbor[5]);
	gc[2][0][2].setWorkingNode(workingNeighbor[6]);
	gc[1][0][2].setWorkingNode(workingNeighbor[7]);
	gc[0][0][2].setWorkingNode(workingNeighbor[8]);
	gc[2][2][1].setWorkingNode(workingNeighbor[9]);
	gc[1][2][1].setWorkingNode(workingNeighbor[10]);
	gc[0][2][1].setWorkingNode(workingNeighbor[11]);
	gc[2][1][1].setWorkingNode(workingNeighbor[12]);
	gc[0][1][1].setWorkingNode(workingNeighbor[13]);
	gc[2][0][1].setWorkingNode(workingNeighbor[14]);
	gc[1][0][1].setWorkingNode(workingNeighbor[15]);
	gc[0][0][1].setWorkingNode(workingNeighbor[16]);
	gc[2][2][0].setWorkingNode(workingNeighbor[17]);
	gc[1][2][0].setWorkingNode(workingNeighbor[18]);
	gc[0][2][0].setWorkingNode(workingNeighbor[19]);
	gc[2][1][0].setWorkingNode(workingNeighbor[20]);
	gc[1][1][0].setWorkingNode(workingNeighbor[21]);
	gc[0][1][0].setWorkingNode(workingNeighbor[22]);
	gc[2][0][0].setWorkingNode(workingNeighbor[23]);
	gc[1][0][0].setWorkingNode(workingNeighbor[24]);
	gc[0][0][0].setWorkingNode(workingNeighbor[25]);
	gc[1][1][1].setWorkingNode(workingNode);
}

octreeNode* const grid::getNode() const {
	return workingNode;
}

double grid::getSideLength() const {
	return sideLength;
}

bool grid::computeNormal(const vector<double*> &v, vector<double*> &normal, double* myb) {
	double myPi = 3.14159265358979323846;
	vector<const double*> tmpPoint;
	const vector<int> & cv = gc[2][2][2].getInCellVertex();
	int cvsize = gc[2][2][2].getNumOfVertex();
	for (int i = 0; i < 3; i++)	{
		for (int j = 0; j < 3; j++)	{
			for (int k = 0; k < 3; k++)	{
				if (gc[i][j][k].getWorkingNode() != NULL) {
					int nvsize = gc[i][j][k].getNumOfVertex();
					if (nvsize > 0) {
						tmpPoint.push_back(gc[i][j][k].getAvgPoint());
					}
				}
			}
		}
	}
	int nvsize = tmpPoint.size();
	if (nvsize < 5)
		return false;
	double centroid[3];
	for (int d = 0; d < 3; d++)	{
		centroid[d] = 0;
		for (int i = 0; i < nvsize; i++) {
			centroid[d] = centroid[d] + tmpPoint[i][d];
		}
		centroid[d] = centroid[d]/nvsize;
	}
	double cov[3][3];
	for (int d1 = 0; d1 < 3; d1++) {
		for (int d2 = 0; d2 < 3; d2++) {
			cov[d1][d2] = 0;
			for (int i = 0; i < nvsize; i++) {
				cov[d1][d2] = cov[d1][d2] + (tmpPoint[i][d1] - centroid[d1])*(tmpPoint[i][d2] - centroid[d2]);
			}
		}
	}
	gsl_matrix * m = gsl_matrix_alloc(3, 3);
	for (int d1 = 0; d1 < 3; d1++) {
		for (int d2 = 0; d2 < 3; d2++) {
			gsl_matrix_set (m, d1, d2, cov[d1][d2]);
		}
	}

	gsl_vector *eval = gsl_vector_alloc (3);
	gsl_matrix *evec = gsl_matrix_alloc (3, 3);
	gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (3);
	gsl_eigen_symmv (m, eval, evec, w);
	gsl_eigen_symmv_free (w);
	double minEigenValue = 100;
	int minIndex = 0;
	double maxEigenValue = -1;
	int maxIndex = 0;
	double midEigenValue = 100;
	int midIndex = 0;
	for (int d = 0; d < 3; d++) {
		if (d == 0 || gsl_vector_get (eval, d) < minEigenValue) {
			minEigenValue = gsl_vector_get (eval, d);
			minIndex = d;
		}
		if (d == 0 || gsl_vector_get (eval, d) > maxEigenValue) {
			maxEigenValue = gsl_vector_get (eval, d);
			maxIndex = d;
		}
	}
	for (int d = 0; d < 3; d++) {
		if (d != minIndex && d != maxIndex) {
			midEigenValue = gsl_vector_get (eval, d);
			midIndex = d;
		}
	}

	if (fabs(maxEigenValue)/fabs(minEigenValue) < 2) {
		return false;
	}

	gsl_vector *smallestEigenVector = gsl_vector_alloc(3);
	gsl_matrix_get_col (smallestEigenVector, evec, minIndex);
	b = 0;
	for (int d = 0; d < 3; d++)
	{
		plane[d] = gsl_vector_get(smallestEigenVector, d);
	}
	for (int i = 0; i < nvsize; i++)
	{
		for (int d = 0; d < 3; d++)
			b = b + plane[d]*tmpPoint[i][d];
	}
	b = b/nvsize;
	b = -1*b;
	for (int cvi = 0; cvi < cvsize; cvi++)
	{
		for (int d = 0; d < 3; d++)
			normal[cv[cvi]][d] = plane[d];
		myb[cv[cvi]] = b;
	}
	gsl_vector_free(eval);
	gsl_matrix_free(evec);
	gsl_matrix_free(m);
	return true;

}

bool grid::hasEmptySmallCell(const vector<double*> &v, const double* rMaxC, const double* rMinC)
{
	double innerLength = sideLength/4;
	double halfLength = sideLength/2;
	double longestDist = sqrt(3)*halfLength;
	for (int i = 0; i < 3; i++)	{
		for (int j = 0; j < 3; j++)	{
			for (int k = 0; k < 3; k++)	{
				if (gc[i][j][k].getNumOfVertex() == 0) {
					const double* nCenter = gc[i][j][k].getCenter();
					double ppd = fabs(plane[0]*nCenter[0] + plane[1]*nCenter[1] + plane[2]*nCenter[2] + b)/sqrt(plane[0]*plane[0] + plane[1]*plane[1] + plane[2]*plane[2]);
					if (ppd < longestDist + 0.00000001) {
						if (gc[i][j][k].setEdgesAndCheckIntersection(innerLength, halfLength, plane, b, gc, i, j, k, rMaxC, rMinC)) {
							return true;
						}
					}
				}
			}
		}
	}
	return false;
}

const double* grid::getNormal() const
{
	return plane;
}

double grid::getb() const
{
	return b;
}

gridCell::gridCell()
{
	numOfVertex = 0;
}

void gridCell::setCenter(double* _center)
{
	for (int d = 0; d < 3; d++)
	{
		center[d] = _center[d];
	}
}

const double* gridCell::getCenter() const
{
	return center;
}

void gridCell::setSideLength(double _sideLength)
{
	sideLength = _sideLength;
}

double gridCell::getSideLength()
{
	return sideLength;
}

//---- Set the edge by considering all case ----
bool gridCell::setEdgesAndCheckIntersection(double innerLength, double halfLength, const double* plane, const double b, gridCell gc[3][3][3], int cx, int cy, int cz, const double* rMaxC, const double* rMinC)
{
	double edges[144][2][3];
	int edgeCells[144][7][3];
	int edgeNeighborCellSize[144];

	int sign1 = 0;
	int sign2 = 0;
	// xy
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			for (int k = 0; k < 2; k++)
			{
				if (i == 0)
				{
					sign1 = -1;
					sign2 = -1;
				}
				else if (i == 1)
				{
					sign1 = 1;
					sign2 = -1;
				}
				else if (i == 2)
				{
					sign1 = -1;
					sign2 = 1;
				}
				else if (i == 3)
				{
					sign1 = 1;
					sign2 = 1;
				}

				edges[3*(i*4 + j*2 + k) + 0][0][0] = center[0] + sign1*halfLength + -1*sign1*k*innerLength;
				edges[3*(i*4 + j*2 + k) + 0][0][1] = center[1] + sign2*halfLength + -1*sign2*j*innerLength;
				edges[3*(i*4 + j*2 + k) + 0][0][2] = center[2] + halfLength;
				edges[3*(i*4 + j*2 + k) + 0][1][0] = center[0] + sign1*halfLength + -1*sign1*k*innerLength;
				edges[3*(i*4 + j*2 + k) + 0][1][1] = center[1] + sign2*halfLength + -1*sign2*j*innerLength;
				edges[3*(i*4 + j*2 + k) + 0][1][2] = center[2] + halfLength - innerLength;
				if (j == 0 && k == 0)
				{
					edgeCells[3*(i*4 + j*2 + k) + 0][0][0] = 0;
					edgeCells[3*(i*4 + j*2 + k) + 0][0][1] = 0;
					edgeCells[3*(i*4 + j*2 + k) + 0][0][2] = 1;

					edgeCells[3*(i*4 + j*2 + k) + 0][1][0] = sign1;
					edgeCells[3*(i*4 + j*2 + k) + 0][1][1] = 0;
					edgeCells[3*(i*4 + j*2 + k) + 0][1][2] = 1;

					edgeCells[3*(i*4 + j*2 + k) + 0][2][0] = 0;
					edgeCells[3*(i*4 + j*2 + k) + 0][2][1] = sign2;
					edgeCells[3*(i*4 + j*2 + k) + 0][2][2] = 1;

					edgeCells[3*(i*4 + j*2 + k) + 0][3][0] = sign1;
					edgeCells[3*(i*4 + j*2 + k) + 0][3][1] = sign2;
					edgeCells[3*(i*4 + j*2 + k) + 0][3][2] = 1;

					edgeCells[3*(i*4 + j*2 + k) + 0][4][0] = sign1;
					edgeCells[3*(i*4 + j*2 + k) + 0][4][1] = 0;
					edgeCells[3*(i*4 + j*2 + k) + 0][4][2] = 0;

					edgeCells[3*(i*4 + j*2 + k) + 0][5][0] = 0;
					edgeCells[3*(i*4 + j*2 + k) + 0][5][1] = sign1;
					edgeCells[3*(i*4 + j*2 + k) + 0][5][2] = 0;

					edgeCells[3*(i*4 + j*2 + k) + 0][6][0] = sign1;
					edgeCells[3*(i*4 + j*2 + k) + 0][6][1] = sign2;
					edgeCells[3*(i*4 + j*2 + k) + 0][6][2] = 0;

					edgeNeighborCellSize[3*(i*4 + j*2 + k) + 0] = 7;
				}
				else if (j == 0 && k == 1)
				{
					edgeCells[3*(i*4 + j*2 + k) + 0][0][0] = 0;
					edgeCells[3*(i*4 + j*2 + k) + 0][0][1] = 0;
					edgeCells[3*(i*4 + j*2 + k) + 0][0][2] = 1;

					edgeCells[3*(i*4 + j*2 + k) + 0][1][0] = 0;
					edgeCells[3*(i*4 + j*2 + k) + 0][1][1] = sign2;
					edgeCells[3*(i*4 + j*2 + k) + 0][1][2] = 1;

					edgeCells[3*(i*4 + j*2 + k) + 0][2][0] = 0;
					edgeCells[3*(i*4 + j*2 + k) + 0][2][1] = sign2;
					edgeCells[3*(i*4 + j*2 + k) + 0][2][2] = 0;

					edgeNeighborCellSize[3*(i*4 + j*2 + k) + 0] = 3;
				}
				else if (j == 1 && k == 0)
				{
					edgeCells[3*(i*4 + j*2 + k) + 0][0][0] = 0;
					edgeCells[3*(i*4 + j*2 + k) + 0][0][1] = 0;
					edgeCells[3*(i*4 + j*2 + k) + 0][0][2] = 1;

					edgeCells[3*(i*4 + j*2 + k) + 0][1][0] = sign1;
					edgeCells[3*(i*4 + j*2 + k) + 0][1][1] = 0;
					edgeCells[3*(i*4 + j*2 + k) + 0][1][2] = 1;

					edgeCells[3*(i*4 + j*2 + k) + 0][2][0] = sign1;
					edgeCells[3*(i*4 + j*2 + k) + 0][2][1] = 0;
					edgeCells[3*(i*4 + j*2 + k) + 0][2][2] = 0;

					edgeNeighborCellSize[3*(i*4 + j*2 + k) + 0] = 3;
				}
				else if (j == 1 && k == 1)
				{
					edgeCells[3*(i*4 + j*2 + k) + 0][0][0] = 0;
					edgeCells[3*(i*4 + j*2 + k) + 0][0][1] = 0;
					edgeCells[3*(i*4 + j*2 + k) + 0][0][2] = 1;

					edgeNeighborCellSize[3*(i*4 + j*2 + k) + 0] = 1;
				}

				edges[3*(i*4 + j*2 + k) + 1][0][0] = center[0] + sign1*halfLength + -1*sign1*k*innerLength;
				edges[3*(i*4 + j*2 + k) + 1][0][1] = center[1] + sign2*halfLength + -1*sign2*j*innerLength;
				edges[3*(i*4 + j*2 + k) + 1][0][2] = center[2] + halfLength - innerLength;
				edges[3*(i*4 + j*2 + k) + 1][1][0] = center[0] + sign1*halfLength + -1*sign1*k*innerLength;
				edges[3*(i*4 + j*2 + k) + 1][1][1] = center[1] + sign2*halfLength + -1*sign2*j*innerLength;
				edges[3*(i*4 + j*2 + k) + 1][1][2] = center[2] - halfLength + innerLength;

				if (j == 0 && k == 0)
				{
					edgeCells[3*(i*4 + j*2 + k) + 1][0][0] = sign1;
					edgeCells[3*(i*4 + j*2 + k) + 1][0][1] = 0;
					edgeCells[3*(i*4 + j*2 + k) + 1][0][2] = 0;

					edgeCells[3*(i*4 + j*2 + k) + 1][1][0] = 0;
					edgeCells[3*(i*4 + j*2 + k) + 1][1][1] = sign2;
					edgeCells[3*(i*4 + j*2 + k) + 1][1][2] = 0;

					edgeCells[3*(i*4 + j*2 + k) + 1][2][0] = sign1;
					edgeCells[3*(i*4 + j*2 + k) + 1][2][1] = sign2;
					edgeCells[3*(i*4 + j*2 + k) + 1][2][2] = 0;

					edgeNeighborCellSize[3*(i*4 + j*2 + k) + 1] = 3;
				}
				else if (j == 0 && k == 1)
				{
					edgeCells[3*(i*4 + j*2 + k) + 1][0][0] = 0;
					edgeCells[3*(i*4 + j*2 + k) + 1][0][1] = sign2;
					edgeCells[3*(i*4 + j*2 + k) + 1][0][2] = 0;

					edgeNeighborCellSize[3*(i*4 + j*2 + k) + 1] = 1;
				}
				else if (j == 1 && k == 0)
				{
					edgeCells[3*(i*4 + j*2 + k) + 1][0][0] = sign1;
					edgeCells[3*(i*4 + j*2 + k) + 1][0][1] = 0;
					edgeCells[3*(i*4 + j*2 + k) + 1][0][2] = 0;

					edgeNeighborCellSize[3*(i*4 + j*2 + k) + 1] = 1;
				}
				else if (j == 1 && k == 1)
				{
					edgeNeighborCellSize[3*(i*4 + j*2 + k) + 1] = 0;
				}

				edges[3*(i*4 + j*2 + k) + 2][0][0] = center[0] + sign1*halfLength + -1*sign1*k*innerLength;
				edges[3*(i*4 + j*2 + k) + 2][0][1] = center[1] + sign2*halfLength + -1*sign2*j*innerLength;
				edges[3*(i*4 + j*2 + k) + 2][0][2] = center[2] - halfLength + innerLength;
				edges[3*(i*4 + j*2 + k) + 2][1][0] = center[0] + sign1*halfLength + -1*sign1*k*innerLength;
				edges[3*(i*4 + j*2 + k) + 2][1][1] = center[1] + sign2*halfLength + -1*sign2*j*innerLength;
				edges[3*(i*4 + j*2 + k) + 2][1][2] = center[2] - halfLength;
				edgeCells[3*(i*4 + j*2 + k) + 2][0][2] = -1;
				if (j == 0 && k == 0)
				{
					edgeCells[3*(i*4 + j*2 + k) + 2][0][0] = 0;
					edgeCells[3*(i*4 + j*2 + k) + 2][0][1] = 0;
					edgeCells[3*(i*4 + j*2 + k) + 2][0][2] = -1;

					edgeCells[3*(i*4 + j*2 + k) + 2][1][0] = sign1;
					edgeCells[3*(i*4 + j*2 + k) + 2][1][1] = 0;
					edgeCells[3*(i*4 + j*2 + k) + 2][1][2] = -1;

					edgeCells[3*(i*4 + j*2 + k) + 2][2][0] = 0;
					edgeCells[3*(i*4 + j*2 + k) + 2][2][1] = sign2;
					edgeCells[3*(i*4 + j*2 + k) + 2][2][2] = -1;

					edgeCells[3*(i*4 + j*2 + k) + 2][3][0] = sign1;
					edgeCells[3*(i*4 + j*2 + k) + 2][3][1] = sign2;
					edgeCells[3*(i*4 + j*2 + k) + 2][3][2] = -1;

					edgeCells[3*(i*4 + j*2 + k) + 2][4][0] = sign1;
					edgeCells[3*(i*4 + j*2 + k) + 2][4][1] = 0;
					edgeCells[3*(i*4 + j*2 + k) + 2][4][2] = 0;

					edgeCells[3*(i*4 + j*2 + k) + 2][5][0] = 0;
					edgeCells[3*(i*4 + j*2 + k) + 2][5][1] = sign2;
					edgeCells[3*(i*4 + j*2 + k) + 2][5][2] = 0;

					edgeCells[3*(i*4 + j*2 + k) + 2][6][0] = sign1;
					edgeCells[3*(i*4 + j*2 + k) + 2][6][1] = sign2;
					edgeCells[3*(i*4 + j*2 + k) + 2][6][2] = 0;

					edgeNeighborCellSize[3*(i*4 + j*2 + k) + 2] = 7;
				}
				else if (j == 0 && k == 1)
				{
					edgeCells[3*(i*4 + j*2 + k) + 2][0][0] = 0;
					edgeCells[3*(i*4 + j*2 + k) + 2][0][1] = 0;
					edgeCells[3*(i*4 + j*2 + k) + 2][0][2] = -1;

					edgeCells[3*(i*4 + j*2 + k) + 2][1][0] = 0;
					edgeCells[3*(i*4 + j*2 + k) + 2][1][1] = sign2;
					edgeCells[3*(i*4 + j*2 + k) + 2][1][2] = -1;

					edgeCells[3*(i*4 + j*2 + k) + 2][2][0] = 0;
					edgeCells[3*(i*4 + j*2 + k) + 2][2][1] = sign2;
					edgeCells[3*(i*4 + j*2 + k) + 2][2][2] = 0;

					edgeNeighborCellSize[3*(i*4 + j*2 + k) + 2] = 3;
				}
				else if (j == 1 && k == 0)
				{
					edgeCells[3*(i*4 + j*2 + k) + 2][0][0] = 0;
					edgeCells[3*(i*4 + j*2 + k) + 2][0][1] = 0;
					edgeCells[3*(i*4 + j*2 + k) + 2][0][2] = -1;

					edgeCells[3*(i*4 + j*2 + k) + 2][1][0] = sign1;
					edgeCells[3*(i*4 + j*2 + k) + 2][1][1] = 0;
					edgeCells[3*(i*4 + j*2 + k) + 2][1][2] = -1;

					edgeCells[3*(i*4 + j*2 + k) + 2][2][0] = sign1;
					edgeCells[3*(i*4 + j*2 + k) + 2][2][1] = 0;
					edgeCells[3*(i*4 + j*2 + k) + 2][2][2] = 0;

					edgeNeighborCellSize[3*(i*4 + j*2 + k) + 2] = 3;
				}
				else if (j == 1 && k == 1)
				{
					edgeCells[3*(i*4 + j*2 + k) + 2][0][0] = 0;
					edgeCells[3*(i*4 + j*2 + k) + 2][0][1] = 0;
					edgeCells[3*(i*4 + j*2 + k) + 2][0][2] = -1;

					edgeNeighborCellSize[3*(i*4 + j*2 + k) + 2] = 1;
				}
			}
		}
	}
	//xz
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			for (int k = 0; k < 2; k++)
			{
				if (i == 0)
				{
					sign1 = -1;
					sign2 = -1;
				}
				else if (i == 1)
				{
					sign1 = 1;
					sign2 = -1;
				}
				else if (i == 2)
				{
					sign1 = -1;
					sign2 = 1;
				}
				else if (i == 3)
				{
					sign1 = 1;
					sign2 = 1;
				}

				edges[48 + 3*(i*4 + j*2 + k) + 0][0][0] = center[0] + sign1*halfLength + -1*sign1*k*innerLength;
				edges[48 + 3*(i*4 + j*2 + k) + 0][0][1] = center[1] + halfLength;
				edges[48 + 3*(i*4 + j*2 + k) + 0][0][2] = center[2] + sign2*halfLength + -1*sign2*j*innerLength;
				edges[48 + 3*(i*4 + j*2 + k) + 0][1][0] = center[0] + sign1*halfLength + -1*sign1*k*innerLength;
				edges[48 + 3*(i*4 + j*2 + k) + 0][1][1] = center[1] + halfLength - innerLength;
				edges[48 + 3*(i*4 + j*2 + k) + 0][1][2] = center[2] + sign2*halfLength + -1*sign2*j*innerLength;
				if (j == 0 && k == 0)
				{
					edgeCells[48 + 3*(i*4 + j*2 + k) + 0][0][0] = 0;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 0][0][1] = 1;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 0][0][2] = 0;

					edgeCells[48 + 3*(i*4 + j*2 + k) + 0][1][0] = sign1;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 0][1][1] = 1;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 0][1][2] = 0;

					edgeCells[48 + 3*(i*4 + j*2 + k) + 0][2][0] = 0;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 0][2][1] = 1;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 0][2][2] = sign2;

					edgeCells[48 + 3*(i*4 + j*2 + k) + 0][3][0] = sign1;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 0][3][1] = 1;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 0][3][2] = sign2;

					edgeCells[48 + 3*(i*4 + j*2 + k) + 0][4][0] = sign1;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 0][4][1] = 0;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 0][4][2] = 0;

					edgeCells[48 + 3*(i*4 + j*2 + k) + 0][5][0] = 0;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 0][5][1] = 0;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 0][5][2] = sign2;

					edgeCells[48 + 3*(i*4 + j*2 + k) + 0][6][0] = sign1;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 0][6][1] = 0;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 0][6][2] = sign2;


					edgeNeighborCellSize[48 + 3*(i*4 + j*2 + k) + 0] = 7;
				}
				else if (j == 0 && k == 1)
				{
					edgeCells[48 + 3*(i*4 + j*2 + k) + 0][0][0] = 0;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 0][0][1] = 1;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 0][0][2] = 0;

					edgeCells[48 + 3*(i*4 + j*2 + k) + 0][1][0] = 0;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 0][1][1] = 1;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 0][1][2] = sign2;

					edgeCells[48 + 3*(i*4 + j*2 + k) + 0][2][0] = 0;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 0][2][1] = 0;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 0][2][2] = sign2;


					edgeNeighborCellSize[48 + 3*(i*4 + j*2 + k) + 0] = 3;
				}
				else if (j == 1 && k == 0)
				{
					edgeCells[48 + 3*(i*4 + j*2 + k) + 0][0][0] = 0;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 0][0][1] = 1;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 0][0][2] = 0;

					edgeCells[48 + 3*(i*4 + j*2 + k) + 0][1][0] = sign1;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 0][1][1] = 1;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 0][1][2] = 0;

					edgeCells[48 + 3*(i*4 + j*2 + k) + 0][2][0] = sign1;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 0][2][1] = 0;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 0][2][2] = 0;


					edgeNeighborCellSize[48 + 3*(i*4 + j*2 + k) + 0] = 3;
				}
				else if (j == 1 && k == 1)
				{
					edgeCells[48 + 3*(i*4 + j*2 + k) + 0][0][0] = 0;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 0][0][1] = 1;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 0][0][2] = 0;

					edgeNeighborCellSize[48 + 3*(i*4 + j*2 + k) + 0] = 1;
				}

				edges[48 + 3*(i*4 + j*2 + k) + 1][0][0] = center[0] + sign1*halfLength + -1*sign1*k*innerLength;
				edges[48 + 3*(i*4 + j*2 + k) + 1][0][1] = center[1] + halfLength - innerLength;
				edges[48 + 3*(i*4 + j*2 + k) + 1][0][2] = center[2] + sign2*halfLength + -1*sign2*j*innerLength;
				edges[48 + 3*(i*4 + j*2 + k) + 1][1][0] = center[0] + sign1*halfLength + -1*sign1*k*innerLength;
				edges[48 + 3*(i*4 + j*2 + k) + 1][1][1] = center[1] - halfLength + innerLength;
				edges[48 + 3*(i*4 + j*2 + k) + 1][1][2] = center[2] + sign2*halfLength + -1*sign2*j*innerLength;

				if (j == 0 && k == 0)
				{
					edgeCells[48 + 3*(i*4 + j*2 + k) + 1][0][0] = sign1;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 1][0][1] = 0;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 1][0][2] = 0;

					edgeCells[48 + 3*(i*4 + j*2 + k) + 1][1][0] = 0;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 1][1][1] = 0;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 1][1][2] = sign2;

					edgeCells[48 + 3*(i*4 + j*2 + k) + 1][2][0] = sign1;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 1][2][1] = 0;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 1][2][2] = sign2;


					edgeNeighborCellSize[48 + 3*(i*4 + j*2 + k) + 1] = 3;
				}
				else if (j == 0 && k == 1)
				{
					edgeCells[48 + 3*(i*4 + j*2 + k) + 1][0][0] = 0;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 1][0][1] = 0;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 1][0][2] = sign2;

					edgeNeighborCellSize[48 + 3*(i*4 + j*2 + k) + 1] = 1;
				}
				else if (j == 1 && k == 0)
				{
					edgeCells[48 + 3*(i*4 + j*2 + k) + 1][0][0] = sign1;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 1][0][1] = 0;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 1][0][2] = 0;

					edgeNeighborCellSize[48 + 3*(i*4 + j*2 + k) + 1] = 1;
				}
				else if (j == 1 && k == 1)
				{
					edgeNeighborCellSize[48 + 3*(i*4 + j*2 + k) + 1] = 0;
				}

				edges[48 + 3*(i*4 + j*2 + k) + 2][0][0] = center[0] + sign1*halfLength + -1*sign1*k*innerLength;
				edges[48 + 3*(i*4 + j*2 + k) + 2][0][1] = center[1] - halfLength + innerLength;
				edges[48 + 3*(i*4 + j*2 + k) + 2][0][2] = center[2] + sign2*halfLength + -1*sign2*j*innerLength;
				edges[48 + 3*(i*4 + j*2 + k) + 2][1][0] = center[0] + sign1*halfLength + -1*sign1*k*innerLength;
				edges[48 + 3*(i*4 + j*2 + k) + 2][1][1] = center[1] - halfLength;
				edges[48 + 3*(i*4 + j*2 + k) + 2][1][2] = center[2] + sign2*halfLength + -1*sign2*j*innerLength;
				if (j == 0 && k == 0)
				{
					edgeCells[48 + 3*(i*4 + j*2 + k) + 2][0][0] = 0;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 2][0][1] = -1;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 2][0][2] = 0;

					edgeCells[48 + 3*(i*4 + j*2 + k) + 2][1][0] = sign1;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 2][1][1] = -1;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 2][1][2] = 0;

					edgeCells[48 + 3*(i*4 + j*2 + k) + 2][2][0] = 0;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 2][2][1] = -1;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 2][2][2] = sign2;

					edgeCells[48 + 3*(i*4 + j*2 + k) + 2][3][0] = sign1;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 2][3][1] = -1;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 2][3][2] = sign2;

					edgeCells[48 + 3*(i*4 + j*2 + k) + 2][4][0] = sign1;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 2][4][1] = 0;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 2][4][2] = 0;

					edgeCells[48 + 3*(i*4 + j*2 + k) + 2][5][0] = 0;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 2][5][1] = 0;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 2][5][2] = sign2;

					edgeCells[48 + 3*(i*4 + j*2 + k) + 2][6][0] = sign1;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 2][6][1] = 0;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 2][6][2] = sign2;

					edgeNeighborCellSize[48 + 3*(i*4 + j*2 + k) + 2] = 7;
				}
				else if (j == 0 && k == 1)
				{
					edgeCells[48 + 3*(i*4 + j*2 + k) + 2][0][0] = 0;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 2][0][1] = -1;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 2][0][2] = 0;

					edgeCells[48 + 3*(i*4 + j*2 + k) + 2][1][0] = 0;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 2][1][1] = -1;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 2][1][2] = sign2;

					edgeCells[48 + 3*(i*4 + j*2 + k) + 2][2][0] = 0;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 2][2][1] = 0;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 2][2][2] = sign2;

					edgeNeighborCellSize[48 + 3*(i*4 + j*2 + k) + 2] = 3;
				}
				else if (j == 1 && k == 0)
				{
					edgeCells[48 + 3*(i*4 + j*2 + k) + 2][0][0] = 0;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 2][0][1] = -1;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 2][0][2] = 0;

					edgeCells[48 + 3*(i*4 + j*2 + k) + 2][1][0] = sign1;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 2][1][1] = -1;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 2][1][2] = 0;

					edgeCells[48 + 3*(i*4 + j*2 + k) + 2][2][0] = sign1;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 2][2][1] = 0;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 2][2][2] = 0;

					edgeNeighborCellSize[48 + 3*(i*4 + j*2 + k) + 2] = 3;
				}
				else if (j == 1 && k == 1)
				{
					edgeCells[48 + 3*(i*4 + j*2 + k) + 2][0][0] = 0;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 2][0][1] = -1;
					edgeCells[48 + 3*(i*4 + j*2 + k) + 2][0][2] = 0;

					edgeNeighborCellSize[48 + 3*(i*4 + j*2 + k) + 2] = 1;
				}
			}
		}
	}
	//yz
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			for (int k = 0; k < 2; k++)
			{
				if (i == 0)
				{
					sign1 = -1;
					sign2 = -1;
				}
				else if (i == 1)
				{
					sign1 = 1;
					sign2 = -1;
				}
				else if (i == 2)
				{
					sign1 = -1;
					sign2 = 1;
				}
				else if (i == 3)
				{
					sign1 = 1;
					sign2 = 1;
				}

				edges[96 + 3*(i*4 + j*2 + k) + 0][0][0] = center[0] + halfLength;
				edges[96 + 3*(i*4 + j*2 + k) + 0][0][1] = center[1] + sign1*halfLength + -1*sign1*k*innerLength;
				edges[96 + 3*(i*4 + j*2 + k) + 0][0][2] = center[2] + sign2*halfLength + -1*sign2*j*innerLength;
				edges[96 + 3*(i*4 + j*2 + k) + 0][1][0] = center[0] + halfLength - innerLength;
				edges[96 + 3*(i*4 + j*2 + k) + 0][1][1] = center[1] + sign1*halfLength + -1*sign1*k*innerLength;
				edges[96 + 3*(i*4 + j*2 + k) + 0][1][2] = center[2] + sign2*halfLength + -1*sign2*j*innerLength;

				if (j == 0 && k == 0)
				{
					edgeCells[96 + 3*(i*4 + j*2 + k) + 0][0][0] = 1;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 0][0][1] = 0;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 0][0][2] = 0;

					edgeCells[96 + 3*(i*4 + j*2 + k) + 0][1][0] = 1;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 0][1][1] = sign1;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 0][1][2] = 0;

					edgeCells[96 + 3*(i*4 + j*2 + k) + 0][2][0] = 1;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 0][2][1] = 0;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 0][2][2] = sign2;

					edgeCells[96 + 3*(i*4 + j*2 + k) + 0][3][0] = 1;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 0][3][1] = sign1;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 0][3][2] = sign2;

					edgeCells[96 + 3*(i*4 + j*2 + k) + 0][4][0] = 0;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 0][4][1] = sign1;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 0][4][2] = 0;

					edgeCells[96 + 3*(i*4 + j*2 + k) + 0][5][0] = 0;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 0][5][1] = 0;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 0][5][2] = sign2;

					edgeCells[96 + 3*(i*4 + j*2 + k) + 0][6][0] = 0;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 0][6][1] = sign1;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 0][6][2] = sign2;


					edgeNeighborCellSize[96 + 3*(i*4 + j*2 + k) + 0] = 7;
				}
				else if (j == 0 && k == 1)
				{
					edgeCells[96 + 3*(i*4 + j*2 + k) + 0][0][0] = 1;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 0][0][1] = 0;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 0][0][2] = 0;

					edgeCells[96 + 3*(i*4 + j*2 + k) + 0][1][0] = 1;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 0][1][1] = 0;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 0][1][2] = sign2;

					edgeCells[96 + 3*(i*4 + j*2 + k) + 0][2][0] = 0;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 0][2][1] = 0;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 0][2][2] = sign2;


					edgeNeighborCellSize[96 + 3*(i*4 + j*2 + k) + 0] = 3;
				}
				else if (j == 1 && k == 0)
				{
					edgeCells[96 + 3*(i*4 + j*2 + k) + 0][0][0] = 1;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 0][0][1] = 0;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 0][0][2] = 0;

					edgeCells[96 + 3*(i*4 + j*2 + k) + 0][1][0] = 1;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 0][1][1] = sign1;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 0][1][2] = 0;

					edgeCells[96 + 3*(i*4 + j*2 + k) + 0][2][0] = 0;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 0][2][1] = sign1;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 0][2][2] = 0;

					edgeNeighborCellSize[96 + 3*(i*4 + j*2 + k) + 0] = 3;
				}
				else if (j == 1 && k == 1)
				{
					edgeCells[96 + 3*(i*4 + j*2 + k) + 0][0][0] = 1;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 0][0][1] = 0;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 0][0][2] = 0;

					edgeNeighborCellSize[96 + 3*(i*4 + j*2 + k) + 0] = 1;
				}

				edges[96 + 3*(i*4 + j*2 + k) + 1][0][0] = center[0] + halfLength - innerLength;
				edges[96 + 3*(i*4 + j*2 + k) + 1][0][1] = center[1] + sign1*halfLength + -1*sign1*k*innerLength;
				edges[96 + 3*(i*4 + j*2 + k) + 1][0][2] = center[2] + sign2*halfLength + -1*sign2*j*innerLength;
				edges[96 + 3*(i*4 + j*2 + k) + 1][1][0] = center[0] - halfLength + innerLength;
				edges[96 + 3*(i*4 + j*2 + k) + 1][1][1] = center[1] + sign1*halfLength + -1*sign1*k*innerLength;
				edges[96 + 3*(i*4 + j*2 + k) + 1][1][2] = center[2] + sign2*halfLength + -1*sign2*j*innerLength;

				if (j == 0 && k == 0)
				{

					edgeCells[96 + 3*(i*4 + j*2 + k) + 1][0][0] = 0;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 1][0][1] = sign1;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 1][0][2] = 0;


					edgeCells[96 + 3*(i*4 + j*2 + k) + 1][1][0] = 0;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 1][1][1] = 0;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 1][1][2] = sign2;


					edgeCells[96 + 3*(i*4 + j*2 + k) + 1][2][0] = 0;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 1][2][1] = sign1;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 1][2][2] = sign2;


					edgeNeighborCellSize[96 + 3*(i*4 + j*2 + k) + 1] = 3;
				}
				else if (j == 0 && k == 1)
				{

					edgeCells[96 + 3*(i*4 + j*2 + k) + 1][0][0] = 0;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 1][0][1] = 0;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 1][0][2] = sign2;

					edgeNeighborCellSize[96 + 3*(i*4 + j*2 + k) + 1] = 1;
				}
				else if (j == 1 && k == 0)
				{
					edgeCells[96 + 3*(i*4 + j*2 + k) + 1][0][0] = 0;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 1][0][1] = sign1;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 1][0][2] = 0;

					edgeNeighborCellSize[96 + 3*(i*4 + j*2 + k) + 1] = 1;
				}
				else if (j == 1 && k == 1)
				{
					edgeNeighborCellSize[96 + 3*(i*4 + j*2 + k) + 1] = 0;
				}

				edges[96 + 3*(i*4 + j*2 + k) + 2][0][0] = center[0] - halfLength + innerLength;
				edges[96 + 3*(i*4 + j*2 + k) + 2][0][1] = center[1] + sign1*halfLength + -1*sign1*k*innerLength;
				edges[96 + 3*(i*4 + j*2 + k) + 2][0][2] = center[2] + sign2*halfLength + -1*sign2*j*innerLength;
				edges[96 + 3*(i*4 + j*2 + k) + 2][1][0] = center[0] - halfLength;
				edges[96 + 3*(i*4 + j*2 + k) + 2][1][1] = center[1] + sign1*halfLength + -1*sign1*k*innerLength;
				edges[96 + 3*(i*4 + j*2 + k) + 2][1][2] = center[2] + sign2*halfLength + -1*sign2*j*innerLength;
				if (j == 0 && k == 0)
				{
					edgeCells[96 + 3*(i*4 + j*2 + k) + 2][0][0] = -1;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 2][0][1] = 0;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 2][0][2] = 0;

					edgeCells[96 + 3*(i*4 + j*2 + k) + 2][1][0] = -1;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 2][1][1] = sign1;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 2][1][2] = 0;

					edgeCells[96 + 3*(i*4 + j*2 + k) + 2][2][0] = -1;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 2][2][1] = 0;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 2][2][2] = sign2;

					edgeCells[96 + 3*(i*4 + j*2 + k) + 2][3][0] = -1;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 2][3][1] = sign1;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 2][3][2] = sign2;

					edgeCells[96 + 3*(i*4 + j*2 + k) + 2][4][0] = 0;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 2][4][1] = sign1;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 2][4][2] = 0;

					edgeCells[96 + 3*(i*4 + j*2 + k) + 2][5][0] = 0;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 2][5][1] = 0;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 2][5][2] = sign2;

					edgeCells[96 + 3*(i*4 + j*2 + k) + 2][6][0] = 0;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 2][6][1] = sign1;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 2][6][2] = sign2;

					edgeNeighborCellSize[96 + 3*(i*4 + j*2 + k) + 2] = 7;
				}
				else if (j == 0 && k == 1)
				{
					edgeCells[96 + 3*(i*4 + j*2 + k) + 2][0][0] = -1;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 2][0][1] = 0;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 2][0][2] = 0;

					edgeCells[96 + 3*(i*4 + j*2 + k) + 2][1][0] = -1;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 2][1][1] = 0;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 2][1][2] = sign2;

					edgeCells[96 + 3*(i*4 + j*2 + k) + 2][2][0] = 0;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 2][2][1] = 0;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 2][2][2] = sign2;

					edgeNeighborCellSize[96 + 3*(i*4 + j*2 + k) + 2] = 3;
				}
				else if (j == 1 && k == 0)
				{
					edgeCells[96 + 3*(i*4 + j*2 + k) + 2][0][0] = -1;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 2][0][1] = 0;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 2][0][2] = 0;

					edgeCells[96 + 3*(i*4 + j*2 + k) + 2][1][0] = -1;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 2][1][1] = sign1;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 2][1][2] = 0;

					edgeCells[96 + 3*(i*4 + j*2 + k) + 2][2][0] = 0;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 2][2][1] = sign1;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 2][2][2] = 0;

					edgeNeighborCellSize[96 + 3*(i*4 + j*2 + k) + 2] = 3;
				}
				else if (j == 1 && k == 1)
				{
					edgeCells[96 + 3*(i*4 + j*2 + k) + 2][0][0] = -1;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 2][0][1] = 0;
					edgeCells[96 + 3*(i*4 + j*2 + k) + 2][0][2] = 0;

					edgeNeighborCellSize[96 + 3*(i*4 + j*2 + k) + 2] = 1;
				}
			}
		}
	}

	double fv1 = 0;
	double fv2 = 0;
	int x = 0;
	int y = 0;
	int z = 0;

	for (int i = 0; i < 144; i++)
	{
		fv1 = edges[i][0][0]*plane[0] + edges[i][0][1]*plane[1] + edges[i][0][2]*plane[2] + b;
		fv2 = edges[i][1][0]*plane[0] + edges[i][1][1]*plane[1] + edges[i][1][2]*plane[2] + b;
		if (fv1*fv2 < 0 || fabs(fv1) < 0.00000001 || fabs(fv2) < 0.00000001)
		{
			bool allEmpty = true;
			for (int nc = 0; nc < edgeNeighborCellSize[i]; nc++)
			{
				x = cx + edgeCells[i][nc][0];
				y = cy + edgeCells[i][nc][1];
				z = cz + edgeCells[i][nc][2];
				if (x > 2 || y > 2 || z > 2 || x < 0 || y < 0 || z < 0 || gc[x][y][z].getNumOfVertex() > 0)
				{
					allEmpty = false;
					break;
				}
				else
				{
					const double* neighborCenter = gc[x][y][z].getCenter();
					if (!allEmpty)
						break;
				}
			}
			if (allEmpty)
				return true;
		}
	}
	return false;
}

const vector<int> & gridCell::getInCellVertex() const
{
	if (workingNode == NULL)
		return *(new vector<int>);
	return workingNode->getVertices();
}

int gridCell::getNumOfVertex() const
{
	if (workingNode == NULL)
		return 0;
	return workingNode->getNumOfVertex();
}

void gridCell::setWorkingNode(octreeNode* n) {
	workingNode = n;
}

octreeNode* gridCell::getWorkingNode() {
	return workingNode;
}

const double* gridCell::getAvgPoint() const {
	return workingNode->getAvgPoint();
}
