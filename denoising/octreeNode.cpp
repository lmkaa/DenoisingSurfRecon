#include "octreeNode.hpp"
#include <sys/time.h>
#include <fstream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <cmath>

void iconvertToOrder2Poly(const double* orgP, double* newP) {
	newP[0] = orgP[0]*orgP[0] + orgP[1]*orgP[1] + orgP[2]*orgP[2];
	newP[1] = orgP[0];
	newP[2] = orgP[1];
	newP[3] = orgP[2];
	newP[4] = 1;
}

void igetJacM(const double* orgP, double jacM[][3]) {
	jacM[0][0] = 2*orgP[1];
	jacM[0][1] = 2*orgP[2];
	jacM[0][2] = 2*orgP[3];
	jacM[1][0] = 1;
	jacM[2][1] = 1;
	jacM[3][2] = 1;
}

bool isLocallyFlat(const vector<const double*>& points, int pointSize, const Vector3d& cp, const Vector3d& normal) {
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

void solvePCA(vector<const double*>& tmpPoint, int nvsize, vector<double>& weight, double totalWeight, gsl_vector* eval, gsl_matrix* evec, int dim) {
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

void solveAD(vector<const double*>& tmpPoint, int nvsize, vector<double>& weight, double totalWeight, gsl_vector* eval, gsl_matrix* evec, int dim) {
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
		igetJacM(tmpPoint[i],cJac);
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

int getSameLevelNeighborIndex(const double* c1, const double* c2, double length)
{
	if (fabs(c1[0] - c2[0]) < length/1000 && fabs(c1[1] - c2[1]) < length/1000 && (c1[2] - c2[2]) > (length - length/1000))
		return X0_Y0_ZP;
	else if (fabs(c1[0] - c2[0]) < length/1000 && fabs(c1[1] - c2[1]) < length/1000 && (c2[2] - c1[2]) > (length - length/1000))
		return X0_Y0_ZN;
	else if (fabs(c1[0] - c2[0]) < length/1000 && (c1[1] - c2[1]) > (length - length/1000) && fabs(c1[2] - c2[2]) < length/1000)
		return X0_YP_Z0;
	else if (fabs(c1[0] - c2[0]) < length/1000 && (c2[1] - c1[1]) > (length - length/1000) && fabs(c1[2] - c2[2]) < length/1000)
		return X0_YN_Z0;
	else if ((c1[0] - c2[0]) > (length - length/1000) && fabs(c1[1] - c2[1]) < length/1000 && fabs(c1[2] - c2[2]) < length/1000)
		return XP_Y0_Z0;
	else if ((c2[0] - c1[0]) > (length - length/1000) && fabs(c1[1] - c2[1]) < length/1000 && fabs(c1[2] - c2[2]) < length/1000)
		return XN_Y0_Z0;
	else if ((c1[0] - c2[0]) > (length - length/1000) && (c1[1] - c2[1]) > (length - length/1000) && fabs(c1[2] - c2[2]) < length/1000)
		return XP_YP_Z0;
	else if ((c1[0] - c2[0]) > (length - length/1000) && (c2[1] - c1[1]) > (length - length/1000) && fabs(c1[2] - c2[2]) < length/1000)
		return XP_YN_Z0;
	else if ((c2[0] - c1[0]) > (length - length/1000) && (c1[1] - c2[1]) > (length - length/1000) && fabs(c1[2] - c2[2]) < length/1000)
		return XN_YP_Z0;
	else if ((c2[0] - c1[0]) > (length - length/1000) && (c2[1] - c1[1]) > (length - length/1000) && fabs(c1[2] - c2[2]) < length/1000)
		return XN_YN_Z0;
	else if ((c1[0] - c2[0]) > (length - length/1000) && fabs(c1[1] - c2[1]) < length/1000 && (c1[2] - c2[2]) > (length - length/1000))
		return XP_Y0_ZP;
	else if ((c1[0] - c2[0]) > (length - length/1000) && fabs(c1[1] - c2[1]) < length/1000 && (c2[2] - c1[2]) > (length - length/1000))
		return XP_Y0_ZN;
	else if ((c2[0] - c1[0]) > (length - length/1000) && fabs(c1[1] - c2[1]) < length/1000 && (c1[2] - c2[2]) > (length - length/1000))
		return XN_Y0_ZP;
	else if ((c2[0] - c1[0]) > (length - length/1000) && fabs(c1[1] - c2[1]) < length/1000 && (c2[2] - c1[2]) > (length - length/1000))
		return XN_Y0_ZN;
	else if (fabs(c1[0] - c2[0]) < length/1000 && (c1[1] - c2[1]) > (length - length/1000) && (c1[2] - c2[2]) > (length - length/1000))
		return X0_YP_ZP;
	else if (fabs(c1[0] - c2[0]) < length/1000 && (c1[1] - c2[1]) > (length - length/1000) && (c2[2] - c1[2]) > (length - length/1000))
		return X0_YP_ZN;
	else if (fabs(c1[0] - c2[0]) < length/1000 && (c2[1] - c1[1]) > (length - length/1000) && (c1[2] - c2[2]) > (length - length/1000))
		return X0_YN_ZP;
	else if (fabs(c1[0] - c2[0]) < length/1000 && (c2[1] - c1[1]) > (length - length/1000) && (c2[2] - c1[2]) > (length - length/1000))
		return X0_YN_ZN;
	else if ((c1[0] - c2[0]) > (length - length/1000) && (c1[1] - c2[1]) > (length - length/1000) && (c1[2] - c2[2]) > (length - length/1000))
		return XP_YP_ZP;
	else if ((c1[0] - c2[0]) > (length - length/1000) && (c1[1] - c2[1]) > (length - length/1000) && (c2[2] - c1[2]) > (length - length/1000))
		return XP_YP_ZN;
	else if ((c1[0] - c2[0]) > (length - length/1000) && (c2[1] - c1[1]) > (length - length/1000) && (c1[2] - c2[2]) > (length - length/1000))
		return XP_YN_ZP;
	else if ((c1[0] - c2[0]) > (length - length/1000) && (c2[1] - c1[1]) > (length - length/1000) && (c2[2] - c1[2]) > (length - length/1000))
		return XP_YN_ZN;
	else if ((c2[0] - c1[0]) > (length - length/1000) && (c1[1] - c2[1]) > (length - length/1000) && (c1[2] - c2[2]) > (length - length/1000))
		return XN_YP_ZP;
	else if ((c2[0] - c1[0]) > (length - length/1000) && (c1[1] - c2[1]) > (length - length/1000) && (c2[2] - c1[2]) > (length - length/1000))
		return XN_YP_ZN;
	else if ((c2[0] - c1[0]) > (length - length/1000) && (c2[1] - c1[1]) > (length - length/1000) && (c1[2] - c2[2]) > (length - length/1000))
		return XN_YN_ZP;
	else if ((c2[0] - c1[0]) > (length - length/1000) && (c2[1] - c1[1]) > (length - length/1000) && (c2[2] - c1[2]) > (length - length/1000))
		return XN_YN_ZN;
}

void getOneLevelChildPoints(octreeNode* currNode, double* refPoint, double sideLength, vector<const double*>& tmpPoint, vector<double>& weight, double& totalWeight) {
	octreeNode* const * currChild = currNode->getChildNodes();
	bool noChild = true;
	for (int ci = 0; ci < 8; ci++) {
		if (currChild[ci] != NULL && !(currChild[ci]->isDeleted()) && currChild[ci]->getNumOfVertex() > 0) {
			double cd = 0;
			const double* cap = currChild[ci]->getAvgPoint();
			for (int d = 0; d < 3; d++)
				cd += ((cap[d] - refPoint[d])*(cap[d] - refPoint[d]));
			if (cd < 3*sideLength) {
				noChild = false;
				double cw = currChild[ci]->getNumOfVertex()*exp(-1*(cd*cd)/(2*sideLength*sideLength));
				weight.push_back(cw);
				totalWeight += cw;
				tmpPoint.push_back(cap);
			}
		}
	}
	if (noChild) {
		double cd = 0;
		const double* cap = currNode->getAvgPoint();
		for (int d = 0; d < 3; d++)
			cd += ((cap[d] - refPoint[d])*(cap[d] - refPoint[d]));
		if (cd < 3*sideLength) {
			noChild = false;
			double cw = currNode->getNumOfVertex()*exp(-1*(cd*cd)/(2*sideLength*sideLength));
			weight.push_back(cw);
			totalWeight += cw;
			tmpPoint.push_back(cap);
		}
	}
}

void getTwoLevelChildPoints(octreeNode* currNode, double* refPoint, double sideLength, vector<const double*>& tmpPoint, vector<double>& weight, double& totalWeight, int dim) {
	octreeNode* const * currChild = currNode->getChildNodes();
	bool noChild = true;
	for (int ci = 0; ci < 8; ci++) {
		if (currChild[ci] != NULL && !(currChild[ci]->isDeleted()) && currChild[ci]->getNumOfVertex() > 0) {
			octreeNode* const * currChildChild = currChild[ci]->getChildNodes();
			bool noChildChild = true;
			for (int cci = 0; cci < 8; cci++) {
				if (currChildChild[cci] != NULL && !(currChildChild[cci]->isDeleted()) && currChildChild[cci]->getNumOfVertex() > 0) {
					double cd = 0;
					const double* cap = currChildChild[cci]->getAvgPoint();
					for (int d = 0; d < 3; d++)
						cd += ((cap[d] - refPoint[d])*(cap[d] - refPoint[d]));
					cd = sqrt(cd);
					if (cd < 4*sideLength) {
						noChildChild = false;
						double cw = currChildChild[cci]->getNumOfVertex()*exp(-1*(cd*cd)/(2*sideLength*sideLength));
						weight.push_back(cw);
						totalWeight += cw;
						double* polyCap = new double[dim];
						iconvertToOrder2Poly(cap, polyCap);
						tmpPoint.push_back(polyCap);
					}
				}
			}
			if (!noChildChild)
				noChild = false;
			else {
				double cd = 0;
				const double* cap = currChild[ci]->getAvgPoint();
				for (int d = 0; d < 3; d++)
					cd += ((cap[d] - refPoint[d])*(cap[d] - refPoint[d]));
				cd = sqrt(cd);
				if (cd < 4*sideLength) {
					noChild = false;
					double cw = currChild[ci]->getNumOfVertex()*exp(-1*(cd*cd)/(2*sideLength*sideLength));
					weight.push_back(cw);
					totalWeight += cw;
					double* polyCap = new double[dim];
					iconvertToOrder2Poly(cap, polyCap);
					tmpPoint.push_back(polyCap);
				}
			}
		}
	}
	if (noChild) {
		double cd = 0;
		const double* cap = currNode->getAvgPoint();
		for (int d = 0; d < 3; d++)
			cd += ((cap[d] - refPoint[d])*(cap[d] - refPoint[d]));
		cd = sqrt(cd);
		if (cd < 4*sideLength) {
			noChild = false;
			double cw = currNode->getNumOfVertex()*exp(-1*(cd*cd)/(1*sideLength*sideLength));
			weight.push_back(cw);
			totalWeight += cw;
			double* polyCap = new double[dim];
			iconvertToOrder2Poly(cap, polyCap);
			tmpPoint.push_back(polyCap);
		}
	}
}

octreeNode::~octreeNode() {
	inNodeVertices.clear();
	vertexPick.clear();
	for (int k = 0; k < 8; k++) {
		if (childNodes[k] != NULL) {
			delete childNodes[k];
		}
	}
}

octreeNode::octreeNode() {
	for (int i = 0; i < 3; i++)
		center[i] = 0;
	sideLength = 0;
	parentNode = NULL;
	nodeIndex = -1;
	for (int i = 0; i < 8; i++)
		childNodes[i] = NULL;
	for (int i = 0; i < 26; i++)
		sameLevelNeighbors[i] = NULL;
	isLeafNode = false;
	isDeletedNode = false;
	isProcessedNode = false;
	numOfVertex = 0;
	numOfKeepPoint = 0;
	for (int d = 0; d < 3; d++) {
		maxCoord[d] = 0;
		minCoord[d] = 0;
	}
}

int octreeNode::getNodeIndex() {
	return nodeIndex;
}

void octreeNode::setNodeIndex(int n) {
	nodeIndex = n;
}

const double* octreeNode::getCenter() const {
	return center;
}

double octreeNode::getSideLength() const {
	return sideLength;
}

octreeNode* const octreeNode::getParentNode() const {
	return parentNode;
}

octreeNode* const * octreeNode::getChildNodes() const {
	return childNodes;
}

octreeNode* const *  octreeNode::getSameLevelNeighbors() const {
	return sameLevelNeighbors;
}

bool octreeNode::isLeaf() const {
	return isLeafNode;
}

void octreeNode::setLeaf(bool _isLeafNode) {
	isLeafNode = _isLeafNode;
}
void octreeNode::setCenter(double _center[3]) {
	for (int i = 0; i < 3; i++)
		center[i] = _center[i];
}
void octreeNode::setSideLength(double _sideLength) {
	sideLength = _sideLength;
}
void octreeNode::setParentNode(octreeNode* _parentNode) {
	parentNode = _parentNode;
}
void octreeNode::setChildNode(int i, octreeNode* _childNode) {
	childNodes[i] = _childNode;
}
void octreeNode::setSameLevelNeighbor(int i, octreeNode* _neighborCell) {
	sameLevelNeighbors[i] = _neighborCell;
}

const vector<int> & octreeNode::getVertices() const {
	return inNodeVertices;
}

void octreeNode::insertVertex(int i, const vector<double*> & v) {
	inNodeVertices.push_back(i);
	vertexPick.push_back(false);
	numOfVertex++;
	numOfKeepPoint++;
	for (int d = 0; d < 3; d++)
	{
		if (numOfVertex == 1 || v[i][d] > maxCoord[d])
			maxCoord[d] = v[i][d];
		if (numOfVertex == 1 || v[i][d] < minCoord[d])
			minCoord[d] = v[i][d];
	}
}

void octreeNode::insertVertexForMesh(int i) {
	inNodeVertices.push_back(i);
	vertexPick.push_back(false);
	numOfVertex++;
	numOfKeepPoint++;
}

int octreeNode::getLevel() const {
	return level;
}
void octreeNode::setLevel(int _level) {
	level = _level;
}

const double* octreeNode::getMaxCoord() const {
	return maxCoord;
}
const double* octreeNode::getMinCoord() const {
	return minCoord;
}

void octreeNode::setRootNode() {
	double maxDist = -1;
	for (int d = 0; d < 3; d++) {
		center[d] = (maxCoord[d] + minCoord[d])/2;
		if (d == 0 || maxDist < fabs(maxCoord[d] - minCoord[d]))
			maxDist = fabs(maxCoord[d] - minCoord[d]);
	}
	sideLength = maxDist;
	isLeafNode = true;
	parentNode = NULL;
	level = 1;
}

void octreeNode::setInternalNode(octreeNode* _parentNode) {
	sideLength = _parentNode->getSideLength()/2;
	isLeafNode = true;
	parentNode = _parentNode;
	level = _parentNode->getLevel() + 1;
}

void octreeNode::setSameLevelNeighbors(int i, octreeNode* _childNodes[]) {
	// Set same parent neighbors
	if (i == 0)
	{
		_childNodes[i]->setSameLevelNeighbor(X0_Y0_ZN, _childNodes[1]);
		_childNodes[i]->setSameLevelNeighbor(X0_YN_Z0, _childNodes[2]);
		_childNodes[i]->setSameLevelNeighbor(X0_YN_ZN, _childNodes[3]);
		_childNodes[i]->setSameLevelNeighbor(XN_Y0_Z0, _childNodes[4]);
		_childNodes[i]->setSameLevelNeighbor(XN_Y0_ZN, _childNodes[5]);
		_childNodes[i]->setSameLevelNeighbor(XN_YN_Z0, _childNodes[6]);
		_childNodes[i]->setSameLevelNeighbor(XN_YN_ZN, _childNodes[7]);
	}
	else if (i == 1)
	{
		_childNodes[i]->setSameLevelNeighbor(X0_Y0_ZP, _childNodes[0]);
		_childNodes[i]->setSameLevelNeighbor(X0_YN_ZP, _childNodes[2]);
		_childNodes[i]->setSameLevelNeighbor(X0_YN_Z0, _childNodes[3]);
		_childNodes[i]->setSameLevelNeighbor(XN_Y0_ZP, _childNodes[4]);
		_childNodes[i]->setSameLevelNeighbor(XN_Y0_Z0, _childNodes[5]);
		_childNodes[i]->setSameLevelNeighbor(XN_YN_ZP, _childNodes[6]);
		_childNodes[i]->setSameLevelNeighbor(XN_YN_Z0, _childNodes[7]);
	}
	else if (i == 2)
	{
		_childNodes[i]->setSameLevelNeighbor(X0_YP_Z0, _childNodes[0]);
		_childNodes[i]->setSameLevelNeighbor(X0_YP_ZN, _childNodes[1]);
		_childNodes[i]->setSameLevelNeighbor(X0_Y0_ZN, _childNodes[3]);
		_childNodes[i]->setSameLevelNeighbor(XN_YP_Z0, _childNodes[4]);
		_childNodes[i]->setSameLevelNeighbor(XN_YP_ZN, _childNodes[5]);
		_childNodes[i]->setSameLevelNeighbor(XN_Y0_Z0, _childNodes[6]);
		_childNodes[i]->setSameLevelNeighbor(XN_Y0_ZN, _childNodes[7]);
	}
	else if (i == 3)
	{
		_childNodes[i]->setSameLevelNeighbor(X0_YP_ZP, _childNodes[0]);
		_childNodes[i]->setSameLevelNeighbor(X0_YP_Z0, _childNodes[1]);
		_childNodes[i]->setSameLevelNeighbor(X0_Y0_ZP, _childNodes[2]);
		_childNodes[i]->setSameLevelNeighbor(XN_YP_ZP, _childNodes[4]);
		_childNodes[i]->setSameLevelNeighbor(XN_YP_Z0, _childNodes[5]);
		_childNodes[i]->setSameLevelNeighbor(XN_Y0_ZP, _childNodes[6]);
		_childNodes[i]->setSameLevelNeighbor(XN_Y0_Z0, _childNodes[7]);
	}
	else if (i == 4)
	{
		_childNodes[i]->setSameLevelNeighbor(XP_Y0_Z0, _childNodes[0]);
		_childNodes[i]->setSameLevelNeighbor(XP_Y0_ZN, _childNodes[1]);
		_childNodes[i]->setSameLevelNeighbor(XP_YN_Z0, _childNodes[2]);
		_childNodes[i]->setSameLevelNeighbor(XP_YN_ZN, _childNodes[3]);
		_childNodes[i]->setSameLevelNeighbor(X0_Y0_ZN, _childNodes[5]);
		_childNodes[i]->setSameLevelNeighbor(X0_YN_Z0, _childNodes[6]);
		_childNodes[i]->setSameLevelNeighbor(X0_YN_ZN, _childNodes[7]);
	}
	else if (i == 5)
	{
		_childNodes[i]->setSameLevelNeighbor(XP_Y0_ZP, _childNodes[0]);
		_childNodes[i]->setSameLevelNeighbor(XP_Y0_Z0, _childNodes[1]);
		_childNodes[i]->setSameLevelNeighbor(XP_YN_ZP, _childNodes[2]);
		_childNodes[i]->setSameLevelNeighbor(XP_YN_Z0, _childNodes[3]);
		_childNodes[i]->setSameLevelNeighbor(X0_Y0_ZP, _childNodes[4]);
		_childNodes[i]->setSameLevelNeighbor(X0_YN_ZP, _childNodes[6]);
		_childNodes[i]->setSameLevelNeighbor(X0_YN_Z0, _childNodes[7]);
	}
	else if (i == 6)
	{
		_childNodes[i]->setSameLevelNeighbor(XP_YP_Z0, _childNodes[0]);
		_childNodes[i]->setSameLevelNeighbor(XP_YP_ZN, _childNodes[1]);
		_childNodes[i]->setSameLevelNeighbor(XP_Y0_Z0, _childNodes[2]);
		_childNodes[i]->setSameLevelNeighbor(XP_Y0_ZN, _childNodes[3]);
		_childNodes[i]->setSameLevelNeighbor(X0_YP_Z0, _childNodes[4]);
		_childNodes[i]->setSameLevelNeighbor(X0_YP_ZN, _childNodes[5]);
		_childNodes[i]->setSameLevelNeighbor(X0_Y0_ZN, _childNodes[7]);
	}
	else if (i == 7)
	{
		_childNodes[i]->setSameLevelNeighbor(XP_YP_ZP, _childNodes[0]);
		_childNodes[i]->setSameLevelNeighbor(XP_YP_Z0, _childNodes[1]);
		_childNodes[i]->setSameLevelNeighbor(XP_Y0_ZP, _childNodes[2]);
		_childNodes[i]->setSameLevelNeighbor(XP_Y0_Z0, _childNodes[3]);
		_childNodes[i]->setSameLevelNeighbor(X0_YP_ZP, _childNodes[4]);
		_childNodes[i]->setSameLevelNeighbor(X0_YP_Z0, _childNodes[5]);
		_childNodes[i]->setSameLevelNeighbor(X0_Y0_ZP, _childNodes[6]);
	}
	// set the rest neighbors
	const double* childCenter = _childNodes[i]->getCenter();
	octreeNode* const * pSameLevelNeighbors = parentNode->getSameLevelNeighbors();
	for (int j = 0; j < 26; j++)
	{
		if (pSameLevelNeighbors[j] != NULL && !(pSameLevelNeighbors[j]->isLeaf()))
		{
			octreeNode* const * neighborChild = pSameLevelNeighbors[j]->getChildNodes();
			for (int k = 0; k < 8; k++)
			{
				if (neighborChild[k] != NULL) {
					const double* ncCenter = neighborChild[k]->getCenter();
					if (fabs(ncCenter[0] - childCenter[0]) < (sideLength + sideLength/2000) && fabs(ncCenter[1] - childCenter[1]) < (sideLength + sideLength/2000) && fabs(ncCenter[2] - childCenter[2]) < (sideLength + sideLength/2000))
					{
						int childToNeighborIndex = getSameLevelNeighborIndex(ncCenter, childCenter, sideLength/2);
						int neighborToChildIndex = getSameLevelNeighborIndex(childCenter, ncCenter, sideLength/2);
						_childNodes[i]->setSameLevelNeighbor(childToNeighborIndex, neighborChild[k]);
						neighborChild[k]->setSameLevelNeighbor(neighborToChildIndex, _childNodes[i]);
					}
				}
			}
		}
	}
}

int octreeNode::getNumOfVertex() const {
	return numOfVertex;
}

void octreeNode::balanceNode(const vector<double*> & v, bool* keepPoint, queue<octreeNode*> & myqueue) {
	if (level > 2) {
		octreeNode* const * pNeighbor = parentNode->getSameLevelNeighbors();
		for (int i = 0; i < 26; i++) {
			if (pNeighbor[i] != NULL && pNeighbor[i]->isLeaf() && pNeighbor[i]->getNumOfVertex() > 0) {
				const double* pnCenter = pNeighbor[i]->getCenter();
				if (fabs(pnCenter[0] - center[0]) < ((3*sideLength)/2 + sideLength/2000) && fabs(pnCenter[1] - center[1]) < ((3*sideLength)/2 + sideLength/2000) && fabs(pnCenter[2] - center[2]) < ((3*sideLength)/2 + sideLength/2000)) {
					pNeighbor[i]->splitNode(v, keepPoint, myqueue);
				}
			}
		}
	}
}

bool octreeNode::checkSplitNode() {
	if (inNodeVertices.size() > 1)
	{
		int minGridIndex[3];
		int maxGridIndex[3];
		bool allEqual = true;
		for (int d = 0; d < 3; d++)
		{
			minGridIndex[d] = ((int) floor(8*fabs(minCoord[d] - (center[d] - sideLength/2))/sideLength));
			maxGridIndex[d] = ((int) floor(8*fabs(maxCoord[d] - (center[d] - sideLength/2))/sideLength));
			if (minGridIndex[d] > 7)
			{
				minGridIndex[d] = 7;
			}
			else if (minGridIndex[d] < 0)
			{
				minGridIndex[d] = 0;
			}
			if (maxGridIndex[d] > 7)
			{
				maxGridIndex[d] = 7;
			}
			else if (maxGridIndex[d] < 0)
			{
				maxGridIndex[d] = 0;
			}
			if (minGridIndex[d] != maxGridIndex[d])
			{
				allEqual = false;
			}
		}
		if (allEqual)
			return false;
		else
			return true;
	}
	else
		return false;

}

void octreeNode::createSingleNode(int ci, const vector<double*>& v, int vi) {
	if (childNodes[ci] == NULL) {
		int i = 0, j = 0, k = 0;
		if (ci == 0) {
			i = 0; j = 0; k = 0;
		}
		else if (ci == 1) {
			i = 0; j = 0; k = 1;
		}
		else if (ci == 2) {
			i = 0; j = 1; k = 0;
		}
		else if (ci == 3) {
			i = 0; j = 1; k = 1;
		}
		else if (ci == 4) {
			i = 1; j = 0; k = 0;
		}
		else if (ci == 5) {
			i = 1; j = 0; k = 1;
		}
		else if (ci == 6) {
			i = 1; j = 1; k = 0;
		}
		else if (ci == 7) {
			i = 1; j = 1; k = 1;
		}
		childNodes[ci] = new octreeNode();
		double tmpChildCenter[3];
		tmpChildCenter[0] = center[0] + (1 - 2*i)*sideLength/4;
		tmpChildCenter[1] = center[1] + (1 - 2*j)*sideLength/4;
		tmpChildCenter[2] = center[2] + (1 - 2*k)*sideLength/4;
		childNodes[ci]->setCenter(tmpChildCenter);
		childNodes[ci]->setInternalNode(this);
		childNodes[ci]->insertVertex(vi, v);
		childNodes[ci]->setSameLevelNeighbors(ci, childNodes);
		childNodes[ci]->computeAvgPoint(v);
	}
	else {
		cout << "You create a node that is already exists" << endl;
		cin.get();
	}
}

void octreeNode::splitNode(const vector<double*> & v, bool* keepPoint, queue<octreeNode*> & myqueue) {
		if (isLeafNode) {
			isLeafNode = false;
			double tmpChildCenter[3];
			for (int i = 0; i < 2 ; i++)
			{
				for (int j = 0; j < 2; j++)
				{
					for (int k = 0; k < 2; k++)
					{
						int ci = 2*2*i + 2*j + k;
						childNodes[ci] = new octreeNode();
						tmpChildCenter[0] = center[0] + (1 - 2*i)*sideLength/4;
						tmpChildCenter[1] = center[1] + (1 - 2*j)*sideLength/4;
						tmpChildCenter[2] = center[2] + (1 - 2*k)*sideLength/4;
						childNodes[ci]->setCenter(tmpChildCenter);
						childNodes[ci]->setInternalNode(this);
					}
				}
			}
			int vsize = inNodeVertices.size();
			for (int i = 0; i < vsize; i++)
			{
				if (keepPoint[inNodeVertices[i]]) {
					int vIndex = inNodeVertices[i];
					if (v[vIndex][0] >= center[0] && v[vIndex][1] >= center[1] && v[vIndex][2] >= center[2])
						childNodes[0]->insertVertex(vIndex, v);
					else if (v[vIndex][0] >= center[0] && v[vIndex][1] >= center[1] && v[vIndex][2] < center[2])
						childNodes[1]->insertVertex(vIndex, v);
					else if (v[vIndex][0] >= center[0] && v[vIndex][1] < center[1] && v[vIndex][2] >= center[2])
						childNodes[2]->insertVertex(vIndex, v);
					else if (v[vIndex][0] >= center[0] && v[vIndex][1] < center[1] && v[vIndex][2] < center[2])
						childNodes[3]->insertVertex(vIndex, v);
					else if (v[vIndex][0] < center[0] && v[vIndex][1] >= center[1] && v[vIndex][2] >= center[2])
						childNodes[4]->insertVertex(vIndex, v);
					else if (v[vIndex][0] < center[0] && v[vIndex][1] >= center[1] && v[vIndex][2] < center[2])
						childNodes[5]->insertVertex(vIndex, v);
					else if (v[vIndex][0] < center[0] && v[vIndex][1] < center[1] && v[vIndex][2] >= center[2])
						childNodes[6]->insertVertex(vIndex, v);
					else if (v[vIndex][0] < center[0] && v[vIndex][1] < center[1] && v[vIndex][2] < center[2])
						childNodes[7]->insertVertex(vIndex, v);
					else {
						cout << "a node does not in any child" << endl;
						cin.get();
					}
				}
			}
			for (int i = 0; i < 8; i++) {
				if (childNodes[i]->getNumOfVertex() == 0) {
					delete childNodes[i];
					childNodes[i] = NULL;
				}
			}
			for (int i = 0; i < 8; i++) {
				if (childNodes[i] != NULL && childNodes[i]->getNumOfVertex() > 0) {
					childNodes[i]->setSameLevelNeighbors(i, childNodes);
					childNodes[i]->computeAvgPoint(v);
					myqueue.push(childNodes[i]);
				}
			}
			balanceNode(v, keepPoint, myqueue);
		}
}

void octreeNode::splitNodeToSameLeafLevel(const vector<double*> & v, queue<octreeNode*> & myqueue, bool& keepSplit) {
	isLeafNode = false;
	double tmpChildCenter[3];
	for (int i = 0; i < 2 ; i++) {
		for (int j = 0; j < 2; j++) {
			for (int k = 0; k < 2; k++) {
				int ci = 2*2*i + 2*j + k;
				childNodes[ci] = new octreeNode();
				tmpChildCenter[0] = center[0] + (1 - 2*i)*sideLength/4;
				tmpChildCenter[1] = center[1] + (1 - 2*j)*sideLength/4;
				tmpChildCenter[2] = center[2] + (1 - 2*k)*sideLength/4;
				childNodes[ci]->setCenter(tmpChildCenter);
				childNodes[ci]->setInternalNode(this);
			}
		}
	}
	int vsize = inNodeVertices.size();
	for (int i = 0; i < vsize; i++) {
		int vIndex = inNodeVertices[i];
		if (v[vIndex][0] >= center[0] && v[vIndex][1] >= center[1] && v[vIndex][2] >= center[2])
			childNodes[0]->insertVertex(vIndex, v);
		else if (v[vIndex][0] >= center[0] && v[vIndex][1] >= center[1] && v[vIndex][2] < center[2])
			childNodes[1]->insertVertex(vIndex, v);
		else if (v[vIndex][0] >= center[0] && v[vIndex][1] < center[1] && v[vIndex][2] >= center[2])
			childNodes[2]->insertVertex(vIndex, v);
		else if (v[vIndex][0] >= center[0] && v[vIndex][1] < center[1] && v[vIndex][2] < center[2])
			childNodes[3]->insertVertex(vIndex, v);
		else if (v[vIndex][0] < center[0] && v[vIndex][1] >= center[1] && v[vIndex][2] >= center[2])
			childNodes[4]->insertVertex(vIndex, v);
		else if (v[vIndex][0] < center[0] && v[vIndex][1] >= center[1] && v[vIndex][2] < center[2])
			childNodes[5]->insertVertex(vIndex, v);
		else if (v[vIndex][0] < center[0] && v[vIndex][1] < center[1] && v[vIndex][2] >= center[2])
			childNodes[6]->insertVertex(vIndex, v);
		else if (v[vIndex][0] < center[0] && v[vIndex][1] < center[1] && v[vIndex][2] < center[2])
			childNodes[7]->insertVertex(vIndex, v);
	}
	for (int i = 0; i < 8; i++) {
		if (childNodes[i]->getNumOfVertex() == 0) {
			delete childNodes[i];
			childNodes[i] = NULL;
		}
	}
	for (int i = 0; i < 8; i++) {
		if (childNodes[i] != NULL && childNodes[i]->getNumOfVertex() > 0) {
			childNodes[i]->setSameLevelNeighbors(i, childNodes);
			childNodes[i]->computeAvgPoint(v);
			myqueue.push(childNodes[i]);
			if (childNodes[i]->checkSplitNode())
				keepSplit = true;
		}
	}
}

void octreeNode::setDelete(bool _isDeletedNode) {
	isDeletedNode = _isDeletedNode;
}

bool octreeNode::isDeleted() const {
	return isDeletedNode;
}

bool octreeNode::isProcessed() {
	return isProcessedNode;
}

void octreeNode::setProcessed(bool _isProcessedNode) {
	isProcessedNode = _isProcessedNode;
}

bool octreeNode::checkSplitNodeForMesh() {
	if (inNodeVertices.size() > 30)
		return true;
	else
		return false;

}

void octreeNode::setRootNodeForMesh() {
	double maxDist = -1;
	for (int d = 0; d < 3; d++) 	{
		center[d] = (maxCoord[d] + minCoord[d])/2;
		if (d == 0 || maxDist < fabs(maxCoord[d] - minCoord[d]))
			maxDist = fabs(maxCoord[d] - minCoord[d]);
	}
	sideLength = maxDist + 0.1;
	isLeafNode = true;
	parentNode = NULL;
	level = 1;
}


void octreeNode::computeAvgPoint(const vector<double*> &v) {
	for (int d = 0; d < 3; d++)
		avgPoint[d] = 0;
	for (int i = 0; i < numOfVertex; i++) {
		for (int d = 0; d < 3; d++)
			avgPoint[d] = avgPoint[d] + v[inNodeVertices[i]][d];
	}
	for (int d = 0; d < 3; d++)
		avgPoint[d] = avgPoint[d]/numOfVertex;
	double minDist = 3*sideLength;
	for (int i = 0; i < numOfVertex; i++) {
		double currDist = 0;
		for (int d = 0; d < 3; d++)
			currDist += ((v[inNodeVertices[i]][d] - avgPoint[d])*(v[inNodeVertices[i]][d] - avgPoint[d]));
		if (currDist < minDist) {
			minDist = currDist;
			nearestToAvgPoint = inNodeVertices[i];
		}
	}
}

void octreeNode::updateAvgPointWithKeepPoint(const vector<double*> &v, bool* keepPoint) {
	for (int d = 0; d < 3; d++)
		avgPoint[d] = 0;
	int numOfCount = 0;
	for (int i = 0; i < numOfVertex; i++) {
		if (keepPoint[inNodeVertices[i]]) {
			numOfCount++;
			for (int d = 0; d < 3; d++)
				avgPoint[d] = avgPoint[d] + v[inNodeVertices[i]][d];
		}
	}
	for (int d = 0; d < 3; d++)
		avgPoint[d] = avgPoint[d]/numOfCount;
	double minDist = 3*sideLength;
	nearestToAvgPoint = -1;
	for (int i = 0; i < numOfVertex; i++) {
		if (keepPoint[inNodeVertices[i]]) {
			double currDist = 0;
			for (int d = 0; d < 3; d++)
				currDist += ((v[inNodeVertices[i]][d] - avgPoint[d])*(v[inNodeVertices[i]][d] - avgPoint[d]));
			if (currDist < minDist) {
				minDist = currDist;
				nearestToAvgPoint = inNodeVertices[i];
			}
		}
	}
}

const double* octreeNode::getAvgPoint() const {
	return avgPoint;
}

double* octreeNode::getAvgPoint() {
	return avgPoint;
}

void octreeNode::setAvgPoint(const double* p) {
	for (int d = 0; d < 3; d++)
		avgPoint[d] = p[d];
}

int octreeNode::getNearestToAvgPoint() const {
	return nearestToAvgPoint;
}

void octreeNode::setOrder1Func(const double* p) {
	for (int d = 0; d < 3; d++)
		order1Func[d] = p[d];
}
void octreeNode::setOrder1b(double pb) {
	order1b = pb;
}

bool octreeNode::computeOrder1FunctionPCA(const vector<double*>& v, bool* keepPoint, int dim) {
	bool keepGoUp = false;
	octreeNode* cNode = this;
	do {
		keepGoUp = false;
		if (cNode->getParentNode() == NULL || cNode->getParentNode()->getParentNode() == NULL)
			return false;
		octreeNode* currNode = cNode->getParentNode()->getParentNode();
		octreeNode* const * cNodeN = cNode->getParentNode()->getSameLevelNeighbors();
		double refPoint[3];
		double nw = cNode->getNumOfVertex();
		for (int d = 0; d < 3; d++) {
			refPoint[d] = cNode->getParentNode()->getNumOfVertex()*((cNode->getParentNode()->getAvgPoint())[d]);
		}
		for (int i = 0; i < 26; i++) {
			if (cNodeN[i] != NULL && !(cNodeN[i]->isDeleted()) && cNodeN[i]->getNumOfVertex() > 0) {
				for (int d = 0; d < 3; d++)
					refPoint[d] += (cNodeN[i]->getNumOfVertex()*((cNodeN[i]->getAvgPoint())[d]));
				nw += cNodeN[i]->getNumOfVertex();
			}
		}
		for (int d = 0; d < 3; d++)
			refPoint[d] /= nw;
		vector<const double*> tmpPoint;
		vector<double> weight;
		double totalWeight = 0;
		octreeNode* const * currNeighbor = currNode->getSameLevelNeighbors();
		getTwoLevelChildPoints(currNode, refPoint, sideLength, tmpPoint, weight, totalWeight, 3);
		for (int i = 0; i < 26; i++) {
			if (currNeighbor[i] != NULL && !(currNeighbor[i]->isDeleted()) && currNeighbor[i]->getNumOfVertex() > 0) {
				getTwoLevelChildPoints(currNeighbor[i], refPoint, sideLength, tmpPoint, weight, totalWeight, 3);
			}
		}

		int nvsize = tmpPoint.size();
		if (nvsize < 10) {
			keepGoUp = true;
		}
		else {
			gsl_vector *eval = gsl_vector_alloc (dim);
			gsl_matrix *evec = gsl_matrix_alloc (dim, dim);

			solvePCA(tmpPoint, nvsize, weight, totalWeight, eval, evec, dim);

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
			if (fabs(maxEigenValue)/fabs(minEigenValue) < 2 && currNode->getParentNode() != NULL) {
				currNode = currNode->getParentNode();
				keepGoUp = true;
			}
			else if (currNode->getParentNode() == NULL) {
				return false;
			}
			else {
				gsl_vector *smallestEigenVector = gsl_vector_alloc(dim);
				gsl_matrix_get_col (smallestEigenVector, evec, minIndex);
				for (int d = 0; d < dim; d++) {
					order1Func[d] = gsl_vector_get(smallestEigenVector, d);
				}
				order1b = 0;
				for (int i = 0; i < nvsize; i++) {
					for (int d = 0; d < dim; d++)
						order1b = order1b + weight[i]*order1Func[d]*tmpPoint[i][d];
				}
				order1b /= totalWeight;
				order1b = -1*order1b;
				gsl_vector_free(eval);
				gsl_matrix_free(evec);
				keepGoUp = false;
				return true;
			}
		}
		if (keepGoUp)
			cNode = cNode->getParentNode();
	} while (keepGoUp);

	return true;
}

double octreeNode::computeNormalPCA(const vector<double*>& v, const vector<octreeNode*>& neighborhoodVector, int dim) {
	int neighborhoodSize = neighborhoodVector.size();
	vector<const double*> tmpPoint;
	vector<double> weight;
	double totalWeight = 0;
	for (int i = 0; i < neighborhoodSize; i++) {
		tmpPoint.push_back(neighborhoodVector[i]->getAvgPoint());
		weight.push_back(1);
		totalWeight += weight[i];
	}
	gsl_vector *eval = gsl_vector_alloc (dim);
	gsl_matrix *evec = gsl_matrix_alloc (dim, dim);
	solvePCA(tmpPoint, neighborhoodSize, weight, totalWeight, eval, evec, dim);

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

	gsl_vector *smallestEigenVector = gsl_vector_alloc(dim);
	gsl_matrix_get_col (smallestEigenVector, evec, minIndex);
	for (int d = 0; d < dim; d++) {
		order1Func[d] = gsl_vector_get(smallestEigenVector, d);
	}
	order1b = 0;
	for (int i = 0; i < neighborhoodSize; i++) {
		for (int d = 0; d < dim; d++)
			order1b = order1b + weight[i]*order1Func[d]*tmpPoint[i][d];
	}
	order1b /= totalWeight;
	order1b = -1*order1b;

	gsl_vector_free(eval);
	gsl_matrix_free(evec);

	return fabs(minEigenValue)/fabs(maxEigenValue);

}

double octreeNode::computeNormalPCA(vector<const double*>& points, int dim, double* normal) {
	int pointsSize = points.size();
	vector<double> weight;
	double totalWeight = 0;
	for (int i = 0; i < pointsSize; i++) {
		weight.push_back(1);
		totalWeight += weight[i];
	}
	gsl_vector *eval = gsl_vector_alloc (dim);
	gsl_matrix *evec = gsl_matrix_alloc (dim, dim);
	solvePCA(points, pointsSize, weight, totalWeight, eval, evec, dim);

	double minEigenValue = 100;
	int minIndex = 0;
	double maxEigenValue = -1;
	int maxIndex = 0;

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

	gsl_vector *smallestEigenVector = gsl_vector_alloc(dim);
	gsl_matrix_get_col (smallestEigenVector, evec, minIndex);
	for (int d = 0; d < dim; d++) {
		normal[d] = gsl_vector_get(smallestEigenVector, d);
	}

	gsl_vector_free(eval);
	gsl_matrix_free(evec);

	return fabs(minEigenValue)/fabs(maxEigenValue);

}

double* octreeNode::getOrder1Func() {
	return order1Func;
}
double octreeNode::getOrder1b() {
	return order1b;
}

int octreeNode::getNumOfKeepPoint() {
	return numOfKeepPoint;
}
void octreeNode::setNumOfKeepPoint(int n) {
	numOfKeepPoint = n;
}

void octreeNode::removeVertex(int vi, bool kp) {

	for (int i = 0; i < numOfVertex; i++) {
		if (inNodeVertices[i] == vi) {
			inNodeVertices.erase(inNodeVertices.begin() + i, inNodeVertices.begin() + i + 1);
			numOfVertex--;
			if (kp)
				numOfKeepPoint--;
			return ;
		}
	}
}

void octreeNode::laplacianProjection(double* projectedPoint, vector<octreeNode*> neighborNodes, double laplacianFactor) {
	for (int d = 0; d < 3; d++)
		projectedPoint[d] = (getAvgPoint())[d];
	octreeNode* currNode = this;
	double laplacianPoint[3];
	double totalWeight = 0;
	for (int d = 0; d < 3; d++) {
		projectedPoint[d] = (1 - laplacianFactor)*(getAvgPoint())[d];
		laplacianPoint[d] = 0;
	}
	//----
	double maxDist = -1;
	for (int i = 0; i < neighborNodes.size(); i++) {
		double cd = 0;
		for (int d = 0; d < 3; d++)
			cd += (((getAvgPoint())[d] - (neighborNodes[i]->getAvgPoint())[d])*((getAvgPoint())[d] - (neighborNodes[i]->getAvgPoint())[d]));
		if (cd > maxDist)
			maxDist = cd;
	}
	//----
	double avgLength = 0;
	for (int i = 0; i < neighborNodes.size(); i++) {
		double cd = 0;
		for (int d = 0; d < 3; d++)
			cd += (((getAvgPoint())[d] - (neighborNodes[i]->getAvgPoint())[d])*((getAvgPoint())[d] - (neighborNodes[i]->getAvgPoint())[d]));
		avgLength += sqrt(cd);
		double currWeight = exp(-1*cd/((1.0/1.0)*maxDist));
		totalWeight = totalWeight + currWeight;
		for (int d = 0; d < 3; d++)
			laplacianPoint[d] = laplacianPoint[d] + currWeight*(neighborNodes[i]->getAvgPoint())[d];
	}
	avgLength /= neighborNodes.size();
	for (int d = 0; d < 3; d++)
		laplacianPoint[d] = laplacianFactor*laplacianPoint[d]/totalWeight;
	for (int d = 0; d < 3; d++)
		projectedPoint[d] = projectedPoint[d] + laplacianPoint[d];
}

void octreeNode::setNS(double d) {
	ns = d;
}

double octreeNode::getNS() {
	return ns;
}

bool octreeNode::checkLocallyFlat(const vector<const double*>& points, int pointSize, const Vector3d& cp, const Vector3d& normal) {
	return isLocallyFlat(points, pointSize, cp, normal);
}

void octreeNode::printNode(bool* keepPoint, const vector<double*> &v, ofstream & ofs, ofstream & eofs) {
		double corners[8][3];
		for (int i = 0; i < 2 ; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				for (int k = 0; k < 2; k++)
				{
					int ci = 2*2*i + 2*j + k;
					const double* tmpCenter = getCenter();
					corners[ci][0] = tmpCenter[0] + (1 - 2*i)*sideLength/2;
					corners[ci][1] = tmpCenter[1] + (1 - 2*j)*sideLength/2;
					corners[ci][2] = tmpCenter[2] + (1 - 2*k)*sideLength/2;
				}
			}
		}
		eofs << corners[0][0] << " " << corners[0][1] << " " << corners[0][2] << " " << corners[1][0] << " " << corners[1][1] << " " << corners[1][2] << endl;
		eofs << corners[0][0] << " " << corners[0][1] << " " << corners[0][2] << " " << corners[2][0] << " " << corners[2][1] << " " << corners[2][2] << endl;
		eofs << corners[0][0] << " " << corners[0][1] << " " << corners[0][2] << " " << corners[4][0] << " " << corners[4][1] << " " << corners[4][2] << endl;
		eofs << corners[1][0] << " " << corners[1][1] << " " << corners[1][2] << " " << corners[3][0] << " " << corners[3][1] << " " << corners[3][2] << endl;
		eofs << corners[2][0] << " " << corners[2][1] << " " << corners[2][2] << " " << corners[3][0] << " " << corners[3][1] << " " << corners[3][2] << endl;
		eofs << corners[2][0] << " " << corners[2][1] << " " << corners[2][2] << " " << corners[6][0] << " " << corners[6][1] << " " << corners[6][2] << endl;
		eofs << corners[5][0] << " " << corners[5][1] << " " << corners[5][2] << " " << corners[7][0] << " " << corners[7][1] << " " << corners[7][2] << endl;
		eofs << corners[6][0] << " " << corners[6][1] << " " << corners[6][2] << " " << corners[7][0] << " " << corners[7][1] << " " << corners[7][2] << endl;
		eofs << corners[3][0] << " " << corners[3][1] << " " << corners[3][2] << " " << corners[7][0] << " " << corners[7][1] << " " << corners[7][2] << endl;
		eofs << corners[4][0] << " " << corners[4][1] << " " << corners[4][2] << " " << corners[5][0] << " " << corners[5][1] << " " << corners[5][2] << endl;
		eofs << corners[4][0] << " " << corners[4][1] << " " << corners[4][2] << " " << corners[6][0] << " " << corners[6][1] << " " << corners[6][2] << endl;
		eofs << corners[1][0] << " " << corners[1][1] << " " << corners[1][2] << " " << corners[5][0] << " " << corners[5][1] << " " << corners[5][2] << endl;
		const vector<int> & nV = getVertices();
		int nVSize = nV.size();

		for (int i = 0; i < nVSize; i++) {
			int vi = nV[i];
			if (keepPoint[vi])
				ofs << v[vi][0] << " " << v[vi][1] << " " << v[vi][2] << endl;
		}
}

