#include <iostream>
#include <fstream>
#include <queue>
#include <ctime>
#include <sys/time.h>
#include <cmath>
#include "octree.hpp"
#include "globalFunction.hpp"
#include "octreeNode.hpp"
#include <time.h>
#include <vector>
#include <map>
#include <algorithm>
#include <stack>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>

bool readPointSet(char* fileName, vector<double*>& v) {
	if (fileName == NULL || strlen(fileName)==0) return false;
	ifstream ifs(fileName);
	if (ifs.fail()) return false;

	while (!ifs.eof()) {
		double* p = new double[3];
		ifs >> p[0] >> p[1] >> p[2];
		if (!ifs.eof()) {
			v.push_back(p);
		}
	}
	ifs.close();
	return true;
}

void buildOctree(octree* myOctree, vector<double*>& v, bool* keepPoint, int op, double minSideLength) {
	myOctree->buildInitOctree(v, keepPoint);
	queue<octreeNode*> myqueue;
	if (op == 0)
		myOctree->splitAndBalanceTree(v, keepPoint, myqueue);
	else
		myOctree->splitTreeToSameLeafLevel(v, myqueue, op, minSideLength);
}

double squarePointDistance(const double* p1, const double* p2) {
	double cd = 0;
	for (int d = 0; d < 3; d++) {
		cd += ((p1[d] - p2[d])*(p1[d] - p2[d]));
	}
	return cd;
}

bool myFunction(const pair<int, octreeNode*>& i, const pair<int, octreeNode*>& j) {
	return (i.first < j.first);
}

void removeOutliersAndHighNoise(octree* myOctree, vector<octreeNode*>& leafNodeVector, int leafNodeVectorSize, vector<double*>& v, int vsize, bool* keepPoint, double refSideLength, double& avgNumOfNeighbor, double& sdNumOfNeighbor, double sdFactor) {
	avgNumOfNeighbor = 0;
	sdNumOfNeighbor = 0;
	int vLimitSize = (int) floor(leafNodeVectorSize/100.0);
	vector<pair<int,octreeNode*> > numOfNeighbors;
	for (int i = 0; i < leafNodeVectorSize; i++) {
		numOfNeighbors.push_back(make_pair(0,leafNodeVector[i]));
	}
	int countNode = 0;
	for (int i = 0; i < leafNodeVectorSize; i++) {
		if (!(leafNodeVector[i]->isDeleted())) {
			octreeNode* currNode = leafNodeVector[i];
			while (currNode->getSideLength() < 2*refSideLength) {
				currNode = currNode->getParentNode();
			}

			numOfNeighbors[i].first = currNode->getNumOfKeepPoint();
			octreeNode* const * currNeighbor = currNode->getSameLevelNeighbors();
			for (int ni = 0; ni < 26; ni++) {
				if (currNeighbor[ni] != NULL && !(currNeighbor[ni]->isDeleted()) && currNeighbor[ni]->getNumOfVertex() > 0) {
					numOfNeighbors[i].first += (currNeighbor[ni]->getNumOfKeepPoint());
				}
			}
			avgNumOfNeighbor += numOfNeighbors[i].first;
			countNode++;
		}
	}
	avgNumOfNeighbor /= countNode;
	countNode = 0;
	for (int i = 0; i < leafNodeVectorSize; i++) {
		if (!(leafNodeVector[i]->isDeleted())) {
			double tmp = numOfNeighbors[i].first - avgNumOfNeighbor;
			sdNumOfNeighbor += (tmp*tmp);
			countNode++;
		}
	}
	sdNumOfNeighbor = sqrt(sdNumOfNeighbor/countNode);

	if (sdFactor*sdNumOfNeighbor > avgNumOfNeighbor) {
		sort(numOfNeighbors.begin(), numOfNeighbors.end(), myFunction);
		int countRN = 0;
		int numRN = 0;
		while (numRN < vLimitSize) {
			if (!(numOfNeighbors[countRN].second->isDeleted())) {
				numOfNeighbors[countRN].second->setDelete(true);
				octreeNode* tmpNode = numOfNeighbors[countRN].second;
				while (tmpNode != NULL) {
					tmpNode->setNumOfKeepPoint(tmpNode->getNumOfKeepPoint() - numOfNeighbors[countRN].second->getNumOfKeepPoint());
					tmpNode = tmpNode->getParentNode();
				}
				const vector<int> & currV = numOfNeighbors[countRN].second->getVertices();
				for (int j = 0; j < currV.size(); j++) {
					keepPoint[currV[j]] = false;
				}
				numRN++;
			}
			countRN++;
		}
	}
}

void computeNeighborhood(octreeNode* leafNode, vector<double*>& v, bool* keepPoint, double refSideLength, vector<int>& neighborhoodVector) {
	octreeNode* currNode = leafNode;
	while (currNode->getSideLength() < 4*refSideLength) {
		currNode = currNode->getParentNode();
	}
	const vector<int>& cv = currNode->getVertices();
	int cvsize = cv.size();
	for (int ci = 0; ci < cvsize; ci++) {
		if (keepPoint[cv[ci]]) {
			double cd = sqrt(squarePointDistance(leafNode->getAvgPoint(), v[cv[ci]]));
			if (cd < 4.001*refSideLength)
				neighborhoodVector.push_back(cv[ci]);
		}
	}
	octreeNode* const * currNeighbor = currNode->getSameLevelNeighbors();
	for (int ni = 0; ni < 26; ni++) {
		if (currNeighbor[ni] != NULL && !(currNeighbor[ni]->isDeleted()) && currNeighbor[ni]->getNumOfVertex() > 0) {
			const vector<int>& nv = currNeighbor[ni]->getVertices();
			int nvsize = nv.size();
			for (int nvi = 0; nvi < nvsize; nvi++) {
				if (keepPoint[nv[nvi]]) {
					double cd = sqrt(squarePointDistance(leafNode->getAvgPoint(), v[nv[nvi]]));
					if (cd < 4.001*refSideLength)
						neighborhoodVector.push_back(nv[nvi]);
				}
			}
		}
	}
}

void updateAvgPoint(octree* myOctree, const vector<double*>& v, bool* keepPoint) {
	queue<octreeNode*> uQueue;
	uQueue.push(myOctree->getRoot());
	while (!(uQueue.empty())) {
		octreeNode* currNode = uQueue.front();
		uQueue.pop();
		currNode->updateAvgPointWithKeepPoint(v, keepPoint);
	}
}

void getLaplacianNeighbor(octreeNode* currNode, vector<octreeNode*>& LNNeighbor) {
	double myPi = 3.14159265358979323846;
	int s1 = 4;
	int s2 = 6;
	bool intersected[s1][s2];
	double currMinDist[s1][s2];
	octreeNode* currPointNeighbor[s1][s2];
	for (int u = 0; u < s1; u++) {
		for (int v = 0; v < s2; v++) {
			intersected[u][v] = false;
			currPointNeighbor[u][v] = NULL;
			currMinDist[u][v] = -1;
		}
	}
	Vector3d currVector(currNode->getAvgPoint());
	const double* currPoint = currNode->getAvgPoint();
	const double* currCenter = currNode->getCenter();
	double currSideLength = currNode->getSideLength();
	queue<octreeNode*> tQueue;
	tQueue.push(currNode->getParentNode()->getParentNode());
	octreeNode* const * currNeighbor = currNode->getParentNode()->getParentNode()->getSameLevelNeighbors();
	for (int i = 0; i < 26; i++) {
		if (currNeighbor[i] != NULL && !(currNeighbor[i]->isDeleted()) && currNeighbor[i]->getNumOfKeepPoint() > 0) {
			tQueue.push(currNeighbor[i]);
		}
	}
	vector<octreeNode*> neighborNodeVector;
	while (!(tQueue.empty())) {
		octreeNode* tNode = tQueue.front();
		tQueue.pop();
		if (tNode->isLeaf() && tNode != currNode) {
			const double* tCenter = tNode->getCenter();
			if (fabs(tCenter[0] - currCenter[0]) < 4.0001*currSideLength && fabs(tCenter[1] - currCenter[1]) < 4.0001*currSideLength && fabs(tCenter[2] - currCenter[2]) < 4.0001*currSideLength) {
				neighborNodeVector.push_back(tNode);
			}
		}
		else {
			octreeNode* const * tChild = tNode->getChildNodes();
			for (int ci = 0; ci < 8; ci++) {
				if (tChild[ci] != NULL && !(tChild[ci]->isDeleted()) && tChild[ci]->getNumOfKeepPoint() > 0) {
					tQueue.push(tChild[ci]);
				}
			}
		}
	}
	int neighborSize = neighborNodeVector.size();
	for (int q = 0; q < neighborSize; q++) {
		//----
		const double* neighborPoint = neighborNodeVector[q]->getAvgPoint();
		int index1 = -1;
		int index2 = -1;
		int intersectDim = -1;
		double projDist[6];
		bool paralle[3];
		double dimDiff[3];
		for (int d = 0; d < 3; d++) {
			dimDiff[d] = neighborPoint[d] - currPoint[d];
			if (fabs(dimDiff[d]) > 0.0000000001*currSideLength)
				paralle[d] = false;
			else
				paralle[d] = true;
		}
		double minDist = 100;
		double projectedPoint[3];
		for (int d = 0; d < 3 && index2 < 0; d++) {
			if (!paralle[d]) {
				projDist[2*d] = (-1*dimDiff[d] + currSideLength/2)/dimDiff[d];
				projDist[2*d + 1] = (-1*dimDiff[d] - currSideLength/2)/dimDiff[d];
				int candidateIndex = -1;
				if (projDist[2*d] < projDist[2*d+1]) {
					candidateIndex = 2*d;
				}
				else {
					candidateIndex = 2*d + 1;
				}
				for (int id = 0; id < 3; id++)
					projectedPoint[id] = neighborPoint[id] + projDist[candidateIndex]*(dimDiff[id]);
				bool outside = false;
				for (int id = 0; id < 3 && !outside; id++) {
					if (id != d) {
						if (fabs(projectedPoint[id] - currPoint[id]) > currSideLength/2)
							outside = true;
					}
				}
				if (!outside) {
					index2 = candidateIndex;
					intersectDim = d;
				}
			}
		}
		double diff1 = 0;
		double diff2 = 0;
		bool diff1Done = false;
		for (int d = 0; d < 3; d++) {
			if (d != intersectDim) {
				if (!diff1Done) {
					diff1 = projectedPoint[d] - currPoint[d];
					diff1Done = true;
				}
				else {
					diff2 = projectedPoint[d] - currPoint[d];
				}
			}
		}
		if (intersectDim >= 0) {
			if (diff1 >= 0 && diff2 >= 0) {
				index1 = 0;
			}
			else if (diff1 >= 0 && diff2 < 0) {
				index1 = 1;
			}
			else if (diff1 < 0 && diff2 >= 0) {
				index1 = 2;
			}
			else if (diff1 < 0 && diff2 < 0) {
				index1 = 3;
			}
			else {
				cout << "unhandle case occur" << endl;
				cout << diff1 << " " << diff2 << endl;
				cout << dimDiff[intersectDim] << " " << paralle[intersectDim] << endl;
				for (int d = 0; d < 3; d++)
					cout << "paralle " << d << " = " << paralle[d] << endl;
				exit(1);
			}
			double diffLength = sqrt(squarePointDistance(currPoint, neighborPoint));
			if (intersected[index1][index2]) {
				if (diffLength < currMinDist[index1][index2] || currMinDist[index1][index2] < 0) {
					currPointNeighbor[index1][index2] = neighborNodeVector[q];
					currMinDist[index1][index2] = diffLength;
				}
			}
			else {
				intersected[index1][index2] = true;
				currPointNeighbor[index1][index2] = neighborNodeVector[q];
				currMinDist[index1][index2] = diffLength;
			}
		}
		//----

	}
	for (int u = 0; u < s1; u++) {
		for (int v = 0; v < s2; v++) {
			if (intersected[u][v])
				LNNeighbor.push_back(currPointNeighbor[u][v]);
		}
	}
}

bool myFunction_b(const pair<int,int>& i, const pair<int,int>& j) {
	return (i.first > j.first);
}

void BFSRemoveOutliers(vector<double*>& globalV, vector<octreeNode*>& leafNodeVector, int leafNodeVectorSize, bool* keepPoint, vector<vector<double*> >& surfaceV, int numOfSurface) {
	int* componentIndex = new int[leafNodeVectorSize];
	int currIndex = 0;
	vector<pair<int,int> > componentSize;
	for (int i = 0; i < leafNodeVectorSize; i++) {
		if (!(leafNodeVector[i]->isProcessed())) {
			currIndex++;
			queue<octreeNode*> traverseQueue;
			traverseQueue.push(leafNodeVector[i]);
			componentIndex[i] = currIndex;
			leafNodeVector[i]->setProcessed(true);
			int currSize = 1;
			while (!(traverseQueue.empty())) {
				octreeNode* currNode = traverseQueue.front();
				traverseQueue.pop();
				octreeNode* const * neighborNode = currNode->getSameLevelNeighbors();
				for (int ni = 0; ni < 26; ni++) {
					if (neighborNode[ni] != NULL && !(neighborNode[ni]->isDeleted()) && neighborNode[ni]->getNodeIndex() >= 0 && !(neighborNode[ni]->isProcessed()) && neighborNode[ni]->getNumOfVertex() > 0) {
						componentIndex[neighborNode[ni]->getNodeIndex()] = currIndex;
						currSize++;
						traverseQueue.push(neighborNode[ni]);
						neighborNode[ni]->setProcessed(true);
					}
				}
			}
			//---- for multi-surface use ----
			componentSize.push_back(make_pair(currSize, currIndex));
			if (componentSize.size() != currIndex) {
				cout << "Error on currIndex" << endl;
				exit(1);
			}
			//---- end of for multi-surface use ----
		}
	}
	sort(componentSize.begin(), componentSize.end(), myFunction_b);
	vector<int> frontComponent;
	for (int i = 0; i < componentSize.size(); i++) {
		frontComponent.push_back(-1);
	}
	for (int i = 0; i < componentSize.size(); i++) {
		frontComponent[componentSize[i].second] = i;
	}

	for (int i = 0; i < leafNodeVectorSize; i++) {
		leafNodeVector[i]->setProcessed(false);
		if (frontComponent[componentIndex[i]] < numOfSurface) {
			const vector<int> & currV = leafNodeVector[i]->getVertices();
			for (int j = 0; j < currV.size(); j++) {
				surfaceV[frontComponent[componentIndex[i]]].push_back(globalV[currV[j]]);
			}
		}
		else {
			leafNodeVector[i]->setDelete(true);
			octreeNode* tmpNode = leafNodeVector[i];
			while (tmpNode != NULL) {
				tmpNode->setNumOfKeepPoint(tmpNode->getNumOfKeepPoint() - leafNodeVector[i]->getNumOfKeepPoint());
				tmpNode = tmpNode->getParentNode();
			}
			const vector<int> & currV = leafNodeVector[i]->getVertices();
			for (int j = 0; j < currV.size(); j++) {
				keepPoint[currV[j]] = false;
			}
		}
	}
	delete[] componentIndex;
}
void solvePCA(const vector<double*>& tmpPoint, int nvsize, const vector<double>& weight, double totalWeight, gsl_vector* eval, gsl_matrix* evec, int dim) {
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

void updateOctreeToSameLevelLeaf(octree* myOctree, vector<double*>& v, double lFactor, bool* keepPoint, double& avgSideLength, int& countNode) {
	vector<octreeNode*> leafNodeVector;
	myOctree->getLeafNodeVector(leafNodeVector);
	int leafNodeVectorSize = leafNodeVector.size();
	for (int i = 0; i < leafNodeVectorSize; i++) {
		if (!(leafNodeVector[i]->isDeleted())) {
			avgSideLength += (leafNodeVector[i]->getSideLength());
			countNode++;
		}
	}
	avgSideLength /= (countNode);

	for (int i = 0; i < leafNodeVectorSize; i++) {
		if (!(leafNodeVector[i]->isDeleted())) {
			if (leafNodeVector[i]->getSideLength() > lFactor*avgSideLength) {
				queue<octreeNode*> sQueue;
				sQueue.push(leafNodeVector[i]);
				while (!(sQueue.empty())) {
					octreeNode* currNode = sQueue.front();
					sQueue.pop();
					if (currNode->getSideLength() > lFactor*avgSideLength)
						currNode->splitNode(v, keepPoint, sQueue);
				}
			}
			else if (leafNodeVector[i]->getParentNode()->getSideLength() <= lFactor*avgSideLength) {
				octreeNode* newLeafNode = leafNodeVector[i];
				while (newLeafNode->getParentNode()->getSideLength() <= lFactor*avgSideLength) {
					newLeafNode = newLeafNode->getParentNode();
				}
				newLeafNode->setLeaf(true);
				octreeNode* const * ncn = newLeafNode->getChildNodes();
				queue<octreeNode*> rQueue;
				for (int ci = 0; ci < 8; ci++) {
					if (ncn[ci] != NULL) {
						rQueue.push(ncn[ci]);
					}
				}
				while (!(rQueue.empty())) {
					octreeNode* rNode = rQueue.front();
					rQueue.pop();
					rNode->setDelete(true);
					octreeNode* const * rcNode = rNode->getChildNodes();
					for (int ci = 0; ci < 8; ci++) {
						if (rcNode[ci] != NULL)
							rQueue.push(rcNode[ci]);
					}
				}
			}
		}
	}
}

double vectorLength(const double* p) {
	return sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
}

void scalePointSet(vector<double*>& v, int vsize, bool* keepPoint, double* centerPoint, double& scale) {
	double maxCoord[3];
	double minCoord[3];
	for (int d = 0; d < 3; d++) {
		maxCoord[d] = 0;
		minCoord[d] = 0;
	}
	for (int i = 0; i < vsize; i++) {
		if (keepPoint[i]) {
			for (int d = 0; d < 3; d++) {
				if (v[i][d] > maxCoord[d] || i == 0)
					maxCoord[d] = v[i][d];
				if (v[i][d] < minCoord[d] || i == 0)
					minCoord[d] = v[i][d];
			}
		}
	}
	double dimPoint[3];
	for (int d = 0; d < 3; d++) {
		centerPoint[d] = (maxCoord[d] + minCoord[d])/2;
		dimPoint[d] = maxCoord[d] - minCoord[d];
	}

	for (int i = 0; i < vsize; i++) {
		for (int d = 0; d < 3; d++)
			v[i][d] = v[i][d] - centerPoint[d];
	}
	double tmpScale = dimPoint[0];
	for (int d = 1; d < 3; d++) {
		if (tmpScale < dimPoint[d])
			tmpScale = dimPoint[d];
	}
	scale = 2.0/tmpScale;

	for (int i = 0; i < vsize; i++) {
		for (int d = 0; d < 3; d++)
			v[i][d] = v[i][d]*scale;
	}
}

void scaleBack(vector<double*>& v, int vsize, bool* keepPoint, double* centerPoint, double scale) {
	for (int i = 0; i < vsize; i++) {
		for (int d = 0; d < 3; d++)
			v[i][d] = v[i][d]/scale;
	}
	for (int i = 0; i < vsize; i++) {
		for (int d = 0; d < 3; d++)
			v[i][d] = v[i][d] + centerPoint[d];
	}
}

int main(int argc, char* argv[]) {
	if (argc < 11) {
		cout << "usage: main inputPointSetFileName laplacianFactor numOfSurface laplacianStopFactor updateBFSSameLevelFactor updateNormalSameLevelFactor laplacianThershold maxIteration sdFactor outputFileName(without extension)" << endl;
		exit(1);
	}
	double laplacianStopFactor = atof(argv[4]);
	double laplacianFactor = atof(argv[2]);

	timeval start, end, beginning;
	gettimeofday(&beginning, NULL);
	gettimeofday(&start, NULL);
	vector<double*> globalV;
	if (!(readPointSet(argv[1], globalV)))
		return false;
	//----
	int globalVSize = globalV.size();
	bool* globalKeepPoint = new bool[globalVSize];
	for (int i = 0; i < globalVSize; i++)
		globalKeepPoint[i] = true;
	double updateBFSSameLevelFactor = atof(argv[5]);
	double updateNormalSameLevelFactor = atof(argv[6]);
	int numOfSurface = atoi(argv[3]);
	double laplacianThershold = atof(argv[7]);
	int maxIteration = atoi(argv[8]);
	double sdFactor = atof(argv[9]);

	//---- build original octree ----
	gettimeofday(&start, NULL);
	octree* myOctree = new octree;
	buildOctree(myOctree, globalV, globalKeepPoint, 0, 0);
	gettimeofday(&end, NULL);
	cout << "Build and split tree: " << end.tv_sec - start.tv_sec + (end.tv_usec - start.tv_usec)/1000000.0 << " seconds" << endl;
	gettimeofday(&start, NULL);

	//---- remove high noise and outliers ----
	double oldAvgNumOfNeighbor = -100;
	double avgNumOfNeighbor = 0;
	double oldSdNumOfNeighbor = -100;
	double sdNumOfNeighbor = 100;

	double avgSideLength = 0;
	int countNode = 0;
	updateOctreeToSameLevelLeaf(myOctree, globalV, updateBFSSameLevelFactor, globalKeepPoint, avgSideLength, countNode);
	vector<octreeNode*> leafNodeVector;
	myOctree->getLeafNodeVector(leafNodeVector);
	int leafNodeVectorSize = leafNodeVector.size();
	//----
	for (int i = 0; i < leafNodeVectorSize; i++) {
		leafNodeVector[i]->setNodeIndex(i);
	}

	vector<vector<double*> > surfaceV;
	for (int i = 0; i < numOfSurface; i++) {
		surfaceV.push_back(vector<double*>());
	}
	BFSRemoveOutliers(globalV, leafNodeVector, leafNodeVectorSize, globalKeepPoint, surfaceV, numOfSurface);

	vector<double*> currPS;
	double outputSideLength = -1;

	for (int si = 0; si < numOfSurface; si++) {
		cout << "working on surface " << si << endl;
		vector<double*> v;
		int vsize = -1;
		bool* keepPoint = NULL;
		if (numOfSurface > 1) {
			vsize = surfaceV[si].size();
			keepPoint = new bool[vsize];
			for (int i = 0; i < vsize; i++) {
				v.push_back(surfaceV[si][i]);
				keepPoint[i] = true;
			}
			if (myOctree != NULL) {
				delete myOctree;
				myOctree = NULL;
			}
			gettimeofday(&start, NULL);
			myOctree = new octree;
			buildOctree(myOctree, v, keepPoint, 0, 0);
			gettimeofday(&end, NULL);
			cout << "Build and split tree: " << end.tv_sec - start.tv_sec + (end.tv_usec - start.tv_usec)/1000000.0 << " seconds" << endl;
			avgSideLength = 0;
			countNode = 0;
			updateOctreeToSameLevelLeaf(myOctree, v, updateNormalSameLevelFactor, keepPoint, avgSideLength, countNode);
		}
		else {
			v = globalV;
			vsize = globalVSize;
			keepPoint = globalKeepPoint;
		}
		updateAvgPoint(myOctree, v, keepPoint);
		leafNodeVector.clear();
		myOctree->getLeafNodeVector(leafNodeVector);
		leafNodeVectorSize = leafNodeVector.size();
		//---- may not be necessary
		for (int i = 0; i < leafNodeVectorSize; i++) {
			leafNodeVector[i]->setNodeIndex(i);
		}
		//----

		while (sdFactor*sdNumOfNeighbor > avgNumOfNeighbor) {
			oldSdNumOfNeighbor = sdNumOfNeighbor - 1;
			while (fabs(sdNumOfNeighbor - oldSdNumOfNeighbor) > 0.0000001 && sdFactor*sdNumOfNeighbor > avgNumOfNeighbor) {
				oldAvgNumOfNeighbor = avgNumOfNeighbor;
				oldSdNumOfNeighbor = sdNumOfNeighbor;
				avgSideLength = 0;
				countNode = 0;
				for (int i = 0; i < leafNodeVectorSize; i++) {
					if (!(leafNodeVector[i]->isDeleted())) {
						avgSideLength += (leafNodeVector[i]->getSideLength());
						countNode++;
					}
				}
				avgSideLength /= (countNode);
				removeOutliersAndHighNoise(myOctree, leafNodeVector, leafNodeVectorSize, v, vsize, keepPoint, avgSideLength, avgNumOfNeighbor, sdNumOfNeighbor, sdFactor);
			}
		}

		gettimeofday(&end, NULL);
		cout << "Remove outliers and high noise: " << end.tv_sec - start.tv_sec + (end.tv_usec - start.tv_usec)/1000000.0 << " seconds" << endl;
		//---- scale point set ----

		//---- double check keepPoint ----
		for (int i = 0; i < leafNodeVectorSize; i++) {
			if (!(leafNodeVector[i]->isDeleted())) {
				const vector<int>& currV = leafNodeVector[i]->getVertices();
				for (int j = 0; j < currV.size(); j++) {
					if (!keepPoint[currV[j]]) {
						cout << "keepPoint error not delete" << endl;
						cin.get();
					}
				}
			}
			else {
				const vector<int>& currV = leafNodeVector[i]->getVertices();
				for (int j = 0; j < currV.size(); j++) {
					if (keepPoint[currV[j]]) {
						cout << "keepPoint error delete" << endl;
						cin.get();
					}
				}
			}
		}

		double center[3];
		double scale;
		scalePointSet(v, vsize, keepPoint, center, scale);

		if (myOctree != NULL) {
			delete myOctree;
			myOctree = NULL;
		}

		gettimeofday(&start, NULL);
		myOctree = new octree;
		buildOctree(myOctree, v, keepPoint, 0, 0);
		gettimeofday(&end, NULL);
		cout << "Build and split tree: " << end.tv_sec - start.tv_sec + (end.tv_usec - start.tv_usec)/1000000.0 << " seconds" << endl;

		avgSideLength = 0;
		countNode = 0;
		updateOctreeToSameLevelLeaf(myOctree, v, updateNormalSameLevelFactor, keepPoint, avgSideLength, countNode);
		updateAvgPoint(myOctree, v, keepPoint);
		leafNodeVector.clear();
		myOctree->getLeafNodeVector(leafNodeVector);
		leafNodeVectorSize = leafNodeVector.size();
		for (int i = 0; i < leafNodeVectorSize; i++) {
			leafNodeVector[i]->setNodeIndex(i);
		}
		avgSideLength = leafNodeVector[0]->getSideLength();

		//---- obtain one ring neighbor ----
		vector<octreeNode*>* LNNeighbor = new vector<octreeNode*>[leafNodeVectorSize];
		for (int i = 0; i < leafNodeVectorSize; i++) {
			if (!(leafNodeVector[i]->isDeleted())) {
				getLaplacianNeighbor(leafNodeVector[i], LNNeighbor[i]);
				//----
				if (LNNeighbor[i].size() < 2) {
					leafNodeVector[i]->setDelete(true);
					octreeNode* tmpNode = leafNodeVector[i];
					while (tmpNode != NULL) {
						tmpNode->setNumOfKeepPoint(tmpNode->getNumOfKeepPoint() - leafNodeVector[i]->getNumOfKeepPoint());
						tmpNode = tmpNode->getParentNode();
					}
					const vector<int> & currV = leafNodeVector[i]->getVertices();
					for (int j = 0; j < currV.size(); j++) {
						keepPoint[currV[j]] = false;
					}
				}
				//----
			}
		}
		for (int i = 0; i < leafNodeVectorSize; i++) {
			if (!(leafNodeVector[i]->isDeleted())) {
				double ns = leafNodeVector[i]->getNumOfKeepPoint();
				octreeNode* const * neighborNode = leafNodeVector[i]->getSameLevelNeighbors();
				for (int ni = 0; ni < 26; ni++) {
					if (neighborNode[ni] != NULL && !(neighborNode[ni]->isDeleted()) && neighborNode[ni]->getNumOfVertex() > 0) {
						ns += neighborNode[ni]->getNumOfKeepPoint();
					}
				}
				leafNodeVector[i]->setNS(ns);
			}
		}

		//---- meshless laplacian smoothing ----
		vector<double*> projectedPoints;
		for (int i = 0; i < leafNodeVectorSize; i++) {
			projectedPoints.push_back(new double[3]);
			for (int d = 0; d < 3; d++)
				projectedPoints[i][d] = 0;
		}
		int q = 0;
		bool keepMove = true;
		double avgMove = 10;
		vector<bool> hasMove;
		int countLN = 0;
		for (int i = 0; i < leafNodeVectorSize; i++) {
			hasMove.push_back(false);
			if (!(leafNodeVector[i]->isDeleted()))
				countLN++;
		}
		//----
		double globalAvgLength = 0;
		countLN = 0;
		for (int i = 0; i < leafNodeVectorSize; i++) {
			if (!(leafNodeVector[i]->isDeleted())) {
				double localAvgLength = 0;
				for (int j = 0; j < LNNeighbor[i].size(); j++) {
					localAvgLength += sqrt(squarePointDistance(leafNodeVector[i]->getAvgPoint(), LNNeighbor[i][j]->getAvgPoint()));
				}
				localAvgLength /= LNNeighbor[i].size();
				globalAvgLength += localAvgLength;
				countLN++;
			}
		}
		globalAvgLength /= countLN;

		double sa = globalAvgLength*globalAvgLength*countLN/2.0;
		int numOfIteration = (int) floor(sa);

		bool hasDeleted = true;
		while (hasDeleted) {
			hasDeleted = false;
			for (int i = 0; i < leafNodeVectorSize; i++) {
				if (!(leafNodeVector[i]->isDeleted())) {
					bool innerDelete = true;
					while (innerDelete) {
						innerDelete = false;
						for (int j = 0; j < LNNeighbor[i].size() && !innerDelete; j++) {
							if (LNNeighbor[i][j]->isDeleted()) {
								innerDelete = true;
								LNNeighbor[i].erase(LNNeighbor[i].begin() + j);
							}
						}
					}
					if (LNNeighbor[i].size() < 4) {
						hasDeleted = true;
						leafNodeVector[i]->setDelete(true);
						octreeNode* tmpNode = leafNodeVector[i];
						while (tmpNode != NULL) {
							tmpNode->setNumOfKeepPoint(tmpNode->getNumOfKeepPoint() - leafNodeVector[i]->getNumOfKeepPoint());
							tmpNode = tmpNode->getParentNode();
						}
						const vector<int> & currV = leafNodeVector[i]->getVertices();
						for (int j = 0; j < currV.size(); j++) {
							keepPoint[currV[j]] = false;
						}
					}
				}
			}
		}
		int countHasMove = 0;
		while (q < numOfIteration) {
			countHasMove = 0;
			for (int i = 0; i < leafNodeVectorSize; i++) {
				if (!(leafNodeVector[i]->isDeleted())) {
					bool hasMove = false;
					leafNodeVector[i]->laplacianProjection(projectedPoints[i], LNNeighbor[i], laplacianFactor);
					if (hasMove)
						countHasMove++;
				}
			}

			for (int i = 0; i < leafNodeVectorSize; i++) {
				if (!(leafNodeVector[i]->isDeleted())) {
					double moveLength = sqrt(squarePointDistance(leafNodeVector[i]->getAvgPoint(), projectedPoints[i]));
					double avgLength = 0;
					for (int j = 0; j < LNNeighbor[i].size(); j++) {
						avgLength += sqrt(squarePointDistance(leafNodeVector[i]->getAvgPoint(), LNNeighbor[i][j]->getAvgPoint()));
					}
					avgLength /= LNNeighbor[i].size();
					if (moveLength >= avgLength/laplacianStopFactor) {
						countHasMove++;
						leafNodeVector[i]->setAvgPoint(projectedPoints[i]);
					}
					else {
						avgLength = 0;
						moveLength = 0;
					}
				}
			}
			q++;
		}
		cout << "done meshless laplacian smoothing" << endl;

		//---- extract point set from octree ----
		vector<double*> tmpPS;
		for (int i = 0; i < leafNodeVectorSize; i++) {
			if (!(leafNodeVector[i]->isDeleted())) {
				double* p = new double[3];
				for (int d = 0; d < 3; d++) {
					p[d] = (leafNodeVector[i]->getAvgPoint())[d];
				}
				tmpPS.push_back(p);
			}
		}
		vsize = tmpPS.size();
		if (myOctree != NULL) {
			delete myOctree;
			myOctree = NULL;
		}
		if (keepPoint != NULL) {
			delete[] keepPoint;
			keepPoint = NULL;
		}
		keepPoint = new bool[vsize];
		for (int i = 0; i < vsize; i++)
			keepPoint[i] = true;

		gettimeofday(&start, NULL);
		myOctree = new octree;
		buildOctree(myOctree, tmpPS, keepPoint, 0, 0);
		gettimeofday(&end, NULL);
		cout << "Build and split tree: " << end.tv_sec - start.tv_sec + (end.tv_usec - start.tv_usec)/1000000.0 << " seconds" << endl;

		avgSideLength = 0;
		countNode = 0;
		updateOctreeToSameLevelLeaf(myOctree, tmpPS, updateNormalSameLevelFactor, keepPoint, avgSideLength, countNode);
		leafNodeVector.clear();
		myOctree->getLeafNodeVector(leafNodeVector);
		leafNodeVectorSize = leafNodeVector.size();
		for (int i = 0; i < leafNodeVectorSize; i++) {
			leafNodeVector[i]->setNodeIndex(i);
		}
		//---- scale back ----
		if (numOfSurface > 1)
			scaleBack(tmpPS, vsize, keepPoint, center, scale);
		//------------
		for (int i = 0; i < leafNodeVectorSize; i++) {
			if (!(leafNodeVector[i]->isDeleted())) {
				const vector<int>& currVertex = leafNodeVector[i]->getVertices();
				for (int j = 0; j < currVertex.size(); j++) {
					double* p = new double[3];
					for (int d = 0; d < 3; d++) {
						p[d] = tmpPS[currVertex[j]][d];
					}
					currPS.push_back(p);
				}
			}
		}
		tmpPS.clear();
		vsize = currPS.size();
		if (leafNodeVector[0]->getSideLength() > outputSideLength)
			outputSideLength = leafNodeVector[0]->getSideLength();
		if (myOctree != NULL) {
			delete myOctree;
			myOctree = NULL;
		}
		if (keepPoint != NULL) {
			delete[] keepPoint;
			keepPoint = NULL;
		}
	}

	vector<double*> outputPointSet;
	for (int i = 0; i < currPS.size(); i++) {
		outputPointSet.push_back(currPS[i]);
	}
	int vsize = outputPointSet.size();
	if (numOfSurface > 1) {
		double center[3];
		double scale;
		bool* keepPoint = new bool[vsize];
		for (int i = 0; i < vsize; i++)
			keepPoint[i] = true;
		scalePointSet(outputPointSet, vsize, keepPoint, center, scale);
		delete[] keepPoint;
		keepPoint = NULL;
	}

	gettimeofday(&end, NULL);
	cout << "Total meshless denoising time: " << end.tv_sec - beginning.tv_sec + (end.tv_usec - beginning.tv_usec)/1000000.0 << " seconds" << endl;

	//---- output result ----
	cout << "outputting the denoised pointset" << endl;
	char* outputFileName = new char[strlen(argv[10])+14];
	strcpy(outputFileName, argv[10]);
	strcat(outputFileName, "_denoised.obj");
	ofstream ofslpp(outputFileName);
	for (int i = 0; i < vsize; i++) {
		ofslpp << outputPointSet[i][0] << " " << outputPointSet[i][1] << " " << outputPointSet[i][2] << endl;
	}
	ofslpp.close();
	return 0;
}
