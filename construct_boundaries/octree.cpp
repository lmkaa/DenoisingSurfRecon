#include "octree.hpp"
#include <sys/time.h>

octree::octree()
{
	root = NULL;
}

octree::~octree() {
	delete root;
}

octreeNode* const octree::getRoot() const
{
	return root;
}

void octree::setRoot(octreeNode* _root)
{
	root = _root;
}

void octree::buildInitOctree(const vector<double*> &v)
{
	root = new octreeNode();
	int vSize = v.size();
	for (int i = 0; i < vSize; i++)
	{
		root->insertVertex(i, v);
	}
	const double* maxCoord = root->getMaxCoord();
	const double* minCoord = root->getMinCoord();
	root->setRootNode();
}

void octree::splitAndBalanceTree(const vector<double*> &v, queue<octreeNode*> & myqueue)
{
	myqueue.push(root);
	while (!myqueue.empty())
	{
		octreeNode* currNode = myqueue.front();
		myqueue.pop();
		if (currNode->checkSplitNode())
		{
			currNode->splitNode(v, myqueue);
		}
	}
}

void octree::setTrimQueue(octreeNode* n, queue<octreeNode*> & trimQueue)
{
	if (n->isLeaf() && n->getNumOfVertex() > 0)
	{
		trimQueue.push(n);
	}
	else
	{
		octreeNode* const * childNodes = n->getChildNodes();
		for (int i = 0; i < 8; i++)
		{
			if (childNodes[i] != NULL)
				setTrimQueue(childNodes[i], trimQueue);
		}
	}
}

void octree::trimOctree(const vector<double*> &v, vector<double>& neighborDist, vector<double*> &normal, double* myb)
{
	queue<octreeNode*> trimQueue;
	setTrimQueue(root, trimQueue);
	while (!(trimQueue.empty()))
	{
		octreeNode* currNode = trimQueue.front();
		trimQueue.pop();
		if (currNode == root)
		{
			return ;
		}
		else
		{
			if (currNode != NULL && !(currNode->isDeleted()) && currNode->isLeaf() && currNode->getNumOfVertex() > 0)
			{
				currNode->trim(v, neighborDist, normal, myb, trimQueue, root->getMaxCoord(), root->getMinCoord());
			}
		}
	}

}

void octree::locallyExtractPoint(const vector<double*> &v, vector<double*> &normal, double* myb, vector<double*> &finalPointSet, bool* vertexPick)
{
	queue<octreeNode*> myQueue;
	myQueue.push(root);
	int currLevel = 1;
	while (!(myQueue.empty()))
	{
		vector<octreeNode*> currLevelNodes;
		octreeNode* currNode = myQueue.front();
		while (!(myQueue.empty()) && currNode->getLevel() == currLevel)
		{
			currLevelNodes.push_back(currNode);
			myQueue.pop();
			if (!(myQueue.empty()))
				currNode = myQueue.front();
		}
		int clnSize = currLevelNodes.size();
		for (int i = 0; i < clnSize; i++)
		{
			if (currLevelNodes[i]->getNumOfVertex() > 0)
			{
				if (currLevelNodes[i]->isLeaf())
				{
					const double* currCenter = currLevelNodes[i]->getCenter();
					double currSideLength = currLevelNodes[i]->getSideLength();
					octreeNode* const * cnNeighbors = currLevelNodes[i]->getSameLevelNeighbors();
					for (int j = 0; j < 26; j++)
					{
						if (cnNeighbors[j] != NULL && !(cnNeighbors[j]->isDeleted()) && !(cnNeighbors[j]->isLeaf()) && !(cnNeighbors[j]->isProcessed()))
						{
							octreeNode* const * cnNeighborChild = cnNeighbors[j]->getChildNodes();
							for (int k = 0; k < 8; k++)
							{
								if (cnNeighborChild[k] != NULL)
								{
										cnNeighborChild[k]->setLeaf(true);
										cnNeighborChild[k]->setProcessed(true);
								}
							}
						}

						if (cnNeighbors[j] != NULL)
						{
							octreeNode* const * cnnNeighbors = cnNeighbors[j]->getSameLevelNeighbors();
							for (int nj = 0; nj < 26; nj++)
							{
								if (cnnNeighbors[nj] != NULL && !(cnnNeighbors[nj]->isDeleted()) && !(cnnNeighbors[nj]->isLeaf()) && !(cnnNeighbors[nj]->isProcessed()))
								{
									octreeNode* const * cnnNeighborChild = cnnNeighbors[nj]->getChildNodes();
									for (int k = 0; k < 8; k++)
									{
										if (cnnNeighborChild[k] != NULL)
										{
											const double* tmpCenter = cnnNeighborChild[k]->getCenter();
											if (fabs(tmpCenter[0] - currCenter[0]) < ((2 + (1/2000))*currSideLength) &&
												fabs(tmpCenter[1] - currCenter[1]) < ((2 + (1/2000))*currSideLength) &&
												fabs(tmpCenter[2] - currCenter[2]) < ((2 + (1/2000))*currSideLength))
											{
												cnnNeighborChild[k]->setLeaf(true);
												cnnNeighborChild[k]->setProcessed(true);
											}
										}
									}
								}
							}
						}
					}
				}
				else
				{
					octreeNode* const * currChild = currLevelNodes[i]->getChildNodes();
					for (int k = 0; k < 8; k++)
					{
						if (currChild[k] != NULL && !(currChild[k]->isDeleted()))
						{
							myQueue.push(currChild[k]);
						}
					}
				}
			}
		}
		currLevel++;
	}

	queue<octreeNode*> finalQueue;
	finalQueue.push(root);
	while (!(finalQueue.empty()))
	{
		octreeNode* currNode = finalQueue.front();
		finalQueue.pop();

		if (currNode->isLeaf() && currNode->getNumOfVertex() > 0)
		{
			octreeNode* const * currChildNodes = currNode->getChildNodes();
			double currSL = currNode->getSideLength();
			int extractedPoint = 0;
			vector<int> currInsertPoint;
			if (!(currNode->isPickedPoint()))
				currNode->setHasPickPoint(true);
			else
			{
				cout << "this node has pick more than 1 point" << endl;
				cin.get();
			}

			for (int ci = 0; ci < 8; ci++)
			{
				currInsertPoint.push_back(-1);
				if (currChildNodes[ci] != NULL)
				{
					const vector<int> & currChildVertices = currChildNodes[ci]->getVertices();
					if (currChildNodes[ci]->getNumOfVertex() > 0)
					{
						const double* ccCenter = currChildNodes[ci]->getCenter();
						double* x = new double[3];
						int sv = currChildNodes[ci]->getNearestToAvgPoint();
						for (int d = 0; d < 3; d++)
						{
							x[d] = v[sv][d];
						}
						if (!(currChildNodes[ci]->isPickedPoint()))
							currChildNodes[ci]->setHasPickPoint(true);
						else
						{
							cout << "this node has pick more than 1 point" << endl;
							cin.get();
						}
						finalPointSet.push_back(x);
						vertexPick[sv] = true;
						extractedPoint++;
					}
				}
			}
			if (extractedPoint == 0)
			{
				const vector<int> & currVertices = currNode->getVertices();
				const double* cCenter = currNode->getCenter();
				double* x = new double[3];
				int sv = currNode->getNearestToAvgPoint();
				for (int d = 0; d < 3; d++)
				{
					x[d] = v[sv][d];
				}
				finalPointSet.push_back(x);
				vertexPick[sv] = true;
			}
		}
		else if (currNode->getNumOfVertex() > 0)
		{
			octreeNode* const * currChild = currNode->getChildNodes();
			for (int k = 0; k < 8; k++)
			{
				if (currChild[k] != NULL && !(currChild[k]->isDeleted()) && currChild[k]->getNumOfVertex() > 0)
					finalQueue.push(currChild[k]);
			}
		}
	}
}

void octree::locallyExtractPoint_onePointEach(const vector<double*> &v, vector<double*> &finalPointSet)
{
	queue<octreeNode*> myQueue;
	myQueue.push(root);
	int currLevel = 1;
	while (!(myQueue.empty()))
	{
		vector<octreeNode*> currLevelNodes;
		octreeNode* currNode = myQueue.front();
		while (!(myQueue.empty()) && currNode->getLevel() == currLevel)
		{
			currLevelNodes.push_back(currNode);
			myQueue.pop();
			if (!(myQueue.empty()))
				currNode = myQueue.front();
		}
		int clnSize = currLevelNodes.size();
		for (int i = 0; i < clnSize; i++)
		{
			if (currLevelNodes[i]->getNumOfVertex() > 0)
			{
				if (currLevelNodes[i]->isLeaf())
				{
					const double* currCenter = currLevelNodes[i]->getCenter();
					double currSideLength = currLevelNodes[i]->getSideLength();
					octreeNode* const * cnNeighbors = currLevelNodes[i]->getSameLevelNeighbors();
					for (int j = 0; j < 26; j++)
					{
						if (cnNeighbors[j] != NULL && !(cnNeighbors[j]->isDeleted()) && !(cnNeighbors[j]->isLeaf()) && !(cnNeighbors[j]->isProcessed()))
						{
							octreeNode* const * cnNeighborChild = cnNeighbors[j]->getChildNodes();
							for (int k = 0; k < 8; k++)
							{
								if (cnNeighborChild[k] != NULL)
								{
									const double* tmpCenter = cnNeighborChild[k]->getCenter();
									if (fabs(tmpCenter[0] - currCenter[0]) < ((1.5 + (1/2000))*currSideLength) &&
										fabs(tmpCenter[1] - currCenter[1]) < ((1.5 + (1/2000))*currSideLength) &&
										fabs(tmpCenter[2] - currCenter[2]) < ((1.5 + (1/2000))*currSideLength))
									{
										cnNeighborChild[k]->setLeaf(true);
										cnNeighborChild[k]->setProcessed(true);
									}
								}
							}
						}
					}
				}
				else
				{
					octreeNode* const * currChild = currLevelNodes[i]->getChildNodes();
					for (int k = 0; k < 8; k++)
					{
						if (currChild[k] != NULL && !(currChild[k]->isDeleted()))
						{
							myQueue.push(currChild[k]);
						}
					}
				}
			}
		}
		currLevel++;
	}
	queue<octreeNode*> finalQueue;
	finalQueue.push(root);
	while (!(finalQueue.empty()))
	{
		octreeNode* currNode = finalQueue.front();
		finalQueue.pop();

		if (currNode->isLeaf() && currNode->getNumOfVertex() > 0)
		{
			const vector<int> & currVertices = currNode->getVertices();
			double currSL = currNode->getSideLength();
			vector<int> currInsertPoint;
			currInsertPoint.push_back(-1);
			const double* cCenter = currNode->getCenter();
			double differenceValue;
			int cvsize = currNode->getNumOfVertex();
			double minDist = 1000;
			int minPointIndex = -1;
			for (int cvi = 0; cvi < cvsize; cvi++)
			{
				double cLength = 0;
				for (int d = 0; d < 3; d++)
				{
					differenceValue = v[currVertices[cvi]][d] - cCenter[d];
					cLength = differenceValue*differenceValue;
				}
				cLength = sqrt(cLength);
				if (cLength < minDist || minPointIndex == -1)
				{
					minDist = cLength;
					minPointIndex = cvi;
				}
			}
			double* x = new double[3];
			for (int d = 0; d < 3; d++)
			{
				x[d] = v[currVertices[minPointIndex]][d];
			}
			finalPointSet.push_back(x);

		}
		else if (currNode->getNumOfVertex() > 0)
		{
			octreeNode* const * currChild = currNode->getChildNodes();
			for (int k = 0; k < 8; k++)
			{
				if (currChild[k] != NULL && !(currChild[k]->isDeleted()) && currChild[k]->getNumOfVertex() > 0)
					finalQueue.push(currChild[k]);
			}
		}
	}
}

void octree::splitAndBalanceTreeForMesh(const vector<double*> &v, queue<octreeNode*> & myqueue)
{
	myqueue.push(root);
	while (!myqueue.empty())
	{
		octreeNode* currNode = myqueue.front();
		myqueue.pop();
		if (currNode->checkSplitNodeForMesh())
		{
			currNode->splitNode(v, myqueue);
		}
	}
}
bool myfunction(const pair<double, int> & p1, const pair<double, int> & p2) {
	return (p1.first < p2.first);
}

double pointToLineDist(const Vector3d& p, const Vector3d& diffV, const Vector3d& orgP) {
	return (orgP + (p - orgP).Dot(diffV)*diffV - p).L2Norm();
}

Vector3d pointToLine(const Vector3d& p, const Vector3d& diffV, const Vector3d& orgP) {
	return orgP + (p - orgP).Dot(diffV)*diffV;
}

double pointDistVar(const vector<double> & vPointDist, int vSize, double avgVal) {
	double squareSum = 0;
	for (int i = 0; i < vSize; i++) {
		double tmp = vPointDist[i] - avgVal;
		squareSum += (tmp*tmp);
	}
	return (squareSum/vSize);
}

void octree::splitAndBalanceTreeForLaplacian(const vector<double*> &v, queue<octreeNode*> & myqueue, bool* keepPoint, double* pointWeight, vector<octreeNode*>& directQueue)
{
	timeval start, end;
	gettimeofday(&start, NULL);
	int currLevel = root->getLevel();
	bool keepSplit = true;
	myqueue.push(root);
	octreeNode* currNode = myqueue.front();
	myqueue.pop();
	while (keepSplit) {
		keepSplit = false;
		while (currNode->getLevel() == currLevel) {
			currLevel = currNode->getLevel();
			currNode->splitNodeForLaplacian(v, myqueue, keepSplit);
			currNode = myqueue.front();
			myqueue.pop();
		}
		if (keepSplit) {
			currLevel = currNode->getLevel();
		}
	}
	gettimeofday(&end, NULL);
	cout << "Split new octree: " << end.tv_sec - start.tv_sec + (end.tv_usec - start.tv_usec)/1000000.0 << " seconds" << endl;

	gettimeofday(&start, NULL);
	avgNeighborDist = 0;
	int numOfCount = 0;
	queue<octreeNode*> checkQueue;
	checkQueue.push(root);
	vector<double> pointDist;
	while (!(checkQueue.empty())) {
		octreeNode* currNode = checkQueue.front();
		checkQueue.pop();
		if (currNode->isLeaf()) {
			directQueue.push_back(currNode);
			double neighborDist = computeSixthNeighborDistance(v, v[(currNode->getVertices())[0]], currNode);
			currNode->setSixthNeighborDist(neighborDist);
			avgNeighborDist += neighborDist;
			pointDist.push_back(neighborDist);
			numOfCount++;
		}
		else if (currNode->getNumOfVertex() > 0) {
			octreeNode* const * currChild = currNode->getChildNodes();
			for (int k = 0; k < 8; k++)	{
				if (currChild[k] != NULL && !(currChild[k]->isDeleted()) && currChild[k]->getNumOfVertex() > 0)
					checkQueue.push(currChild[k]);
			}
		}
	}
	gettimeofday(&end, NULL);
	cout << "Compute sixth neighbor distance: " << end.tv_sec - start.tv_sec + (end.tv_usec - start.tv_usec)/1000000.0 << " seconds" << endl;
	gettimeofday(&start, NULL);
	avgNeighborDist = avgNeighborDist/((double) numOfCount);
	double varNeighborDist = pointDistVar(pointDist, numOfCount, avgNeighborDist);
	double sdNeighborDist = sqrt(varNeighborDist);
	int leafNodeSize = directQueue.size();
	for (int i = 0; i < leafNodeSize; i++) {
		octreeNode* currNode = directQueue[i];
		if (!(currNode->isDeleted())) {
			double neighborDist = currNode->getSixthNeighborDist();
		}
		else {
			keepPoint[(currNode->getVertices())[0]] = false;
		}
	}
	gettimeofday(&end, NULL);
	cout << "remove some outliers: " << end.tv_sec - start.tv_sec + (end.tv_usec - start.tv_usec)/1000000.0 << " seconds" << endl;
	gettimeofday(&start, NULL);
	vector<pair<double,int> >* neighborPoint = new vector<pair<double,int> >[leafNodeSize];
	for (int i = 0; i < leafNodeSize; i++) {
		octreeNode* currNode = directQueue[i];
		if (currNode->isLeaf() && currNode->getNumOfVertex() > 0 && !(currNode->isDeleted())) {
			int myPointIndex = (currNode->getVertices())[0];
			double neighborDist = currNode->getSixthNeighborDist();
			octreeNode* rangeNode = currNode;
			while (rangeNode->getSideLength() < 3*neighborDist)
				rangeNode = rangeNode->getParentNode();
			const vector<int> & cv = rangeNode->getVertices();
			int cvSize = cv.size();
			for (int q = 0; q < cvSize; q++) {
				if (cv[q] != myPointIndex && keepPoint[cv[q]]) {
					double currDist = 0;
					for (int d = 0; d < 3; d++)
						currDist += (v[myPointIndex][d] - v[cv[q]][d])*(v[myPointIndex][d] - v[cv[q]][d]);
					currDist = sqrt(currDist);
					neighborPoint[i].push_back(make_pair(currDist, cv[q]));
					pointWeight[myPointIndex] = pointWeight[myPointIndex] + 1;
				}
			}
			octreeNode* const * neighborNodes = rangeNode->getSameLevelNeighbors();
			for (int cn = 0; cn < 26; cn++) {
				if (neighborNodes[cn] != NULL && neighborNodes[cn]->getNumOfVertex() > 0) {
					int nvSize = neighborNodes[cn]->getNumOfVertex();
					const vector<int> & nv = neighborNodes[cn]->getVertices();
					for (int j = 0; j < nvSize; j++) {
						if (keepPoint[nv[j]]) {
							double currDist = 0;
							for (int d = 0; d < 3; d++)
								currDist += (v[myPointIndex][d] - v[nv[j]][d])*(v[myPointIndex][d] - v[nv[j]][d]);
							currDist = sqrt(currDist);
							if (currDist < 3*neighborDist) {
								neighborPoint[i].push_back(make_pair(currDist, nv[j]));
								pointWeight[myPointIndex] = pointWeight[myPointIndex] + 1;
							}
						}
					}
				}
			}
		}
	}
	gettimeofday(&end, NULL);
	cout << "get neighborhood: " << end.tv_sec - start.tv_sec + (end.tv_usec - start.tv_usec)/1000000.0 << " seconds" << endl;

//--------------
	gettimeofday(&start, NULL);
	double avgNP = 0;
	double sdNP = 0;
	for (int i = 0; i < leafNodeSize; i++) {
		octreeNode* currNode = directQueue[i];
		int numOfCount = 0;
		if (currNode->isLeaf() && currNode->getNumOfVertex() > 0 && !(currNode->isDeleted())) {
			avgNP += (neighborPoint[i].size());
			numOfCount++;
		}
	}
	avgNP = avgNP/numOfCount;
	for (int i = 0; i < leafNodeSize; i++) {
		octreeNode* currNode = directQueue[i];
		if (currNode->isLeaf() && currNode->getNumOfVertex() > 0 && !(currNode->isDeleted())) {
			double tmp = neighborPoint[i].size() - avgNP;
			sdNP += (tmp*tmp);
		}
	}
	sdNP = sqrt(sdNP/numOfCount);
	cout << "avgNP = " << avgNP << " sdNP = " << sdNP << endl;
	for (int i = 0; i < leafNodeSize; i++) {
		octreeNode* currNode = directQueue[i];
		if (currNode->isLeaf() && currNode->getNumOfVertex() > 0 && !(currNode->isDeleted())) {
			int myPointIndex = (currNode->getVertices())[0];
			if (neighborPoint[i].size() > (avgNP + 0.5*sdNP)) {
				keepPoint[myPointIndex] = false;
				currNode->setDelete(true);
			}
		}
	}

	gettimeofday(&end, NULL);
	cout << "remove outliers by num of neighbor: " << end.tv_sec - start.tv_sec + (end.tv_usec - start.tv_usec)/1000000.0 << " seconds" << endl;
//--------------


	gettimeofday(&start, NULL);

	double myPi = 3.14159265358979323846;
	for (int i = 0; i < leafNodeSize; i++) {
		octreeNode* currNode = directQueue[i];
		if (currNode->isLeaf() && currNode->getNumOfVertex() > 0 && !(currNode->isDeleted())) {
			int myPointIndex = (currNode->getVertices())[0];
			int neighborSize = neighborPoint[i].size();
			int s1 = 8;
			int s2 = 4;
			bool intersected[s1][s2];
			double currMinDist[s1][s2];
			int currPointIndex[s1][s2];
			for (int u = 0; u < s1; u++) {
				for (int v = 0; v < s2; v++) {
					intersected[u][v] = false;
					currPointIndex[u][v] = -1;
					currMinDist[u][v] = -1;
				}
			}
			Vector3d currVector(v[myPointIndex]);
			for (int q = 0; q < neighborSize; q++) {
				if (keepPoint[(neighborPoint[i])[q].second]) {
					int neighborIndex = (neighborPoint[i])[q].second;
					Vector3d neighborVector(v[neighborIndex]);
					Vector3d diffVector = neighborVector - currVector;
					double diffLength = diffVector.L2Norm();
					double xzvLength = sqrt(diffVector.X()*diffVector.X() + diffVector.Z()*diffVector.Z());
					double angleX;
					if (diffVector.Z() < 0)
						angleX = 2*myPi - acos(diffVector.X()/xzvLength);
					else
						angleX = acos(diffVector.X()/xzvLength);
					diffVector.normalize();
					double angleY = acos(diffVector.Y());
					int indexX;
					if (isnan(angleX) && fabs(xzvLength) < 0.000000001)
						indexX = 0;
					else
						indexX = (int) floor(angleX/(myPi/s2));
					int indexY = (int) floor(angleY/(myPi/s2));
					if (indexX == s1)
						indexX--;
					if (indexY == s2)
						indexY--;
					if (intersected[indexX][indexY]) {
						if (diffLength < currMinDist[indexX][indexY]) {
							currPointIndex[indexX][indexY] = neighborIndex;
							currMinDist[indexX][indexY] = diffLength;
						}
					}
					else {
						intersected[indexX][indexY] = true;
						currPointIndex[indexX][indexY] = neighborIndex;
						currMinDist[indexX][indexY] = diffLength;
					}
				}
			}
			for (int u = 0; u < s1; u++) {
				for (int v = 0; v < s2; v++) {
					if (intersected[u][v])
						currNode->insertLN(currPointIndex[u][v]);
				}
			}
			const vector<int> & currLN = currNode->getLN();
			int lnSize = currLN.size();
			double avgLNDist = 0;
			vector<double> LNDist;
			for (int k = 0; k < lnSize; k++) {
				(currNode->deleteLapN).push_back(false);
				double currDist = 0;
				for (int d = 0; d < 3; d++) {
					currDist += ((v[myPointIndex][d] - v[currLN[k]][d])*(v[myPointIndex][d] - v[currLN[k]][d]));
				}
				currDist = sqrt(currDist);
				LNDist.push_back(currDist);
				avgLNDist += currDist;
			}
			int countDel = 0;
			avgLNDist = avgLNDist/((double) lnSize);
			double varLNDist = pointDistVar(LNDist, lnSize, avgLNDist);
			double sdLNDist = sqrt(varLNDist);
			for (int k = 0; k < lnSize; k++) {
				if (LNDist[k] > (avgLNDist + 200*sdLNDist)) {
					countDel++;
					(currNode->deleteLapN)[k] = true;
				}
			}

			if (lnSize - countDel < 3) {
				currNode->setDelete(true);
				keepPoint[myPointIndex] = false;
			}
		}
	}

	gettimeofday(&end, NULL);
	cout << "get Laplacian neighborhood: " << end.tv_sec - start.tv_sec + (end.tv_usec - start.tv_usec)/1000000.0 << " seconds" << endl;
	gettimeofday(&start, NULL);
	delete[] neighborPoint;
}

void octree::postDel(const vector<double*> &v, const vector<octreeNode*>& directQueue, bool* keepPoint) {
	int leafNodeSize = directQueue.size();
	double* maxLNDist = new double[leafNodeSize];
	bool* keepLNDist = new bool[leafNodeSize];
	int numOfCount = 0;
	double avgMaxDist = 0;
	for (int i = 0; i < leafNodeSize; i++) {
		maxLNDist[i] = 0;
		keepLNDist[i] = false;
		octreeNode* currNode = directQueue[i];
		if (currNode->isLeaf() && currNode->getNumOfVertex() > 0 && !(currNode->isDeleted())) {
			int myPointIndex = (currNode->getVertices())[0];
			if (keepPoint[myPointIndex]) {
				int countNumOfLNPoint = 0;
				const vector<int> & currLN = currNode->getLN();
				int lnSize = currLN.size();
				vector<double> LNDist;
				for (int k = 0; k < lnSize; k++) {
					if (!((currNode->deleteLapN)[k])) {
						double currDist = 0;
						for (int d = 0; d < 3; d++) {
							currDist += ((v[myPointIndex][d] - v[currLN[k]][d])*(v[myPointIndex][d] - v[currLN[k]][d]));
						}
						currDist = sqrt(currDist);

						if (currDist > maxLNDist[i]) {
							maxLNDist[i] = currDist;
						}
					}
				}

				if (maxLNDist[i] > 0) {
					avgMaxDist += maxLNDist[i];
					keepLNDist[i] = true;
					numOfCount++;
				}
			}
		}
	}
	avgMaxDist = avgMaxDist/numOfCount;
	double squareSum = 0;
	for (int i = 0; i < leafNodeSize; i++) {
		if (keepLNDist[i]) {
			double tmp = maxLNDist[i] - avgMaxDist;
			squareSum += (tmp*tmp);
		}
	}
	double sd = sqrt(squareSum/numOfCount);
	cout << "avgMaxDist = " << avgMaxDist << " sd = " << sd << endl;

	for (int i = 0; i < leafNodeSize; i++) {
		octreeNode* currNode = directQueue[i];
		if (currNode->isLeaf() && currNode->getNumOfVertex() > 0 && !(currNode->isDeleted())) {
			int myPointIndex = (currNode->getVertices())[0];
			if (keepPoint[myPointIndex]) {
				if (maxLNDist[i] > (avgMaxDist + 4*sd)) {
					currNode->setDelete(true);
					keepPoint[myPointIndex] = false;
				}
			}
		}
	}

	delete[] maxLNDist;
	delete[] keepLNDist;
}

void octree::buildInitOctreeForMesh(const vector<double*> &v)
{
	root = new octreeNode();
	int vSize = v.size();
	for (int i = 0; i < vSize; i++)
	{
		root->insertVertex(i, v);
	}
	const double* maxCoord = root->getMaxCoord();
	const double* minCoord = root->getMinCoord();
	root->setRootNodeForMesh();
}

octreeNode* octree::findLeafNode(const double* v) const {
	octreeNode* currNode = root;
	while (!currNode->isLeaf()) {
		octreeNode* const * childNodes = currNode->getChildNodes();
		double currSideLength = currNode->getSideLength();
		double childSideLength = currSideLength/2;
		double childHalfSideLength = childSideLength/2;
		int i = 0;
		bool finished = false;
		while (i < 8 && !finished) {
			if (childNodes[i] != NULL) {
				const double* childCenter = childNodes[i]->getCenter();
				if (fabs(v[0] - childCenter[0]) <= childHalfSideLength + 0.000000001 && fabs(v[1] - childCenter[1]) <= childHalfSideLength + 0.000000001 && fabs(v[2] - childCenter[2]) <= childHalfSideLength + 0.000000001) {
					currNode = childNodes[i];
					finished = true;
					break;
				}
			}
			i++;
		}
	}
	return currNode;
}

double octree::computeSixthNeighborDistance(const vector<double*> &v, const double* currPoint) const {
	vector<int> neighborPoints;
	octreeNode* currNode = findLeafNode(currPoint);
	Vector3d currP(currPoint[0],currPoint[1],currPoint[2]);
	int nodeSize = currNode->getNumOfVertex();
	while (nodeSize < 6 && currNode->getParentNode() != NULL) {
		currNode = currNode->getParentNode();
		nodeSize = currNode->getNumOfVertex();
	}
	if (nodeSize > 0) {
		const vector<int>& nodeVertices = currNode->getVertices();
		for (int i = 0; i < nodeSize; i++) {
			neighborPoints.push_back(nodeVertices[i]);
		}
	}
	octreeNode* const * neighborNodes = currNode->getSameLevelNeighbors();
	for (int i = 0; i < 26; i++) {
		if (neighborNodes[i] != NULL && neighborNodes[i]->getNumOfVertex() > 0) {
			nodeSize = neighborNodes[i]->getNumOfVertex();
			const vector<int>& nodeVertices = neighborNodes[i]->getVertices();
			for (int j = 0; j < nodeSize; j++) {
				neighborPoints.push_back(nodeVertices[j]);
			}
		}
	}

	vector<double> neighborDist;
	for (int i = 0; i < neighborPoints.size(); i++) {
		double currDist = (currP-Vector3d(v[neighborPoints[i]])).L2Norm();
		neighborDist.push_back(currDist);
	}
	sort(neighborDist.begin(), neighborDist.end());
	return neighborDist[5];
}

double octree::computeSixthNeighborDistance(const vector<double*> &v, const double* currPoint, octreeNode* currNode) const {
	vector<int> neighborPoints;
	Vector3d currP(currPoint[0],currPoint[1],currPoint[2]);
	int nodeSize = currNode->getNumOfVertex();
	while (nodeSize < 6 && currNode->getParentNode() != NULL) {
		currNode = currNode->getParentNode();
		nodeSize = currNode->getNumOfVertex();
	}
	if (nodeSize > 0) {
		const vector<int>& nodeVertices = currNode->getVertices();
		for (int i = 0; i < nodeSize; i++) {
			neighborPoints.push_back(nodeVertices[i]);
		}
	}
	octreeNode* const * neighborNodes = currNode->getSameLevelNeighbors();
	for (int i = 0; i < 26; i++) {
		if (neighborNodes[i] != NULL && neighborNodes[i]->getNumOfVertex() > 0) {
			nodeSize = neighborNodes[i]->getNumOfVertex();
			const vector<int>& nodeVertices = neighborNodes[i]->getVertices();
			for (int j = 0; j < nodeSize; j++) {
				neighborPoints.push_back(nodeVertices[j]);
			}
		}
	}
	vector<double> neighborDist;
	for (int i = 0; i < neighborPoints.size(); i++) {
		double currDist = (currP-Vector3d(v[neighborPoints[i]])).L2Norm();
		neighborDist.push_back(currDist);
	}
	sort(neighborDist.begin(), neighborDist.end());
	return neighborDist[5];
}

void octree::insertVertexForMesh(int vIndex, const double* v) {
	queue<octreeNode*> myQueue;
	myQueue.push(root);
	octreeNode* currNode = root;
	while (!(currNode->isLeaf())) {
		const double* currCenter = currNode->getCenter();
		octreeNode* const * childNodes = currNode->getChildNodes();
		if (v[0] >= currCenter[0] && v[1] >= currCenter[1] && v[2] >= currCenter[2])
			currNode = childNodes[0];
		else if (v[0] >= currCenter[0] && v[1] >= currCenter[1] && v[2] < currCenter[2])
			currNode = childNodes[1];
		else if (v[0] >= currCenter[0] && v[1] < currCenter[1] && v[2] >= currCenter[2])
			currNode = childNodes[2];
		else if (v[0] >= currCenter[0] && v[1] < currCenter[1] && v[2] < currCenter[2])
			currNode = childNodes[3];
		else if (v[0] < currCenter[0] && v[1] >= currCenter[1] && v[2] >= currCenter[2])
			currNode = childNodes[4];
		else if (v[0] < currCenter[0] && v[1] >= currCenter[1] && v[2] < currCenter[2])
			currNode = childNodes[5];
		else if (v[0] < currCenter[0] && v[1] < currCenter[1] && v[2] >= currCenter[2])
			currNode = childNodes[6];
		else if (v[0] < currCenter[0] && v[1] < currCenter[1] && v[2] < currCenter[2])
			currNode = childNodes[7];
		else {
			cout << "A vertex cannot fall into any child node" << endl;
			exit(0);
		}
	}
	currNode->insertVertexForMesh(vIndex);
}

void octree::computeAvgSideLength() {
	queue<octreeNode*> myQueue;
	myQueue.push(root);
	avgSideLength = 0;
	int numOfLeafNode = 0;
	while (!(myQueue.empty())) {
		octreeNode* currNode = myQueue.front();
		myQueue.pop();
		if (currNode->isLeaf() && currNode->getNumOfVertex() > 0) {
			avgSideLength = avgSideLength + currNode->getSideLength();
			numOfLeafNode++;
		}
		else if (currNode->getNumOfVertex() > 0) {
			octreeNode* const * currChild = currNode->getChildNodes();
			for (int k = 0; k < 8; k++) {
				if (currChild[k] != NULL && !(currChild[k]->isDeleted()) && currChild[k]->getNumOfVertex() > 0)
					myQueue.push(currChild[k]);
			}
		}
	}
	avgSideLength = avgSideLength/numOfLeafNode;
}

double octree::getAvgSideLength() const {
	return avgSideLength;
}

void octree::myProjection(const vector<double*>& v) {
	cout << "start of my projection" << endl;
	cout << "computing neighborhood" << endl;
	computeAvgSideLength();
	queue<octreeNode*> finalQueue;
	finalQueue.push(root);
	while (!(finalQueue.empty())) {
		octreeNode* currNode = finalQueue.front();
		finalQueue.pop();

		if (currNode->isLeaf() && currNode->getNumOfVertex() > 0) {
			octreeNode* rangeNode = currNode;
			while (rangeNode->getSideLength() < avgSideLength)
				rangeNode = rangeNode->getParentNode();
			rangeNode = rangeNode->getParentNode();
			int sv = currNode->getNearestToAvgPoint();
			Vector3d cp(v[sv]);
			currNode->insertSvNeighbor(sv);
			octreeNode* const * cnNeighbors = rangeNode->getSameLevelNeighbors();

			int cc = 0;
			queue<octreeNode*> neighborQueue;
			neighborQueue.push(rangeNode);
			for (int j = 0; j < 26; j++) {
				if (cnNeighbors[j] != NULL && cnNeighbors[j]->getNumOfVertex() > 0) {
					neighborQueue.push(cnNeighbors[j]);
				}
			}
			while (!(neighborQueue.empty())) {
				octreeNode* currNeighbor = neighborQueue.front();
				neighborQueue.pop();
				if (currNeighbor->isLeaf() && currNeighbor->getNumOfVertex() > 0) {
					int nsv = currNeighbor->getNearestToAvgPoint();
					Vector3d np(v[nsv]);
					if ((cp - np).L2Norm() < 3*avgSideLength) {
						currNode->insertSvNeighbor(nsv);
						cc++;
					}
				}
				else if (currNeighbor->getNumOfVertex() > 0) {
					octreeNode* const * currChild = currNeighbor->getChildNodes();
					for (int k = 0; k < 8; k++) {
						if (currChild[k] != NULL && !(currChild[k]->isDeleted()) && currChild[k]->getNumOfVertex() > 0)
							neighborQueue.push(currChild[k]);
					}
				}
			}
			currNode->computeAvgDist(v);
			if (cc < 7)
				currNode->setDelete(true);

		}
		else if (currNode->getNumOfVertex() > 0) {
			octreeNode* const * currChild = currNode->getChildNodes();
			for (int k = 0; k < 8; k++)
			{
				if (currChild[k] != NULL && !(currChild[k]->isDeleted()) && currChild[k]->getNumOfVertex() > 0)
					finalQueue.push(currChild[k]);
			}
		}
	}

	for (int q = 0; q < 100; q++) {
		cout << "my iteration " << q << endl;
		double laplacianFactor = 0.1;
		queue<octreeNode*> finalQueueq;
		finalQueueq.push(root);
		while (!(finalQueueq.empty())) {
			octreeNode* currNode = finalQueueq.front();
			finalQueueq.pop();

			if (currNode->isLeaf() && currNode->getNumOfVertex() > 0 && !(currNode->isDeleted())) {
				// Laplacian smoothing

				const vector<int>& svNeighbor = currNode->getSvNeighbor();
				int svSize = svNeighbor.size();
				double projectedPoint[3];
				double laplacianPoint[3];
				double totalWeight = 0;
				int sv = currNode->getNearestToAvgPoint();
				for (int d = 0; d < 3; d++) {
					projectedPoint[d] = (1 - laplacianFactor)*v[sv][d];
					laplacianPoint[d] = 0;
				}
				double avgDist = currNode->getAvgDist();
				for (int j = 1; j < svSize; j++) {
					int nsv = svNeighbor[j];
					double currDist = (v[sv][0] - v[nsv][0])*(v[sv][0] - v[nsv][0]) + (v[sv][1] - v[nsv][1])*(v[sv][1] - v[nsv][1]) + (v[sv][2] - v[nsv][2])*(v[sv][2] - v[nsv][2]);
					double currWeight = 1/(avgSideLength/2 + currDist);
					totalWeight = totalWeight + currWeight;
					for (int d = 0; d < 3; d++)
						laplacianPoint[d] = laplacianPoint[d] + currWeight*v[nsv][d];
				}

				for (int d = 0; d < 3; d++)
					laplacianPoint[d] = laplacianFactor*laplacianPoint[d]/totalWeight;
				for (int d = 0; d < 3; d++)
					v[sv][d] = projectedPoint[d] + laplacianPoint[d];
			}
			else if (currNode->getNumOfVertex() > 0) {
				octreeNode* const * currChild = currNode->getChildNodes();
				for (int k = 0; k < 8; k++)
				{
					if (currChild[k] != NULL && !(currChild[k]->isDeleted()) && currChild[k]->getNumOfVertex() > 0)
						finalQueueq.push(currChild[k]);
				}
			}
		}
	}
}

void octree::planeProjection(const vector<double*>& v) {
	cout << "start plane projection" << endl;
	int vsize = v.size();
	vector<double*> projectedPoints;
	for (int i = 0; i < vsize; i++) {
		projectedPoints.push_back(new double[3]);
	}
	queue<octreeNode*> finalQueue;
	finalQueue.push(root);
	while (!(finalQueue.empty())) {
		octreeNode* currNode = finalQueue.front();
		finalQueue.pop();
		if (currNode->isLeaf() && currNode->getNumOfVertex() > 0) {
			octreeNode* rangeNode = currNode;
			double normal[3];
			double b;
			while (rangeNode->getSideLength() < 4*avgSideLength)
				rangeNode = rangeNode->getParentNode();
			bool smallEigenRadio = true;
			int sv = currNode->getNearestToAvgPoint();
			int cc = 0;

			rangeNode = rangeNode->getParentNode();
			cc = 0;
			vector<double*> tmpPoint;
			Vector3d cp(v[sv]);

			octreeNode* const * cnNeighbors = rangeNode->getSameLevelNeighbors();
			queue<octreeNode*> neighborQueue;
			neighborQueue.push(rangeNode);
			for (int j = 0; j < 26; j++) {
				if (cnNeighbors[j] != NULL && cnNeighbors[j]->getNumOfVertex() > 0) {
					neighborQueue.push(cnNeighbors[j]);
				}
			}
			while (!(neighborQueue.empty())) {
				octreeNode* currNeighbor = neighborQueue.front();
				neighborQueue.pop();
				if (currNeighbor->isLeaf() && currNeighbor->getNumOfVertex() > 0) {
					int nsv = currNeighbor->getNearestToAvgPoint();
					Vector3d np(v[nsv]);
					if ((cp - np).L2Norm() < 6*avgSideLength) {
						tmpPoint.push_back(currNeighbor->getAvgPoint());
						cc++;
					}
				}
				else if (currNeighbor->getNumOfVertex() > 0) {
					octreeNode* const * currChild = currNeighbor->getChildNodes();
					for (int k = 0; k < 8; k++) {
						if (currChild[k] != NULL && !(currChild[k]->isDeleted()) && currChild[k]->getNumOfVertex() > 0)
							neighborQueue.push(currChild[k]);
					}
				}
			}

			if (cc < 7 || rangeNode == NULL)
				currNode->setDelete(true);

			if (!(currNode->isDeleted())) {
				computeNormal(tmpPoint, normal, b);

				Vector3d currNormal(normal);
				Vector3d currPoint(v[sv]);
				double myt = (-1*(currNormal.Dot(currPoint) + b))/(currNormal.Dot(currNormal));
				for (int d = 0; d < 3; d++)
					projectedPoints[sv][d] = v[sv][d] + myt*normal[d];
				Vector3d projectedP(projectedPoints[sv]);
				if (fabs(projectedP.Dot(currNormal) + b) > 0.0000000001) {
					cout << "a point not on plane" << endl;
					cin.get();
				}
			}
		}
		else if (currNode->getNumOfVertex() > 0) {
			octreeNode* const * currChild = currNode->getChildNodes();
			for (int k = 0; k < 8; k++)
			{
				if (currChild[k] != NULL && !(currChild[k]->isDeleted()) && currChild[k]->getNumOfVertex() > 0)
					finalQueue.push(currChild[k]);
			}
		}
	}
	for (int i = 0; i < vsize; i++) {
		for (int d = 0; d < 3; d++)
			v[i][d] = projectedPoints[i][d];
	}
}

void simpleComputeNormal(const vector<double*> &tmpPoint, int nvsize, double* normal)
{
	double myPi = 3.14159265358979323846;
	double centroid[3];
	for (int d = 0; d < 3; d++) {
		centroid[d] = 0;
		for (int i = 0; i < nvsize; i++)
			centroid[d] = centroid[d] + tmpPoint[i][d];
		centroid[d] = centroid[d]/nvsize;
	}
	double cov[3][3];
	for (int d1 = 0; d1 < 3; d1++) {
		for (int d2 = 0; d2 < 3; d2++) {
			cov[d1][d2] = 0;
			for (int i = 0; i < nvsize; i++)
				cov[d1][d2] = cov[d1][d2] + (tmpPoint[i][d1] - centroid[d1])*(tmpPoint[i][d2] - centroid[d2]);
		}
	}
	gsl_matrix * m = gsl_matrix_alloc(3, 3);
	for (int d1 = 0; d1 < 3; d1++) {
		for (int d2 = 0; d2 < 3; d2++)
			gsl_matrix_set (m, d1, d2, cov[d1][d2]);
	}

	gsl_vector *eval = gsl_vector_alloc (3);
	gsl_matrix *evec = gsl_matrix_alloc (3, 3);
	gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (3);
	gsl_eigen_symmv (m, eval, evec, w);
	gsl_eigen_symmv_free (w);
	double minEigenValue = 100;
	int minIndex = 0;
	for (int d = 0; d < 3; d++) {
		if (d == 0 || gsl_vector_get (eval, d) < minEigenValue) {
			minEigenValue = gsl_vector_get (eval, d);
			minIndex = d;
		}
	}

	gsl_vector *smallestEigenVector = gsl_vector_alloc(3);
	gsl_matrix_get_col (smallestEigenVector, evec, minIndex);
	for (int d = 0; d < 3; d++)
		normal[d] = gsl_vector_get(smallestEigenVector, d);

	gsl_vector_free(eval);
	gsl_matrix_free(evec);
	gsl_matrix_free(m);

}

void octree::laplacianProjectPoints(const vector<double*> &v, vector<double*> & newv, bool& keepSmoothing, double* pointWeight, bool* keepPoint, double laplacianFactor, int weightApp) {
	double myPi = 3.14159265358979323846;
	queue<octreeNode*> finalQueue;
	finalQueue.push(root);
	while (!(finalQueue.empty())) {
		octreeNode* currNode = finalQueue.front();
		finalQueue.pop();

		if (currNode->isLeaf() && currNode->getNumOfVertex() > 0) {
			// Laplacian smoothing
			double projectedPoint[3];
			double laplacianPoint[3];
			double totalWeight = 0;
			int sv = (currNode->getVertices())[0];
			for (int d = 0; d < 3; d++) {
				projectedPoint[d] = (1 - laplacianFactor)*v[sv][d];
				laplacianPoint[d] = 0;
			}
			const vector<int> & neighborPoints = currNode->getLN();

			int neighborSize = neighborPoints.size();
			currNode->computeAvgDist(v);
			double avgDist = currNode->getAvgDist();
			for (int i = 0; i < neighborSize; i++) {
				int nsv = neighborPoints[i];
				if(!((currNode->deleteLapN)[i]) && keepPoint[nsv]) {
					double currDist = (v[sv][0] - v[nsv][0])*(v[sv][0] - v[nsv][0]) + (v[sv][1] - v[nsv][1])*(v[sv][1] - v[nsv][1]) + (v[sv][2] - v[nsv][2])*(v[sv][2] - v[nsv][2]);
					double currWeight;
					if (weightApp == 0)
						currWeight = sqrt(currDist);
					else
						currWeight = 1/(avgDist/2 + currDist);
					totalWeight = totalWeight + currWeight;
					for (int d = 0; d < 3; d++)
						laplacianPoint[d] = laplacianPoint[d] + currWeight*v[nsv][d];
				}
			}

				keepSmoothing = true;
				for (int d = 0; d < 3; d++)
					laplacianPoint[d] = laplacianFactor*laplacianPoint[d]/totalWeight;
				for (int d = 0; d < 3; d++)
					newv[sv][d] = projectedPoint[d] + laplacianPoint[d];
		}
		else if (currNode->getNumOfVertex() > 0) {
			octreeNode* const * currChild = currNode->getChildNodes();
			for (int k = 0; k < 8; k++) {
				if (currChild[k] != NULL && !(currChild[k]->isDeleted()) && currChild[k]->getNumOfVertex() > 0)
					finalQueue.push(currChild[k]);
			}
		}
	}
}

bool octree::computeNormal(const vector<double*> &tmpPoint, double* normal, double& myb)
{
	double myPi = 3.14159265358979323846;
	int nvsize = tmpPoint.size();
	double* weights = new double[nvsize];
	double* nDists = new double[nvsize];
	double avgDist = 0;
	for (int i = 0; i < nvsize; i++) {
		double currDist = 0;
		for (int d = 0; d < 3; d++)
			currDist += ((tmpPoint[i][d] - tmpPoint[0][d])*(tmpPoint[i][d] - tmpPoint[0][d]));
		nDists[i] = currDist;
		avgDist += sqrt(currDist);
	}
	avgDist = avgDist/nvsize;
	double totalWeight = 0;
	for (int i = 0; i < nvsize; i++) {
		weights[i] = exp(-1*nDists[i]/(0.34*avgDist*avgDist));
		totalWeight += weights[i];
	}

	double centroid[3];
	for (int d = 0; d < 3; d++) {
		centroid[d] = 0;
		for (int i = 0; i < nvsize; i++)
			centroid[d] = centroid[d] + weights[i]*tmpPoint[i][d];
		centroid[d] = centroid[d]/totalWeight;
	}
	double cov[3][3];
	for (int d1 = 0; d1 < 3; d1++) {
		for (int d2 = 0; d2 < 3; d2++) {
			cov[d1][d2] = 0;
			for (int i = 0; i < nvsize; i++)
				cov[d1][d2] = cov[d1][d2] + weights[i]*(tmpPoint[i][d1] - centroid[d1])*(tmpPoint[i][d2] - centroid[d2]);
		}
	}
	gsl_matrix * m = gsl_matrix_alloc(3, 3);
	for (int d1 = 0; d1 < 3; d1++) {
		for (int d2 = 0; d2 < 3; d2++)
			gsl_matrix_set (m, d1, d2, cov[d1][d2]);
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

	gsl_vector *smallestEigenVector = gsl_vector_alloc(3);
	gsl_matrix_get_col (smallestEigenVector, evec, minIndex);
	myb = 0;
	for (int d = 0; d < 3; d++)
		normal[d] = gsl_vector_get(smallestEigenVector, d);
	for (int i = 0; i < nvsize; i++) {
		for (int d = 0; d < 3; d++)
			myb = myb + normal[d]*weights[i]*tmpPoint[i][d];
	}
	myb = myb/totalWeight;
	myb = -1*myb;

	gsl_vector_free(eval);
	gsl_matrix_free(evec);
	gsl_matrix_free(m);

	return false;
}

void octree::MLSProjection(const vector<double*>& v) {
	cout << "start MLS projection" << endl;
	int vsize = v.size();
	vector<double*> projectedPoints;
	vector<double*> vNormal;
	bool* hasNormal = new bool[vsize];
	for (int i = 0; i < vsize; i++) {
		projectedPoints.push_back(new double[3]);
		vNormal.push_back(new double[3]);
		hasNormal[i] = false;
	}
	queue<octreeNode*> finalQueue;
	finalQueue.push(root);
	while (!(finalQueue.empty())) {
		octreeNode* currNode = finalQueue.front();
		finalQueue.pop();
		if (currNode->isLeaf() && currNode->getNumOfVertex() > 0) {
			octreeNode* rangeNode = currNode;
			double normal[3];
			double b;
			while (rangeNode->getSideLength() < 4*avgSideLength)
				rangeNode = rangeNode->getParentNode();
			int sv = currNode->getNearestToAvgPoint();
			int cc = 0;
			rangeNode = rangeNode->getParentNode();
			cc = 0;
			vector<double*> tmpPoint;
			Vector3d cp(v[sv]);
			octreeNode* const * cnNeighbors = rangeNode->getSameLevelNeighbors();
			queue<octreeNode*> neighborQueue;
			neighborQueue.push(rangeNode);
			for (int j = 0; j < 26; j++) {
				if (cnNeighbors[j] != NULL && cnNeighbors[j]->getNumOfVertex() > 0) {
					neighborQueue.push(cnNeighbors[j]);
				}
			}
			while (!(neighborQueue.empty())) {
				octreeNode* currNeighbor = neighborQueue.front();
				neighborQueue.pop();
				if (currNeighbor->isLeaf() && currNeighbor->getNumOfVertex() > 0) {
					int nsv = currNeighbor->getNearestToAvgPoint();
					Vector3d np(v[nsv]);
					if ((cp - np).L2Norm() < 6*avgSideLength) {
						tmpPoint.push_back(v[nsv]);
						currNode->insertSvNeighbor(nsv);
						cc++;
					}
				}
				else if (currNeighbor->getNumOfVertex() > 0) {
					octreeNode* const * currChild = currNeighbor->getChildNodes();
					for (int k = 0; k < 8; k++) {
						if (currChild[k] != NULL && !(currChild[k]->isDeleted()) && currChild[k]->getNumOfVertex() > 0)
							neighborQueue.push(currChild[k]);
					}
				}
			}
			if (cc < 7 || rangeNode == NULL) {
				currNode->setDelete(true);
				currNode->setHasNormal(false);
			}

			if (!(currNode->isDeleted())) {
				computeNormal(tmpPoint, normal, b);
				double currLength = 0;
				for (int d = 0; d < 3; d++)
					currLength += (normal[d]*normal[d]);
				for (int d = 0; d < 3; d++) {
					normal[d] = normal[d]/sqrt(currLength);
				}
				currNode->setSvNormal(normal);
				for (int d = 0; d < 3; d++)
					vNormal[sv][d] = normal[d];
				hasNormal[sv] = true;
			}
		}
		else if (currNode->getNumOfVertex() > 0) {
			octreeNode* const * currChild = currNode->getChildNodes();
			for (int k = 0; k < 8; k++)
			{
				if (currChild[k] != NULL && !(currChild[k]->isDeleted()) && currChild[k]->getNumOfVertex() > 0)
					finalQueue.push(currChild[k]);
			}
		}
	}
		cout << "finish compute normal" << endl;
		queue<octreeNode*> finalQueueq;
		finalQueueq.push(root);
		while (!(finalQueueq.empty())) {
			octreeNode* currNode = finalQueueq.front();
			finalQueueq.pop();

			if (currNode->isLeaf() && currNode->getNumOfVertex() > 0 && !(currNode->isDeleted())) {

				const vector<int>& svNeighbor = currNode->getSvNeighbor();
				int svSize = svNeighbor.size();
				int sv = currNode->getNearestToAvgPoint();
				double oldfv = 10000;
				int sign = 1;
				currNode->computeAvgDist(v);
				double avgDist = currNode->getAvgDist();
				double fv = 0;
				do {
					double totalWeight = 0;
					fv = 0;
					Vector3d currPoint(v[sv]);
					for (int j = 0; j < svSize; j++) {
						if (hasNormal[svNeighbor[j]]) {
							Vector3d neighborPoint(v[svNeighbor[j]]);
							Vector3d neighborNormal(vNormal[svNeighbor[j]]);
							double currDist = (currPoint - neighborPoint).SquareLength();
							double currWeight = exp(-1*(currDist)/(0.3*avgDist*avgDist));
							fv += (((currPoint - neighborPoint).Dot(neighborNormal))*currWeight);
							totalWeight += currWeight;
						}
					}
					fv = fv/totalWeight;
					if (fabs(fv) > fabs(oldfv))
						sign = -1;
					for (int d = 0; d < 3; d++)
						v[sv][d] = v[sv][d] + sign*fv*vNormal[sv][d];
				} while (fabs(fv) > 0.0001);
			}
			else if (!(currNode->isDeleted()) && currNode->getNumOfVertex() > 0) {
				octreeNode* const * currChild = currNode->getChildNodes();
				for (int k = 0; k < 8; k++)
				{
					if (currChild[k] != NULL && !(currChild[k]->isDeleted()) && currChild[k]->getNumOfVertex() > 0)
						finalQueueq.push(currChild[k]);
				}
			}
		}
}
