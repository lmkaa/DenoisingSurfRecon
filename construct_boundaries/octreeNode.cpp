#include "octreeNode.hpp"
#include <cmath>

octreeNode::~octreeNode() {
		inNodeVertices.clear();
		vertexPick.clear();
		svNeighbor.clear();
		laplacianNeighbor.clear();
		parentNode = NULL;
		for (int k = 0; k < 8; k++) {
			if (childNodes[k] != NULL) {
				delete childNodes[k];
				childNodes[k] = NULL;
			}
		}
		for (int k = 0; k > 26; k++) {
			sameLevelNeighbors[k] = NULL;
		}
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

octreeNode::octreeNode()
{
	for (int i = 0; i < 3; i++)
		center[i] = 0;
	sideLength = 0;
	parentNode = NULL;
	for (int i = 0; i < 8; i++)
		childNodes[i] = NULL;
	for (int i = 0; i < 26; i++)
		sameLevelNeighbors[i] = NULL;
	isLeafNode = false;
	isDeletedNode = false;
	isProcessedNode = false;
	hasPickPoint = false;
	_hasPlane = false;
	numOfVertex = 0;
	for (int d = 0; d < 3; d++)
	{
		maxCoord[d] = 0;
		minCoord[d] = 0;
	}
}

const double* octreeNode::getCenter() const
{
	return center;
}

double octreeNode::getSideLength() const
{
	return sideLength;
}

octreeNode* const octreeNode::getParentNode() const
{
	return parentNode;
}

octreeNode* const * octreeNode::getChildNodes() const
{
	return childNodes;
}

octreeNode* const *  octreeNode::getSameLevelNeighbors() const
{
	return sameLevelNeighbors;
}

bool octreeNode::isLeaf() const
{
	return isLeafNode;
}

void octreeNode::setLeaf(bool _isLeafNode)
{
	isLeafNode = _isLeafNode;
}
void octreeNode::setCenter(double _center[3])
{
	for (int i = 0; i < 3; i++)
		center[i] = _center[i];
}
void octreeNode::setSideLength(double _sideLength)
{
	sideLength = _sideLength;
}
void octreeNode::setParentNode(octreeNode* _parentNode)
{
	parentNode = _parentNode;
}
void octreeNode::setChildNode(int i, octreeNode* _childNode)
{
	childNodes[i] = _childNode;
}
void octreeNode::setSameLevelNeighbor(int i, octreeNode* _neighborCell)
{
	sameLevelNeighbors[i] = _neighborCell;
}

const vector<int> & octreeNode::getVertices() const
{
	return inNodeVertices;
}

void octreeNode::insertVertex(int i, const vector<double*> & v)
{
	inNodeVertices.push_back(i);
	vertexPick.push_back(false);
	numOfVertex++;
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
}

int octreeNode::getLevel() const
{
	return level;
}
void octreeNode::setLevel(int _level)
{
	level = _level;
}

const double* octreeNode::getMaxCoord() const
{
	return maxCoord;
}
const double* octreeNode::getMinCoord() const
{
	return minCoord;
}

void octreeNode::setRootNode()
{
	double maxDist = -1;
	for (int d = 0; d < 3; d++)
	{
		center[d] = (maxCoord[d] + minCoord[d])/2;
		if (d == 0 || maxDist < fabs(maxCoord[d] - minCoord[d]))
			maxDist = fabs(maxCoord[d] - minCoord[d]);
	}
	sideLength = maxDist;
	isLeafNode = true;
	parentNode = NULL;
	level = 1;
}

void octreeNode::setInternalNode(octreeNode* _parentNode)
{
	sideLength = _parentNode->getSideLength()/2;
	isLeafNode = true;
	parentNode = _parentNode;
	level = _parentNode->getLevel() + 1;
}

void octreeNode::setSameLevelNeighbors(int i, octreeNode* _childNodes[])
{
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

int octreeNode::getNumOfVertex() const
{
	return numOfVertex;
}

void octreeNode::balanceNode(const vector<double*> & v, queue<octreeNode*> & myqueue)
{
	if (level > 2)
	{
		octreeNode* const * pNeighbor = parentNode->getSameLevelNeighbors();
		for (int i = 0; i < 26; i++)
		{
			if (pNeighbor[i] != NULL && pNeighbor[i]->isLeaf() && pNeighbor[i]->getNumOfVertex() > 0)
			{
				const double* pnCenter = pNeighbor[i]->getCenter();
				if (fabs(pnCenter[0] - center[0]) < ((3*sideLength)/2 + sideLength/2000) && fabs(pnCenter[1] - center[1]) < ((3*sideLength)/2 + sideLength/2000) && fabs(pnCenter[2] - center[2]) < ((3*sideLength)/2 + sideLength/2000))
				{
					pNeighbor[i]->splitNode(v, myqueue);
				}
			}
		}
	}
}

bool octreeNode::checkSplitNode()
{
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

void octreeNode::splitNode(const vector<double*> & v, queue<octreeNode*> & myqueue)
{
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
				}
			}
			balanceNode(v, myqueue);
}

bool octreeNode::trim(const vector<double*> &v, vector<double>& neighborDist, vector<double*> &normal, double* myb, queue<octreeNode*> & trimQueue, const double* rMaxC, const double* rMinC)
{
	double myPi = 3.14159265358979323846;
	grid mygrid(this, center, sideLength);

	bool normalFound = mygrid.computeNormal(v, normal, myb);
	bool hasEmptySCell = false;

	if ((!normalFound) || hasEmptySCell)
	{
		octreeNode* const* pcNode = parentNode->getChildNodes();
		for (int i = 0; i < 8; i++)
		{
			if (pcNode[i] != NULL)
			{
				pcNode[i]->setLeaf(false);
				pcNode[i]->setDelete(true);
			}
		}
		isLeafNode = false;
		isDeletedNode = true;
		parentNode->setLeaf(true);
		trimQueue.push(parentNode);
		return false;
	}
	else
		return true;
}

void octreeNode::setDelete(bool _isDeletedNode)
{
	isDeletedNode = _isDeletedNode;
}

bool octreeNode::isDeleted() const
{
	return isDeletedNode;
}

bool octreeNode::isProcessed()
{
	return isProcessedNode;
}

void octreeNode::setProcessed(bool _isProcessedNode)
{
	isProcessedNode = _isProcessedNode;
}

void octreeNode::setHasPickPoint(bool _p)
{
	hasPickPoint = _p;
}

bool octreeNode::isPickedPoint() const
{
	return hasPickPoint;
}

bool octreeNode::checkSplitNodeForMesh()
{
	if (inNodeVertices.size() > 30)
		return true;
	else
		return false;

}

void octreeNode::setRootNodeForMesh()
{
	double maxDist = -1;
	for (int d = 0; d < 3; d++)
	{
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

void octreeNode::setSvNormal(const double* n) {
	for (int d = 0; d < 3; d++)
		svNormal[d] = n[d];
	hNormal = true;
}
const double* octreeNode::getSvNormal() const {
	return svNormal;
}
void octreeNode::insertSvNeighbor(int i) {
	svNeighbor.push_back(i);
}
const vector<int>& octreeNode::getSvNeighbor() const {
	return svNeighbor;
}
double octreeNode::getAvgDist() {
	return avgDist;
}
void octreeNode::computeAvgDist(const vector<double*> & v) {
	avgDist = 0;
	int svSize = laplacianNeighbor.size();
	for (int i = 0; i < svSize; i++) {
		double currDist = 0;
		for (int d = 0; d < 3; d++)
			currDist += ((v[laplacianNeighbor[i]][d] - v[inNodeVertices[0]][d])*(v[laplacianNeighbor[i]][d] - v[inNodeVertices[0]][d]));
		avgDist += sqrt(currDist);
	}
	avgDist = avgDist/svSize;
}
void octreeNode::setHasNormal(bool b) {
	hNormal = b;
}
bool octreeNode::hasNormal() {
	return hNormal;
}

void octreeNode::splitNodeForLaplacian(const vector<double*> & v, queue<octreeNode*> & myqueue, bool& keepSplit)
{
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
			myqueue.push(childNodes[i]);
			if (childNodes[i]->getNumOfVertex() > 1)
				keepSplit = true;
		}
	}
}

void octreeNode::insertLN(int i) {
	laplacianNeighbor.push_back(i);
}
const vector<int> & octreeNode::getLN() const {
	return laplacianNeighbor;
}

double octreeNode::getSixthNeighborDist() {
	return sixthNeighborDist;
}
void octreeNode::setSixthNeighborDist(double d) {
	sixthNeighborDist = d;
}
