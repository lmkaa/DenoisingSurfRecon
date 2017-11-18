#include "octree.hpp"

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

void octree::buildInitOctree(const vector<double*> &v, const bool* keepPoint)
{
	root = new octreeNode();
	int vSize = v.size();
	for (int i = 0; i < vSize; i++)
	{
		if (keepPoint[i])
			root->insertVertex(i, v);
	}
	const double* maxCoord = root->getMaxCoord();
	const double* minCoord = root->getMinCoord();
	root->setRootNode();
}

void octree::splitAndBalanceTree(const vector<double*> &v, bool* keepPoint, queue<octreeNode*> & myqueue)
{
	myqueue.push(root);
	while (!myqueue.empty())
	{
		octreeNode* currNode = myqueue.front();
		myqueue.pop();
		if (currNode->checkSplitNode())
		{
			currNode->splitNode(v, keepPoint, myqueue);
		}
	}
}

void octree::splitTreeToSameLeafLevel(const vector<double*> & v, queue<octreeNode*> & myqueue, int op, double minSideLength) {
	int currLevel = root->getLevel();
	bool keepSplit = true;
	myqueue.push(root);
	octreeNode* currNode = myqueue.front();
	myqueue.pop();
	while (keepSplit) {
		keepSplit = false;
		while (currNode->getLevel() == currLevel) {
			currLevel = currNode->getLevel();
			currNode->splitNodeToSameLeafLevel(v, myqueue, keepSplit);
			currNode = myqueue.front();
			myqueue.pop();
		}
		if (op == 1) {
			if (currNode->getSideLength() < 1.99*minSideLength)
				keepSplit = false;
		}
		if (keepSplit) {
			currLevel = currNode->getLevel();
		}
	}
}

void octree::getLeafNodeVector(vector<octreeNode*>& leafNodeVector) {
	queue<octreeNode*> myQueue;
	myQueue.push(root);
	while (!(myQueue.empty())) {
		octreeNode* currNode = myQueue.front();
		myQueue.pop();
		if (currNode->isLeaf()) {
			leafNodeVector.push_back(currNode);
		}
		else {
			octreeNode* const * childNodes = currNode->getChildNodes();
			for (int i = 0; i < 8; i++) {
				if (childNodes[i] != NULL && !(childNodes[i]->isDeleted()) && childNodes[i]->getNumOfVertex() > 0)
					myQueue.push(childNodes[i]);
			}
		}
	}
}

octreeNode* octree::findLeafNode(const double* v) const {
	octreeNode* currNode = root;
	while (currNode != NULL && !currNode->isLeaf()) {
		octreeNode* const * childNodes = currNode->getChildNodes();
		double currSideLength = currNode->getSideLength();
		double childSideLength = currSideLength/2;
		double childHalfSideLength = childSideLength/2;
		int i = 0;
		bool finished = false;
		while (i < 8 && !finished) {
			if (childNodes[i] != NULL && !(childNodes[i]->isDeleted())) {
				const double* childCenter = childNodes[i]->getCenter();
				if (fabs(v[0] - childCenter[0]) <= childHalfSideLength + 0.000000001 && fabs(v[1] - childCenter[1]) <= childHalfSideLength + 0.000000001 && fabs(v[2] - childCenter[2]) <= childHalfSideLength + 0.000000001) {
					currNode = childNodes[i];
					finished = true;
					break;
				}
			}
			i++;
		}
		if (!finished)
			return currNode;
	}
	return currNode;
}
