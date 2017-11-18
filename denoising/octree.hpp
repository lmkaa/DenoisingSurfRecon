#ifndef OCTREE_HPP
#define OCTREE_HPP
#include <iostream>
#include <vector>
#include <queue>
#include <cmath>
#include "vector3.h"
#include "octreeNode.hpp"

using namespace std;

class octree
{
	private:
		octreeNode* root;
	public:
		octree();
		~octree();
		octreeNode* const getRoot() const;
		void setRoot(octreeNode* _root);
		void buildInitOctree(const vector<double*> &v, const bool* keepPoint);
		void splitAndBalanceTree(const vector<double*> &v, bool* keepPoint, queue<octreeNode*> & myqueue);
		void splitTreeToSameLeafLevel(const vector<double*> & v, queue<octreeNode*> & myqueue, int op, double minSideLength);
		octreeNode* findLeafNode(const double* v) const;
		void getLeafNodeVector(vector<octreeNode*>& leafNodeVector);
};
#endif
