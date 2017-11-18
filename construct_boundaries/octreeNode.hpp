#include <iostream>
#include <vector>
#include <queue>
#include "grid.hpp"
using namespace std;

#ifndef OCTREENODE_HPP
#define OCTREENODE_HPP

const int XP_YP_ZP = 0;
const int X0_YP_ZP = 1;
const int XN_YP_ZP = 2;
const int XP_Y0_ZP = 3;
const int X0_Y0_ZP = 4;
const int XN_Y0_ZP = 5;
const int XP_YN_ZP = 6;
const int X0_YN_ZP = 7;
const int XN_YN_ZP = 8;
const int XP_YP_Z0 = 9;
const int X0_YP_Z0 = 10;
const int XN_YP_Z0 = 11;
const int XP_Y0_Z0 = 12;
const int XN_Y0_Z0 = 13;
const int XP_YN_Z0 = 14;
const int X0_YN_Z0 = 15;
const int XN_YN_Z0 = 16;
const int XP_YP_ZN = 17;
const int X0_YP_ZN = 18;
const int XN_YP_ZN = 19;
const int XP_Y0_ZN = 20;
const int X0_Y0_ZN = 21;
const int XN_Y0_ZN = 22;
const int XP_YN_ZN = 23;
const int X0_YN_ZN = 24;
const int XN_YN_ZN = 25;

class octreeNode
{
	private:
		double center[3];
		double avgPoint[3];
		int nearestToAvgPoint;
		double sideLength;
		octreeNode* parentNode;
		octreeNode* childNodes[8];
		octreeNode* sameLevelNeighbors[26];
		vector<int> inNodeVertices;
		vector<bool> vertexPick;
		bool isLeafNode;
		int numOfVertex;
		double maxCoord[3];
		double minCoord[3];
		bool _hasPlane;
		int level;
		bool isDeletedNode;
		bool isProcessedNode;
		bool hasPickPoint;
		double svNormal[3];
		vector<int> svNeighbor;
		double avgDist;
		double normal[3];
		bool hNormal;
		vector<int> laplacianNeighbor;
		double sixthNeighborDist;
	public:
		octreeNode();
		~octreeNode();
		const double* getCenter() const;
		double getSideLength() const;
		octreeNode* const getParentNode() const;
		octreeNode* const * getChildNodes() const;
		octreeNode* const * getSameLevelNeighbors() const;
		const vector<int> & getVertices() const;
		int getNumOfVertex() const;
		int numberOfVertex() const;
		int getLevel() const;
		bool isDeleted() const;
		const double* getMaxCoord() const;
		const double* getMinCoord() const;
		bool isLeaf() const;
		void setLeaf(bool _isLeafNode);
		void setDelete(bool _isDeletedNode);
		void setCenter(double _center[3]);
		void setSideLength(double _sideLength);
		void setParentNode(octreeNode* _parentNode);
		void setChildNode(int i, octreeNode* _childNode);
		void setSameLevelNeighbor(int i, octreeNode* _neighborCell);
		void insertVertex(int i, const vector<double*> & v);
		void setLevel(int _level);
		void setRootNode();
		void setRootNodeForMesh();
		void setInternalNode(octreeNode* _parentNode);
		void setSameLevelNeighbors(int i, octreeNode* _childNodes[]);
		void balanceNode(const vector<double*> & v, queue<octreeNode*> & myqueue);
		bool checkSplitNode();
		void splitNode(const vector<double*> & v, queue<octreeNode*> & myqueue);
		bool trim(const vector<double*> &v, vector<double>& neighborDist, vector<double*> &normal, double* myb, queue<octreeNode*> & trimQueue, const double* rMaxC, const double* rMinC);
		bool isProcessed();
		void setProcessed(bool _isProcessedNode);
		void setHasPickPoint(bool _p);
		bool isPickedPoint() const;
		bool checkSplitNodeForMesh();
		void computeAvgPoint(const vector<double*> &v);
		const double* getAvgPoint() const;
		double* getAvgPoint();
		void setAvgPoint(const double* p);
		int getNearestToAvgPoint() const;
		void insertVertexForMesh(int i);
		void setSvNormal(const double* n);
		const double* getSvNormal() const;
		void insertSvNeighbor(int i);
		const vector<int>& getSvNeighbor() const;
		void computeAvgDist(const vector<double*> & v);
		double getAvgDist();
		void setHasNormal(bool b);
		bool hasNormal();
		void splitNodeForLaplacian(const vector<double*> & v, queue<octreeNode*> & myqueue, bool& keepSplit);
		void insertLN(int i);
		const vector<int> & getLN() const;
		double getSixthNeighborDist();
		void setSixthNeighborDist(double d);
		vector<bool> deleteLapN;
};
#endif
