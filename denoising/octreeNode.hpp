#ifndef OCTREENODE_HPP
#define OCTREENODE_HPP

#include <iostream>
#include <vector>
#include <queue>
#include "vector3.h"
using namespace std;

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
		double order1Func[3];
		double order1b;
		int nearestToAvgPoint;
		double sideLength;
		octreeNode* parentNode;
		octreeNode* childNodes[8];
		octreeNode* sameLevelNeighbors[26];
		vector<int> inNodeVertices;
		vector<bool> vertexPick;
		bool isLeafNode;
		int numOfVertex;
		int numOfKeepPoint;
		double maxCoord[3];
		double minCoord[3];
		int level;
		bool isDeletedNode;
		bool isProcessedNode;
		double ns;
		int nodeIndex;
	public:
		octreeNode();
		~octreeNode();
		int getNodeIndex();
		void setNodeIndex(int n);
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
		void balanceNode(const vector<double*> & v, bool* keepPoint, queue<octreeNode*> & myqueue);
		bool checkSplitNode();
		void splitNode(const vector<double*> & v, bool* keepPoint, queue<octreeNode*> & myqueue);
		void splitNodeToSameLeafLevel(const vector<double*> & v, queue<octreeNode*> & myqueue, bool& keepSplit);
		bool isProcessed();
		void setProcessed(bool _isProcessedNode);
		bool checkSplitNodeForMesh();
		void computeAvgPoint(const vector<double*> &v);
		void updateAvgPointWithKeepPoint(const vector<double*> &v, bool* keepPoint);
		const double* getAvgPoint() const;
		double* getAvgPoint();
		void setAvgPoint(const double* p);
		int getNearestToAvgPoint() const;
		void insertVertexForMesh(int i);
		bool computeOrder1FunctionPCA(const vector<double*>& v, bool* keepPoint, int dim);
		double computeNormalPCA(const vector<double*>& v, const vector<octreeNode*>& neighborhoodVector, int dim);
		double computeNormalPCA(vector<const double*>& points, int dim, double* normal);
		void setOrder1Func(const double* p);
		void setOrder1b(double pb);
		double* getOrder1Func();
		double getOrder1b();
		void removeVertex(int vi, bool kp);
		void createSingleNode(int ci, const vector<double*>& v, int vi);
		int getNumOfKeepPoint();
		void setNumOfKeepPoint(int n);
		void laplacianProjection(double* projectedPoint, vector<octreeNode*> neighborNodes, double laplacianFactor);
		bool checkLocallyFlat(const vector<const double*>& points, int pointSize, const Vector3d& cp, const Vector3d& normal);
		void printNode(bool* keepPoint, const vector<double*> &v, ofstream & ofs, ofstream & eofs);
		void setNS(double d);
		double getNS();
};
#endif
