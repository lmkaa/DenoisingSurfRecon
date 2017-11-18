#ifndef MESH_H
#define MESH_H

#include <cstdlib>
#include <cmath>
#include <vector>
#include <set>
#include <assert.h>
#include <map>
#include <list>
#include "vector3.h"
#include "octree.hpp"

// classes

class Vertex;
class Face;
class Mesh;

// types
typedef std::vector<Vertex*> VertexList;
typedef std::vector<Face*> FaceList;

class Vertex {
private:
	Vector3d position;
	FaceList f_list;
	set<Vertex*> neighbor_vertex;
	int index;
	bool visited;
public:
	Vertex();
	Vertex(const Vector3d & _p);
	void setPosition(const Vector3d & _p);
	const Vector3d & getPosition() const;
	void insert_n_vert(Vertex* v) { neighbor_vertex.insert(v); };
	set<Vertex*> & neighborVert_list() { return neighbor_vertex; };
	void setIndex(int i) { index = i; };
	int getIndex() const { return index; };
	void insertFace(Face* f) { f_list.push_back(f); };
	void eraseFace(Face* f);
	void erase_n_vert(Vertex* v) {
		neighbor_vertex.erase(v);
	}
	const FaceList& getFaceList() const { return f_list; };
	bool isVisited() const { return visited; };
	void setVisited(bool b) { visited = b; };
};

class Face {
private:
	int vIndex[3];
	int fIndex;
	bool isDeleted;
	bool flag;
	double circumRadius;
	double area;
	double maxLength;
public:
	Face();
	Face(int _v[3]);
	Face(int _v1, int _v2, int _v3);
	void setVerticesIndex(int _v[3]);
	void setVerticesIndex(int _v1, int _v2, int _v3);
	void setVertexIndex(int _v, int _i);
	const int* getVerticesIndex() const;
	int getVertexIndex(int _i) const;
	void setFaceIndex(int fi) { fIndex = fi; };
	int getFaceIndex() const { return fIndex; };
	void setDeleted(bool b) { isDeleted = b; };
	bool deleted() const { return isDeleted; };
	void setFlag(bool b) { flag = b; };
	bool isFlag() const { return flag; };
	void computeCircumRadius(const VertexList& vList);
	double getCircumRadius() const;
	void computeFaceArea(const VertexList& v);
	double getArea() const;
	void computeMaxLength(const VertexList& v);
	double getMaxLength() const;
};

class Mesh {
private:
	VertexList vList;
	vector<double*> vListDouble;
	octree myOctree;
	FaceList fList;
	int numOfFace, numOfVertex;
	vector<int> grid[42][42][42];
public:
	Mesh();
	void insertVertex(Vertex* v);
	void insertFace(Face* f);
	const VertexList& getVertexList() const;
	const FaceList& getFaceList() const;
	void AddVertex(Vertex* v) { vList.push_back(v); numOfVertex++; }
	Face* AddFace(int v1, int v2, int v3);
	void readOFF(const char* fileName);
	Face* findNearestFace(const double* newPoint, bool& inside);
	void insertVertexToGrid(int i, int j, int k, int vIndex) { grid[i][j][k].push_back(vIndex); };
	void createOctree();
	octreeNode* findNode(const double* v);
	Vertex* findNearestVertex(const double* v);
	void splitFace(const double* v);
	bool flippable(int eL, int eR, int vU, int vD);
	void recursiveFlipEdge(int eL, int eR, int vU, int vD, queue<int*> &flipQueue);
	void flipEdge(int cv, int v1, int v2, int v3, Face* f1, Face* f2, Face* f3);
	void flipEdge(int cv, int eL, int eR, int vU, int vD, Face* f1, Face* f2, Face* f3, Face* f4);
	Face* findFace(const FaceList& fl, int v1, int v2, int v3);
	int findFaceVI(const FaceList& fl, int v1, int v2, int v3);
	void flipOneEdge(int eL, int eR, int vU, int vD, Face* f1, Face* f2);
	void laplacianSmoothing();
	void meanSdFaceArea(double& mean, double& sd, int operation);
	void removeLargeFace(double mean, double sd, int operation);
};
#endif
