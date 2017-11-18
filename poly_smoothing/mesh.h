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
#include "order2Projection.hpp"

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
	Vector3d & getPosition();
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
};

class Mesh {
private:
	VertexList vList;
	vector<double*> vListDouble;
	FaceList fList;
	int numOfFace, numOfVertex;
	vector<int> grid[42][42][42];
public:
	Mesh();
	void insertVertex(Vertex* v);
	void insertFace(Face* f);
	const VertexList& getVertexList() const;
	const FaceList& getFaceList() const;
	VertexList& getVertexList() ;
	FaceList& getFaceList() ;
	void AddVertex(Vertex* v) { vList.push_back(v); numOfVertex++; }
	Face* AddFace(int v1, int v2, int v3);
	void readOFF(const char* fileName);
	Face* findNearestFace(const double* newPoint, bool& inside);
	void insertVertexToGrid(int i, int j, int k, int vIndex) { grid[i][j][k].push_back(vIndex); };
	double laplacianSmoothing(double avgFactor, int& numOfMove);
};
#endif
