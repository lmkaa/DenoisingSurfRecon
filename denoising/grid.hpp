#ifndef GRID_HPP
#define GRID_HPP
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "vector3.h"
#include "octreeNode.hpp"

using namespace std;
class octreeNode;

#ifndef GRIDCELL_HPP
#define GRIDCELL_HPP

class gridCell
{
	private:
		double center[3];
		double sideLength;
		int numOfVertex;
		octreeNode* workingNode;
	public:
		gridCell();
		void setCenter(double* _center);
		const double* getCenter() const;
		void setSideLength(double _sideLength);
		void setWorkingNode(octreeNode* n);
		octreeNode* getWorkingNode();
		int getNumOfVertex() const;
		double getSideLength();
		const vector<int> & getInCellVertex() const;
		int* getEdgeNeighborCellSize();
		bool setEdgesAndCheckIntersection(double innerLength, double halfLength, const double* plane, const double b, gridCell gc[3][3][3], int cx, int cy, int cz, const double* rMaxC, const double* rMinC);
		const double* getAvgPoint() const;
};
#endif

class grid
{
	private:
		octreeNode* workingNode;
		double center[3];
		double sideLength;
		gridCell gc[3][3][3];
		double minCoord[3];
		double plane[3];
		double b;
	public:
		grid(octreeNode* _node, const double* _center, double _sideLength);
		octreeNode* const getNode() const;
		double getSideLength() const;
		bool computeNormal(const vector<double*> &v, vector<double*> &normal, double* myb);
		const double* getNormal() const;
		double getb() const;
		bool hasEmptySmallCell(const vector<double*> &v, const double* rMaxC, const double* rMinC);
};
#endif
