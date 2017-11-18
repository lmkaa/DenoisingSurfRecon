#ifndef GLOBALFUNCTION_HPP
#define GLOBALFUNCTION_HPP

void convertToOrder2Poly(double* orgP, double* newP) {
	newP[0] = orgP[0]*orgP[0] + orgP[1]*orgP[1] + orgP[2]*orgP[2];
	newP[1] = orgP[0];
	newP[2] = orgP[1];
	newP[3] = orgP[2];
	newP[4] = 1;
}

void normalOfOrder2Poly(double* polyFunc, double* p, double* normal) {
	normal[0] = 2*polyFunc[0]*p[0] + polyFunc[1];
	normal[1] = 2*polyFunc[0]*p[1] + polyFunc[2];
	normal[2] = 2*polyFunc[0]*p[2] + polyFunc[3];
}

bool projectToOrder2Poly(const vector<double*>& v, double* polyFunc, double b, double* normal, double* orgP, double* projectedP1, double* projectedP2) {
	double A = polyFunc[0]*(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
	double B = 2*polyFunc[0]*(orgP[0]*normal[0] + orgP[1]*normal[1] + orgP[2]*normal[2]) + polyFunc[1]*normal[0] + polyFunc[2]*normal[1] + polyFunc[3]*normal[2];
	double C = polyFunc[0]*(orgP[0]*orgP[0] + orgP[1]*orgP[1] + orgP[2]*orgP[2]) + polyFunc[1]*orgP[0] + polyFunc[2]*orgP[1] + polyFunc[3]*orgP[2] + b;

	double discrim = B*B - 4*A*C;

	double ip1[3];
	double ip2[3];

	if (discrim >= 0) {
		double sol1 = (-B + sqrt(discrim))/(2*A);
		double sol2 = (-B - sqrt(discrim))/(2*A);

		for (int d = 0; d < 3; d++) {
			ip1[d] = orgP[d] + sol1*normal[d];
			ip2[d] = orgP[d] + sol2*normal[d];
		}
		for (int d = 0; d < 3; d++) {
			projectedP1[d] = ip1[d];
			projectedP2[d] = ip2[d];
		}
		return true;
	}
	return false;
}

bool projectToOrder1Poly(double* polyFunc, double b, double* orgP, double* projectedP) {
	Vector3d currNormal(polyFunc);
	Vector3d currPoint(orgP);
	double myt = (-1*(currNormal.Dot(currPoint) + b))/(currNormal.Dot(currNormal));
	for (int d = 0; d < 3; d++)
		projectedP[d] = orgP[d] + myt*polyFunc[d];
	Vector3d projectedV(projectedP);
	if (fabs(projectedV.Dot(currNormal) + b) > 0.0000000001) {
		cout << "a point not on plane" << endl;
		for (int d = 0; d < 3; d++)
			cout << polyFunc[d] << " ";
		cout << b << endl;
		for (int d = 0; d < 3; d++)
			cout << orgP[d] << " ";
		cout << endl;
		for (int d = 0; d < 3; d++)
			cout << projectedP[d] << " ";
		cout << endl;
		cin.get();
	}
	return true;
}

void updateOctree(const vector<double*>& v, octreeNode* leafNode, int vi, double* p, double refSideLength, bool* keepPoint) {
	octreeNode* currNode = leafNode;
	const double* currCenter = currNode->getCenter();
	double currHalfSideLength = (currNode->getSideLength())/2;
	while (fabs(p[0] - currCenter[0]) > currHalfSideLength || fabs(p[1] - currCenter[1]) > currHalfSideLength || fabs(p[2] - currCenter[2]) > currHalfSideLength) {
		if (currNode->getParentNode() == NULL) {
			currNode->removeVertex(vi, keepPoint[vi]);
			keepPoint[vi] = false;
			return ;
		}
		currNode->removeVertex(vi, keepPoint[vi]);
		currNode = currNode->getParentNode();
		currCenter = currNode->getCenter();
		currHalfSideLength = (currNode->getSideLength())/2;
	}
	if (currNode != leafNode && !(currNode->isLeaf())) {
		bool done = false;
		bool needNewNode = false;
		while (!done) {
			currCenter = currNode->getCenter();
			octreeNode* const * currChild = currNode->getChildNodes();
			int nextNodeIndex = -1;
			if (p[0] >= currCenter[0] && p[1] >= currCenter[1] && p[2] >= currCenter[2]) {
				nextNodeIndex = 0;
				if (currChild[0] != NULL)
					currChild[0]->insertVertex(vi, v);
				else
					needNewNode = true;
			}
			else if (p[0] >= currCenter[0] && p[1] >= currCenter[1] && p[2] < currCenter[2]) {
				nextNodeIndex = 1;
				if (currChild[1] != NULL)
					currChild[1]->insertVertex(vi, v);
				else
					needNewNode = true;
			}
			else if (p[0] >= currCenter[0] && p[1] < currCenter[1] && p[2] >= currCenter[2]) {
				nextNodeIndex = 2;
				if (currChild[2] != NULL)
					currChild[2]->insertVertex(vi, v);
				else
					needNewNode = true;
			}
			else if (p[0] >= currCenter[0] && p[1] < currCenter[1] && p[2] < currCenter[2]) {
				nextNodeIndex = 3;
				if (currChild[3] != NULL)
					currChild[3]->insertVertex(vi, v);
				else
					needNewNode = true;
			}
			else if (p[0] < currCenter[0] && p[1] >= currCenter[1] && p[2] >= currCenter[2]) {
				nextNodeIndex = 4;
				if (currChild[4] != NULL)
					currChild[4]->insertVertex(vi, v);
				else
					needNewNode = true;
			}
			else if (p[0] < currCenter[0] && p[1] >= currCenter[1] && p[2] < currCenter[2]) {
				nextNodeIndex = 5;
				if (currChild[5] != NULL)
					currChild[5]->insertVertex(vi, v);
				else
					needNewNode = true;
			}
			else if (p[0] < currCenter[0] && p[1] < currCenter[1] && p[2] >= currCenter[2]) {
				nextNodeIndex = 6;
				if (currChild[6] != NULL)
					currChild[6]->insertVertex(vi, v);
				else
					needNewNode = true;
			}
			else if (p[0] < currCenter[0] && p[1] < currCenter[1] && p[2] < currCenter[2]) {
				nextNodeIndex = 7;
				if (currChild[7] != NULL)
					currChild[7]->insertVertex(vi, v);
				else
					needNewNode = true;
			}
			if (needNewNode) {
				currNode->createSingleNode(nextNodeIndex, v, vi);
				currNode = currChild[nextNodeIndex];
				if (currNode->getSideLength() < refSideLength)
					currNode->setLeaf(true);
				else
					currNode->setLeaf(false);
			}
			if (currNode->isLeaf())
				done = true;
			else
				currNode = currChild[nextNodeIndex];
		}
	}
}

#endif
