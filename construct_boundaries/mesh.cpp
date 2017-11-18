#include "mesh.h"
#include <fstream>
#include <queue>

static inline double PointToSegmentDistance
(const Vector3d& p, const Vector3d& a, const Vector3d& b)
{
    Vector3d B = b - a, P = p - a;
    double t = P.Dot(B) / B.Dot(B);
    if(t < 0) return (p - a).L2Norm();
    if(t > 1) return (p - b).L2Norm();
    return (p - (a + t * B)).L2Norm();
}

static inline double PointToTriangleDistance
(const Vector3d& p, const Vector3d& a, const Vector3d& b, const Vector3d& c, bool& inside)
{
    Vector3d B = b - a, C = c - a, P = p - a;
    double BP = B.Dot(P), BB = B.Dot(B), CP = C.Dot(P), CC = C.Dot(C), BC = B.Dot(C);
    double div = BB * CC - BC * BC;
    double u = (BP * CC - CP * BC) / div;
    double v = (CP * BB - BP * BC) / div;
    if(u > 0 && v > 0 && u + v < 1)
    {
		inside = true;
        return (p - (a + u * B + v * C)).L2Norm();
    }
    else
    {
		inside = false;
        double d1 = PointToSegmentDistance(p, a, b);
        double d2 = PointToSegmentDistance(p, a, c);
        double d3 = PointToSegmentDistance(p, b, c);
        return std::min(std::min(d1, d2), d3);
    }
}

static inline bool OnSameSide(const Vector3d& e1,
                              const Vector3d& e2,
                              const Vector3d& u,
                              const Vector3d& v)
{
    return (((e1 - e2).unitcross(u - e2)).Dot((e1 - e2).unitcross(v - e2)) > 0);
}

static inline double Area(const Vector3d& a, const Vector3d& b, const Vector3d& c) {
	Vector3d ab = a-b;
	Vector3d cb = c-b;
	return ( ab.L2Norm()*cb.L2Norm()*sin(ab.angle(cb))/2  );

}

Vertex::Vertex() {
	position = Vector3d(0,0,0);
	index = -1;
	visited = false;
}

Vertex::Vertex(const Vector3d & _p) {
	position = _p;
	visited = false;
}

void Vertex::setPosition(const Vector3d & _p) {
	position = _p;
}

const Vector3d & Vertex::getPosition() const {
	return position;
}
void Vertex::eraseFace(Face* f) {
	int f_list_size = f_list.size();
	for (int i = 0; i < f_list_size; i++) {
		if (f_list[i] == f) {
			f_list.erase(f_list.begin() + i);
			return ;
		}
	}
}

Face::Face() {
	isDeleted = false;
	flag = false;
	circumRadius = -1;
	for (int t = 0; t < 3; t++)
		vIndex[t] = -1;
}

Face::Face(int _v[3]) {
	isDeleted = false;
	flag = false;
	circumRadius = -1;
	for (int t = 0; t < 3; t++)
		vIndex[t] = _v[t];
}

Face::Face(int _v1, int _v2, int _v3) {
	isDeleted = false;
	flag = false;
	circumRadius = -1;
	vIndex[0] = _v1;
	vIndex[1] = _v2;
	vIndex[2] = _v3;
}

void Face::computeCircumRadius(const VertexList& vList) {
	// va - vb
	Vector3d a = vList[vIndex[0]]->getPosition() - vList[vIndex[1]]->getPosition();
	// va - vc
	Vector3d b = vList[vIndex[0]]->getPosition() - vList[vIndex[2]]->getPosition();
	// vb - vc
	Vector3d c = vList[vIndex[1]]->getPosition() - vList[vIndex[2]]->getPosition();
	double lenA = sqrt(a.Dot(a));
	double lenB = sqrt(b.Dot(b));
	double lenC = sqrt(c.Dot(c));

	circumRadius = (lenA*lenB*lenC)/sqrt((lenA + lenB + lenC)*(lenB + lenC - lenA)*(lenC + lenA - lenB)*(lenA + lenB - lenC));
	if (isnan(circumRadius)) {
		cout << lenA << endl;
		cout << lenB << endl;
		cout << lenC << endl;
		cin.get();
	}
}

double Face::getCircumRadius() const {
	return circumRadius;
}

void Face::computeFaceArea(const VertexList& v) {
	area = Area(v[vIndex[0]]->getPosition(), v[vIndex[1]]->getPosition(), v[vIndex[2]]->getPosition());
}

double Face::getArea() const {
	return area;
}

void Face::computeMaxLength(const VertexList& v) {
	Vector3d a = v[vIndex[0]]->getPosition();
	Vector3d b = v[vIndex[1]]->getPosition();
	Vector3d c = v[vIndex[2]]->getPosition();
	Vector3d ab = b - a;
	Vector3d bc = c - b;
	Vector3d ac = c - a;
	double abL = ab.L2Norm();
	double bcL = bc.L2Norm();
	double acL = ac.L2Norm();
	if (abL > bcL && abL > acL)
		maxLength = abL;
	else if (bcL > abL && bcL > acL)
		maxLength = bcL;
	else
		maxLength = acL;
}
double Face::getMaxLength() const {
	return maxLength;
}

void Face::setVerticesIndex(int _v[3]) {
	for (int t = 0; t < 3; t++)
		vIndex[t] = _v[t];
}

void Face::setVerticesIndex(int _v1, int _v2, int _v3) {
	vIndex[0] = _v1;
	vIndex[1] = _v2;
	vIndex[2] = _v3;
}

void Face::setVertexIndex(int _v, int _i) {
	vIndex[_i] = _v;
}

const int* Face::getVerticesIndex() const {
	return vIndex;
}

int Face::getVertexIndex(int _i) const {
	return vIndex[_i];
}

Mesh::Mesh() {
	numOfFace = 0;
	numOfVertex = 0;
}

void Mesh::insertVertex(Vertex* v) {
	vList.push_back(v);
}

void Mesh::insertFace(Face* f) {
	fList.push_back(f);
}

const VertexList& Mesh::getVertexList() const {
	return vList;
}

const FaceList& Mesh::getFaceList() const {
	return fList;
}
int k = 0;
Face* Mesh::AddFace(int v1, int v2, int v3) {
	vList[v1]->setIndex(v1);
	vList[v2]->setIndex(v2);
	vList[v3]->setIndex(v3);

	Face* f = new Face();

	f->setVerticesIndex(v1, v2, v3);
	vList[v1]->insertFace(f);
	vList[v2]->insertFace(f);
	vList[v3]->insertFace(f);
	vList[v1]->insert_n_vert(vList[v2]);
	vList[v1]->insert_n_vert(vList[v3]);
	vList[v2]->insert_n_vert(vList[v1]);
	vList[v2]->insert_n_vert(vList[v3]);
	vList[v3]->insert_n_vert(vList[v1]);
	vList[v3]->insert_n_vert(vList[v2]);

	Vector3d p1 = vList[v1]->getPosition();
	Vector3d p2 = vList[v2]->getPosition();
	Vector3d p3 = vList[v3]->getPosition();

	f->setFaceIndex(numOfFace);
	fList.push_back(f);
	numOfFace++;
	return f;
}

void Mesh::readOFF(const char* fileName) {
	if (fileName==NULL || strlen(fileName)==0) {
		cerr << "file not found" << endl;
		exit(1);
	}
	ifstream ifs(fileName);
	if (ifs.fail()) {
		cerr << "file not found" << endl;
		exit(1);
	}

	char buf[1024];
	int num_vert = 0;
	int num_tri = 0;
	int dummy;
	ifs >> buf >> num_vert >> num_tri >> dummy;
	for (int j = 0; j < num_vert; j++) {
		double* v = new double[3];
		ifs >> v[0] >> v[1] >> v[2];
		AddVertex(new Vertex(Vector3d(v[0],v[1],v[2])));
		vListDouble.push_back(v);
	}
	int countLine = 0;
	while (!(ifs.eof())) {
		char first_char = '0';
		ifs >> first_char;
		if (first_char == '3') {
			int v[3];
			ifs >> v[0] >> v[1] >> v[2];
			AddFace(v[0], v[1], v[2]);
			ifs.getline(buf, 1024);
		}
		else {
			ifs.getline(buf, 1024);
		}
	}

	ifs.close();

	int vSize = vList.size();
	int fSize = fList.size();
}

Face* Mesh::findNearestFace(const double* newPoint, bool& inside) {
	Vertex* v = findNearestVertex(newPoint);
	queue<Vertex*> myqueue;
	queue<Vertex*> visitedQueue;
	queue<Face*> visitedFaceQueue;
	myqueue.push(v);
	visitedQueue.push(v);
	v->setVisited(true);
	double minFaceDist = -1;
	Face* minFace = NULL;
	bool tmp = false;
	while (!myqueue.empty()) {
		Vertex* currV = myqueue.front();
		myqueue.pop();
		const FaceList& vfList = currV->getFaceList();
		int vfListSize = vfList.size();
		for (int i = 0; i < vfListSize; i++) {
			if (!(vfList[i]->deleted()) && !(vfList[i]->isFlag())) {
				vfList[i]->setFlag(true);
				visitedFaceQueue.push(vfList[i]);
				const int* faceVI = vfList[i]->getVerticesIndex();
				double pfDist = PointToTriangleDistance(Vector3d(newPoint[0], newPoint[1], newPoint[2]), vList[faceVI[0]]->getPosition(), vList[faceVI[1]]->getPosition(), vList[faceVI[2]]->getPosition(), tmp);
				if (pfDist < minFaceDist || minFaceDist < 0) {
					inside = tmp;
					minFaceDist = pfDist;
					minFace = vfList[i];
				}
			}
		}
		const int* minFaceVI = minFace->getVerticesIndex();
		for (int d = 0; d < 3; d++) {
			if (!(vList[minFaceVI[d]]->isVisited())) {
				myqueue.push(vList[minFaceVI[d]]);
				visitedQueue.push(vList[minFaceVI[d]]);
				vList[minFaceVI[d]]->setVisited(true);
			}
		}
	}

	while (!visitedQueue.empty()) {
		(visitedQueue.front())->setVisited(false);
		visitedQueue.pop();
	}
	while (!visitedFaceQueue.empty()) {
		(visitedFaceQueue.front())->setFlag(false);
		visitedFaceQueue.pop();
	}
	return minFace;
}

void Mesh::createOctree() {
	myOctree.buildInitOctreeForMesh(vListDouble);
	queue<octreeNode*> myqueue;
	myOctree.splitAndBalanceTreeForMesh(vListDouble, myqueue);
}

octreeNode* Mesh::findNode(const double* v) {
	octreeNode* currNode = myOctree.getRoot();
	while (!currNode->isLeaf()) {
		octreeNode* const * childNodes = currNode->getChildNodes();
		double currSideLength = currNode->getSideLength();
		double childSideLength = currSideLength/2;
		double childHalfSideLength = childSideLength/2;
		int i = 0;
		bool finished = false;
		while (i < 8 && !finished) {
			const double* childCenter = childNodes[i]->getCenter();
			if (fabs(v[0] - childCenter[0]) <= childHalfSideLength + 0.000000001 && fabs(v[1] - childCenter[1]) <= childHalfSideLength + 0.000000001 && fabs(v[2] - childCenter[2]) <= childHalfSideLength + 0.000000001) {
				currNode = childNodes[i];
				finished = true;
				break;
			}
			i++;
		}
	}
	return currNode;
}

Vertex* Mesh::findNearestVertex(const double* v) {
	octreeNode* currNode = findNode(v);
	Vector3d currP(v[0],v[1],v[2]);
	bool vertexFound = false;
	int nodeSize = currNode->getNumOfVertex();
	vector<int> visitedVertex;
	double minDist = -1;
	Vertex* minVertex = NULL;
	if (nodeSize > 0) {
		const vector<int>& nodeVertices = currNode->getVertices();
		for (int i = 0; i < nodeSize; i++) {
			if (!(vList[nodeVertices[i]]->isVisited()) && vList[nodeVertices[i]]->getFaceList().size() > 0) {
				double dist = (currP - vList[nodeVertices[i]]->getPosition()).L2Norm();
				if (dist < minDist || minDist < 0) {
					minDist = dist;
					minVertex = vList[nodeVertices[i]];
					vList[nodeVertices[i]]->setVisited(true);
					visitedVertex.push_back(nodeVertices[i]);
				}
			}
		}
	}
	bool noUpdate = false;
	while (!noUpdate) {
		noUpdate = true;
		octreeNode* const * neighborNodes = currNode->getSameLevelNeighbors();
		for (int i = 0; i < 26; i++) {
			if (neighborNodes[i] != NULL && neighborNodes[i]->getNumOfVertex() > 0) {
				nodeSize = neighborNodes[i]->getNumOfVertex();
				const vector<int>& nodeVertices = neighborNodes[i]->getVertices();
				for (int j = 0; j < nodeSize; j++) {
					if (!(vList[nodeVertices[j]]->isVisited()) && vList[nodeVertices[j]]->getFaceList().size() > 0) {
						double dist = (currP - vList[nodeVertices[j]]->getPosition()).L2Norm();
						if (dist < minDist || minDist < 0) {
							minDist = dist;
							minVertex = vList[nodeVertices[j]];
							vList[nodeVertices[j]]->setVisited(true);
							visitedVertex.push_back(nodeVertices[j]);
							noUpdate = false;
						}
					}
				}
			}
		}

		if (minVertex == NULL)
			noUpdate = false;
		currNode = currNode->getParentNode();
		if (currNode == NULL) {
			cerr << "cannot find neighbor up to parent node" << endl;
			exit(1);
		}
	}
	int visitedSize = visitedVertex.size();
	for (int i = 0; i < visitedSize; i++)
		vList[visitedVertex[i]]->setVisited(false);
	if (minVertex == NULL)
	{ cout << "minVertex is null" << endl; exit(1); }
	return minVertex;
}

void Mesh::splitFace(const double* v) {
	bool inside = true;
	Face* nearestFace = findNearestFace(v, inside);
	const int* fi = nearestFace->getVerticesIndex();
	Vector3d newP(Vector3d(v[0],v[1],v[2]));
	Vector3d a = vList[fi[0]]->getPosition();
	Vector3d b = vList[fi[1]]->getPosition();
	Vector3d c = vList[fi[2]]->getPosition();

	double dAB = PointToSegmentDistance(newP, a, b);
	double dAC = PointToSegmentDistance(newP, a, c);
	double dBC = PointToSegmentDistance(newP, b, c);

	Vector3d faceNormal = (a-c).unitcross(b-c);
	double facePAB = ((a - newP).unitcross(b - newP)).Dot(faceNormal/faceNormal.L2Norm());
	double facePAC = ((a - newP).unitcross(c - newP)).Dot(faceNormal/faceNormal.L2Norm());
	double facePBC = ((b - newP).unitcross(c - newP)).Dot(faceNormal/faceNormal.L2Norm());
	bool org_flip = true;
		int eL = -1, eR = -1, vU = -1, vD = -1;
		if (dAB <= dAC && dAB <= dBC) {
			const FaceList& afList = vList[fi[0]]->getFaceList();
			vD = findFaceVI(afList, fi[0], fi[1], fi[2]);
			eL = fi[0]; eR = fi[1]; vU = fi[2];
		}
		else if (dAC <= dAB && dAC <= dBC) {
			const FaceList& afList = vList[fi[0]]->getFaceList();
			vD = findFaceVI(afList, fi[0], fi[2], fi[1]);
			eL = fi[0]; eR = fi[2]; vU = fi[1];
		}
		else {
			const FaceList& afList = vList[fi[1]]->getFaceList();
			vD = findFaceVI(afList, fi[1], fi[2], fi[0]);
			eL = fi[1]; eR = fi[2]; vU = fi[0];
		}
		if (vD != -1) {
			org_flip = false;
			Vertex* newVertex = new Vertex(newP);
			newVertex->setIndex(numOfVertex);
			vList.push_back(newVertex);

			numOfVertex++;
			Face* f1 = AddFace(eL, newVertex->getIndex(), vU);
			Face* f2 = AddFace(eR, newVertex->getIndex(), vU);
			Face* f3 = AddFace(eL, newVertex->getIndex(), vD);
			Face* f4 = AddFace(eR, newVertex->getIndex(), vD);
			f1->setFlag(true);
			f2->setFlag(true);
			f3->setFlag(true);
			f4->setFlag(true);
			nearestFace->setDeleted(true);
			const FaceList& eLfList = vList[eL]->getFaceList();
			Face* otherFace = findFace(eLfList, eL, eR, vD);
			if (otherFace == NULL) {
				cerr << "cannot find other face" << endl;
				exit(1);
			}
			otherFace->setDeleted(true);
			flipEdge(newVertex->getIndex(), eL, eR, vU, vD, f1, f2, f3, f4);
		}
	else if (org_flip) {
		Vertex* newVertex = new Vertex(newP);
		newVertex->setIndex(numOfVertex);
		vList.push_back(newVertex);
		numOfVertex++;
		Face* f1 = AddFace(fi[0], fi[1], newVertex->getIndex());
		Face* f2 = AddFace(fi[1], fi[2], newVertex->getIndex());
		Face* f3 = AddFace(fi[2], fi[0], newVertex->getIndex());
		f1->setFlag(true);
		f2->setFlag(true);
		f3->setFlag(true);
		nearestFace->setDeleted(true);
		flipEdge(newVertex->getIndex(), fi[0], fi[1], fi[2], f1, f2, f3);
	}
}

void Mesh::flipEdge(int cv, int eL, int eR, int vU, int vD, Face* f1, Face* f2, Face* f3, Face* f4) {
	queue<int*> flipQueue;
	const FaceList& eLfList = vList[eL]->getFaceList();
	const FaceList& eRfList = vList[eR]->getFaceList();
	int newV;
	newV = findFaceVI(eLfList, eL, vU, cv);
	int* newFlip;
	if (newV != -1) {
		newFlip = new int[4];
		newFlip[0] = eL; newFlip[1] = vU; newFlip[2] = cv; newFlip[3] = newV;
		flipQueue.push(newFlip);
	}

	newV = findFaceVI(eRfList, eR, vU, cv);
	if (newV != -1) {
		newFlip = new int[4];
		newFlip[0] = eR; newFlip[1] = vU; newFlip[2] = cv; newFlip[3] = newV;
		flipQueue.push(newFlip);
	}

	newV = findFaceVI(eLfList, eL, vD, cv);
	if (newV != -1) {
		newFlip = new int[4];
		newFlip[0] = eL; newFlip[1] = vD; newFlip[2] = cv; newFlip[3] = newV;
		flipQueue.push(newFlip);
	}

	newV = findFaceVI(eRfList, eR, vD, cv);
	if (newV != -1) {
		newFlip = new int[4];
		newFlip[0] = eR; newFlip[1] = vD; newFlip[2] = cv; newFlip[3] = newV;
		flipQueue.push(newFlip);
	}

	f1->setFlag(false);
	f2->setFlag(false);
	f3->setFlag(false);
	f4->setFlag(false);
	while (!(flipQueue.empty())) {
		int* currFlip = flipQueue.front();
		flipQueue.pop();
		recursiveFlipEdge(currFlip[0], currFlip[1], currFlip[2], currFlip[3], flipQueue);
	}
}

Face* Mesh::findFace(const FaceList& fl, int v1, int v2, int v3) {
	if (v1 == v2 || v2 == v3 || v1 == v3) {
		cerr << "Invalid triangle" << endl;
		exit(1);
	}
	int flSize = fl.size();
	for (int i = 0; i < flSize; i++) {
		if (!(fl[i]->isFlag()) && !(fl[i]->deleted())) {
			const int* vI = fl[i]->getVerticesIndex();
			int count = 0;
			for (int d = 0; d < 3; d++) {
				if (vI[d] == v1 || vI[d] == v2 || vI[d] == v3)
					count++;
			}
			if (count == 3)
				return fl[i];
		}
	}
	return NULL;
}

int Mesh::findFaceVI(const FaceList& fl, int v1, int v2, int v3) {
	if (v1 == v2 || v2 == v3 || v1 == v3) {
		cerr << "Invalid triangle" << endl;
		exit(1);
	}
	int flSize = fl.size();
	for (int i = 0; i < flSize; i++) {
		if (!(fl[i]->isFlag()) && !(fl[i]->deleted())) {
			const int* vI = fl[i]->getVerticesIndex();
			int count = 0;
			int countv3 = 0;
			for (int d = 0; d < 3; d++) {
				if (vI[d] == v1 || vI[d] == v2)
					count++;
				if (vI[d] == v3)
					countv3++;
			}
			if (count == 2 && countv3 == 0) {
				for (int d = 0; d < 3; d++) {
					if (vI[d] != v1 && vI[d] != v2)
						return vI[d];
				}
			}
		}
	}
	return -1;
}

bool Mesh::flippable(int eL, int eR, int vU, int vD) {
	Vector3d eLp = vList[eL]->getPosition();
	Vector3d eRp = vList[eR]->getPosition();
	Vector3d vUp = vList[vU]->getPosition();
	Vector3d vDp = vList[vD]->getPosition();

	if (OnSameSide(vUp, vDp, eLp, eRp))
		return false;

	const FaceList& eLfList = vList[eL]->getFaceList();
	const FaceList& eRfList = vList[eR]->getFaceList();
	for (int i = 0; i < eLfList.size(); i++) {
		if (!(eLfList[i]->deleted())) {
			const int* faceVI = eLfList[i]->getVerticesIndex();
			if ((faceVI[0] == vU && faceVI[1] == vD) ||
				(faceVI[0] == vD && faceVI[1] == vU) ||
				(faceVI[0] == vU && faceVI[2] == vD) ||
				(faceVI[0] == vD && faceVI[2] == vU) ||
				(faceVI[1] == vU && faceVI[2] == vD) ||
				(faceVI[1] == vD && faceVI[2] == vU))
				return false;
		}
	}
	for (int i = 0; i < eRfList.size(); i++) {
		if (!(eRfList[i]->deleted())) {
			const int* faceVI = eRfList[i]->getVerticesIndex();
			if ((faceVI[0] == vU && faceVI[1] == vD) ||
				(faceVI[0] == vD && faceVI[1] == vU) ||
				(faceVI[0] == vU && faceVI[2] == vD) ||
				(faceVI[0] == vD && faceVI[2] == vU) ||
				(faceVI[1] == vU && faceVI[2] == vD) ||
				(faceVI[1] == vD && faceVI[2] == vU))
				return false;
		}
	}

	double face1Angles[3];
	double face2Angles[3];
	face1Angles[0] = (eLp-vDp).angle(eRp-vDp);
	face1Angles[1] = (eLp-eRp).angle(eLp-vDp);
	face1Angles[2] = (eRp-eLp).angle(eRp-vDp);
	face2Angles[0] = (eLp-vUp).angle(eRp-vUp);
	face2Angles[1] = (eLp-eRp).angle(eLp-vUp);
	face2Angles[2] = (eRp-eLp).angle(eRp-vUp);

	double minAngle = -1;
	for (int i = 0; i < 3; i++) {
		if (face1Angles[i] < minAngle || minAngle < 0) {
			minAngle = face1Angles[i];
		}
	}
	for (int i = 0; i < 3; i++) {
		if (face2Angles[i] < minAngle) {
			minAngle = face2Angles[i];
		}
	}

	if (isnan(minAngle) || minAngle <= 0 || minAngle > 1.57) {
		return true;
	}

	Vector3d& new_eLp = vUp;
	Vector3d& new_eRp = vDp;
	Vector3d& new_vUp = eRp;
	Vector3d& new_vDp = eLp;

	double new_face1Angles[3];
	double new_face2Angles[3];
	new_face1Angles[0] = (new_eLp-new_vDp).angle(new_eRp-new_vDp);
	new_face1Angles[1] = (new_eLp-new_eRp).angle(new_eLp-new_vDp);
	new_face1Angles[2] = (new_eRp-new_eLp).angle(new_eRp-new_vDp);
	new_face2Angles[0] = (new_eLp-new_vUp).angle(new_eRp-new_vUp);
	new_face2Angles[1] = (new_eLp-new_eRp).angle(new_eLp-new_vUp);
	new_face2Angles[2] = (new_eRp-new_eLp).angle(new_eRp-new_vUp);

	double new_minAngle = -1;
	for (int i = 0; i < 3; i++) {
		if (new_face1Angles[i] < new_minAngle || new_minAngle < 0) {
			new_minAngle = new_face1Angles[i];
		}
	}
	for (int i = 0; i < 3; i++) {
		if (new_face2Angles[i] < new_minAngle) {
			new_minAngle = new_face2Angles[i];
		}
	}

	if (new_minAngle > minAngle + 0.0000001)
		return true;
	if (isnan(new_minAngle) && isnan(minAngle)) {
		cerr << "org is nan and become nan after flip" << endl;
		cin.get();
	}
	return false;
}

void Mesh::recursiveFlipEdge(int eL, int eR, int vU, int vD, queue<int*> &flipQueue) {
	const FaceList& eLfList = vList[eL]->getFaceList();
	Face* f1 = findFace(eLfList, eL, eR, vU);
	if (f1 == NULL)
		return ;
	f1->setFlag(true);
	const FaceList& eRfList = vList[eR]->getFaceList();
	Face* f2 = findFace(eRfList, eL, eR, vD);
	if (f2 == NULL) {
		f1->setFlag(false);
		return ;
	}
	f2->setFlag(true);
	if (flippable(eL, eR, vU, vD)) {
		f1->setVerticesIndex(vU, eR, vD);
		f2->setVerticesIndex(vU, vD, eL);
		vList[eR]->eraseFace(f2);
		vList[eL]->eraseFace(f1);
		vList[eL]->erase_n_vert(vList[eR]);
		vList[eR]->erase_n_vert(vList[eL]);
		vList[vU]->insertFace(f2);
		vList[vD]->insertFace(f1);
		vList[vU]->insert_n_vert(vList[vD]);
		vList[vD]->insert_n_vert(vList[vU]);
		int newV = findFaceVI(eLfList, vU, eL, vD);
		int* newFlip;
		if (newV != -1) {
			newFlip = new int[4];
			newFlip[0] = vU; newFlip[1] = eL; newFlip[2] = vD; newFlip[3] = newV;
			flipQueue.push(newFlip);
		}
		newV = findFaceVI(eRfList, vU, eR, vD);
		if (newV != -1) {
			newFlip = new int[4];
			newFlip[0] = vU; newFlip[1] = eR; newFlip[2] = vD; newFlip[3] = newV;
			flipQueue.push(newFlip);
		}
		newV = findFaceVI(eLfList, vD, eL, vU);
		if (newV != -1) {
			newFlip = new int[4];
			newFlip[0] = vD; newFlip[1] = eL; newFlip[2] = vU; newFlip[3] = newV;
			flipQueue.push(newFlip);
		}
		newV = findFaceVI(eRfList, vD, eR, vU);
		if (newV != -1) {
			newFlip = new int[4];
			newFlip[0] = vD; newFlip[1] = eR; newFlip[2] = vU; newFlip[3] = newV;
			flipQueue.push(newFlip);
		}
	}
	f1->setFlag(false);
	f2->setFlag(false);
}

void Mesh::flipEdge(int cv, int v1, int v2, int v3, Face* f1, Face* f2, Face* f3) {
	queue<int*> flipQueue;
	const FaceList& v1fList = vList[v1]->getFaceList();
	int newV = findFaceVI(v1fList, v1, v2, cv);
	int* newFlip = new int[4];
	newFlip[0] = v1; newFlip[1] = v2; newFlip[2] = cv; newFlip[3] = newV;
	flipQueue.push(newFlip);

	newV = findFaceVI(v1fList, v1, v3, cv);
	newFlip = new int[4];
	newFlip[0] = v1; newFlip[1] = v3; newFlip[2] = cv; newFlip[3] = newV;
	flipQueue.push(newFlip);

	const FaceList& v3fList = vList[v3]->getFaceList();
	newV = findFaceVI(v3fList, v2, v3, cv);
	newFlip = new int[4];
	newFlip[0] = v2; newFlip[1] = v3; newFlip[2] = cv; newFlip[3] = newV;
	flipQueue.push(newFlip);

	f1->setFlag(false);
	f2->setFlag(false);
	f3->setFlag(false);
	while (!(flipQueue.empty())) {
		int* currFlip = flipQueue.front();
		flipQueue.pop();
		recursiveFlipEdge(currFlip[0], currFlip[1], currFlip[2], currFlip[3], flipQueue);
	}
}

void Mesh::flipOneEdge(int eL, int eR, int vU, int vD, Face* f1, Face* f2) {
	f1->setVerticesIndex(vU, eR, vD);
	f2->setVerticesIndex(vU, vD, eL);
	vList[eR]->eraseFace(f2);
	vList[eL]->eraseFace(f1);
	vList[eL]->erase_n_vert(vList[eR]);
	vList[eR]->erase_n_vert(vList[eL]);
	vList[vU]->insertFace(f2);
	vList[vD]->insertFace(f1);
	vList[vU]->insert_n_vert(vList[vD]);
	vList[vD]->insert_n_vert(vList[vU]);
}

void Mesh::laplacianSmoothing() {
	Vector3d* newPos = new Vector3d[numOfVertex];
	double laplacianFactor = 0.1;
	for (int i = 0; i < numOfVertex; i++) {
		newPos[i] = vList[i]->getPosition();
		set<Vertex*> & currNeighbor = vList[i]->neighborVert_list();
		if (currNeighbor.size() > 2) {
			const Vector3d& currPos = vList[i]->getPosition();
			double totalWeight = 0;
			Vector3d laplacianPos(0,0,0);
			set<Vertex*>::iterator iter;
			for (iter = currNeighbor.begin(); iter != currNeighbor.end(); iter++) {
				if ((*iter)->getIndex() != vList[i]->getIndex()) {
					double currWeight = 1;
					totalWeight += currWeight;
					laplacianPos = laplacianPos + currWeight*(*iter)->getPosition();
				}
			}
			if (totalWeight > 0) {
				laplacianPos = laplacianPos/totalWeight;
				newPos[i] = (1 - laplacianFactor)*currPos + laplacianFactor*laplacianPos;
			}
			else {
				const FaceList& vfL = vList[i]->getFaceList();
				for (int fi = 0; fi < vfL.size(); fi++) {
					vfL[fi]->setDeleted(true);
				}
			}
		}
		else {
			const FaceList& vfL = vList[i]->getFaceList();
			for (int fi = 0; fi < vfL.size(); fi++) {
				vfL[fi]->setDeleted(true);
			}
		}
	}
	for (int i = 0; i < numOfVertex; i++) {
		vList[i]->setPosition(newPos[i]);
	}
}

void Mesh::meanSdFaceArea(double& mean, double& sd, int operation) {
	int nf = fList.size();
	mean = 0;
	sd = 0;
	if (operation == 0) {
		for (int i = 0; i < nf; i++) {
			fList[i]->computeMaxLength(vList);
		}
		for (int i = 0; i < nf; i++) {
			mean += (fList[i]->getMaxLength());
		}
		mean = mean/nf;
		for (int i = 0; i < nf; i++) {
			double tmp = fList[i]->getMaxLength() - mean;
			sd += (tmp*tmp);
		}
		sd = sqrt(sd/nf);
	}
	else if (operation == 1) {
		for (int i = 0; i < nf; i++) {
			fList[i]->computeFaceArea(vList);
		}
		for (int i = 0; i < nf; i++) {
			mean += (fList[i]->getArea());
		}
		mean = mean/nf;
		for (int i = 0; i < nf; i++) {
			double tmp = fList[i]->getArea() - mean;
			sd += (tmp*tmp);
		}
		sd = sqrt(sd/nf);
	}
	else {
		for (int i = 0; i < nf; i++) {
			fList[i]->computeCircumRadius(vList);
		}
		for (int i = 0; i < nf; i++) {
			mean += (fList[i]->getCircumRadius());
		}
		mean = mean/nf;
		for (int i = 0; i < nf; i++) {
			double tmp = fList[i]->getCircumRadius() - mean;
			sd += (tmp*tmp);
		}
		sd = sqrt(sd/nf);
	}
}

void Mesh::removeLargeFace(double mean, double sd, int operation) {
	int nf = fList.size();
	if (operation == 0) {
		for (int i = 0; i < nf; i++) {
			if (fList[i]->getMaxLength() > (mean + 4*sd)) {
				fList[i]->setDeleted(true);
			}
		}
	}
	else if (operation == 1) {
		for (int i = 0; i < nf; i++) {
			if (fList[i]->getArea() > (mean + 4*sd)) {
				fList[i]->setDeleted(true);
			}
		}
	}
	else {
		for (int i = 0; i < nf; i++) {
			if (fList[i]->getCircumRadius() > (mean + 4*sd)) {
				fList[i]->setDeleted(true);
			}
		}
	}
}


