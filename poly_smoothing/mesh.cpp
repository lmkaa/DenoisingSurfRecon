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
Vector3d & Vertex::getPosition() {
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

VertexList& Mesh::getVertexList() {
	return vList;
}

FaceList& Mesh::getFaceList() {
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

void polySmoothing(VertexList& vList, int vsize) {
	vector<double*> polyP;
	vector<double*> projectedP;
	for (int i = 0; i < vsize; i++) {
		polyP.push_back(new double[5]);
		projectedP.push_back(new double[3]);
		convertToOrder2Poly(vList[i]->getPosition().ToArray(), polyP[i]);
	}
	for (int i = 0; i < vsize; i++) {
		if (vList[i]->getFaceList().size() > 0) {
			set<Vertex*> & nvList = vList[i]->neighborVert_list();
			set<Vertex*>::iterator it;
			vector<Vertex*> neighborVertex;
			queue<Vertex*> traverseQueue;
			neighborVertex.push_back(vList[i]);
			vList[i]->setVisited(true);
			for (it = nvList.begin(); it != nvList.end(); it++) {
				traverseQueue.push(*it);
				neighborVertex.push_back(*it);
				(*it)->setVisited(true);
			}
			while (!(traverseQueue.empty())) {
				Vertex* currVertex = traverseQueue.front();
				traverseQueue.pop();
				set<Vertex*>& nnv = currVertex->neighborVert_list();
				for (it = nnv.begin(); it != nnv.end(); it++) {
					if (!(*it)->isVisited()) {
						neighborVertex.push_back(*it);
						(*it)->setVisited(true);
					}
				}
			}
			int nvsize = neighborVertex.size();
			vector<const double*> neighborPoints;
			vector<double*> neighborPolyPoints;
			for (int ni = 0; ni < nvsize; ni++) {
				neighborVertex[ni]->setVisited(false);
				neighborPoints.push_back(neighborVertex[ni]->getPosition().ToArray());
				neighborPolyPoints.push_back(polyP[neighborVertex[ni]->getIndex()]);
			}
			double polyFunc[4];
			double polyb;
			computePolyFunction(neighborPoints, neighborPolyPoints, vList[i]->getPosition().ToArray(), 1, 5, polyFunc, polyb);
			double normal[3];
			normalOfOrder2Poly(polyFunc, vList[i]->getPosition().ToArray(), normal);
			if (projectToOrder2Poly(polyFunc, polyb, normal, vList[i]->getPosition().ToArray(), projectedP[i])) {
				vList[i]->setPosition(Vector3d(projectedP[i]));
			}
		}
	}
}

double Mesh::laplacianSmoothing(double avgFactor, int& numOfMove) {
	Vector3d* newPos = new Vector3d[numOfVertex];
	double laplacianFactor = 0.1;
	int numOfCount = 0;
	double totalMove = 0;
	vector<bool> needMove;
	for (int i = 0; i < numOfVertex; i++) {
		needMove.push_back(true);
	}
	for (int i = 0; i < numOfVertex; i++) {
		newPos[i] = vList[i]->getPosition();
		set<Vertex*> & currNeighbor = vList[i]->neighborVert_list();
		if (currNeighbor.size() > 2) {
			numOfCount++;
			const Vector3d& currPos = vList[i]->getPosition();
			double totalWeight = 0;
			Vector3d laplacianPos(0,0,0);
			set<Vertex*>::iterator iter;
			double avgLength = 0;
			int countNeighbor = 0;
			for (iter = currNeighbor.begin(); iter != currNeighbor.end(); iter++) {
				if ((*iter)->getIndex() != vList[i]->getIndex()) {
					double currWeight = 1;
					totalWeight += currWeight;
					laplacianPos = laplacianPos + currWeight*(*iter)->getPosition();
					avgLength += ((*iter)->getPosition() - vList[i]->getPosition()).L2Norm();
					countNeighbor++;
				}
			}
			avgLength /= countNeighbor;
			if (totalWeight > 0) {
				laplacianPos = laplacianPos/totalWeight;
				newPos[i] = (1 - laplacianFactor)*currPos + laplacianFactor*laplacianPos;
				if ((newPos[i] - vList[i]->getPosition()).L2Norm() < avgLength/avgFactor)
					needMove[i] = false;
			}
			else {
				const FaceList& vfL = vList[i]->getFaceList();
				for (int fi = 0; fi < vfL.size(); fi++) {
					vfL[fi]->setDeleted(true);
				}
			}
			totalMove += (vList[i]->getPosition() - newPos[i]).length();
		}
		else {
			const FaceList& vfL = vList[i]->getFaceList();
			for (int fi = 0; fi < vfL.size(); fi++) {
				vfL[fi]->setDeleted(true);
			}
		}
	}
	numOfMove = 0;
	for (int i = 0; i < numOfVertex; i++) {
		if (needMove[i]) {
			vList[i]->setPosition(newPos[i]);
			numOfMove++;
		}
	}
	delete[] newPos;
	return totalMove/numOfVertex;
}

int main(int argc, char* argv[]) {
	if (argc < 6) {
		cerr << "Usage: polySmooth meshFileName outputFileName sideLength maxIteration avgFactor" << endl;
		exit(1);
	}
        timeval end, beginning;
        gettimeofday(&beginning, NULL);

	Mesh mesh;
	mesh.readOFF(argv[1]);
	VertexList& vList = mesh.getVertexList();
	FaceList& fList = mesh.getFaceList();
	int vsize = vList.size();
	int fsize = fList.size();
	double avgMove = 10;
	double sideLength = atof(argv[3]);
	double prevAvgMove = 100;
	int maxIteration = atoi(argv[4]);
	double avgFactor = atof(argv[5]);
	int numOfMove = 10;
	int q = 0;
	while (numOfMove > 0 && q < maxIteration) {
		numOfMove = 0;
		prevAvgMove = avgMove;
		avgMove = mesh.laplacianSmoothing(avgFactor, numOfMove);
		q++;
	}
        gettimeofday(&end, NULL);
        cout << "Final laplacian smoothing: " << end.tv_sec - beginning.tv_sec + (end.tv_usec - beginning.tv_usec)/1000000.0 << " seconds" << endl;

	ofstream ofs(argv[2]);
	ofs << "OFF" << endl;
	ofs << vsize << " " << fsize << " " << 0 << endl;
	for (int i = 0; i < vsize; i++) {
		Vector3d p = vList[i]->getPosition();
		ofs << p.X() << " " << p.Y() << " " << p.Z() << endl;
	}
	for (int i = 0; i < fsize; i++) {
		const int* _vi = fList[i]->getVerticesIndex();
		ofs << "3 " << _vi[0] << " " << _vi[1] << " " << _vi[2] << endl;
	}
	ofs.close();
}
