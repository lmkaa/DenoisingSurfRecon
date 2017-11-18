#include <iostream>
#include <fstream>
#include <queue>
#include <ctime>
#include <sys/time.h>
#include <cmath>
#include "mesh.h"
#include <time.h>
#include <vector>

void removeLargeCircumRadiusFace(const FaceList& fList, double thershold) {
	int nf = fList.size();
	for (int i = 0; i < nf; i++) {
		if (fList[i]->getCircumRadius() > thershold) {
			fList[i]->setDeleted(true);
		}
	}
}

int main(int argc, char* argv[])
{
	if (argc < 5) {
		cerr << "usage: meshmain inputModel sideLength factor outputFileName" << endl;
		exit(1);
	}
	timeval end, beginning;
	gettimeofday(&beginning, NULL);
	Mesh mesh;
	cout << "reading mesh file" << endl;
	if (argv[1] == NULL || strlen(argv[1])==0) return false;
	mesh.readOFF(argv[1]);
	const VertexList& vList = mesh.getVertexList();
	const FaceList& fList = mesh.getFaceList();
	double sideLength = atof(argv[2]);
	double factor = atof(argv[3]);
	for (int i = 0; i < fList.size(); i++) {
		fList[i]->computeCircumRadius(vList);
	}
	double avgCircumRadius = 0;
	for (int i = 0; i < fList.size(); i++) {
		avgCircumRadius += fList[i]->getCircumRadius();
	}
	avgCircumRadius /= fList.size();
	double sdCircumRadius = 0;
	for (int i = 0; i < fList.size(); i++) {
		sdCircumRadius += ((fList[i]->getCircumRadius() - avgCircumRadius)*(fList[i]->getCircumRadius() - avgCircumRadius));
	}
	sdCircumRadius /= fList.size();
	sdCircumRadius = sqrt(sdCircumRadius);

	if (sdCircumRadius > 0) {
		removeLargeCircumRadiusFace(fList, avgCircumRadius + factor*sdCircumRadius);
	}
	int mvsize = vList.size();
	int mfsize = 0;
	for (int i = 0; i < fList.size(); i++) {
		if (!(fList[i]->deleted()))
			mfsize++;
	}
	ofstream ofss(argv[4]);
	ofss << "OFF" << endl;
	ofss << mvsize << " " << mfsize << " " << 0 << endl;
	for (int i = 0; i < vList.size(); i++) {
		Vector3d p = vList[i]->getPosition();
		ofss << p.X() << " " << p.Y() << " " << p.Z() << endl;
	}
	for (int i = 0; i < fList.size(); i++) {
		if (!(fList[i]->deleted())) {
			const int* _vi = fList[i]->getVerticesIndex();
			ofss << "3 " << _vi[0] << " " << _vi[1] << " " << _vi[2] << endl;
		}
	}
	ofss.close();
	gettimeofday(&end, NULL);
	cout << "Remove large triangle: " << end.tv_sec - beginning.tv_sec + (end.tv_usec - beginning.tv_usec)/1000000.0 << " seconds" << endl;
	return 0;
}
