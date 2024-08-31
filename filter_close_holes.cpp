#include "filter_close_holes.h"
#include <D:\Learning\CGIS\meshlab\src\external\downloads\libigl-2.4.0\include\igl\harmonic.h>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCore>
#include <cmath>
#include <vcg/complex/algorithms/hole.h>
#include <vcg/complex/algorithms/smooth.h>
#include <vcg/complex/complex.h>

using namespace std;
using namespace vcg;
using namespace Eigen;

const int    INF      = 1e9;
const double _epsilon = 1e-10;
const double PI       = 3.14159265359;

void Freeze(MeshModel* m)
{
	tri::UpdatePosition<CMeshO>::Matrix(m->cm, m->cm.Tr, true);
	tri::UpdateBounding<CMeshO>::Box(m->cm);
	m->cm.shot.ApplyRigidTransformation(m->cm.Tr);
	m->cm.Tr.SetIdentity();
}

void ApplyTransform(
	MeshDocument&    md,
	const Matrix44m& tr,
	bool             toAllFlag,
	bool             freeze,
	bool             invertFlag   = false,
	bool             composeFlage = true)
{
	if (toAllFlag) {
		MeshModel* m = nullptr;
		while ((m = md.nextVisibleMesh(m))) {
			if (invertFlag)
				m->cm.Tr = vcg::Inverse(m->cm.Tr);
			if (composeFlage)
				m->cm.Tr = tr * m->cm.Tr;
			else
				m->cm.Tr = tr;
			if (freeze)
				Freeze(m);
		}

		for (RasterModel& rm : md.rasterIterator())
			if (rm.isVisible())
				rm.shot.ApplyRigidTransformation(tr);
	}
	else {
		MeshModel* m = md.mm();
		if (invertFlag)
			m->cm.Tr = vcg::Inverse(m->cm.Tr);
		if (composeFlage)
			m->cm.Tr = tr * m->cm.Tr;
		else
			m->cm.Tr = tr;
		if (freeze)
			Freeze(md.mm());
	}
}

Matrix44m rotateHoleCenter(MeshDocument& md, Point3m holeCenter)
{
	Matrix44m trRot, trTran, trTranInv, transfM;
	// Point3m axis, tranVec;

	Point3m tranVec(0, 0, 0); // suppose holeCenter is P(x, y, z) tranVec is vector OP
	tranVec = -holeCenter;
	Point3m zAxis(0, 0, 1);
	Point3m tranAxis = holeCenter ^ zAxis; // tranAxis is cross product of OP and Oz axis

	Scalarm angleRad = Angle(zAxis, holeCenter);
	Scalarm angleDeg = angleRad * (180.0 / 3.141592653589793238463);

	trRot.SetRotateDeg(angleDeg, tranAxis);
	trTran.SetTranslate(tranVec);
	trTranInv.SetTranslate(-tranVec);
	transfM = trRot * trTran;

	ApplyTransform(md, transfM, false, true);
	return transfM;
};

void rotateInverse(MeshDocument& md, Matrix44m mt)
{
	Matrix44m imt = vcg::Inverse(mt);

	ApplyTransform(md, imt, false, true);
}

struct Hole
{
	vector<int> vertIdx;
	vector<int> faceIdx;
};

struct Vector
{
	double x, y, z;

	Vector operator+(const Vector& other) const { return {x + other.x, y + other.y, z + other.z}; }
	Vector operator*(const double& other) const { return {x * other, y * other, z * other}; }
	Vector operator/(const double& other) const { return {x / other, y / other, z / other}; }

	Vector cross(const Vector& other) const
	{
		return {y * other.z - z * other.y, z * other.x - x * other.z, x * other.y - y * other.x};
	}
	double dot(const Vector& other) const { return x * other.x + y * other.y + z * other.z; }
	double norm() const { return sqrt(x * x + y * y + z * z); }
};

struct Triangle
{
	int       a, b, c;
	Triangle& operator=(const Triangle& other)
	{
		a = other.a;
		b = other.b;
		c = other.c;
		return *this;
	};
};

QhullPlugin::QhullPlugin()
{
	typeList = {FP_CLOSEHOLES};

	for (ActionIDType tt : types())
		actionList.push_back(new QAction(filterName(tt), this));
}

QhullPlugin::~QhullPlugin()
{
}

QString QhullPlugin::pluginName() const
{
	return "FilterCloseHoles2";
}

QString QhullPlugin::filterName(ActionIDType f) const
{
	switch (f) {
	case FP_CLOSEHOLES: return QString("Close holes 2");
	default: assert(0); return QString();
	}
}

QString QhullPlugin::pythonFilterName(ActionIDType filterId) const
{
	switch (filterId) {
	case FP_CLOSEHOLES: return QString("Method to close holes");
	default: assert(0); return QString();
	}
}

QString QhullPlugin::filterInfo(ActionIDType filterId) const
{
	switch (filterId) {
	case FP_CLOSEHOLES: return QString("Method to close holes");
	default: assert(0);
	}
	return QString("Error: Unknown Filter");
}

QhullPlugin::FilterClass QhullPlugin::getClass(const QAction* a) const
{
	switch (ID(a)) {
	case FP_CLOSEHOLES: return FilterClass(FilterPlugin::Remeshing);
	default: assert(0);
	}
	return FilterClass(0);
}

RichParameterList QhullPlugin::initParameterList(const QAction* action, const MeshModel& m)
{
	RichParameterList parlst;
	switch (ID(action)) {
	case FP_CLOSEHOLES:
		parlst.addParam(RichFloat(
			"edgeLength",
			(float) -1.0,
			"Edge length ",
			"Chose the edge lengths of the new patch, at normal it will predict the edge length"));
		parlst.addParam(RichFloat(
			"densityFactor",
			(float) 1.41,
			"Density factor ",
			"The value to control to length of each edges"));
		parlst.addParam(RichInt(
			"MaxHoleSize",
			(int) 50,
			"Max size to be closed ",
			"The size is expressed as number of edges composing the hole boundary"));
		parlst.addParam(RichInt(
			"itFair",
			(int) 0,
			"The number of fairing times",
			"The number of iterations using harmonic weight function "));
		parlst.addParam(RichBool(
			"holeRefine",
			false,
			"Refine the hole",
			"Add more triangle to mark the hole more regular"));
		parlst.addParam(RichBool(
			"fairing",
			false,
			"Fairing the mesh",
			"Make the hole mor smooth and adapt the surrounding curvature"));
		parlst.addParam(RichBool(
			"openSurface",
			false,
			"Open Surface",
			"Choose this if your input is the open surface mesh"));
		break;
	default: break; // do not add any parameter for the other filters
	}
	return parlst;
}

double weightFunction(Point3m a, Point3m b, Point3m c)
{
	Vector first  = {b.X() - a.X(), b.Y() - a.Y(), b.Z() - a.Z()};
	Vector second = {c.X() - a.X(), c.Y() - a.Y(), c.Z() - a.Z()};

	Vector ans = first.cross(second);
	return 0.5 * sqrt(ans.x * ans.x + ans.y * ans.y + ans.z * ans.z);
}

void trace(vector<int> subhole, CMeshO& cm, vector<vector<int>> O, vector<Triangle>& S, int i, int k)
{
	if (i + 2 == k)
		S.push_back({subhole[i], subhole[i + 1], subhole[k]});
	else {
		int curO = O[i][k];
		if (curO != i + 1)
			trace(subhole, cm, O, S, i, curO);
		S.push_back({subhole[i], subhole[curO], subhole[k]});
		if (curO != k - 1)
			trace(subhole, cm, O, S, curO, k);
	}
}

float calcSigma(CMeshO& cm, int vertIdx)
{
	vector<CMeshO::VertexPointer>      vecvt;
	face::VFIterator<CMeshO::FaceType> vfi(&cm.vert[vertIdx]);
	face::Pos<CMeshO::FaceType>        fpos;
	fpos.Set(vfi.f, &cm.vert[vertIdx]);
	face::VVOrderedStarFF<CMeshO::FaceType>(fpos, vecvt);
	float cnt      = 0.0;
	float curSigma = 0.0;
	for (auto&& vp : vecvt) {
		cnt = cnt + 1.0;
		curSigma += Distance(cm.vert[vertIdx].P(), vp->P());
	}
	curSigma /= cnt;
	return curSigma;
}

bool isIn(Point3m a, Point3m b, Point3m c, Point3m d)
{
	Eigen::Vector3d v1;
	v1 << a.X(), a.Y(), a.Z();
	Eigen::Vector3d v2;
	v2 << b.X(), b.Y(), b.Z();
	Eigen::Vector3d v3;
	v3 << c.X(), c.Y(), c.Z();
	Eigen::Vector3d op;
	op << d.X(), d.Y(), d.Z();

	if ((v1 - v2).norm() < (v3 - op).norm())
		return false;

	Eigen::Matrix<double, 2, 3> projection; // 3D to 2D projection
	Eigen::Vector3d             v10 = (v2 - v1).normalized();
	Eigen::Vector3d             v   = v3 - v1;
	Eigen::Vector3d             n   = v.cross(v10);
	Eigen::Vector3d             v20 = v10.cross(n).normalized();

	projection.row(0) = v10;
	projection.row(1) = v20;

	Eigen::Vector2d A;
	A.fill(0);

	Eigen::Matrix<double, 2, 1> B   = projection * (v2 - v1);
	Eigen::Matrix<double, 2, 1> C   = projection * (v3 - v1);
	Eigen::Matrix<double, 2, 1> v11 = projection * (op - v1);

	Eigen::Matrix4d M;
	M(0, 0) = v11.squaredNorm();
	M(1, 0) = A.squaredNorm();
	M(2, 0) = B.col(0).squaredNorm();
	M(3, 0) = C.col(0).squaredNorm();

	M(0, 1) = v11(0);
	M(1, 1) = A(0);
	M(2, 1) = B(0);
	M(3, 1) = C(0);

	M(0, 2) = v11(1);
	M(1, 2) = A(1);
	M(2, 2) = B(1);
	M(3, 2) = C(1);

	M(0, 3) = 1;
	M(1, 3) = 1;
	M(2, 3) = 1;
	M(3, 3) = 1;

	if (M.determinant() > 0)
		return false;
	else
		return true;
}

// Make sure that we not create 2 same triangles
bool check(vector<Triangle>& triangles, int vIdx1, int vIdx2, int vOpIdx2, int vOpIdx1)
{
	unordered_set<int> tri1, tri2;
	tri1.insert(vIdx1);
	tri1.insert(vOpIdx2);
	tri1.insert(vOpIdx1);
	tri2.insert(vIdx2);
	tri2.insert(vOpIdx1);
	tri2.insert(vOpIdx2);

	for (int tri_id = 0; tri_id < triangles.size(); tri_id++) {
		Triangle curTri = triangles[tri_id];

		if (tri1.count(curTri.a) > 0 && tri1.count(curTri.b) > 0 && tri1.count(curTri.c) > 0)
			return false;
		if (tri2.count(curTri.a) > 0 && tri2.count(curTri.b) > 0 && tri2.count(curTri.c) > 0)
			return false;
	}

	return true;
}

bool relaxEdge(
	CMeshO&                    cm,
	vector<Triangle>&          triangles,
	int                        vIdx1,
	int                        vIdx2,
	vector<Point3m>&           centroid,
	unordered_map<int, float>& sigma,
	vector<float>&             centroidSigma)
{
	vector<int> swapTri;
	for (int tri_id = 0; tri_id < triangles.size(); tri_id++) {
		unordered_set<int> tri;
		tri.insert(triangles[tri_id].a);
		tri.insert(triangles[tri_id].b);
		tri.insert(triangles[tri_id].c);

		if (tri.count(vIdx1) && tri.count(vIdx2))
			swapTri.push_back(tri_id);
	}

	// To make sure that it just has two triangles that share an edge v1, v2 and adjacent
	if (swapTri.size() != 2)
		return false;

	// Identify the edge that will be replaced the v1v2 edge
	int vOpIdx1 = -1;
	for (int v_id = 0; v_id < 3; v_id++) {
		int triIdx = swapTri[0];
		int curV;
		if (!v_id)
			curV = triangles[triIdx].a;
		else if (v_id == 1)
			curV = triangles[triIdx].b;
		else
			curV = triangles[triIdx].c;

		if (curV != vIdx1 && curV != vIdx2) {
			vOpIdx1 = curV;
			break;
		}
	}

	int vOpIdx2 = -1;
	for (int v_id = 0; v_id < 3; v_id++) {
		int triIdx = swapTri[1];
		int curV;
		if (!v_id)
			curV = triangles[triIdx].a;
		else if (v_id == 1)
			curV = triangles[triIdx].b;
		else
			curV = triangles[triIdx].c;

		if (curV != vIdx1 && curV != vIdx2) {
			vOpIdx2 = curV;
			break;
		}
	}

	if (vOpIdx1 == -1 || vOpIdx2 == -1)
		return false;

	if (isIn(cm.vert[vIdx1].P(), cm.vert[vIdx2].P(), cm.vert[vOpIdx1].P(), cm.vert[vOpIdx2].P()) ||
		isIn(cm.vert[vIdx1].P(), cm.vert[vIdx2].P(), cm.vert[vOpIdx2].P(), cm.vert[vOpIdx1].P())) {
		if (check(triangles, vIdx1, vIdx2, vOpIdx2, vOpIdx1)) {
			// Update two new swap triangles
			triangles[swapTri[0]] = {vOpIdx1, vOpIdx2, vIdx1};
			triangles[swapTri[1]] = {vOpIdx1, vOpIdx2, vIdx2};

			// Calculate new centroids of two new swap triangles
			centroid[swapTri[0]] =
				(cm.vert[vOpIdx1].P() + cm.vert[vOpIdx2].P() + cm.vert[vIdx1].P()) / 3.0;
			centroid[swapTri[1]] =
				(cm.vert[vOpIdx1].P() + cm.vert[vOpIdx2].P() + cm.vert[vIdx2].P()) / 3.0;

			// Calculate new centroid sigmas of two new swap triangles
			centroidSigma[swapTri[0]] = (sigma[vOpIdx1] + sigma[vOpIdx2] + sigma[vIdx1]) / 3.0;
			centroidSigma[swapTri[1]] = (sigma[vOpIdx1] + sigma[vOpIdx2] + sigma[vIdx2]) / 3.0;
			return true;
		}
	}
	return false;
}

void hashFunc(CMeshO& cm, vector<string> edges_hash)
{
	for (int i = 0; i < cm.face.size(); i++) {
		for (int j = 0; j < 3; j++) {
			int v0 = cm.face[i].V(j)->Index();
			int v1 = cm.face[i].V((j + 1) % 3)->Index();
			if (v0 > v1)
				swap(v0, v1);
			edges_hash.push_back(to_string(v0) + "_" + to_string(v1));
		}
	}
}

void closeHoleByMinimumArea(
	vector<int> subhole,
	CMeshO&           cm,
	vector<Triangle>& triangles)
{
	subhole.push_back(subhole[0]);

	int                   sz = subhole.size();
	vector<vector<float>> f(sz, vector<float>(sz, 1000000000.0));
	vector<vector<int>>   O(sz, vector<int>(sz, 0));

	// Call f[i][j] is the min area of all triangle can form by the vertices from i to j
	for (int i = 0; i <= sz - 2; i++)
		f[i][i + 1] = 0;
	for (int i = 0; i <= sz - 3; i++)
		f[i][i + 2] = weightFunction(
			cm.vert[subhole[i]].P(),
			cm.vert[subhole[i + 1]].P(),
			cm.vert[subhole[i + 2]].P());

	for (int j = 2; j <= sz - 1; j++) {
		for (int i = 0; i <= sz - j - 1; i++) {
			int k = i + j;
			for (int m = i + 1; m < k; m++) {
				float val =
					f[i][m] + f[m][k] +
					weightFunction(
						cm.vert[subhole[i]].P(), cm.vert[subhole[m]].P(), cm.vert[subhole[k]].P());
				if (val < f[i][k]) {
					f[i][k] = val;
					O[i][k] = m;
				}
			}
		}
	}
	trace(subhole, cm, O, triangles, 0, sz - 1);
	subhole.pop_back();
}

Point3m calcHoleCenter(CMeshO& cm, vector<int> vertIdx)
{
	Point3m center(0, 0, 0);
	for (int i = 0; i < vertIdx.size() - 1; i++) {
		center += cm.vert[vertIdx[i]].N();
	}
	center /= (vertIdx.size());
	return center;
}

int closeHoleByCenter2(CMeshO& cm, vector<int> vertIdx, vector<Triangle>& triangles)
{
	vertIdx.push_back(vertIdx[0]);
	Point3m center(0, 0, 0);
	for (int i = 0; i < vertIdx.size() - 1; i++) {
		center += cm.vert[vertIdx[i]].P();
	}
	center /= (vertIdx.size() - 1);

	tri::Allocator<CMeshO>::AddVertex(cm, center);
	cm.vert.back().C() = Color4b::Green;
	int centerIdx      = cm.vert.back().Index();

	for (int i = 0; i < vertIdx.size() - 1; i++) {
		triangles.push_back({vertIdx[i], vertIdx[i + 1], centerIdx});
		/*tri::Allocator<CMeshO>::AddFace(cm, vertIdx[i], vertIdx[i + 1], centerIdx);
		cm.face.back().C() = Color4b::Magenta;*/
	}

	// Update topology, normal. Delete duplicate vertices, faces
	tri::UpdateBounding<CMeshO>::Box(cm);
	tri::UpdateNormal<CMeshO>::NormalizePerVertex(cm);
	tri::UpdateTopology<CMeshO>::FaceFace(cm);
	tri::UpdateTopology<CMeshO>::VertexFace(cm);
	tri::Clean<CMeshO>::RemoveDuplicateFace(cm);
	tri::Clean<CMeshO>::RemoveDuplicateVertex(cm);

	return centerIdx;
}

void closeHoleByCenter(CMeshO& cm, vector<int> vertIdx)
{
	vertIdx.push_back(vertIdx[0]);
	Point3m center(0, 0, 0);
	for (int i = 0; i < vertIdx.size() - 1; i++) {
		center += cm.vert[vertIdx[i]].P();
	}
	center /= (vertIdx.size() - 1);

	tri::Allocator<CMeshO>::AddVertex(cm, center);
	cm.vert.back().C() = Color4b::Green;
	int centerIdx      = cm.vert.back().Index();

	for (int i = 0; i < vertIdx.size() - 1; i++) {
		tri::Allocator<CMeshO>::AddFace(cm, vertIdx[i], vertIdx[i + 1], centerIdx);
		cm.face.back().C() = Color4b::Magenta;
	}

	// Update topology, normal. Delete duplicate vertices, faces
	tri::UpdateBounding<CMeshO>::Box(cm);
	tri::UpdateNormal<CMeshO>::NormalizePerVertex(cm);
	tri::UpdateTopology<CMeshO>::FaceFace(cm);
	tri::UpdateTopology<CMeshO>::VertexFace(cm);
	tri::Clean<CMeshO>::RemoveDuplicateFace(cm);
	tri::Clean<CMeshO>::RemoveDuplicateVertex(cm);
}

void divideCentroid(
	CMeshO& cm,
	vector<Triangle>& triangles,
	vector<int> dividedTriangles,
	vector<Point3m>&                       centroid,
	vector<float>&                         centroidSigma,
	unordered_map<int, float>&             sigma,
	vector <pair<int, int>>& swapEdges,
	int patchBit)
{
	for (int i = 0; i < dividedTriangles.size(); i++) {
		int      curDividedTriIdx = dividedTriangles[i];
		Point3m  curCentroid      = centroid[curDividedTriIdx];
		Triangle curTri           = triangles[curDividedTriIdx];

		// Add centroid to the mesh
		tri::Allocator<CMeshO>::AddVertex(cm, curCentroid);
		cm.vert.back().C() = Color4b::Green;
		int curCentroidIdx = cm.vert.back().Index();
		cm.vert.back().SetUserBit(patchBit);

		// Add edge to relax
		swapEdges.push_back({curTri.a, curTri.b});
		swapEdges.push_back({curTri.b, curTri.c});
		swapEdges.push_back({curTri.c, curTri.a});

		// Create new triangles with the new centroid
		Triangle tri1 = {
			curCentroidIdx, triangles[curDividedTriIdx].a, triangles[curDividedTriIdx].b};
		Triangle tri2 = {
			curCentroidIdx, triangles[curDividedTriIdx].b, triangles[curDividedTriIdx].c};
		Triangle tri3 = {
			curCentroidIdx, triangles[curDividedTriIdx].c, triangles[curDividedTriIdx].a};

		triangles[curDividedTriIdx] = tri1;
		triangles.push_back(tri2);
		triangles.push_back(tri3);

		sigma[curCentroidIdx] = centroidSigma[curDividedTriIdx];

		// Calculate centroids and centroid sigmas for the new triangles
		Point3m centroid1 = (cm.vert[tri1.a].P() + cm.vert[tri1.b].P() + cm.vert[tri1.c].P()) / 3.0;
		Point3m centroid2 = (cm.vert[tri2.a].P() + cm.vert[tri2.b].P() + cm.vert[tri2.c].P()) / 3.0;
		Point3m centroid3 = (cm.vert[tri3.a].P() + cm.vert[tri3.b].P() + cm.vert[tri3.c].P()) / 3.0;

		float sigma1 = (sigma[tri1.a] + sigma[tri1.b] + sigma[tri1.c]) / 3.0;
		float sigma2 = (sigma[tri2.a] + sigma[tri2.b] + sigma[tri2.c]) / 3.0;
		float sigma3 = (sigma[tri3.a] + sigma[tri3.b] + sigma[tri3.c]) / 3.0;

		centroid[curDividedTriIdx]      = centroid1;
		centroidSigma[curDividedTriIdx] = sigma1;

		centroid.push_back(centroid2);
		centroid.push_back(centroid3);

		centroidSigma.push_back(sigma2);
		centroidSigma.push_back(sigma3);
	}
}

vector<int> fillInSegment(MeshModel& m, CMeshO& cm, Point3m a, Point3m b, vector<Point3m> boundaryPoints, int bitCheck, int idxMin, int idxMax, int idxSegment) 
{
	vector<int> newPoints;
	float disAB = Distance(a, b);
	float avgLen = 0;
	for (int i = 0; i < boundaryPoints.size()-1; i++) {
		avgLen += Distance(boundaryPoints[i], boundaryPoints[i + 1]);
	}
	avgLen /= (boundaryPoints.size()-1);
	float i;
	for (i = 1; i < round(disAB / avgLen); i++) {
		float ratio = (float) i / round(disAB / avgLen);
		Point3m newPoint = a + (b - a) * ratio;
		tri::Allocator<CMeshO>::AddVertex(cm, newPoint);
		newPoints.push_back(cm.vert.back().Index());
		cm.vert.back().SetUserBit(bitCheck);
	}

	int cnt1 = 0, cnt2 = 0;
	for (face::VFIterator<CMeshO::FaceType> vfi(&(cm.vert[idxMin])); !vfi.End(); ++vfi) {
		int boundaryVCnt = 0;
		if (vfi.f->V(0)->IsUserBit(bitCheck))
			boundaryVCnt++;
		if (vfi.f->V(1)->IsUserBit(bitCheck))
			boundaryVCnt++;
		if (vfi.f->V(2)->IsUserBit(bitCheck))
			boundaryVCnt++;
		if (boundaryVCnt < 2)
			cnt1++;
	}
	if (cnt1) {
		for (face::VFIterator<CMeshO::FaceType> vfi(&(cm.vert[idxMin])); !vfi.End(); ++vfi) {
			int boundaryVCnt = 0;
			if (vfi.f->V(0)->IsUserBit(bitCheck))
				boundaryVCnt++;
			if (vfi.f->V(1)->IsUserBit(bitCheck))
				boundaryVCnt++;
			if (vfi.f->V(2)->IsUserBit(bitCheck))
				boundaryVCnt++;
			if (boundaryVCnt < 2)
				if (!vfi.f->V(0)->IsUserBit(bitCheck))
					vfi.f->V(0)->C() = Color4b::Red;
				else if (!vfi.f->V(1)->IsUserBit(bitCheck))
					vfi.f->V(1)->C() = Color4b::Red;
				else if (!vfi.f->V(2)->IsUserBit(bitCheck))
					vfi.f->V(2)->C() = Color4b::Red;
		}
	}
	else {
		for (face::VFIterator<CMeshO::FaceType> vfi(&(cm.vert[idxMin])); !vfi.End(); ++vfi) {
			if (!vfi.f->V(0)->IsUserBit(bitCheck))
				vfi.f->V(0)->C() = Color4b::Red;
			else if (!vfi.f->V(1)->IsUserBit(bitCheck))
				vfi.f->V(1)->C() = Color4b::Red;
			else if (!vfi.f->V(2)->IsUserBit(bitCheck))
				vfi.f->V(2)->C() = Color4b::Red;
			break;
		}
	}
	for (face::VFIterator<CMeshO::FaceType> vfi(&(cm.vert[idxSegment])); !vfi.End(); ++vfi) {
		int boundaryVCnt = 0;
		if (vfi.f->V(0)->IsUserBit(bitCheck))
			boundaryVCnt++;
		if (vfi.f->V(1)->IsUserBit(bitCheck))
			boundaryVCnt++;
		if (vfi.f->V(2)->IsUserBit(bitCheck))
			boundaryVCnt++;
		if (boundaryVCnt < 2)
			cnt2++;
	}
	if (cnt2 == 1) {
		for (face::VFIterator<CMeshO::FaceType> vfi(&(cm.vert[idxSegment])); !vfi.End(); ++vfi) {
			int boundaryVCnt = 0;
			if (vfi.f->V(0)->IsUserBit(bitCheck))
				boundaryVCnt++;
			if (vfi.f->V(1)->IsUserBit(bitCheck))
				boundaryVCnt++;
			if (vfi.f->V(2)->IsUserBit(bitCheck))
				boundaryVCnt++;
			if (boundaryVCnt < 2)
				if (!vfi.f->V(0)->IsUserBit(bitCheck))
					vfi.f->V(0)->C() = Color4b::Red;
				else if (!vfi.f->V(1)->IsUserBit(bitCheck))
					vfi.f->V(1)->C() = Color4b::Red;
				else if (!vfi.f->V(2)->IsUserBit(bitCheck))
					vfi.f->V(2)->C() = Color4b::Red;
		}
	}
	else {
		for (face::VFIterator<CMeshO::FaceType> vfi(&(cm.vert[idxSegment])); !vfi.End(); ++vfi) {
			if (!vfi.f->V(0)->IsUserBit(bitCheck))
				vfi.f->V(0)->C() = Color4b::Red;
			else if (!vfi.f->V(1)->IsUserBit(bitCheck))
				vfi.f->V(1)->C() = Color4b::Red;
			else if (!vfi.f->V(2)->IsUserBit(bitCheck))
				vfi.f->V(2)->C() = Color4b::Red;
			break;
		}
	}

	return newPoints;
}

void refine(MeshModel&                m,
			CMeshO&                   cm,
			float                     densityFactor,
			vector<Triangle>          triangles,
			vector<int>               subhole1,
			Hole                      curHole,
			int                       patchBit,
			unordered_map<int, float>& sigma,
			unordered_map<int, int>&   revPos,
			float                     total_len)
{
	float delta = densityFactor; // Density factor

	vector<float>   centroidSigma;
	vector<Point3m> centroid;
	for (int tri_id = 0; tri_id < triangles.size(); tri_id++) {
		Triangle curTri = triangles[tri_id];
		Point3m  curCentroid =
			(cm.vert[curTri.a].P() + cm.vert[curTri.b].P() + cm.vert[curTri.c].P()) / 3.0;
		float curCentroidSigma = (sigma[curTri.a] + sigma[curTri.b] + sigma[curTri.c]) / 3.0;
		centroidSigma.push_back(curCentroidSigma);
		centroid.push_back(curCentroid);
	}

	while (true) {
		// Identify divided triangle
		bool                   created = false;
		vector<pair<int, int>> swapEdges;
		vector<int>            dividedTriangles;
		for (int tri_id = 0; tri_id < triangles.size(); tri_id++) {
			Triangle curTri           = triangles[tri_id];
			float    curCentroidSigma = centroidSigma[tri_id];
			Point3m  curCentroid      = centroid[tri_id];

			// Check 3 vertices of current triangle
			bool flag = true;
			for (int m1 = 0; m1 < 3; m1++) {
				Point3m curV;
				float   curSigma;
				if (!m1)
					curSigma = sigma[curTri.a], curV = cm.vert[curTri.a].P();
				else if (m1 == 1)
					curSigma = sigma[curTri.b], curV = cm.vert[curTri.b].P();
				else
					curSigma = sigma[curTri.c], curV = cm.vert[curTri.c].P();

				if (!(delta * Distance(curCentroid, curV) > curCentroidSigma) ||
					!(delta * Distance(curCentroid, curV) > curSigma)) {
					flag = false;
					break;
				}
			}
			if (flag ||
				((delta * Distance(cm.vert[curTri.a].P(), cm.vert[curTri.b].P()) > total_len) &&
				 (delta * Distance(cm.vert[curTri.a].P(), cm.vert[curTri.c].P()) > total_len) &&
				 (delta * Distance(cm.vert[curTri.b].P(), cm.vert[curTri.b].P()) > total_len))) {
				created = true;
				dividedTriangles.push_back(tri_id);
			}
		}
		divideCentroid(
			cm, triangles, dividedTriangles, centroid, centroidSigma, sigma, swapEdges, patchBit);

		// No triangle can be divided => The patching mesh is completed
		if (!created)
			break;

		for (int e_id = 0; e_id < swapEdges.size(); e_id++) {
			relaxEdge(
				cm,
				triangles,
				swapEdges[e_id].first,
				swapEdges[e_id].second,
				centroid,
				sigma,
				centroidSigma);
		}

		// Relax all interior edges
		int  cnt       = 0;
		bool isSwapped = false;
		while (!isSwapped && cnt <= 100) {
			isSwapped = false;
			// relax(cm, triangles, cntSwap, cntSame, centroid, sigma, centroidSigma);
			for (int tri_id = 0; tri_id < triangles.size(); tri_id++) {
				Triangle curTri = triangles[tri_id];

				// Make sure this is not boundary edge
				int cntCheck = 0;
				for (int v_id = 0; v_id < curHole.vertIdx.size(); v_id++) {
					if (curTri.a == curHole.vertIdx[v_id])
						cntCheck++;
					if (curTri.b == curHole.vertIdx[v_id])
						cntCheck++;
				}
				if (cntCheck != 2)
					isSwapped = relaxEdge(
						cm, triangles, curTri.a, curTri.b, centroid, sigma, centroidSigma);
				else if (
					abs(revPos[cm.vert[curTri.a].Index()] - revPos[cm.vert[curTri.b].Index()]) !=
					1) {
					isSwapped = relaxEdge(
						cm, triangles, curTri.a, curTri.b, centroid, sigma, centroidSigma);
				}

				cntCheck = 0;
				for (int v_id = 0; v_id < curHole.vertIdx.size(); v_id++) {
					if (curTri.b == curHole.vertIdx[v_id])
						cntCheck++;
					if (curTri.c == curHole.vertIdx[v_id])
						cntCheck++;
				}
				if (cntCheck != 2)
					isSwapped = relaxEdge(
						cm, triangles, curTri.b, curTri.c, centroid, sigma, centroidSigma);
				else if (
					abs(revPos[cm.vert[curTri.c].Index()] - revPos[cm.vert[curTri.b].Index()]) !=
					1) {
					isSwapped = relaxEdge(
						cm, triangles, curTri.b, curTri.c, centroid, sigma, centroidSigma);
				}

				cntCheck = 0;
				for (int v_id = 0; v_id < curHole.vertIdx.size(); v_id++) {
					if (curTri.a == curHole.vertIdx[v_id])
						cntCheck++;
					if (curTri.c == curHole.vertIdx[v_id])
						cntCheck++;
				}
				if (cntCheck != 2)
					isSwapped = relaxEdge(
						cm, triangles, curTri.c, curTri.a, centroid, sigma, centroidSigma);
				else if (
					abs(revPos[cm.vert[curTri.c].Index()] - revPos[cm.vert[curTri.a].Index()]) !=
					1) {
					isSwapped = relaxEdge(
						cm, triangles, curTri.c, curTri.a, centroid, sigma, centroidSigma);
				}
			}
			cnt++;
		}
	}

	// Add triangles into mesh
	for (int tri_id = 0; tri_id < triangles.size(); tri_id++) {
		Triangle curTri = triangles[tri_id];
		tri::Allocator<CMeshO>::AddFace(cm, curTri.a, curTri.b, curTri.c);
		cm.face.back().SetV();
		cm.vert[curTri.a].SetV();
		cm.vert[curTri.b].SetV();
		cm.vert[curTri.c].SetV();
		m.cm.face.back().C() = Color4b::Magenta;
	}

	// Update topology, normal. Delete duplicate vertices, faces
	tri::UpdateBounding<CMeshO>::Box(cm);
	tri::UpdateNormal<CMeshO>::NormalizePerVertex(cm);
	tri::UpdateTopology<CMeshO>::FaceFace(cm);
	tri::UpdateTopology<CMeshO>::VertexFace(cm);
	tri::Clean<CMeshO>::RemoveDuplicateFace(cm);
	tri::Clean<CMeshO>::RemoveDuplicateVertex(cm);
}

std::map<std::string, QVariant> QhullPlugin::applyFilter(
	const QAction*           filter,
	const RichParameterList& par,
	MeshDocument&            md,
	unsigned int& /*postConditionMask*/,
	vcg::CallBackPos* /* cb*/)
{
	switch (ID(filter)) {
	case FP_CLOSEHOLES: {
		MeshModel& m  = *md.mm();
		//MeshModel& nm = *md.addNewMesh("", "test");
		CMeshO&    cm = m.cm;
		m.updateDataMask(MeshModel::MM_FACECOLOR);
		m.updateDataMask(MeshModel::MM_FACEFACETOPO);
		m.updateDataMask(MeshModel::MM_VERTCOLOR);
		m.updateDataMask(MeshModel::MM_VERTFACETOPO);
		m.cm.face.EnableFFAdjacency();
		m.cm.vert.EnableVFAdjacency();
		m.cm.face.EnableVFAdjacency();

		/* Precondition */
		if (tri::Clean<CMeshO>::CountNonManifoldEdgeFF(m.cm) > 0) {
			log("Mesh is two manifold, cannot process");
			break;
		}

		// User bit define
		int patchBit    = CVertexO::NewBitFlag();
		int boundaryBit = CVertexO::NewBitFlag();

		// Original value
		int originalNumbVertices = cm.vert.size();

		// Initial parameter
		int  MaxHoleSize = par.getInt("MaxHoleSize");
		float edgeLength    = par.getFloat("edgeLength");
		float  densityFactor = par.getFloat("densityFactor");
		int  itFair      = par.getInt("itFair");
		bool holeRefine  = par.getBool("holeRefine");
		bool fairing     = par.getBool("fairing");
		bool openSurface = par.getBool("openSurface");

		/* Clear bit */
		tri::UpdateFlags<CMeshO>::FaceClearV(cm);
		tri::UpdateFlags<CMeshO>::FaceClearS(cm);

		/* Extract the hole boundary */
		vector<Hole> holeList;
		int          maxHoleSize         = -1;
		int          exteriorBoundaryIdx = -1;
		for (CMeshO::FaceIterator fi = cm.face.begin(); fi != cm.face.end(); fi++) {
			// Ignore face has deleted flag
			if ((*fi).IsD()) {
				log("This face has deleted flag, skipped!");
				continue;
			}

			// Check 3 edge of each face
			for (int i = 0; i < 3; i++) {
				// Found the boundary face that not yet visited
				if (face::IsBorder(*fi, i) && !(*fi).IsV()) {
					(*fi).SetV();
					// Starpos to iterate through face in cw order
					face::Pos<CMeshO::FaceType> sp(&*fi, i, (*fi).V(i));
					face::Pos<CMeshO::FaceType> fp = sp;

					Hole curHole;

					sp.f->SetV();
					do {
						sp.f->SetV();
						sp.NextB();
						sp.f->SetV();

						// Mark color of boundary vertices & faces
						curHole.vertIdx.push_back(sp.v->Index());
						sp.v->C() = Color4b::Blue;
						curHole.faceIdx.push_back(sp.f->Index());
						sp.f->C() = Color4b::Green;
						sp.v->SetUserBit(boundaryBit);
					} while (sp != fp);

					if (maxHoleSize < (int) curHole.vertIdx.size()) {
						maxHoleSize = (int) curHole.vertIdx.size();
						exteriorBoundaryIdx = (int) holeList.size();
					}
					holeList.push_back(curHole);
				}
			}
		}

		if (openSurface) {
			for (int i = 0; i < holeList[exteriorBoundaryIdx].faceIdx.size(); i++) {
				cm.face[holeList[exteriorBoundaryIdx].faceIdx[i]].C() = Color4b::Red;
			}
		}

		/* Hole filling process */
		// Fill hole with ONLY boundary vertices
		for (int hole_id = 0; hole_id < holeList.size(); hole_id++) {
			Hole curHole = holeList[hole_id];
			if (curHole.vertIdx.size() > MaxHoleSize)
				continue;
			else if (curHole.vertIdx.size() <= 10) {
				closeHoleByCenter(cm, curHole.vertIdx);
				continue;
			}

			for (int i = 0; i < curHole.vertIdx.size(); i++)
				if (!cm.vert[curHole.vertIdx[i]].IsUserBit(boundaryBit))
					cm.vert[curHole.vertIdx[i]].SetUserBit(patchBit);

			vector<Triangle> triangles;

			//Color again the hole
			for (int i : curHole.vertIdx) {
				cm.vert[i].C() = Color4b(
					255, i / curHole.vertIdx.size() * 255, i / curHole.vertIdx.size() * 255, 125);
			}

			// Rotate z-axis up
			Point3m center = calcHoleCenter(cm, curHole.vertIdx);
			Matrix44m rotateM = rotateHoleCenter(md, center);

			// Z_min, Z_max, Z_avg
			float minZ = 1000000000.0;
			float maxZ = -1000000000.0;
			int   idxMax  = -1;
			int   idxMin  = -1;
			int   idxMin2; // Index in vertIdx
			float avgZ = 0;
			for (int i = 0; i < curHole.vertIdx.size(); i++) {
				avgZ += cm.vert[curHole.vertIdx[i]].P().Z();
				int nxt = (i == curHole.vertIdx.size() - 1 ? 0 : i + 1);
				if (minZ > cm.vert[curHole.vertIdx[i]].P().Z()) {
					idxMin = curHole.vertIdx[i];
					idxMin2 = i;
					minZ   = cm.vert[curHole.vertIdx[i]].P().Z();
				}
				if (maxZ < cm.vert[curHole.vertIdx[i]].P().Z()) {
					idxMax = curHole.vertIdx[i];
					maxZ   = cm.vert[curHole.vertIdx[i]].P().Z();
				}
			}
			avgZ /= curHole.vertIdx.size();

			float minDiffZ = 1e9;
			int   minDiffZIdx;
			int   minDiffZIdx2;
			for (int i = 0; i < curHole.vertIdx.size(); i++) {
				if (curHole.vertIdx[i] == idxMin || avgZ < cm.vert[curHole.vertIdx[i]].P().Z() || abs(idxMin2 - i) < curHole.vertIdx.size()/3)
					continue;

				float diffZ = abs(cm.vert[curHole.vertIdx[i]].P().Z() - cm.vert[idxMin].P().Z());
				float dis   = Distance(cm.vert[curHole.vertIdx[i]].P(), cm.vert[idxMin].P());
				log("DiffZ: %f", diffZ);
				if (minDiffZ > diffZ) {
					minDiffZ = diffZ;
					minDiffZIdx = curHole.vertIdx[i];
					minDiffZIdx2 = i;
				}
				log("Distance: %f", dis);
			}
			cm.vert[minDiffZIdx].C() = Color4b::Black;
			vector<Point3m> boundaryPoints;
			for (int i = 0; i < curHole.vertIdx.size(); i++) {
				boundaryPoints.push_back(cm.vert[curHole.vertIdx[i]].P());
			}
			vector<int> newPoints = fillInSegment(
				m,
				cm,
				cm.vert[idxMin].P(),
				cm.vert[minDiffZIdx].P(),
				boundaryPoints,
				boundaryBit,
				idxMin,
				idxMax,
				minDiffZIdx);
			vector<int> subhole1, subhole2;
			//First subhole
			subhole1.push_back(curHole.vertIdx[idxMin2]);
			int nxtIdx = idxMin2;
			while (true) {
				nxtIdx = (nxtIdx == 0 ? curHole.vertIdx.size() - 1 : nxtIdx - 1);
				if (nxtIdx == minDiffZIdx2) {
					subhole1.push_back(curHole.vertIdx[nxtIdx]);
					break;
				}
				subhole1.push_back(curHole.vertIdx[nxtIdx]);
			}
			reverse(newPoints.begin(), newPoints.end());
			for (auto i : newPoints) {
				subhole1.push_back(i);
			}
			int idxCenter = closeHoleByCenter2(cm, subhole1, triangles);
			subhole1.push_back(idxCenter);

			// Second subhole
			vector<Triangle> triangles2;
			subhole2.push_back(curHole.vertIdx[idxMin2]);
			nxtIdx = idxMin2;
			while (true) {
				nxtIdx = (nxtIdx == curHole.vertIdx.size() - 1 ? 0 : nxtIdx + 1);
				if (nxtIdx == minDiffZIdx2) {
					subhole2.push_back(curHole.vertIdx[nxtIdx]);
					break;
				}
				subhole2.push_back(curHole.vertIdx[nxtIdx]);
			}
			for (auto i : newPoints) {
				subhole2.push_back(i);
			}
			idxCenter = closeHoleByCenter2(cm, subhole2, triangles2);
			subhole2.push_back(idxCenter);
			rotateInverse(md, rotateM);

			// Mesh refinement
			if (holeRefine) {
				unordered_map<int, float> sigma;
				unordered_map<int, int>   revPos;
				float                     total_len = 0;
				for (int v_id = 0; v_id < curHole.vertIdx.size(); v_id++) {
					revPos[curHole.vertIdx[v_id]] = v_id;
					int idx                       = curHole.vertIdx[v_id];
					sigma[idx]                    = (edgeLength == -1.0 ? calcSigma(cm, idx) : edgeLength);
					total_len += calcSigma(cm, idx);
				}
				total_len /= curHole.vertIdx.size();
				for (int v_id = 0; v_id < subhole1.size(); v_id++) {
					sigma[subhole1[v_id]] = (edgeLength == -1.0 ? total_len : edgeLength);
				}
				refine(m, cm, densityFactor, triangles, subhole1, curHole, patchBit, sigma, revPos, total_len);
				sigma.clear();
				revPos.clear();
				for (int v_id = 0; v_id < curHole.vertIdx.size(); v_id++) {
					revPos[curHole.vertIdx[v_id]] = v_id;
					int idx                       = curHole.vertIdx[v_id];
					sigma[idx]                    = (edgeLength == -1.0 ? total_len : edgeLength);
				}
				for (int v_id = 0; v_id < subhole2.size(); v_id++) {
					sigma[subhole2[v_id]] = (edgeLength == -1.0 ? total_len : edgeLength);
				}
				refine(m, cm, densityFactor, triangles2, subhole2, curHole, patchBit, sigma, revPos, total_len);
				if (fairing) {
					for (int it = 0; it < itFair; it++) {
						int      vNum   = cm.vert.size();
						int      fNum   = cm.face.size();
						int      inVNum = 0;
						MatrixXd V(vNum, 3);
						MatrixXi F(fNum, 3);
						for (int i = 0; i < vNum; i++) {
							if (cm.vert[i].IsUserBit(patchBit))
								inVNum++;
							Point3m curPoint = cm.vert[i].P();
							V(i, 0)          = curPoint.X();
							V(i, 1)          = curPoint.Y();
							V(i, 2)          = curPoint.Z();
						}
						for (int i = 0; i < fNum; i++) {
							F(i, 0) = cm.face[i].V(0)->Index();
							F(i, 1) = cm.face[i].V(1)->Index();
							F(i, 2) = cm.face[i].V(2)->Index();
						}
						MatrixXd bc(vNum - inVNum, 3);
						VectorXi b(vNum - inVNum);
						int      idx = 0;
						for (int i = 0; i < vNum; i++) {
							if (!cm.vert[i].IsUserBit(patchBit)) {
								b(idx)     = i;
								bc(idx, 0) = cm.vert[i].P().X();
								bc(idx, 1) = cm.vert[i].P().Y();
								bc(idx, 2) = cm.vert[i].P().Z();
								idx++;
							}
						}
						MatrixXd update;
						igl::harmonic(V, F, b, bc, 2, update);
						for (int i = 0; i < vNum; i++) {
							if (cm.vert[i].IsUserBit(patchBit))
								if (!std::isnan(update(i, 0)) && !std::isnan(update(i, 1)) &&
									!std::isnan(update(i, 2))) {
									cm.vert[i].P() =
										Point3m(update(i, 0), update(i, 1), update(i, 2));
								}
						}
					}
					for (int i = 0; i < curHole.vertIdx.size(); i++)
						if (cm.vert[curHole.vertIdx[i]].IsUserBit(patchBit))
							cm.vert[curHole.vertIdx[i]].ClearUserBit(patchBit);
				}
			}
		}

	} break;
	default: wrongActionCalled(filter);
	}
	return std::map<std::string, QVariant>();
}
MESHLAB_PLUGIN_NAME_EXPORTER(QhullPlugin)
