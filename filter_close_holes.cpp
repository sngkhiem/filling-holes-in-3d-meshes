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
	case FP_CLOSEHOLES: return QString("Close holes 2 (Liepa)");
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

void trace(Hole hole, CMeshO& cm, vector<vector<int>> O, vector<Triangle>& S, int i, int k)
{
	if (i + 2 == k)
		S.push_back({hole.vertIdx[i], hole.vertIdx[i + 1], hole.vertIdx[k]});
	else {
		int curO = O[i][k];
		if (curO != i + 1)
			trace(hole, cm, O, S, i, curO);
		S.push_back({hole.vertIdx[i], hole.vertIdx[curO], hole.vertIdx[k]});
		if (curO != k - 1)
			trace(hole, cm, O, S, curO, k);
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
	Hole&              curHole,
	CMeshO&           cm,
	vector<Triangle>& triangles)
{
	curHole.vertIdx.push_back(curHole.vertIdx[0]);

	int                   sz = curHole.vertIdx.size();
	vector<vector<float>> f(sz, vector<float>(sz, 1000000000.0));
	vector<vector<int>>   O(sz, vector<int>(sz, 0));

	// Call f[i][j] is the min area of all triangle can form by the vertices from i to j
	for (int i = 0; i <= sz - 2; i++)
		f[i][i + 1] = 0;
	for (int i = 0; i <= sz - 3; i++)
		f[i][i + 2] = weightFunction(
			cm.vert[curHole.vertIdx[i]].P(), cm.vert[curHole.vertIdx[i + 1]].P(), cm.vert[curHole.vertIdx[i + 2]].P());

	for (int j = 2; j <= sz - 1; j++) {
		for (int i = 0; i <= sz - j - 1; i++) {
			int k = i + j;
			for (int m = i + 1; m < k; m++) {
				float val =
					f[i][m] + f[m][k] +
					weightFunction(
						cm.vert[curHole.vertIdx[i]].P(), cm.vert[curHole.vertIdx[m]].P(), cm.vert[curHole.vertIdx[k]].P());
				if (val < f[i][k]) {
					f[i][k] = val;
					O[i][k] = m;
				}
			}
		}
	}

	trace(curHole, cm, O, triangles, 0, sz - 1);
	curHole.vertIdx.pop_back();
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


void closeHoleByConcavity(CMeshO& cm, vector<int> vertIdx, vector<Triangle>& triangles)
{
	vector<string> edges_hash;
	hashFunc(cm, edges_hash);

	if (vertIdx.size() == 3) {
		triangles.push_back({vertIdx[0], vertIdx[1], vertIdx[2]});
		return;
	}

	vector<vector<int>> queue = {vector<int>(vertIdx.rbegin(), vertIdx.rend())};
	while (!queue.empty()) {
		vector<int> curV = queue.back();
		queue.pop_back();
		if (curV.size() == 3) {
			triangles.push_back({curV[0], curV[1], curV[2]});
			continue;
		}

		vector<double> edge_len(curV.size());
		for (int i = 0; i < curV.size(); i++) {
			edge_len[i] = (cm.vert[curV[(i + 1) % curV.size()]].P() - cm.vert[curV[i]].P()).Norm();
		}

		double hole_len           = accumulate(edge_len.begin(), edge_len.end(), 0.0);
		double min_concave_degree = 1000000000.0;
		int    target_i = -1, target_j = -1;
		for (int i = 0; i < curV.size(); i++) {
			vector<double> eu_dists(curV.size());
			for (int j = 0; j < curV.size(); j++) {
				eu_dists[j] = (cm.vert[curV[i]].P() - cm.vert[curV[j]].P()).Norm();
			}

			for (int j = 0; j < curV.size(); j++) {
				int    v0        = min(curV[i], curV[j]);
				int    v1        = max(curV[i], curV[j]);
				string edge_hash = to_string(v0) + "_" + to_string(v1);
				if (find(edges_hash.begin(), edges_hash.end(), edge_hash) != edges_hash.end()) {
					eu_dists[j] = 1000000000.0;
				}
			}

			vector<double> geo_dists(curV.size());
			partial_sum(edge_len.begin(), edge_len.end(), geo_dists.begin());
			rotate(geo_dists.begin(), geo_dists.begin() + i, geo_dists.end());
			transform(geo_dists.begin(), geo_dists.end(), geo_dists.begin(), [hole_len](double d) {
				return min(d, hole_len - d);
			});
			rotate(geo_dists.begin(), geo_dists.begin() + i, geo_dists.end());
			vector<double> concave_degree(curV.size());
			transform(
				eu_dists.begin(),
				eu_dists.end(),
				geo_dists.begin(),
				concave_degree.begin(),
				[](double eu, double geo) { return eu / (geo * geo + _epsilon); });
			concave_degree[i] = -numeric_limits<double>::infinity(); // There may exist two duplicate vertices

			vector<int> sorted_indices(concave_degree.size());
			iota(sorted_indices.begin(), sorted_indices.end(), 0);
			sort(sorted_indices.begin(), sorted_indices.end(), [&concave_degree](int a, int b) {
				return concave_degree[a] < concave_degree[b];
			});

			int idx = 1;
			int j   = sorted_indices[idx];
			while (min((j + curV.size() - i) % curV.size(), (i + curV.size() - j) % curV.size()) <=
				   1) {
				idx++;
				j = sorted_indices[idx];
			}

			if (concave_degree[j] < min_concave_degree) {
				min_concave_degree = concave_degree[j];
				target_i           = min(i, (int) j);
				target_j           = max(i, (int) j);
			}
		}

		queue.push_back(vector<int>(curV.begin() + target_i, curV.begin() + target_j + 1));
		queue.push_back(vector<int>(curV.begin() + target_j, curV.end()));
		queue.back().insert(queue.back().end(), curV.begin(), curV.begin() + target_i + 1);
	}
}

double calculateAngle(const Point3f& p1, const Point3f& p2, const Point3f& p3)
{
	Point3f v1 = p1 - p2;
	Point3f v2 = p3 - p2;
	double dotProduct = v1.dot(v2);
	return acos(dotProduct);
}


void closeHoleByConcavity2(CMeshO& cm, vector<int> vertIdx, vector<Triangle>& triangles)
{
	vector<string> edges_hash;
	hashFunc(cm, edges_hash);

	if (vertIdx.size() == 3) {
		triangles.push_back({vertIdx[0], vertIdx[1], vertIdx[2]});
		return;
	}

	vector<vector<int>> queue = {vector<int>(vertIdx.rbegin(), vertIdx.rend())};
	while (!queue.empty()) {
		vector<int> curV = queue.back();
		queue.pop_back();
		if (curV.size() == 3) {
			triangles.push_back({curV[0], curV[1], curV[2]});
			continue;
		}

		// Calculate angles
		double min_angle = numeric_limits<double>::max();
		int    min_idx   = -1;
		for (int i = 0; i < curV.size(); i++) {
			int    prev  = (i - 1 + curV.size()) % curV.size();
			int    next  = (i + 1) % curV.size();
			double angle = calculateAngle(
				cm.vert[curV[prev]].P(), cm.vert[curV[i]].P(), cm.vert[curV[next]].P());
			if (angle < min_angle) {
				min_angle = angle;
				min_idx   = i;
			}
		}

		// Split the hole at the minimum angle
		int prev = (min_idx - 1 + curV.size()) % curV.size();
		int next = (min_idx + 1) % curV.size();
		queue.push_back(vector<int>(curV.begin() + prev, curV.begin() + next + 1));
		queue.push_back(vector<int>(curV.begin() + next, curV.end()));
		queue.back().insert(queue.back().end(), curV.begin(), curV.begin() + prev + 1);
	}
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

void fillInSegment(MeshModel& m, CMeshO& cm, Point3m a, Point3m b, vector<Point3m> boundaryPoints) 
{
	float disAB = Distance(a, b);
	float avgLen = 0;
	for (int i = 0; i < boundaryPoints.size()-1; i++) {
		avgLen += Distance(boundaryPoints[i], boundaryPoints[i + 1]);
	}
	avgLen /= (boundaryPoints.size()-1);

	qDebug("%f, %f", disAB, avgLen);
	qDebug("%f", round(disAB / avgLen));

	for (float i = 1; i < round(disAB / avgLen); i++) {
		float ratio = (float) i / round(disAB / avgLen);
		Point3m newPoint = a + (b - a) * ratio;
		tri::Allocator<CMeshO>::AddVertex(cm, newPoint);
		//m.updateDataMask();
	}
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
			closeHoleByConcavity(cm, curHole.vertIdx, triangles);

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
			float avgZ = 0;
			vector<Point3m> boundaryPoints;
			for (int i = 0; i < curHole.vertIdx.size(); i++) {
				boundaryPoints.push_back(cm.vert[curHole.vertIdx[i]].P());
				avgZ += cm.vert[curHole.vertIdx[i]].P().Z();

				if (maxZ < cm.vert[curHole.vertIdx[i]].P().Z()) {
					maxZ = cm.vert[curHole.vertIdx[i]].P().Z();
					idxMax  = curHole.vertIdx[i];
				}
				if (minZ > cm.vert[curHole.vertIdx[i]].P().Z()) {
					minZ   = cm.vert[curHole.vertIdx[i]].P().Z();
					idxMin = curHole.vertIdx[i];
				}
			}
			avgZ /= curHole.vertIdx.size();
			cm.vert[idxMin].C() = Color4b::Yellow;

			cm.vert[188].C() = Color4b::Yellow;
			fillInSegment(m, cm, cm.vert[idxMin].P(), cm.vert[188].P(), boundaryPoints);

			//rotateInverse(md, rotateM);


			// Mesh refinement
			if (holeRefine) {
				float delta = 1.55; // Density factor

				unordered_map<int, float> sigma;
				unordered_map<int, int>   revPos;
				float                     total_len = 0;
				for (int v_id = 0; v_id < curHole.vertIdx.size(); v_id++) {
					revPos[curHole.vertIdx[v_id]] = v_id;
					int idx    = curHole.vertIdx[v_id];
					sigma[idx] = calcSigma(cm, idx);
					total_len += sigma[idx];
				}
				total_len /= curHole.vertIdx.size();
				for (int v_id = 0; v_id < curHole.vertIdx.size(); v_id++) {
					int idx    = curHole.vertIdx[v_id];
					sigma[idx] = total_len;
				}


				vector<float>   centroidSigma;
				vector<Point3m> centroid;
				for (int tri_id = 0; tri_id < triangles.size(); tri_id++) {
					Triangle curTri = triangles[tri_id];
					Point3m  curCentroid =
						(cm.vert[curTri.a].P() + cm.vert[curTri.b].P() + cm.vert[curTri.c].P()) /
						3.0;
					float curCentroidSigma =
						(sigma[curTri.a] + sigma[curTri.b] + sigma[curTri.c]) / 3.0;
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
							((delta * Distance(cm.vert[curTri.a].P(), cm.vert[curTri.b].P()) >
							  total_len) &&
							 (delta * Distance(cm.vert[curTri.a].P(), cm.vert[curTri.c].P()) >
							  total_len) &&
							 (delta * Distance(cm.vert[curTri.b].P(), cm.vert[curTri.b].P()) >
							  total_len))) {
							created                   = true;
							dividedTriangles.push_back(tri_id);
						}
					}
					divideCentroid(cm, triangles, dividedTriangles, centroid, centroidSigma, sigma, swapEdges, patchBit);

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
									cm,
									triangles,
									curTri.a,
									curTri.b,
									centroid,
									sigma,
									centroidSigma);
							else if (abs(revPos[cm.vert[curTri.a].Index()] - revPos[cm.vert[curTri.b].Index()]) != 1) {
								isSwapped = relaxEdge(
									cm,
									triangles,
									curTri.a,
									curTri.b,
									centroid,
									sigma,
									centroidSigma);
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
									cm,
									triangles,
									curTri.b,
									curTri.c,
									centroid,
									sigma,
									centroidSigma);
							else if (abs(revPos[cm.vert[curTri.c].Index()] - revPos[cm.vert[curTri.b].Index()]) != 1) {
								isSwapped = relaxEdge(
									cm,
									triangles,
									curTri.b,
									curTri.c,
									centroid,
									sigma,
									centroidSigma);
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
									cm,
									triangles,
									curTri.c,
									curTri.a,
									centroid,
									sigma,
									centroidSigma);
							else if (abs(revPos[cm.vert[curTri.c].Index()] - revPos[cm.vert[curTri.a].Index()]) != 1) {
								isSwapped = relaxEdge(
									cm,
									triangles,
									curTri.c,
									curTri.a,
									centroid,
									sigma,
									centroidSigma);
							}
						}
						cnt++;
					}
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
								cm.vert[i].P() = Point3m(update(i, 0), update(i, 1), update(i, 2));
							}
					}
				}
				for (int i = 0; i < curHole.vertIdx.size(); i++)
					if (cm.vert[curHole.vertIdx[i]].IsUserBit(patchBit))
						cm.vert[curHole.vertIdx[i]].ClearUserBit(patchBit);
			}
		}

	} break;
	default: wrongActionCalled(filter);
	}
	return std::map<std::string, QVariant>();
}
MESHLAB_PLUGIN_NAME_EXPORTER(QhullPlugin)
