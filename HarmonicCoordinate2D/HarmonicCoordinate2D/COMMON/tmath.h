#pragma once

#pragma unmanaged

//Eigen 
#include <Dense>
#include <Geometry> 

//stl
#include <iostream>
#include <vector>
#include <list>
#include <deque>

#pragma unmanaged

typedef Eigen::Vector2i EVec2i;
typedef Eigen::Vector2d EVec2d;
typedef Eigen::Vector2f EVec2f;

typedef Eigen::Vector3i EVec3i;
typedef Eigen::Vector3d EVec3d;
typedef Eigen::Vector3f EVec3f;

typedef Eigen::Vector4i EVec4i;
typedef Eigen::Vector4d EVec4d;
typedef Eigen::Vector4f EVec4f;

typedef Eigen::Matrix2d EMat2d;
typedef Eigen::Matrix3f EMat3f;
typedef Eigen::Matrix4f EMat4f;
typedef Eigen::Matrix3d EMat3d;
typedef Eigen::Matrix4d EMat4d;

typedef Eigen::VectorXf EVecXf;
typedef Eigen::MatrixXf EMatXf;

#ifndef byte
typedef unsigned char byte;
#endif

#ifndef max3
#define max3(a,b,c)      std::max((a), std::max( (b),(c) ))
#endif

#ifndef min3
#define min3(a,b,c)      std::min((a), std::min( (b),(c) ))
#endif


#ifndef max3id
#define max3id(a,b,c)    ((a>=b&&a>=c)?0:(b>=c)?1:2)
#endif




inline void Trace(const EVec3d &v)
{
  std::cout << v[0] << " " << v[1] << " " << v[2] << "\n";
}


inline void Trace(const EVec3f &v)
{
  std::cout << v[0] << " " << v[1] << " " << v[2] << "\n";
}
inline void Trace(const EVec2f &v)
{
  std::cout << v[0] << " " << v[1] << "\n";
}


inline void Trace(const EVec2d &v)
{
  std::cout << v[0] << " " << v[1] << "\n";
}

inline void Trace(const EMat4d &m)
{
  for (int i = 0; i < 4; i++) 
    std::cout << m(i,0) << " " << m(i,1) << " " << m(i,2) << " " << m(i,3) << "\n";
}

inline void Trace(const EMat3d &m)
{
  for (int i = 0; i < 3; i++)
    std::cout << m(i,0) << " " << m(i,1) << " " << m(i,2) << "\n";

}

inline void Trace(const EMat3f &m)
{
  for (int i = 0; i < 3; i++) 
    std::cout << m(i,0) << " " << m(i,1) << " " << m(i,2) << "\n";
}






inline float t_distRayToPoint(const EVec3f &rayP, const EVec3f &rayD, const EVec3f &p)
{
  float t = (p - rayP).dot(rayD) / rayD.dot(rayD);
  return (rayP + t * rayD - p).norm();
}


inline float t_distRayToPoint(const EVec3f &rayP, const EVec3f &rayD, const EVec3f &p, float &t)
{
  t = (p - rayP).dot(rayD) / rayD.dot(rayD);
  return (rayP + t * rayD - p).norm();
}



inline float t_distPointToLineSegment_sq(const EVec3f &p,
  const EVec3f &lineP0,
  const EVec3f &lineP1)
{
  EVec3f dir = (lineP1 - lineP0);
  float len = dir.norm();
  if (len == 0) return (p - lineP0).squaredNorm();

  dir /= len;
  float t = (p - lineP0).dot(dir);

  if (t < 0) return (p - lineP0).squaredNorm();
  if (t > 1) return (p - lineP1).squaredNorm();

  EVec3f v = lineP0 + t*dir - p;
  return v.squaredNorm();
}


inline float t_distPointToLineSegment(const EVec3f &p,
  const EVec3f &lineP0,
  const EVec3f &lineP1)
{
  return sqrt(t_distPointToLineSegment_sq(p, lineP0, lineP1));
}





inline float t_distPointToLineSegment(
  const EVec2f &p     ,
  const EVec2f &lineP0,
  const EVec2f &lineP1,
  float &t)
{
  EVec2f dir = (lineP1 - lineP0);
  float len = dir.norm();

  if (len == 0)
  {
    t = 0;
    return (p - lineP0).norm();
  }

  t = (p - lineP0).dot(dir) / (len*len);

  if (t < 0) return (p - lineP0).norm();
  if (t > 1) return (p - lineP1).norm();

  return (lineP0 + t*dir - p).norm();
}






/*

inline double t_distRayToPointD(const EVec3d &rayP, const EVec3d &rayD, const EVec3d &p)
{
  double t = (p - rayP).dot(rayD) / rayD.dot(rayD);
  return (rayP + t * rayD - p).norm();
}

inline bool t_intersectRayToTriangle
(
  const EVec3d &rayP,
  const EVec3d &rayD,
  const EVec3d &x0,
  const EVec3d &x1,
  const EVec3d &x2,
  EVec3d &pos
  )
{
  Eigen::Matrix3d A;
  A << x1[0] - x0[0], x2[0] - x0[0], -rayD[0],
    x1[1] - x0[1], x2[1] - x0[1], -rayD[1],
    x1[2] - x0[2], x2[2] - x0[2], -rayD[2];

  EVec3d stu = A.inverse()*(rayP - x0);

  if (0 <= stu[0] && stu[0] <= 1 &&
    0 <= stu[1] && stu[1] <= 1 &&
    0 <= stu[0] + stu[1] && stu[0] + stu[1] <= 1)
  {
    pos = rayP + stu[2] * rayD;
    return true;
  }

  return false;
}
*/

inline bool t_intersectRayToTriangle
(
  const EVec3f &rayP,
  const EVec3f &rayD,
  const EVec3f &x0,
  const EVec3f &x1,
  const EVec3f &x2,
  EVec3f &pos
)
{
  Eigen::Matrix3f A;
  A << x1[0] - x0[0], x2[0] - x0[0], -rayD[0],
    x1[1] - x0[1], x2[1] - x0[1], -rayD[1],
    x1[2] - x0[2], x2[2] - x0[2], -rayD[2];

  EVec3f stu = A.inverse()*(rayP - x0);

  if (0 <= stu[0] && stu[0] <= 1 &&
    0 <= stu[1] && stu[1] <= 1 &&
    0 <= stu[0] + stu[1] && stu[0] + stu[1] <= 1)
  {
    pos = rayP + stu[2] * rayD;
    return true;
  }

  return false;
}



inline bool t_intersectRayToQuad
(
  const EVec3f &rayP,
  const EVec3f &rayD,

  const EVec3f &x0,
  const EVec3f &x1,
  const EVec3f &x2,
  const EVec3f &x3,

  EVec3f &pos
)
{
  if (t_intersectRayToTriangle(rayP, rayD, x0, x1, x2, pos)) return true;
  if (t_intersectRayToTriangle(rayP, rayD, x0, x2, x3, pos)) return true;
  return false;
}


template<class T>
void t_getMaxMinOfArray(const int N, const T* src, T& minV, T& maxV)
{
  if (N == 0) return;

  minV = src[0];
  maxV = src[0];
  for (int i = 0; i < N; ++i)
  {
    if (src[i] < minV) minV = src[i];
    if (src[i] > maxV) maxV = src[i];
  }

}




template<class T>
inline T t_crop(const T &minV, const T &maxV, const T &v)
{
  return std::min(maxV, std::max(minV, v));
}




inline bool t_bInWindow3D(const EVec3f &minW, const EVec3f &maxW, const EVec3f &p, float offset = 0)
{
  return  minW[0] - offset <= p[0] && p[0] <= maxW[0] + offset &&
    minW[1] - offset <= p[1] && p[1] <= maxW[1] + offset &&
    minW[2] - offset <= p[2] && p[2] <= maxW[2] + offset;

}





inline float t_dist(const EVec2f &p1, const EVec2f &p2)
{
  return sqrt((p1[0] - p2[0]) * (p1[0] - p2[0]) +
              (p1[1] - p2[1]) * (p1[1] - p2[1]));
}


inline double t_dist(const EVec3d &p1, const EVec3d &p2)
{
  return sqrt((p1[0] - p2[0]) * (p1[0] - p2[0]) +
    (p1[1] - p2[1]) * (p1[1] - p2[1]) +
    (p1[2] - p2[2]) * (p1[2] - p2[2]));
}


inline float t_dist(const EVec3f &p1, const EVec3f &p2)
{
  return sqrt((p1[0] - p2[0]) * (p1[0] - p2[0]) +
    (p1[1] - p2[1]) * (p1[1] - p2[1]) +
    (p1[2] - p2[2]) * (p1[2] - p2[2]));
}


inline float t_dist_sq(const EVec3f &p1, const EVec3f &p2)
{
  return (p1[0] - p2[0]) * (p1[0] - p2[0]) +
         (p1[1] - p2[1]) * (p1[1] - p2[1]) +
         (p1[2] - p2[2]) * (p1[2] - p2[2]);
}




template<class T>
inline T t_calcGravityCenter(const std::vector<T> &vs)
{
  T gc(0, 0, 0);
  for (const auto &p : vs) gc += p;
  gc /= (float)vs.size();
  return gc;

}


template<class T>
inline void t_curveSmoothing(int TIME, std::vector<T> &curve)
{
  const int N = (int)curve.size();
  if (N < 3) return;

  for (int time = 0; time < TIME; ++time)
  {
    std::vector<T> c(curve.size());

    c[0] = 0.25*curve[N - 1] + 0.5*curve[0] + 0.25*curve[1];
    for (int i = 1; i < N - 1; ++i) c[i] = 0.25*curve[i - 1] + 0.5*curve[i] + 0.25*curve[i + 1];
    c[N - 1] = 0.25*curve[N - 2] + 0.5*curve[N - 1] + 0.25*curve[0];
    curve = c;

  }
}






//calc length of a polyline
inline float t_verts_Length(const std::vector<EVec3f> &verts, bool bClosed = false)
{
  const int N = (int)verts.size();
  float d = 0;
  if (bClosed && N >= 2) d += t_dist(verts.back(), verts.front());

  for (int i = 1; i < N; ++i) d += t_dist(verts[i], verts[i - 1]);

  return d;
}




inline void t_verts_Translate(const EVec3f vec, std::vector<EVec3f> &verts)
{
  for (auto& it : verts) it += vec;
}
//
//inline float t_verts_MultMat( const EMat4f &M, std::vector<EVec3f> &verts)
//{
//	for( auto& it : verts ) it *= M;
//}
//

inline void t_verts_Rotate(const Eigen::AngleAxis<float> &R, const EVec3f &center, std::vector<EVec3f> &verts)
{
  for (auto& it : verts) it -= center;
  for (auto& it : verts) it = R * it;
  for (auto& it : verts) it += center;
}




//----------------------------------------------------------------------
// resampling a polyline into "n" section with equal interval
// This function returns n+1 points
// ---------------------------------------------------------------------
inline void t_verts_ResampleEqualInterval
(
  const int n,
  const std::vector<EVec3f> &verts,
  std::vector<EVec3f> &result
)
{
  result.clear();
  if (verts.size() < 2) { result.resize(n + 1); return; }

  const float stepD = t_verts_Length(verts) / (float)n;

  if (stepD == 0) {
    result.resize(n + 1);
    return;
  }

  result.push_back(verts[0]);
  float distance = 0;

  EVec3f vec, pivot = verts[0];

  for (int index = 1; index < (int)verts.size();)
  {
    distance += t_dist(verts[index], pivot);

    if (distance >= stepD)//steo over
    {
      vec = pivot - verts[index];
      vec *= (distance - stepD) / vec.norm();

      pivot = verts[index] + vec;
      result.push_back(pivot);
      distance = 0;
    }
    else
    {
      pivot = verts[index];
      ++index;
    }
  }
  if (result.size() != n + 1) result.push_back(verts.back());
}



inline void t_verts_Smoothing(std::vector<EVec3f> &verts)
{
  const int N = (int)verts.size();
  if (N < 2) return;

  std::vector<EVec3f> result = verts;

  for (int i = 1; i < N - 1; ++i) result[i] = 0.5 * verts[i] + 0.25 * verts[i - 1] + 0.25 * verts[i + 1];
  verts = result;
}


inline void t_verts_Smoothing(const int times, std::vector<EVec3f> &verts)
{
  for (int i = 0; i < times; ++i) t_verts_Smoothing(verts);
}




//	  | a b | |s|    w1
//    | c d | |t|  = w2
inline bool t_solve2by2LinearEquation(const double a, const double b,
  const double c, const double d, const double w1, const double w2,
  double &s, double &t)
{
  double det = (a*d - b*c);
  if (det == 0) return false;
  det = 1.0 / det;
  s = (d*w1 - b*w2) * det;
  t = (-c*w1 + a*w2) * det;
  return true;
}



static void t_calcBoundBox2D(const std::vector<EVec2f> &verts, EVec2f &BBmin, EVec2f &BBmax)
{
  BBmin << FLT_MAX, FLT_MAX;
  BBmax <<-FLT_MAX,-FLT_MAX;

  for (int i = 0, s = (int)verts.size(); i < s; ++i)
  {
    BBmin[0] = std::min(BBmin[0], verts[i][0]);
    BBmin[1] = std::min(BBmin[1], verts[i][1]);
    BBmax[0] = std::max(BBmax[0], verts[i][0]);
    BBmax[1] = std::max(BBmax[1], verts[i][1]);
  }
}



static void t_calcBoundBox2D(const std::vector<EVec2i> &verts, EVec2i &BBmin, EVec2i &BBmax)
{
  BBmin << INT_MAX, INT_MAX;
  BBmax << INT_MIN, INT_MIN;

  for (int i = 0; i < (int)verts.size(); ++i)
  {
    BBmin[0] = std::min(BBmin[0], verts[i][0]);
    BBmin[1] = std::min(BBmin[1], verts[i][1]);
    BBmax[0] = std::max(BBmax[0], verts[i][0]);
    BBmax[1] = std::max(BBmax[1], verts[i][1]);
  }
}



static void t_calcBoundBox3D(const int N, const EVec3f* verts, EVec3f &BBmin, EVec3f &BBmax)
{
  BBmin << FLT_MAX, FLT_MAX, FLT_MAX;
  BBmax << -FLT_MAX, -FLT_MAX, -FLT_MAX;

  for (int i = 0; i < N; ++i)
  {
    BBmin[0] = std::min(BBmin[0], (float)verts[i][0]);
    BBmin[1] = std::min(BBmin[1], (float)verts[i][1]);
    BBmin[2] = std::min(BBmin[2], (float)verts[i][2]);
    BBmax[0] = std::max(BBmax[0], (float)verts[i][0]);
    BBmax[1] = std::max(BBmax[1], (float)verts[i][1]);
    BBmax[2] = std::max(BBmax[2], (float)verts[i][2]);
  }
}







inline double t_calcAngleTwoVec(
  EVec3f v1,
  EVec3f v2,
  EVec3f axis
)
{
  double v1_len = v1.norm();
  double v2_len = v2.norm();
  if (v1_len == 0 || v2_len == 0) return 0;
  double cosT = v1.dot(v2) / v1_len / v2_len;
  if (cosT >= 1.0) cosT = 1;
  if (cosT <= -1.0) cosT = -1;

  if (axis.dot(v1.cross(v2)) >= 0) return acos(cosT);
  else return -acos(cosT);
}




static double t_calcTriangleArea(
  const EVec3f &x0,
  const EVec3f &x1,
  const EVec3f &x2
)
{
  EVec3f v1 = x1 - x0;
  EVec3f v2 = x2 - x0;
  double v1_len = v1.norm();
  double v2_len = v2.norm();
  if (v1_len == 0 || v2_len == 0) return 0;
  double cosT = v1.dot(v2) / v1_len / v2_len;
  return v1.norm() * v2.norm() * sqrt(1 - cosT*cosT);
}




//-------------------------------------------------------------------//
//---CONTOUR MATCHING ALGORITHM--------------------------------------//
// paper : Optimal surface reconstruction from planar contours
// implemented by Chika Tomiyama and Takashi Ijiri Oct. 2017
//-------------------------------------------------------------------//



// get shortest path in (2*pN x (qN+1) ) graph  
// from (pivPi, 0) to (pivPi + pN, qN)
// 
// edgeCostV and edgeCostH store edge capacity
//
static double c_findShortestPathByDP(
  const int pivPi,
  const int pN,
  const int qN,
  double **edgeCostV, //2pN x (qN+1) NotChanged
  double **edgeCostH, //2pN x (qN+1) NotChanged
  double **costArray, //2pN x (qN+1)
  int    **fromArray, //2pN x (qN+1)
  std::list<std::pair<int, int>> &path
)
{
  for (int p = 0; p < 2 * pN; ++p)
    for (int q = 0; q <= qN; ++q)
    {
      costArray[p][q] = DBL_MAX;
      fromArray[p][q] = -1;
    }

  path.clear();
  costArray[pivPi][0] = 0.0;
  fromArray[pivPi][0] = 0;

  for (int p = pivPi; p <= pivPi + pN; ++p)
  {
    for (int q = 0; q <= qN; ++q)
    {
      if (p == pivPi && q == 0) continue;
      double costH = (q == 0) ? DBL_MAX : costArray[p][q - 1] + edgeCostH[p][q - 1];
      double costV = (p == pivPi) ? DBL_MAX : costArray[p - 1][q] + edgeCostV[p - 1][q];
      fromArray[p][q] = (costV < costH) ? 1 : 2;
      costArray[p][q] = (costV < costH) ? costV : costH;
    }
  }

  int p = pivPi + pN;
  int q = qN;
  path.push_front(std::make_pair(p, q));

  while (true)
  {
    if (fromArray[p][q] == 1) p--;
    else q--;
    path.push_front(std::make_pair(p, q));

    if (p == pivPi && q == 0) break;
  }

  return costArray[pivPi + pN][qN];
}


enum POLY_MACH_METRIC
{
  PMM_AREA, //triangle area
  PMM_ARCLEN   //distance between matched vertices of P and Q
};


inline void c_polylineMatching
(
  const std::vector<EVec3f>   &P,
  const std::vector<EVec3f>   &Q,
  const POLY_MACH_METRIC metric,
  std::list<std::pair<int, int>>   &matchedIds
)
{
  matchedIds.clear();
  const int pNum = (int)P.size();
  const int qNum = (int)Q.size();
  if (qNum < 3 || pNum < 3) return;

  //initialzie edgeCost
  double  **edgeCostV = new double*[2 * pNum];
  double  **edgeCostH = new double*[2 * pNum]; //-1:init, 0:leftmost, 1:from Left, 2:from top
  double  **costArray = new double*[2 * pNum];
  int     **fromArray = new    int*[2 * pNum]; //-1:init, 0:start, 1:top, 2:from left  

  for (int p = 0; p < 2 * pNum; ++p)
  {
    edgeCostV[p] = new double[qNum + 1];
    edgeCostH[p] = new double[qNum + 1];
    costArray[p] = new double[qNum + 1];
    fromArray[p] = new int[qNum + 1];

    for (int q = 0; q <= qNum; ++q)
    {
      if (metric == PMM_AREA) {
        edgeCostV[p][q] = t_calcTriangleArea(P[p % pNum], Q[q % qNum], P[(p + 1) % pNum]);
        edgeCostH[p][q] = t_calcTriangleArea(P[p % pNum], Q[q % qNum], Q[(q + 1) % qNum]);
      }
      else
      {
        //edgeCostV[p][q] = (Q[q % qNum] - 0.5*P[p % pNum] - 0.5*P[(p + 1) % pNum]).norm(); 
        //edgeCostH[p][q] = (P[p % pNum] - 0.5*Q[q % qNum] - 0.5*Q[(q + 1) % qNum]).norm();
        edgeCostV[p][q] = (P[p % pNum] - Q[q % qNum]).norm() + (P[(p + 1) % pNum] - Q[q % qNum]).norm();
        edgeCostH[p][q] = (P[p % pNum] - Q[q % qNum]).norm() + (P[p % pNum] - Q[(q + 1) % qNum]).norm();
      }
      costArray[p][q] = DBL_MAX;
      fromArray[p][q] = -1;
    }
  }

  double minCost = DBL_MAX;

  for (int pivPi = 0; pivPi < pNum; ++pivPi)
  {
    std::list< std::pair<int, int> > pivPiPath;
    double pivPiCost = c_findShortestPathByDP(pivPi, pNum, qNum, edgeCostV, edgeCostH, costArray, fromArray, pivPiPath);
    if (pivPiCost < minCost)
    {
      matchedIds = pivPiPath;
      minCost = pivPiCost;
    }
  }

  for (int i = 0; i < pNum * 2; i++)
  {
    delete[] edgeCostV[i];
    delete[] edgeCostH[i];
    delete[] costArray[i];
    delete[] fromArray[i];
  }
  delete[] edgeCostV;
  delete[] edgeCostH;
  delete[] costArray;
  delete[] fromArray;
}




static void t_CalcBoundingBox(const std::vector<EVec2i> &verts, EVec2i &BBmin, EVec2i &BBmax)
{
  BBmin << INT_MAX, INT_MAX;
  BBmax << INT_MIN, INT_MIN;

  for (int i = 0; i < (int)verts.size(); ++i)
  {
    BBmin[0] = std::min(BBmin[0], verts[i][0]);
    BBmin[1] = std::min(BBmin[1], verts[i][1]);
    BBmax[0] = std::max(BBmax[0], verts[i][0]);
    BBmax[1] = std::max(BBmax[1], verts[i][1]);
  }
}



static void t_CalcBoundingBox(const int N, const EVec3f* verts, EVec3f &BBmin, EVec3f &BBmax)
{
  BBmin <<  FLT_MAX,  FLT_MAX,  FLT_MAX;
  BBmax << -FLT_MAX, -FLT_MAX, -FLT_MAX;

  for (int i = 0; i < N; ++i)
  {
    BBmin[0] = std::min(BBmin[0], (float)verts[i][0]);
    BBmin[1] = std::min(BBmin[1], (float)verts[i][1]);
    BBmin[2] = std::min(BBmin[2], (float)verts[i][2]);
    BBmax[0] = std::max(BBmax[0], (float)verts[i][0]);
    BBmax[1] = std::max(BBmax[1], (float)verts[i][1]);
    BBmax[2] = std::max(BBmax[2], (float)verts[i][2]);
  }
}

inline void t_CalcBoundingBox
(
  const std::vector<EVec3f> &verts,
  EVec3f &BBmin,
  EVec3f &BBmax
)
{
  BBmin << FLT_MAX, FLT_MAX, FLT_MAX;
  BBmax <<-FLT_MAX,-FLT_MAX,-FLT_MAX;
  
  for( const auto &v : verts ){
    BBmin[0] = std::min(BBmin[0], v[0]);
    BBmin[1] = std::min(BBmin[1], v[1]);
    BBmin[2] = std::min(BBmin[2], v[2]);
    BBmax[0] = std::max(BBmax[0], v[0]);
    BBmax[1] = std::max(BBmax[1], v[1]);
    BBmax[2] = std::max(BBmax[2], v[2]);
  }
}






//calculate the angle of two 3D vectors
//axisに対して右ねじ方向が正
inline double t_CalcAngle(
  const EVec3f &v1,
  const EVec3f &v2,
  const EVec3f &axis
)
{
  double v1_len = v1.norm();
  double v2_len = v2.norm();
  if (v1_len == 0 || v2_len == 0) return 0;
  double cosT = v1.dot(v2) / v1_len / v2_len;
  if (cosT >=  1.0) cosT = 1;
  if (cosT <= -1.0) cosT = -1;

  if (axis.dot(v1.cross(v2)) >= 0) return acos(cosT);
  else return -acos(cosT);
}



//2D座標系において左回りが正（右ねじ座標系でz(0,0,1)を軸にして右回りが正）
inline double t_CalcAngle(const EVec2d &d1, const EVec2d &d2)
{
	double l = d1.norm() * d2.norm();
	if( l == 0 ) return 0;

	double cosT = t_crop<double>(-1., 1., ( d1.dot(d2) ) / l);

	if( d1[0] * d2[1] - d1[1] * d2[0] >= 0)
    return  acos(cosT);
	else
    return -acos(cosT);
}







// grow distans transform volume (vol_dist)
// return minimum value of grown voxel
inline static float __DistTransGrowFronteer
(
    const int W, const int H, const int D,
    const EVec3f &pitch, 

    std::deque<EVec4i> &Q1, //input : current fronteer 
    std::deque<EVec4i> &Q2, //output: new fronteer

    float *vol_dist,
    byte  *vol_flg  
)
{
	const int WH = W*H, WHD = W*H*D;
	float pX   = pitch[0];
  float pY   = pitch[1];
  float pZ   = pitch[2];
	float pXY  = sqrt(pX*pX + pY*pY        );
	float pYZ  = sqrt(        pY*pY + pZ*pZ);
	float pXZ  = sqrt(pX*pX         + pZ*pZ);
	float pXYZ = sqrt(pX*pX + pY*pY + pZ*pZ);
	
	float min_dist = FLT_MAX; 

	while(!Q1.empty())
	{
		const int x = Q1.front()[0];
		const int y = Q1.front()[1];
		const int z = Q1.front()[2];
		const int I = Q1.front()[3];
		Q1.pop_front();

		int K ;
	  
    //1. vol[I]の26近傍を見て，距離が最小となるものを入れる
    //2. one ring分growさせる
		//calc Dist for voxel (x,y,z) 
    //unisotropic volumeで、seedが曲線状なのでこの実装のほうが多少正確
		K = I-1-W-WH; if(x>0  &&y>0  &&z>0	&& vol_flg[K]) vol_dist[I] = std::min( vol_dist[I], vol_dist[K] + pXYZ );
		K = I  -W-WH; if(       y>0  &&z>0  && vol_flg[K]) vol_dist[I] = std::min( vol_dist[I], vol_dist[K] + pYZ  );
		K = I+1-W-WH; if(x<W-1&&y>0  &&z>0  && vol_flg[K]) vol_dist[I] = std::min( vol_dist[I], vol_dist[K] + pXYZ );
		K = I-1  -WH; if(x>0         &&z>0  && vol_flg[K]) vol_dist[I] = std::min( vol_dist[I], vol_dist[K] + pXZ  );
		K = I    -WH; if(              z>0  && vol_flg[K]) vol_dist[I] = std::min( vol_dist[I], vol_dist[K] + pZ   );
		K = I+1  -WH; if(x<W-1       &&z>0  && vol_flg[K]) vol_dist[I] = std::min( vol_dist[I], vol_dist[K] + pXZ  );
		K = I-1+W-WH; if(x>0  &&y<H-1&&z>0  && vol_flg[K]) vol_dist[I] = std::min( vol_dist[I], vol_dist[K] + pXYZ );
		K = I  +W-WH; if(       y<H-1&&z>0  && vol_flg[K]) vol_dist[I] = std::min( vol_dist[I], vol_dist[K] + pYZ  );
		K = I+1+W-WH; if(x<W-1&&y<H-1&&z>0  && vol_flg[K]) vol_dist[I] = std::min( vol_dist[I], vol_dist[K] + pXYZ );
		
		K = I-1-W   ; if(x>0  &&y>0         && vol_flg[K]) vol_dist[I] = std::min( vol_dist[I], vol_dist[K] + pXY );
		K = I  -W   ; if(       y>0         && vol_flg[K]) vol_dist[I] = std::min( vol_dist[I], vol_dist[K] + pY  );
		K = I+1-W   ; if(x<W-1&&y>0         && vol_flg[K]) vol_dist[I] = std::min( vol_dist[I], vol_dist[K] + pXY );
		K = I-1     ; if(x>0                && vol_flg[K]) vol_dist[I] = std::min( vol_dist[I], vol_dist[K] + pX  );
		K = I+1     ; if(x<W-1              && vol_flg[K]) vol_dist[I] = std::min( vol_dist[I], vol_dist[K] + pX  );
		K = I-1+W   ; if(x>0  &&y<H-1       && vol_flg[K]) vol_dist[I] = std::min( vol_dist[I], vol_dist[K] + pXY );
		K = I  +W   ; if(       y<H-1       && vol_flg[K]) vol_dist[I] = std::min( vol_dist[I], vol_dist[K] + pY  );
		K = I+1+W   ; if(x<W-1&&y<H-1       && vol_flg[K]) vol_dist[I] = std::min( vol_dist[I], vol_dist[K] + pXY );
                                                                     
		K = I-1-W+WH; if(x>0  &&y>0  &&z<D-1&& vol_flg[K]) vol_dist[I] = std::min( vol_dist[I], vol_dist[K] + pXYZ );
		K = I  -W+WH; if(       y>0  &&z<D-1&& vol_flg[K]) vol_dist[I] = std::min( vol_dist[I], vol_dist[K] + pYZ  );
		K = I+1-W+WH; if(x<W-1&&y>0  &&z<D-1&& vol_flg[K]) vol_dist[I] = std::min( vol_dist[I], vol_dist[K] + pXYZ );
		K = I-1  +WH; if(x>0         &&z<D-1&& vol_flg[K]) vol_dist[I] = std::min( vol_dist[I], vol_dist[K] + pXZ  );
		K = I    +WH; if(              z<D-1&& vol_flg[K]) vol_dist[I] = std::min( vol_dist[I], vol_dist[K] + pZ   );
		K = I+1  +WH; if(x<W-1       &&z<D-1&& vol_flg[K]) vol_dist[I] = std::min( vol_dist[I], vol_dist[K] + pXZ  );
		K = I-1+W+WH; if(x>0  &&y<H-1&&z<D-1&& vol_flg[K]) vol_dist[I] = std::min( vol_dist[I], vol_dist[K] + pXYZ );
		K = I  +W+WH; if(       y<H-1&&z<D-1&& vol_flg[K]) vol_dist[I] = std::min( vol_dist[I], vol_dist[K] + pYZ  );
		K = I+1+W+WH; if(x<W-1&&y<H-1&&z<D-1&& vol_flg[K]) vol_dist[I] = std::min( vol_dist[I], vol_dist[K] + pXYZ );

		//grow fronteer to deque
		K = I-1-W-WH; if(x>0  &&y>0  &&z>0	&&!vol_flg[K]) { Q2.push_back(EVec4i(x-1, y-1, z-1, K)); vol_flg[K] = 1;}	
		K = I  -W-WH; if(       y>0  &&z>0  &&!vol_flg[K]) { Q2.push_back(EVec4i(x  , y-1, z-1, K)); vol_flg[K] = 1;}	
		K = I+1-W-WH; if(x<W-1&&y>0  &&z>0  &&!vol_flg[K]) { Q2.push_back(EVec4i(x+1, y-1, z-1, K)); vol_flg[K] = 1;}	
		K = I-1  -WH; if(x>0         &&z>0  &&!vol_flg[K]) { Q2.push_back(EVec4i(x-1, y  , z-1, K)); vol_flg[K] = 1;}	
		K = I    -WH; if(              z>0  &&!vol_flg[K]) { Q2.push_back(EVec4i(x  , y  , z-1, K)); vol_flg[K] = 1;}
		K = I+1  -WH; if(x<W-1       &&z>0  &&!vol_flg[K]) { Q2.push_back(EVec4i(x+1, y  , z-1, K)); vol_flg[K] = 1;}
		K = I-1+W-WH; if(x>0  &&y<H-1&&z>0  &&!vol_flg[K]) { Q2.push_back(EVec4i(x-1, y+1, z-1, K)); vol_flg[K] = 1;}
		K = I  +W-WH; if(       y<H-1&&z>0  &&!vol_flg[K]) { Q2.push_back(EVec4i(x  , y+1, z-1, K)); vol_flg[K] = 1;}
		K = I+1+W-WH; if(x<W-1&&y<H-1&&z>0  &&!vol_flg[K]) { Q2.push_back(EVec4i(x+1, y+1, z-1, K)); vol_flg[K] = 1;}
		
		K = I-1-W   ; if(x>0  &&y>0         &&!vol_flg[K]) { Q2.push_back(EVec4i(x-1, y-1, z  , K)); vol_flg[K] = 1;}
		K = I  -W   ; if(       y>0         &&!vol_flg[K]) { Q2.push_back(EVec4i(x  , y-1, z  , K)); vol_flg[K] = 1;}
		K = I+1-W   ; if(x<W-1&&y>0         &&!vol_flg[K]) { Q2.push_back(EVec4i(x+1, y-1, z  , K)); vol_flg[K] = 1;}
		K = I-1     ; if(x>0                &&!vol_flg[K]) { Q2.push_back(EVec4i(x-1, y  , z  , K)); vol_flg[K] = 1;}
		K = I+1     ; if(x<W-1              &&!vol_flg[K]) { Q2.push_back(EVec4i(x+1, y  , z  , K)); vol_flg[K] = 1;}
		K = I-1+W   ; if(x>0  &&y<H-1       &&!vol_flg[K]) { Q2.push_back(EVec4i(x-1, y+1, z  , K)); vol_flg[K] = 1;}
		K = I  +W   ; if(       y<H-1       &&!vol_flg[K]) { Q2.push_back(EVec4i(x  , y+1, z  , K)); vol_flg[K] = 1;}
		K = I+1+W   ; if(x<W-1&&y<H-1       &&!vol_flg[K]) { Q2.push_back(EVec4i(x+1, y+1, z  , K)); vol_flg[K] = 1;}

		K = I-1-W+WH; if(x>0  &&y>0  &&z<D-1&&!vol_flg[K]) { Q2.push_back(EVec4i(x-1, y-1, z+1, K)); vol_flg[K] = 1;}
		K = I  -W+WH; if(       y>0  &&z<D-1&&!vol_flg[K]) { Q2.push_back(EVec4i(x  , y-1, z+1, K)); vol_flg[K] = 1;}
		K = I+1-W+WH; if(x<W-1&&y>0  &&z<D-1&&!vol_flg[K]) { Q2.push_back(EVec4i(x+1, y-1, z+1, K)); vol_flg[K] = 1;}
		K = I-1  +WH; if(x>0         &&z<D-1&&!vol_flg[K]) { Q2.push_back(EVec4i(x-1, y  , z+1, K)); vol_flg[K] = 1;}
		K = I    +WH; if(              z<D-1&&!vol_flg[K]) { Q2.push_back(EVec4i(x  , y  , z+1, K)); vol_flg[K] = 1;}
		K = I+1  +WH; if(x<W-1       &&z<D-1&&!vol_flg[K]) { Q2.push_back(EVec4i(x+1, y  , z+1, K)); vol_flg[K] = 1;}
		K = I-1+W+WH; if(x>0  &&y<H-1&&z<D-1&&!vol_flg[K]) { Q2.push_back(EVec4i(x-1, y+1, z+1, K)); vol_flg[K] = 1;}
		K = I  +W+WH; if(       y<H-1&&z<D-1&&!vol_flg[K]) { Q2.push_back(EVec4i(x  , y+1, z+1, K)); vol_flg[K] = 1;}
		K = I+1+W+WH; if(x<W-1&&y<H-1&&z<D-1&&!vol_flg[K]) { Q2.push_back(EVec4i(x+1, y+1, z+1, K)); vol_flg[K] = 1;}

		min_dist = std::min( min_dist, vol_dist[I] );
	}
	return min_dist;
}



inline void CalcDistanceTransform(
  const EVec3f &pitch,
  const EVec3i &reso ,
  const float  max_dist, //この値以上の領域は計算しない
  const std::vector<EVec4i> &seed_voxels, //xyzI
  
  float* vol_dt //should be allocated
)
{
  const int W = reso[0];
  const int H = reso[1];
  const int D = reso[2];
	const int WHD = reso[0] * reso[1] * reso[2];
	
  //vol_flg -- 0:yet, 1:visited. 
	byte *vol_flg = new byte[WHD];

	for(int i = 0; i < WHD; ++i)
	{
		vol_dt [i] = FLT_MAX;
		vol_flg[i] = 0;
	}
	
	std::deque<EVec4i> Q1, Q2;

	for (const auto &p : seed_voxels)
	{
		Q1.push_back(p);
		vol_dt [ p[3] ] = 0.0;
		vol_flg[ p[3] ] = 1  ;
	}

	while( !Q1.empty() || !Q2.empty() )
	{
		float min_dist;
		min_dist = __DistTransGrowFronteer(W,H,D,pitch, Q1,Q2, vol_dt, vol_flg);
		if( min_dist > max_dist) break;
		min_dist = __DistTransGrowFronteer(W,H,D,pitch, Q2,Q1, vol_dt, vol_flg);			
		if( min_dist > max_dist) break;
	}

	delete[] vol_flg;
  std::cout << "...Distance Transform Done\n";
}




#pragma managed

