#pragma once

#pragma unmanaged 

//Eigen 
#include "./3rdparty/Eigen/Dense"
#include "./3rdparty/Eigen/Geometry"
#include "./3rdparty/Eigen/Sparse"

typedef Eigen::Vector2i EVec2i;
typedef Eigen::Vector2d EVec2d;

#include <vector>
using namespace std;



//実数解を一つ返す
//solve x^3 + a x^2 + b x + c
inline double getRealSolutionOfCubicFunc
(
	const double a,
	const double b,
	const double c
)
{
	const double p = b - a * a / 3;
	const double q = 2 * a*a*a / 27 - a*b/3 + c;

	const double tmp = q*q/4 + p*p*p / 27;
	if( tmp < 0) return 0;

	const double mTmp = - q / 2 + sqrt(tmp);
	const double nTmp = - q / 2 - sqrt(tmp);
	const double m = (mTmp < 0) ?  - pow( - mTmp, 1/3.0) :  pow(mTmp, 1/3.0);
	const double n = (nTmp < 0) ?  - pow( - nTmp, 1/3.0) :  pow(nTmp, 1/3.0);

	return - a / 3 + m + n;
}

inline void Trace(EVec2d &p)
{
	printf( "%f %f\n", p[0], p[1]);
}

inline double TriangleArea(const EVec2d &p0, const EVec2d &p1, const EVec2d &p2)
{
	EVec2d v1 = p1-p0;
	EVec2d v2 = p2-p0;
	double v1_len = v1.norm() ;
	double v2_len = v2.norm() ;
	double cosT = v1.dot(v2) / (v1_len * v2_len);
	double area = v1_len * v2_len * sqrt( 1.0 - cosT * cosT);


	//printf( "area== %f %f\n", cosT, area);
	//Trace(v1);
	//Trace(v2);
	return area;
}






typedef Eigen::SparseMatrix<double>    ESpMat;
typedef Eigen::Triplet<double>         ETrip ;



//   solve 
//   2  3  0  0  0     x1       8 
//   3  0  4  0  6     x2      45
//   0 -1 -3  2  0     x3    = -3
//   0  0  1  0  0     x4       3
//   0  4  2  0  1     x5      19
inline void EigenSparseMatPractice()
{
	//prepare field 
	ESpMat A(5,5);
	Eigen::VectorXd b(5);

	//fill A
	vector< Eigen::Triplet<double> > entries; //{row, col, val}
	entries.push_back( Eigen::Triplet<double>(0,0, 2) );
	entries.push_back( Eigen::Triplet<double>(1,0, 3) );

	entries.push_back( Eigen::Triplet<double>(0,1, 3) );
	entries.push_back( Eigen::Triplet<double>(2,1,-1) );
	entries.push_back( Eigen::Triplet<double>(4,1, 4) );

	entries.push_back( Eigen::Triplet<double>(2,3, 2) );
	
	entries.push_back( Eigen::Triplet<double>(1,4, 6) );
	entries.push_back( Eigen::Triplet<double>(4,4, 1) );

	entries.push_back( Eigen::Triplet<double>(1,2, 4) );
	entries.push_back( Eigen::Triplet<double>(2,2,-3) );
	entries.push_back( Eigen::Triplet<double>(3,2, 1) );
	entries.push_back( Eigen::Triplet<double>(4,2, 2) );


	A.setFromTriplets( entries.begin(), entries.end() );

	// fill b
	b[0] = 8;
	b[1] =45;
	b[2] =-3;
	b[3] = 3;
	b[4] =19;

	// solve Ax = b
	Eigen::SparseLU<ESpMat> LU(A);  
	Eigen::VectorXd x = LU.solve(b);
	
	//printf("%f %f %f %f %f\n", x[0], x[1], x[2], x[3], x[4]);
	return;
}






















// input : cps (2D control points,        size:N)
// output: C1  (2D Bezier control points, size:N)
// output:λi  (weight to compute Bezier cp, size:N)
//
// note
// 各制御点 cps[i] に対して、2次ベジエ制御点 Ci0, Ci1, Ci2を生成する
//
// ただし，
// C_i1 は上記の C1[i]に対応
// C_i0, C_i+1,0 については、『C_i0 = C_i+1,0 = (1-λi) Ci + λi Ci+1 』と計算できるのでλiとC1[i]のみ保持する


const int ITER_NUM = 15;	

inline void compute_kCurves
(
	const vector<EVec2d> &cps,

	int  sampleN, //sampling number of each curve segment
	vector<EVec2d> &controlPoints, 	
	vector<EVec2d> &points	
)
{
	if( cps.size() < 3 ) return;
	if( sampleN < 3) sampleN = 3;

	const int N = (int)cps.size();
	vector<double> lambda(N, 0.5), ti(N);
	vector<EVec2d> Ci0(N), Ci1 = cps, Ci2(N);

	for( int iter = 0; iter < ITER_NUM; ++iter)
	{
		//1. update lambda (use 0.5 for the first iteration (更新してもいい気もするけど..))
		if( iter != 0 ) 
		{
			for( int i = 0; i < N; ++i)
			{
				int nextI = (i+1)%N;
				double A = sqrt( TriangleArea(Ci0[i], Ci1[i    ], Ci1[nextI] ) );
				double B = sqrt( TriangleArea(Ci0[i], Ci1[i    ], Ci1[nextI] ) );
				double C = sqrt( TriangleArea(Ci1[i], Ci1[nextI], Ci2[nextI] ) );
				lambda[i] = A / (B + C);
			}
		}

		//2. update Ci0, Ci2
		for( int i = 0; i < N; ++i)
		{
			int nextI = (i+1)%N;
			Ci2[i] = Ci0[nextI] = (1 - lambda[i]) * Ci1[  i  ]  +  
				                       lambda[i]  * Ci1[nextI];
		}

		//3. update ti
		for( int i = 0; i < N; ++i)
		{
			//solve eq(4)
			EVec2d Ci2_Ci0 = Ci2[i] - Ci0[i];
			EVec2d Ci0_pi  = Ci0[i] - cps[i];
			double a = Ci2_Ci0.squaredNorm();
			double b = 3 * Ci2_Ci0.dot( Ci0_pi );
			double c = (3 * Ci0[i] - 2 * cps[i] - Ci2[i]).dot( Ci0_pi );
			double d = - Ci0_pi.squaredNorm();

			ti[i] = (a==0) ? 0.5 : getRealSolutionOfCubicFunc(b/a, c/a, d/a);
			//printf( " ti %f (%f %f %f %f)\n", ti[i], a,b,c,d);
		}

		//4. update C1
		ESpMat A(N,N);
		Eigen::VectorXd b1(N), b2(N);
		vector< Eigen::Triplet<double> > entries; //{row, col, val}

		for( int i = 0; i < N; ++i)
		{
			int nextI = (i+1) % N;
			int prevI = (i-1) < 0 ? N - 1 : i - 1;
			double t = ti[i];
			double a = (1-lambda[prevI]) * (1 - t) * (1 - t);
			double b = lambda[i] * t * t;
			double c =   lambda[prevI] * (1 - t) * (1 - t) + (2-(1+lambda[i]) * t) * t;

			entries.push_back( ETrip(i, prevI, a) );
			entries.push_back( ETrip(i, i    , c) );
			entries.push_back( ETrip(i, nextI, b) );
			b1[i] = cps[i][0];
			b2[i] = cps[i][1];
		}


		A.setFromTriplets(entries.begin(), entries.end());
		// solve Ax = b
		Eigen::SparseLU<ESpMat> LU(A);  
		Eigen::VectorXd x = LU.solve(b1);
		Eigen::VectorXd y = LU.solve(b2);

		for( int i = 0; i < N; ++i) Ci1[i] << x[i], y[i];

	}

	controlPoints.clear();
	for( int i = 0; i < N; ++i)
	{
		controlPoints.push_back( Ci0[i] );
		controlPoints.push_back( Ci1[i] );
	}


	points.clear();
	for( int i = 0; i < N; ++i)
	{
		const EVec2d &C0 = Ci0[i];
		const EVec2d &C1 = Ci1[i];
		const EVec2d &C2 = Ci2[i];
		for( int j=0; j < sampleN; ++j)
		{
			double t = j * 1.0 / sampleN;
			points.push_back(  (1-t)*(1-t)*C0 + 2*t*(1-t)*C1 +  t*t*C2 ); 
		}
	}
}

#pragma managed 



