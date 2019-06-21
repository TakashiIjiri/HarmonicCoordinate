#include "stdafx.h"
#include "HarmCoord3D.h"


#pragma unmanaged


using namespace std;

//todo generate map
//todo generate constratin 
//todo perform smoothing 
//todo perform smoothing in hierarchical order 
//todo compute hc for all vertices






//note : map[i,j,k] corresponds to [ m_map_pos + m_map_pitch * (i+0.5,j+0.5,k+0.5) ]
//note : (x,y,z) in 3d space corresponds to map[i,j,k] 
//              (i.j,k) = (int)( ((x,y,z) - m_map_pos) / m_map_pitch ); 

//init map from cage
void HarmCoord3D::InitMap()
{
  //calc map configuration
  EVec3f bbmin, bbmax;
  m_cage.getBoundBox(bbmin, bbmax);
  EVec3f bbsize = bbmax - bbmin;

  int MAX_RESO = 128;
  m_map_pitch = max3( bbsize[0], bbsize[1], bbsize[2]) / (MAX_RESO - 2);
  m_map_pos  << bbmin[0] - m_map_pitch, bbmin[1] - m_map_pitch, bbmin[2] - m_map_pitch; 
  
  // resolutionÇ≈ÇÕÅAêÿÇËè„Ç∞Ç≈âëúìxÇåàíËÇ∑ÇÈ
  m_map_reso << (int)(bbsize[0] / m_map_pitch + 0.9999f) + 2,
                (int)(bbsize[1] / m_map_pitch + 0.9999f) + 2,
                (int)(bbsize[2] / m_map_pitch + 0.9999f) + 2 ;

  std::cout << "compute Harmonic coordinate volume configuration\n";

  Trace(m_map_reso);
  Trace(m_map_pos );
  std::cout << "pitch  =  " << m_map_pitch << "\n";

  //calc flag image (0:out,1:inside,2:boundary)
  m_map_flag = new byte[m_map_reso[0]*m_map_reso[1]*m_map_reso[2]];
  
  // cast ray in Y axis ( divide ZX plane )  
  TMesh cage = m_cage;
  cage.Translate( -m_map_pos );
  EVec3f pitch(m_map_pitch, m_map_pitch,m_map_pitch);
  genBinaryVolumeInTriangleMeshY(m_map_reso, pitch, cage, m_map_flag);

  const int W = m_map_reso[0];
  const int H = m_map_reso[1];
  const int D = m_map_reso[2], WH = W * H;

  for ( int z = 1; z < D-1; ++z)
  {
    for ( int y = 1; y < H-1; ++y)
    {
      for ( int x = 1; x < W-1; ++x)
      {
        int i= x + y * W + z * WH;
        if ( m_map_flag[i] != 1 ) continue;

        if( m_map_flag[i-1 ] == 0 || m_map_flag[i+1 ] == 0 || 
            m_map_flag[i-W ] == 0 || m_map_flag[i+W ] == 0 ||
            m_map_flag[i-WH] == 0 || m_map_flag[i+WH] == 0 )
        {
          m_map_flag[i] = 2; // set as boundary 
        }
      }
    }
  }
}



// h = x0 + s(x1-x0) + t(x2-x0)
// (h-p) * (x1-x0) = 0  --> d1d1 s + d1d2 t = d1 (p-x0)
// (h-p) * (x2-x0) = 0  --> d1d2 s + d2d2 t = d2 (p-x0)    
//ÇΩÇæÇµ d1 = (x1-x0)


static bool CalcBarycentricCoordinate(
  const EVec3f &p, 
  const EVec3f &x0, 
  const EVec3f &x1, 
  const EVec3f &x2, 
  float &s, float &t)
{
  EVec3f d1 = x1-x0;
  EVec3f d2 = x2-x0;
  return t_solve2by2LinearEquationF( 
    d1.dot(d1), d1.dot(d2), 
    d1.dot(d2), d2.dot(d2),  d1.dot(p-x0), 
                             d2.dot(p-x0), s, t);
}




static void CalcBoundBox3D(
  const EVec3f &x0, 
  const EVec3f &x1, 
  const EVec3f &x2, 
  float  offset, 
  EVec3f &bbmin, 
  EVec3f &bbmax)
{

  bbmin << min3(x0[0], x1[0], x2[0]) - offset, 
           min3(x0[1], x1[1], x2[1]) - offset, 
           min3(x0[2], x1[2], x2[2]) - offset;

  bbmax << max3(x0[0], x1[0], x2[0]) + offset, 
           max3(x0[1], x1[1], x2[1]) + offset, 
           max3(x0[2], x1[2], x2[2]) + offset;

}

static bool IsInBoundingBox(
  const EVec3f &p, 
  const EVec3f &bbmin, 
  const EVec3f &bbmax 
)
{
  return bbmin[0] <= p[0] && p[0] <= bbmax[0] && 
         bbmin[1] <= p[1] && p[1] <= bbmax[1] && 
         bbmin[2] <= p[2] && p[2] <= bbmax[2] ;
}






static void SetBoundaryCondition
(
  const int W, const int H, const int D,
  const float   pitch,
  const EVec3f &map_pos  ,
  const byte   *map_flg, // 0:out, 1:in, 2:boundary
  
  const EVec3f &x0, const EVec3f &x1, const EVec3f &x2,
  float *map_hc
)
{
  EVec3f bbmin, bbmax;
  CalcBoundBox3D(x0,x1,x2, pitch*1.5f, bbmin, bbmax);


  //set boundary condition value
  const int WH = W*H;

  for ( int z = 0; z < D; ++z)
  {
    for ( int y = 0; y < H; ++y)
    {
      for ( int x = 0; x < W; ++x)
      {
        int i = x + y * W + z * WH;
        if ( map_flg[i] != 2 ) continue;

        const EVec3f p( (x+0.5f)*pitch + map_pos[0], 
                        (y+0.5f)*pitch + map_pos[1], 
                        (z+0.5f)*pitch + map_pos[2] );
      
        if ( !IsInBoundingBox(p, bbmin, bbmax) ) continue;

        //bacycentric coordinateÇåvéZÇµéOäpå`ÇÃíÜÇ…ì¸Ç¡ÇƒÇ¢ÇΩÇÁíl 1-s-tÇìoò^
        float s,t;
        if( CalcBarycentricCoordinate(p,x0,x1,x2,s,t) && 
            -0.0f <= s && s <= 1.0f && -0.0f <= t && t <= 1.0f && 
            s + t <= 1.00f) 
        {
          map_hc[i] = 1-s-t;
        }
      }
    }
  }

}

static void SetBoundaryCondition
(
  const int W, const int H, const int D,
  const float   pitch,
  const EVec3f &map_pos  ,
  const byte   *map_flg, // 0:out, 1:in, 2:boundary
  
  const TMesh &cage,
  const int piv_cage_vi,
  float *map_hc
)
{
  //init map_hc
  memset( map_hc, 0, sizeof(float) * W*H*D );

  //ëŒè€Ç∆Ç»ÇÈéOäpå`Çã´äEêßñÒÇ∆ÇµÇƒï`âÊ
  for( int pi = 0; pi < cage.m_pSize; ++pi )
  {
    int* poly = cage.m_pPolys[pi].idx;
    EVec3f &x0 = cage.m_vVerts[poly[0]]; 
    EVec3f &x1 = cage.m_vVerts[poly[1]]; 
    EVec3f &x2 = cage.m_vVerts[poly[2]]; 

    if ( poly[0] == piv_cage_vi ) {
      SetBoundaryCondition(W,H,D, pitch, map_pos, map_flg, x0,x1,x2, map_hc);
    }else if( poly[1] == piv_cage_vi ) {
      SetBoundaryCondition(W,H,D, pitch, map_pos, map_flg, x1,x2,x0, map_hc);
    }else if( poly[2] == piv_cage_vi ) {
      SetBoundaryCondition(W,H,D, pitch, map_pos, map_flg, x2,x0,x1, map_hc);
    }else {
      continue;
    }
  }
}




static void LaplacianSmoothing(
  const int W, const int H, const int D, 
  const byte *map_flg , // 0:out, 1:in, 2:boundary
  const int NUN_ITER, 
  float *map_hc         //contains values only at boundary
)
{
  const int WH = W*H;

  float *tmp = new float[W*H*D];
  memcpy(tmp, map_hc, sizeof(float)*W*H*D);

  for ( int iter = 0; iter < NUN_ITER; ++iter)
  {
    for ( int z = 1; z < D-1; ++z)
      for ( int y = 1; y < H-1; ++y)
        for ( int x = 1; x < W-1; ++x)
        {
          int i = x + y * W + z * WH;
          if ( map_flg[i] != 1 ) continue;
          
          tmp[i] =  (map_hc[i- 1] + map_hc[i+ 1] + 
                     map_hc[i- W] + map_hc[i+ W] + 
                     map_hc[i-WH] + map_hc[i+WH]) / 6.0f;

        }
    memcpy(map_hc, tmp, sizeof(float)*W*H*D);
  }
  
  std::cout << "smoothign done--";
}


void HarmCoord3D::CalcHarmonicCoordinateMap()
{
  const int N = m_cage.m_vSize;
  const int W = m_map_reso[0];
  const int H = m_map_reso[1];
  const int D = m_map_reso[2], WH = W * H;

  for( int ci = 0; ci < N; ++ci)
  {
    float *map_hc = new float[W*H*D];
    m_map_hc.push_back(map_hc);
    SetBoundaryCondition(W,H,D, m_map_pitch, m_map_pos, m_map_flag, m_cage, ci, map_hc);
  }

#pragma omp parallel for
  for( int ci = 0; ci < N; ++ci)
    LaplacianSmoothing( W,H,D, m_map_flag, 5000,  m_map_hc[ci]);

}



// note : âÊëfÇÃíÜêS (0.5)ÇÃà íuÇ…âÊëfílÇ™Ç†ÇÈÇ∆Ç∑ÇÈ
// 
// posx =  (v[0] - m_map_pos[0]) / m_map_pitch Ç≈ÅAâÊëfäiéqè„Ç…Ç®ÇØÇÈ âÊëfà íuÇ™ìæÇÁÇÍÇÈ
// posx = posx-0.5 Ç≈â∫ê}ÇÃÇÊÇ§Ç…ÅCÇøÇÂÇ§Ç«0,1,2Ç…âÊëfílÇ™Ç†ÇÈãÛä‘Ç…Ç®ÇØÇÈà íuÇ…ïœä∑Ç≈Ç´ÇÈÅ@
//
//              posx
//  0    1    2 |  3    4
//  +----+----+----+----+


std::vector<float> HarmCoord3D::GetHarmonicCoordinate(const EVec3f &v)
{
  const int N = m_cage.m_vSize;
  const int W = m_map_reso[0];
  const int H = m_map_reso[1];
  const int D = m_map_reso[2], WH = W * H;

  const float posx = (v[0] - m_map_pos[0]) / m_map_pitch - 0.5f;
  const float posy = (v[1] - m_map_pos[1]) / m_map_pitch - 0.5f;
  const float posz = (v[2] - m_map_pos[2]) / m_map_pitch - 0.5f;
  const int   x    = (int) posx;
  const int   y    = (int) posy;
  const int   z    = (int) posz;
  const float xt = posx - x;
  const float yt = posy - y;
  const float zt = posz - z;

  float sum = 0;
  vector<float> hc(N, 0);

  //tri linearï‚ä‘Çé¿é{Ç∑ÇÈÅicageÇÕè\ï™ëÂÇ´Ç≠ÅCÇÕÇ›èoÇ∑éñÇÕñ≥Ç¢ëOíÒÇ≈é¿ëïÇ∑ÇÈÅj

  if ( x < 0 || y < 0 || z < 0 ||  W-1 <= x || H-1 <= y || D-1 <= z ) {
    std::cout << "cage is too small\n";
    return hc;
  }

  const int map_idx = x + y * W + z * W * H;
  
  for( int i=0; i < N; ++i)
  {
    float v0 = (1-yt) * ( (1-xt) * m_map_hc[i][map_idx     ] + xt * m_map_hc[i][map_idx+1     ] ) +
                 yt   * ( (1-xt) * m_map_hc[i][map_idx+W   ] + xt * m_map_hc[i][map_idx+1+W   ] );
    float v1 = (1-yt) * ( (1-xt) * m_map_hc[i][map_idx  +WH] + xt * m_map_hc[i][map_idx+1  +WH] ) + 
                 yt   * ( (1-xt) * m_map_hc[i][map_idx+W+WH] + xt * m_map_hc[i][map_idx+1+W+WH] );

    hc[i] = (1-zt) * v0 + zt * v1;
    sum += hc[i];
  }

  //std::cout << "sum of hc value (should be 1.0)" << sum << "\n";
  return hc;

}


HarmCoord3D::HarmCoord3D(const TMesh &cage, const int num_verts, const EVec3f *verts) : m_cage(cage)
{
  InitMap();
  CalcHarmonicCoordinateMap();

  m_verts_hc.clear();
  for ( int i=0; i < num_verts; ++i)
    m_verts_hc.push_back( GetHarmonicCoordinate(verts[i]) );


}

#pragma managed