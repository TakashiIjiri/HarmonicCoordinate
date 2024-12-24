#include "stdafx.h"
#include "HarmCoord3D.h"


#pragma unmanaged


using namespace std;

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
  
  // resolutionでは、切り上げで解像度を決定する
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
  EVec3f pitch(m_map_pitch, m_map_pitch, m_map_pitch);
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
//ただし d1 = (x1-x0)


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
  CalcBoundBox3D(x0,x1,x2, pitch * 2.5f, bbmin, bbmax); //modified


  //set boundary condition value
  const int WH = W*H;

  for ( int z = 0; z < D; ++z)
  {
    for ( int y = 0; y < H; ++y)
    {
      for ( int x = 0; x < W; ++x)
      {
        int i = x + y * W + z * WH;

        //consider boundary and background (bkg is necessary for target close to boundary)
        if ( map_flg[i] == 1 ) continue; 

        const EVec3f p( (x + 0.5f) * pitch + map_pos[0], 
                        (y + 0.5f) * pitch + map_pos[1], 
                        (z + 0.5f) * pitch + map_pos[2] );
      
        if ( !IsInBoundingBox(p, bbmin, bbmax) ) continue;

        //bacycentric coordinateを計算し三角形の中に入っていたら値 1-s-tを登録
        float s,t;
        if( CalcBarycentricCoordinate(p,x0,x1,x2,s,t) && 
            (x0 + s*(x1 - x0) + t*(x2 - x0) - p).norm() < 3 * pitch && 
            -0.01f <= s && s <= 1.01f && -0.01f <= t && t <= 1.01f &&
            s + t <= 1.0f) 
        {
          map_hc[i] = 1-s-t;
        }
      }
    }
  }
}

//
//+ todo bounding conditionの計算方法を改める
//+ todo 舌モデルでもちゃんと動くか確かめる --> roipainterに移動する
//+ で初期変形がちゃんとできるか確かめる



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

  //対象となる三角形を境界制約として描画
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
    const int W, 
    const int H, 
    const int D, 
    const byte *map_flg, // 0:out, 1:in, 2:boundary
    const int NUN_ITER , 
    float *map_hc        
)
{
  const int WH = W*H;
  float *tmp = new float[W*H*D];
  memcpy(tmp, map_hc, sizeof(float)*W*H*D);

  for ( int iter = 0; iter < NUN_ITER; ++iter)
  {
    for ( int z = 1; z < D - 1; ++z)
    {
      for ( int y = 1; y < H - 1; ++y)
      { 
        for ( int x = 1; x < W - 1; ++x)
        {
          int i = x + y * W + z * WH;
          if ( map_flg[i] != 1 ) continue;
          
          tmp[i] =  (map_hc[i- 1] + map_hc[i+ 1] + 
                     map_hc[i- W] + map_hc[i+ W] + 
                     map_hc[i-WH] + map_hc[i+WH]) / 6.0f;

        }
      }
    }
    memcpy(map_hc, tmp, sizeof(float)*W*H*D);
  }
}


static void LaplacianSmoothing_multilevel(
    const int  W,
    const int  H,
    const int  D,
    const int smoothing_iter_n,
    const int max_lv,
    const byte* map_flg, // 0:out, 1:in, 2:boundary
    float* map_hc        //contains values only at boundary
)
{
  std::vector<EVec3i> maps_reso;
  std::vector<byte* > maps_flag;
  std::vector<float*> maps_val;
  maps_reso.push_back(EVec3i(W, H, D));
  maps_flag.push_back((byte*)map_flg);
  maps_val .push_back(map_hc);

  //down sampling
  for (int lv = 1; lv < max_lv; ++lv)
  {
    const EVec3i p_reso = maps_reso[lv - 1];
    const byte*  p_flg  = maps_flag[lv - 1];
    const float* p_val  = maps_val [lv - 1];

    //calc resolution
    int pW = p_reso[0], pH = p_reso[1], pD = p_reso[2], pWH = pW * pH;
    int w = (pW + 1) / 2, h = (pH + 1) / 2, d = (pD + 1) / 2, wh = w * h;

    //gen flg map
    byte* flg = new byte[w * h * d];
    float* val = new float[w * h * d];
    byte* sum = new byte[w * h * d];
    memset(flg, 0, sizeof(byte) * w * h * d);
    memset(val, 0, sizeof(float) * w * h * d);
    memset(sum, 0, sizeof(byte) * w * h * d);
    maps_reso.push_back(EVec3i(w, h, d));
    maps_flag.push_back(flg);
    maps_val.push_back(val);

    //8個のchild cellsをまわる...
    //flag: ひとつでもboundary ならboundaryに
    //val : boundary cellの平均を入れる（in/out cellには0を入れる）

    //pre lv --> new lv
    for (int z = 0; z < pD; ++z)
      for (int y = 0; y < pH; ++y)
        for (int x = 0; x < pW; ++x)
        {
          int p_idx = x + y * pW + z * pWH;
          if (p_flg[p_idx] == 2) {
            int i = (x / 2) + (y / 2) * w + (z / 2) * wh;
            flg[i] = 2;
            sum[i] += 1;
            val[i] += p_val[p_idx];
          }
        }

    //boundary valueを平均値に
    for (int z = 0; z < d; ++z)
      for (int y = 0; y < h; ++y)
        for (int x = 0; x < w; ++x)
        {
          int i = x + y * w + z * w * h;
          if (flg[i] == 2) {
            val[i] /= (float)sum[i];
          }
          else {
            flg[i] = p_flg[2 * x + 2 * y * pW + 2 * z * pWH];
          }
        }

    delete[] sum;
  }

  // smooting and upsampling
  for (int lv = max_lv - 1; lv >= 0; --lv)
  {
    EVec3i reso = maps_reso[lv];
    byte* flg = maps_flag[lv];
    float* val = maps_val[lv];

    const int w = reso[0], h = reso[1], d = reso[2];

    //updampling copy value from coarse layer
    if (lv != max_lv - 1)
    {
      const int pW = maps_reso[lv + 1][0], pH = maps_reso[lv + 1][1], pD = maps_reso[lv + 1][2];
      const float* p_val = maps_val[lv + 1];

      for (int z = 0; z < d; ++z)
        for (int y = 0; y < h; ++y)
          for (int x = 0; x < w; ++x)
          {
            int i = x + y * w + z * w * h;
            if (flg[i] != 1) continue;
            val[i] = p_val[(x / 2) + (y / 2) * pW + (z / 2) * pW * pH];
          }
    }

    //smoothing 
    LaplacianSmoothing(reso[0], reso[1], reso[2], flg, smoothing_iter_n, val);
  }

  //creanup 
  for (int lv = 1; lv < max_lv; ++lv)
  {
    delete[] maps_val[lv];
    delete[] maps_flag[lv];
  }

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
  { 
    std::cout << "start LaplacianSmoothing_multilevel -- " << ci << "\n";
    //LaplacianSmoothing(W, H, D, m_map_flag, 5000, m_map_hc[ci]);
    LaplacianSmoothing_multilevel(W,H,D, 400, 5, m_map_flag, m_map_hc[ci]);
    std::cout << " -- " << ci << "DONE \n";
  }

}



// note : 画素の中心 (0.5)の位置に画素値があるとする
// 
// posx =  (v[0] - m_map_pos[0]) / m_map_pitch で、画素格子上における 画素位置が得られる
// posx = posx-0.5 で下図のように，ちょうど0,1,2に画素値がある空間における位置に変換できる　
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

  //tri linear補間を実施する（cageは十分大きく，はみ出す事は無い前提で実装する）

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

  ////debug 
  //int f1 = m_map_flag[map_idx];
  //int f2 = m_map_flag[map_idx + 1];
  //int f3 = m_map_flag[map_idx + W];
  //int f4 = m_map_flag[map_idx + 1 + W];
  //int f5 = m_map_flag[map_idx + WH];
  //int f6 = m_map_flag[map_idx + 1 + WH];
  //int f7 = m_map_flag[map_idx + W + WH];
  //int f8 = m_map_flag[map_idx + 1 + W + WH];

  //bool tf = 
  //     f1 == 0 || f2 == 0 || f3 == 0 || f4 == 0 || f5 == 0 || f6 == 0 || f7 == 0 || f8 == 0 ||
  //     f1 == 3 || f2 == 3 || f3 == 3 || f4 == 3 || f5 == 3 || f6 == 3 || f7 == 3 || f8 == 3;

  ////debug 
  //std::cout << "sum of hc value (should be 1.0)" << sum ;
  //std::cout << (tf?" outside " : " ") << f1 << f2 << f3 << f4 << f5 << f6 << f7 << f8 <<"\n";

  return hc;

}


HarmCoord3D::HarmCoord3D(const TMesh &cage, const int num_verts, const EVec3f *verts) : m_cage(cage)
{ 
  //initialzie map
  //m_map_pos
  //m_map_reso
  //m_map_flag :0 outside, 1 inside, 2: boundary
  InitMap();
  
  //compute harmonic coodinate map (volume)
  //m_map_hc : harmonic coordianate for each vertex
  CalcHarmonicCoordinateMap();

  //compute harmonic coordinate for each vertex
  m_verts_hc.clear();
  for ( int i = 0; i < num_verts; ++i)
  {
    m_verts_hc.push_back(GetHarmonicCoordinate(verts[i]));
  }
}

#pragma managed