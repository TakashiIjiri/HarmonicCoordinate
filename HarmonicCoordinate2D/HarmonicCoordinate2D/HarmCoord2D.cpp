#include "stdafx.h"
#include "HarmCoord2D.h"


#include <map>

using namespace std;


#pragma unmanaged




// cast ray in Y axis   
static void PaintInCageByRayCasting
(
  const EVec2i &map_reso ,
  const float  &map_pitch,
  const EVec2f &map_pos  ,
  const std::vector<EVec2f> &cageverts, 

  byte *map_flg //allocated[WxH], 0:out, 1:in
)
{
  const int N = (int)cageverts.size();
  const int W = map_reso[0];
  const int H = map_reso[1];

  memset( map_flg, 0, sizeof(byte)*W*H);

  // ray casting along x axis to fill inside the mesh 
#pragma omp parallel for
  for (int xI = 0; xI < W; ++xI) 
  {
    //cast ray along y axis (x=this value and direction is (0,1,0) )
    
    float ray_x = (0.5f + xI) * map_pitch + map_pos[0];
    std::multimap<double, double> blist;// (yPos, normInYir);

    for ( int i=0; i < N; ++i) 
    {
      const EVec2f &c0 = cageverts[i      ];
      const EVec2f &c1 = cageverts[(i+1)%N];

      //calc intersect line secment (c0,c1) and ray (x = x)
      bool isIntersect = (c0[0] <= ray_x && ray_x <= c1[0]) || (c0[0] >= ray_x && ray_x >= c1[0]);
      if( c0[0] == c1[0] || !isIntersect ) continue;

      float t = (ray_x - c0[0]) / (c1[0]-c0[0]); 
      float y = c0[1] + (c1[1]-c0[1]) * t;
      float normalx = -(c1[1]-c0[1]); //not used 
      float normaly =   c1[0]-c0[0];
      blist.insert( std::make_pair(y, normaly)); //(y座標, normal in y)
    }

    //clean up blist (頂点上で起こった交差重複を削除)
    while (blist.size() != 0)
    {
      if (blist.size() == 1) { blist.clear(); break; }

      bool found = false;
      auto it0 = blist.begin();
      auto it1 = blist.begin(); it1++;

      for (; it1 != blist.end(); ++it0, ++it1) if (it0->second * it1->second > 0)
      {
        blist.erase(it1);
        found = true;
        break;
      }
      if (!found) break;
    }

    bool flag = false;
    int yI = 0;

    //int pivIdx = xI ;
    for (auto it = blist.begin(); it != blist.end(); ++it)
    {
      int piv_yi = (int)( (it->first - map_pos[1]) / map_pitch);
      for (; yI <= piv_yi && yI < H; ++yI) map_flg[ xI + yI * W] = flag ? 1 : 0;
      flag = !flag;
    }
    if (flag == true) std::cout << "error double check here!";
  }
}









//note : map[i,j] corresponds to [ m_map_pos + m_map_pitch * (i+0.5,j+0.5) ]
//note : (x,y) in 2d space corresponds to map[i,j] 
//              (i.j) = ((x,y) - m_map_pos) / m_map_pitch - 0.5; 
void HarmCoord2D::InitMap()
{
  //calc map configuration
  EVec2f bb_min, bb_max;
  t_calcBoundBox2D( m_cage, bb_min, bb_max);

  m_map_pitch = std::max(bb_max[0]-bb_min[0], bb_max[1]-bb_min[1] ) / (m_map_reso - 2);
  m_map_pos   << bb_min[0] - m_map_pitch, bb_min[1] - m_map_pitch; 

  //calc flag image (0:out,1:inside,2:boundary)
  m_map_flag = new byte[m_map_reso*m_map_reso];
  PaintInCageByRayCasting( EVec2i(m_map_reso,m_map_reso), m_map_pitch, m_map_pos, m_cage, m_map_flag);

  for ( int y = 1; y < m_map_reso-1; ++y)
  {
    for ( int x = 1; x < m_map_reso-1; ++x)
    {
      int i= x + y * m_map_reso;
      if ( m_map_flag[i] != 1 ) continue;

      if( m_map_flag[i-1] == 0 || m_map_flag[i+1] == 0 || 
          m_map_flag[i-m_map_reso] == 0 || m_map_flag[i+m_map_reso] == 0)
      {
        m_map_flag[i] = 2; // set as boundary 
      }
    }
  }

}



static void CalcBoundingBox(const EVec2f &p0, const EVec2f &p1, EVec2f &bbmin, EVec2f &bbmax)
{
  bbmin << std::min(p0[0], p1[0]), std::min(p0[1], p1[1]);

  bbmax << std::max(p0[0], p1[0]), std::max(p0[1], p1[1]);
}

static bool IsInBoundingBox(const EVec2f &p, const EVec2f &bbmin, const EVec2f &bbmax)
{
  return bbmin[0] <= p[0]     && bbmin[1] <= p[1] && 
             p[0] <= bbmax[0] &&     p[1] <= bbmax[1];
}


// cagepoint  --> 1
// cagepoint1 --> 0
// cagepoint2 --> 0
static void SetBoundaryCondition
(
  const int    &map_reso ,
  const float  &map_pitch,
  const EVec2f &map_pos  ,
  const byte   *map_flg, // 0:out, 1:in, 2:boundary

  const EVec2f &c0, //cage point (pivot) 
  const EVec2f &c1, //cage point (neighbor)
  const EVec2f &c2, //cage point (neighbor)
  
  float *map_hc
)
{
  //init map_hc
  memset( map_hc, 0, sizeof(float) * map_reso * map_reso );

  EVec2f bb1_min, bb1_max, bb2_min, bb2_max;
  CalcBoundingBox(c0, c1, bb1_min, bb1_max);
  CalcBoundingBox(c0, c2, bb2_min, bb2_max);
  bb1_min[0] -= 2.5f * map_pitch; bb1_min[1] -= 2.5f * map_pitch;
  bb2_min[0] -= 2.5f * map_pitch; bb2_min[1] -= 2.5f * map_pitch;
  bb1_max[0] += 2.5f * map_pitch; bb1_max[1] += 2.5f * map_pitch;
  bb2_max[0] += 2.5f * map_pitch; bb2_max[1] += 2.5f * map_pitch;


  //set boundary condition value
  int counter = 0;
  for ( int y = 0; y < map_reso; ++y)
  {
    for ( int x = 0; x < map_reso; ++x)
    {
      int i = x + y * map_reso;
      if ( map_flg[i] != 2 ) continue;

      const EVec2f p( (x+0.5f)*map_pitch + map_pos[0], 
                      (y+0.5f)*map_pitch + map_pos[1]);
      
      //check  cagepoint0-1
      if ( IsInBoundingBox(p, bb1_min, bb1_max) )
      {
        float t1;
        if( t_distPointToLineSegment( p, c0, c1,t1) < map_pitch * 2.0 ) {
          map_hc[i] = 1.0f - t1;
          counter++;
        }
      }

      //check  cagepoint0-2
      if ( IsInBoundingBox(p, bb2_min, bb2_max) )
      {
        float t2;
        if( t_distPointToLineSegment( p, c0, c2, t2) < map_pitch * 2.0) {
          map_hc[i] = 1.0f - t2;
          counter++;
        }
      }
    }
  }

  std::cout << " NUM OF BOUNDARY VTX " << counter << "\n";
}

static void LaplacianSmoothing(
  const int    &map_reso ,
  const byte   *map_flg, // 0:out, 1:in, 2:boundary
  float *map_hc //contains values only at boundary
)
{
  //todo here
  float *tmp = new float[map_reso*map_reso];
  memcpy(tmp, map_hc, sizeof(float)*map_reso*map_reso);

  for ( int iter = 0; iter < 12000; ++iter)
  {
    float diff = 0;
    for ( int y = 1; y < map_reso-1; ++y)
    {
      for ( int x = 1; x < map_reso-1; ++x)
      {
        int i = x + y * map_reso;
        if ( map_flg[i] != 1 ) continue;
        tmp[i] = 0.25f * (map_hc[i-1] + map_hc[i+1] + map_hc[i-map_reso] + map_hc[i+map_reso]);

        diff += (tmp[i] - map_hc[i])*(tmp[i] - map_hc[i]);
      }
    }
    memcpy(map_hc, tmp, sizeof(float)*map_reso*map_reso);
    //std::cout << iter << " " << diff << "\n";
  }
}




void HarmCoord2D::InitHarmonicCoordinate()
{
  const int N = (int)m_cage.size();

  for( int ci = 0; ci < N; ++ci)
  {
    float *map_hc = new float[m_map_reso*m_map_reso];
    const EVec2f &c0 = m_cage[ci        ];
    const EVec2f &c1 = m_cage[(ci+1  )%N];
    const EVec2f &c2 = m_cage[(ci-1+N)%N];

    SetBoundaryCondition(m_map_reso, m_map_pitch, m_map_pos, m_map_flag, c0,c1,c2, map_hc);
    m_map_hc.push_back(map_hc);
  }

#pragma omp parallel for
  for( int ci = 0; ci < N; ++ci)
    LaplacianSmoothing(m_map_reso, m_map_flag, m_map_hc[ci]);
}



HarmCoord2D::HarmCoord2D( const std::vector<EVec2f> &cage, int reso = 256 ) :
  m_map_reso(reso), m_cage(cage)
{
  //compute m_flg_map
  InitMap();

  //compute harmonic coordinate for all cage points
  InitHarmonicCoordinate();

}


std::vector<float> HarmCoord2D::GetHarmonicCoodinateForPoint(const EVec2f p)
{
  const int N = (int) m_cage.size();
  int x = (int)( (p[0] - m_map_pos[0]) / m_map_pitch );
  int y = (int)( (p[1] - m_map_pos[1]) / m_map_pitch );

  float sum = 0;
  vector<float> hc(N);
  for( int i=0; i < N; ++i)
  {
    hc[i] = m_map_hc[i][x + y * m_map_reso];
    sum += hc[i];
  }

  std::cout << "aaaaa  " << sum << "\n";
  return hc;
}




#pragma managed
