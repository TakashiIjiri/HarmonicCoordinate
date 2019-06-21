#pragma once

#include "COMMON/tmesh.h"

class HarmCoord3D
{
public:
  //original cage 
  const TMesh m_cage  ;
 
  // -- 3D map configuration -- 
  //等方性ボクセル，解像度（辺の長さは） 2のべき，
  //少なくとも周囲1ボクセル分は背景に
  //cage vertexごとの計算結果を保持するのは現実的ではないので頂点ごとに計算したあと破棄する
  EVec3i  m_map_reso ; // 3d map resolution 
  float   m_map_pitch; // cubic voxel pitch
  EVec3f  m_map_pos  ; // origin position of the map
  byte   *m_map_flag ; // 0:boundary, 1:inside, 2:outside

  //harmonic coordinate mapを保持する実装（メモリ的にはしんどいか）
  std::vector<float*> m_map_hc; //harmonic coordinate map for each cage vertex

  //harmonic coordinate for internal vertices
  std::vector<std::vector<float>> m_verts_hc;

public:

  HarmCoord3D(const TMesh &cage, const int num_verts, const EVec3f *verts);


private:
  void InitMap();
  void CalcHarmonicCoordinateMap();
  std::vector<float> GetHarmonicCoordinate(const EVec3f &v);
};

