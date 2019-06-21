#pragma once

#include "COMMON/tmesh.h"

class HarmCoord3D
{
public:
  //original cage 
  const TMesh m_cage  ;
 
  // -- 3D map configuration -- 
  //�������{�N�Z���C�𑜓x�i�ӂ̒����́j 2�ׂ̂��C
  //���Ȃ��Ƃ�����1�{�N�Z�����͔w�i��
  //cage vertex���Ƃ̌v�Z���ʂ�ێ�����̂͌����I�ł͂Ȃ��̂Œ��_���ƂɌv�Z�������Ɣj������
  EVec3i  m_map_reso ; // 3d map resolution 
  float   m_map_pitch; // cubic voxel pitch
  EVec3f  m_map_pos  ; // origin position of the map
  byte   *m_map_flag ; // 0:boundary, 1:inside, 2:outside

  //harmonic coordinate map��ێ���������i�������I�ɂ͂���ǂ����j
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

