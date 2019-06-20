#pragma once

#include "tmath.h"
#include <vector>


#pragma unmanaged

class HarmCoord2D
{ 
  //ƒJƒvƒZƒ‹‰»‚µ‚Ä‚à‚¢‚¢‚ñ‚¾‚¯‚Ç‚Æ‚è‚ ‚¦‚¸public‚Å
public: 

  //original cage 
  const std::vector<EVec2f> m_cage  ;
 
  //2D map resolution for computing harmonic coordinates (should be 2^x)
  const int m_map_reso ; // resolution of the map 
  float     m_map_pitch; // edge length of one pixel in map
  EVec2f    m_map_pos  ; // origin position of the map
  byte     *m_map_flag ; // 0:boundary, 1:inside, 2:outside
  std::vector<float*> m_map_hc; 



public:
  HarmCoord2D(const std::vector<EVec2f> &cage, int reso );

  std::vector<float> GetHarmonicCoodinateForPoint(const EVec2f p);


private:
  void InitMap();
  void InitHarmonicCoordinate();

};

#pragma managed
