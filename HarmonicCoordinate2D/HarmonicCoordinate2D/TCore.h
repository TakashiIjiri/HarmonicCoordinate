#pragma once


#include "tmath.h"
#include "OglForCLI.h"
#include "HarmCoord2D.h"
#include <vector>

#pragma unmanaged

class TCore
{

  std::vector<EVec2f> m_cage ; //cage
  std::vector<EVec2f> m_verts; //“à•”‚Ì’¸“_

  std::vector<std::vector<float>> m_verts_in_hc; //vertices in harmonic coordinates

  HarmCoord2D *m_harmcoord;

  bool m_vis_hc_weight;
  bool m_vis_hc_flag   ;
  int  m_vis_hc_idx;

  //mouse manipulation
  bool m_bR, m_bL, m_bM;
  int  m_drag_cagept_id; // -1 when false
  std::vector<EVec2f> m_stroke; //“à•”‚Ì’¸“_
  


private:
  TCore();
  ~TCore();
public:
  static TCore* getInst(){ 
    static TCore p;
    return &p;
  }

  void LBtnDown(EVec2i p, OglForCLI *ogl);
  void RBtnDown(EVec2i p, OglForCLI *ogl);
  void MBtnDown(EVec2i p, OglForCLI *ogl);
  void LBtnUp(EVec2i p, OglForCLI *ogl);
  void RBtnUp(EVec2i p, OglForCLI *ogl);
  void MBtnUp(EVec2i p, OglForCLI *ogl);
  void MouseMove(EVec2i p, OglForCLI *ogl);

  void DrawScene();

  void KeyDown(int k);

private:
  void DrawHarmonicCoordinateFlag  ();
  void DrawHarmonicCoordinateWeight();

  void UpdateVertsByHarmonicCoord();

  int  PickCagePoint(EVec2f p);
};

inline bool IsCtrKeyOn  (){ return GetKeyState( VK_CONTROL ) < 0 ; }
inline bool IsSpaceKeyOn(){ return GetKeyState( VK_SPACE   ) < 0 ; }
inline bool IsShiftKeyOn(){ return GetKeyState( VK_SHIFT   ) < 0 ; }
inline bool IsAltKeyOn  (){ return GetKeyState( VK_MENU    ) < 0 ; }


#pragma managed
