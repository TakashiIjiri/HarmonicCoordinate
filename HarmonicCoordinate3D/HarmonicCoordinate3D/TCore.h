#pragma once

#include "COMMON/tmath.h"
#include "COMMON/tmesh.h"
#include "COMMON/OglForCLI.h"
#include "HarmCoord3D.h"



#pragma unmanaged


class TCore
{
private:
  //mouse listener 
  bool m_bL, m_bR, m_bM;
  int  m_drag_cage_id  ; //-1 when false
  TMesh m_cp_mesh;


  //vis flg
 
  int m_vis_hc_idx;
  bool m_vis_harmcoord_flg;
  bool m_vis_harmcoord_value;

  TMesh m_mesh, m_cage;

  HarmCoord3D *m_harmcoord;


  TCore();
public:

  static TCore* GetInst(){
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
  int PickCageVertex(EVec3f ray_pos, EVec3f ray_dir );
  void DrawHarmCoordFlag();
  void DrawHarmCoordValue();

  void UpdateShapeByHC();
};



inline bool IsCtrKeyOn  (){ return GetKeyState( VK_CONTROL ) < 0 ; }
inline bool IsSpaceKeyOn(){ return GetKeyState( VK_SPACE   ) < 0 ; }
inline bool IsShiftKeyOn(){ return GetKeyState( VK_SHIFT   ) < 0 ; }
inline bool IsAltKeyOn  (){ return GetKeyState( VK_MENU    ) < 0 ; }


#pragma managed
