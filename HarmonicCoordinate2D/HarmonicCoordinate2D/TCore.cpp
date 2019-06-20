#include "stdafx.h"
#include <iostream>
#include "TCore.h"
#include "MainForm.h"
#include <vector>

using namespace std;


TCore::~TCore()
{
  delete m_harmcoord;
}


TCore::TCore()
{  
  //generate cage model
  m_cage = {EVec2f( 70,200),EVec2f(76,168), EVec2f(50,160), EVec2f(17,106), EVec2f(40,90), 
            EVec2f( 75,140),EVec2f(66,89 ), EVec2f(55,48 ), EVec2f(49,  9), EVec2f(84,9  ), EVec2f(87,43),
            EVec2f(100,82),
            EVec2f(113,43), EVec2f(116,9  ), EVec2f(151, 9), EVec2f(145,48 ),EVec2f(134,89 ), EVec2f( 125,140), 
            EVec2f(160,90), EVec2f(183,106), EVec2f(150,160),EVec2f(124,168), EVec2f(130,200)};

  //normalize the model: note screen size is [-10,10]x[-10,10]
  float TRGT_SIZE = 10.f + 10.f - 2.f;

  EVec2f bb_min, bb_max;
  t_calcBoundBox2D( m_cage, bb_min, bb_max);
  float scale = std::max(bb_max[0]-bb_min[0], bb_max[1]-bb_min[1]);

  for ( auto& c : m_cage) {
    c -= bb_min;
    c /= scale;
    c *= TRGT_SIZE;
    c[0] -= TRGT_SIZE/2;
    c[1] -= TRGT_SIZE/2;
    std::cout << c[0] << " " << c[1];
  }
  
  //generate vertex in the cage!!
  m_harmcoord = new HarmCoord2D( m_cage, 256 );

  //init mouse state
  m_bR = m_bL = m_bM = false;
  m_drag_cagept_id = false;

  //init vis flg
  m_vis_hc_weight = true;
  m_vis_hc_flag   = false;
  m_vis_hc_idx = 0;

}




//Mouse listeners

void TCore::LBtnDown(EVec2i p, OglForCLI *ogl)
{
  m_bL = true;

  EVec3f rayp, rayd;
  ogl->GetCursorRay(p, rayp, rayd);
  EVec2f rayp2f(rayp[0], rayp[1]);
  m_drag_cagept_id = PickCagePoint(rayp2f);
}

void TCore::RBtnDown(EVec2i p, OglForCLI *ogl)
{
  m_bR = true;
}

void TCore::MBtnDown(EVec2i p, OglForCLI *ogl)
{
  m_bM = true;
}

void TCore::LBtnUp(EVec2i p, OglForCLI *ogl)
{
  //update m_verts_in_hc
  for (const auto p : m_stroke)
  {
    auto hc = m_harmcoord->GetHarmonicCoodinateForPoint(p);
    m_verts.push_back(p);
    m_verts_in_hc.push_back(hc);
  }

  UpdateVertsByHarmonicCoord();
  m_stroke.clear();
  m_drag_cagept_id = -1;
  m_bL = false;
}


void TCore::RBtnUp(EVec2i p, OglForCLI *ogl)
{
  m_bR = false;
}
void TCore::MBtnUp(EVec2i p, OglForCLI *ogl)
{
  m_bM = false;
}


void TCore::MouseMove(EVec2i p, OglForCLI *ogl)
{
  if( !m_bR && !m_bM && !m_bL) return ;

  EVec3f rayPos, rayDir;
  ogl->GetCursorRay(p, rayPos, rayDir);

  if( m_bL && IsShiftKeyOn() )
  {
    EVec2f p( rayPos[0], rayPos[1] );
    if ( m_stroke.size() ==0 || (m_stroke.back() - p).norm() > 0.7 )
      m_stroke.push_back(p);
  }

  if ( m_drag_cagept_id >= 0 )
  {
    m_cage[m_drag_cagept_id] << rayPos[0], rayPos[1];
    UpdateVertsByHarmonicCoord();
  }

  HarmonicCoordinate2D::MainForm_RedrawMainPanel();
}


//key listener
void TCore::KeyDown(int k)
{
  if ( k == 37 ) //left
  {
    --m_vis_hc_idx;
    if( m_vis_hc_idx < 0  ) m_vis_hc_idx = m_cage.size() - 1;
  }

  if ( k == 39 ) //right
  {
    ++m_vis_hc_idx;
    if( m_vis_hc_idx >= m_cage.size() ) m_vis_hc_idx = 0;
  }

  if ( k == 49 ) //1 key
  {
    m_vis_hc_flag = !m_vis_hc_flag; 
  }
  if ( k == 50 ) //2 key
  {
    m_vis_hc_weight = !m_vis_hc_weight; 
  }
  if ( k == 51 ) //3 key
  {
    //init cage 
    m_cage = m_harmcoord->m_cage;
    UpdateVertsByHarmonicCoord();
  }

  HarmonicCoordinate2D::MainForm_RedrawMainPanel();
}




void TCore::DrawScene(){

  glDisable(GL_LIGHTING);
  glLineWidth(3);

  if( m_vis_hc_flag  ) DrawHarmonicCoordinateFlag();
  if( m_vis_hc_weight) DrawHarmonicCoordinateWeight();

  glColor3d(1,1,0);
  glBegin(GL_LINE_STRIP);
  for ( const auto& c : m_cage) glVertex2fv(c.data());
  if ( m_cage.size() > 0 ) glVertex2fv(m_cage.front().data() );
  glEnd();

  // render cage points
  glPointSize(14);
  glColor3d(1,0,0);
  glBegin(GL_POINTS);
  for ( const auto& c : m_cage) glVertex2fv(c.data());
  glEnd();

  // render internal points 
  glPointSize(4);

  glBegin(GL_POINTS);
    glColor3d(1,0,1);
    for ( const auto& c : m_stroke) glVertex2fv(c.data());
    glColor3d(0,1,1);
    for ( const auto& c : m_verts) glVertex2fv(c.data());
  glEnd();





}


void TCore::DrawHarmonicCoordinateFlag  ()
{
  const int    reso  = m_harmcoord->m_map_reso ;
  const float  pitch = m_harmcoord->m_map_pitch;
  const EVec2f posi  = m_harmcoord->m_map_pos  ;
  const byte   *flg  = m_harmcoord->m_map_flag;

  glPointSize(4);
  glDisable(GL_LIGHTING);
  glBegin(GL_POINTS);
  for ( int y = 0; y < reso; ++y)
  {
    for ( int x = 0; x < reso; ++x)
    {
      int i= x + y*reso;
      if (     flg[i] == 0 ) glColor3f(0,0,1);
      else if (flg[i] == 1 ) glColor3f(1,0,0);
      else glColor3f(0,1,0);    

      glVertex2f((x+0.5f)*pitch + posi[0], (y+0.5f)*pitch + posi[1]);
    }
  }
  glEnd();


}



void TCore::DrawHarmonicCoordinateWeight()
{
  const int    reso  = m_harmcoord->m_map_reso ;
  const float  pitch = m_harmcoord->m_map_pitch;
  const EVec2f posi  = m_harmcoord->m_map_pos  ;
  const float *hc    = m_harmcoord->m_map_hc[m_vis_hc_idx];

  glPointSize(4);
  glDisable(GL_LIGHTING);
  glBegin(GL_POINTS);
  for ( int y = 0; y < reso; ++y)
  {
    for ( int x = 0; x < reso; ++x)
    {
      int i= x + y*reso;
      float c = hc[i] * 2;
      glColor3f(c,c,0);    
      glVertex2f((x+0.5f)*pitch + posi[0], (y+0.5f)*pitch + posi[1]);
    }
  }
  glEnd();


}




void TCore::UpdateVertsByHarmonicCoord()
{
  const int cN = (int) m_cage.size();

  for( int i = 0; i < (int) m_verts.size(); ++i)
  {
    EVec2f p(0,0);
    for ( int j = 0; j < cN; ++j)
      p += m_verts_in_hc[i][j] * m_cage[j];
    m_verts[i] = p;
  }

}



int  TCore::PickCagePoint(EVec2f p)
{
  const int N = (int) m_cage.size();

  for ( int i=0; i < N; ++i )
    if ( t_dist(m_cage[i], p) < 0.2) //閾値適当。。。
      return i;

  return -1;
}
