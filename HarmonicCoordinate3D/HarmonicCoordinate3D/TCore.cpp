#include "stdafx.h"
#include "TCore.h"
#include "Mainform.h"



#pragma unmanaged

//cage
//     6 7 8  14 15 16   23 24 25
// ^   3 4 5  12    13   20 21 22
// |   0 1 2  09 10 11   17 18 19
// +-->x
//
//

using namespace HarmonicCoordinate3D;
using namespace std;

static float CP_RADIUS = 0.15f;


TCore::TCore()
{
  //load 3d mesh
  
  if (true)
  { 
    m_mesh.initialize("bunny.obj");
    EVec3f bbmin, bbmax;
    m_mesh.getBoundBox(bbmin, bbmax);
    float length = max3(bbmax[0] - bbmin[0], bbmax[1] - bbmin[1], bbmax[2] - bbmin[2]);
    m_mesh.Translate( -0.5 * (bbmin + bbmax) );
    m_mesh.Scale(4.0f / length);

    m_cage.initialize("bunnycage.obj");
  }
  else
  {
    //for different obj and cage
    m_mesh.initialize("bunny.obj");
    m_cage.initialize("bunnycage.obj");

    EVec3f cuboid(10, 10, 10);
    EVec3f bbmin, bbmax;
    m_cage.getBoundBox(bbmin, bbmax);
    int idx = max3id(bbmax[0] - bbmin[0], bbmax[1] - bbmin[1], bbmax[2] - bbmin[2]);
    m_cage.Scale(0.5f * cuboid[idx] / (bbmax[idx] - bbmin[idx]));
    m_mesh.Scale(0.5f * cuboid[idx] / (bbmax[idx] - bbmin[idx]));
    EVec3f gc = m_cage.getGravityCenter();
    m_cage.Translate(-gc + cuboid / 2);
    m_mesh.Translate(-gc + cuboid / 2);
  }


  m_cp_mesh.initializeIcosaHedron( CP_RADIUS );


  m_vis_harmcoord_flg   = true;
  m_vis_harmcoord_value = true;
  m_vis_hc_idx = 0;
  m_drag_cage_id = -1;


  m_harmcoord = new HarmCoord3D(m_cage, m_mesh.m_vSize, m_mesh.m_vVerts);

}


/* 
  //gen cage
  //直方体ケージ生成（最初は使ってた）
  m_mesh.getBoundBox(bbmin, bbmax);
  float x0 = bbmin[0]-0.5f, x2 = bbmax[0] + 0.2f, x1 = 0.5f * (x0 + x2); 
  float y0 = bbmin[1]-0.5f, y2 = bbmax[1] + 0.2f, y1 = 0.5f * (y0 + y2); 
  float z0 = bbmin[2]-0.5f, z2 = bbmax[2] + 0.2f, z1 = 0.5f * (z0 + z2); 

  vector<EVec3f> Vs = {
    EVec3f(x0,y0,z0), EVec3f(x1,y0,z0), EVec3f(x2,y0,z0), 
    EVec3f(x0,y1,z0), EVec3f(x1,y1,z0), EVec3f(x2,y1,z0), 
    EVec3f(x0,y2,z0), EVec3f(x1,y2,z0), EVec3f(x2,y2,z0), 
    EVec3f(x0,y0,z1), EVec3f(x1,y0,z1), EVec3f(x2,y0,z1), 
    EVec3f(x0,y1,z1),                   EVec3f(x2,y1,z1), 
    EVec3f(x0,y2,z1), EVec3f(x1,y2,z1), EVec3f(x2,y2,z1), 
    EVec3f(x0,y0,z2), EVec3f(x1,y0,z2), EVec3f(x2,y0,z2), 
    EVec3f(x0,y1,z2), EVec3f(x1,y1,z2), EVec3f(x2,y1,z2), 
    EVec3f(x0,y2,z2), EVec3f(x1,y2,z2), EVec3f(x2,y2,z2)};
  vector<TPoly > Ps = {
    TPoly(0,3,4), TPoly(0,4,1), TPoly(1,4,5), TPoly(1,5,2), 
    TPoly(3,6,7), TPoly(3,7,4), TPoly(4,7,8), TPoly(4,8,5), 

    TPoly(0,10, 9), TPoly(0,1,10), TPoly(1,11,10), TPoly(1,2,11), 
    TPoly(2,13,11), TPoly(2,5,13), TPoly(5,16,13), TPoly(5,8,16), 
    TPoly(8,15,16), TPoly(8,7,15), TPoly(7,14,15), TPoly(7,6,14), 
    TPoly(6,12,14), TPoly(6,3,12), TPoly(3, 9,12), TPoly(3,0, 9),

    TPoly( 9,18,17), TPoly( 9,10,18), TPoly(10,19,18), TPoly(10,11,19), 
    TPoly(11,22,19), TPoly(11,13,22), TPoly(13,25,22), TPoly(13,16,25), 
    TPoly(16,24,25), TPoly(16,15,24), TPoly(15,23,24), TPoly(15,14,23), 
    TPoly(14,20,23), TPoly(14,12,20), TPoly(12,17,20), TPoly(12,9,17),

    TPoly(17,21,20), TPoly(17,18,21), TPoly(18,22,21), TPoly(18,19,22), 
    TPoly(20,24,23), TPoly(20,21,24), TPoly(21,25,24), TPoly(21,22,25) 
  };
  
  m_cage.initialize(Vs,Ps);

*/




void TCore::LBtnDown(EVec2i p, OglForCLI *ogl)
{
  m_bL = true;
  if( IsShiftKeyOn() )
  {
    EVec3f ray_pos, ray_dir;
    ogl->GetCursorRay(p, ray_pos, ray_dir);
    m_drag_cage_id = PickCageVertex(ray_pos, ray_dir);
  }
  else
  {
    ogl->BtnDown_Trans(p);
  }
}

void TCore::RBtnDown(EVec2i p, OglForCLI *ogl)
{
  m_bR = true;
  ogl->BtnDown_Rot(p);
}

void TCore::MBtnDown(EVec2i p, OglForCLI *ogl)
{
  m_bM = true;
  ogl->BtnDown_Zoom(p);
}

void TCore::LBtnUp(EVec2i p, OglForCLI *ogl)
{
  if( m_drag_cage_id != -1) UpdateShapeByHC();

  m_bL = false;
  m_drag_cage_id = -1;
  ogl->BtnUp();
}

void TCore::RBtnUp(EVec2i p, OglForCLI *ogl)
{
  m_bR = false;
  ogl->BtnUp();
}

void TCore::MBtnUp(EVec2i p, OglForCLI *ogl)
{
  m_bM = false;
  ogl->BtnUp();
}

void TCore::MouseMove(EVec2i p, OglForCLI *ogl)
{
  if ( m_bL && m_bR && m_bM ) return;

  if( m_drag_cage_id != -1)
  {
    EVec3f ray_pos, ray_dir;
    ogl->GetCursorRay(p, ray_pos, ray_dir);

    float d = (m_cage.m_vVerts[m_drag_cage_id] - ray_pos).norm();
    m_cage.m_vVerts[m_drag_cage_id] = ray_pos + d * ray_dir;
    UpdateShapeByHC();


  }
  else
  {
    ogl->MouseMove(p);
  }

  MainForm_RedrawMainPanel();
}



void TCore::KeyDown(int k)
{
  std::cout << k;
  if( k == 83 ) m_cage.exportObjNoTexCd("test.obj");
  if( k == 49 ) m_vis_harmcoord_flg   = !m_vis_harmcoord_flg  ;
  if( k == 50 ) m_vis_harmcoord_value = !m_vis_harmcoord_value;

   if ( k == 37 ) //left
  {
    --m_vis_hc_idx;
    if( m_vis_hc_idx < 0  ) m_vis_hc_idx = m_cage.m_vSize - 1;
  }

  if ( k == 39 ) //right
  {
    ++m_vis_hc_idx;
    if( m_vis_hc_idx >= m_cage.m_vSize ) m_vis_hc_idx = 0;
  }
  MainForm_RedrawMainPanel();

}


static float diff[4] = {0.5f, 0.5f, 0.2f, 0.5f};
static float ambi[4] = {0.5f, 0.5f, 0.2f, 0.5f};
static float diffR[4] ={1.0f, 0.2f, 0.2f, 0.5f};
static float ambiR[4] ={1.0f, 0.2f, 0.2f, 0.5f};
static float spec[4] = {1.0f, 1.0f, 1.0f, 0.5f};
static float shin[1] = {64};


void TCore::DrawScene()
{
  glEnable(GL_DEPTH_TEST);
  glDisable(GL_BLEND);
  //draw world coordinate frame
  glDisable(GL_LIGHTING );

  glLineWidth(3);
  glBegin(GL_LINES);
  glColor3d(1,0,0); glVertex3d(0,0,0); glVertex3d(2,0,0);
  glColor3d(0,1,0); glVertex3d(0,0,0); glVertex3d(0,2,0);
  glColor3d(0,0,1); glVertex3d(0,0,0); glVertex3d(0,0,2);
  glEnd();


  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHT1);
  glEnable(GL_LIGHT2);
  if (!IsSpaceKeyOn())
    m_mesh.draw(diff, ambi, spec, shin);
  //m_cage.draw(diff, ambi, spec, shin);

  for( int i=0; i < m_cage.m_vSize; ++i) {
    glTranslatef( m_cage.m_vVerts[i][0], m_cage.m_vVerts[i][1], m_cage.m_vVerts[i][2]);
    m_cp_mesh.draw(diffR,ambiR,spec,shin);
    glTranslatef(-m_cage.m_vVerts[i][0],-m_cage.m_vVerts[i][1],-m_cage.m_vVerts[i][2]);
    
  }

  glDisable(GL_LIGHTING);
  m_cage.drawEdges(2,0.5,0.5,0.5);


  if( m_vis_harmcoord_flg   ) DrawHarmCoordFlag ();
  if( m_vis_harmcoord_value ) DrawHarmCoordValue();


}




void TCore::DrawHarmCoordValue()
{
  const int   W = m_harmcoord->m_map_reso[0];
  const int   H = m_harmcoord->m_map_reso[1];
  const int   D = m_harmcoord->m_map_reso[2], WH = W*H;
  const float p = m_harmcoord->m_map_pitch  ;
  const EVec3f orig = m_harmcoord->m_map_pos;
  const byte* flg = m_harmcoord->m_map_flag ;
  const float* hc = m_harmcoord->m_map_hc[m_vis_hc_idx];

  glDisable(GL_LIGHTING );
  glPointSize(3);
  glBegin( GL_POINTS);
  for ( int z = 0; z < D; ++z) 
  {
    for ( int y = 0; y < H; ++y) 
    {
      for ( int x = 0; x < W; ++x) 
      {
        const int i = x + y * W + z * WH;
        if (flg[i] == 2 || (flg[i] == 0 && hc[i] > 0))
        {
          float xx = orig[0] + (x + 0.5f) * p;
          float yy = orig[1] + (y + 0.5f) * p;
          float zz = orig[2] + (z + 0.5f) * p;
          float c=  hc[i] * 1.5f;
          glColor3f(c, c, 0);
          glVertex3f(xx,yy,zz);
        }
      }
    }
  }
  glEnd();

}

void TCore::DrawHarmCoordFlag()
{
  const int   W = m_harmcoord->m_map_reso[0];
  const int   H = m_harmcoord->m_map_reso[1];
  const int   D = m_harmcoord->m_map_reso[2], WH = W*H;
  const float p = m_harmcoord->m_map_pitch  ;
  const EVec3f orig = m_harmcoord->m_map_pos;
  const byte* flg = m_harmcoord->m_map_flag ;

  glDisable(GL_LIGHTING );
  glPointSize(3);
  glBegin( GL_POINTS);
  for ( int z = 0; z < D; ++z) 
  {
    for ( int y = 0; y < H; ++y) 
    {
      for ( int x = 0; x < W; ++x) 
      {
        const int i = x + y * W + z * WH;
        if( flg[i] != 2) continue;
        float xx = orig[0] + (x + 0.5f) * p;
        float yy = orig[1] + (y + 0.5f) * p;
        float zz = orig[2] + (z + 0.5f) * p;
        //if( flg[i] == 0) glColor3d(0.2, 1.0, 0.2);
        //else 
        glColor3d(1.0, 1.0, 0.2);
        glVertex3f(xx,yy,zz);
      }
    }
  }
  glEnd();

}





int TCore::PickCageVertex( EVec3f ray_pos, EVec3f ray_dir )
{
  for( int i=0; i < m_cage.m_vSize; ++i) 
    if( t_distRayToPoint(ray_pos, ray_dir, m_cage.m_vVerts[i]) < CP_RADIUS ) return i;
  return -1;
}


void TCore::UpdateShapeByHC()
{
  const int cN = m_cage.m_vSize;

  for ( int i = 0; i < m_mesh.m_vSize; ++i)
  {
    const vector<float> &hc = m_harmcoord->m_verts_hc[i];
    EVec3f v(0,0,0);
    for ( int j = 0; j < cN; ++j) v += hc[j] * m_cage.m_vVerts[j];
    m_mesh.m_vVerts[i] = v;
  }

  m_mesh.updateNormal();
}



#pragma managed
