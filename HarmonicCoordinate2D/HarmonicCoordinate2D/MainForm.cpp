#include "stdafx.h"
#include "MainForm.h"
#include "TCore.h"


#pragma comment( lib, "opengl32.lib" )
#pragma comment( lib, "glu32.lib" )
#pragma comment( lib, "gdi32.lib" )
#pragma comment( lib, "User32.lib" )


using namespace HarmonicCoordinate2D;

MainForm::MainForm()
{
  m_ogl = 0;
  InitializeComponent();
  m_ogl = new OglForCLI(GetDC((HWND)mainpanel->Handle.ToPointer()));
  m_ogl->SetCam(EVec3f(0,0,10), EVec3f(0,0,0), EVec3f(0,1,0));
}



System::Void MainForm::mainpanel_MouseUp(System::Object^  sender, System::Windows::Forms::MouseEventArgs^  e) 
 {
  if (e->Button == System::Windows::Forms::MouseButtons::Left  ) TCore::getInst()->LBtnUp(EVec2i(e->X, e->Y), m_ogl);
  if (e->Button == System::Windows::Forms::MouseButtons::Middle) TCore::getInst()->MBtnUp(EVec2i(e->X, e->Y), m_ogl);
  if (e->Button == System::Windows::Forms::MouseButtons::Right ) TCore::getInst()->RBtnUp(EVec2i(e->X, e->Y), m_ogl);
}


System::Void MainForm::mainpanel_MouseMove(System::Object^  sender, System::Windows::Forms::MouseEventArgs^  e) 
{
  TCore::getInst()->MouseMove(EVec2i(e->X, e->Y), m_ogl);
}

System::Void MainForm::mainpanel_MouseDown(System::Object^  sender, System::Windows::Forms::MouseEventArgs^  e) 
{
  if (e->Button == System::Windows::Forms::MouseButtons::Left  ) TCore::getInst()->LBtnDown(EVec2i(e->X, e->Y), m_ogl);
  if (e->Button == System::Windows::Forms::MouseButtons::Middle) TCore::getInst()->MBtnDown(EVec2i(e->X, e->Y), m_ogl);
  if (e->Button == System::Windows::Forms::MouseButtons::Right ) TCore::getInst()->RBtnDown(EVec2i(e->X, e->Y), m_ogl);
}

System::Void MainForm::MainForm_KeyDown(System::Object^  sender, System::Windows::Forms::KeyEventArgs^  e) 
{
  TCore::getInst()->KeyDown(e->KeyValue);
}







void MainForm::RedrawMainPanel()
{
  if( m_ogl == 0 )return;

  m_ogl->SetBgColor(0.0f, 0.3f, 0.3f, 0.5f);


  const float scene_size = 10;
  m_ogl->OnDrawBegin_ortho( mainpanel->Width, mainpanel->Height, scene_size );
 
  TCore::getInst()->DrawScene();

  m_ogl->OnDrawEnd();
}
