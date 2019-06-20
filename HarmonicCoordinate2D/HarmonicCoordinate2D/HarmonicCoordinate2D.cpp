#include "stdafx.h"

#include "MainForm.h"

using namespace System;



//TODO LIST
//1. gen form               OK
//2. attach OpenGL to form  OK     
//3. mouse event            OK
//4. gen 2D object and cage 
//6. compute binding 
//7. compute binding in hierarchical structure
//8. deformation

using namespace HarmonicCoordinate2D;



[STAThreadAttribute]
int main(array<System::String ^> ^args)
{

    MainForm::getInst()->ShowDialog();

    return 0;
}
