#include "stdafx.h"

#include "MainForm.h"

using namespace System;

using namespace HarmonicCoordinate3D;


[STAThreadAttribute]
int main(array<System::String ^> ^args)
{
    MainForm::getInst()->ShowDialog();

    return 0;
}

//TODO LIST
//1. gen form                OK
//2. attach OpenGL to form   OK  
//3. mouse event             OK 
//4. load 3D object and cage TODO 
//5. cage deformation ui     TODO
//6. compute binding         TODO
//7. compute binding in hierarchical structure  TODO
//8. deformation  TODO


