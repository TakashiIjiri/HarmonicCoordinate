#pragma once

#include "COMMON/tmath.h"
#include "COMMON/OglForCLI.h"



namespace HarmonicCoordinate3D {

	using namespace System;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Drawing;

	/// <summary>
	/// MainForm の概要
	/// </summary>
	public ref class MainForm : public System::Windows::Forms::Form
	{
  private:
    private: System::Windows::Forms::Panel^  mainpanel;


  private:
    OglForCLI *m_ogl;
    static MainForm^ m_singleton;
		MainForm(void);

  public:
    static MainForm^ getInst()
    {
      if (m_singleton == nullptr) m_singleton = gcnew MainForm();
      return m_singleton;
    }

    void RedrawMainPanel();


	public:

	protected:
		/// <summary>
		/// 使用中のリソースをすべてクリーンアップします。
		/// </summary>
		~MainForm()
		{
			if (components)
			{
				delete components;
			}
		}
  protected:

	private:
		/// <summary>
		/// 必要なデザイナー変数です。
		/// </summary>
		System::ComponentModel::Container ^components;

#pragma region Windows Form Designer generated code
		/// <summary>
		/// デザイナー サポートに必要なメソッドです。このメソッドの内容を
		/// コード エディターで変更しないでください。
		/// </summary>
		void InitializeComponent(void)
		{
      this->mainpanel = (gcnew System::Windows::Forms::Panel());
      this->SuspendLayout();
      // 
      // mainpanel
      // 
      this->mainpanel->Anchor = static_cast<System::Windows::Forms::AnchorStyles>((((System::Windows::Forms::AnchorStyles::Top | System::Windows::Forms::AnchorStyles::Bottom)
        | System::Windows::Forms::AnchorStyles::Left)
        | System::Windows::Forms::AnchorStyles::Right));
      this->mainpanel->Location = System::Drawing::Point(2, 12);
      this->mainpanel->Name = L"mainpanel";
      this->mainpanel->Size = System::Drawing::Size(754, 603);
      this->mainpanel->TabIndex = 0;
      this->mainpanel->Paint += gcnew System::Windows::Forms::PaintEventHandler(this, &MainForm::mainpanel_Paint);
      this->mainpanel->MouseDown += gcnew System::Windows::Forms::MouseEventHandler(this, &MainForm::mainpanel_MouseDown);
      this->mainpanel->MouseMove += gcnew System::Windows::Forms::MouseEventHandler(this, &MainForm::mainpanel_MouseMove);
      this->mainpanel->MouseUp += gcnew System::Windows::Forms::MouseEventHandler(this, &MainForm::mainpanel_MouseUp);
      this->mainpanel->Resize += gcnew System::EventHandler(this, &MainForm::mainpanel_Resize);
      // 
      // MainForm
      // 
      this->AutoScaleDimensions = System::Drawing::SizeF(6, 12);
      this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
      this->ClientSize = System::Drawing::Size(758, 622);
      this->Controls->Add(this->mainpanel);
      this->Name = L"MainForm";
      this->Text = L"MainForm";
      this->KeyDown += gcnew System::Windows::Forms::KeyEventHandler(this, &MainForm::MainForm_KeyDown);
      this->ResumeLayout(false);

    }
#pragma endregion
  private: System::Void mainpanel_MouseDown(System::Object^  sender, System::Windows::Forms::MouseEventArgs^  e) ;
  private: System::Void mainpanel_MouseUp  (System::Object^  sender, System::Windows::Forms::MouseEventArgs^  e) ;
  private: System::Void mainpanel_MouseMove(System::Object^  sender, System::Windows::Forms::MouseEventArgs^  e) ;
  private: System::Void mainpanel_Resize(System::Object^  sender, System::EventArgs^  e) {
    if ( m_ogl != 0) RedrawMainPanel();
  }
  private: System::Void mainpanel_Paint(System::Object^  sender, System::Windows::Forms::PaintEventArgs^  e) {
    if ( m_ogl != 0) RedrawMainPanel();
  }
  private: System::Void MainForm_KeyDown(System::Object^  sender, System::Windows::Forms::KeyEventArgs^  e) ;
  };

  //interfaces
  inline void MainForm_RedrawMainPanel(){ MainForm::getInst()->RedrawMainPanel(); }
}
