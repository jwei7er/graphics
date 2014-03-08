/*=========================================================================

 Author: Jordan Weiler
 Date:   May 21, 2013

 This program produces a dog written with OpenGL. The dog is supposed to
 be a Cairn Terrier. In the past, it has been mistaken for a cat.

=========================================================================*/

#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkInteractorStyle.h"
#include "vtkObjectFactory.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty.h"
#include "vtkCamera.h"
#include "vtkLight.h"
#include "vtkOpenGLPolyDataMapper.h"
#include "vtkJPEGReader.h"
#include "vtkImageData.h"

#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>

#include <vector>

class Triangle
{
	public:
		double         X[3];
		double         Y[3];
		double         Z[3];
};

class vtk441Mapper : public vtkOpenGLPolyDataMapper
{
	protected:
		GLuint displayList;
		bool   initialized;

	public:
		vtk441Mapper()
		{
			initialized = false;
		}

		void RemoveVTKOpenGLStateSideEffects()
		{
			float Info[4] = { 0, 0, 0, 1 };
			glLightModelfv(GL_LIGHT_MODEL_AMBIENT, Info);
			float ambient[4] = { 1,1, 1, 1.0 };
			glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
			float diffuse[4] = { 1, 1, 1, 1.0 };
			glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse);
			float specular[4] = { 1, 1, 1, 1.0 };
			glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
		}

		void SetupLight(void)
		{
			glEnable(GL_LIGHTING);
			glEnable(GL_LIGHT0);
			GLfloat diffuse0[4] = { 0.6, 0.6, 0.6, 1 };
			GLfloat ambient0[4] = { 0.2, 0.2, 0.2, 1 };
			GLfloat specular0[4] = { 0.0, 0.0, 0.0, 1 };
			GLfloat pos0[4] = { 0, .707, 0.707, 0 };
			glLightfv(GL_LIGHT0, GL_POSITION, pos0);
			glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse0);
			glLightfv(GL_LIGHT0, GL_AMBIENT, ambient0);
			glLightfv(GL_LIGHT0, GL_SPECULAR, specular0);
			glEnable(GL_LIGHT1);
			GLfloat pos1[4] = { .707, -.707, 0 };
			glLightfv(GL_LIGHT1, GL_POSITION, pos1);
			glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse0);
			glLightfv(GL_LIGHT1, GL_AMBIENT, ambient0);
			glLightfv(GL_LIGHT1, GL_SPECULAR, specular0);
			glDisable(GL_LIGHT2);
			glDisable(GL_LIGHT3);
			glDisable(GL_LIGHT5);
			glDisable(GL_LIGHT6);
			glDisable(GL_LIGHT7);
		}
};

class vtk441MapperPart3 : public vtk441Mapper
{
	public:
		static vtk441MapperPart3 *New();

		GLuint displayList;
		bool   initialized;
		double headSize; 	// absolute terms
		double eyeballSize; // as proportion of headSize
		double eyeHeight; 	// as proportion of radius above center of sphere
		double pupilSize; 	// as proportion of eyeballSize
		double eyeAngle; 	// from center-line, in degrees.  

		vtk441MapperPart3()
		{
			initialized = false;
			headSize = 4;
			eyeballSize = 0.15;
			eyeHeight = 0.25;
			pupilSize = 0.25;
			eyeAngle = 22;
		}
		
		void DrawCone()
		{
			int nfacets = 30;
			glBegin(GL_TRIANGLES);
			for (int i = 0; i < nfacets; i++)
			{
				double angle = 3.14159 * 2.0 * i / nfacets;
				double nextAngle = (i == nfacets-1 ? 0 : 3.14159 * 2.0 * (i + 1) / nfacets);
				glNormal3f(0, 0, 1);
				glVertex3f(0, 1, 1);
				glNormal3f(cos(nextAngle), -sin(nextAngle), 0);
				glVertex3f(cos(nextAngle), 1-sin(nextAngle), 0);
				glNormal3f(cos(angle), -sin(angle), 0);
				glVertex3f(cos(angle), 1-sin(angle), 0);
				glNormal3f(0,0,-1);
				glVertex3f(0,1,0);
				glVertex3f(cos(angle), 1-sin(angle), 0);
				glVertex3f(cos(nextAngle), 1-sin(nextAngle), 0);
			}
			glEnd();
		}

		void DrawCylinder()
		{
			int nfacets = 30;
			glBegin(GL_TRIANGLES);
			for (int i = 0 ; i < nfacets ; i++)
			{
				double angle = 3.14159*2.0*i/nfacets;
				double nextAngle = (i == nfacets-1 ? 0 : 3.14159*2.0*(i+1)/nfacets);
				glNormal3f(0,0,1);
				glVertex3f(0, 0, 1);
				glVertex3f(cos(angle), sin(angle), 1);
				glVertex3f(cos(nextAngle), sin(nextAngle), 1);
				glNormal3f(0,0,-1);
				glVertex3f(0, 0, 0);
				glVertex3f(cos(angle), -sin(angle), 0);
				glVertex3f(cos(nextAngle), -sin(nextAngle), 0);
			}
			glEnd();
			glBegin(GL_QUADS);
			for (int i = 0 ; i < nfacets ; i++)
			{
				double angle = 3.14159*2.0*i/nfacets;
				double nextAngle = (i == nfacets-1 ? 0 : 3.14159*2.0*(i+1)/nfacets);
				glNormal3f(cos(angle), sin(angle), 0);
				glVertex3f(cos(angle), sin(angle), 1);
				glVertex3f(cos(angle), sin(angle), 0);
				glNormal3f(cos(nextAngle), sin(nextAngle), 0);
				glVertex3f(cos(nextAngle), sin(nextAngle), 0);
				glVertex3f(cos(nextAngle), sin(nextAngle), 1);
			}
			glEnd();
		}
		
		void DrawEar(void)
		{
			int ntriangles = 15;
			glBegin(GL_TRIANGLES);
			for (int i = 0 ; i < ntriangles ; i++)
			{
				double angle = M_PI * i / ((double) ntriangles-1);
				double nextAngle = M_PI * (i+1) / ((double) ntriangles-1);
				glNormal3f(0, 1, 0);
				glVertex3f(0, 1, 0);
				glNormal3f(cos(angle), 0, -sin(angle));
				glVertex3f(cos(angle), 0, -sin(angle));
				glNormal3f(cos(nextAngle), 0, -sin(nextAngle));
				glVertex3f(cos(nextAngle), 0, -sin(nextAngle));
			}
			glEnd();
		}
		
		// used to recursively split sphere triangles into 4 smaller triangles so that the sphere is more realistic
		std::vector<Triangle> SplitTriangle(std::vector<Triangle> &list)
		{
			std::vector<Triangle> output(4*list.size());
			for (unsigned int i = 0 ; i < list.size() ; i++)
			{
				double mid1[3], mid2[3], mid3[3];
				mid1[0] = (list[i].X[0]+list[i].X[1])/2;
				mid1[1] = (list[i].Y[0]+list[i].Y[1])/2;
				mid1[2] = (list[i].Z[0]+list[i].Z[1])/2;
				mid2[0] = (list[i].X[1]+list[i].X[2])/2;
				mid2[1] = (list[i].Y[1]+list[i].Y[2])/2;
				mid2[2] = (list[i].Z[1]+list[i].Z[2])/2;
				mid3[0] = (list[i].X[0]+list[i].X[2])/2;
				mid3[1] = (list[i].Y[0]+list[i].Y[2])/2;
				mid3[2] = (list[i].Z[0]+list[i].Z[2])/2;
				output[4*i+0].X[0] = list[i].X[0];
				output[4*i+0].Y[0] = list[i].Y[0];
				output[4*i+0].Z[0] = list[i].Z[0];
				output[4*i+0].X[1] = mid1[0];
				output[4*i+0].Y[1] = mid1[1];
				output[4*i+0].Z[1] = mid1[2];
				output[4*i+0].X[2] = mid3[0];
				output[4*i+0].Y[2] = mid3[1];
				output[4*i+0].Z[2] = mid3[2];
				output[4*i+1].X[0] = list[i].X[1];
				output[4*i+1].Y[0] = list[i].Y[1];
				output[4*i+1].Z[0] = list[i].Z[1];
				output[4*i+1].X[1] = mid2[0];
				output[4*i+1].Y[1] = mid2[1];
				output[4*i+1].Z[1] = mid2[2];
				output[4*i+1].X[2] = mid1[0];
				output[4*i+1].Y[2] = mid1[1];
				output[4*i+1].Z[2] = mid1[2];
				output[4*i+2].X[0] = list[i].X[2];
				output[4*i+2].Y[0] = list[i].Y[2];
				output[4*i+2].Z[0] = list[i].Z[2];
				output[4*i+2].X[1] = mid3[0];
				output[4*i+2].Y[1] = mid3[1];
				output[4*i+2].Z[1] = mid3[2];
				output[4*i+2].X[2] = mid2[0];
				output[4*i+2].Y[2] = mid2[1];
				output[4*i+2].Z[2] = mid2[2];
				output[4*i+3].X[0] = mid1[0];
				output[4*i+3].Y[0] = mid1[1];
				output[4*i+3].Z[0] = mid1[2];
				output[4*i+3].X[1] = mid2[0];
				output[4*i+3].Y[1] = mid2[1];
				output[4*i+3].Z[1] = mid2[2];
				output[4*i+3].X[2] = mid3[0];
				output[4*i+3].Y[2] = mid3[1];
				output[4*i+3].Z[2] = mid3[2];
			}
			return output;
		}
		
		void DrawSphere()
		{
			int recursionLevel = 3;
			Triangle t;
			t.X[0] = 1;
			t.Y[0] = 0;
			t.Z[0] = 0;
			t.X[1] = 0;
			t.Y[1] = 1;
			t.Z[1] = 0;
			t.X[2] = 0;
			t.Y[2] = 0;
			t.Z[2] = 1;
			std::vector<Triangle> list;
			list.push_back(t);
			for (int r = 0 ; r < recursionLevel ; r++)
			{
				list = SplitTriangle(list);
			}

			// really draw
			for (int octent = 0 ; octent < 8 ; octent++)
			{
				glPushMatrix();
				glRotatef(90*(octent%4), 1, 0, 0);
				if (octent >= 4)
					glRotatef(180, 0, 0, 1);
				glBegin(GL_TRIANGLES);
				for (unsigned int i = 0 ; i < list.size() ; i++)
				{
					for (int j = 0 ; j < 3 ; j++)
					{
						double ptMag = sqrt(list[i].X[j]*list[i].X[j] + list[i].Y[j]*list[i].Y[j] + list[i].Z[j]*list[i].Z[j]);
						glNormal3f(list[i].X[j]/ptMag, list[i].Y[j]/ptMag, list[i].Z[j]/ptMag);
						glVertex3f(list[i].X[j]/ptMag, list[i].Y[j]/ptMag, list[i].Z[j]/ptMag);
					}
				}
				glEnd();
				glPopMatrix();
			}
		}

		void Pink(void) { glColor3ub(250, 128, 114); };
		void White(void) { glColor3ub(255, 255, 255); };
		void Black(void) { glColor3ub(0, 0, 0); };
		void DimGray(void) { glColor3ub(105, 105, 105); };

		void Brown1(void) { glColor3ub(185,157,130); }; //light brown
		void Brown2(void) { glColor3ub(173,141,109); };
		void Brown3(void) { glColor3ub(161,125,88); };
		void Brown4(void) { glColor3ub(150,108,67); };
		void Brown5(void) { glColor3ub(138,92,46); }; //dark brown
		
		void DrawHead(void)
		{
			glPushMatrix();
			Brown3();
			glScalef(headSize, headSize, headSize);
			DrawSphere();

			// draw eyeballs
			double angleRadians = eyeAngle/360.0 * 2 * 3.14159;
			double mag = sqrt(1*1-eyeHeight*eyeHeight)*(1-eyeballSize/2);

			glPushMatrix();
			glTranslatef(mag*cos(angleRadians), mag*sin(angleRadians), eyeHeight);
			DrawEyeball();
			glPopMatrix();

			glPushMatrix();
			glTranslatef(mag*cos(-angleRadians), mag*sin(-angleRadians), eyeHeight);
			DrawEyeball();
			glPopMatrix();

			glPopMatrix();
		}
		
		void DrawEyeball(void)
		{
			glPushMatrix();
			White();
			glScalef(eyeballSize, eyeballSize, eyeballSize);
			DrawSphere();
			Black();
			glTranslatef(1-pupilSize/2, 0, 0);
			glScalef(pupilSize, pupilSize, pupilSize);
			DrawSphere();
			glPopMatrix();
		}
		
		void DrawNeck(void)
		{
			glPushMatrix();
			Brown2();
			glRotatef(225,0,1,0);
			glTranslatef(0,0,2);
			glScalef(2,2,4);
			DrawCylinder();
			glPopMatrix();
		}
		
		void DrawBody(void)
		{
			glPushMatrix();
			Brown3();
			glTranslatef(-10,0,-4);
			glScalef(8,4,4);
			DrawSphere();
			glPopMatrix();
		}
		
		void DrawTail(void)
		{
			glPushMatrix();
			Brown2();
			glTranslatef(-17,0,-4);
			glRotatef(-30,0,1,0);
			glScalef(0.8,0.8,5);
			DrawCylinder();
			glPopMatrix();
			
			// draw pointy tip of tail
			glPushMatrix();
			Brown1();
			glTranslatef(-17,0,-4);
			glRotatef(-30,0,1,0);
			glTranslatef(0,-0.8,5);
			glScalef(0.8,0.8,2);
			DrawCone();
			glPopMatrix();
		}
		
		void DrawEars(void)
		{
			glPushMatrix();
			Brown4();
			glRotatef(-15,1,0,0);
			glTranslatef(0,1.5,3);
			glRotatef(150,0,0,1);
			glRotatef(90,1,0,0);
			glScalef(1.5,4,1.5);
			DrawEar();
			glPopMatrix();
			
			glPushMatrix();
			glRotatef(15,1,0,0);
			glTranslatef(0,-1.5,3);
			glRotatef(30,0,0,1);
			glRotatef(90,1,0,0);
			glScalef(1.5,4,1.5);
			DrawEar();
			
			glPopMatrix();
		}
		
		void DrawLegs(void)
		{
			glPushMatrix();
			Brown4();
			glTranslatef(-6,2,-11);
			glScalef(0.75,0.75,6);
			DrawCylinder();
			
			glPushMatrix();
			glTranslatef(0,-5.33,0);
			DrawCylinder();
			
			glTranslatef(-12,0,0);
			DrawCylinder();
			glPopMatrix();
			
			glTranslatef(-12,0,0);
			DrawCylinder();
			glPopMatrix();
			
		}
		
		void DrawFeet(void)
		{
			glPushMatrix();
			Brown5();
			glTranslatef(-5.5,2,-11);
			glScalef(1.5,1,0.5);
			DrawSphere();
			
			glPushMatrix();
			glTranslatef(0,-4,0);
			DrawSphere();
			
			glTranslatef(-6,0,0);
			DrawSphere();
			
			glPopMatrix();
			glTranslatef(-6,0,0);
			DrawSphere();
			
			glPopMatrix();
			
		}
		
		void DrawMouth(void)
		{
			glPushMatrix();
			DimGray();
			glTranslatef(3.2,0,-0.6);
			glScalef(1.7,1.5,0.5);
			DrawSphere();
			glPopMatrix();
			
			glPushMatrix();
			glTranslatef(3.2,0,-1.6);
			glScalef(1.4,1.2,0.5);
			DrawSphere();
			glPopMatrix();
		}
			
		void DrawTongue(void)
		{
			glPushMatrix();
			Pink();
			glTranslatef(2.5,0,-1.1);
			glScalef(2,0.8,0.4);
			DrawSphere();
			glPopMatrix();
		}
			
		void DrawNose(void)
		{
			glPushMatrix();
			Black();
			glTranslatef(4.6,0,-0.2);
			glScalef(0.4,0.4,0.4);
			DrawSphere();
			glPopMatrix();
		}
		
		virtual void RenderPiece(vtkRenderer *ren, vtkActor *act)
		{
			RemoveVTKOpenGLStateSideEffects();
			SetupLight();
			glEnable(GL_COLOR_MATERIAL);

			glMatrixMode(GL_MODELVIEW);

			DrawHead();
			
			DrawNeck();
			
			DrawBody();
			
			DrawTail();
			
			DrawEars();
			
			DrawLegs();
			
			DrawFeet();
			
			DrawMouth();
			
			DrawTongue();
			
			DrawNose();
		}
};

vtkStandardNewMacro(vtk441MapperPart3);


int main()
{
	// Dummy input so VTK pipeline mojo is happy.
	vtkSmartPointer<vtkSphereSource> sphere = vtkSmartPointer<vtkSphereSource>::New();
	sphere->SetThetaResolution(100);
	sphere->SetPhiResolution(50);

	vtkSmartPointer<vtk441MapperPart3> win3Mapper = vtkSmartPointer<vtk441MapperPart3>::New();
	win3Mapper->SetInputConnection(sphere->GetOutputPort());

	vtkSmartPointer<vtkActor> win3Actor = vtkSmartPointer<vtkActor>::New();
	win3Actor->SetMapper(win3Mapper);

	vtkSmartPointer<vtkRenderer> ren3 = vtkSmartPointer<vtkRenderer>::New();

	vtkSmartPointer<vtkRenderWindow> renWin = vtkSmartPointer<vtkRenderWindow>::New();
	renWin->AddRenderer(ren3);
	ren3->SetViewport(0, 0, 1, 1);

	vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	iren->SetRenderWindow(renWin);

	// Add the actors to the renderer, set the background and size.
	//
	bool doWindow3 = true;
	if (doWindow3)
		ren3->AddActor(win3Actor);
	ren3->SetBackground(0.3, 0.3, 0.3);
	renWin->SetSize(500, 500);

	// Set up the lighting.
	//
	vtkSmartPointer<vtkLight> light = vtkSmartPointer<vtkLight>::New();
	light->SetFocalPoint(1.875,0.6125,0);
	light->SetPosition(0.875,1.6125,1);
	ren3->AddLight(light);

	ren3->GetActiveCamera()->SetFocalPoint(-10,0,-5);
	ren3->GetActiveCamera()->SetPosition(0,0,70);
	ren3->GetActiveCamera()->SetViewUp(0,1,0);
	ren3->GetActiveCamera()->SetClippingRange(20, 120);
	ren3->GetActiveCamera()->SetDistance(70);

	// This starts the event loop and invokes an initial render.
	//
	((vtkInteractorStyle *)iren->GetInteractorStyle())->SetAutoAdjustCameraClippingRange(0);
	iren->Initialize();
	iren->Start();

	return EXIT_SUCCESS;
}
