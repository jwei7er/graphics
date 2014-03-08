/*=========================================================================

 Author: Jordan Weiler
 Date:   April 10, 2013

 This program reads triangles from a file and creates an image. For each
 triangle, the scanline algorithm is applied to fill up the color buffer.

=========================================================================*/

#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>

using std::cerr;
using std::endl;

double ceil441(double f)
{
    return ceil(f-0.00001);
}

double floor441(double f)
{
    return floor(f+0.00001);
}
	
double max(double a, double b)
{
	return (b > a) ? b: a;
}

double min(double a, double b)
{
	return (a > b) ? b: a;
}

double max(double a, double b, double c)
{
	return (b > a) ? max(b, c) : max(a, c);
}

double min(double a, double b, double c)
{
	return (a > b) ? min(b, c) : min(a, c);
}

vtkImageData *
NewImage(int width, int height)
{
	vtkImageData *image = vtkImageData::New();
	image->SetDimensions(width, height, 1);
	image->SetWholeExtent(0, width-1, 0, height-1, 0, 0);
	image->SetUpdateExtent(0, width-1, 0, height-1, 0, 0);
	image->SetNumberOfScalarComponents(3);
	image->SetScalarType(VTK_UNSIGNED_CHAR);
	image->AllocateScalars();

	return image;
}

void
WriteImage(vtkImageData *img, const char *filename)
{
	std::string full_filename = filename;
	full_filename += ".png";
	vtkPNGWriter *writer = vtkPNGWriter::New();
	writer->SetInput(img);
	writer->SetFileName(full_filename.c_str());
	writer->Write();
	writer->Delete();
}

class Screen
{
	public:
		unsigned char   *buffer;
		int width, height;

	// would some methods for accessing and setting pixels be helpful?
};

class TriangleSide
{
	public:
		double x1, y1;
		double x2, y2;
		double slope;
		double yIntercept;
		bool vertical, ignore;
		double minY, maxY;
		
		TriangleSide(double nX1, double nY1, double nX2, double nY2)
		{
			x1 = nX1;
			y1 = nY1;
			x2 = nX2;
			y2 = nY2;
			
			minY = min(y1, y2);
			maxY = max(y1, y2);
			
			if (x2 - x1 == 0) {
				vertical = true;
				slope = x1;
			}
			else {
				vertical = false;
				slope = (y2 - y1) / (x2 - x1);
			}
			yIntercept = y1 - slope * x1;
			if (y2 - y1 == 0)
				ignore = true;
			else
				ignore = false;
		}
		
		double getXValueGivenY(double y)
		{
			if (vertical == true)
				return slope;
			else {
				if (slope == 0)
					return 0;
				else
					return (y - yIntercept) / slope;
			}
		}
		
		bool withinBounds(int y)
		{
			if (minY <= y and y <= maxY)
				return true;
			else
				return false;
		}
};

class Triangle
{
	public:
		double X[3];
		double Y[3];
		unsigned char color[3];
		Screen screen;

	void runScanline()
	{
		double minX = max(ceil441(min(X[0], X[1], X[2])),0);
		double maxX = min(floor441(max(X[0], X[1], X[2])),screen.width-1);
		
		double minY = max(ceil441(min(Y[0], Y[1], Y[2])),0);
		double maxY = min(floor441(max(Y[0], Y[1], Y[2])),screen.height-1);
		
		//cout << "(" << X[0] << "," << Y[0] << ") (" << X[1] << "," << Y[1] << ") (" << X[2] << "," << Y[2] << ")\n";
		
		//cout << "minX: " << minX << " maxX: " << maxX << " minY: " << minY << " maxY: " << maxY << "\n";
		
		TriangleSide ts1 = TriangleSide(X[0], Y[0], X[1], Y[1]);
		//cout << "side 1: y = " << ts1.slope << "x + " << ts1.yIntercept << "\n";
		TriangleSide ts2 = TriangleSide(X[1], Y[1], X[2], Y[2]);
		//cout << "side 2: y = " << ts2.slope << "x + " << ts2.yIntercept << "\n";
		TriangleSide ts3 = TriangleSide(X[2], Y[2], X[0], Y[0]);
		//cout << "side 3: y = " << ts3.slope << "x + " << ts3.yIntercept << "\n";
		
		int offset;
		double minXGivenY, maxXGivenY, tempX;
		for (int y = minY; y <= maxY; y++)
		{
			offset = screen.width * 3 * y;
			minXGivenY = screen.width;
			maxXGivenY = 0;
			if (!ts1.ignore and ts1.withinBounds(y)) {
				tempX = ts1.getXValueGivenY(y);
				//cout << "ts1: " << tempX << " ";
				minXGivenY = ts1.getXValueGivenY(y);
				maxXGivenY = minXGivenY;
			}
			if (!ts2.ignore and ts2.withinBounds(y)) {
				tempX = ts2.getXValueGivenY(y);
				//cout << "ts2: " << tempX << " ";
				if (tempX < minXGivenY) {
					minXGivenY = tempX;
				}
				if (tempX > maxXGivenY) {
					maxXGivenY = tempX;
				}
			}
			if (!ts3.ignore and ts3.withinBounds(y)) {
				tempX = ts3.getXValueGivenY(y);
				//cout << "ts3: " << tempX << " ";
				if (tempX < minXGivenY) {
					minXGivenY = tempX;
				}
				if (tempX > maxXGivenY) {
					maxXGivenY = tempX;
				}
			}
			
			//cout << "y: " << y << " min x: " << minXGivenY << " max x: " << maxXGivenY << "\n";
			minXGivenY = max(ceil441(minXGivenY), minX);
			maxXGivenY = min(floor441(maxXGivenY), maxX);
			
			for (int x = minXGivenY; x <= maxXGivenY; x++)
			{
				screen.buffer[offset + (x * 3)] = color[0];
				screen.buffer[offset + (x * 3) + 1] = color[1];
				screen.buffer[offset + (x * 3) + 2] = color[2];
			}
		}
		//cout << "\n\n";
	}
};

std::vector<Triangle>
GetTriangles(void)
{
    vtkPolyDataReader *rdr = vtkPolyDataReader::New();
    rdr->SetFileName("proj1c_geometry.vtk");
    cerr << "Reading" << endl;
    rdr->Update();
    cerr << "Done reading" << endl;
    if (rdr->GetOutput()->GetNumberOfCells() == 0)
    {
        cerr << "Unable to open file!!" << endl;
        exit(EXIT_FAILURE);
    }
    vtkPolyData *pd = rdr->GetOutput();
    int numTris = pd->GetNumberOfCells();
    vtkPoints *pts = pd->GetPoints();
    vtkCellArray *cells = pd->GetPolys();
    vtkFloatArray *colors = (vtkFloatArray *) pd->GetPointData()->GetArray("color_nodal");
    float *color_ptr = colors->GetPointer(0);
    std::vector<Triangle> tris(numTris);
    vtkIdType npts;
    vtkIdType *ptIds;
    int idx;
    for (idx = 0, cells->InitTraversal() ; cells->GetNextCell(npts, ptIds) ; idx++)
    {
        if (npts != 3)
        {
            cerr << "Non-triangles!! ???" << endl;
            exit(EXIT_FAILURE);
        }
        tris[idx].X[0] = pts->GetPoint(ptIds[0])[0];
        tris[idx].X[1] = pts->GetPoint(ptIds[1])[0];
        tris[idx].X[2] = pts->GetPoint(ptIds[2])[0];
        tris[idx].Y[0] = pts->GetPoint(ptIds[0])[1];
        tris[idx].Y[1] = pts->GetPoint(ptIds[1])[1];
        tris[idx].Y[2] = pts->GetPoint(ptIds[2])[1];
        tris[idx].color[0] = (unsigned char) color_ptr[4*ptIds[0]+0];
        tris[idx].color[1] = (unsigned char) color_ptr[4*ptIds[0]+1];
        tris[idx].color[2] = (unsigned char) color_ptr[4*ptIds[0]+2];
    }

    return tris;
}

int main()
{
	int imageWidth = 1786;
	int imageHeight = 1344;
	vtkImageData *image = NewImage(imageWidth, imageHeight);
	unsigned char *buffer = 
		(unsigned char *) image->GetScalarPointer(0,0,0);
	int npixels = imageWidth*imageHeight;
	for (int i = 0 ; i < npixels*3 ; i++)
		buffer[i] = 0;

	std::vector<Triangle> triangles = GetTriangles();

	Screen screen;
	screen.buffer = buffer;
	screen.width = imageWidth;
	screen.height = imageHeight;

	for (int i = 0; i < triangles.size(); i++) //triangles.size(); i++)
	{
		Triangle t = triangles[i];
		t.screen = screen;
		t.runScanline();
	}

	WriteImage(image, "allTriangles");
}
