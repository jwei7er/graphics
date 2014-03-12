/*=========================================================================

 Author: Jordan Weiler
 Date:   April 9, 2013

 This program creates an image with 100 triangles. For each triangle, 
 the scanline algorithm is used to fill up the image buffer with colors.

=========================================================================*/

#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>

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

vtkImageData * NewImage(int width, int height)
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

void WriteImage(vtkImageData *img, const char *filename)
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
        unsigned char *buffer;
        int width, height;
};

class TriangleSide
{
    // class to solve the y = mx + b equation for a side of a triangle
    public:
    double x1, y1;
    double x2, y2;
    double slope;
    double yIntercept;
    bool vertical, ignore;

    TriangleSide(double nX1, double nY1, double nX2, double nY2)
    {
        x1 = nX1;
        y1 = nY1;
        x2 = nX2;
        y2 = nY2;
        
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
        if (vertical == true) {
            return slope;
        }
        else {
            if (slope == 0) {
                return 0;
            }
            else {
                return (y - yIntercept) / slope;
            }
        }
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
            // Scanline algorithm to find the x coordinates of a triangle
            double minX = max(min(X[0], X[1], X[2]),0);
            double maxX = min(max(X[0], X[1], X[2]),screen.width);
    
            double minY = max(min(Y[0], Y[1], Y[2]),0);
            double maxY = min(max(Y[0], Y[1], Y[2]),screen.height);
    
            TriangleSide ts1 = TriangleSide(X[0], Y[0], X[1], Y[1]);
            TriangleSide ts2 = TriangleSide(X[1], Y[1], X[2], Y[2]);
            TriangleSide ts3 = TriangleSide(X[2], Y[2], X[0], Y[0]);
    
            int offset;
            double minXGivenY, maxXGivenY, tempX;
            for (int y = minY; y <= maxY; y++)
            {
                offset = screen.width * 3 * y;
                minXGivenY = screen.width;
                maxXGivenY = 0;
                if (!ts1.ignore) {
                    tempX = ts1.getXValueGivenY(y);
                    minXGivenY = tempX;
                    maxXGivenY = tempX;
                }
                if (!ts2.ignore) {
                    tempX = ts2.getXValueGivenY(y);
                    if (tempX < minXGivenY)
                        minXGivenY = tempX;
                    if (tempX > maxXGivenY)
                        maxXGivenY = tempX;
                }
                if (!ts3.ignore) {
                    tempX = ts3.getXValueGivenY(y);
                    if (tempX < minXGivenY)
                        minXGivenY = tempX;
                    if (tempX > maxXGivenY)
                        maxXGivenY = tempX;
                }
        
                minXGivenY = max(ceil441(minXGivenY), minX);
                maxXGivenY = min(floor441(maxXGivenY), maxX);
        
                // only loop over minX and maxX range for given Y value to fill in buffer
                for (int x = minXGivenY; x <= maxXGivenY; x++) {
                    screen.buffer[offset + (x * 3)] = color[0];
                    screen.buffer[offset + (x * 3) + 1] = color[1];
                    screen.buffer[offset + (x * 3) + 2] = color[2];
                }
            }
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
};

std::vector<Triangle> GetTriangles(void)
{
    std::vector<Triangle> rv(100);

    unsigned char colors[6][3] = { {255,128,0}, {255, 0, 127}, {0,204,204}, 
                                    {76,153,0}, {255, 204, 204}, {204, 204, 0}};
    for (int i = 0 ; i < 100 ; i++)
    {
        int idxI = i%10;
        int posI = idxI*100;
        int idxJ = i/10;
        int posJ = idxJ*100;
        int firstPt = (i%3);
        rv[i].X[firstPt] = posI;
        if (i == 50)
            rv[i].X[firstPt] = -10;
        rv[i].Y[firstPt] = posJ;
        rv[i].X[(firstPt+1)%3] = posI+99;
        rv[i].Y[(firstPt+1)%3] = posJ;
        rv[i].X[(firstPt+2)%3] = posI+i;
        rv[i].Y[(firstPt+2)%3] = posJ+10*(idxJ+1);
        if (i == 95)
            rv[i].Y[(firstPt+2)%3] = 1050;
        rv[i].color[0] = colors[i%6][0];
        rv[i].color[1] = colors[i%6][1];
        rv[i].color[2] = colors[i%6][2];
    }

    return rv;
}

int main()
{
    vtkImageData *image = NewImage(1000, 1000);
    unsigned char *buffer = (unsigned char *) image->GetScalarPointer(0,0,0);
    int npixels = 1000*1000;
    for (int i = 0 ; i < npixels*3 ; i++)
        buffer[i] = 0;

    std::vector<Triangle> triangles = GetTriangles();

    Screen screen;
    screen.buffer = buffer;
    screen.width = 1000;
    screen.height = 1000;

    for (int i = 0; i < 100; i++)
    {
        Triangle t = triangles[i];
        t.screen = screen;
        t.runScanline();
    }

    WriteImage(image, "allTriangles");
}
