/*=========================================================================

 Author: Jordan Weiler
 Date:   April 19, 2013

 This program reads in triangles from a file and produces an image. Each
 vertex of each triangle has a specific color. Interpolation occurs for
 each triangle to determine the color of every pixel in the triangle.

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
        double *zBuffer;
        int width, height;

        void setZBuffer()
        {
            zBuffer = new double[width*height];
    
            for (int i = 0; i < width * height; i++)
            {
                zBuffer[i] = -1.0;
            }
        }
};

class TriangleSide
{
    public:
        double x1, y1, z1, r1, g1, b1;
        double x2, y2, z2, r2, g2, b2;
        double slope;
        double yIntercept;
        bool vertical, ignore;
        double minY, maxY;
    
        TriangleSide()
        {
            // Empty constructor
        }
    
        TriangleSide(double nX1, double nY1, double nZ1, double nR1, double nG1, double nB1, double nX2, double nY2, double nZ2, double nR2, double nG2, double nB2)
        {
            x1 = nX1;
            y1 = nY1;
            z1 = nZ1;
            r1 = nR1;
            g1 = nG1;
            b1 = nB1;
            x2 = nX2;
            y2 = nY2;
            z2 = nZ2;
            r2 = nR2;
            g2 = nG2;
            b2 = nB2;
        
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

        double interpolation(double v1, double v2, double A, double B, double v)
        {
            return A + ((v - v1) / (v2 - v1)) * (B - A);
        }
    
        double getZValueGivenY(double y)
        {
            return interpolation(y1, y2, z1, z2, y);
        }
    
        double getRValueGivenY(double y)
        {
            return interpolation(y1, y2, r1, r2, y);
        }
    
        double getGValueGivenY(double y)
        {
            return interpolation(y1, y2, g1, g2, y);
        }
    
        double getBValueGivenY(double y)
        {
            return interpolation(y1, y2, b1, b2, y);
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
        double Z[3];
        double colors[3][3];
        Screen screen;

        void runScanline(bool debug)
        {
            double minX = max(ceil441(min(X[0], X[1], X[2])),0);
            double maxX = min(floor441(max(X[0], X[1], X[2])),screen.width-1);
    
            double minY = max(ceil441(min(Y[0], Y[1], Y[2])),0);
            double maxY = min(floor441(max(Y[0], Y[1], Y[2])),screen.height-1);
    
            if (debug && (((abs(X[0] - 551.02) < 0.00001) && (abs(Y[0] - 877.551) < 0.00001)) ||
                          ((abs(X[1] - 551.02) < 0.00001) && (abs(Y[1] - 877.551) < 0.00001)) ||
                          ((abs(X[2] - 551.02) < 0.00001) && (abs(Y[2] - 877.551) < 0.00001))))
            {
                cout << "(" << X[0] << "," << Y[0] << "," << Z[0] << ") (" << X[1] << "," << Y[1] << "," << Z[1] << ") (" << X[2] << "," << Y[2] << "," << Z[2] << ")\n";
                cout << "(" << colors[0][0] << "," << colors[0][1] << "," << colors[0][2] << ")";
                cout << "(" << colors[1][0] << "," << colors[1][1] << "," << colors[1][2] << ")";
                cout << "(" << colors[2][0] << "," << colors[2][1] << "," << colors[2][2] << ")\n\n";
            }
            
            TriangleSide ts1 = TriangleSide(X[0], Y[0], Z[0], colors[0][0], colors[0][1], colors[0][2], X[1], Y[1], Z[1], colors[1][0], colors[1][1], colors[1][2]);
            TriangleSide ts2 = TriangleSide(X[1], Y[1], Z[1], colors[1][0], colors[1][1], colors[1][2], X[2], Y[2], Z[2], colors[2][0], colors[2][1], colors[2][2]);
            TriangleSide ts3 = TriangleSide(X[2], Y[2], Z[2], colors[2][0], colors[2][1], colors[2][2], X[0], Y[0], Z[0], colors[0][0], colors[0][1], colors[0][2]);
            
            int offset, zOffset;
            double minXGivenY, maxXGivenY, tempX;
            TriangleSide tsLeft, tsRight;
            for (int y = minY; y <= maxY; y++)
            {
                offset = screen.width * 3 * y;
                zOffset = screen.width * y;
                minXGivenY = screen.width;
                maxXGivenY = 0;
                if (!ts1.ignore and ts1.withinBounds(y)) {
                    tempX = ts1.getXValueGivenY(y);
                    minXGivenY = ts1.getXValueGivenY(y);
                    maxXGivenY = minXGivenY;
                    tsLeft = ts1;
                    tsRight = ts1;
                }
                if (!ts2.ignore and ts2.withinBounds(y)) {
                    tempX = ts2.getXValueGivenY(y);
                    if (tempX < minXGivenY) {
                        minXGivenY = tempX;
                        tsLeft = ts2;
                    }
                    if (tempX > maxXGivenY) {
                        maxXGivenY = tempX;
                        tsRight = ts2;
                    }
                }
                if (!ts3.ignore and ts3.withinBounds(y)) {
                    tempX = ts3.getXValueGivenY(y);
                    if (tempX < minXGivenY) {
                        minXGivenY = tempX;
                        tsLeft = ts3;
                    }
                    if (tempX > maxXGivenY) {
                        maxXGivenY = tempX;
                        tsRight = ts3;
                    }
                }
                
                minXGivenY = max(ceil441(minXGivenY), minX);
                maxXGivenY = min(floor441(maxXGivenY), maxX);
        
                double leftZ = tsLeft.getZValueGivenY(y);
                double leftR = tsLeft.getRValueGivenY(y);
                double leftG = tsLeft.getGValueGivenY(y);
                double leftB = tsLeft.getBValueGivenY(y);
                double rightZ = tsRight.getZValueGivenY(y);
                double rightR = tsRight.getRValueGivenY(y);
                double rightG = tsRight.getGValueGivenY(y);
                double rightB = tsRight.getBValueGivenY(y);
        
                for (int x = minXGivenY; x <= maxXGivenY; x++)
                {
                    double proportion;
                    if (minXGivenY != maxXGivenY)
                        proportion = (((double) x) - minXGivenY) / (maxXGivenY - minXGivenY);
                    else
                        proportion = 1.0;
                    
                    double z = leftZ + proportion * (rightZ - leftZ);
                    if (screen.zBuffer[zOffset + x] <= z)
                    {
                        if (debug) {
                            cout << "assigning pixel (" << x << "," << y << ")";
                            cout << "R:" << ceil441(255.0 * (leftR + proportion * (rightR - leftR))) << " G:" << ceil441(255.0 * (leftG + proportion * (rightG - leftG))) << " B:" << ceil441(255.0 * (leftB + proportion * (rightB - leftB))) << " Z:" << z << "\n";
                        }
                        screen.zBuffer[zOffset + x] = z;
                        screen.buffer[offset + (x * 3)] = (unsigned char) ceil441(255.0 * (leftR + proportion * (rightR - leftR)));
                        screen.buffer[offset + (x * 3) + 1] = (unsigned char) ceil441(255.0 * (leftG + proportion * (rightG - leftG)));
                        screen.buffer[offset + (x * 3) + 2] = (unsigned char) ceil441(255.0 * (leftB + proportion * (rightB - leftB)));
                    }
                }
            }
        }
};

std::vector<Triangle> GetTriangles(void)
{
    vtkPolyDataReader *rdr = vtkPolyDataReader::New();
    rdr->SetFileName("proj1d_geometry.vtk");
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
    vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray("hardyglobal");
    float *color_ptr = var->GetPointer(0);
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
        tris[idx].Z[0] = pts->GetPoint(ptIds[0])[2];
        tris[idx].Z[1] = pts->GetPoint(ptIds[1])[2];
        tris[idx].Z[2] = pts->GetPoint(ptIds[2])[2];
        
        // 1->2 interpolate between light blue, dark blue
        // 2->2.5 interpolate between dark blue, cyan
        // 2.5->3 interpolate between cyan, green
        // 3->3.5 interpolate between green, yellow
        // 3.5->4 interpolate between yellow, orange
        // 4->5 interpolate between orange, brick
        // 5->6 interpolate between brick, salmon
        double mins[7] = { 1, 2, 2.5, 3, 3.5, 4, 5 };
        double maxs[7] = { 2, 2.5, 3, 3.5, 4, 5, 6 };
        unsigned char RGB[8][3] = { { 71, 71, 219 }, 
                                    { 0, 0, 91 },
                                    { 0, 255, 255 },
                                    { 0, 128, 0 },
                                    { 255, 255, 0 },
                                    { 255, 96, 0 },
                                    { 107, 0, 0 },
                                    { 224, 76, 76 } 
                                  };
        for (int j = 0 ; j < 3 ; j++)
        {
            float val = color_ptr[ptIds[j]];
            int r;
            for (r = 0 ; r < 7 ; r++)
            {
                if (mins[r] <= val && val < maxs[r])
                    break;
            }
            if (r == 7)
            {
                cerr << "Could not interpolate color for " << val << endl;
                exit(EXIT_FAILURE);
            }
            double proportion = (val-mins[r]) / (maxs[r]-mins[r]);
            tris[idx].colors[j][0] = (RGB[r][0] + proportion * (RGB[r+1][0] - RGB[r][0])) / 255.0;
            tris[idx].colors[j][1] = (RGB[r][1] + proportion * (RGB[r+1][1] - RGB[r][1])) / 255.0;
            tris[idx].colors[j][2] = (RGB[r][2] + proportion * (RGB[r+1][2] - RGB[r][2])) / 255.0;
        }
    }

    return tris;
}

int main()
{
    int imageWidth = 1000;
    int imageHeight = 1000;
    vtkImageData *image = NewImage(imageWidth, imageHeight);
    unsigned char *buffer = (unsigned char *) image->GetScalarPointer(0,0,0);
    int npixels = imageWidth*imageHeight;
    for (int i = 0 ; i < npixels*3 ; i++) 
        buffer[i] = 0;

    std::vector<Triangle> triangles = GetTriangles();

    Screen screen;
    screen.buffer = buffer;
    screen.width = imageWidth;
    screen.height = imageHeight;
    screen.setZBuffer();

    for (int i = 0; i < triangles.size(); i++) //triangles.size(); i++)
    {
        Triangle t = triangles[i];
        t.screen = screen;
        if (i == 2)
            t.runScanline(true);
        else
            t.runScanline(false);
    }

    WriteImage(image, "allTriangles");
}
