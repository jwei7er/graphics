/*=========================================================================

 Author: Jordan Weiler
 Date:   May 7, 2013

 This program produces 1000 images from rotating the camera around a 3D
 isosurface. Triangles are loaded from a VTK file. The triangles are 
 transformed to device space. Lighting and Phong shading are applied. 

=========================================================================*/

#include <iostream>
#include <sstream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkDataSetWriter.h>

using std::cerr;
using std::endl;

double ceil441(double f) {
    return ceil(f-0.00001);
}

double floor441(double f) {
    return floor(f+0.00001);
}
	
double max(double a, double b) {
	return (b > a) ? b: a;
}

double min(double a, double b) {
	return (a > b) ? b: a;
}

double max(double a, double b, double c) {
	return (b > a) ? max(b, c) : max(a, c);
}

double min(double a, double b, double c) {
	return (a > b) ? min(b, c) : min(a, c);
}

double dotProduct(double array1[], double array2[], int length)
{
	//cout << "cal dot prod\n";
	double product = 0;
	for (int i=0; i < length; i++) {
		//cout << "a1: " << array1[i] << ", a2: " << array2[i];
		product += (array1[i] * array2[i]);
		//cout << " = " << product << "\n";
	}
	//cout << "total product: " << product << "\n";
	return product;
}

class Screen
{
	public:
		unsigned char *buffer;
		double *zBuffer;
		int width, height;

	// would some methods for accessing and setting pixels be helpful?
	void setZBuffer()
	{
		zBuffer = new double[width*height];
		
		for (int i = 0; i < width * height; i++)
		{
			zBuffer[i] = -1.0;
		}
	}
};

std::vector<double> crossProduct(double *a, double *b)
{
	std::vector<double> product(3);
	product[0] = (a[1] * b[2]) - (a[2] * b[1]);
	product[1] = (a[2] * b[0]) - (a[0] * b[2]);
	product[2] = (a[0] * b[1]) - (a[1] * b[0]);
	
	return product;
}

struct LightingParameters
{
    LightingParameters(void)
    {
         lightDir[0] = +0.6;
         lightDir[1] = 0;
         lightDir[2] = +0.8;
         Ka = 0.3;
         Kd = 0.7;
         Ks = 5.3;
         gamma = 7.5;
		 S = 1.0;
    };

    double lightDir[3];  // The direction of the light source
    double Ka;           // The coefficient for ambient lighting.
    double Kd;           // The coefficient for diffuse lighting.
    double Ks;           // The coefficient for specular lighting.
    double gamma;        // The exponent term for specular lighting.
    double S;		 // Another coefficient used with specular lighting.
};

class Matrix
{
	public:
		double A[4][4];

		void TransformPoint(const double *ptIn, double *ptOut);
		static Matrix ComposeMatrices(const Matrix &, const Matrix &);
		void Print(ostream &o);
};

void Matrix::Print(ostream &o)
{
    for (int i = 0 ; i < 4 ; i++)
    {
        char str[256];
        sprintf(str, "(%.7f %.7f %.7f %.7f)\n", A[i][0], A[i][1], A[i][2], A[i][3]);
        o << str;
    }
}

Matrix Matrix::ComposeMatrices(const Matrix &M1, const Matrix &M2)
{
    Matrix rv;
    for (int i = 0 ; i < 4 ; i++)
        for (int j = 0 ; j < 4 ; j++)
        {
            rv.A[i][j] = 0;
            for (int k = 0 ; k < 4 ; k++)
                rv.A[i][j] += M1.A[i][k]*M2.A[k][j];
        }

    return rv;
}

void Matrix::TransformPoint(const double *ptIn, double *ptOut)
{
    ptOut[0] = ptIn[0]*A[0][0]
             + ptIn[1]*A[1][0]
             + ptIn[2]*A[2][0]
             + ptIn[3]*A[3][0];
    ptOut[1] = ptIn[0]*A[0][1]
             + ptIn[1]*A[1][1]
             + ptIn[2]*A[2][1]
             + ptIn[3]*A[3][1];
    ptOut[2] = ptIn[0]*A[0][2]
             + ptIn[1]*A[1][2]
             + ptIn[2]*A[2][2]
             + ptIn[3]*A[3][2];
    ptOut[3] = ptIn[0]*A[0][3]
             + ptIn[1]*A[1][3]
             + ptIn[2]*A[2][3]
             + ptIn[3]*A[3][3];
}

class Camera
{
	public:
		double near, far;
		double angle;
		double position[3];
		double focus[3];
		double up[3];

		double O[3];
		double v1[3];
		double v2[3];
		double v3[3];
			
		Matrix CameraTransform;
		Matrix ViewTransform;
		Matrix DeviceTransform;
	
		void printCamera()
		{
			cout << "position: (" << position[0] << " " << position[1] << " " << position[2] << ")\n";
			cout << "focus: (" << focus[0] << " " << focus[1] << " " << focus[2] << ")\n";
			cout << "up: (" << up[0] << " " << up[1] << " " << up[2] << ")\n\n";
		}
		
		void calculateCameraTransform()
		{
			// O
			for (int i = 0; i < 3; i++) {
				O[i] = position[i];
			}
			//cout << "O: " << O[0] << " " << O[1] << " " << O[2] << "\n";
			
			double OmF[3];
			for (int i = 0; i < 3; i++) {
				OmF[i] = O[i] - focus[i];
			}
			//cout << "OmF: " << OmF[0] << " " << OmF[1] << " " << OmF[2] << "\n";
			
			// V1
			std::vector<double> vect1 = crossProduct(up, OmF);
			for (int i = 0; i < 3; i++) {
				v1[i] = vect1[i];
			}
			double total = 0.0;
			for (int i = 0; i < 3; i++) {
				double temp = v1[i];
				total += temp * temp;
			}
			double v1Sq = total;
			for (int i = 0; i < 3; i++) {
				if (fabs(v1Sq) < 0.00001)
					v1[i] = 0;
				else
					v1[i] = v1[i] / sqrt(v1Sq);
			}
			//cout << "v1 val: " << v1[0] << " " << v1[1] << " " << v1[2] << "\n\n";
			
			// V2
			std::vector<double> vect2 = crossProduct(OmF, v1);
			for (int i = 0; i < 3; i++) {
				v2[i] = vect2[i];
			}
			total = 0.0;
			for (int i = 0; i < 3; i++) {
				double temp = v2[i];
				total += temp * temp;
			}
			double v2Sq = total;
			for (int i = 0; i < 3; i++) {
				if (fabs(v2Sq) < 0.00001)
					v2[i] = 0;
				else
					v2[i] = v2[i] / sqrt(v2Sq);
			}
			//cout << "v2 val: " << v2[0] << " " << v2[1] << " " << v2[2] << "\n\n";
			
			// V3
			for (int i = 0; i < 3; i++) {
				v3[i] = OmF[i];
			}
			total = 0.0;
			for (int i = 0; i < 3; i++) {
				double temp = v3[i];
				total += temp * temp;
			}
			double v3Sq = total;
			for (int i = 0; i < 3; i++) {
				if (fabs(v3Sq) < 0.00001) {
					v3[i] = 0;
				}
				else {
					v3[i] = v3[i] / sqrt(v3Sq);
				}
			}
			//cout << "v3 val: " << v3[0] << " " << v3[1] << " " << v3[2] << "\n\n";
			
			double t[3];
			t[0] = 0 - O[0];
			t[1] = 0 - O[1];
			t[2] = 0 - O[2];
			
			CameraTransform.A[0][0] = v1[0];
			CameraTransform.A[0][1] = v2[0];
			CameraTransform.A[0][2] = v3[0]; 
			CameraTransform.A[0][3] = 0;
			CameraTransform.A[1][0] = v1[1];
			CameraTransform.A[1][1] = v2[1];
			CameraTransform.A[1][2] = v3[1];
			CameraTransform.A[1][3] = 0;
			CameraTransform.A[2][0] = v1[2];
			CameraTransform.A[2][1] = v2[2];
			CameraTransform.A[2][2] = v3[2];
			CameraTransform.A[2][3] = 0;
			CameraTransform.A[3][0] = dotProduct(v1, t, 3);
			CameraTransform.A[3][1] = dotProduct(v2, t, 3);
			CameraTransform.A[3][2] = dotProduct(v3, t, 3);
			CameraTransform.A[3][3] = 1;
			
			//cout << "printing Camera Transform \n";
			//CameraTransform.Print(cout);
		}
		
		void calculateViewTransform()
		{
			ViewTransform.A[0][0] = 1 / tan(angle/2);
			ViewTransform.A[0][1] = 0;
			ViewTransform.A[0][2] = 0;
			ViewTransform.A[0][3] = 0;
			ViewTransform.A[1][0] = 0;
			ViewTransform.A[1][1] = 1 / tan(angle/2);
			ViewTransform.A[1][2] = 0;
			ViewTransform.A[1][3] = 0;
			ViewTransform.A[2][0] = 0;
			ViewTransform.A[2][1] = 0;
			ViewTransform.A[2][2] = (far + near) / (far - near);
			ViewTransform.A[2][3] = -1;
			ViewTransform.A[3][0] = 0;
			ViewTransform.A[3][1] = 0;
			ViewTransform.A[3][2] = (2 * far * near) / (far - near);
			ViewTransform.A[3][3] = 0;
			
			//cout << "printing View Transform \n";
			//ViewTransform.Print(cout);
		}
		
		void calculateDeviceTransform(Screen s)
		{
			DeviceTransform.A[0][0] = s.width / 2;
			DeviceTransform.A[0][1] = 0;
			DeviceTransform.A[0][2] = 0;
			DeviceTransform.A[0][3] = 0;
			DeviceTransform.A[1][0] = 0;
			DeviceTransform.A[1][1] = s.width / 2;
			DeviceTransform.A[1][2] = 0;
			DeviceTransform.A[1][3] = 0;
			DeviceTransform.A[2][0] = 0;
			DeviceTransform.A[2][1] = 0;
			DeviceTransform.A[2][2] = 1;
			DeviceTransform.A[2][3] = 0;
			DeviceTransform.A[3][0] = s.width / 2;
			DeviceTransform.A[3][1] = s.width / 2;
			DeviceTransform.A[3][2] = 0;
			DeviceTransform.A[3][3] = 1;
			
			//cout << "printing Device Transform \n";
			//DeviceTransform.Print(cout);
		}
};

double SineParameterize(int curFrame, int nFrames, int ramp)
{
    int nNonRamp = nFrames-2*ramp;
    double height = 1./(nNonRamp + 4*ramp/M_PI);
    if (curFrame < ramp)
    {
        double factor = 2*height*ramp/M_PI;
        double eval = cos(M_PI/2*((double)curFrame)/ramp);
        return (1.-eval)*factor;
    }
    else if (curFrame > nFrames-ramp)
    {
        int amount_left = nFrames-curFrame;
        double factor = 2*height*ramp/M_PI;
        double eval = cos(M_PI/2*((double)amount_left/ramp));
        return 1. - (1-eval)*factor;
    }
    double amount_in_quad = ((double)curFrame-ramp);
    double quad_part = amount_in_quad*height;
    double curve_part = height*(2*ramp)/M_PI;
    return quad_part+curve_part;
}

Camera GetCamera(int frame, int nframes)
{
    double t = SineParameterize(frame, nframes, nframes/10);
    Camera c;
    c.near = 5;
    c.far = 200;
    c.angle = M_PI/6;
    c.position[0] = 40*sin(2*M_PI*t);
    c.position[1] = 40*cos(2*M_PI*t);
    c.position[2] = 40;
    c.focus[0] = 0;
    c.focus[1] = 0;
    c.focus[2] = 0;
    c.up[0] = 0;
    c.up[1] = 1;
    c.up[2] = 0;
    return c;
}

double calculatePhongShading(LightingParameters &lp, double *viewDirection, double *normal)
{
	double shading_amount = 0, diffuse = 0, specular = 0;
	
	diffuse = fabs(dotProduct(lp.lightDir, normal, 3));
	double rTemp = 2 * dotProduct(lp.lightDir, normal, 3);
	double R[] = {rTemp * normal[0] - lp.lightDir[0], rTemp * normal[1] - lp.lightDir[1], rTemp * normal[2] - lp.lightDir[2]};
	specular = max(0, lp.S * pow(dotProduct(R, viewDirection, 3), lp.gamma));
	shading_amount = lp.Ka + lp.Kd * diffuse + lp.Ks * specular;
	
	return shading_amount;
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

class TriangleSide
{
	public:
		double x1, y1, z1, rgb1[3], normal1[3];
		double x2, y2, z2, rgb2[3], normal2[3];
		double slope;
		double yIntercept;
		bool vertical, ignore;
		double minY, maxY;
		
		TriangleSide() {
			// default constructor
		}
		
		TriangleSide(double nX1, double nY1, double nZ1, double nRGB1[], double nNormal1[], double nX2, double nY2, double nZ2, double nRGB2[], double nNormal2[])
		{
			x1 = nX1;
			y1 = nY1;
			z1 = nZ1;
			x2 = nX2;
			y2 = nY2;
			z2 = nZ2;
			for (int i = 0; i < 3; i++) {
				rgb1[i] = nRGB1[i];
				normal1[i] = nNormal1[i];
				rgb2[i] = nRGB2[i];
				normal2[i] = nNormal2[i];
			}
			
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
		
		double getXValueGivenY(double y) {
			if (vertical == true)
				return slope;
			else {
				if (slope == 0)
					return 0;
				else
					return (y - yIntercept) / slope;
			}
		}

		double interpolation(double v1, double v2, double A, double B, double v) {
			return A + ((v - v1) / (v2 - v1)) * (B - A);
		}
		
		double getZValueGivenY(double y) {
			return interpolation(y1, y2, z1, z2, y);
		}
		
		double getRValueGivenY(double y) {
			return interpolation(y1, y2, rgb1[0], rgb2[0], y);
		}
		
		double getGValueGivenY(double y) {
			return interpolation(y1, y2, rgb1[1], rgb2[1], y);
		}
		
		double getBValueGivenY(double y) {
			return interpolation(y1, y2, rgb1[2], rgb2[2], y);
		}
		
		double getNormXValueGivenY(double y) {
			return interpolation(y1, y2, normal1[0], normal2[0], y);
		}
		
		double getNormYValueGivenY(double y) {
			return interpolation(y1, y2, normal1[1], normal2[1], y);
		}
		
		double getNormZValueGivenY(double y) {
			return interpolation(y1, y2, normal1[2], normal2[2], y);
		}
		
		bool withinBounds(int y) {
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
		double normals[3][3];
		double imageView[3]; // view vector for image
		Screen screen;
		LightingParameters lp;

	void print()
	{
		cout << "Triangle\nX: (" << X[0] << "," << X[1] << "," << X[2] << ")\n";
		cout << "Y: (" << Y[0] << "," << Y[1] << "," << Y[2] << ")\n";
		cout << "Z: (" << Z[0] << "," << Z[1] << "," << Z[2] << ")\n";
	}
	
	void runScanline(bool debug)
	{
		double minX = max(ceil441(min(X[0], X[1], X[2])),0);
		double maxX = min(floor441(max(X[0], X[1], X[2])),screen.width-1);
		
		double minY = max(ceil441(min(Y[0], Y[1], Y[2])),0);
		double maxY = min(floor441(max(Y[0], Y[1], Y[2])),screen.height-1);
		
        if (debug) {
            cout << "X,Y,Z (" << X[0] << "," << Y[0] << "," << Z[0] << ") (" << X[1] << "," << Y[1] << "," << Z[1] << ") (" << X[2] << "," << Y[2] << "," << Z[2] << ")\n";
            cout << "RGB (" << colors[0][0] << "," << colors[0][1] << "," << colors[0][2] << ") ";
            cout << "(" << colors[1][0] << "," << colors[1][1] << "," << colors[1][2] << ") ";
            cout << "(" << colors[2][0] << "," << colors[2][1] << "," << colors[2][2] << ")\n";
			cout << "Normals (" << normals[0][0] << "," << normals[0][1] << "," << normals[0][2] << ") ";
			cout << "(" << normals[1][0] << "," << normals[1][1] << "," << normals[1][2] << ") ";
			cout << "(" << normals[2][0] << "," << normals[2][1] << "," << normals[2][2] << ")\n";
        }
		
		TriangleSide ts1 = TriangleSide(X[0], Y[0], Z[0], colors[0], normals[0], X[1], Y[1], Z[1], colors[1], normals[1]);
		//cout << "side 1: y = " << ts1.slope << "x + " << ts1.yIntercept << "\n";
		TriangleSide ts2 = TriangleSide(X[1], Y[1], Z[1], colors[1], normals[1], X[2], Y[2], Z[2], colors[2], normals[2]);
		//cout << "side 2: y = " << ts2.slope << "x + " << ts2.yIntercept << "\n";
		TriangleSide ts3 = TriangleSide(X[2], Y[2], Z[2], colors[2], normals[2], X[0], Y[0], Z[0], colors[0], normals[0]);
		//cout << "side 3: y = " << ts3.slope << "x + " << ts3.yIntercept << "\n";
		
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
				//cout << "ts1: " << tempX << " ";
				minXGivenY = ts1.getXValueGivenY(y);
				maxXGivenY = minXGivenY;
				tsLeft = ts1;
				tsRight = ts1;
			}
			if (!ts2.ignore and ts2.withinBounds(y)) {
				tempX = ts2.getXValueGivenY(y);
				//cout << "ts2: " << tempX << " ";
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
				//cout << "ts3: " << tempX << " ";
				if (tempX < minXGivenY) {
					minXGivenY = tempX;
					tsLeft = ts3;
				}
				if (tempX > maxXGivenY) {
					maxXGivenY = tempX;
					tsRight = ts3;
				}
			}
			
			double leftZ = tsLeft.getZValueGivenY(y);
			double leftRGB[] = {tsLeft.getRValueGivenY(y), tsLeft.getGValueGivenY(y), tsLeft.getBValueGivenY(y)};
			double leftNormal[] = {tsLeft.getNormXValueGivenY(y), tsLeft.getNormYValueGivenY(y), tsLeft.getNormZValueGivenY(y)};
			
			double rightZ = tsRight.getZValueGivenY(y);
			double rightRGB[] = {tsRight.getRValueGivenY(y), tsRight.getGValueGivenY(y), tsRight.getBValueGivenY(y)};
			double rightNormal[] = {tsRight.getNormXValueGivenY(y), tsRight.getNormYValueGivenY(y), tsRight.getNormZValueGivenY(y)};
			
			if (debug) {
				cout << "\nScanline " << y << " left: " << minXGivenY << ", depth: " << leftZ;
				cout << ", color: " << leftRGB[0] << "/" << leftRGB[1] << "/" << leftRGB[2];
				cout << ", normal: " << leftNormal[0] << "/" << leftNormal[1] << "/" << leftNormal[2];
				cout << " right: " << maxXGivenY << ", depth: " << rightZ;
				cout << ", color: " << rightRGB[0] << "/" << rightRGB[1] << "/" << rightRGB[2];
				cout << ", normal: " << rightNormal[0] << "/" << rightNormal[1] << "/" << rightNormal[2] << "\n";
			}
			
			int lineMinX = max(ceil441(minXGivenY), minX);
			int lineMaxX = min(floor441(maxXGivenY), maxX);
			for (int x = lineMinX; x <= lineMaxX; x++)
			{
				double proportion;
				if (minXGivenY != maxXGivenY)
					proportion = (((double) x) - minXGivenY) / (maxXGivenY - minXGivenY);
				else
					proportion = 1.0;
			
				double z = leftZ + proportion * (rightZ - leftZ);
				if (screen.zBuffer[zOffset + x] <= z)
				{
                	double pixelNormal[3];
					pixelNormal[0] = leftNormal[0] + proportion * (rightNormal[0] - leftNormal[0]);
					pixelNormal[1] = leftNormal[1] + proportion * (rightNormal[1] - leftNormal[1]);
					pixelNormal[2] = leftNormal[2] + proportion * (rightNormal[2] - leftNormal[2]);
					
					double pixelShading = calculatePhongShading(lp, imageView, pixelNormal);
					
					// weight each pixel color by the phong shading value
					double pixelRGB[3];
					pixelRGB[0] = min(1.0, (leftRGB[0] + proportion * (rightRGB[0] - leftRGB[0])) * pixelShading);
					pixelRGB[1] = min(1.0, (leftRGB[1] + proportion * (rightRGB[1] - leftRGB[1])) * pixelShading);
					pixelRGB[2] = min(1.0, (leftRGB[2] + proportion * (rightRGB[2] - leftRGB[2])) * pixelShading);
					
					if (debug) {
                        cout << "assigning pixel (" << x << "," << y << ")";
                        cout << "R:" << ceil441(255.0 * pixelRGB[0]) << " G:" << ceil441(255.0 * pixelRGB[1]) << " B:" << ceil441(255.0 * pixelRGB[2]) << " Z:" << z << "\n";
                    }
					
					screen.zBuffer[zOffset + x] = z;
					screen.buffer[offset + (x * 3)] = (unsigned char) ceil441(255.0 * pixelRGB[0]);
					screen.buffer[offset + (x * 3) + 1] = (unsigned char) ceil441(255.0 * pixelRGB[1]);
					screen.buffer[offset + (x * 3) + 2] = (unsigned char) ceil441(255.0 * pixelRGB[2]);
				}
			}
		}
		
		if (debug) {
			cout << "\n\n";		
		}
	}
};

std::vector<Triangle> GetTriangles(void)
{
    vtkPolyDataReader *rdr = vtkPolyDataReader::New();
    rdr->SetFileName("proj1f_geometry.vtk");
    cerr << "Reading" << endl;
    rdr->Update();
    cerr << "Done reading" << endl;
    if (rdr->GetOutput()->GetNumberOfCells() == 0)
    {
        cerr << "Unable to open file!!" << endl;
        exit(EXIT_FAILURE);
    }
    vtkPolyData *pd = rdr->GetOutput();
/*
vtkDataSetWriter *writer = vtkDataSetWriter::New();
writer->SetInput(pd);
writer->SetFileName("hrc.vtk");
writer->Write();
 */

    int numTris = pd->GetNumberOfCells();
    vtkPoints *pts = pd->GetPoints();
    vtkCellArray *cells = pd->GetPolys();
    vtkDoubleArray *var = (vtkDoubleArray *) pd->GetPointData()->GetArray("hardyglobal");
    double *color_ptr = var->GetPointer(0);
    //vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray("hardyglobal");
    //float *color_ptr = var->GetPointer(0);
    vtkFloatArray *n = (vtkFloatArray *) pd->GetPointData()->GetNormals();
    float *normals = n->GetPointer(0);
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
        double *pt = NULL;
        pt = pts->GetPoint(ptIds[0]);
        tris[idx].X[0] = pt[0];
        tris[idx].Y[0] = pt[1];
        tris[idx].Z[0] = pt[2];
        tris[idx].normals[0][0] = normals[3*ptIds[0]+0];
        tris[idx].normals[0][1] = normals[3*ptIds[0]+1];
        tris[idx].normals[0][2] = normals[3*ptIds[0]+2];
        pt = pts->GetPoint(ptIds[1]);
        tris[idx].X[1] = pt[0];
        tris[idx].Y[1] = pt[1];
        tris[idx].Z[1] = pt[2];
        tris[idx].normals[1][0] = normals[3*ptIds[1]+0];
        tris[idx].normals[1][1] = normals[3*ptIds[1]+1];
        tris[idx].normals[1][2] = normals[3*ptIds[1]+2];
        pt = pts->GetPoint(ptIds[2]);
        tris[idx].X[2] = pt[0];
        tris[idx].Y[2] = pt[1];
        tris[idx].Z[2] = pt[2];
        tris[idx].normals[2][0] = normals[3*ptIds[2]+0];
        tris[idx].normals[2][1] = normals[3*ptIds[2]+1];
        tris[idx].normals[2][2] = normals[3*ptIds[2]+2];

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
            tris[idx].colors[j][0] = (RGB[r][0]+proportion*(RGB[r+1][0]-RGB[r][0]))/255.0;
            tris[idx].colors[j][1] = (RGB[r][1]+proportion*(RGB[r+1][1]-RGB[r][1]))/255.0;
            tris[idx].colors[j][2] = (RGB[r][2]+proportion*(RGB[r+1][2]-RGB[r][2]))/255.0;
        }
    }

    return tris;
}

void SaveImage(vtkImageData *image, int count)
{
	std::string iName;
	std::ostringstream convert;
	convert << "image" << count;
	iName = convert.str();
	const char * c = iName.c_str();
	WriteImage(image, c);
}

Screen InitializeScreen(unsigned char *buffer, int imageWidth, int imageHeight)
{
	Screen screen;
	screen.buffer = buffer;
	screen.width = imageWidth;
	screen.height = imageHeight;
	screen.setZBuffer();
	
	// initialize the buffer to black
	int npixels = imageWidth*imageHeight;
	for (int i = 0 ; i < npixels*3 ; i++) 
		screen.buffer[i] = 0;
	
	return screen;
}

Triangle transformTriangleToDeviceSpace(Triangle t, Matrix m, bool debug)
{
	if (debug) {
		t.print();
	}
	
	Triangle newT;
	for (int i = 0; i < 3; i++) 
	{
		double pointIn[4];
		pointIn[0] = t.X[i];
		pointIn[1] = t.Y[i];
		pointIn[2] = t.Z[i];
		pointIn[3] = 1;
		double pointOut[4];
		m.TransformPoint(pointIn, pointOut);
		newT.X[i] = pointOut[0] / pointOut[3];
		newT.Y[i] = pointOut[1] / pointOut[3];
		newT.Z[i] = pointOut[2] / pointOut[3];
		if (debug) {
			cout << "Transformed (" << t.X[i] << "," << t.Y[i] << "," << t.Z[i] << ")";
			cout << " to (" << newT.X[i] << "," << newT.Y[i] << "," << newT.Z[i] << ")\n";
		}
	}
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			newT.colors[i][j] = t.colors[i][j];
			newT.normals[i][j] = t.normals[i][j];
		}
	}
	
	if (debug) {t.print();};
	return newT;
}

int main()
{
	std::vector<Triangle> triangles = GetTriangles();
	int imageWidth = 1000;
	int imageHeight = 1000;
	vtkImageData *image = NewImage(imageWidth, imageHeight);
	unsigned char *buffer = 
		(unsigned char *) image->GetScalarPointer(0,0,0);
	Screen s;
	Camera c;

	for (int i = 0; i < 1000; i++)
	{
		s = InitializeScreen(buffer, imageWidth, imageHeight);
		c = GetCamera(i, 1000);
		c.calculateCameraTransform();
		c.calculateViewTransform();
		c.calculateDeviceTransform(s);
		
		Matrix m2 = Matrix::ComposeMatrices(Matrix::ComposeMatrices(c.CameraTransform, c.ViewTransform), c.DeviceTransform);
		
		for (int j = 0; j < triangles.size(); j++)
		{
			Triangle t = triangles[j];
			Triangle newT = transformTriangleToDeviceSpace(t, m2, false);
			for (int k = 0; k < 3; k++) {
				newT.imageView[k] = c.v3[k];
			}
			newT.screen = s;
			newT.runScanline(false);
		}

		SaveImage(image, i);
		
		delete [] s.zBuffer;
	}
}

