/*=========================================================================

 Author: Jordan Weiler
 Date:   April 4, 2013

 This program makes an image with 27 horizontal stripes, each with 50 pixels.
 The color for the Xth strip should be:
 	X % 3 = 0 		-> B=0
 	X % 3 = 1 		-> B=126
 	X % 3 = 2 		-> B=255
 	(X/3) % 3 = 0 	-> G=0
 	(X/3) % 3 = 1 	-> G=128
 	(X/3) % 3 = 2 	-> G=255
 	X/9 = 0 		-> R=0
 	X/9 = 1 		-> R=128
 	X/9 = 2 		-> R=255

=========================================================================*/

#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>

using std::cerr;
using std::endl;

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


int main()
{
   std::cerr << "In main!" << endl;
   vtkImageData *image = NewImage(1024, 1350);
   image->Print(cerr);
   unsigned char *buffer = 
     (unsigned char *) image->GetScalarPointer(0,0,0);
   cerr << "buffer is " << (void *) buffer << endl;

   // 27 strips each with 50 pixels   
   int i, j, offset;

   for(i = 0; i < 27; i++)
   {
      offset = i * 1024 * 50 * 3;

      for(j = 0; j < 1024 * 50 * 3; j++)
      {
         if (j % 3 == 0) {
            // Red pixels
            if (i / 9 == 0) {
               buffer[offset + j] = 0;
            }
            else if (i / 9 == 1) {
               buffer[offset + j] = 128;
            }
            else {
               buffer[offset + j] = 255;
            }
         }
         else if (j % 3 == 1) {
            // Green pixels
            if ((i / 3) % 3 == 0) {
               buffer[offset + j] = 0;
            }
            else if ((i / 3) % 3 == 1) {
               buffer[offset + j] = 128;
            }
            else {
               buffer[offset + j] = 255;
            }
         }
         else {
            // Blue pixels
            if (i % 3 == 0) {
               buffer[offset + j] = 0;
            }
            else if (i % 3 == 1) {
               buffer[offset + j] = 128;
            }
            else {
               buffer[offset + j] = 255;
            }
         }
      }
   }
   
   WriteImage(image, "p1Image");
}
