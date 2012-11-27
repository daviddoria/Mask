/*=========================================================================
 *
 *  Copyright David Doria 2012 daviddoria@gmail.com
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#include "MaskVTK.h"

// VTK
#include <vtkImageData.h>

void SetMaskTransparency(const Mask* const input, vtkImageData* outputImage)
{
  assert(input);

  // Setup and allocate the VTK image
  outputImage->SetDimensions(input->GetLargestPossibleRegion().GetSize()[0],
                             input->GetLargestPossibleRegion().GetSize()[1],
                             1);

  outputImage->AllocateScalars(VTK_UNSIGNED_CHAR, 4);

  // Copy all of the pixels to the output
  itk::ImageRegionConstIteratorWithIndex<Mask> imageIterator(input, input->GetLargestPossibleRegion());
  imageIterator.GoToBegin();

  while(!imageIterator.IsAtEnd())
  {
    unsigned char* pixel = static_cast<unsigned char*>(outputImage->GetScalarPointer(imageIterator.GetIndex()[0],
                                                       imageIterator.GetIndex()[1],0));
    /*
    // Set masked pixels to bright green and opaque. Set non-masked pixels to black and fully transparent.
    pixel[0] = 0;
    pixel[1] = imageIterator.Get();
    pixel[2] = 0;
    */

    // Set masked pixels to bright red and opaque. Set non-masked pixels to black and fully transparent.
    pixel[0] = 255;
    pixel[1] = 0;
    pixel[2] = 0;

    if(input->IsHole(imageIterator.GetIndex()))
    {
      pixel[3] = 255;
    }
    else
    {
      pixel[3] = 0;
    }

    ++imageIterator;
  }

  outputImage->Modified();
}
