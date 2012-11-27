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

#ifndef MaskVTK_HPP
#define MaskVTK_HPP

// STL
#include <stdexcept>

// Custom
#include "Mask.h"

// VTK
#include <vtkImageData.h>

namespace MaskVTK
{

  template <typename TImage>
  void ITKImageToVTKImageMasked(const TImage* const image, const Mask* const mask,
                                                            vtkImageData* const outputImage, const unsigned char maskColor[3])
  {
    assert(mask);
    // This function assumes an ND (with N>3) image has the first 3 channels as RGB and extra
    // information in the remaining channels.

    //std::cout << "ITKImagetoVTKRGBImage()" << std::endl;
    if(image->GetNumberOfComponentsPerPixel() < 3)
    {
      std::cerr << "The input image has " << image->GetNumberOfComponentsPerPixel()
                << " components, but at least 3 are required." << std::endl;
      return;
    }

    // Setup and allocate the image data
    //outputImage->SetNumberOfScalarComponents(3);
    //outputImage->SetScalarTypeToUnsignedChar();
    outputImage->SetDimensions(image->GetLargestPossibleRegion().GetSize()[0],
                                                           image->GetLargestPossibleRegion().GetSize()[1],
                                                           1);
    //outputImage->AllocateScalars();
    outputImage->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

    // Copy all of the input image pixels to the output image
    itk::ImageRegionConstIteratorWithIndex<TImage>
                   imageIterator(image,image->GetLargestPossibleRegion());
    imageIterator.GoToBegin();

    while(!imageIterator.IsAtEnd())
    {
      unsigned char* VTKPixel = static_cast<unsigned char*>(
            outputImage->GetScalarPointer(imageIterator.GetIndex()[0], imageIterator.GetIndex()[1],0));
      if(mask->IsValid(imageIterator.GetIndex()))
      {
        for(unsigned int component = 0; component < 3; component++)
        {
          VTKPixel[component] = static_cast<unsigned char>(imageIterator.Get()[component]);
        }
      }
      else
      {
        for(unsigned int component = 0; component < 3; component++)
        {
          VTKPixel[component] = maskColor[component];
        }
      }

      ++imageIterator;
    }

    outputImage->Modified();
  }

  template <typename TPixel>
  void ITKImageToVTKImageMasked(const typename itk::VectorImage<TPixel, 2>* const image, const Mask* const mask,
                                                            vtkImageData* const outputImage, const unsigned char maskColor[3])
  {
    assert(mask);

    typedef typename itk::VectorImage<TPixel, 2> VectorImageType;
    // This function assumes an ND (with N>3) image has the first 3 channels as RGB and extra
    // information in the remaining channels.

    //std::cout << "ITKImagetoVTKRGBImage()" << std::endl;
    if(image->GetNumberOfComponentsPerPixel() < 3)
    {
      std::cerr << "The input image has " << image->GetNumberOfComponentsPerPixel()
                << " components, but at least 3 are required." << std::endl;
      return;
    }

    // Setup and allocate the image data
    //outputImage->SetNumberOfScalarComponents(3);
    //outputImage->SetScalarTypeToUnsignedChar();
    outputImage->SetDimensions(image->GetLargestPossibleRegion().GetSize()[0],
                                                           image->GetLargestPossibleRegion().GetSize()[1],
                                                           1);
    //outputImage->AllocateScalars();
    outputImage->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

    // Copy all of the input image pixels to the output image
    itk::ImageRegionConstIteratorWithIndex<VectorImageType>
                   imageIterator(image,image->GetLargestPossibleRegion());
    imageIterator.GoToBegin();

    while(!imageIterator.IsAtEnd())
    {
      unsigned char* VTKPixel = static_cast<unsigned char*>(
            outputImage->GetScalarPointer(imageIterator.GetIndex()[0], imageIterator.GetIndex()[1],0));
      if(mask->IsValid(imageIterator.GetIndex()))
      {
        for(unsigned int component = 0; component < 3; component++)
        {
          VTKPixel[component] = static_cast<unsigned char>(imageIterator.Get()[component]);
        }
      }
      else
      {
        for(unsigned int component = 0; component < 3; component++)
        {
          VTKPixel[component] = maskColor[component];
        }
      }

      ++imageIterator;
    }

    outputImage->Modified();
  }

} // end namespace

#endif // MaskVTK
