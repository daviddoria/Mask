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

#ifndef MaskVTK_H
#define MaskVTK_H

#include "Mask.h"

// Submodules
#include <ITKHelpers/ITKHelpers.h>

// VTK
class vtkImageData;

namespace MaskVTK
{
  ////////////////// Functions with VTK ////////////////
  void SetMaskTransparency(const Mask* const input, vtkImageData* outputImage);

  ////////////////// Function templates with VTK ////////////////
  template <typename TImage>
  void ITKImageToVTKImageMasked(const ITKHelpers::FloatVectorImageType* const image, const Mask* const mask,
                                vtkImageData* const outputImage, const unsigned char maskColor[3]);

  template <typename TPixel>
  void ITKImageToVTKImageMasked(const typename itk::VectorImage<TPixel, 2>* const image, const Mask* const mask,
                                vtkImageData* const outputImage, const unsigned char maskColor[3]);

}

#include "MaskVTK.hpp"

#endif
