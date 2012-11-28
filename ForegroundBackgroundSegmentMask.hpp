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

#ifndef ForegroundBackgroundSegmentMask_HPP
#ifndef ForegroundBackgroundSegmentMask_HPP

#include "ForegroundBackgroundSegmentMask.h"

// Submodules
#include <ITKHelpers/ITKHelpers.h>


template <typename TPixel>
void Mask::ReadFromImage(const std::string& filename,
                         const ForegroundPixelValueWrapper<TPixel>& foregroundValue,
                         const BackgroundPixelValueWrapper<TPixel>& backgroundValue)
{
  std::cout << "Reading mask from image: " << filename << std::endl;

  // Ensure the input image can be interpreted as a mask.
  unsigned int numberOfComponents =
      ITKHelpers::GetNumberOfComponentsPerPixelInFile(filename);

  if(!(numberOfComponents == 1 || numberOfComponents == 3))
  {
    std::stringstream ss;
    ss << "Number of components for a mask must be 1 or 3! (" << filename
       << " is " << numberOfComponents << ")";
    throw std::runtime_error(ss.str());
  }

  // Read the image
  typedef int ReadPixelType;
  typedef itk::Image<ReadPixelType, 2> ImageType;
  typedef  itk::ImageFileReader<ImageType> ImageReaderType;
  ImageReaderType::Pointer imageReader = ImageReaderType::New();
  imageReader->SetFileName(filename);
  imageReader->Update();

  CreateHolesFromValue(imageReader->GetOutput(),
                       static_cast<ReadPixelType>(holeValue.Value));
  CreateValidPixelsFromValue(imageReader->GetOutput(),
                             static_cast<ReadPixelType>(validValue.Value));
}

#endif
