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

#ifndef StrokeMask_HPP
#define StrokeMask_HPP

#include "StrokeMask.h"

// Submodules
#include <ITKHelpers/ITKHelpers.h>

template <typename TPixel>
void StrokeMask::
ReadFromImage(const std::string& filename,
              const TPixel& strokeValue)
{
  std::cout << "Reading stroke mask from image: " << filename << std::endl;

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

  this->SetRegions(imageReader->GetOutput()->GetLargestPossibleRegion());
  this->Allocate();

  itk::ImageRegionConstIteratorWithIndex<ImageType>
      imageIterator(imageReader->GetOutput(),
                    imageReader->GetOutput()->GetLargestPossibleRegion());
  while(!imageIterator.IsAtEnd())
  {
    if(imageIterator.Get() == strokeValue)
    {
      this->SetPixel(imageIterator.GetIndex(),
                     StrokeMaskPixelTypeEnum::STROKE);
    }
    else
    {
      this->SetPixel(imageIterator.GetIndex(),
                     StrokeMaskPixelTypeEnum::NOTSTROKE);
    }
    ++imageIterator;
  }

}

template <typename TPixel>
void StrokeMask::
Write(const std::string& filename,
      const TPixel& strokeValue)
{
  typedef itk::Image<TPixel, 2> ImageType;
  typename ImageType::Pointer image = ImageType::New();

  image->SetRegions(this->GetLargestPossibleRegion());
  image->Allocate();

  itk::ImageRegionConstIteratorWithIndex<StrokeMask>
      maskIterator(this,
                   this->GetLargestPossibleRegion());
  while(!maskIterator.IsAtEnd())
  {
    if(maskIterator.Get() == StrokeMaskPixelTypeEnum::STROKE)
    {
      image->SetPixel(maskIterator.GetIndex(),
                      strokeValue);
    }
    else
    {
      image->SetPixel(maskIterator.GetIndex(), itk::NumericTraits<TPixel>::Zero);
    }

    ++maskIterator;
  }

  typedef  itk::ImageFileWriter<ImageType> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(filename);
  writer->SetInput(image);
  writer->Update();

}

#endif
