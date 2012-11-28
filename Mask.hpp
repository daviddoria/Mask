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

#ifndef Mask_HPP
#define Mask_HPP

#include "Mask.h" // Appease syntax parser

// ITK
#include "itkImageRegionIterator.h"

// Submodules
#include <ITKHelpers/ITKHelpers.h>

template<typename TImage, typename TColor>
void Mask::ApplyToRGBImage(TImage* const image, const TColor& color) const
{
  // Using generics, we allow any Color class that has .red(), .green(), and .blue() member functions
  // to be used to specify the color.
  if(image->GetLargestPossibleRegion() != this->GetLargestPossibleRegion())
  {
    std::cerr << "Image and mask must be the same size!" << std::endl
              << "Image region: " << image->GetLargestPossibleRegion() << std::endl
              << "Mask region: " << this->GetLargestPossibleRegion() << std::endl;
    return;
  }

  // Color the hole pixels in the image.
  typename TImage::PixelType holeValue;
  holeValue.SetSize(image->GetNumberOfComponentsPerPixel());
  holeValue.Fill(0);
  if(image->GetNumberOfComponentsPerPixel() >= 3)
  {
    holeValue[0] = color.red();
    holeValue[1] = color.green();
    holeValue[2] = color.blue();
  }

  itk::ImageRegionConstIterator<Mask> maskIterator(this, this->GetLargestPossibleRegion());

  while(!maskIterator.IsAtEnd())
  {
    if(this->IsHole(maskIterator.GetIndex()))
    {
      image->SetPixel(maskIterator.GetIndex(), holeValue);
    }

    ++maskIterator;
  }
}

template<typename TImage>
void Mask::ApplyToImage(TImage* const image, const typename TImage::PixelType& color) const
{
  if(image->GetLargestPossibleRegion() != this->GetLargestPossibleRegion())
  {
    std::cerr << "Image and mask must be the same size!" << std::endl
              << "Image region: " << image->GetLargestPossibleRegion() << std::endl
              << "Mask region: " << this->GetLargestPossibleRegion() << std::endl;
    return;
  }
  ApplyRegionToImageRegion(this->GetLargestPossibleRegion(), image, image->GetLargestPossibleRegion(), color);
}

template<typename TImage>
void Mask::ApplyRegionToImageRegion(const itk::ImageRegion<2>& maskRegion, TImage* const image,
                                    const itk::ImageRegion<2>& imageRegion, const typename TImage::PixelType& color) const
{
  if(maskRegion.GetSize() != imageRegion.GetSize())
    {
    std::cerr << "imageRegion and maskRegion must be the same size!" << std::endl
              << "Image region: " << imageRegion << std::endl
              << "Mask region: " << maskRegion << std::endl;
    return;
    }

  itk::ImageRegionConstIterator<Mask> maskIterator(this, maskRegion);
  itk::ImageRegionIterator<TImage> imageIterator(image, imageRegion);

  while(!maskIterator.IsAtEnd())
  {
    if(this->IsHole(maskIterator.GetIndex()))
    {
      imageIterator.Set(color);
    }

    ++maskIterator;
    ++imageIterator;
  }
}

template<typename TImage>
void Mask::ApplyToScalarImage(TImage* const image, typename TImage::PixelType holeValue) const
{
  if(image->GetLargestPossibleRegion() != this->GetLargestPossibleRegion())
  {
    std::cerr << "Image and mask must be the same size!" << std::endl
              << "Image region: " << image->GetLargestPossibleRegion() << std::endl
              << "Mask region: " << this->GetLargestPossibleRegion() << std::endl;
    return;
  }

  itk::ImageRegionConstIterator<Mask> maskIterator(this, this->GetLargestPossibleRegion());

  while(!maskIterator.IsAtEnd())
  {
    if(this->IsHole(maskIterator.GetIndex()))
    {
      image->SetPixel(maskIterator.GetIndex(), holeValue);
    }

    ++maskIterator;
  }
}

template<typename TImage>
void Mask::CreateFromImage(const TImage* const image, const HolePixelValueWrapper<typename TImage::PixelType>& holeValue,
                           const ValidPixelValueWrapper<typename TImage::PixelType>& validValue)
{
  this->SetRegions(image->GetLargestPossibleRegion());
  this->Allocate();

  itk::ImageRegionConstIterator<TImage> imageIterator(image, image->GetLargestPossibleRegion());

  unsigned int holeCounter = 0;
  unsigned int validCounter = 0;
  unsigned int undeterminedCounter = 0;
  while(!imageIterator.IsAtEnd())
  {
    typename TImage::PixelType currentPixel = imageIterator.Get();
    //std::cout << "Current color: " << currentPixel << std::endl;
    if(currentPixel == holeValue)
    {
      this->SetPixel(imageIterator.GetIndex(), HoleMaskPixelTypeEnum::HOLE);
      holeCounter++;
    }
    else if(currentPixel == validValue)
    {
      this->SetPixel(imageIterator.GetIndex(), HoleMaskPixelTypeEnum::VALID);
      validCounter++;
    }
    else
    {
      this->SetPixel(imageIterator.GetIndex(), HoleMaskPixelTypeEnum::UNDETERMINED);
      undeterminedCounter++;
    }

    ++imageIterator;
  }
  std::cout << "Mask::CreateFromImage: There were "
               << holeCounter << " hole pixels." << std::endl
               << validCounter << " valid pixels." << std::endl
               << undeterminedCounter << " undetermined pixels." << std::endl;
}

template <typename TImage>
void Mask::CreateHolesFromValue(const TImage* const inputImage,
                                const typename TImage::PixelType value)
{
  assert(inputImage->GetLargestPossibleRegion() == this->GetLargestPossibleRegion());

  itk::ImageRegionConstIterator<TImage> inputIterator(inputImage, inputImage->GetLargestPossibleRegion());
  itk::ImageRegionIterator<Mask> thisIterator(this, this->GetLargestPossibleRegion());

  while(!inputIterator.IsAtEnd())
  {
    if(inputIterator.Get() == value)
    {
      thisIterator.Set(HoleMaskPixelTypeEnum::HOLE);
    }
    ++inputIterator;
    ++thisIterator;
  }
}

template <typename TImage>
void Mask::CreateValidPixelsFromValue(const TImage* const inputImage,
                                      const typename TImage::PixelType value)
{
  assert(inputImage->GetLargestPossibleRegion() == this->GetLargestPossibleRegion());

  itk::ImageRegionConstIterator<TImage> inputIterator(inputImage, inputImage->GetLargestPossibleRegion());
  itk::ImageRegionIterator<Mask> thisIterator(this, this->GetLargestPossibleRegion());

  while(!inputIterator.IsAtEnd())
  {
    if(inputIterator.Get() == value)
    {
      thisIterator.Set(HoleMaskPixelTypeEnum::VALID);
    }
    ++inputIterator;
    ++thisIterator;
  }
}

template <typename TPixel>
void Mask::ReadFromImage(const std::string& filename,
                         const HolePixelValueWrapper<TPixel>& holeValue,
                         const ValidPixelValueWrapper<TPixel>& validValue)
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

  this->SetRegions(imageReader->GetOutput()->GetLargestPossibleRegion());
  this->Allocate();

  CreateHolesFromValue(imageReader->GetOutput(),
                       static_cast<ReadPixelType>(holeValue.Value));
  CreateValidPixelsFromValue(imageReader->GetOutput(),
                             static_cast<ReadPixelType>(validValue.Value));
}

#endif // Mask_HPP
