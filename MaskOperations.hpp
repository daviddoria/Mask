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

// STL
#include <stdexcept>

// Custom
#include "Mask.h"
#include <ITKHelpers/ITKHelpers.h>

// ITK
#include "itkBresenhamLine.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkGaussianOperator.h"
#include "itkImageRegionIterator.h"
#include "itkLaplacianOperator.h"
#include "itkMedianImageFilter.h"
#include "itkSignedDanielssonDistanceMapImageFilter.h"
#include "itkSimpleFastMutexLock.h"

// VTK
#ifdef MaskUseVTK
  #include <vtkImageData.h>
#endif

namespace MaskOperations
{
  
template <class TImage>
void CopySelfPatchIntoHoleOfTargetRegion(TImage* const image, const Mask* const mask,
                                  const itk::ImageRegion<2>& sourceRegionInput,
                                  const itk::ImageRegion<2>& destinationRegionInput)
{
  CopySourcePatchIntoHoleOfTargetRegion(image, image, mask, sourceRegionInput, destinationRegionInput);
}

template <class TImage>
void CopyRegionIntoHolePortionOfTargetRegion(const TImage* const sourceImage, TImage* const targetImage,
                                           const Mask* const mask,
                                           itk::ImageRegion<2> sourceRegion,
                                           itk::ImageRegion<2> destinationRegion)
{
  itk::ImageRegion<2> fullImageRegion = sourceImage->GetLargestPossibleRegion();

  sourceRegion = ITKHelpers::CropRegionAtPosition(sourceRegion, fullImageRegion, destinationRegion);
  destinationRegion.Crop(fullImageRegion);

  itk::ImageRegionConstIterator<TImage> sourceIterator(sourceImage, sourceRegion);
  itk::ImageRegionIterator<TImage> destinationIterator(targetImage, destinationRegion);
  itk::ImageRegionConstIterator<Mask> maskIterator(mask, destinationRegion);

  while(!sourceIterator.IsAtEnd())
  {
    if(maskIterator.Get() == mask->GetHoleValue()) // we are in the target region
    {
      destinationIterator.Set(sourceIterator.Get());
    }
    ++sourceIterator;
    ++maskIterator;
    ++destinationIterator;
  }
}

template<typename TImage>
void CreatePatchImage(const TImage* const image, const itk::ImageRegion<2>& sourceRegion,
                      const itk::ImageRegion<2>& targetRegion,
                      const Mask* const mask, TImage* const result)
{
  // The input 'result' is expected to already be sized and initialized.

  itk::ImageRegionConstIterator<TImage> sourceRegionIterator(image, sourceRegion);
  itk::ImageRegionConstIterator<TImage> targetRegionIterator(image, targetRegion);

  itk::ImageRegionIterator<TImage> resultIterator(result, result->GetLargestPossibleRegion());

  while(!sourceRegionIterator.IsAtEnd())
    {

    if(mask->IsHole(targetRegionIterator.GetIndex()))
      {
      resultIterator.Set(sourceRegionIterator.Get());
      }
    else
      {
      resultIterator.Set(targetRegionIterator.Get());
      }

    ++sourceRegionIterator;
    ++targetRegionIterator;
    ++resultIterator;
    }
}

template<typename TImage>
void FindMaximumValueInMaskedRegion(const TImage* const image, const Mask* const mask,
                                    const itk::ImageRegion<2>& region, const Mask::PixelType maskValue,
                                    typename TImage::PixelType& maxValue)
{
  // Return the location of the highest pixel in 'image' out of the pixels with 'maskValue' in the 'mask'.
  // Return the value of that pixel by reference.

  std::vector<itk::Index<2> > maskedPixelLocations = ITKHelpers::GetPixelsWithValueInRegion(mask, region, maskValue);

  if(maskedPixelLocations.size() <= 0)
  {
    throw std::runtime_error("FindMaximumValueInMaskedRegion(): No boundary pixels!");
  }

  std::vector<typename TImage::PixelType> maskedPixels = ITKHelpers::GetPixelValues(image, maskedPixelLocations);

  // Initialize
  Helpers::MaxOfAllIndices(maskedPixels, maxValue);
}

template<typename TImage>
void FindMinimumValueInMaskedRegion(const TImage* const image, const Mask* const mask, const itk::ImageRegion<2>& region,
                                    const Mask::PixelType maskValue, typename TImage::PixelType& minValue)
{
  // Return the location of the lowest pixel in 'image' out of the pixels with 'maskValue' in the 'mask'.
  // Return the value of that pixel by reference.

  std::vector<itk::Index<2> > maskedPixelLocations = ITKHelpers::GetPixelsWithValueInRegion(mask, region, maskValue);

  if(maskedPixelLocations.size() <= 0)
  {
    throw std::runtime_error("FindMinimumValueInMaskedRegion(): No masked pixels!");
  }

  std::vector<typename TImage::PixelType> maskedPixels = ITKHelpers::GetPixelValues(image, maskedPixelLocations);

  // Initialize
  Helpers::MinOfAllIndices(maskedPixels, minValue);
}

template<typename TImage, typename TRegionIndicatorImage>
itk::Index<2> FindHighestValueInNonZeroRegion(const TImage* const image, float& maxValue)
{
  // Create a mask from the indicator image
  Mask::Pointer mask = Mask::New();
  mask->CreateFromImage(image, itk::NumericTraits<typename TRegionIndicatorImage::PixelType>::Zero);
  return FindHighestValueInMaskedRegion(image, maxValue, mask);
}

template<typename TImage>
std::vector<typename TImage::PixelType> GetValidPixelsInRegion(const TImage* const image, const Mask* const mask,
                                                               const itk::ImageRegion<2>& region)
{
  std::vector<typename TImage::PixelType> validPixels;
  
  itk::ImageRegionConstIteratorWithIndex<TImage> imageIterator(image, region);

  while(!imageIterator.IsAtEnd())
    {
    if(mask->IsValid(imageIterator.GetIndex()))
      {
      validPixels.push_back(imageIterator.Get());
      }
    ++imageIterator;
    }

  return validPixels;
}

template<typename TImage>
void AddNoiseInHole(TImage* const image, const Mask* const mask, const float noiseVariance)
{
  itk::ImageRegionIterator<TImage> imageIterator(image, image->GetLargestPossibleRegion());

  while(!imageIterator.IsAtEnd())
    {
    if(mask->IsHole(imageIterator.GetIndex()))
      {
      float randomNumber = drand48();
      //std::cout << "randomNumber: " << randomNumber << std::endl;
      float noise = (randomNumber - .5) * noiseVariance;
      //std::cout << "noise: " << noise << std::endl;
      imageIterator.Set(imageIterator.Get() + noise);
      }
    ++imageIterator;
    }
}

template<typename TImage>
void InterpolateHole(TImage* const image, const Mask* const mask)
{
  struct WeightedPixel
  {
    float Weight;
    itk::Index<2> Pixel;
    
    WeightedPixel(const itk::Index<2>& pixel, const float weight)
    {
      Weight = weight;
      Pixel = pixel;
    }
    
    bool operator<(const WeightedPixel &other) const
    {
      if(Weight < other.Weight)
      {
        return true;
      }
      else
      {
        return false;
      }
    }
  };

  // Compute the bounding box of the hole
  itk::ImageRegion<2> boundingBox = MaskOperations::ComputeHoleBoundingBox(mask);

  // Compute the distance function in the bounding box
  typedef itk::SignedDanielssonDistanceMapImageFilter<TImage, ITKHelpersTypes::FloatScalarImageType>
          SignedDanielssonDistanceMapImageFilterType;
  typename SignedDanielssonDistanceMapImageFilterType::Pointer distanceMapFilter =
           SignedDanielssonDistanceMapImageFilterType::New();
  distanceMapFilter->SetInput(image);
  distanceMapFilter->SetInsideIsPositive(true);
  distanceMapFilter->GetOutput()->SetRequestedRegion(boundingBox);
  distanceMapFilter->Update();

  // We want to process pixels with the smallest distance first (closer to the boundary)
  std::priority_queue<WeightedPixel> queue;
  std::vector<float> weights;
  
  itk::ImageRegionConstIteratorWithIndex<Mask> maskIterator(mask, mask->GetLargestPossibleRegion());

  while(!maskIterator.IsAtEnd())
    {
    if(mask->IsHole(maskIterator.GetIndex()))
      {
      WeightedPixel weightedPixel(maskIterator.GetIndex(),
                                  distanceMapFilter->GetOutput()->GetPixel(maskIterator.GetIndex()));
      queue.push(weightedPixel);
      weights.push_back(weightedPixel.Weight);
      }

    ++maskIterator;
    }

  // Here we need to find the distance from each hole pixel to *every* boundary pixel. Not sure that we need the
  // distance map computed above, as that is just the distance to the *closest* boundary pixel.

//  // Find the sum of the weights
//  float weightSum = Helpers::Sum(weights);

//  // "Inverse Distance Weighting Interpolation"
  
//  while (!queue.empty())
//  {
//    WeightedPixel p = queue.top();   //print out the highest priority element

//    image->SetPixel(p.Pixel, value);

//    queue.pop(); // remove the highest priority element
//  }
}

template<typename TImage>
void InterpolateThroughHole(TImage* const image, Mask* const mask, const itk::Index<2>& p0,
                            const itk::Index<2>& p1, const unsigned int lineThickness)
{
  // This function sets the pixels on the line in the mask to valid and sets the corresponding pixels
  // in the image to the interpolated values.
  if(mask->IsHole(p0) || mask->IsHole(p1))
  {
    throw std::runtime_error("Both p0 and p1 must be valid (not holes)!");
  }
  
  itk::BresenhamLine<2> line;

  std::vector< itk::Index<2> > pixels = line.BuildLine(p0, p1);
  std::cout << "Line contains " << pixels.size() << " pixels." << std::endl;
  
  std::vector< itk::Index<2> > holePixels;

  // Find the start and end of the line
  typename TImage::PixelType value0;
  unsigned int firstHolePixelIndex = 0;

  // Look for the first hole pixel, and set value0 to the one before it (the pixel on the valid side of the hole boundary)
  unsigned int pixelId = 0;
  for(; pixelId < pixels.size(); ++pixelId)
    {
    // std::cout << "pixel " << i << " " << pixels[i] << std::endl;
    if(mask->IsHole(pixels[pixelId]))
      {
      firstHolePixelIndex = pixelId;
      value0 = image->GetPixel(pixels[pixelId-1]);
      break;
      }
    else
      {
      //std::cout << "Skipping pixel " << pixels[pixelId] << " while looking for beginning of hole." << std::endl;
      }
    }

  std::cout << "First pixel in hole ID " << firstHolePixelIndex << std::endl;
    
  typename TImage::PixelType value1;
  // Look for the last hole pixel, and set value1 to the one after it
  // (the pixel on the valid side of the hole boundary)
  unsigned int lastHolePixelIndex = 0;
  for(; pixelId < pixels.size(); ++pixelId)
    {
    // std::cout << "pixel " << i << " " << pixels[i] << std::endl;
    if(mask->IsValid(pixels[pixelId])) // We found a pixel on the other side of the hole.
      {
      value1 = image->GetPixel(pixels[pixelId]);
      lastHolePixelIndex = pixelId - 1; // The previous pixel is the last one in the hole.
      break;
      }
    else
      {
      //std::cout << "Skipping pixel " << pixels[pixelId] << " while looking for end of hole." << std::endl;
      }
    }

  float difference = value1 - value0;

  std::cout << "Last pixel in hole ID " << lastHolePixelIndex << std::endl;

  unsigned int numberOfPixelsInHole = lastHolePixelIndex - firstHolePixelIndex;
  float step = difference / static_cast<float>(numberOfPixelsInHole);
  std::cout << "There are " << numberOfPixelsInHole << " pixels in the hole." << std::endl;
  
  if(lastHolePixelIndex - firstHolePixelIndex == 0)
  {
    throw std::runtime_error("Something is wrong, there are zero pixels to change!");
  }
  
  for(unsigned int holePixelId = firstHolePixelIndex; holePixelId <= lastHolePixelIndex; ++holePixelId)
    {
    std::cout << "Changing pixel " << holePixelId << " " << pixels[holePixelId] << std::endl;
//     if(!mask->IsHole(pixels[holePixelId]))
//       {
//       throw std::runtime_error("Something went wrong, we should only have hole pixels!");
//       }
    
  // For a line thickness of 0
    image->SetPixel(pixels[holePixelId], value0 + holePixelId * step);
    mask->SetPixel(pixels[holePixelId], mask->GetValidValue());

    // For a hacky thicker line
    std::vector<itk::Index<2> > indices = ITKHelpers::GetIndicesInRegion(ITKHelpers::GetRegionInRadiusAroundPixel(pixels[holePixelId], lineThickness));
    for(unsigned int neighborId = 0; neighborId < indices.size(); ++neighborId)
      {
      if(mask->IsHole(indices[neighborId]))
        {
        image->SetPixel(indices[neighborId], value0 + holePixelId * step);
        mask->SetPixel(indices[neighborId], mask->GetValidValue());
        }
      }
    }
}

template<typename TImage>
void InteroplateLineBetweenPointsWithFilling(TImage* const image, Mask* const mask,
                                             const itk::Index<2>& p0, const itk::Index<2>& p1)
{
  itk::BresenhamLine<2> line;

  std::vector< itk::Index<2> > pixels = line.BuildLine(p0, p1);

  typename TImage::PixelType value0 = image->GetPixel(p0);
  typename TImage::PixelType value1 = image->GetPixel(p1);

  float difference = value1 - value0;
  float step = difference / static_cast<float>(pixels.size());

  for(unsigned int i = 0; i < pixels.size(); i++)
    {
    //std::cout << "pixel " << i << " " << pixels[i] << std::endl;
    image->SetPixel(pixels[i], value0 + i * step);
    mask->SetPixel(pixels[i], mask->GetValidValue());
    }
}

template<typename TImage>
void BlurInHole(TImage* const image, const Mask* const mask, const float kernelVariance)
{
  typedef itk::DiscreteGaussianImageFilter<TImage, TImage> DiscreteGaussianImageFilterType;

  // Create and setup a Gaussian filter
  typename DiscreteGaussianImageFilterType::Pointer gaussianFilter = DiscreteGaussianImageFilterType::New();
  gaussianFilter->SetInput(image);

  typename DiscreteGaussianImageFilterType::ArrayType varianceArray;
  varianceArray.Fill(kernelVariance);
  gaussianFilter->SetVariance(varianceArray);
  gaussianFilter->Update();

  itk::ImageRegionIteratorWithIndex<TImage> imageIterator(image, image->GetLargestPossibleRegion());

  while(!imageIterator.IsAtEnd())
    {
    if(mask->IsHole(imageIterator.GetIndex()))
      {
      imageIterator.Set(gaussianFilter->GetOutput()->GetPixel(imageIterator.GetIndex()));
      }

    ++imageIterator;
    }
}

template<typename TImage>
void MedianFilterInHole(TImage* const image, const Mask* const mask, const unsigned int kernelRadius)
{
  std::cout << "Median filtering with radius " << kernelRadius << std::endl;
  typedef itk::MedianImageFilter<TImage, TImage> MedianFilterType;
  typename MedianFilterType::Pointer medianFilter = MedianFilterType::New();
  typename MedianFilterType::InputSizeType radius;
  radius.Fill(kernelRadius);
  medianFilter->SetRadius(radius);
  medianFilter->SetInput(image);
  medianFilter->Update();

  itk::ImageRegionIteratorWithIndex<TImage> imageIterator(image, image->GetLargestPossibleRegion());

  while(!imageIterator.IsAtEnd())
    {
    if(mask->IsHole(imageIterator.GetIndex()))
      {
//       std::cout << "Changing " << image->GetPixel(imageIterator.GetIndex()) << " to "
//                 << medianFilter->GetOutput()->GetPixel(imageIterator.GetIndex()) << std::endl;
      imageIterator.Set(medianFilter->GetOutput()->GetPixel(imageIterator.GetIndex()));
      }

    ++imageIterator;
    }
}

/** Clip the values in the image inside the hole. */
template<typename TImage>
void ClipInHole(TImage* const image, const Mask* const mask, const float min, const float max)
{
  itk::ImageRegionIteratorWithIndex<TImage> imageIterator(image, image->GetLargestPossibleRegion());

  while(!imageIterator.IsAtEnd())
    {
    if(mask->IsHole(imageIterator.GetIndex()))
      {
      //std::cout << "Changing " << image->GetPixel(imageIterator.GetIndex()) << " to "
      //          << medianFilter->GetOutput()->GetPixel(imageIterator.GetIndex()) << std::endl;
      if(imageIterator.Get() < min)
      {
        imageIterator.Set(min);
      }
      if(imageIterator.Get() > max)
      {
        imageIterator.Set(max);
      }
      }

    ++imageIterator;
    }
}

template<typename TImage>
void AddConstantInHole(TImage* const image, const float value, const Mask* const mask)
{
  itk::ImageRegionIteratorWithIndex<TImage> imageIterator(image, image->GetLargestPossibleRegion());

  while(!imageIterator.IsAtEnd())
    {
    if(mask->IsHole(imageIterator.GetIndex()))
      {
      imageIterator.Set(value);
      }

    ++imageIterator;
    }
}

template<typename TImage>
typename TImage::PixelType AverageHoleValue(const TImage* const image, const Mask* const mask)
{
  itk::ImageRegionConstIteratorWithIndex<TImage> imageIterator(image, image->GetLargestPossibleRegion());

  typename TImage::PixelType sum(image->GetNumberOfComponentsPerPixel());
  sum.Fill(0);
  unsigned int numberOfHolePixels = 0;
  while(!imageIterator.IsAtEnd())
    {
    if(mask->IsHole(imageIterator.GetIndex()))
      {
      sum += imageIterator.Get();
      numberOfHolePixels++;
      }

    ++imageIterator;
    }
  std::cout << "Sum: " << sum << std::endl;
  std::cout << "numberOfHolePixels: " << numberOfHolePixels << std::endl;
  
  return sum/static_cast<float>(numberOfHolePixels);
}


template <class TImage>
void CopyPatchIntoImage(const TImage* const patch, TImage* const image, const Mask* const mask,
                        const itk::Index<2>& position)
{
  // This function copies 'patch' into 'image' centered at 'position' only where the 'mask' is non-zero

  // 'Mask' must be the same size as 'image'
  if(mask->GetLargestPossibleRegion().GetSize() != image->GetLargestPossibleRegion().GetSize())
    {
    throw std::runtime_error("mask and image must be the same size!");
    }

  // The PasteFilter expects the lower left corner of the destination position, but we have passed the center pixel.
  position[0] -= patch->GetLargestPossibleRegion().GetSize()[0]/2;
  position[1] -= patch->GetLargestPossibleRegion().GetSize()[1]/2;

  itk::ImageRegion<2> region = GetRegionInRadiusAroundPixel(position,
                                                            patch->GetLargestPossibleRegion().GetSize()[0]/2);

  itk::ImageRegionConstIterator<TImage> patchIterator(patch,patch->GetLargestPossibleRegion());
  itk::ImageRegionConstIterator<Mask> maskIterator(mask,region);
  itk::ImageRegionIterator<TImage> imageIterator(image, region);

  while(!patchIterator.IsAtEnd())
    {
    if(mask->IsHole(maskIterator.GetIndex())) // we are in the target region
      {
      imageIterator.Set(patchIterator.Get());
      }
    ++imageIterator;
    ++maskIterator;
    ++patchIterator;
    }
}


/** Compute the average of all unmasked pixels in a region.*/
template<typename TImage>
typename TypeTraits<typename TImage::PixelType>::LargerType AverageInRegionMasked(const TImage* const image,
                                                                                  const Mask* const mask,
                                                                                  const itk::ImageRegion<2>& region)
{
  typename itk::ImageRegionConstIteratorWithIndex<Mask> maskIterator(mask, region);

  std::vector<typename TImage::PixelType> pixels;

  while(!maskIterator.IsAtEnd())
    {
    if(mask->IsValid(maskIterator.GetIndex()))
      {
      pixels.push_back(image->GetPixel(maskIterator.GetIndex()));
      }
    ++maskIterator;
    }

  using Statistics::Average;
  //using ITKHelpers::Average;
  return Average(pixels);
}

template<typename TImage>
typename TImage::PixelType AverageValidNeighborValue(const TImage* const image, const Mask* const mask,
                                                     const itk::Index<2>& pixel)
{
  std::vector<itk::Index<2> > neighborIndices = ITKHelpers::Get8NeighborsInRegion(image->GetLargestPossibleRegion(),
                                                                                  pixel);

  std::vector<typename TImage::PixelType> pixels;

  for(unsigned int i = 0; i < neighborIndices.size(); ++i)
    {
    if(mask->IsValid(neighborIndices[i]))
      {
      pixels.push_back(image->GetPixel(neighborIndices[i]));
      }
    }

  return Statistics::Average(pixels);
}

template<typename TImage>
typename TImage::PixelType AverageHoleNeighborValue(const TImage* const image, const Mask* const mask,
                                                      const itk::Index<2>& pixel)
{
  std::vector<itk::Index<2> > neighborIndices = ITKHelpers::Get8NeighborsInRegion(image->GetLargestPossibleRegion(),
                                                                                  pixel);

  std::vector<typename TImage::PixelType> pixels;

  for(unsigned int i = 0; i < neighborIndices.size(); ++i)
    {
    if(mask->IsHole(neighborIndices[i]))
      {
      pixels.push_back(image->GetPixel(neighborIndices[i]));
      }
    }

  return Statistics::Average(pixels);
}

/** Compute the average of all unmasked pixels in a region.*/
template<typename TImage>
typename TypeTraits<typename TImage::PixelType>::LargerType VarianceInRegionMasked(const TImage* const image,
                                                                                   const Mask* const mask,
                                                                             const itk::ImageRegion<2>& region)
{
  typename itk::ImageRegionConstIterator<Mask> maskIterator(mask, region);
  std::vector<typename TImage::PixelType> pixels;
  while(!maskIterator.IsAtEnd())
    {
    if(mask->IsValid(maskIterator.GetIndex()))
      {
      pixels.push_back(image->GetPixel(maskIterator.GetIndex()));
      }
    ++maskIterator;
    }

  return Statistics::Variance(pixels);
}

template<typename TImage>
void WriteMaskedRegion(const TImage* const image, const Mask* mask, const itk::ImageRegion<2>& region,
                       const std::string& filename, const typename TImage::PixelType& holeColor)
{
  typedef itk::RegionOfInterestImageFilter<TImage, TImage> RegionOfInterestImageFilterType;
  typename RegionOfInterestImageFilterType::Pointer regionOfInterestImageFilter =
            RegionOfInterestImageFilterType::New();
  regionOfInterestImageFilter->SetRegionOfInterest(region);
  regionOfInterestImageFilter->SetInput(image);
  regionOfInterestImageFilter->Update();

  typedef itk::RegionOfInterestImageFilter<Mask, Mask> RegionOfInterestMaskFilterType;
  typename RegionOfInterestMaskFilterType::Pointer regionOfInterestMaskFilter =
            RegionOfInterestMaskFilterType::New();
  regionOfInterestMaskFilter->SetRegionOfInterest(region);
  regionOfInterestMaskFilter->SetInput(mask);
  regionOfInterestMaskFilter->Update();

  itk::ImageRegionIterator<TImage> imageIterator(regionOfInterestImageFilter->GetOutput(),
                                                 regionOfInterestImageFilter->GetOutput()->
                                                 GetLargestPossibleRegion());

  while(!imageIterator.IsAtEnd())
    {
    typename TImage::PixelType pixel = imageIterator.Get();

    itk::Index<2> index = imageIterator.GetIndex();

    if(regionOfInterestMaskFilter->GetOutput()->IsHole(imageIterator.GetIndex()))
      {
      regionOfInterestImageFilter->GetOutput()->SetPixel(index, holeColor);
      }

    ++imageIterator;
    }

  typename itk::ImageFileWriter<TImage>::Pointer writer = itk::ImageFileWriter<TImage>::New();
  writer->SetFileName(filename);
  writer->SetInput(regionOfInterestImageFilter->GetOutput());
  writer->Update();
}

template<typename TImage>
void WriteMaskedRegionPNG(const TImage* const image, const Mask* mask,
                          itk::ImageRegion<2> region, const std::string& filename,
                          const typename TImage::PixelType& holeColor)
{
  region.Crop(image->GetLargestPossibleRegion());

  if(region.GetSize()[0] == 0 || region.GetSize()[1] == 0 )
  {
    throw std::runtime_error("Cropped region is 0 in at least one dimension!");
  }

  typedef itk::RegionOfInterestImageFilter<TImage, TImage> RegionOfInterestImageFilterType;
  typename RegionOfInterestImageFilterType::Pointer regionOfInterestImageFilter =
            RegionOfInterestImageFilterType::New();
  regionOfInterestImageFilter->SetRegionOfInterest(region);
  regionOfInterestImageFilter->SetInput(image);
  regionOfInterestImageFilter->Update();

  typedef itk::RegionOfInterestImageFilter<Mask, Mask> RegionOfInterestMaskFilterType;
  typename RegionOfInterestMaskFilterType::Pointer regionOfInterestMaskFilter =
            RegionOfInterestMaskFilterType::New();
  regionOfInterestMaskFilter->SetRegionOfInterest(region);
  regionOfInterestMaskFilter->SetInput(mask);
  regionOfInterestMaskFilter->Update();

  itk::ImageRegionIterator<TImage> imageIterator(regionOfInterestImageFilter->GetOutput(),
                                                 regionOfInterestImageFilter->GetOutput()->
                                                 GetLargestPossibleRegion());

  while(!imageIterator.IsAtEnd())
  {
//    typename TImage::PixelType pixel = imageIterator.Get();

    itk::Index<2> index = imageIterator.GetIndex();

    if(regionOfInterestMaskFilter->GetOutput()->IsHole(imageIterator.GetIndex()))
    {
      regionOfInterestImageFilter->GetOutput()->SetPixel(index, holeColor);
    }

    ++imageIterator;
  }

  ITKHelpers::RGBImageType::Pointer rgbImage = ITKHelpers::RGBImageType::New();
  ITKHelpers::VectorImageToRGBImage(regionOfInterestImageFilter->GetOutput(), rgbImage.GetPointer());

  ITKHelpers::WriteImage(rgbImage.GetPointer(), filename);
}


// This struct is used inside MaskedBlur()
template <typename TPixel>
struct Contribution
{
  Contribution() : Weight(0.0f), Value(itk::NumericTraits<TPixel>::Zero){}
  float Weight;
  TPixel Value;
  itk::Offset<2> Offset;
};

template <typename TImage>
void MaskedBlur(const TImage* const inputImage, const Mask* const mask, const float blurVariance,
                TImage* const output)
{
//  itk::SimpleFastMutexLock mutex;

//  mutex.Lock();
  MaskedBlurInRegion(inputImage, mask, inputImage->GetLargestPossibleRegion(), blurVariance, output);
//  mutex.Unlock();
}

/** Blur the 'image' only where 'mask' is valid, and only using pixels where 'mask' is valid. */
template <typename TImage>
void MaskedBlurInRegion(const TImage* const inputImage, const Mask* const mask,
                        const itk::ImageRegion<2>& region,
                        const float blurVariance, TImage* const output)

{
//  std::cout << "MaskedBlurInRegion()" << region << std::endl;
//  std::cout << "image region: " << inputImage->GetLargestPossibleRegion()
//            << " mask region: " << mask->GetLargestPossibleRegion() << std::endl;
  assert(inputImage->GetLargestPossibleRegion() == mask->GetLargestPossibleRegion());
  assert(mask->GetLargestPossibleRegion().IsInside(region));

  if(!mask->GetLargestPossibleRegion().IsInside(region))
  {
    std::cout << "MaskedBlurInRegion: region is not inside image!" << region << std::endl;
  }

  // Create a Gaussian kernel
  typedef itk::GaussianOperator<float, 1> GaussianOperatorType;

  // Make a (2*kernelRadius+1)x1 kernel
  itk::Size<1> radius;
  radius.Fill(20); // Make a length 41 kernel

  GaussianOperatorType gaussianOperator;
  // It doesn't matter which direction we set - we will be interpreting the kernel as 1D (no direction)
  gaussianOperator.SetDirection(0); 
  gaussianOperator.SetVariance(blurVariance);
  gaussianOperator.CreateToRadius(radius);

//   {
//   // Debugging only
//   std::cout << "gaussianOperator: " << gaussianOperator << std::endl;
//   for(unsigned int i = 0; i < gaussianOperator.Size(); i++)
//     {
//     //std::cout << i << " : " << gaussianOperator.GetOffset(i) << std::endl;
//     std::cout << i << " : " << gaussianOperator.GetElement(i) << std::endl;
//     }
//   }

  // Create the output image - data will be deep copied into it
//  std::cout << "MaskedBlurInRegion: creating blurredImage" << std::endl;
//  std::cout << "MaskedBlurInRegion: inputImage region " << inputImage->GetLargestPossibleRegion() << std::endl;
  typename TImage::Pointer blurredImage = TImage::New();
  ITKHelpers::InitializeImage(blurredImage.GetPointer(), inputImage->GetLargestPossibleRegion());
//  std::cout << "MaskedBlurInRegion: after init" << std::endl;
  // We apply the filter to the same image, one dimension at a time. To do this,
  // instead of applying it to the 'input' the first pass, then needing to copy
  // the output to the input before running the next pass, we instead create
  // an intermediate image that is used each pass.
  typename TImage::Pointer operatingImage = TImage::New();
//  operatingImage->SetRegions(region);
  operatingImage->SetRegions(inputImage->GetLargestPossibleRegion());
  operatingImage->Allocate();
//  std::cout << "MaskedBlurInRegion: deep copy" << std::endl;
  ITKHelpers::DeepCopyInRegion(inputImage, region, operatingImage.GetPointer());
//  std::cout << "MaskedBlurInRegion: after deep copy" << std::endl;
  for(unsigned int dimensionPass = 0; dimensionPass < 2; dimensionPass++) // The image is 2D
  {
//    std::cout << "MaskedBlurInRegion: inner iter" << std::endl;
    itk::ImageRegionIterator<TImage> imageIterator(operatingImage, region);
//    std::cout << "MaskedBlurInRegion: after inner iter" << std::endl;
    while(!imageIterator.IsAtEnd())
    {
      itk::Index<2> centerPixel = imageIterator.GetIndex();

      // We should not compute derivatives for pixels in the hole.
//      std::cout << "MaskedBlurInRegion: checking IsHole()" << std::endl;
      if(mask->IsHole(centerPixel))
      {
        ++imageIterator;
        continue;
      }

      // Loop over all of the pixels in the kernel and use the ones that fit a criteria
      typedef Contribution<typename TImage::PixelType> ContributionType;
      std::vector<ContributionType> contributions;
      for(unsigned int i = 0; i < gaussianOperator.Size(); i++)
      {
        // Since we use 1D kernels, we must manually construct a 2D offset with 0 in all
        // dimensions except the dimension of the current pass
        itk::Offset<2> offset = ITKHelpers::OffsetFrom1DOffset(gaussianOperator.GetOffset(i), dimensionPass);

        itk::Index<2> index = centerPixel + offset;
        if(blurredImage->GetLargestPossibleRegion().IsInside(index) && mask->IsValid(index))
        {
          ContributionType contribution;
          contribution.Weight = gaussianOperator.GetElement(i);
          contribution.Value = operatingImage->GetPixel(index);
          contribution.Offset = offset;
          contributions.push_back(contribution);
        }
      }

      if(contributions.size() == 0)
      {
        std::stringstream ss;
        ss << "Pixel " << centerPixel << " does not have any valid neighbors!";
        throw std::runtime_error(ss.str());
      }

      float totalWeight = 0.0f;
      for(unsigned int i = 0; i < contributions.size(); i++)
      {
        totalWeight += contributions[i].Weight;
      }

      // Determine the new pixel value
      typedef typename TypeTraits<typename TImage::PixelType>::LargerType LargerType;
      LargerType newPixelValue = itk::NumericTraits<LargerType>::Zero;

      for(unsigned int i = 0; i < contributions.size(); i++)
      {
        itk::Index<2> contributionIndex = centerPixel + contributions[i].Offset;
        float contributionWeight = contributions[i].Weight/totalWeight;
        newPixelValue += static_cast<LargerType>(operatingImage->GetPixel(contributionIndex))
                         * contributionWeight;
      }

      blurredImage->SetPixel(centerPixel, newPixelValue);
      ++imageIterator;
    }

    // For the separable Gaussian filtering concept to work,
    // the next pass must operate on the output of the current pass.
//    std::cout << "MaskedBlurInRegion: copying into operating image" << std::endl;
    ITKHelpers::DeepCopy(blurredImage.GetPointer(), operatingImage.GetPointer());
//    std::cout << "MaskedBlurInRegion: after copying into operating image" << std::endl;
  }

  // Copy the final image to the output.
//  std::cout << "MaskedBlurInRegion: copying into output" << std::endl;
  ITKHelpers::DeepCopy(blurredImage.GetPointer(), output);
//  std::cout << "MaskedBlurInRegion: after copying into output" << std::endl;
}


template<typename TImage>
void CopyAtValues(const TImage* const input, const Mask::PixelType& value,
                  const Mask* const mask, TImage* const output)
{
  // It is sometimes desired to copy the hole pixels from an image into the corresponding pixels
  // in a larger image, so we do not do the following.
//   if(input->GetLargestPossibleRegion().GetSize() != output->GetLargestPossibleRegion().GetSize())
//   {
//     std::stringstream ss;
//     ss << "Input size (" << input->GetLargestPossibleRegion().GetSize() << ") must match output size ("
//        << output->GetLargestPossibleRegion().GetSize() << ")";
//     throw std::runtime_error(ss.str());
//   }

  if(!output->GetLargestPossibleRegion().IsInside(input->GetLargestPossibleRegion()))
  {
    std::stringstream ss;
    ss << "Input is not smaller than output! Input: " << input->GetLargestPossibleRegion().GetSize()
       << " output: " << output->GetLargestPossibleRegion().GetSize();
    throw std::runtime_error(ss.str());
  }

  if(input->GetLargestPossibleRegion().GetSize() != output->GetLargestPossibleRegion().GetSize())
  {
    std::stringstream ss;
    ss << "Input size (" << input->GetLargestPossibleRegion().GetSize() << ") must match output size ("
       << output->GetLargestPossibleRegion().GetSize() << ")";
    throw std::runtime_error(ss.str());
  }

  itk::ImageRegionConstIterator<Mask> maskIterator(mask, mask->GetLargestPossibleRegion());

  while(!maskIterator.IsAtEnd())
  {
    if(maskIterator.Get() == value)
    {
      output->SetPixel(maskIterator.GetIndex(), input->GetPixel(maskIterator.GetIndex()));
    }
    ++maskIterator;
  }
}

template <class TImage>
void CopyInHoleRegion(const TImage* const input, TImage* const output, const Mask* const mask)
{
  CopyAtValues(input, mask->GetHoleValue(), mask, output);
}

template <class TImage>
void CopyInValidRegion(const TImage* const input, TImage* const output, const Mask* const mask)
{
  CopyAtValues(input, mask->GetValidValue(), mask, output);
}

template<typename TImage>
void SetHolePixelsToConstant(TImage* const image, const typename TImage::PixelType& value,
                             const Mask* const mask)
{
  itk::ImageRegionIterator<TImage> imageIterator(image, image->GetLargestPossibleRegion());
  
  while(!imageIterator.IsAtEnd())
  {
    if(mask->IsHole(imageIterator.GetIndex()))
    {
      imageIterator.Set(value);
    }
    ++imageIterator;
  }
}

#ifdef MaskUseVTK
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
#endif // #if MaskUseVTK

} // end namespace
