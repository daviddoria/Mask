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
#include "ITKHelpers/ITKHelpers.h"

// ITK
#include "itkBresenhamLine.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkGaussianOperator.h"
#include "itkImageRegionIterator.h"
#include "itkLaplacianOperator.h"
#include "itkMedianImageFilter.h"
#include "itkSignedDanielssonDistanceMapImageFilter.h"

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
void CopySourcePatchIntoHoleOfTargetRegion(const TImage* const sourceImage, TImage* const targetImage,
                                           const Mask* const mask,
                                           const itk::ImageRegion<2>& sourceRegionInput,
                                           const itk::ImageRegion<2>& destinationRegionInput)
{
  itk::ImageRegion<2> fullImageRegion = sourceImage->GetLargestPossibleRegion();

  // We pass the regions by const reference, so copy them here before they are mutated
  itk::ImageRegion<2> sourceRegion = sourceRegionInput;
  itk::ImageRegion<2> destinationRegion = destinationRegionInput;

  // Move the source region to the desintation region
  itk::Offset<2> offset = destinationRegion.GetIndex() - sourceRegion.GetIndex();
  sourceRegion.SetIndex(sourceRegion.GetIndex() + offset);

  // Make the destination be entirely inside the image
  destinationRegion.Crop(fullImageRegion);
  sourceRegion.Crop(fullImageRegion);

  // Move the source region back
  sourceRegion.SetIndex(sourceRegion.GetIndex() - offset);

  itk::ImageRegionConstIterator<TImage> sourceIterator(sourceImage, sourceRegion);
  itk::ImageRegionIterator<TImage> destinationIterator(targetImage, destinationRegion);
  itk::ImageRegionConstIterator<Mask> maskIterator(mask, destinationRegion);

  while(!sourceIterator.IsAtEnd())
    {
    if(mask->IsHole(maskIterator.GetIndex())) // we are in the target region
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
itk::Index<2> FindHighestValueInMaskedRegion(const TImage* const image, float& maxValue, const Mask* const maskImage)
{
  //EnterFunction("FindHighestValueOnBoundary()");
  // Return the location of the highest pixel in 'image' out of the non-zero pixels
  // in 'boundaryImage'. Return the value of that pixel by reference.

  // Explicity find the maximum on the boundary
  maxValue = 0.0f; // priorities are non-negative, so anything better than 0 will win

  std::vector<itk::Index<2> > boundaryPixels = ITKHelpers::GetNonZeroPixels(maskImage);

  if(boundaryPixels.size() <= 0)
    {
    throw std::runtime_error("FindHighestValueOnBoundary(): No boundary pixels!");
    }

  itk::Index<2> locationOfMaxValue = boundaryPixels[0];

  for(unsigned int i = 0; i < boundaryPixels.size(); ++i)
    {
    if(image->GetPixel(boundaryPixels[i]) > maxValue)
      {
      maxValue = image->GetPixel(boundaryPixels[i]);
      locationOfMaxValue = boundaryPixels[i];
      }
    }
  //DebugMessage<float>("Highest value: ", maxValue);
  //LeaveFunction("FindHighestValueOnBoundary()");
  return locationOfMaxValue;
}

template<typename TImage, typename TRegionIndicatorImage>
itk::Index<2> FindHighestValueInNonZeroRegion(const TImage* const image, float& maxValue,
                                              const TRegionIndicatorImage* const indicatorImage)
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

  itk::ImageRegion<2> boundingBox = MaskOperations::ComputeHoleBoundingBox(mask);

  typedef itk::SignedDanielssonDistanceMapImageFilter<TImage, ITKHelpersTypes::FloatScalarImageType>
          SignedDanielssonDistanceMapImageFilterType;
  typename SignedDanielssonDistanceMapImageFilterType::Pointer distanceMapFilter =
           SignedDanielssonDistanceMapImageFilterType::New();
  distanceMapFilter->SetInput(image);
  distanceMapFilter->SetInsideIsPositive(true);
  distanceMapFilter->GetOutput()->SetRequestedRegion(boundingBox);
  distanceMapFilter->Update();

  // the first element will be the one with the smallest value
  std::priority_queue <WeightedPixel, std::vector<WeightedPixel> > queue;
  
  itk::ImageRegionConstIteratorWithIndex<Mask> maskIterator(mask, mask->GetLargestPossibleRegion());

  while(!maskIterator.IsAtEnd())
    {
    if(mask->IsHole(maskIterator.GetIndex()))
      {
      WeightedPixel weightedPixel(maskIterator.GetIndex(),
                                  distanceMapFilter->GetOutput()->GetPixel(maskIterator.GetIndex()));
      queue.push(weightedPixel);
      }

    ++maskIterator;
    }

  // TODO: Do "Inverse Distance Weighting Interpolation" here
  
  while (!queue.empty())
  {
    WeightedPixel p = queue.top();   //print out the highest priority element
    queue.pop();                   //remove the highest priority element
  }
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
      std::cout << "Changing " << image->GetPixel(imageIterator.GetIndex()) << " to "
                << medianFilter->GetOutput()->GetPixel(imageIterator.GetIndex()) << std::endl;
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
                          const itk::ImageRegion<2>& region, const std::string& filename,
                          const typename TImage::PixelType& holeColor)
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

  ITKHelpers::RGBImageType::Pointer rgbImage = ITKHelpers::RGBImageType::New();
  //ITKHelpers::VectorImageToRGBImage(regionOfInterestImageFilter->GetOutput(), rgbImage.GetPointer());

  typename itk::ImageFileWriter<ITKHelpers::RGBImageType>::Pointer writer =
             itk::ImageFileWriter<ITKHelpers::RGBImageType>::New();
  writer->SetFileName(filename);
  writer->SetInput(rgbImage.GetPointer());
  writer->Update();
}


// This struct is used inside MaskedBlur()
struct Contribution
{
  float weight;
  unsigned char value;
  itk::Offset<2> offset;
};

template <typename TImage>
void MaskedBlur(const TImage* const inputImage, const Mask* const mask, const float blurVariance,
                TImage* const output)
{
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
  typename TImage::Pointer blurredImage = TImage::New();
  ITKHelpers::InitializeImage<TImage>(blurredImage, inputImage->GetLargestPossibleRegion());

  // Initialize
  typename TImage::Pointer operatingImage = TImage::New();
  ITKHelpers::DeepCopy(inputImage, operatingImage.GetPointer());

  for(unsigned int dimensionPass = 0; dimensionPass < 2; dimensionPass++) // The image is 2D
    {
    itk::ImageRegionIterator<TImage> imageIterator(operatingImage, operatingImage->GetLargestPossibleRegion());

    while(!imageIterator.IsAtEnd())
      {
      itk::Index<2> centerPixel = imageIterator.GetIndex();

      // We should not compute derivatives for pixels in the hole.
      if(mask->IsHole(centerPixel))
        {
        ++imageIterator;
        continue;
        }

      // Loop over all of the pixels in the kernel and use the ones that fit a criteria
      std::vector<Contribution> contributions;
      for(unsigned int i = 0; i < gaussianOperator.Size(); i++)
        {
        // Since we use 1D kernels, we must manually construct a 2D offset with 0 in all
        // dimensions except the dimension of the current pass
        itk::Offset<2> offset = ITKHelpers::OffsetFrom1DOffset(gaussianOperator.GetOffset(i), dimensionPass);

        itk::Index<2> pixel = centerPixel + offset;
        if(blurredImage->GetLargestPossibleRegion().IsInside(pixel) && mask->IsValid(pixel))
          {
          Contribution contribution;
          contribution.weight = gaussianOperator.GetElement(i);
          contribution.value = operatingImage->GetPixel(pixel);
          contribution.offset = ITKHelpers::OffsetFrom1DOffset(gaussianOperator.GetOffset(i), dimensionPass);
          contributions.push_back(contribution);
          }
        }

      float total = 0.0f;
      for(unsigned int i = 0; i < contributions.size(); i++)
        {
        total += contributions[i].weight;
        }

      // Determine the new pixel value
      float newPixelValue = 0.0f;
      for(unsigned int i = 0; i < contributions.size(); i++)
        {
        itk::Index<2> pixel = centerPixel + contributions[i].offset;
        newPixelValue += contributions[i].weight/total * operatingImage->GetPixel(pixel);
        }

      blurredImage->SetPixel(centerPixel, newPixelValue);
      ++imageIterator;
      }

    // For the separable Gaussian filtering concept to work,
    // the next pass must operate on the output of the current pass.
    ITKHelpers::DeepCopy(blurredImage.GetPointer(), operatingImage.GetPointer());
    }

  // Copy the final image to the output.
  ITKHelpers::DeepCopy(blurredImage.GetPointer(), output);
}

template <class TImage>
void CopyInHoleRegion(const TImage* const input, TImage* const output, const Mask* const mask)
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

  if(input->GetLargestPossibleRegion().GetSize() != mask->GetLargestPossibleRegion().GetSize())
  {
    std::stringstream ss;
    ss << "Input size (" << input->GetLargestPossibleRegion().GetSize() << ") must match mask size ("
       << mask->GetLargestPossibleRegion().GetSize() << ")";
    throw std::runtime_error(ss.str());
  }
  
  itk::ImageRegionConstIterator<TImage> imageIterator(input, input->GetLargestPossibleRegion());
  //std::cout << "Hole value is " << static_cast<int>(mask->GetHoleValue()) << std::endl;
  while(!imageIterator.IsAtEnd())
  {
    if(mask->IsHole(imageIterator.GetIndex()))
    {
      output->SetPixel(imageIterator.GetIndex(), imageIterator.Get());
    }
    ++imageIterator;
  }
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

} // end namespace
