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

#ifndef MaskOperations_H
#define MaskOperations_H

// Custom
class Mask;
#include "ITKHelpers/ITKHelpers.h"

// ITK
#include "itkIndex.h"
#include "itkImageRegion.h"

// VTK
class vtkImageData;

namespace MaskOperations
{

/** Write a 'region' of an 'image' to 'filename', coloring any invalid pixels in 'mask' the color 'holeColor'. */
template<typename TImage>
void WriteMaskedRegion(const TImage* const image, const Mask* mask, const itk::ImageRegion<2>& region,
                       const std::string& filename,
                       const typename TImage::PixelType& holeColor);

/** Write a 'region' of an 'image' to 'filename', coloring any invalid pixels in 'mask' the color 'holeColor'. */
template<typename TImage>
void WriteMaskedRegionPNG(const TImage* const image, const Mask* mask, const itk::ImageRegion<2>& region,
                          const std::string& filename, const typename TImage::PixelType& holeColor);

void ITKImageToVTKImageMasked(const ITKHelpers::FloatVectorImageType* const image, const Mask* const mask,
                              vtkImageData* const outputImage, const unsigned char maskColor[3]);

void SetMaskTransparency(const Mask* const input, vtkImageData* outputImage);

/** Return a random region that is entirely inside the hole. */
itk::ImageRegion<2> RandomRegionInsideHole(const Mask* const mask, const unsigned int halfWidth);

/** Return a random region that is entirely valid. */
itk::ImageRegion<2> RandomValidRegion(const Mask* const mask, const unsigned int halfWidth);

/** Compute the bounding box of the mask. */
itk::ImageRegion<2> ComputeHoleBoundingBox(const Mask* const mask);

/** Look from a pixel across the hole in a specified direction and return the
  * pixel that exists on the other side of the hole. */
itk::Index<2> FindPixelAcrossHole(const itk::Index<2>& queryPixel,
                                  const ITKHelpers::FloatVector2Type& direction, const Mask* const mask);


/** Get all regions of a particular size that contain only valid pixels.*/
std::vector<itk::ImageRegion<2> > GetAllFullyValidRegions(const Mask* const mask, const unsigned int patchRadius);

/** Get all regions of a particular size inside searchRegion that contain only valid pixels.*/
std::vector<itk::ImageRegion<2> > GetAllFullyValidRegions(const Mask* const mask,
                                                          const itk::ImageRegion<2>& searchRegion,
                                                          const unsigned int patchRadius);

/** Get a random fully valid patch in the specified region, if it exists.*/
itk::ImageRegion<2> GetRandomValidPatchInRegion(const Mask* const mask,
                                                const itk::ImageRegion<2>& searchRegion,
                                                const unsigned int patchRadius,
                                                const unsigned int maxNumberOfAttempts);

/** Get a random fully valid patch in the specified region. Try harder than the above function - takes much longer.*/
itk::ImageRegion<2> GetRandomValidPatchInRegion(const Mask* const mask,
                                                const itk::ImageRegion<2>& searchRegion,
                                                const unsigned int patchRadius);

////////////////// Templates ////////////////
template <typename TImage>
void MaskedBlur(const TImage* const inputImage, const Mask* const mask, const float blurVariance,
                TImage* const output);

template <class TImage>
void CopySelfPatchIntoHoleOfTargetRegion(TImage* const image, const Mask* const mask,
                                         const itk::ImageRegion<2>& sourceRegionInput,
                                         const itk::ImageRegion<2>& destinationRegionInput);

template <class TImage>
void CopyInHoleRegion(const TImage* const input, TImage* const output, const Mask* const mask);

template <class TImage>
void CopySourcePatchIntoHoleOfTargetRegion(const TImage* const sourceImage, TImage* const targetImage,
                                           const Mask* const mask,
                                           const itk::ImageRegion<2>& sourceRegionInput,
                                           const itk::ImageRegion<2>& destinationRegionInput);

template<typename TImage>
void CreatePatchImage(const TImage* const image, const itk::ImageRegion<2>& sourceRegion,
                      const itk::ImageRegion<2>& targetRegion, const Mask* const mask, TImage* const result);

template<typename TImage>
void AddConstantInHole(TImage* const image, const float value, const Mask* const maskImage);

template<typename TImage>
void SetHolePixelsToConstant(TImage* const image, const typename TImage::PixelType& value,
                             const Mask* const maskImage);

/** Return the highest value of the specified image out of the pixels under a specified BoundaryImage. */
template<typename TImage>
itk::Index<2> FindHighestValueInMaskedRegion(const TImage* const image, float& maxValue,
                                             const Mask* const maskImage);

template<typename TImage, typename TRegionIndicatorImage>
itk::Index<2> FindHighestValueInNonZero(const TImage* const image, float& maxValue,
                                        const TRegionIndicatorImage* const maskImage);

/** Get the average value of the masked pixels. */
template<typename TImage>
typename TImage::PixelType AverageHoleValue(const TImage* const image, const Mask* const mask);

/** Get the average value of the non-masked neighbors of a pixel. */
template<typename TImage>
typename TImage::PixelType AverageValidNeighborValue(const TImage* const image, const Mask* const mask,
                                                     const itk::Index<2>& pixel);

/** Get the average value of the masked neighbors of a pixel. */
template<typename TImage>
typename TImage::PixelType AverageHoleNeighborValue(const TImage* const image, const Mask* const mask,
                                                    const itk::Index<2>& pixel);

/** Get the average value of the masked neighbors of a pixel. */
template<typename TImage>
std::vector<typename TImage::PixelType> GetValidPixelsInRegion(const TImage* const image, const Mask* const mask,
                                                               const itk::ImageRegion<2>& region);

/** Get the average value of the masked neighbors of a pixel. */
template<typename TImage>
void AddNoiseInHole(TImage* const image, const Mask* const mask, const float noiseVariance);

/** Interpolate values through a hole, filling only pixels that are in the hole.
  * This function assumes one hole entry and one hole exit (i.e. the line between
  * p0 and p1 only intersects the hole twice). */
template<typename TImage>
void InterpolateThroughHole(TImage* const image, const Mask* const mask, const itk::Index<2>& p0,
                            const itk::Index<2>& p1, const unsigned int lineThickness = 0);

/** Interpolate values in a hole. */
template<typename TImage>
void InterpolateHole(TImage* const image, const Mask* const mask);

/** Blur an image using all of its values but only replaced the pixel values with
  * the blurred values inside the hole. */
template<typename TImage>
void BlurInHole(TImage* const image, const Mask* const mask, const float kernelVariance = 1.0f);

/** Median filter an image using all of its values but only replaced the pixel values
 * with the blurred values inside the hole. */
template<typename TImage>
void MedianFilterInHole(TImage* const image, const Mask* const mask, const unsigned int kernelRadius = 1);

/** Clip the values in the image inside the hole. */
template<typename TImage>
void ClipInHole(TImage* const image, const Mask* const mask, const float min, const float max);


/** Intersect a line with a hole. Return a pair containing the starting and ending
  * points of the line intersection. */
std::pair<itk::Index<2>, itk::Index<2> > IntersectLineWithHole(const std::vector<itk::Index<2> >& line,
                                                               const Mask* const mask,
                                                               bool &hasInteriorLine);

} // end namespace

#include "MaskOperations.hpp"

#endif
