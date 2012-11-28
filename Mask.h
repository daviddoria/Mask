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

/**
\class Mask
\brief This class is a subclass of itkImage that provides the concept of "valid" pixels
       and hole pixels. Pixels that are any other value are never used in computations.
       Using itkImageFileReader, the first channel of any input image will be attempted
       to be converted to a Mask. If the image is a 4 channel image where the 4th
       channel represents alpha (or a 2 channel image where the 2nd channel represents alpha)
       the reader sometimes produces a blank image. We throw an exception if this is the case.
       Your mask should be a 1 or 3 channel image.
*/

#ifndef MASK_H
#define MASK_H

// ITK
#include "itkImage.h"

/** The pixels in the mask have only these possible values. */
enum class HoleMaskPixelTypeEnum {HOLE, VALID, UNDETERMINED};

/** This must be defined in order to create an itk::Image<HoleMaskPixelTypeEnum> because
  * The Set/Get macros require a way to output the pixel type. */
std::ostream& operator<<(std::ostream& output, const HoleMaskPixelTypeEnum &pixelType);

/** This class forces us to pass functions values as HoleValueWrapper(0) instead of just "0"
  * so that we can be sure that a hole value is getting passed where a hole value is expected,
  * and not accidentally confuse the order of hole/valid arguments silently. */
template <typename T>
struct HolePixelValueWrapper
{
  HolePixelValueWrapper(T value) : Value(value){}

  operator T()
  {
    return this->Value;
  }

  T Value;
};

/** This class forces us to pass functions values as ValidValueWrapper(0) instead of just "0"
  * so that we can be sure that a value value is getting passed where a value value is expected,
  * and not accidentally confuse the order of hole/valid arguments silently. */
template <typename T>
struct ValidPixelValueWrapper
{
  ValidPixelValueWrapper(T value) : Value(value){}

  operator T()
  {
    return this->Value;
  }

  T Value;
};

class Mask : public itk::Image<HoleMaskPixelTypeEnum, 2>
{
public:
  /** Standard typedefs. */
  typedef Mask                       Self;
  typedef itk::Image< unsigned char, 2> Superclass;
  typedef itk::SmartPointer< Self >              Pointer;
  typedef itk::SmartPointer< const Self >        ConstPointer;
  typedef itk::WeakPointer< const Self >         ConstWeakPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(Mask, Image);

  /** Dimension of the image. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      Superclass::ImageDimension);

  /** Types derived from the Superclass */
  typedef Superclass::IndexType IndexType;

  typedef Superclass::IOPixelType IOPixelType;

  /** Tyepdef for the functor used to access a neighborhood of pixel
  * pointers. */
  typedef itk::NeighborhoodAccessorFunctor< Self >
  NeighborhoodAccessorFunctorType;

  /** Return the NeighborhoodAccessor functor. This method is called by the
   * neighborhood iterators. */
  NeighborhoodAccessorFunctorType GetNeighborhoodAccessor()
  { return NeighborhoodAccessorFunctorType(); }

  /** Return the NeighborhoodAccessor functor. This method is called by the
   * neighborhood iterators. */
  const NeighborhoodAccessorFunctorType GetNeighborhoodAccessor() const
  { return NeighborhoodAccessorFunctorType(); }

  /** Determine if a pixel is a hole pixel.*/
  bool IsHole(const itk::Index<2>& index) const;

  /** Determine if an entire region consists of hole pixels.*/
  bool IsHole(const itk::ImageRegion<2>& region) const;

  /** Determine if an entire region is valid.*/
  bool IsValid(const itk::ImageRegion<2>& region) const;

  /** Determine if a pixel is valid.*/
  bool IsValid(const itk::Index<2>& index) const;

  /** Create a binary image of holes and valid pixels.*/
  typedef itk::Image<unsigned char, 2> UnsignedCharImageType;
  void CreateBinaryImage(UnsignedCharImageType* const image, const unsigned char holeColor,
                         const unsigned char validColor);

  /** Create an image of holes, valid pixels, and undetermined pixels.*/
  void CreateImage(UnsignedCharImageType* const image, const unsigned char holeColor,
                   const unsigned char validColor, const unsigned char undeterminedColor);

  /** Invert the mask by setting all hole pixels to ValidValue and all valid pixels to HoleValue.*/
  void InvertData();

  /** Invert the mask by switching HoleValue and ValidValue.*/
  void InvertInterpretation();

  /** Snap the pixel values to either 'hole' or 'valid'.*/
  void Cleanup();

  /** Only keep the largest separate hole.*/
  void KeepLargestHole();

  /** Increase the size of the hole.*/
  void ExpandHole(const unsigned int kernelRadius);

  /** Decrease the size of the hole.*/
  void ShrinkHole(const unsigned int kernelRadius);

  /** Mark the pixel as a hole.*/
  void SetHole(const itk::Index<2>& index);

  /** Mark the pixel as valid.*/
  void SetValid(const itk::Index<2>& index);

  /** Mark the region as valid.*/
  void SetValid(const itk::ImageRegion<2>& region);

  /** Copy the holes from a mask.*/
  void CopyHolesFrom(const Mask* const inputMask);

  /** Create holes from specified pixels in an image.*/
  template <typename TImage>
  void CreateHolesFromValue(const TImage* const inputImage,
                            const typename TImage::PixelType value);

  /** Create valid pixels from specified pixels in an image.*/
  template <typename TImage>
  void CreateValidPixelsFromValue(const TImage* const inputImage,
                                  const typename TImage::PixelType value);

  /** Find the boundary of the Mask. The output image is black with boundary pixels of value
    * 'outputBoundaryPixelValue'. */
  typedef itk::Image<unsigned char, 2> BoundaryImageType;
  void CreateBoundaryImage(BoundaryImageType* const boundary, const HoleMaskPixelTypeEnum& whichSideOfBoundary,
                           const BoundaryImageType::PixelType& outputBoundaryPixelValue = 255) const;
  void CreateBoundaryImageInRegion(const itk::ImageRegion<2>& region, BoundaryImageType* const boundaryImage,
                                   const HoleMaskPixelTypeEnum& whichSideOfBoundary,
                                   const BoundaryImageType::PixelType& outputBoundaryPixelValue = 255) const;

  /** Recolor the hole pixels in 'image' a specified 'color'. 'color' cannot have a
    * default value (even itk::NumericTraits<T>::ZeroValue())
    * because for a itk::VectorImage<ScalarType, 2>, the ::PixelType is a
    * itk::VariableLengthVector<ScalarType>, which does not have a size,
    * so therefore the ZeroValue() function will not do what we'd expect. */
  template<typename TImage>
  void ApplyToImage(TImage* const image, const typename TImage::PixelType& color) const;

  /** Recolor the hole pixels in 'maskRegion' in 'imageRegion' in 'image' a specified 'color'.*/
  template<typename TImage>
  void ApplyRegionToImageRegion(const itk::ImageRegion<2>& maskRegion, TImage* const image,
                                const itk::ImageRegion<2>& imageRegion,
                                const typename TImage::PixelType& color) const;

  /** Change the hole pixels in 'image' to a specified 'holeValue'. 'holeValue' is not const because it might
    * need to be modified if it is not provided or is invalid. */
  template<typename TImage>
  void ApplyToScalarImage(TImage* const image,
                          typename TImage::PixelType holeValue =
                            itk::NumericTraits<typename TImage::PixelType>::ZeroValue()) const;

  /** Recolor the hole pixels in 'image' a specified 'color'.
   * Here 'TColor' must have .red(), .green(), and .blue() functions.*/
  template<typename TImage, typename TColor>
  void ApplyToRGBImage(TImage* const image, const TColor& color) const;

  /** Create a mask from a mask image. That is, take a binary image (or grayscale) and convert it to a Mask. */
  template<typename TImage>
  void CreateFromImage(const TImage* const image,
                       const HolePixelValueWrapper<typename TImage::PixelType>& holeValue,
                       const ValidPixelValueWrapper<typename TImage::PixelType>& validValue);

  /** Get a list of the valid neighbors of a pixel.*/
  std::vector<itk::Index<2> > GetValid8Neighbors(const itk::Index<2>& pixel) const;

  /** Get a list of the valid neighbors of a 'pixel' in 'region'.*/
  std::vector<itk::Index<2> > GetValid8NeighborsInRegion(const itk::Index<2>& pixel,
                                                        const itk::ImageRegion<2>& region) const;

  /** Determine if a pixel has at least 1 hole 8-neighbor.*/
  bool HasHole8Neighbor(const itk::Index<2>& pixel) const;

  /** Determine if a pixel has at least 1 hole 8-neighbor in 'region'.*/
  bool HasHole8NeighborInRegion(const itk::Index<2>& pixel, const itk::ImageRegion<2>& region) const;

  /** Determine if a pixel has at least 1 valid 8-neighbor.*/
  bool HasValid8Neighbor(const itk::Index<2>& pixel) const;

  /** Determine if a pixel has at least 1 valid 4-neighbor.*/
  bool HasValid4Neighbor(const itk::Index<2>& pixel);

  /** Determine which of the 4-neighbors of a 'pixel' are valid and inside 'region'.*/
  std::vector<itk::Index<2> > GetValid4NeighborIndices(const itk::Index<2>& pixel,
                                                       const itk::ImageRegion<2>& region);

  /** Determine which of the 4-neighbors of a 'pixel' are valid.*/
  std::vector<itk::Index<2> > GetValid4Neighbors(const itk::Index<2>& pixel);

  /** Get a list of the hole neighbors of a pixel.*/
  std::vector<itk::Index<2> > GetHole8Neighbors(const itk::Index<2>& pixel) const;

  /** Get a list of the hole neighbors of a 'pixel' in 'region'.*/
  std::vector<itk::Index<2> > GetHole8NeighborsInRegion(const itk::Index<2>& pixel,
                                                       const itk::ImageRegion<2>& region) const;

  /** Get a list of the offsets of the valid neighbors of a pixel.*/
  std::vector<itk::Offset<2> > GetValid8NeighborOffsets(const itk::Index<2>& pixel) const;

  /** Get a list of the offsets of the hole neighbors of a pixel.*/
  std::vector<itk::Offset<2> > GetHole8NeighborOffsets(const itk::Index<2>& pixel) const;

  /** Get a list of the valid pixels in a region. 'region' is not passed by reference because
    * it is cropped by the image before computing the valid pixels. */
  std::vector<itk::Index<2> > GetValidPixelsInRegion(itk::ImageRegion<2> region,
                                                     const bool forward = true) const;

  /** Get a list of the valid pixels in the mask.*/
  std::vector<itk::Index<2> > GetValidPixels(const bool forward = true) const;

  /** Get a list of the hole pixels in a region. 'region' is not passed by reference because
    * it is cropped by the image before computing the valid pixels. */
  std::vector<itk::Index<2> > GetHolePixelsInRegion(itk::ImageRegion<2> region) const;

  /** Get a list of the hole pixels in the mask.*/
  std::vector<itk::Index<2> > GetHolePixels() const;

  /** Get a list of the offsets of the valid pixels in a region.*/
  std::vector<itk::Offset<2> > GetValidOffsetsInRegion(itk::ImageRegion<2> region) const;

  /** Get a list of the offsets of the hole pixels in a region.*/
  std::vector<itk::Offset<2> > GetHoleOffsetsInRegion(itk::ImageRegion<2> region) const;

  /** Count hole pixels in a region.*/
  unsigned int CountHolePixels(const itk::ImageRegion<2>& region) const;

  /** Find hole pixels that are touching valid pixels.*/
  std::vector<itk::Index<2> > FindBoundaryPixelsInRegion(const itk::ImageRegion<2>& region,
                                                         const HoleMaskPixelTypeEnum& whichSideOfBoundary) const;

  /** Find hole pixels that are touching valid pixels.*/
  std::vector<itk::Index<2> > FindBoundaryPixels(const Mask::PixelType& whichSideOfBoundary) const;

  /** Count hole pixels that are touching valid pixels.*/
  unsigned int CountBoundaryPixels(const itk::ImageRegion<2>& region,
                                   const HoleMaskPixelTypeEnum& whichSideOfBoundary) const;

  /** Count hole pixels that are touching valid pixels.*/
  unsigned int CountBoundaryPixels(const Mask::PixelType& whichSideOfBoundary) const;

  /** Count hole pixels in the whole mask.*/
  unsigned int CountHolePixels() const;

  /** Determine if the mask has any valid pixels.*/
  bool HasValidPixels() const;

  /** Determine if the mask has any valid pixels in 'region'.*/
  bool HasValidPixels(const itk::ImageRegion<2>& region) const;

  /** Determine if the mask has any hole pixels.*/
  bool HasHolePixels() const;

  /** Determine if the mask has any hole pixels in 'region'.*/
  bool HasHolePixels(const itk::ImageRegion<2>& region) const;

  /** Count valid pixels in a region.*/
  unsigned int CountValidPixels(const itk::ImageRegion<2>& region) const;

  /** Count valid pixels in a region.*/
  unsigned int CountValidPatches(const unsigned int patchRadius) const;

  /** Find the first valid patch of radius 'patchRadius' in raster scan order.*/
  itk::ImageRegion<2> FindFirstValidPatch(const unsigned int patchRadius);

  /** Count valid pixels in the whole mask.*/
  unsigned int CountValidPixels() const;

  /** Read the mask from a .mask file.*/
  void Read(const std::string& filename);

  /** Read the mask from an image file.*/
  template <typename TPixel>
  void ReadFromImage(const std::string& filename, const HolePixelValueWrapper<TPixel>& holeValue,
                     const ValidPixelValueWrapper<TPixel>& validValue);

  /** Mark the pixel as a hole.*/
  void MarkAsHole(const itk::Index<2>& pixel);

  /** Mark the pixel as a valid pixel.*/
  void MarkAsValid(const itk::Index<2>& pixel);

private:

  Mask(const Self &);    //purposely not implemented
  void operator=(const Self &); //purposely not implemented

  /** Constructor. */
  Mask();
};

#include "Mask.hpp"

#endif
