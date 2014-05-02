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
\class ForegroundBackground
\brief This class is a subclass of itkImage that provides the concept of "valid" pixels
       and hole pixels. Pixels that are any other value are never used in computations.
       Using itkImageFileReader, the first channel of any input image will be attempted
       to be converted to a Mask. If the image is a 4 channel image where the 4th
       channel represents alpha (or a 2 channel image where the 2nd channel represents alpha)
       the reader sometimes produces a blank image. We throw an exception if this is the case.
       Your mask should be a 1 or 3 channel image.
*/

#ifndef ForegroundBackgroundSegmentMask_H
#define ForegroundBackgroundSegmentMask_H

// ITK
#include "itkImage.h"

/** The pixels in the mask have only these possible values. */
enum class ForegroundBackgroundSegmentMaskPixelTypeEnum {FOREGROUND, BACKGROUND};

/** This must be defined in order to create an itk::Image<HoleMaskPixelTypeEnum> because
  * The Set/Get macros require a way to output the pixel type. */
std::ostream& operator<<(std::ostream& output, const ForegroundBackgroundSegmentMaskPixelTypeEnum &pixelType);


/** This class forces us to pass functions values as ForegroundValueWrapper(0) instead of just "0"
  * so that we can be sure that a foreground value is getting passed where a foreground value is expected,
  * and not accidentally confuse the order of foreground/background arguments silently. */
template <typename T>
struct ForegroundPixelValueWrapper
{
  ForegroundPixelValueWrapper(T value) : Value(value){}

  operator T()
  {
    return this->Value;
  }

  T Value;
};

/** This class forces us to pass functions values as BackgroundValueWrapper(0) instead of just "0"
  * so that we can be sure that a background value is getting passed where a background value is expected,
  * and not accidentally confuse the order of foreground/background arguments silently. */
template <typename T>
struct BackgroundPixelValueWrapper
{
  BackgroundPixelValueWrapper(T value) : Value(value){}

  operator T()
  {
    return this->Value;
  }

  T Value;
};


class ForegroundBackgroundSegmentMask : public itk::Image<ForegroundBackgroundSegmentMaskPixelTypeEnum, 2>
{
public:
  /** Standard typedefs. */
  typedef ForegroundBackgroundSegmentMask                       Self;
  typedef itk::Image< unsigned char, 2> Superclass;
  typedef itk::SmartPointer< Self >              Pointer;
  typedef itk::SmartPointer< const Self >        ConstPointer;
  typedef itk::WeakPointer< const Self >         ConstWeakPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ForegroundBackgroundSegmentMask, Image);

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

  /** Determine if a pixel is a foreground pixel.*/
  bool IsForeground(const itk::Index<2>& index) const;

  /** Determine if a pixel is a background pixel.*/
  bool IsBackground(const itk::Index<2>& index) const;

  template <typename TPixel>
  void ReadFromImage(const std::string& filename, const ForegroundPixelValueWrapper<TPixel>& foregroundValue,
                     const BackgroundPixelValueWrapper<TPixel>& backgroundValue);

  template <typename TPixel>
  void Write(const std::string& filename, const ForegroundPixelValueWrapper<TPixel>& foregroundValue,
             const BackgroundPixelValueWrapper<TPixel>& backgroundValue);

  template <typename TImage>
  void ApplyToImage(TImage* image,
                    const typename TImage::PixelType& backgroundValue);

  /** Read the mask from a .mask file.*/
  void Read(const std::string& filename);

  /** Count foreground pixels in the whole mask.*/
  unsigned int CountForegroundPixels() const;

  /** Count background pixels in the whole mask.*/
  unsigned int CountBackgroundPixels() const;

private:

  ForegroundBackgroundSegmentMask(const Self &);    //purposely not implemented
  void operator=(const Self &); //purposely not implemented

  ForegroundBackgroundSegmentMask(){} // required by itkNewMacro
};

#include "ForegroundBackgroundSegmentMask.hpp"

#endif
