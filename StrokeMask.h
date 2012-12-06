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
\class StrokeMask
\brief This class is a subclass of itkImage that provides the concept of a "stroke"
       drawn on a background.
*/

#ifndef StrokeMask_H
#define StrokeMask_H

// ITK
#include "itkImage.h"

/** The pixels in the mask have only these possible values. */
enum class StrokeMaskPixelTypeEnum {STROKE, NOTSTROKE};

/** This must be defined in order to create an itk::Image<StrokeMaskPixelTypeEnum> because
  * The Set/Get macros require a way to output the pixel type. */
std::ostream& operator<<(std::ostream& output, const StrokeMaskPixelTypeEnum &pixelType);

class StrokeMask : public itk::Image<StrokeMaskPixelTypeEnum, 2>
{
public:
  /** Standard typedefs. */
  typedef StrokeMask                       Self;
  typedef itk::Image< unsigned char, 2> Superclass;
  typedef itk::SmartPointer< Self >              Pointer;
  typedef itk::SmartPointer< const Self >        ConstPointer;
  typedef itk::WeakPointer< const Self >         ConstWeakPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(StrokeMask, Image);

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

  /** Determine if a pixel is a stroke pixel.*/
  bool IsStroke(const itk::Index<2>& index) const;

  template <typename TPixel>
  void ReadFromImage(const std::string& filename, const TPixel& strokevalue);

  template <typename TPixel>
  void Write(const std::string& filename, const TPixel& strokeValue);

  /** Read the mask from a .mask file.*/
  void Read(const std::string& filename);

  /** Count foreground pixels in the whole mask.*/
  unsigned int CountStrokePixels() const;

private:

  StrokeMask(const Self &);    //purposely not implemented
  void operator=(const Self &); //purposely not implemented

  StrokeMask(){} // required by itkNewMacro
};

#include "StrokeMask.hpp"

#endif
