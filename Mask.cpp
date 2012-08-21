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

#include "Mask.h"

// Submodules
#include <Helpers/Helpers.h>
#include <ITKHelpers/ITKHelpers.h>

// ITK
#include "itkBinaryContourImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkFlatStructuringElement.h"
#include "itkInvertIntensityImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkLabelShapeKeepNObjectsImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"

Mask::Mask() : HoleValue(255), ValidValue(0)
{

}

void Mask::Read(const std::string& filename)
{
  std::string extension = Helpers::GetFileExtension(filename);
  if(extension != "mask")
  {
    std::stringstream ss;
    ss << "Cannot read any file except .mask! Specified file was ." << extension;
    throw std::runtime_error(ss.str());
  }

  //Create an input stream for file
  std::ifstream fin(filename.c_str());

  if(!fin )
    {
    throw std::runtime_error("File not found!");
    }

  std::string line;
  std::stringstream linestream;

  int holeValue;
  int validValue;

  std::string imageFileName;
  getline(fin, line);

  linestream.clear();
  linestream << line;
  // Can't do this directly because HoleValue and ValidValue are unsigned char,
  // which will only read one character.
//   linestream >> this->HoleValue;
//   linestream >> this->ValidValue;
  linestream >> holeValue;
  linestream >> validValue;
  linestream >> imageFileName;

  this->HoleValue = holeValue;
  this->ValidValue = validValue;

  std::string path = Helpers::GetPath(filename);

  std::cout << "Reading mask: HoleValue " << static_cast<int>(this->HoleValue)
            << " ValidValue: " << static_cast<int>(this->ValidValue) << std::endl;

  std::string fullImageFileName = path + imageFileName;

  ReadFromImage(fullImageFileName);
}

void Mask::ReadFromImage(const std::string& filename)
{
  std::cout << "Reading mask from image: " << filename << std::endl;

  // Ensure the input image can be interpreted as a mask.
  {
  typedef itk::VectorImage<float, 2> TestImageType;
  typedef  itk::ImageFileReader<TestImageType> ImageReaderType;
  ImageReaderType::Pointer imageReader = ImageReaderType::New();
  imageReader->SetFileName(filename);
  imageReader->Update();

  unsigned int numberOfComponents = imageReader->GetOutput()->GetNumberOfComponentsPerPixel();
  if(!(numberOfComponents == 1 || numberOfComponents == 3))
    {
    std::stringstream ss;
    ss << "Number of components for a mask must be 1 or 3! (" << filename << " is " << numberOfComponents << ")";
    throw std::runtime_error(ss.str());
    }
  }

  // Should probably check that all 3 components are the same for all pixels (if numberOfComponents == 3)
  typedef itk::ImageFileReader<Mask> ImageReaderType;
  ImageReaderType::Pointer imageReader = ImageReaderType::New();
  imageReader->SetFileName(filename);
  imageReader->Update();

  unsigned char tempHoleValue = this->HoleValue;
  unsigned char tempValidValue = this->ValidValue;

  DeepCopyFrom(imageReader->GetOutput());

  this->HoleValue = tempHoleValue;
  this->ValidValue = tempValidValue;
}

unsigned int Mask::CountBoundaryPixels(const itk::ImageRegion<2>& region) const
{
  itk::ImageRegionConstIteratorWithIndex<Mask> maskIterator(this, region);
  unsigned int numberOfBoundaryPixels = 0;
  while(!maskIterator.IsAtEnd())
    {
    //if(this->IsHole(maskIterator.GetIndex()))
    if(maskIterator.Get() == this->HoleValue)
    {
      if(this->HasValidNeighbor(maskIterator.GetIndex()))
        {
        numberOfBoundaryPixels++;
        }
    }

    ++maskIterator;
    }
  return numberOfBoundaryPixels;
}

unsigned int Mask::CountBoundaryPixels() const
{
  return CountBoundaryPixels(this->GetLargestPossibleRegion());
}

unsigned int Mask::CountHolePixels(const itk::ImageRegion<2>& region) const
{
  return GetHolePixelsInRegion(region).size();
}

bool Mask::HasValidPixels() const
{
  return HasValidPixels(this->GetLargestPossibleRegion());
}

bool Mask::HasValidPixels(const itk::ImageRegion<2>& region) const
{
  if(CountValidPixels(region) > 0)
  {
    return true;
  }
  else
  {
    return false;
  }
}

bool Mask::HasHolePixels() const
{
  return HasHolePixels(this->GetLargestPossibleRegion());
}

bool Mask::HasHolePixels(const itk::ImageRegion<2>& region) const
{
  if(CountHolePixels(region) > 0)
  {
    return true;
  }
  else
  {
    return false;
  }
}

std::vector<itk::Index<2> > Mask::GetHolePixels() const
{
  return GetHolePixelsInRegion(this->GetLargestPossibleRegion());
}

unsigned int Mask::CountHolePixels() const
{
  return CountHolePixels(this->GetLargestPossibleRegion());
}

unsigned int Mask::CountValidPixels(const itk::ImageRegion<2>& region) const
{
  return GetValidPixelsInRegion(region).size();
}

unsigned int Mask::CountValidPixels() const
{
  return CountValidPixels(this->GetLargestPossibleRegion());
}

std::vector<itk::Offset<2> > Mask::GetValidOffsetsInRegion(itk::ImageRegion<2> region) const
{
  region.Crop(this->GetLargestPossibleRegion());

  std::vector<itk::Offset<2> > validOffsets;

  itk::ImageRegionConstIterator<Mask> iterator(this, region);

  while(!iterator.IsAtEnd())
    {
    //if(this->IsValid(iterator.GetIndex()))
    if(iterator.Get() == this->ValidValue)
      {
      validOffsets.push_back(iterator.GetIndex() - region.GetIndex());
      }

    ++iterator;
    }
  return validOffsets;
}

std::vector<itk::Offset<2> > Mask::GetHoleOffsetsInRegion(itk::ImageRegion<2> region) const
{
  region.Crop(this->GetLargestPossibleRegion());

  std::vector<itk::Offset<2> > holeOffsets;

  itk::ImageRegionConstIterator<Mask> iterator(this, region);

  while(!iterator.IsAtEnd())
    {
    //if(this->IsHole(iterator.GetIndex()))
    if(iterator.Get() == this->HoleValue)
      {
      holeOffsets.push_back(iterator.GetIndex() - region.GetIndex());
      }

    ++iterator;
    }
  return holeOffsets;
}

std::vector<itk::Index<2> > Mask::GetValidPixels(const bool forward) const
{
  return GetValidPixelsInRegion(this->GetLargestPossibleRegion(), forward);
}

std::vector<itk::Index<2> > Mask::GetValidPixelsInRegion(itk::ImageRegion<2> region,
                                                         const bool forward) const
{
  region.Crop(this->GetLargestPossibleRegion());

  std::vector<itk::Index<2> > validPixels;

  itk::ImageRegionConstIterator<Mask> iterator(this, region);

  while(!iterator.IsAtEnd())
    {
//    if(this->IsValid(iterator.GetIndex()))
//      {
//      validPixels.push_back(iterator.GetIndex());
//      }

    if(iterator.Get() == this->ValidValue)
      {
      validPixels.push_back(iterator.GetIndex());
      }

    ++iterator;
    }

  if(!forward)
  {
    std::reverse(validPixels.begin( ), validPixels.end( ) );
  }

  return validPixels;
}

std::vector<itk::Index<2> > Mask::GetHolePixelsInRegion(itk::ImageRegion<2> region) const
{
  region.Crop(this->GetLargestPossibleRegion());

  std::vector<itk::Index<2> > holePixels;

  itk::ImageRegionConstIterator<Mask> iterator(this, region);

  while(!iterator.IsAtEnd())
    {
    // Don't do this, because it loses the coherency of using an iterator.
//    if(this->IsHole(iterator.GetIndex()))
//      {
//      holePixels.push_back(iterator.GetIndex());
//      }
    if(iterator.Get() == this->HoleValue)
      {
      holePixels.push_back(iterator.GetIndex());
      }
    ++iterator;
    }
  return holePixels;
}

bool Mask::IsHole(const itk::Index<2>& index) const
{
  if(this->GetPixel(index) == this->HoleValue)
    {
    return true;
    }
  return false;
}

bool Mask::IsHole(const itk::ImageRegion<2>& region) const
{
  // If any of the pixels in the region are not hole pixels, the region is not entirely hole pixels.

  itk::ImageRegionConstIterator<Mask> maskIterator(this, region);

  while(!maskIterator.IsAtEnd())
    {
    //if(!this->IsHole(maskIterator.GetIndex()))
    if(maskIterator.Get() != this->HoleValue)
      {
      return false;
      }

    ++maskIterator;
    }
  return true;
}

bool Mask::IsValid(const itk::ImageRegion<2>& region) const
{
  // If any of the pixels in the region are invalid, the region is invalid.

  itk::ImageRegionConstIteratorWithIndex<Mask> maskIterator(this, region);

  while(!maskIterator.IsAtEnd())
    {
    //if(!this->IsValid(maskIterator.GetIndex()))
    if(maskIterator.Get() != this->ValidValue)
      {
      //std::cout << "Mask::IsValid - Pixel " << maskIterator.GetIndex() << " has value " << static_cast<unsigned int>(maskIterator.Get())
      //          << " which makes the region invalid because Mask::ValidValue = " << static_cast<unsigned int>(this->ValidValue) << std::endl;
      return false;
      }

    ++maskIterator;
    }
  return true;
}

bool Mask::IsValid(const itk::Index<2>& index) const
{
  if(this->GetPixel(index) == this->ValidValue)
    {
    return true;
    }
  return false;
}

void Mask::InvertInterpretation()
{
  unsigned char oldHoleValue = this->HoleValue;

  this->HoleValue = this->ValidValue;

  this->ValidValue = oldHoleValue;
}

void Mask::InvertData()
{
  // Exchange HoleValue and ValidValue, but leave everything else alone.
  itk::ImageRegionIterator<Mask> maskIterator(this, this->GetLargestPossibleRegion());
  unsigned int invertedCounter = 0;
  while(!maskIterator.IsAtEnd())
    {
    if(this->IsValid(maskIterator.GetIndex()))
      {
      maskIterator.Set(this->HoleValue);
      invertedCounter++;
      }
    else if(this->IsHole(maskIterator.GetIndex()))
      {
      maskIterator.Set(this->ValidValue);
      invertedCounter++;
      }
    ++maskIterator;
    }
  //std::cout << "Inverted " << invertedCounter << " in the mask." << std::endl;
}

void Mask::Cleanup()
{
  // We want to interpret pixels that are "pretty much hole value" as holes, and pixels that
  // are "pretty much valid value" as valid. The "do not use" pixels must be very far away from both of these values.
  itk::ImageRegionIterator<Mask> maskIterator(this, this->GetLargestPossibleRegion());

  float tolerance = 4;
  while(!maskIterator.IsAtEnd())
    {
    if(abs(maskIterator.Get() - this->ValidValue) < tolerance)
      {
      //std::cout << "Setting valid pixel to " << static_cast<unsigned int>(this->ValidValue) << std::endl;
      maskIterator.Set(this->ValidValue);
      }
    else if(abs(maskIterator.Get() - this->HoleValue) < tolerance)
      {
      //std::cout << "Setting hole pixel to " << static_cast<unsigned int>(this->HoleValue) << std::endl;
      maskIterator.Set(this->HoleValue);
      }
    ++maskIterator;
    }

}

void Mask::SetHoleValue(const unsigned char value)
{
  this->HoleValue = value;
}

void Mask::SetValidValue(const unsigned char value)
{
  this->ValidValue = value;
}

unsigned char Mask::GetHoleValue() const
{
  return this->HoleValue;
}

unsigned char Mask::GetValidValue() const
{
  return this->ValidValue;
}

void Mask::OutputMembers() const
{
  std::cout << "HoleValue: " << static_cast<unsigned int>(this->HoleValue) << std::endl;
  std::cout << "ValidValue: " << static_cast<unsigned int>(this->ValidValue) << std::endl;
}

void Mask::CopyInformationFrom(const Mask* const inputMask)
{
  assert(inputMask);
  this->SetHoleValue(inputMask->GetHoleValue());
  this->SetValidValue(inputMask->GetValidValue());
}

void Mask::DeepCopyFrom(const Mask* const inputMask)
{
  assert(inputMask);
  this->SetRegions(inputMask->GetLargestPossibleRegion());
  this->Allocate();

  itk::ImageRegionConstIterator<Mask> inputIterator(inputMask, inputMask->GetLargestPossibleRegion());
  itk::ImageRegionIterator<Mask> thisIterator(this, this->GetLargestPossibleRegion());

  while(!inputIterator.IsAtEnd())
    {
    thisIterator.Set(inputIterator.Get());
    ++inputIterator;
    ++thisIterator;
    }
  this->SetHoleValue(inputMask->GetHoleValue());
  this->SetValidValue(inputMask->GetValidValue());
}

void Mask::CopyHolesFrom(const Mask* const inputMask)
{
  itk::ImageRegionConstIterator<Mask> inputIterator(inputMask, inputMask->GetLargestPossibleRegion());
  itk::ImageRegionIterator<Mask> thisIterator(this, this->GetLargestPossibleRegion());

  while(!inputIterator.IsAtEnd())
    {
    if(inputMask->IsHole(inputIterator.GetIndex()))
      {
      thisIterator.Set(this->HoleValue);
      }
    ++inputIterator;
    ++thisIterator;
    }
}

void Mask::ExpandHole(const unsigned int kernelRadius)
{
  UnsignedCharImageType::Pointer binaryHoleImage = UnsignedCharImageType::New();
  this->CreateBinaryImage(binaryHoleImage, 255, 0);

//   std::cout << "binaryHoleImage: " << std::endl;
//   ITKHelpers::PrintImage(binaryHoleImage.GetPointer());

  typedef itk::FlatStructuringElement<2> StructuringElementType;
  StructuringElementType::RadiusType radius;
  radius.Fill(kernelRadius); // This is correct that the RadiusType expects the region radius, not the side length.

  StructuringElementType structuringElement = StructuringElementType::Box(radius);
  typedef itk::BinaryDilateImageFilter<UnsignedCharImageType, UnsignedCharImageType, StructuringElementType> BinaryDilateImageFilterType;
  BinaryDilateImageFilterType::Pointer dilateFilter = BinaryDilateImageFilterType::New();
  dilateFilter->SetInput(binaryHoleImage);
  dilateFilter->SetKernel(structuringElement);
  dilateFilter->Update();

//   std::cout << "dilateFilter output: " << std::endl;
//   ITKHelpers::PrintImage(dilateFilter->GetOutput());

  // There will now be more hole pixels than there were previously. Copy them into the mask.
  this->CopyHolesFromValue(dilateFilter->GetOutput(), 255);
}

void Mask::ShrinkHole(const unsigned int kernelRadius)
{
  UnsignedCharImageType::Pointer binaryHoleImage = UnsignedCharImageType::New();
  this->CreateBinaryImage(binaryHoleImage, 255, 0);

//   std::cout << "binaryHoleImage: " << std::endl;
//   ITKHelpers::PrintImage(binaryHoleImage.GetPointer());

  typedef itk::FlatStructuringElement<2> StructuringElementType;
  StructuringElementType::RadiusType radius;
  radius.Fill(kernelRadius); // This is correct that the RadiusType expects the region radius, not the side length.

  StructuringElementType structuringElement = StructuringElementType::Box(radius);
  typedef itk::BinaryErodeImageFilter<UnsignedCharImageType, UnsignedCharImageType, StructuringElementType> BinaryErodeImageFilterType;
  BinaryErodeImageFilterType::Pointer erodeFilter = BinaryErodeImageFilterType::New();
  erodeFilter->SetInput(binaryHoleImage);
  erodeFilter->SetKernel(structuringElement);
  erodeFilter->Update();

//   std::cout << "erodeFilter output: " << std::endl;
//   ITKHelpers::PrintImage(erodeFilter->GetOutput());

  // There will now be more valid pixels than there were previously. Copy them into the mask.
  this->CopyValidPixelsFromValue(erodeFilter->GetOutput(), 0);
}

void Mask::CreateImage(UnsignedCharImageType* const image, const unsigned char holeColor,
                            const unsigned char validColor, const unsigned char otherColor)
{
  image->SetRegions(this->GetLargestPossibleRegion());
  image->Allocate();

  itk::ImageRegionIterator<UnsignedCharImageType> binaryImageIterator(image,
                                                                      image->GetLargestPossibleRegion());

  while(!binaryImageIterator.IsAtEnd())
    {
    if(this->IsHole(binaryImageIterator.GetIndex()))
      {
      binaryImageIterator.Set(holeColor);
      }
    else if(this->IsValid(binaryImageIterator.GetIndex()))
      {
      binaryImageIterator.Set(validColor);
      }
    else
      {
      binaryImageIterator.Set(otherColor);
      }
    ++binaryImageIterator;
    }
}

void Mask::CreateBinaryImage(UnsignedCharImageType* const image, const unsigned char holeColor,
                            const unsigned char validColor)
{
  CreateImage(image, holeColor, validColor, validColor);
}

void Mask::FindBoundaryInRegion(const itk::ImageRegion<2>& region, BoundaryImageType* const boundaryImage,
                                const PixelTypeEnum& whichSideOfBoundary,
                                const BoundaryImageType::PixelType& outputBoundaryPixelValue) const
{
  Mask::Pointer extractedRegionMask = Mask::New();
  extractedRegionMask->SetRegions(region);
  extractedRegionMask->Allocate();
  ITKHelpers::CopyRegion(this, extractedRegionMask.GetPointer(), region, region);

  //HelpersOutput::WriteImageConditional<Mask>(holeOnly, "Debug/FindBoundary.HoleOnly.mha", this->DebugImages);
  //HelpersOutput::WriteImageConditional<Mask>(holeOnly, "Debug/FindBoundary.HoleOnly.png", this->DebugImages);

  // Since the hole is white, we want the foreground value of the contour filter to be black.
  // This means that the boundary will
  // be detected in the black pixel region, which is on the outside edge of the hole like we want. However,
  // The BinaryContourImageFilter will change all non-boundary pixels to the background color,
  // so the resulting output will be inverted - the boundary pixels will be black and the
  // non-boundary pixels will be white.

  // Find the boundary
  typedef itk::BinaryContourImageFilter<Mask, Mask> binaryContourImageFilterType;
  binaryContourImageFilterType::Pointer binaryContourFilter = binaryContourImageFilterType::New();
  binaryContourFilter->SetInput(this);
  binaryContourFilter->SetFullyConnected(true);

  unsigned char foregroundValue;
  unsigned char backgroundValue;
  if(whichSideOfBoundary == VALID)
  {
    // we want the boundary pixels to be in the valid region.
    foregroundValue = this->GetValidValue();
    backgroundValue = this->GetHoleValue();
  }
  else if(whichSideOfBoundary == HOLE)
  {
    // we want the boundary pixels to be in the hole region.
    foregroundValue = this->GetHoleValue();
    backgroundValue = this->GetValidValue();
  }
  else
  {
    throw std::runtime_error("An invalid side of the boundary was requested.");
  }

  binaryContourFilter->SetForegroundValue(foregroundValue);
  binaryContourFilter->SetBackgroundValue(backgroundValue);
  binaryContourFilter->Update();

  binaryContourFilter->GetOutput()->CopyInformationFrom(this);


  if(whichSideOfBoundary == VALID)
  {
    // CreateBinaryImage(holeColor, validColor)
    binaryContourFilter->GetOutput()->CreateBinaryImage(boundaryImage, 0, outputBoundaryPixelValue);
  }
  else if(whichSideOfBoundary == HOLE)
  {
    binaryContourFilter->GetOutput()->CreateBinaryImage(boundaryImage, outputBoundaryPixelValue, 0);
  }

}

void Mask::FindBoundary(itk::Image<unsigned char, 2>* const boundaryImage, const PixelTypeEnum& whichSideOfBoundary,
                        const BoundaryImageType::PixelType& outputBoundaryPixelValue) const
{
  FindBoundaryInRegion(this->GetLargestPossibleRegion(), boundaryImage,
                       whichSideOfBoundary, outputBoundaryPixelValue);
}

/** Get a list of the valid neighbors of a pixel.*/
std::vector<itk::Index<2> > Mask::GetValidNeighbors(const itk::Index<2>& pixel) const
{
  return ITKHelpers::Get8NeighborsWithValue(pixel, this, this->ValidValue);
}

bool Mask::HasHoleNeighbor(const itk::Index<2>& pixel) const
{
  if(GetHoleNeighbors(pixel).size() > 0)
    {
    return true;
    }
  return false;
}

bool Mask::HasValidNeighbor(const itk::Index<2>& pixel) const
{
  if(GetValidNeighbors(pixel).size() > 0)
    {
    return true;
    }
  return false;
}

/** Get a list of the hole neighbors of a pixel.*/
std::vector<itk::Index<2> > Mask::GetHoleNeighbors(const itk::Index<2>& pixel) const
{
  return ITKHelpers::Get8NeighborsWithValue(pixel, this, this->HoleValue);
}

std::vector<itk::Offset<2> > Mask::GetValidNeighborOffsets(const itk::Index<2>& pixel) const
{
  std::vector<itk::Index<2> > indices = ITKHelpers::Get8NeighborsWithValue(pixel, this, this->ValidValue);
  std::vector<itk::Offset<2> > offsets;
  for(unsigned int i = 0; i < indices.size(); ++i)
  {
    offsets.push_back(indices[i] - pixel);
  }
  return offsets;
}

std::vector<itk::Offset<2> > Mask::GetHoleNeighborOffsets(const itk::Index<2>& pixel) const
{
  std::vector<itk::Index<2> > indices = ITKHelpers::Get8NeighborsWithValue(pixel, this, this->HoleValue);
  std::vector<itk::Offset<2> > offsets;
  for(unsigned int i = 0; i < indices.size(); ++i)
  {
    offsets.push_back(indices[i] - pixel);
  }
  return offsets;
}

void Mask::MarkAsHole(const itk::Index<2>& pixel)
{
  this->SetPixel(pixel, this->HoleValue);
}

void Mask::MarkAsValid(const itk::Index<2>& pixel)
{
  this->SetPixel(pixel, this->ValidValue);
}

bool Mask::HasValid4Neighbor(const itk::Index<2>& pixel)
{
  std::vector<itk::Index<2> > neighbors =
        ITKHelpers::Get4NeighborIndicesInsideRegion(pixel, this->GetLargestPossibleRegion());

  for(unsigned int i = 0; i < neighbors.size(); ++i)
    {
    if(this->IsValid(neighbors[i]))
      {
      return true;
      }
    }

  return false;
}

std::vector<itk::Index<2> > Mask::GetValid4Neighbors(const itk::Index<2>& pixel)
{
  std::vector<itk::Index<2> > neighborhood =
          ITKHelpers::Get4NeighborIndicesInsideRegion(pixel, this->GetLargestPossibleRegion());

  std::vector<itk::Index<2> > validNeighbors;

  for(unsigned int i = 0; i < neighborhood.size(); ++i)
    {
    if(this->IsValid(neighborhood[i]))
      {
      validNeighbors.push_back(neighborhood[i]);
      }
    }

  return validNeighbors;
}

bool Mask::IsHoleValue(const unsigned char value) const
{
  if(value == this->HoleValue)
  {
    return true;
  }
  else
  {
    return false;
  }
}


bool Mask::IsValidValue(const unsigned char value) const
{
  if(value == this->ValidValue)
  {
    return true;
  }
  else
  {
    return false;
  }
}

void Mask::KeepLargestHole()
{
  // Only keep the largest segment
  typedef itk::ConnectedComponentImageFilter<Mask, Mask> ConnectedComponentImageFilterType;
  ConnectedComponentImageFilterType::Pointer connectedComponentFilter = ConnectedComponentImageFilterType::New ();
  connectedComponentFilter->SetInput(this);
  connectedComponentFilter->Update();

  typedef itk::LabelShapeKeepNObjectsImageFilter<Mask> LabelShapeKeepNObjectsImageFilterType;
  LabelShapeKeepNObjectsImageFilterType::Pointer labelShapeKeepNObjectsImageFilter =
           LabelShapeKeepNObjectsImageFilterType::New();
  labelShapeKeepNObjectsImageFilter->SetInput(connectedComponentFilter->GetOutput());
  labelShapeKeepNObjectsImageFilter->SetBackgroundValue(0);
  labelShapeKeepNObjectsImageFilter->SetNumberOfObjects(1);
  labelShapeKeepNObjectsImageFilter
            ->SetAttribute(LabelShapeKeepNObjectsImageFilterType::LabelObjectType::NUMBER_OF_PIXELS);
  labelShapeKeepNObjectsImageFilter->Update();

  typedef itk::RescaleIntensityImageFilter<Mask, Mask> RescaleFilterType;
  RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
  rescaleFilter->SetOutputMinimum(0);
  rescaleFilter->SetOutputMaximum(255);
  rescaleFilter->SetInput(labelShapeKeepNObjectsImageFilter->GetOutput());
  rescaleFilter->Update();

  ITKHelpers::DeepCopy(rescaleFilter->GetOutput(), this);
}

unsigned int Mask::CountValidPatches(const unsigned int patchRadius) const
{

  itk::ImageRegionConstIteratorWithIndex<Mask> maskIterator(this, this->GetLargestPossibleRegion());

  unsigned int counter = 0;
  // std::cout << "CountValidPatches (patch radius " << patchRadius << ")..." << std::endl;
  while(!maskIterator.IsAtEnd())
    {
    itk::ImageRegion<2> region = ITKHelpers::GetRegionInRadiusAroundPixel(maskIterator.GetIndex(), patchRadius);

    if(this->IsValid(region))
      {
      counter++;
      }
    ++maskIterator;
    }

  //std::cout << "There were " << counter << " valid patches." << std::endl;
  return counter;
}

itk::ImageRegion<2> Mask::FindFirstValidPatch(const unsigned int patchRadius)
{
  itk::ImageRegionConstIteratorWithIndex<Mask> maskIterator(this, this->GetLargestPossibleRegion());

  while(!maskIterator.IsAtEnd())
    {
    itk::ImageRegion<2> region = ITKHelpers::GetRegionInRadiusAroundPixel(maskIterator.GetIndex(), patchRadius);

    if(this->IsValid(region))
      {
      return region;
      }

    ++maskIterator;
    }

  throw std::runtime_error("No valid patches found!");

  // We should never reach this point
  itk::ImageRegion<2> dummyRegion;
  return dummyRegion;
}

void Mask::SetHole(const itk::Index<2>& index)
{
  this->SetPixel(index, this->HoleValue);
}

void Mask::SetValid(const itk::Index<2>& index)
{
  this->SetPixel(index, this->ValidValue);
}

void Mask::SetValid(const itk::ImageRegion<2>& region)
{
  itk::ImageRegionIteratorWithIndex<Mask> maskIterator(this, region);

  while(!maskIterator.IsAtEnd())
    {
    maskIterator.Set(this->ValidValue);
    ++maskIterator;
    }
}
