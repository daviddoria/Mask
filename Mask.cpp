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

std::string Mask::GetFilenameFromMaskFile(const std::string& maskFileName)
{
  // Currently this is massive duplication between this function and Read().
  // We should change Read to take advantage of a function to similar to this.
  /**
   * The format of the .mask file is:
   * hole 0
   * valid 255
   * Mask.png
   *
   * OR
   *
   * valid 255
   * hole 0
   * Mask.png
   *
   * That is, the "valid VALUE" line can be either on the first or second line.
   * Note that the 0 and 255 here are arbitrary and can be anything.
   */
  std::string extension = Helpers::GetFileExtension(maskFileName);
  if(extension != "mask")
  {
    std::stringstream ss;
    ss << "Cannot read any file except .mask! Specified file was ." << extension
       << " You might want ReadFromImage instead.";
    throw std::runtime_error(ss.str());
  }

  //Create an input stream for file
  std::ifstream fin(maskFileName.c_str());

  if(!fin )
  {
    throw std::runtime_error("File not found!");
  }

  std::string line;

  getline(fin, line);

  getline(fin, line);

  std::string imageFileName;

  getline(fin, imageFileName);
  std::cout << "Mask image file: " << imageFileName << std::endl;
  return imageFileName;
}

void Mask::Read(const std::string& filename)
{
  /**
   * The format of the .mask file is:
   * hole 0
   * valid 255
   * Mask.png
   *
   * OR
   *
   * valid 255
   * hole 0
   * Mask.png
   *
   * That is, the "valid VALUE" line can be either on the first or second line.
   * Note that the 0 and 255 here are arbitrary and can be anything.
   */
  std::string extension = Helpers::GetFileExtension(filename);
  if(extension != "mask")
  {
    std::stringstream ss;
    ss << "Cannot read any file except .mask! Specified file was ." << extension
       << " You might want ReadFromImage instead.";
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
  int holeValue = 0;
  int validValue = 0;

  std::string type1;
  std::string type2;
  int value1 = 0;
  int value2 = 0;

  getline(fin, line);
  linestream.clear();
  linestream << line;
  linestream >> type1 >> value1;

  if(type1 == "hole")
  {
    holeValue = value1;
  }
  else if(type1 == "valid")
  {
    validValue = value1;
  }
  else
  {
    throw std::runtime_error("Invalid .mask file!");
  }

  getline(fin, line);
  linestream.clear();
  linestream << line;
  linestream >> type2 >> value2;

  if(type1 == type2)
  {
    throw std::runtime_error("Invalid .mask file! Valid or hole value listed twice!");
  }

  if(type2 == "hole")
  {
    holeValue = value2;
  }
  else if(type2 == "valid")
  {
    validValue = value2;
  }
  else
  {
    throw std::runtime_error("Invalid .mask file!");
  }

  std::string imageFileName;
//  linestream >> imageFileName;
  getline(fin, imageFileName);
  std::cout << "Mask image file: " << imageFileName << std::endl;

  std::string path = Helpers::GetPath(filename);

  std::string fullImageFileName = path + imageFileName;

  std::cout << "Full mask image file: " << fullImageFileName << std::endl;

  ReadFromImage(fullImageFileName, HolePixelValueWrapper<int>(holeValue),
                ValidPixelValueWrapper<int>(validValue));
}

unsigned int Mask::CountBoundaryPixels(const itk::ImageRegion<2>& region,
                                       const Mask::PixelType& whichSideOfBoundary) const
{
  return FindBoundaryPixelsInRegion(region, whichSideOfBoundary).size();
}

unsigned int Mask::CountBoundaryPixels(const Mask::PixelType& whichSideOfBoundary) const
{
  return CountBoundaryPixels(this->GetLargestPossibleRegion(), whichSideOfBoundary);
}

/** Find hole pixels that are touching valid pixels.*/
std::vector<itk::Index<2> > Mask::FindBoundaryPixelsInRegion(const itk::ImageRegion<2>& region,
                                                             const Mask::PixelType& whichSideOfBoundary) const
{
  std::vector<itk::Index<2> > boundaryPixels;

  itk::ImageRegionConstIteratorWithIndex<Mask> maskIterator(this, region);

  while(!maskIterator.IsAtEnd())
  {
    if(maskIterator.Get() == whichSideOfBoundary)
    {
      if(ITKHelpers::HasNeighborWithValueOtherThan(maskIterator.GetIndex(), this,
                                                   whichSideOfBoundary))
      {
        boundaryPixels.push_back(maskIterator.GetIndex());
      }
    }

    ++maskIterator;
  }
  return boundaryPixels;
}

/** Find hole pixels that are touching valid pixels.*/
std::vector<itk::Index<2> > Mask::FindBoundaryPixels(const Mask::PixelType& whichSideOfBoundary) const
{
  return FindBoundaryPixelsInRegion(this->GetLargestPossibleRegion(), whichSideOfBoundary);
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
  // Ensure the region is inside the image
  region.Crop(this->GetLargestPossibleRegion());

  std::vector<itk::Index<2> > indices =
      ITKHelpers::GetPixelsWithValueInRegion(this, region, HoleMaskPixelTypeEnum::VALID);

  std::vector<itk::Offset<2> > validOffsets =
      ITKHelpers::IndicesToOffsets(indices, region.GetIndex());

  return validOffsets;
}

std::vector<itk::Offset<2> > Mask::GetHoleOffsetsInRegion(itk::ImageRegion<2> region) const
{
  // Ensure the region is inside the image
  region.Crop(this->GetLargestPossibleRegion());

  std::vector<itk::Index<2> > indices =
      ITKHelpers::GetPixelsWithValueInRegion(this, region, HoleMaskPixelTypeEnum::HOLE);

  std::vector<itk::Offset<2> > holeOffsets =
      ITKHelpers::IndicesToOffsets(indices, region.GetIndex());

  return holeOffsets;
}

std::vector<itk::Index<2> > Mask::GetValidPixels(const bool forward) const
{
  return GetValidPixelsInRegion(this->GetLargestPossibleRegion(), forward);
}

std::vector<itk::Index<2> > Mask::GetValidPixelsInRegion(itk::ImageRegion<2> region,
                                                         const bool forward) const
{
  // Ensure the region is inside the image
  region.Crop(this->GetLargestPossibleRegion());

  std::vector<itk::Index<2> > validPixels =
      ITKHelpers::GetPixelsWithValueInRegion(this, region, HoleMaskPixelTypeEnum::VALID);

  if(!forward)
  {
    std::reverse(validPixels.begin( ), validPixels.end( ) );
  }

  return validPixels;
}

std::vector<itk::Index<2> > Mask::GetHolePixelsInRegion(itk::ImageRegion<2> region) const
{
  // Ensure the region is inside the image
  region.Crop(this->GetLargestPossibleRegion());

  std::vector<itk::Index<2> > holePixels =
      ITKHelpers::GetPixelsWithValueInRegion(this, region, HoleMaskPixelTypeEnum::HOLE);

  return holePixels;
}

bool Mask::IsHole(const itk::Index<2>& index) const
{
  if(this->GetPixel(index) == HoleMaskPixelTypeEnum::HOLE)
  {
    return true;
  }
  return false;
}

bool Mask::IsHole(const itk::ImageRegion<2>& region) const
{
  return ITKHelpers::AllPixelsEqualTo(this, region,
                                      HoleMaskPixelTypeEnum::HOLE);
}

bool Mask::IsValid(const itk::ImageRegion<2>& region) const
{
  return ITKHelpers::AllPixelsEqualTo(this, region,
                                      HoleMaskPixelTypeEnum::VALID);
}

bool Mask::IsValid(const itk::Index<2>& index) const
{
  if(this->GetPixel(index) == HoleMaskPixelTypeEnum::VALID)
  {
    return true;
  }
  return false;
}

void Mask::InvertData()
{
  // Exchange HoleValue and ValidValue, but leave everything else alone.
  itk::ImageRegionIterator<Mask>
      maskIterator(this, this->GetLargestPossibleRegion());
  unsigned int invertedCounter = 0;
  while(!maskIterator.IsAtEnd())
  {
    if(this->IsValid(maskIterator.GetIndex()))
    {
      maskIterator.Set(HoleMaskPixelTypeEnum::HOLE);
      invertedCounter++;
    }
    else if(this->IsHole(maskIterator.GetIndex()))
    {
      maskIterator.Set(HoleMaskPixelTypeEnum::VALID);
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

  throw std::runtime_error("Cleanup() is not yet implemeneted!");
//  itk::ImageRegionIterator<Mask> maskIterator(this, this->GetLargestPossibleRegion());

//  float tolerance = 4;
//  while(!maskIterator.IsAtEnd())
//  {
//    if(abs(maskIterator.Get() - this->ValidValue) < tolerance)
//    {
//      //std::cout << "Setting valid pixel to " << static_cast<unsigned int>(this->ValidValue) << std::endl;
//      maskIterator.Set(this->ValidValue);
//    }
//    else if(abs(maskIterator.Get() - this->HoleValue) < tolerance)
//    {
//      //std::cout << "Setting hole pixel to " << static_cast<unsigned int>(this->HoleValue) << std::endl;
//      maskIterator.Set(this->HoleValue);
//    }
//    ++maskIterator;
//  }

}

void Mask::CopyHolesFrom(const Mask* const inputMask)
{
  itk::ImageRegionConstIterator<Mask> inputIterator(inputMask, inputMask->GetLargestPossibleRegion());
  itk::ImageRegionIterator<Mask> thisIterator(this, this->GetLargestPossibleRegion());

  while(!inputIterator.IsAtEnd())
  {
    if(inputMask->IsHole(inputIterator.GetIndex()))
    {
      thisIterator.Set(HoleMaskPixelTypeEnum::HOLE);
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
  this->CreateHolesFromValue(dilateFilter->GetOutput(), 255);
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
  this->CreateValidPixelsFromValue(erodeFilter->GetOutput(), 0);
}

void Mask::CreateImageInRegion(const itk::ImageRegion<2>& region, UnsignedCharImageType* const image,
                               const unsigned char holeColor,
                               const unsigned char validColor, const unsigned char undeterminedColor) const
{
    image->SetRegions(this->GetLargestPossibleRegion());
    image->Allocate();

    itk::ImageRegionIterator<UnsignedCharImageType>
        binaryImageIterator(image, region);

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
        binaryImageIterator.Set(undeterminedColor);
      }
      ++binaryImageIterator;
    }
}

void Mask::CreateImage(UnsignedCharImageType* const image, const unsigned char holeColor,
                       const unsigned char validColor, const unsigned char undeterminedColor) const
{
  CreateImageInRegion(image->GetLargestPossibleRegion(), image, holeColor,
                      validColor, undeterminedColor);
}

void Mask::CreateBinaryImageInRegion(const itk::ImageRegion<2>& region, UnsignedCharImageType* const image,
                                     const unsigned char holeColor,
                                     const unsigned char validColor) const
{
  CreateImageInRegion(region, image, holeColor, validColor, validColor);
}

void Mask::CreateBinaryImage(UnsignedCharImageType* const image, const unsigned char holeColor,
                             const unsigned char validColor) const
{
  CreateImage(image, holeColor, validColor, validColor);
}

void Mask::CreateBoundaryImageInRegion(const itk::ImageRegion<2>& region, BoundaryImageType* const boundaryImage,
                                       const HoleMaskPixelTypeEnum& whichSideOfBoundary) const
{
  // Create a binary image of the mask
  unsigned char holeColor = 255;
  unsigned char validColor = 0;
  UnsignedCharImageType::Pointer fullBinaryImage = UnsignedCharImageType::New();
  CreateBinaryImageInRegion(region, fullBinaryImage, holeColor, validColor);

  // Extract the relevant region from the binary image
  UnsignedCharImageType::Pointer binaryImage = UnsignedCharImageType::New();
  binaryImage->SetRegions(region);
  binaryImage->Allocate();
  ITKHelpers::CopyRegion(fullBinaryImage.GetPointer(), binaryImage.GetPointer(), region, region);

  // Extract the relevant region from the mask
  Mask::Pointer extractedRegionMask = Mask::New();
  extractedRegionMask->SetRegions(region);
  extractedRegionMask->Allocate();
  ITKHelpers::CopyRegion(this, extractedRegionMask.GetPointer(), region, region);

  // Since the hole is white (we have specified this at the beginning of this function),
  // we want the foreground value of the contour filter to be black.
  // This means that the boundary will be detected in the black pixel region,
  // which is on the outside edge of the hole like we want. However,
  // The BinaryContourImageFilter will change all non-boundary pixels to the background color,
  // so the resulting output will be inverted - the boundary pixels will be black and the
  // non-boundary pixels will be white.

  // Find the boundary
  typedef itk::BinaryContourImageFilter<UnsignedCharImageType, UnsignedCharImageType> binaryContourImageFilterType;
  binaryContourImageFilterType::Pointer binaryContourFilter = binaryContourImageFilterType::New();
  binaryContourFilter->SetInput(binaryImage);
  binaryContourFilter->SetFullyConnected(true);

  if(whichSideOfBoundary == HoleMaskPixelTypeEnum::VALID)
  {
    // we want the boundary pixels to be in the valid region.
      binaryContourFilter->SetForegroundValue(validColor);
      binaryContourFilter->SetBackgroundValue(holeColor);
  }
  else if(whichSideOfBoundary == HoleMaskPixelTypeEnum::HOLE)
  {
    // we want the boundary pixels to be in the hole region.
      binaryContourFilter->SetForegroundValue(holeColor);
      binaryContourFilter->SetBackgroundValue(validColor);
  }
  else
  {
    throw std::runtime_error("An invalid side of the boundary was requested.");
  }

  binaryContourFilter->Update();

  ITKHelpers::DeepCopy(binaryContourFilter->GetOutput(), boundaryImage);
}

void Mask::CreateBoundaryImage(itk::Image<unsigned char, 2>* const boundaryImage,
                               const HoleMaskPixelTypeEnum& whichSideOfBoundary) const
{
  CreateBoundaryImageInRegion(this->GetLargestPossibleRegion(), boundaryImage,
                              whichSideOfBoundary);
}

/** Get a list of the valid neighbors of a pixel.*/
std::vector<itk::Index<2> > Mask::GetValid8Neighbors(const itk::Index<2>& pixel) const
{
  return ITKHelpers::Get8NeighborsWithValue(pixel, this, HoleMaskPixelTypeEnum::VALID);
}

std::vector<itk::Index<2> > Mask::GetValid8NeighborsInRegion(const itk::Index<2>& pixel, const itk::ImageRegion<2>& region) const
{
  return ITKHelpers::Get8NeighborsInRegionWithValue(pixel, this, region,
                                                    HoleMaskPixelTypeEnum::VALID);
}

bool Mask::HasHole8NeighborInRegion(const itk::Index<2>& pixel,
                                   const itk::ImageRegion<2>& region) const
{
  // Return true of there are more than 0 neighbors
  return GetHole8NeighborsInRegion(pixel, region).size() > 0;
}

bool Mask::HasHole8Neighbor(const itk::Index<2>& pixel) const
{
  // Return true of there are more than 0 neighbors
  return GetHole8Neighbors(pixel).size() > 0;
}

bool Mask::HasValid8Neighbor(const itk::Index<2>& pixel) const
{
  // Return true of there are more than 0 neighbors
  return GetValid8Neighbors(pixel).size() > 0;
}

/** Get a list of the hole neighbors of a pixel.*/
std::vector<itk::Index<2> > Mask::GetHole8Neighbors(const itk::Index<2>& pixel) const
{
  return ITKHelpers::Get8NeighborsWithValue(pixel, this, HoleMaskPixelTypeEnum::HOLE);
}

std::vector<itk::Index<2> > Mask::GetHole8NeighborsInRegion(const itk::Index<2>& pixel,
                                                            const itk::ImageRegion<2>& region) const
{
  return ITKHelpers::Get8NeighborsInRegionWithValue(pixel, this, region,
                                                    HoleMaskPixelTypeEnum::HOLE);
}

std::vector<itk::Offset<2> > Mask::GetValid8NeighborOffsets(const itk::Index<2>& pixel) const
{
  std::vector<itk::Index<2> > indices =
      ITKHelpers::Get8NeighborsWithValue(pixel, this, HoleMaskPixelTypeEnum::VALID);

  std::vector<itk::Offset<2> > offsets = ITKHelpers::IndicesToOffsets(indices, pixel);
  return offsets;
}

std::vector<itk::Offset<2> > Mask::GetHole8NeighborOffsets(const itk::Index<2>& pixel) const
{
  std::vector<itk::Index<2> > indices =
      ITKHelpers::Get8NeighborsWithValue(pixel, this, HoleMaskPixelTypeEnum::HOLE);

  std::vector<itk::Offset<2> > offsets = ITKHelpers::IndicesToOffsets(indices, pixel);

  return offsets;
}

void Mask::MarkAsHole(const itk::Index<2>& pixel)
{
  this->SetPixel(pixel, HoleMaskPixelTypeEnum::HOLE);
}

void Mask::MarkAsValid(const itk::Index<2>& pixel)
{
  this->SetPixel(pixel, HoleMaskPixelTypeEnum::VALID);
}

bool Mask::HasValid4Neighbor(const itk::Index<2>& pixel)
{
  return ITKHelpers::Has4NeighborsWithValue(this, pixel, HoleMaskPixelTypeEnum::VALID);
}

std::vector<itk::Index<2> > Mask::GetValid4Neighbors(const itk::Index<2>& pixel)
{
  return ITKHelpers::Get4NeighborsWithValue(this, pixel, HoleMaskPixelTypeEnum::VALID);
}

void Mask::KeepLargestHole()
{
  throw std::runtime_error("Mask::KeepLargestHole() needs to be tested!");

  typedef itk::Image<unsigned int> LabelImageType;
  // Only keep the largest segment
  typedef itk::ConnectedComponentImageFilter<Mask, LabelImageType>
      ConnectedComponentImageFilterType;
  ConnectedComponentImageFilterType::Pointer connectedComponentFilter =
      ConnectedComponentImageFilterType::New();
  connectedComponentFilter->SetInput(this);
  connectedComponentFilter->Update();

  ITKHelpers::WriteImage(connectedComponentFilter->GetOutput(), "ConnectedComponents.mha");

  typedef itk::LabelShapeKeepNObjectsImageFilter<LabelImageType>
      LabelShapeKeepNObjectsImageFilterType;
  LabelShapeKeepNObjectsImageFilterType::Pointer labelShapeKeepNObjectsImageFilter =
           LabelShapeKeepNObjectsImageFilterType::New();
  labelShapeKeepNObjectsImageFilter->SetInput(connectedComponentFilter->GetOutput());
  labelShapeKeepNObjectsImageFilter->SetBackgroundValue(0); // TODO: why zero?
  labelShapeKeepNObjectsImageFilter->SetNumberOfObjects(1);
  labelShapeKeepNObjectsImageFilter
            ->SetAttribute(LabelShapeKeepNObjectsImageFilterType::LabelObjectType::NUMBER_OF_PIXELS);
  labelShapeKeepNObjectsImageFilter->Update();

  ITKHelpers::WriteImage(labelShapeKeepNObjectsImageFilter->GetOutput(),
                         "LargestComponent.mha");

  // TODO: replace whatever value is the foreground (by looking at the output above)
  // with HOLE
//  ITKHelpers::DeepCopy(rescaleFilter->GetOutput(), this);
}

unsigned int Mask::CountValidPatches(const unsigned int patchRadius) const
{

  itk::ImageRegionConstIteratorWithIndex<Mask> maskIterator(this, this->GetLargestPossibleRegion());

  unsigned int counter = 0;
  // std::cout << "CountValidPatches (patch radius " << patchRadius << ")..." << std::endl;
  while(!maskIterator.IsAtEnd())
  {
    itk::ImageRegion<2> region =
        ITKHelpers::GetRegionInRadiusAroundPixel(maskIterator.GetIndex(), patchRadius);

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
  itk::ImageRegionConstIteratorWithIndex<Mask>
      maskIterator(this, this->GetLargestPossibleRegion());

  while(!maskIterator.IsAtEnd())
  {
    itk::ImageRegion<2> region =
        ITKHelpers::GetRegionInRadiusAroundPixel(maskIterator.GetIndex(), patchRadius);

    if(this->IsValid(region))
    {
      return region;
    }

    ++maskIterator;
  }

  throw std::runtime_error("No valid patches found!");
}

void Mask::SetHole(const itk::Index<2>& index)
{
  this->SetPixel(index, HoleMaskPixelTypeEnum::HOLE);
}

void Mask::SetValid(const itk::Index<2>& index)
{
  this->SetPixel(index, HoleMaskPixelTypeEnum::VALID);
}

void Mask::SetValid(const itk::ImageRegion<2>& region)
{
  ITKHelpers::SetRegionToConstant(this, region, HoleMaskPixelTypeEnum::VALID);
}

std::ostream& operator<<(std::ostream& output, const HoleMaskPixelTypeEnum &pixelType)
{
  if(pixelType == HoleMaskPixelTypeEnum::HOLE)
  {
    output << "HoleMaskPixelTypeEnum::HOLE" << std::endl;
  }
  else if(pixelType == HoleMaskPixelTypeEnum::VALID)
  {
    output << "HoleMaskPixelTypeEnum::VALID" << std::endl;
  }
  else
  {
    output << "HoleMaskPixelTypeEnum::UNDETERMINED" << std::endl;
  }

  return output;
}
