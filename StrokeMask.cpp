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

#include "StrokeMask.h"

// Submodules
#include <Helpers/Helpers.h>
#include <ITKHelpers/ITKHelpers.h>

void StrokeMask::Read(const std::string& filename)
{
  /**
   * The format of the .stroke file is:
   * stroke 255
   * Mask.png
   */
  std::string extension = Helpers::GetFileExtension(filename);
  if(extension != "stroke")
  {
    std::stringstream ss;
    ss << "StrokeMask cannot read files with extension other than .stroke! Specified file had extension ." << extension
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
  std::string strokeString;
  int strokeValue;

  getline(fin, line);
  linestream.clear();
  linestream << line;
  linestream >> strokeString >> strokeValue;

  if(strokeString != "stroke")
  {
    throw std::runtime_error("Invalid .stroke file!");
  }

  std::cout << "strokeValue: " << strokeValue << std::endl;

  std::string imageFileName;
  getline(fin, imageFileName);

  if(imageFileName.length() == 0)
  {
    throw std::runtime_error("Image file name was empty!");
  }

  std::string path = Helpers::GetPath(filename);

  std::string fullImageFileName = path + imageFileName;

  ReadFromImage(fullImageFileName, strokeValue);
}


bool StrokeMask::IsStroke(const itk::Index<2>& index) const
{
  if(this->GetPixel(index) == StrokeMaskPixelTypeEnum::STROKE)
  {
    return true;
  }
  return false;
}

unsigned int StrokeMask::CountStrokePixels() const
{
  std::vector<itk::Index<2> > foregroundPixels =
      ITKHelpers::GetPixelsWithValueInRegion(this, this->GetLargestPossibleRegion(),
                                             StrokeMaskPixelTypeEnum::STROKE);

  return foregroundPixels.size();
}

std::ostream& operator<<(std::ostream& output,
                         const StrokeMaskPixelTypeEnum &pixelType)
{
  if(pixelType == StrokeMaskPixelTypeEnum::STROKE)
  {
    output << "StrokeMaskPixelTypeEnum::STROKE" << std::endl;
  }
  else if(pixelType == StrokeMaskPixelTypeEnum::NOTSTROKE)
  {
    output << "StrokeMaskPixelTypeEnum::NOTSTROKE" << std::endl;
  }

  return output;
}
