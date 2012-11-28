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

#include "ForegroundBackgroundSegmentMask.h"

// Submodules
#include <Helpers/Helpers.h>
#include <ITKHelpers/ITKHelpers.h>

void ForegroundBackgroundSegmentMask::Read(const std::string& filename)
{
  /**
   * The format of the .fbmask (foreground/background mask) file is:
   * foreground 0
   * background 255
   * Mask.png
   *
   * OR
   *
   * background 255
   * foreground 0
   * Mask.png
   *
   * That is, the "foreground [VALUE]" line can be either on the first or second line.
   * Note that the 0 and 255 here are arbitrary and can be anything.
   */
  std::string extension = Helpers::GetFileExtension(filename);
  if(extension != "mask")
  {
    std::stringstream ss;
    ss << "Cannot read any file except .fbmask! Specified file was ." << extension
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
  int foregroundValue;
  int backgroundValue;

  std::string type1;
  std::string type2;
  int value1;
  int value2;

  getline(fin, line);
  linestream.clear();
  linestream << line;
  linestream >> type1 >> value1;

  if(type1 == "foreground")
  {
    foregroundValue = value1;
  }
  else if(type1 == "background")
  {
    backgroundValue = value1;
  }
  else
  {
    throw std::runtime_error("Invalid .fbmask file!");
  }

  getline(fin, line);
  linestream.clear();
  linestream << line;
  linestream >> type2 >> value2;

  if(type1 == type2)
  {
    throw std::runtime_error("Invalid .fbmask file! Foreground or background value listed twice!");
  }

  if(type2 == "foreground")
  {
    foregroundValue = value2;
  }
  else if(type2 == "background")
  {
    backgroundValue = value2;
  }
  else
  {
    throw std::runtime_error("Invalid .fbmask file!");
  }

  std::string imageFileName;
  linestream >> imageFileName;

  std::string path = Helpers::GetPath(filename);

  std::string fullImageFileName = path + imageFileName;

  ReadFromImage(fullImageFileName, ForegroundPixelValueWrapper<int>(foregroundValue),
                BackgroundPixelValueWrapper<int>(backgroundValue));
}


bool ForegroundBackgroundSegmentMask::IsForeground(const itk::Index<2>& index) const
{
  if(this->GetPixel(index) == ForegroundBackgroundSegmentMaskPixelTypeEnum::FOREGROUND)
  {
    return true;
  }
  return false;
}

bool ForegroundBackgroundSegmentMask::IsBackground(const itk::Index<2>& index) const
{
  if(this->GetPixel(index) == ForegroundBackgroundSegmentMaskPixelTypeEnum::BACKGROUND)
  {
    return true;
  }
  return false;
}


std::ostream& operator<<(std::ostream& output,
                         const ForegroundBackgroundSegmentMaskPixelTypeEnum &pixelType)
{
  if(pixelType == ForegroundBackgroundSegmentMaskPixelTypeEnum::FOREGROUND)
  {
    output << "HoleMaskPixelTypeEnum::FOREGROUND" << std::endl;
  }
  else if(pixelType == ForegroundBackgroundSegmentMaskPixelTypeEnum::BACKGROUND)
  {
    output << "HoleMaskPixelTypeEnum::BACKGROUND" << std::endl;
  }

  return output;
}
