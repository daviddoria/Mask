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

#include "MaskQt.h"

// ITK
#include "itkRegionOfInterestImageFilter.h"

// Qt
#include <QColor>

namespace MaskQt
{
  
QImage GetQtImage(const Mask* const mask)
{
  return GetQtImage(mask, mask->GetLargestPossibleRegion());
}

QImage GetQtImage(const Mask* const mask, const itk::ImageRegion<2>& region)
{
  QImage qimage(region.GetSize()[0], region.GetSize()[1], QImage::Format_ARGB32);

  typedef itk::RegionOfInterestImageFilter<Mask, Mask> RegionOfInterestImageFilterType;
  typename RegionOfInterestImageFilterType::Pointer regionOfInterestImageFilter =
           RegionOfInterestImageFilterType::New();
  regionOfInterestImageFilter->SetRegionOfInterest(region);
  regionOfInterestImageFilter->SetInput(mask);
  regionOfInterestImageFilter->Update();

  itk::ImageRegionConstIterator<Mask> imageIterator(regionOfInterestImageFilter->GetOutput(),
                                                    regionOfInterestImageFilter->GetOutput()->GetLargestPossibleRegion());

  while(!imageIterator.IsAtEnd())
  {
    typename Mask::PixelType pixel = imageIterator.Get();

    itk::Index<2> index = imageIterator.GetIndex();
    // In a binary mask, R, G, and B are set to the same value.

    int r = static_cast<int>(pixel);
    int g = static_cast<int>(pixel);
    int b = static_cast<int>(pixel);

    unsigned int alpha = 0;

    if(mask->IsHole(index))
    {
      alpha = 255; // opaque
    }

    QColor pixelColor(r,g,b,alpha);

    qimage.setPixel(index[0], index[1], pixelColor.rgba());

    ++imageIterator;
  }

  return qimage; // The actual image region
}

QImage SetPixelsToTransparent(QImage image, const Mask* const mask,
                              HoleMaskPixelTypeEnum pixelValue)
{
  itk::ImageRegionConstIteratorWithIndex<Mask>
      maskIterator(mask, mask->GetLargestPossibleRegion());

  while(!maskIterator.IsAtEnd())
  {
    itk::Index<2> index = maskIterator.GetIndex();

    int alpha = 255; // opaque
    if(maskIterator.Get() == pixelValue)
    {
      alpha = 0; // transparent
    }

    QColor pixelColor = image.pixel(index[0], index[1]);
    pixelColor.setAlpha(alpha);

    image.setPixel(index[0], index[1], pixelColor.rgba());

    ++maskIterator;
  }

  return image;
}

} // end namespace
