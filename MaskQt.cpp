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

  itk::ImageRegionIterator<Mask> imageIterator(regionOfInterestImageFilter->GetOutput(),
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

  //return qimage; // The actual image region
  return qimage.mirrored(false, true); // The flipped image region
}

} // end namespace