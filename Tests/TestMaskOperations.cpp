#include "Mask.h"
#include "MaskOperations.h"

#include "ITKHelpers/ITKHelpers.h"

void TestComputeHoleBoundingBox();

int main( int argc, char ** argv )
{
  TestComputeHoleBoundingBox();
  return 0;
}

void TestComputeHoleBoundingBox()
{
  Mask::Pointer mask = Mask::New();
  itk::Index<2> corner = {{0,0}};
  itk::Size<2> size = {{100,100}};
  itk::ImageRegion<2> imageRegion(corner, size);
  mask->SetRegions(imageRegion);
  mask->Allocate();

  mask->SetValidValue(255);
  mask->SetHoleValue(0);
  mask->FillBuffer(mask->GetValidValue());

  itk::ImageRegionIterator<Mask> maskIterator(mask, mask->GetLargestPossibleRegion());

  while(!maskIterator.IsAtEnd())
    {
    if(maskIterator.GetIndex()[0] > 50 && maskIterator.GetIndex()[0] < 70 &&
      maskIterator.GetIndex()[1] > 50 && maskIterator.GetIndex()[1] < 70)
      {
      maskIterator.Set(mask->GetHoleValue());
      }

    ++maskIterator;
    }

  itk::ImageRegion<2> boundingBox = MaskOperations::ComputeHoleBoundingBox(mask);

  std::cout << "Bounding box: " << boundingBox << std::endl;
}
