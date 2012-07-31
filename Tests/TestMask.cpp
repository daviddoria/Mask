#include "Mask.h"

// Submodules
#include <ITKHelpers/ITKHelpers.h>

void TestFindBoundaryInRegion();

int main( int argc, char ** argv )
{
  TestFindBoundaryInRegion();
  return 0;
}

void TestFindBoundaryInRegion()
{
  // Test with 255=valid
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

    ITKHelpers::WriteImage(mask.GetPointer(), "mask_A.png");

    itk::Index<2> regionCorner = {{40,40}};
    itk::Size<2> regionSize = {{20,20}};
    itk::ImageRegion<2> queryRegion(regionCorner, regionSize);

    {
    Mask::BoundaryImageType::Pointer boundaryImage = Mask::BoundaryImageType::New();
    mask->FindBoundaryInRegion(queryRegion, boundaryImage, Mask::VALID, 255);
    ITKHelpers::WriteImage(boundaryImage.GetPointer(), "boundary_A_valid.png");
    }

    {
    Mask::BoundaryImageType::Pointer boundaryImage = Mask::BoundaryImageType::New();
    mask->FindBoundaryInRegion(queryRegion, boundaryImage, Mask::HOLE, 255);
    ITKHelpers::WriteImage(boundaryImage.GetPointer(), "boundary_A_hole.png");
    }
  }

  // Test with 255=hole
  {
    Mask::Pointer mask = Mask::New();
    itk::Index<2> corner = {{0,0}};
    itk::Size<2> size = {{100,100}};
    itk::ImageRegion<2> imageRegion(corner, size);
    mask->SetRegions(imageRegion);
    mask->Allocate();

    mask->SetValidValue(0);
    mask->SetHoleValue(255);
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

    ITKHelpers::WriteImage(mask.GetPointer(), "mask_B.png");

    itk::Index<2> regionCorner = {{40,40}};
    itk::Size<2> regionSize = {{20,20}};
    itk::ImageRegion<2> queryRegion(regionCorner, regionSize);

    {
    Mask::BoundaryImageType::Pointer boundaryImage = Mask::BoundaryImageType::New();
    mask->FindBoundaryInRegion(queryRegion, boundaryImage, Mask::VALID, 255);
    ITKHelpers::WriteImage(boundaryImage.GetPointer(), "boundary_B_valid.png");
    }

    {
    Mask::BoundaryImageType::Pointer boundaryImage = Mask::BoundaryImageType::New();
    mask->FindBoundaryInRegion(queryRegion, boundaryImage, Mask::HOLE, 255);
    ITKHelpers::WriteImage(boundaryImage.GetPointer(), "boundary_B_hole.png");
    }
  }
}
