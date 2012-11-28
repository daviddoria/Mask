#include "Mask.h"

// Submodules
#include <ITKHelpers/ITKHelpers.h>

static bool TestFindBoundaryInRegion();

static bool TestRead();

int main()
{
  bool allPass = true;
  allPass &= TestFindBoundaryInRegion();

  allPass &= TestRead();

  if(allPass)
  {
    return EXIT_SUCCESS;
  }
  else
  {
    return EXIT_FAILURE;
  }
}

bool TestFindBoundaryInRegion()
{
//  // Test with 255=valid
//  {
//    Mask::Pointer mask = Mask::New();
//    itk::Index<2> corner = {{0,0}};
//    itk::Size<2> size = {{100,100}};
//    itk::ImageRegion<2> imageRegion(corner, size);
//    mask->SetRegions(imageRegion);
//    mask->Allocate();

//    mask->FillBuffer(HoleMaskPixelTypeEnum::VALID);

//    itk::ImageRegionIterator<Mask> maskIterator(mask, mask->GetLargestPossibleRegion());

//    while(!maskIterator.IsAtEnd())
//    {
//      if(maskIterator.GetIndex()[0] > 50 && maskIterator.GetIndex()[0] < 70 &&
//         maskIterator.GetIndex()[1] > 50 && maskIterator.GetIndex()[1] < 70)
//      {
//        maskIterator.Set(HoleMaskPixelTypeEnum::HOLE);
//      }

//      ++maskIterator;
//    }

//    ITKHelpers::WriteImage(mask.GetPointer(), "mask_A.png");

//    itk::Index<2> regionCorner = {{40,40}};
//    itk::Size<2> regionSize = {{20,20}};
//    itk::ImageRegion<2> queryRegion(regionCorner, regionSize);

//    {
//    Mask::BoundaryImageType::Pointer boundaryImage = Mask::BoundaryImageType::New();
//    mask->CreateBoundaryImageInRegion(queryRegion, boundaryImage, HoleMaskPixelTypeEnum::VALID, 255);
//    ITKHelpers::WriteImage(boundaryImage.GetPointer(), "boundary_A_valid.png");
//    }

//    {
//    Mask::BoundaryImageType::Pointer boundaryImage = Mask::BoundaryImageType::New();
//    mask->CreateBoundaryImageInRegion(queryRegion, boundaryImage, HoleMaskPixelTypeEnum::HOLE, 255);
//    ITKHelpers::WriteImage(boundaryImage.GetPointer(), "boundary_A_hole.png");
//    }
//  }

//  // Test with 255=hole
//  {
//    Mask::Pointer mask = Mask::New();
//    itk::Index<2> corner = {{0,0}};
//    itk::Size<2> size = {{100,100}};
//    itk::ImageRegion<2> imageRegion(corner, size);
//    mask->SetRegions(imageRegion);
//    mask->Allocate();

//    mask->FillBuffer(HoleMaskPixelTypeEnum::VALID);

//    itk::ImageRegionIterator<Mask> maskIterator(mask, mask->GetLargestPossibleRegion());

//    while(!maskIterator.IsAtEnd())
//    {
//      if(maskIterator.GetIndex()[0] > 50 && maskIterator.GetIndex()[0] < 70 &&
//         maskIterator.GetIndex()[1] > 50 && maskIterator.GetIndex()[1] < 70)
//      {
//        maskIterator.Set(HoleMaskPixelTypeEnum::HOLE);
//      }

//      ++maskIterator;
//    }

//    ITKHelpers::WriteImage(mask.GetPointer(), "mask_B.png");

//    itk::Index<2> regionCorner = {{40,40}};
//    itk::Size<2> regionSize = {{20,20}};
//    itk::ImageRegion<2> queryRegion(regionCorner, regionSize);

//    {
//    Mask::BoundaryImageType::Pointer boundaryImage = Mask::BoundaryImageType::New();
//    mask->CreateBoundaryImageInRegion(queryRegion, boundaryImage, HoleMaskPixelTypeEnum::VALID, 255);
//    ITKHelpers::WriteImage(boundaryImage.GetPointer(), "boundary_B_valid.png");
//    }

//    {
//    Mask::BoundaryImageType::Pointer boundaryImage = Mask::BoundaryImageType::New();
//    mask->CreateBoundaryImageInRegion(queryRegion, boundaryImage, HoleMaskPixelTypeEnum::HOLE, 255);
//    ITKHelpers::WriteImage(boundaryImage.GetPointer(), "boundary_B_hole.png");
//    }
//  }

//  return true;

  return false;
}

bool TestRead()
{
  Mask::Pointer mask = Mask::New();
  mask->Read("/media/portable/Projects/ImageGraphCutSegmentation/Mask/Tests/data/TestMask.mask");
  unsigned int numberOfHoles = ITKHelpers::CountPixelsWithValue(mask.GetPointer(), HoleMaskPixelTypeEnum::HOLE);
  unsigned int numberOfValid = ITKHelpers::CountPixelsWithValue(mask.GetPointer(), HoleMaskPixelTypeEnum::VALID);

  std::cout << "numberOfHoles: " << numberOfHoles << " numberOfValid: " << numberOfValid << std::endl;

  return true;
}
