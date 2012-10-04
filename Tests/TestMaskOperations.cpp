#include "Mask.h"
#include "MaskOperations.h"

// Submodules
#include <ITKHelpers/ITKHelpers.h>

static void TestComputeHoleBoundingBox();
static void TestInterpolateHole();
static void TestMaskedBlur();
static void TestFindMinimumValueInMaskedRegion();
static void TestFindMaximumValueInMaskedRegion();

template <typename TImage>
static void CreateImage(TImage* const image);

static void CreateMask(Mask* const mask);

int main(int, char**)
{
  TestInterpolateHole();
  TestComputeHoleBoundingBox();

  TestMaskedBlur();

  TestFindMinimumValueInMaskedRegion();
  TestFindMaximumValueInMaskedRegion();

  return 0;
}

void TestInterpolateHole()
{
  itk::Index<2> corner = {{0,0}};
  itk::Size<2> size = {{100,100}};
  itk::ImageRegion<2> imageRegion(corner, size);
  
  Mask::Pointer mask = Mask::New();
  mask->SetRegions(imageRegion);
  mask->Allocate();

  typedef itk::Image<float, 2> ImageType;
  ImageType::Pointer image = ImageType::New();
  image->SetRegions(imageRegion);
  image->Allocate();
  
  MaskOperations::InterpolateHole(image.GetPointer(), mask);
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


template <typename TImage>
void CreateImage(TImage* const image)
{
  typename TImage::IndexType corner;
  corner.Fill(0);

  typename TImage::SizeType size;
  size.Fill(100);

  typename TImage::RegionType region(corner, size);

  image->SetRegions(region);
  image->Allocate();

  itk::ImageRegionIterator<TImage> imageIterator(image,region);

  while(!imageIterator.IsAtEnd())
  {
//    if(imageIterator.GetIndex()[0] < 70)
//    {
//      imageIterator.Set(255);
//    }
//    else
//    {
//      imageIterator.Set(0);
//    }

    imageIterator.Set(rand() % 255);
    ++imageIterator;
  }

}

void CreateMask(Mask* const mask)
{
  typename Mask::IndexType corner;
  corner.Fill(0);

  typename Mask::SizeType size;
  size.Fill(100);

  typename Mask::RegionType region(corner, size);

  mask->SetRegions(region);
  mask->Allocate();

  itk::ImageRegionIterator<Mask> maskIterator(mask, region);

  while(!maskIterator.IsAtEnd())
  {
    if(maskIterator.GetIndex()[0] < 70)
    {
      maskIterator.Set(mask->GetValidValue());
    }
    else
    {
      maskIterator.Set(mask->GetHoleValue());
    }

    ++maskIterator;
  }

}

void TestMaskedBlur()
{
  // Scalar
  {
  typedef itk::Image<unsigned char, 2> ImageType;
  ImageType::Pointer image = ImageType::New();
  CreateImage(image.GetPointer());
  ITKHelpers::WriteImage(image.GetPointer(), "ScalarImage.png");

  Mask::Pointer mask = Mask::New();
  CreateMask(mask);
  std::cout << "Mask hole value: " << static_cast<int>(mask->GetHoleValue()) << std::endl;
  ITKHelpers::WriteImage(mask.GetPointer(), "Mask.png");

  // Test with a normal variance
  float blurVariance = 2.0f;
  ImageType::Pointer output_2 = ImageType::New();
  MaskOperations::MaskedBlur(image.GetPointer(), mask, blurVariance, output_2.GetPointer());
  ITKHelpers::WriteImage(output_2.GetPointer(), "ScalarBlurred_2.png");

  // Test with zero variance (i.e. don't blur the image)
  blurVariance = 0.0f;
  ImageType::Pointer output_0 = ImageType::New();
  MaskOperations::MaskedBlur(image.GetPointer(), mask, blurVariance, output_0.GetPointer());
  ITKHelpers::WriteImage(output_0.GetPointer(), "ScalarBlurred_0.png");
  }

  // Vector
  {
  typedef itk::Image<itk::CovariantVector<unsigned char, 3>, 2> ImageType;
  ImageType::Pointer image = ImageType::New();
  CreateImage(image.GetPointer());
  ITKHelpers::WriteImage(image.GetPointer(), "VectorImage.png");

  Mask::Pointer mask = Mask::New();
  CreateMask(mask);
  std::cout << "Mask hole value: " << static_cast<int>(mask->GetHoleValue()) << std::endl;
  ITKHelpers::WriteImage(mask.GetPointer(), "Mask.png");

  float blurVariance = 2.0f;
  ImageType::Pointer output = ImageType::New();
  MaskOperations::MaskedBlur(image.GetPointer(), mask, blurVariance, output.GetPointer());

  ITKHelpers::WriteImage(output.GetPointer(), "VectorBlurred.png");
  }
}

void TestFindMaximumValueInMaskedRegion()
{
  // Scalar
  {
  typedef itk::Image<int, 2> ImageType;
  ImageType::Pointer image = ImageType::New();
  CreateImage(image.GetPointer());
//  ITKHelpers::WriteImage(image.GetPointer(), "ScalarImage.png");

  Mask::Pointer mask = Mask::New();
  CreateMask(mask);
  std::cout << "Mask hole value: " << static_cast<int>(mask->GetHoleValue()) << std::endl;
//  ITKHelpers::WriteImage(mask.GetPointer(), "Mask.png");

  itk::Index<2> regionCorner = {{0,0}};
  itk::Size<2> regionSize = {{10,10}};
  itk::ImageRegion<2> region(regionCorner, regionSize);

  ImageType::PixelType maxValue;
  MaskOperations::FindMaximumValueInMaskedRegion(image.GetPointer(), mask,
                                                 region, mask->GetValidValue(),
                                                 maxValue);

  std::cout << "Scalar max: " << maxValue << std::endl;
  }

  // Vector
  {
  typedef itk::Image<itk::CovariantVector<unsigned char, 3>, 2> ImageType;
  ImageType::Pointer image = ImageType::New();
  CreateImage(image.GetPointer());
//  ITKHelpers::WriteImage(image.GetPointer(), "VectorImage.png");

  Mask::Pointer mask = Mask::New();
  CreateMask(mask);
  std::cout << "Mask hole value: " << static_cast<int>(mask->GetHoleValue()) << std::endl;
//  ITKHelpers::WriteImage(mask.GetPointer(), "Mask.png");

  itk::Index<2> regionCorner = {{0,0}};
  itk::Size<2> regionSize = {{10,10}};
  itk::ImageRegion<2> region(regionCorner, regionSize);

  ImageType::PixelType maxValue;
  MaskOperations::FindMaximumValueInMaskedRegion(image.GetPointer(), mask,
                                                 region, mask->GetValidValue(),
                                                 maxValue);

  std::cout << "Vector max: " << maxValue << std::endl;
  }

}

void TestFindMinimumValueInMaskedRegion()
{
  // Scalar
  {
  typedef itk::Image<int, 2> ImageType;
  ImageType::Pointer image = ImageType::New();
  CreateImage(image.GetPointer());
//  ITKHelpers::WriteImage(image.GetPointer(), "ScalarImage.png");

  Mask::Pointer mask = Mask::New();
  CreateMask(mask);
  std::cout << "Mask hole value: " << static_cast<int>(mask->GetHoleValue()) << std::endl;
//  ITKHelpers::WriteImage(mask.GetPointer(), "Mask.png");

  itk::Index<2> regionCorner = {{0,0}};
  itk::Size<2> regionSize = {{10,10}};
  itk::ImageRegion<2> region(regionCorner, regionSize);

  ImageType::PixelType minValue;
  MaskOperations::FindMinimumValueInMaskedRegion(image.GetPointer(), mask,
                                                 region, mask->GetValidValue(),
                                                 minValue);

  std::cout << "Scalar min: " << minValue << std::endl;
  }

  // Vector
  {
  typedef itk::Image<itk::CovariantVector<unsigned char, 3>, 2> ImageType;
  ImageType::Pointer image = ImageType::New();
  CreateImage(image.GetPointer());
//  ITKHelpers::WriteImage(image.GetPointer(), "VectorImage.png");

  Mask::Pointer mask = Mask::New();
  CreateMask(mask);
  std::cout << "Mask hole value: " << static_cast<int>(mask->GetHoleValue()) << std::endl;
//  ITKHelpers::WriteImage(mask.GetPointer(), "Mask.png");

  itk::Index<2> regionCorner = {{0,0}};
  itk::Size<2> regionSize = {{10,10}};
  itk::ImageRegion<2> region(regionCorner, regionSize);

  ImageType::PixelType minValue;
  MaskOperations::FindMinimumValueInMaskedRegion(image.GetPointer(), mask,
                                                 region, mask->GetValidValue(),
                                                 minValue);

  std::cout << "Vector min: " << minValue << std::endl;
  }
}

