#include "ForegroundBackgroundSegmentMask.h"

// Submodules
#include <ITKHelpers/ITKHelpers.h>

int main(int argc, char*argv[])
{
  if(argc < 2)
  {
    std::cerr << "Usage: file.mask" << std::endl;
    return EXIT_FAILURE;
  }

  std::string fileName = argv[1];

  ForegroundBackgroundSegmentMask::Pointer mask = ForegroundBackgroundSegmentMask::New();
  mask->Read(fileName);
  unsigned int numberOfForegroundPixels = mask->CountForegroundPixels();
  unsigned int numberOfBackgroundPixels = mask->CountBackgroundPixels();

  std::cout << "numberOfForegroundPixels: " << numberOfForegroundPixels
            << " numberOfBackgroundPixels: " << numberOfBackgroundPixels << std::endl;

  // There should be 100 foreground pixels (10 x 10)
  // There should be (100*100) - (10*10), a 100x100 image with a 10x10 foreground region
  if(numberOfForegroundPixels == 100 &&
     numberOfBackgroundPixels == 100*100 - 100)
  {
    return EXIT_SUCCESS;
  }
  else
  {
    return EXIT_FAILURE;
  }
}
