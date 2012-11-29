#include "Mask.h"

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

  Mask::Pointer mask = Mask::New();
  mask->Read(fileName);
  unsigned int numberOfHolePixels = mask->CountHolePixels();
  unsigned int numberOfValidPixels = mask->CountValidPixels();

  std::cout << "numberOfHolePixels: " << numberOfHolePixels
            << " numberOfValidPixels: " << numberOfValidPixels << std::endl;

  // There should be 100 hole pixels (10 x 10)
  // There should be (100*100) - (10*10), a 100x100 image with a 10x10 hole in it
  if(numberOfHolePixels == 100 && numberOfValidPixels == 100*100 - 100)
  {
    return EXIT_SUCCESS;
  }
  else
  {
    return EXIT_FAILURE;
  }
}

