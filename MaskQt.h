#ifndef MaskQt_H
#define MaskQt_H

#include "Mask.h"

// Qt
#include <QImage>

namespace MaskQt
{
QImage GetQtImage(const Mask* const mask);

QImage GetQtImage(const Mask* const mask, const itk::ImageRegion<2>& region);
}

#endif
