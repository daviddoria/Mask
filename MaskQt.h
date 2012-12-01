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

#ifndef MaskQt_H
#define MaskQt_H

#include "Mask.h"

// Qt
#include <QImage>

namespace MaskQt
{
  QImage GetQtImage(const Mask* const mask);

  QImage GetQtImage(const Mask* const mask, const itk::ImageRegion<2>& region);

  QImage SetPixelsToTransparent(QImage image, const Mask* const mask,
                                HoleMaskPixelTypeEnum pixelValue);
}

#endif
