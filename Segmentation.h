/*
Copyright (C) 2006 Pedro Felzenszwalb

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
*/

#ifndef SEGMENTATION_H
#define SEGMENTATION_H

#include <image.h>
#include <misc.h>
#include <pnmfile.h>

#include <QImage>
#include "disjoint-set.h"

#include "mw3.h"

class Segmentation
{
public:
    Segmentation(const QImage & input_image);
    ~Segmentation();
    QImage segment(float sigma, float c, int min_size);
    QImage removeRegionUnder(Point_plane ArrayPixelsMask, int len_ArrayPixelsMask);
    QImage getResult();
protected:
    image<rgb> *input;
    image<rgb> *segmentation;
    Universe *u ;
    edge *edges;
    QImage result;
    int num;
    std::map<int, QColor> component_colors;
    std::map<int, bool> removed_regions;
};

#endif // SEGMENTATION_H
