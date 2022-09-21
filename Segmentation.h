/*
 Module to implement segmentation on input image.

Copyright (C) 2022, noura.faraj@umontpellier.fr, lucia.bouza-heguerte@u-paris.fr, julie.delon@u-paris.fr

This program is free software: you can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License along with this program. If not, see <http://www.gnu.org/licenses/>
*/

#ifndef SEGMENTATION_H
#define SEGMENTATION_H

#include <image.h>
#include <misc.h>
#include <pnmfile.h>
#include <map>
#include <QImage>
#include "disjoint-set.h"

class Segmentation
{
public:
    Segmentation(const QImage & input_image);
    ~Segmentation();
    QImage segment(float sigma, float c, int min_size);
    QImage getResult( );
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
