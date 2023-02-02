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
};

#endif // SEGMENTATION_H
