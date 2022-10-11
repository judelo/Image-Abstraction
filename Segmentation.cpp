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

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include "Segmentation.h"
#include <QColor>
#include <image.h>
#include <misc.h>
#include <filter.h>
#include "segment-graph.h"


// dissimilarity measure between pixels
static inline float diff(image<float> *r, image<float> *g, image<float> *b,
                         int x1, int y1, int x2, int y2) {
    return sqrt(square(imRef(r, x1, y1)-imRef(r, x2, y2)) +
                square(imRef(g, x1, y1)-imRef(g, x2, y2)) +
                square(imRef(b, x1, y1)-imRef(b, x2, y2)));
}


Segmentation::Segmentation(const QImage & input_image )
{

    int width = input_image.width();
    int height = input_image.height();

    input = new image<rgb> ( width, height );

    // smooth each color channel
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {

            QColor color =input_image.pixel( x, y );
            imRef(input, x, y).r = (uchar)color.red();
            imRef(input, x, y).g = (uchar)color.green();
            imRef(input, x, y).b = (uchar)color.blue();
        }
    }

    result = QImage(input_image);
}

Segmentation::~Segmentation(){
    delete [] edges;
    delete [] u;
}

/*
 * Segment an image
 *
 * Returns a color image representing the segmentation.
 *
 * im: image to segment.
 * sigma: to smooth the image.
 * c: constant for treshold function.
 * min_size: minimum component size (enforced by post-processing stage).
 */
QImage Segmentation::segment(float sigma, float c, int min_size) {
    int width = input->width();
    int height = input->height();

    image<float> *r = new image<float>(width, height);
    image<float> *g = new image<float>(width, height);
    image<float> *b = new image<float>(width, height);

    // smooth each color channel
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            imRef(r, x, y) = imRef(input, x, y).r;
            imRef(g, x, y) = imRef(input, x, y).g;
            imRef(b, x, y) = imRef(input, x, y).b;
        }
    }
    image<float> *smooth_r = smooth(r, sigma);
    image<float> *smooth_g = smooth(g, sigma);
    image<float> *smooth_b = smooth(b, sigma);
    delete r;
    delete g;
    delete b;

    // build graph
    edges = new edge[width*height*4];
    num = 0;
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            if (x < width-1) {
                edges[num].a = y * width + x;
                edges[num].b = y * width + (x+1);
                edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x+1, y);
                num++;
            }

            if (y < height-1) {
                edges[num].a = y * width + x;
                edges[num].b = (y+1) * width + x;
                edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x, y+1);
                num++;
            }

            if ((x < width-1) && (y < height-1)) {
                edges[num].a = y * width + x;
                edges[num].b = (y+1) * width + (x+1);
                edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x+1, y+1);
                num++;
            }

            if ((x < width-1) && (y > 0)) {
                edges[num].a = y * width + x;
                edges[num].b = (y-1) * width + (x+1);
                edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x+1, y-1);
                num++;
            }
        }
    }
    delete smooth_r;
    delete smooth_g;
    delete smooth_b;

    // segment
    u = segment_graph(width*height, num, edges, c);

    // post process small components
    for (int i = 0; i < num; i++) {
        int a = u->find(edges[i].a);
        int b = u->find(edges[i].b);
        if ((a != b) && ((u->size(a) < min_size) || (u->size(b) < min_size)))
            u->join(a, b);
    } 

    int num_ccs = u->num_sets();

    segmentation = new image<rgb>(width, height);

    //colors for each component
    std::map<int, float> average_r;
    std::map<int, float> average_g;
    std::map<int, float> average_b;

    std::map<int, int> nb_pixels;

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int comp = u->find(y * width + x);

            average_r[comp] += imRef(input, x, y).r;
            average_g[comp] += imRef(input, x, y).g;
            average_b[comp] += imRef(input, x, y).b;
            nb_pixels[comp]++;
        }
    }

    std::cout << "Num ccs "<< num_ccs << " and " << nb_pixels.size() << std::endl;

    for( std::map<int, int>::iterator it = nb_pixels.begin() ; it != nb_pixels.end() ; it++ ){
        int comp = it->first;
        component_colors[comp] = QColor(int(float(average_r[comp])/nb_pixels[comp]), int(float(average_g[comp])/nb_pixels[comp]), int(float(average_b[comp])/nb_pixels[comp]));       
    }

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int comp = u->find(y * width + x);
            result.setPixel(x, y , qRgb(component_colors[comp].red(), component_colors[comp].green(), component_colors[comp].blue()));
        }
    }

    return result;
}