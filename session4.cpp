//
// Created by Lander on 5/4/18.
//
#include "session2.cpp"

void
draw_zbuff_line(ZBuffer &zBuffer, img::EasyImage &image, unsigned int x0, unsigned int y0, double z0, unsigned int x1,
                unsigned int y1, double z1, img::Color &color) {
    assert(x0 < image.get_width() && y0 < image.get_height());
    assert(x1 < image.get_width() && y1 < image.get_height());
    double invzvalue;
    double a;
    if (x0 == x1) {
        //special case for x0 == x1
        a = std::max(y0, y1) - std::min(y0, y1);
        for (unsigned int i = std::min(y0, y1); i <= std::max(y0, y1); i++) {
            invzvalue = ((i - std::min(y0, y1)) / a) / z0 + (1 - (i - std::min(y0, y1)) / a) / z1;
            if (invzvalue < zBuffer.zBuffer[i][x0]) {
                (image)(x0, i) = color;
                zBuffer.zBuffer[i][x0] = invzvalue;
            }

        }
    } else if (y0 == y1) {
        //special case for y0 == y1
        a = std::max(x0, x1) - std::min(x0, x1);
        for (unsigned int i = std::min(x0, x1); i <= std::max(x0, x1); i++) {
            invzvalue = ((i - std::min(x0, x1)) / a) / z0 + (1 - (i - std::min(x0, x1)) / a) / z1;
            if (invzvalue < zBuffer.zBuffer[y0][i]) {
                (image)(i, y0) = color;
                zBuffer.zBuffer[y0][i] = invzvalue;
            }
        }
    } else {
        if (x0 > x1) {
            //flip points if x1>x0: we want x0 to have the lowest value
            std::swap(x0, x1);
            std::swap(y0, y1);
        }
        double m = ((double) y1 - (double) y0) / ((double) x1 - (double) x0);
        if (-1.0 <= m && m <= 1.0) {
            for (unsigned int i = 0; i <= (x1 - x0); i++) {
                invzvalue = (double) i / (x1 - x0) / z0 + (1 - ((double) i / (x1 - x0))) / z1;
                if (invzvalue < zBuffer.zBuffer[(unsigned int) round(y0 + m * i)][x0 + i]) {
                    (image)(x0 + i, (unsigned int) round(y0 + m * i)) = color;
                    zBuffer.zBuffer[(unsigned int) round(y0 + m * i)][x0 + i] = invzvalue;
                }

            }
        } else if (m > 1.0) {
            for (unsigned int i = 0; i <= (y1 - y0); i++) {
                invzvalue = (double) i / (y1 - y0) / z0 + (1 - ((double) i / (y1 - y0))) / z1;
                if (invzvalue < zBuffer.zBuffer[y0 + i][(unsigned int) round(x0 + (i / m))]) {
                    (image)((unsigned int) round(x0 + (i / m)), y0 + i) = color;
                    zBuffer.zBuffer[y0 + i][(unsigned int) round(x0 + (i / m))] = invzvalue;
                }
            }
        } else if (m < -1.0) {
            for (unsigned int i = 0; i <= (y0 - y1); i++) {
                invzvalue = (double) i / (y0 - y1) / z0 + (1 - ((double) i / (y0 - y1))) / z1;
                if (invzvalue < zBuffer.zBuffer[y0 - i][(unsigned int) round(x0 - (i / m))]) {
                    (image)((unsigned int) round(x0 - (i / m)), y0 - i) = color;
                    zBuffer.zBuffer[y0 - i][(unsigned int) round(x0 - (i / m))] = invzvalue;
                }
            }
        }
    }
}

img::EasyImage drawLines_zBuff(const Lines2D &lines, const int size, const Color &backgroundColor) {
//    std::chrono::high_resolution_clock::time_point begin = std::chrono::high_resolution_clock::now();

    double xmin = 0, xmax = 0, ymin = 0, ymax = 0;
    double maxRange, Imagex, Imagey, d, DCx, DCy, dx, dy;

    // finding min and max
    for (Line2D l : lines) {
        if (l.p1.x < xmin) {
            xmin = l.p1.x;
        }
        if (l.p2.x < xmin) {
            xmin = l.p2.x;
        }
        if (l.p1.x > xmax) {
            xmax = l.p1.x;
        }
        if (l.p2.x > xmax) {
            xmax = l.p2.x;
        }
        if (l.p1.y < ymin) {
            ymin = l.p1.y;
        }
        if (l.p2.y < ymin) {
            ymin = l.p2.y;
        }
        if (l.p1.y > ymax) {
            ymax = l.p1.y;
        }
        if (l.p2.y > ymax) {
            ymax = l.p2.y;
        }
    }

    //calculating image size
    maxRange = std::max(xmax - xmin, ymax - ymin);
    Imagex = size * ((xmax - xmin) / maxRange);
    Imagey = size * ((ymax - ymin) / maxRange);
    //calculating scaling coefficient
    d = 0.95 * (Imagex / (xmax - xmin));

    //calculating center
    DCx = d / 2 * (xmin + xmax);
    DCy = d / 2 * (ymin + ymax);

    //calculating moving factor
    dx = Imagex / 2 - DCx;
    dy = Imagey / 2 - DCy;

    //making base image
    img::Color bgcolor((uint8_t) convertToColor(backgroundColor.red), (uint8_t) convertToColor(backgroundColor.green),
                       (uint8_t) convertToColor(backgroundColor.blue));
    img::EasyImage image((unsigned) roundToInt(Imagex), (unsigned) roundToInt(Imagey), bgcolor);
    ZBuffer zBuffer((unsigned) roundToInt(Imagex), (unsigned) roundToInt(Imagey));

    //drawing lines from vector
    for (Line2D l : lines) {
        img::Color linecolor((uint8_t) convertToColor(l.color.red), (uint8_t) convertToColor(l.color.green),
                             (uint8_t) convertToColor(l.color.blue));
        draw_zbuff_line(zBuffer, image, (unsigned) roundToInt(l.p1.x * d + dx), (unsigned) roundToInt(l.p1.y * d + dy),
                        l.z1, (unsigned) roundToInt(l.p2.x * d + dx), (unsigned) roundToInt(l.p2.y * d + dy), l.z2,
                        linecolor);
//        image.draw_line((unsigned)roundToInt(l.p1.x * d + dx) , (unsigned)roundToInt(l.p1.y * d + dy), (unsigned)roundToInt(l.p2.x * d + dx), (unsigned)roundToInt(l.p2.y * d + dy), linecolor);
//        draw_zbuff_line(zBuffer, image, (unsigned)std::min(abs(roundToInt(l.p1.x * d + dx)), (int)image.get_width()-1) , (unsigned)std::min(abs(roundToInt(l.p1.y * d + dy)), (int)image.get_height()-1), l.z1, (unsigned)std::min(abs(roundToInt(l.p2.x * d + dx)), (int)image.get_width()-1), (unsigned)std::min(abs(roundToInt(l.p2.y * d + dy)), (int)image.get_height()-1), l.z2, linecolor);
    }


//    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
//    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
//    std::cout << "ZBuffer drawing took " << duration << "µs to complete." << std::endl;
    return image;
};

img::EasyImage zbuffdrawLines3D(const ini::Configuration &configuration) {
    Figures3D ImageFigures;
    const unsigned int size = (unsigned) configuration["General"]["size"].as_int_or_die();
    const std::vector<double> bg = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
    Color bgcolor;
    bgcolor.red = bg[0];
    bgcolor.green = bg[1];
    bgcolor.blue = bg[2];

    generateImage(ImageFigures, configuration);

    if (ImageFigures.Figures.empty()) {
        return img::EasyImage();
    }

    return drawLines_zBuff(ImageFigures.doProjection(), size, bgcolor);
};