//
// Created by Lander on 4/4/18.
//

#include "easy_image.h"
#include "ini_configuration.h"
#include "l_parser.h"
#include "Intro.cpp"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <cmath>
#include <list>
#include <stack>

//TODO: Implement stochastic replacement rules

std::string replaceChar(const LParser::LSystem2D &l_system, const char c, unsigned int amount) {
    if (c == '+' || c == '-' || c == '(' || c == ')') {
        std::string w;
        return w + c;
    }
    if (amount == 1) {
        return l_system.get_replacement(c);
    }
    std::string replaced;
    for (char w : l_system.get_replacement(c)) {
        replaced.append(replaceChar(l_system, w, amount - 1));
    }
    return replaced;
}

Lines2D generateLSystem(const LParser::LSystem2D &l_system, const Color &lc) {
    std::stack<Point2D> positions;
    std::stack<double> angles;
    std::vector<Line2D> lines;
    std::string drawstring;
    double currangle = l_system.get_starting_angle();
    double angle = l_system.get_angle();
    Point2D current;
    current.x = 0;
    current.y = 0;

    if (l_system.get_nr_iterations() == 0) {
        drawstring = l_system.get_initiator();
    } else {
        for (char c : l_system.get_initiator()) {
            drawstring.append(replaceChar(l_system, c, l_system.get_nr_iterations()));
        }
    }

    for (char c : drawstring) {
        if (c == '+') {
            currangle += angle;
            continue;
        }
        if (c == '-') {
            currangle -= angle;
            continue;
        }
        if (c == '(') {
            positions.push(current);
            angles.push(currangle);
            continue;
        }
        if (c == ')') {
            current = positions.top();
            positions.pop();
            currangle = angles.top();
            angles.pop();
            continue;
        }

        Point2D next;
        next.x = current.x + std::cos(convertToRadiant(currangle));
        next.y = current.y + std::sin(convertToRadiant(currangle));

        if (l_system.draw(c)) {
            Line2D currline;
            currline.color = lc;
            currline.p1 = current;
            currline.p2 = next;
            lines.push_back(currline);
        }

        current = next;
    }
    return lines;
};

img::EasyImage drawLines(const Lines2D &lines, const int size, const Color &backgroundColor) {

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

    //drawing lines from vector
    for (Line2D l : lines) {
        img::Color linecolor((uint8_t) convertToColor(l.color.red), (uint8_t) convertToColor(l.color.green),
                             (uint8_t) convertToColor(l.color.blue));
        image.draw_line((unsigned) roundToInt(l.p1.x * d + dx), (unsigned) roundToInt(l.p1.y * d + dy),
                        (unsigned) roundToInt(l.p2.x * d + dx), (unsigned) roundToInt(l.p2.y * d + dy), linecolor);
//        image.draw_line((unsigned)std::min(abs(roundToInt(l.p1.x * d + dx)), (int)image.get_width()-1) , (unsigned)std::min(abs(roundToInt(l.p1.y * d + dy)), (int)image.get_height()-1), (unsigned)std::min(abs(roundToInt(l.p2.x * d + dx)), (int)image.get_width()-1), (unsigned)std::min(abs(roundToInt(l.p2.y * d + dy)), (int)image.get_height()-1), linecolor);
    }

    return image;
};

img::EasyImage drawLSystem(const ini::Configuration &configuration) {
    LParser::LSystem2D l_system;

    const unsigned int size = (unsigned) configuration["General"]["size"].as_int_or_die();
    const std::vector<double> bg = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
    const std::vector<double> lc = configuration["2DLSystem"]["color"].as_double_tuple_or_die();

    Color linecolor;
    linecolor.red = lc[0];
    linecolor.green = lc[1];
    linecolor.blue = lc[2];
    Color bgcolor;
    bgcolor.red = bg[0];
    bgcolor.green = bg[1];
    bgcolor.blue = bg[2];

    std::ifstream input_stream(configuration["2DLSystem"]["inputfile"].as_string_or_die());
    input_stream >> l_system;
    input_stream.close();

    return drawLines(generateLSystem(l_system, linecolor), size, bgcolor);
};









//img::EasyImage Lines(const ini::Configuration &configuration) {
//
//    const unsigned int size = (unsigned)configuration["General"]["size"].as_int_or_die();
//    const std::vector<double> bg = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
//    const std::vector<double> lc = configuration["General"]["color"].as_double_tuple_or_die();
//    const std::vector<int> lines = configuration["General"]["lines"].as_int_tuple_or_die();
//
//    Color linecolor; linecolor.red = lc[0]; linecolor.green = lc[1]; linecolor.blue = lc[2];
//    Color bgcolor; bgcolor.red = bg[0]; bgcolor.green = bg[1]; bgcolor.blue = bg[2];
//
//    Lines2D ProcessedLines;
//    int i=0;
//    while ((unsigned)i+1<lines.size()) {
//        Point2D p1; p1.x = lines[i]; p1.y = lines[i+1];
//        Point2D p2; p2.x = lines[i+2]; p2.y = lines[i+3];
//        Line2D l; l.color = linecolor; l.p1 = p1; l.p2 = p2;
//        ProcessedLines.push_back(l);
//        i += 4;
//    }
//
//    return drawLines(ProcessedLines, size, bgcolor);
//};