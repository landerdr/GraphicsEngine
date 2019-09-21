//
// Created by lander on 2/24/18.
//
#include "easy_image.h"
#include "ini_configuration.h"
#include "l_parser.h"
#include "EngineUtilities.cpp"


#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <cmath>
#include <list>

img::EasyImage generate_IntroColorRectangle(const ini::Configuration &configuration) {
    const unsigned int width = (unsigned)configuration["ImageProperties"]["width"].as_int_or_die();
    const unsigned int height = (unsigned)configuration["ImageProperties"]["height"].as_int_or_die();

    img::EasyImage image(width,height);
    for(unsigned int i = 0; i < width; i++)
    {
        for(unsigned int j = 0; j < height; j++)
        {
            int new_i = roundToInt(i * (255.0/width));
            int new_j = roundToInt(j * (255.0/height));
            image(i,j).red = new_i;
            image(i,j).green = new_j;
            image(i,j).blue = (new_i+new_j)%256;
        }
    }
    return image;
}

img::EasyImage generate_IntroBlocks(const ini::Configuration &configuration) {
    const unsigned int width = (unsigned)configuration["ImageProperties"]["width"].as_int_or_die();
    const unsigned int height = (unsigned)configuration["ImageProperties"]["height"].as_int_or_die();
    const std::vector<double> white = configuration["BlockProperties"]["colorWhite"].as_double_tuple_or_die();
    const std::vector<double> black = configuration["BlockProperties"]["colorBlack"].as_double_tuple_or_die();
    const unsigned int nrXBlocks = (unsigned)configuration["BlockProperties"]["nrXBlocks"].as_int_or_die();
    const unsigned int nrYBlocks = (unsigned)configuration["BlockProperties"]["nrYBlocks"].as_int_or_die();

    const int blockHeight = height/nrXBlocks;
    const int blockWidth = width/nrYBlocks;
    img::EasyImage image(width,height);
    for(unsigned int i = 0; i < width; i++)
    {
        for(unsigned int j = 0; j < height; j++)
        {
            if((i/blockWidth + j/blockHeight)%2 == 0) {
                image(i,j).red = convertToColor(white[0]);
                image(i,j).green = convertToColor(white[1]);
                image(i,j).blue = convertToColor(white[2]);
            } else {
                image(i,j).red = convertToColor(black[0]);
                image(i,j).green = convertToColor(black[1]);
                image(i,j).blue = convertToColor(black[2]);
            }
        }
    }
    return image;
}

img::EasyImage generate_IntroLines(const ini::Configuration &configuration) {
    const unsigned int width = (unsigned)configuration["ImageProperties"]["width"].as_int_or_die();
    const unsigned int height = (unsigned)configuration["ImageProperties"]["height"].as_int_or_die();
    const std::string type = configuration["LineProperties"]["figure"].as_string_or_die();
    const std::vector<double> bg = configuration["LineProperties"]["backgroundcolor"].as_double_tuple_or_die();
    const std::vector<double> line = configuration["LineProperties"]["lineColor"].as_double_tuple_or_die();
    const unsigned int nrLines = (unsigned)configuration["LineProperties"]["nrLines"].as_int_or_die();
    img::Color bgc((uint8_t)convertToColor(bg[0]), (uint8_t)convertToColor(bg[1]), (uint8_t)convertToColor(bg[2]));
    img::EasyImage image(width, height, bgc);

    if(type == "QuarterCircle") {
        const double Hs = (double)height/(nrLines-1);
        const double Ws = static_cast<double>(width)/(nrLines-1);
        img::Color linecolor((uint8_t)convertToColor(line[0]), (uint8_t)convertToColor(line[1]), (uint8_t)convertToColor(line[2]));
        for (unsigned int i=0; i<nrLines; i++) {
            image.draw_line(0, (uint8_t)std::min(roundToInt(i*Hs), (int)height-1), (uint8_t)std::min(roundToInt(i*Ws), (int)width-1), height-1, linecolor);
        }
        return image;
    }
    if(type == "Eye") {
        const double Hs = (double)height/(nrLines/2-1);
        const double Ws = static_cast<double>(width)/(nrLines/2-1);
        img::Color linecolor((uint8_t)convertToColor(line[0]), (uint8_t)convertToColor(line[1]), (uint8_t)convertToColor(line[2]));
        for (unsigned int i=0; i<nrLines/2; i++) {
            image.draw_line(0, (unsigned)std::min(roundToInt(i*Hs), (int)height-1), (unsigned)std::min(roundToInt(i*Ws), (int)width-1), height-1, linecolor);
            image.draw_line(width-1, (unsigned)std::min(roundToInt(i*Hs), (int)height-1), (unsigned)std::min(roundToInt(i*Ws), (int)width-1), 0, linecolor);
        }
        return image;
    }
    if(type == "Diamond") {
        const double Hs = (double)height/(nrLines/2-1);
        const double Ws = static_cast<double>(width)/(nrLines/2-1);
        img::Color linecolor((uint8_t)convertToColor(line[0]), (uint8_t)convertToColor(line[1]), (uint8_t)convertToColor(line[2]));
        for (unsigned int i=0; i<nrLines/4; i++) {
            image.draw_line(width/2, height/2 - std::min(roundToInt((nrLines/4-1-i)*Hs), (int)height-1), width/2 - std::min(roundToInt(i*Ws), (int)width-1), height/2, linecolor);
            image.draw_line(width/2, height/2 + std::min(roundToInt((nrLines/4-1-i)*Hs), (int)height-1), width/2 + std::min(roundToInt(i*Ws), (int)width-1), height/2, linecolor);
            image.draw_line(width/2, height/2 + std::min(roundToInt((nrLines/4-1-i)*Hs), (int)height-1), width/2 - std::min(roundToInt(i*Ws), (int)width-1), height/2, linecolor);
            image.draw_line(width/2, height/2 - std::min(roundToInt((nrLines/4-1-i)*Hs), (int)height-1), width/2 + std::min(roundToInt(i*Ws), (int)width-1), height/2, linecolor);

        }
        return image;
    }
    return img::EasyImage();
}