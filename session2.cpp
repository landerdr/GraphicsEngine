//
// Created by Lander on 4/4/18.
//

#include "easy_image.h"
#include "ini_configuration.h"
#include "l_parser.h"
#include "vector3d.h"
#include "session3.cpp"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <cmath>
#include <list>

void generateImage(Figures3D &ImageFigures, const ini::Configuration &configuration) {
    const unsigned int nrFigures = (unsigned) configuration["General"]["nrFigures"].as_int_or_die();
    const std::vector<double> eyep = configuration["General"]["eye"].as_double_tuple_or_die();
    Vector3D eye;
    eye.x = eyep[0];
    eye.y = eyep[1];
    eye.z = eyep[2];

    Matrix eyetrans = generateEyePerspectiveMatrix(eye);

    for (unsigned int i = 0; i < nrFigures; ++i) {
        Figure singleFigure;
        const std::string type = configuration["Figure" + std::to_string(i)]["type"].as_string_or_die();
        const double rotateX = configuration["Figure" + std::to_string(i)]["rotateX"].as_double_or_die();
        const double rotateY = configuration["Figure" + std::to_string(i)]["rotateY"].as_double_or_die();
        const double rotateZ = configuration["Figure" + std::to_string(i)]["rotateZ"].as_double_or_die();
        const double scale = configuration["Figure" + std::to_string(i)]["scale"].as_double_or_die();
        const std::vector<double> center = configuration["Figure" +
                                                         std::to_string(i)]["center"].as_double_tuple_or_die();
        Vector3D centerp;
        centerp.x = center[0];
        centerp.y = center[1];
        centerp.z = center[2];
        std::vector<double> lc;

        if (configuration["General"]["type"].as_string_or_die() != "LightedZBuffering") {
            lc = configuration["Figure" + std::to_string(i)]["color"].as_double_tuple_or_die();
        } else {
            lc = configuration["Figure" + std::to_string(i)]["ambientReflection"].as_double_tuple_or_die();
        }

        singleFigure.ambientReflection.red = lc[0];
        singleFigure.ambientReflection.green = lc[1];
        singleFigure.ambientReflection.blue = lc[2];

        if (configuration["Figure" + std::to_string(i)]["diffuseReflection"].as_double_tuple_if_exists(lc)) {
            singleFigure.diffuseReflection.red = lc[0];
            singleFigure.diffuseReflection.green = lc[1];
            singleFigure.diffuseReflection.blue = lc[2];
        } else {
            singleFigure.diffuseReflection.red = 0;
            singleFigure.diffuseReflection.green = 0;
            singleFigure.diffuseReflection.blue = 0;
        }

        if (configuration["Figure" + std::to_string(i)]["specularReflection"].as_double_tuple_if_exists(lc)) {
            singleFigure.specularReflection.red = lc[0];
            singleFigure.specularReflection.green = lc[1];
            singleFigure.specularReflection.blue = lc[2];
            singleFigure.reflectionCoefficient = configuration["Figure" + std::to_string(i)]["reflectionCoefficient"].as_double_or_die();
        } else {
            singleFigure.specularReflection.red = 0;
            singleFigure.specularReflection.green = 0;
            singleFigure.specularReflection.blue = 0;
            singleFigure.reflectionCoefficient = 0;
        }

        Matrix adj = generateScaleMatrix(scale) * generateXRotationMatrix(convertToRadiant(rotateX)) *
                     generateYRotationMatrix(convertToRadiant(rotateY)) *
                     generateZRotationMatrix(convertToRadiant(rotateZ)) * generateTranslationMatrix(centerp) * eyetrans;

        if (type == "LineDrawing") {
            const unsigned int nrPoints = (unsigned) configuration["Figure" +
                                                                   std::to_string(i)]["nrPoints"].as_int_or_die();
            const unsigned int nrLines = (unsigned) configuration["Figure" +
                                                                  std::to_string(i)]["nrLines"].as_int_or_die();
            for (unsigned int j = 0; j < nrPoints; ++j) {
                const std::vector<double> vpoint = configuration["Figure" + std::to_string(i)]["point" + std::to_string(
                        j)].as_double_tuple_or_die();
                Vector3D point;
                point.x = vpoint[0];
                point.y = vpoint[1];
                point.z = vpoint[2];
                singleFigure.points.push_back(point);
            }
            for (unsigned int j = 0; j < nrLines; ++j) {
                const std::vector<int> vline = configuration["Figure" + std::to_string(i)]["line" + std::to_string(
                        j)].as_int_tuple_or_die();
                Face face;
                face.addIndexes(vline);
                singleFigure.faces.push_back(face);
            }
            singleFigure.applyTransformation(adj);
            ImageFigures.Figures.push_back(singleFigure);
        } else if (type == "Cube") {
            singleFigure.generateCube();
            singleFigure.applyTransformation(adj);
            ImageFigures.Figures.push_back(singleFigure);
        } else if (type == "FractalCube") {
            const double fractscale = configuration["Figure" + std::to_string(i)]["fractalScale"].as_double_or_die();
            const unsigned int nr_it = (unsigned) configuration["Figure" +
                                                                std::to_string(i)]["nrIterations"].as_int_or_die();
            singleFigure.generateCube();
            singleFigure.applyTransformation(adj);
            generateFractal(singleFigure, ImageFigures, nr_it, fractscale);
        } else if (type == "Tetrahedron") {
            singleFigure.generateTetrahedron();
            singleFigure.applyTransformation(adj);
            ImageFigures.Figures.push_back(singleFigure);
        } else if (type == "FractalTetrahedron") {
            const double fractscale = configuration["Figure" + std::to_string(i)]["fractalScale"].as_double_or_die();
            const unsigned int nr_it = (unsigned) configuration["Figure" +
                                                                std::to_string(i)]["nrIterations"].as_int_or_die();
            singleFigure.generateTetrahedron();
            singleFigure.applyTransformation(adj);
            generateFractal(singleFigure, ImageFigures, nr_it, fractscale);
        } else if (type == "Icosahedron") {
            singleFigure.generateIcosahedron();
            singleFigure.applyTransformation(adj);
            ImageFigures.Figures.push_back(singleFigure);
        } else if (type == "FractalIcosahedron") {
            const double fractscale = configuration["Figure" + std::to_string(i)]["fractalScale"].as_double_or_die();
            const unsigned int nr_it = (unsigned) configuration["Figure" +
                                                                std::to_string(i)]["nrIterations"].as_int_or_die();
            singleFigure.generateIcosahedron();
            singleFigure.applyTransformation(adj);
            generateFractal(singleFigure, ImageFigures, nr_it, fractscale);
        } else if (type == "Octahedron") {
            singleFigure.generateOctahedron();
            singleFigure.applyTransformation(adj);
            ImageFigures.Figures.push_back(singleFigure);
        } else if (type == "FractalOctahedron") {
            const double fractscale = configuration["Figure" + std::to_string(i)]["fractalScale"].as_double_or_die();
            const unsigned int nr_it = (unsigned) configuration["Figure" +
                                                                std::to_string(i)]["nrIterations"].as_int_or_die();
            singleFigure.generateOctahedron();
            singleFigure.applyTransformation(adj);
            generateFractal(singleFigure, ImageFigures, nr_it, fractscale);
        } else if (type == "Dodecahedron") {
            singleFigure.generateDodecahedron();
            singleFigure.applyTransformation(adj);
            ImageFigures.Figures.push_back(singleFigure);
        } else if (type == "FractalDodecahedron") {
            const double fractscale = configuration["Figure" + std::to_string(i)]["fractalScale"].as_double_or_die();
            const unsigned int nr_it = (unsigned) configuration["Figure" +
                                                                std::to_string(i)]["nrIterations"].as_int_or_die();
            singleFigure.generateDodecahedron();
            singleFigure.applyTransformation(adj);
            generateFractal(singleFigure, ImageFigures, nr_it, fractscale);
        } else if (type == "Cone") {
            const unsigned int n = (unsigned) configuration["Figure" + std::to_string(i)]["n"].as_int_or_die();
            const double h = configuration["Figure" + std::to_string(i)]["height"].as_double_or_die();
            singleFigure.generateCone(n, h);
            singleFigure.applyTransformation(adj);
            ImageFigures.Figures.push_back(singleFigure);
        } else if (type == "Cylinder") {
            const unsigned int n = (unsigned) configuration["Figure" + std::to_string(i)]["n"].as_int_or_die();
            const double h = configuration["Figure" + std::to_string(i)]["height"].as_double_or_die();
            singleFigure.generateCylinder(n, h);
            singleFigure.applyTransformation(adj);
            ImageFigures.Figures.push_back(singleFigure);
        } else if (type == "Torus") {
            const double r = configuration["Figure" + std::to_string(i)]["r"].as_double_or_die();
            const double R = configuration["Figure" + std::to_string(i)]["R"].as_double_or_die();
            const unsigned int n = (unsigned) configuration["Figure" + std::to_string(i)]["n"].as_int_or_die();
            const unsigned int m = (unsigned) configuration["Figure" + std::to_string(i)]["m"].as_int_or_die();
            singleFigure.generateTorus(r, R, n, m);
            singleFigure.applyTransformation(adj);
            ImageFigures.Figures.push_back(singleFigure);
        } else if (type == "Sphere") {
            const unsigned int n = (unsigned) configuration["Figure" + std::to_string(i)]["n"].as_int_or_die();
            singleFigure.generateSphere(n);
            singleFigure.applyTransformation(adj);
            ImageFigures.Figures.push_back(singleFigure);
        } else if (type == "3DLSystem") {
            LParser::LSystem3D l_system;
            std::ifstream input_stream(configuration["Figure" + std::to_string(i)]["inputfile"].as_string_or_die());
            input_stream >> l_system;
            input_stream.close();

            generateLSystem3D(singleFigure, l_system);
            singleFigure.applyTransformation(adj);
            ImageFigures.Figures.push_back(singleFigure);
        } else {
            return;
        }
    }

}

img::EasyImage drawLines3D(const ini::Configuration &configuration) {
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

    return drawLines(ImageFigures.doProjection(), size, bgcolor);
};