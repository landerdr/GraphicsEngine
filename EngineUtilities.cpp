//
// Created by Lander on 4/4/18.
//

#include <cmath>
#include <chrono>
#include <vector>
#include <array>
#include <limits>
#include <cassert>
#include "vector3d.h"
#include "easy_image.h"

//intro

inline int roundToInt(double d) {
    return static_cast<int> (std::round(d));
};

inline int convertToColor(double d) {
    if (d < 0) return 0;
    if (d > 1) return 255;
    return roundToInt(d * 255);
};

inline double convertToRadiant(double degrees) {
    return degrees * M_PI / 180;
};

//session 1

class Color {
public:
    double red;
    double green;
    double blue;
};

class Point2D {
public:
    double x;
    double y;
};

class Line2D {
public:
    Point2D p1;
    Point2D p2;
    Color color;
    double z1;
    double z2;
};

typedef std::vector<Line2D> Lines2D;


//session 2

const Matrix generateScaleMatrix(const double scale) {
    Matrix m;
    for (int i = 1; i < 4; i++) {
        m(i, i) = scale;
    }
    return m;
};

const Matrix generateXRotationMatrix(const double angle) {
    Matrix m;
    m(2, 2) = std::cos(angle);
    m(2, 3) = std::sin(angle);
    m(3, 2) = -std::sin(angle);
    m(3, 3) = std::cos(angle);
    return m;
};

const Matrix generateYRotationMatrix(const double angle) {
    Matrix m;
    m(1, 1) = std::cos(angle);
    m(1, 3) = -std::sin(angle);
    m(3, 1) = std::sin(angle);
    m(3, 3) = std::cos(angle);
    return m;
};

const Matrix generateZRotationMatrix(const double angle) {
    Matrix m;
    m(1, 1) = std::cos(angle);
    m(1, 2) = std::sin(angle);
    m(2, 1) = -std::sin(angle);
    m(2, 2) = std::cos(angle);
    return m;
};

const Matrix generateTranslationMatrix(const Vector3D &translation) {
    Matrix m;
    m(4, 1) = translation.x;
    m(4, 2) = translation.y;
    m(4, 3) = translation.z;
    return m;
};

const Matrix generateEyePerspectiveMatrix(const Vector3D &eye) {
    Matrix m;
    double r = eye.length();
    double theta = std::atan2(eye.y, eye.x);
    double phi = std::acos(eye.z / r);

    m(1, 1) = -std::sin(theta);
    m(1, 2) = -std::cos(theta) * std::cos(phi);
    m(1, 3) = std::cos(theta) * std::sin(phi);
    m(2, 1) = std::cos(theta);
    m(2, 2) = -std::sin(theta) * std::cos(phi);
    m(2, 3) = std::sin(theta) * std::sin(phi);
    m(3, 2) = std::sin(phi);
    m(3, 3) = std::cos(phi);
    m(4, 3) = -r;

    return m;
};

Point2D doPProjection(const Vector3D &point, const double d) {
    Point2D p;

    p.x = -(d * point.x) / point.z;
    p.y = -(d * point.y) / point.z;

    return p;
};

class Face {
public:
    std::vector<int> index;

    void addIndexes(const std::vector<int> &v) {
        for (int i : v) {
            index.push_back(i);
        }
    }
};

class Figure {
public:
    std::vector<Vector3D> points;
    std::vector<Face> faces;
    //Color color;
    Color ambientReflection;
    Color diffuseReflection;
    Color specularReflection;
    double reflectionCoefficient;


    Lines2D doProjection() {
        std::vector<Line2D> Lines;

        for (Face &f : faces) {
            for (unsigned int i = 0; i < f.index.size(); i++) {
                Line2D l;
                l.p1 = doPProjection(points[f.index[i]], 1);
                l.z1 = points[f.index[i]].z;
                l.p2 = doPProjection(points[f.index[(i + 1) % f.index.size()]], 1);
                l.z2 = points[f.index[(i + 1) % f.index.size()]].z;

                l.color = ambientReflection;
                Lines.push_back(l);
            }
        }
        return Lines;
    }

    void applyTransformation(const Matrix &m) {
        for (Vector3D &v : points) {
            v *= m;
        }
    }

    void generateCube() {
        Vector3D point;
        point.x = 1;
        point.y = -1;
        point.z = -1;
        points.push_back(point);
        point.x = -1;
        point.y = 1;
        points.push_back(point);
        point.x = 1;
        point.z = 1;
        points.push_back(point);
        point.x = -1;
        point.y = -1;
        points.push_back(point);
        point.x = 1;
        point.y = 1;
        point.z = -1;
        points.push_back(point);
        point.x = -1;
        point.y = -1;
        points.push_back(point);
        point.x = 1;
        point.z = 1;
        points.push_back(point);
        point.x = -1;
        point.y = 1;
        points.push_back(point);

        Face g;
        g.index = {0, 4, 2, 6};
        faces.push_back(g);
        g.index = {4, 1, 7, 2};
        faces.push_back(g);
        g.index = {1, 5, 3, 7};
        faces.push_back(g);
        g.index = {5, 0, 6, 3};
        faces.push_back(g);
        g.index = {6, 2, 7, 3};
        faces.push_back(g);
        g.index = {0, 5, 1, 4};
        faces.push_back(g);
    };

    void generateTetrahedron() {
        Vector3D point;
        point.x = 1;
        point.y = -1;
        point.z = -1;
        points.push_back(point);
        point.x = -1;
        point.y = 1;
        points.push_back(point);
        point.x = 1;
        point.z = 1;
        points.push_back(point);
        point.x = -1;
        point.y = -1;
        points.push_back(point);

        Face g;
        g.index = {0, 1, 2};
        faces.push_back(g);
        g.index = {1, 3, 2};
        faces.push_back(g);
        g.index = {0, 3, 1};
        faces.push_back(g);
        g.index = {0, 2, 3};
        faces.push_back(g);
    };

    void generateOctahedron() {
        Vector3D point;
        point.x = 1;
        point.y = 0;
        point.z = 0;
        points.push_back(point);
        point.x = 0;
        point.y = 1;
        points.push_back(point);
        point.x = -1;
        point.y = 0;
        points.push_back(point);
        point.x = 0;
        point.y = -1;
        points.push_back(point);
        point.y = 0;
        point.z = -1;
        points.push_back(point);
        point.z = 1;
        points.push_back(point);

        Face g;
        g.index = {0, 1, 5};
        faces.push_back(g);
        g.index = {1, 2, 5};
        faces.push_back(g);
        g.index = {2, 3, 5};
        faces.push_back(g);
        g.index = {3, 0, 5};
        faces.push_back(g);
        g.index = {1, 0, 4};
        faces.push_back(g);
        g.index = {2, 1, 4};
        faces.push_back(g);
        g.index = {3, 2, 4};
        faces.push_back(g);
        g.index = {0, 3, 4};
        faces.push_back(g);
    };

    void generateIcosahedron() {
        Vector3D point;
        point.x = 0;
        point.y = 0;
        point.z = sqrt(5) / 2;
        points.push_back(point);
        point.z = 0.5;
        for (unsigned int i = 2; i < 7; i++) {
            point.x = std::cos((i - 2) * 2 * M_PI / 5);
            point.y = std::sin((i - 2) * 2 * M_PI / 5);
            points.push_back(point);
        }
        point.z = -0.5;
        for (unsigned int i = 7; i < 12; i++) {
            point.x = std::cos(M_PI / 5 * (1 + (i - 7) * 2));
            point.y = std::sin(M_PI / 5 * (1 + (i - 7) * 2));
            points.push_back(point);
        }
        point.x = 0;
        point.y = 0;
        point.z = -sqrt(5) / 2;
        points.push_back(point);

        Face g;
        g.index = {0, 1, 2};
        faces.push_back(g);
        g.index = {0, 2, 3};
        faces.push_back(g);
        g.index = {0, 3, 4};
        faces.push_back(g);
        g.index = {0, 4, 5};
        faces.push_back(g);
        g.index = {0, 5, 1};
        faces.push_back(g);
        g.index = {1, 6, 2};
        faces.push_back(g);
        g.index = {2, 6, 7};
        faces.push_back(g);
        g.index = {2, 7, 3};
        faces.push_back(g);
        g.index = {3, 7, 8};
        faces.push_back(g);
        g.index = {3, 8, 4};
        faces.push_back(g);
        g.index = {4, 8, 9};
        faces.push_back(g);
        g.index = {4, 9, 5};
        faces.push_back(g);
        g.index = {5, 9, 10};
        faces.push_back(g);
        g.index = {5, 10, 1};
        faces.push_back(g);
        g.index = {1, 10, 6};
        faces.push_back(g);
        g.index = {11, 7, 6};
        faces.push_back(g);
        g.index = {11, 8, 7};
        faces.push_back(g);
        g.index = {11, 9, 8};
        faces.push_back(g);
        g.index = {11, 10, 9};
        faces.push_back(g);
        g.index = {11, 6, 10};
        faces.push_back(g);
    };

    void generateDodecahedron() {
        Figure g;
        g.generateIcosahedron();
        Vector3D point;
        for (unsigned int i = 0; i < 20; i++) {
            point = (g.points[g.faces[i].index[0]] + g.points[g.faces[i].index[1]] + g.points[g.faces[i].index[2]]) / 3;
            points.push_back(point);
        }

        Face face;
        face.index = {0, 1, 2, 3, 4};
        faces.push_back(face);
        face.index = {0, 5, 6, 7, 1};
        faces.push_back(face);
        face.index = {1, 7, 8, 9, 2};
        faces.push_back(face);
        face.index = {2, 9, 10, 11, 3};
        faces.push_back(face);
        face.index = {3, 11, 12, 13, 4};
        faces.push_back(face);
        face.index = {4, 13, 14, 5, 0};
        faces.push_back(face);
        face.index = {19, 18, 17, 16, 15};
        faces.push_back(face);
        face.index = {19, 14, 13, 12, 18};
        faces.push_back(face);
        face.index = {18, 12, 11, 10, 17};
        faces.push_back(face);
        face.index = {17, 10, 9, 8, 16};
        faces.push_back(face);
        face.index = {16, 8, 7, 6, 15};
        faces.push_back(face);
        face.index = {15, 6, 5, 14, 19};
        faces.push_back(face);

    };

    void generateCone(const unsigned int n, const double h) {
        Vector3D point;
        point.z = 0;
        for (unsigned int i = 0; i < n; i++) {
            point.x = std::cos(2 * M_PI * i / n);
            point.y = std::sin(2 * M_PI * i / n);
            points.push_back(point);
        }
        point.x = 0;
        point.y = 0;
        point.z = h;
        points.push_back(point);

        Face g;
        for (unsigned int i = 0; i < n; i++) {
            g.index = {i, (i + 1) % n, n};
            faces.push_back(g);
        }
        g.index = {};
        for (int i = n - 1; i >= 0; i--) {
            g.index.push_back(i);
        }
        faces.push_back(g);
    };

    void generateCylinder(const unsigned int n, const double h) {
        Vector3D point;
        std::vector<Vector3D> points_back;
        for (unsigned int i = 0; i < n; i++) {
            point.x = std::cos(2 * M_PI * i / n);
            point.y = std::sin(2 * M_PI * i / n);
            point.z = 0;
            points.push_back(point);
            point.z = h;
            points_back.push_back(point);
        }
        points.insert(points.end(), points_back.begin(), points_back.end());

        Face g;
        Face k;
        for (unsigned int i = 0; i < n; i++) {
            g.index = {i, (i + 1) % n, n + (i + 1) % n, i + n};
            faces.push_back(g);
        }
        g.index = {};
        for (int i = n - 1; i >= 0; i--) {
            g.index.push_back(i);
            k.index.push_back(i + n);
        }
        faces.push_back(g);
        faces.push_back(k);

    };

    void generateSphere(const int n) {
        std::vector<Face> a;
        Face g;
        Vector3D point;
        generateIcosahedron();
        for (unsigned int i = 0; i < n; i++) {
            a = {};
            for (Face &f : faces) {
                for (unsigned int j = 0; j < 3; ++j) {
                    point = (points[f.index[j]] + points[f.index[(j + 1) % 3]]) / 2;
                    points.push_back(point);
                }
                g.index = {f.index[0], points.size() - 3, points.size() - 1};
                a.push_back(g);
                g.index = {f.index[1], points.size() - 2, points.size() - 3};
                a.push_back(g);
                g.index = {f.index[2], points.size() - 1, points.size() - 2};
                a.push_back(g);
                g.index = {points.size() - 3, points.size() - 2, points.size() - 1};
                a.push_back(g);
            }
            faces = a;
        }
        for (Vector3D &p : points) {
            p.normalise();
        }
    };

    void generateTorus(const double r, const double R, const int n, const int m) {
        Vector3D point;
        Face g;
        for (unsigned int i = 0; i < n; i++) {
            for (unsigned int j = 0; j < m; j++) {
                point.x = (R + r * std::cos(2 * M_PI * j / m)) * std::cos(2 * M_PI * i / n);
                point.y = (R + r * std::cos(2 * M_PI * j / m)) * std::sin(2 * M_PI * i / n);
                point.z = r * std::sin(2 * M_PI * j / m);
                points.push_back(point);

                g.index = {m * i + j, m * ((i + 1) % n) + j, m * ((i + 1) % n) + (j + 1) % m, m * i + (j + 1) % m};
                faces.push_back(g);
            }
        }
    };

    void triangulate() {
        std::vector<Face> a;
        Face g;
        for (Face &f : faces) {
            if (f.index.size() > 3) {
                for (unsigned int i = 1; i < f.index.size() - 1; i++) {
                    g.index = {f.index[0], f.index[i], f.index[i + 1]};
                    a.push_back(g);
                }
            } else {
                a.push_back(f);
            }
        }
        faces = a;
    };

};

class Figures3D {
public:
    std::vector<Figure> Figures;

    Lines2D doProjection() {
        std::vector<Line2D> Lines;
        for (Figure f : Figures) {
            for (Line2D &l : f.doProjection()) {
                Lines.push_back(l);
            }
        }

        return Lines;
    }

    void applyTransformation(const Matrix &m) {
        for (Figure &f : Figures) {
            f.applyTransformation(m);
        }
    }
};

// Session 4

class ZBuffer {
public:
    std::vector<std::vector<double>> zBuffer;

    ZBuffer(unsigned int width, unsigned int height) {
        std::chrono::high_resolution_clock::time_point begin = std::chrono::high_resolution_clock::now();
        std::vector<double> row;
        row.assign(width, std::numeric_limits<double>::infinity());
        zBuffer.assign(height, row);
        std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
        std::cout << "ZBuffer generation took " << duration << "µs to complete." << std::endl;
    }
};

// Session 6

void generateFractal(Figure &fig, Figures3D &fractal, const int nr_iterations, const double scale) {
    Matrix scalem = generateScaleMatrix(1 / scale);
    std::vector<Figure> baseFig = {fig};
    std::vector<Figure> appliedFigs;
    for (unsigned int i = 0; i < nr_iterations; ++i) {
        appliedFigs = {};
        for (Figure &f : baseFig) {
            for (unsigned int p = 0; p < f.points.size(); ++p) {
                Figure newFig = f;
                newFig.applyTransformation(scalem);
                Matrix trans = generateTranslationMatrix(f.points[p] - newFig.points[p]);
                newFig.applyTransformation(trans);
                appliedFigs.push_back(newFig);
            }
        }
        baseFig = appliedFigs;
    }
    fractal.Figures.insert(std::end(fractal.Figures), std::begin(baseFig), std::end(baseFig));
};

// Session 7

class Light {
public:
    //de ambiente licht component
    Color ambientLight;
    //de diffuse licht component
    Color diffuseLight;
    //de diffuse licht component
    Color specularLight;
    //de richting waarin het
    //licht schijnt
    Vector3D ldVector;
    bool infLight;
    //de locatie van de puntbron
    Vector3D location;
    bool posLight;
    bool specularLighting;
};

typedef std::vector<Light> Lights3D;
