//
// Created by Lander on 5/5/18.
//

#include "session4.cpp"

void
draw_zbuff_triangle(ZBuffer &zBuffer, img::EasyImage &image, Vector3D &A, Vector3D &B, Vector3D &C, double d, double dx,
                    double dy, Color ambientReflection, Color diffuseReflection, Color specularReflection,
                    double reflectionCoeff, Lights3D &lights) {
    int ymax, ymin;
    double ZG_inv, dxdz, dydz, k, z_inv;

    // AB, AC, (a,b,c)
    Vector3D u, v, w, n;

    // Projected points
    Point2D A_proj = doPProjection(A, d);
    A_proj.x += dx;
    A_proj.y += dy;
    Point2D B_proj = doPProjection(B, d);
    B_proj.x += dx;
    B_proj.y += dy;
    Point2D C_proj = doPProjection(C, d);
    C_proj.x += dx;
    C_proj.y += dy;
    Point2D G_proj;
    G_proj.x = (A_proj.x + B_proj.x + C_proj.x) / 3;
    G_proj.y = (A_proj.y + B_proj.y + C_proj.y) / 3;

    // Extrema
    ymax = roundToInt(std::max(std::max(A_proj.y, B_proj.y), C_proj.y) - 0.5);
    ymin = roundToInt(std::min(std::min(A_proj.y, B_proj.y), C_proj.y) + 0.5);

    u = B - A;
    v = C - A;

//    w = Vector3D::vector(u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x);
    w = Vector3D::cross(u, v);

    k = w.x * A.x + w.y * A.y + w.z * A.z;

    n = w/k;
    n.normalise();

    // Calculating ambient/diffuse colour
    Color ambientColor;
    ambientColor.red = 0;
    ambientColor.green = 0;
    ambientColor.blue = 0;

    Color diffuseColor;
    diffuseColor.red = 0;
    diffuseColor.green = 0;
    diffuseColor.blue = 0;

    Color specularColor;
    specularColor.red = 0;
    specularColor.green = 0;
    specularColor.blue = 0;

    for (unsigned int i = 0; i < lights.size(); i++) {
        ambientColor.red += lights[i].ambientLight.red;
        ambientColor.green += lights[i].ambientLight.green;
        ambientColor.blue += lights[i].ambientLight.blue;

        if (lights[i].infLight) {
            Vector3D ld = lights[i].ldVector;
            ld.normalise();
            double a = ld.x*n.x + ld.y*n.y + ld.z*n.z;
//            a *= -1;
//            std::cout << a << std::endl;
            if (a>0) {
                diffuseColor.red += lights[i].diffuseLight.red*a;
                diffuseColor.green += lights[i].diffuseLight.green*a;
                diffuseColor.blue += lights[i].diffuseLight.blue*a;


            }
        }

    }

    Color basecolor;
    basecolor.red = ambientReflection.red * ambientColor.red + diffuseReflection.red * diffuseColor.red;
    basecolor.green = ambientReflection.green * ambientColor.green + diffuseReflection.green * diffuseColor.green;
    basecolor.blue = ambientReflection.blue * ambientColor.blue + diffuseReflection.blue * diffuseColor.blue;

    dxdz = -w.x / d / k;
    dydz = -w.y / d / k;

    ZG_inv = (1 / A.z + 1 / B.z + 1 / C.z) / 3;


    for (int yi = ymin; yi <= ymax; ++yi) {
        double xl_AB = std::numeric_limits<double>::infinity(), xl_BC = std::numeric_limits<double>::infinity(), xl_AC = std::numeric_limits<double>::infinity(),
                xr_AB = -std::numeric_limits<double>::infinity(), xr_BC = -std::numeric_limits<double>::infinity(), xr_AC = -std::numeric_limits<double>::infinity();
        int xl = 0, xr = 0;

//        std::cout << (yi-A_proj.y)*(yi-B_proj.y) << std::endl;
        if (A_proj.y != B_proj.y && (yi - A_proj.y) * (yi - B_proj.y) <= 0) {
            xl_AB = B_proj.x + (A_proj.x - B_proj.x) * (yi - B_proj.y) / (A_proj.y - B_proj.y);
            xr_AB = xl_AB;
        }
//        std::cout << (yi-B_proj.y)*(yi-C_proj.y) << std::endl;
        if (B_proj.y != C_proj.y && (yi - B_proj.y) * (yi - C_proj.y) <= 0) {
            xl_BC = C_proj.x + (B_proj.x - C_proj.x) * (yi - C_proj.y) / (B_proj.y - C_proj.y);
            xr_BC = xl_BC;
        }
//        std::cout << (yi-A_proj.y)*(yi-C_proj.y) << std::endl;
        if (A_proj.y != C_proj.y && (yi - A_proj.y) * (yi - C_proj.y) <= 0) {
            xl_AC = C_proj.x + (A_proj.x - C_proj.x) * (yi - C_proj.y) / (A_proj.y - C_proj.y);
            xr_AC = xl_AC;
        }
//        std::cout << xl_AB << " : " << xl_BC << " : " << xl_AC << std::endl;
//        std::cout << xr_AB << " : " << xr_BC << " : " << xr_AC << std::endl;
        xl = roundToInt(std::min(std::min(xl_AB, xl_BC), xl_AC) + 0.5);
        xr = roundToInt(std::max(std::max(xr_AB, xr_BC), xr_AC) - 0.5);

        for (int xi = xl; xi <= xr; ++xi) {
            z_inv = 1.0001 * ZG_inv + (xi - G_proj.x) * dxdz + (yi - G_proj.y) * dydz;
            if (z_inv < zBuffer.zBuffer[yi][xi]) {
                zBuffer.zBuffer[yi][xi] = z_inv;

                Vector3D p = Vector3D::point((xi-dx)*-1/z_inv/d, (yi-dy)*-1/z_inv/d, 1/z_inv);

                Color diffusePointColor;
                diffusePointColor.red = 0;
                diffusePointColor.green = 0;
                diffusePointColor.blue = 0;

                Color specularPointColor;
                specularPointColor.red = 0;
                specularPointColor.green = 0;
                specularPointColor.blue = 0;

                for (unsigned int i = 0; i < lights.size(); i++) {
                    if (lights[i].posLight) {
                        Vector3D ld = Vector3D::normalise(lights[i].location - p);
                        double a = Vector3D::dot(n, ld);

                        a *= -1;
                        if (a>0) {
                            diffusePointColor.red += lights[i].diffuseLight.red*a;
                            diffusePointColor.green += lights[i].diffuseLight.green*a;
                            diffusePointColor.blue += lights[i].diffuseLight.blue*a;

                            if (lights[i].specularLighting) {
                                Vector3D r = 2*a*n+ld;
                                r.normalise();
                                Vector3D pn = Vector3D::vector(p);
                                pn.normalise();
                                double b = Vector3D::dot(pn, r);
                                if (b>0) {
                                    specularPointColor.red += lights[i].specularLight.red*pow(b, reflectionCoeff);
                                    specularPointColor.green += lights[i].specularLight.green*pow(b, reflectionCoeff);
                                    specularPointColor.blue += lights[i].specularLight.blue*pow(b, reflectionCoeff);
                                }
                            }
                        }
                    } else if (lights[i].specularLighting) {
                        Vector3D ld = lights[i].ldVector;
                        ld.normalise();
                        double a = Vector3D::dot(n, ld);

                        if (a>0) {
                            Vector3D r = 2*a*n-ld;
                            r.normalise();
                            Vector3D pn = Vector3D::vector(p);
                            pn.normalise();
                            double b = Vector3D::dot(pn, r);
                            if (b>0) {
                                specularPointColor.red += lights[i].specularLight.red*pow(b, reflectionCoeff);
                                specularPointColor.green += lights[i].specularLight.green*pow(b, reflectionCoeff);
                                specularPointColor.blue += lights[i].specularLight.blue*pow(b, reflectionCoeff);
                            }
                        }
                    }
                }
                img::Color color((uint8_t) convertToColor(basecolor.red + diffuseReflection.red * diffusePointColor.red + specularReflection.red * specularPointColor.red),
                                 (uint8_t) convertToColor(basecolor.green + diffuseReflection.green * diffusePointColor.green + specularReflection.green * specularPointColor.green),
                                 (uint8_t) convertToColor(basecolor.blue + diffuseReflection.blue * diffusePointColor.blue + specularReflection.blue * specularPointColor.blue));
                image((unsigned) xi, (unsigned) yi) = color;
            }
        }
    }
}

img::EasyImage
drawFigures3D_zBuff_triangle(Figures3D &figures, Lights3D &lights, const int size, const Color &backgroundColor) {
//    std::chrono::high_resolution_clock::time_point begin = std::chrono::high_resolution_clock::now();

    double xmin = 0, xmax = 0, ymin = 0, ymax = 0;
    double maxRange, Imagex, Imagey, d, DCx, DCy, dx, dy;

    // finding min and max
    for (Figure &f : figures.Figures) {
        for (Vector3D &p : f.points) {
            Point2D pt = doPProjection(p, 1.0);
            if (pt.x < xmin) {
                xmin = pt.x;
            }
            if (pt.x > xmax) {
                xmax = pt.x;
            }
            if (pt.y < ymin) {
                ymin = pt.y;
            }
            if (pt.y > ymax) {
                ymax = pt.y;
            }
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
    for (Figure &f : figures.Figures) {
        f.triangulate();
        for (Face &face: f.faces) {
            draw_zbuff_triangle(zBuffer, image, f.points[face.index[0]], f.points[face.index[1]],
                                f.points[face.index[2]], d, dx, dy, f.ambientReflection, f.diffuseReflection,
                                f.specularReflection, f.reflectionCoefficient, lights);
        }
    }


//    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
//    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
//    std::cout << "ZBuffer drawing took " << duration << "µs to complete." << std::endl;
    return image;
};

img::EasyImage drawLines3D_zbuff_triangle(const ini::Configuration &configuration) {
    Figures3D ImageFigures;
    Lights3D lights;
    const unsigned int size = (unsigned) configuration["General"]["size"].as_int_or_die();
    const std::vector<double> bg = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
    Color bgcolor;
    bgcolor.red = bg[0];
    bgcolor.green = bg[1];
    bgcolor.blue = bg[2];

    generateImage(ImageFigures, configuration);
    Light light;
    light.ambientLight = {1, 1, 1};
    lights.push_back(light);

    if (ImageFigures.Figures.empty()) {
        return img::EasyImage();
    }

    return drawFigures3D_zBuff_triangle(ImageFigures, lights, size, bgcolor);
};