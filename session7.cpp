//
// Created by Lander on 5/7/18.
//

#include "session5.cpp"

void generateLights(Lights3D &lights, const ini::Configuration &configuration) {
    const unsigned int amount = (unsigned) configuration["General"]["nrLights"].as_int_or_die();
    const std::vector<double> eyep = configuration["General"]["eye"].as_double_tuple_or_die();
    Vector3D eye;
    eye.x = eyep[0];
    eye.y = eyep[1];
    eye.z = eyep[2];

    Matrix eyetrans = generateEyePerspectiveMatrix(eye);

    for (unsigned int i = 0; i < amount; ++i) {
        bool infLight;
        Light light;
        if (configuration["Light" + std::to_string(i)]["infinity"].as_bool_if_exists(infLight)) {
            if (infLight) {
                light.infLight = true;
                light.posLight = false;
                const std::vector<double> orientation = configuration["Light" + std::to_string(i)]["direction"].as_double_tuple_or_die();
                light.ldVector = Vector3D::vector(orientation[0], orientation[1], orientation[2]);
                light.ldVector *= eyetrans;
            } else {
                light.posLight = true;
                light.infLight = false;
                const std::vector<double> orientation = configuration["Light" + std::to_string(i)]["location"].as_double_tuple_or_die();
                light.location = Vector3D::point(orientation[0], orientation[1], orientation[2]);
                light.location *= eyetrans;
            }
        }
        std::vector<double> diffuseLight;
        const std::vector<double> ambientLight = configuration["Light" + std::to_string(i)]["ambientLight"].as_double_tuple_or_die();
        light.ambientLight.red = ambientLight[0];
        light.ambientLight.green = ambientLight[1];
        light.ambientLight.blue = ambientLight[2];

        if (configuration["Light" + std::to_string(i)]["diffuseLight"].as_double_tuple_if_exists(diffuseLight)) {
            light.diffuseLight.red = diffuseLight[0];
            light.diffuseLight.green = diffuseLight[1];
            light.diffuseLight.blue = diffuseLight[2];
        } else {
            light.diffuseLight.red = 0;
            light.diffuseLight.green = 0;
            light.diffuseLight.blue = 0;
        }
        std::vector<double> specularLight;
        if (configuration["Light" + std::to_string(i)]["specularLight"].as_double_tuple_if_exists(specularLight)) {
            light.specularLighting = true;
            light.specularLight.red = specularLight[0];
            light.specularLight.green = specularLight[1];
            light.specularLight.blue = specularLight[2];
        } else {
            light.specularLighting = false;
            light.specularLight.red = 0;
            light.specularLight.green = 0;
            light.specularLight.blue = 0;
        }

        lights.push_back(light);
    }
}

img::EasyImage drawLines3D_zbuff_triangle_light(const ini::Configuration &configuration) {
    Figures3D ImageFigures;
    Lights3D lights;
    const unsigned int size = (unsigned) configuration["General"]["size"].as_int_or_die();
    const std::vector<double> bg = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
    Color bgcolor;
    bgcolor.red = bg[0];
    bgcolor.green = bg[1];
    bgcolor.blue = bg[2];

    generateImage(ImageFigures, configuration);
//    std::cout << "figures" << std::endl;
    generateLights(lights, configuration);
//    std::cout << "lights" << std::endl;

    if (ImageFigures.Figures.empty()) {
        return img::EasyImage();
    }

    return drawFigures3D_zBuff_triangle(ImageFigures, lights, size, bgcolor);
};