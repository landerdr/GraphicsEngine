#include "easy_image.h"
#include "ini_configuration.h"
#include "session7.cpp"


img::EasyImage generate_image(const ini::Configuration &configuration) {
    const std::string type = configuration["General"]["type"].as_string_or_die();
    if (type == "IntroColorRectangle") {
        return generate_IntroColorRectangle(configuration);
    }
    if (type == "IntroBlocks") {
        return generate_IntroBlocks(configuration);
    }
    if (type == "IntroLines") {
        return generate_IntroLines(configuration);
    }
    if (type == "2DLSystem") {
        return drawLSystem(configuration);
    }
//    if (type == "Lines") {
//        return Lines(configuration);
//    }
    if (type == "Wireframe") {
        return drawLines3D(configuration);
    }
    if (type == "ZBufferedWireframe") {
        return zbuffdrawLines3D(configuration);
    }
    if (type == "ZBuffering") {
        return drawLines3D_zbuff_triangle(configuration);
    }
    if (type == "LightedZBuffering") {
        return drawLines3D_zbuff_triangle_light(configuration);
    }
    return img::EasyImage();
}

int main(int argc, char const *argv[]) {
    int retVal = 0;
    try {
        for (int i = 1; i < argc; ++i) {
            ini::Configuration conf;
            try {
                std::ifstream fin(argv[i]);
                fin >> conf;
                fin.close();
            }
            catch (ini::ParseException &ex) {
                std::cerr << "Error parsing file: " << argv[i] << ": " << ex.what() << std::endl;
                retVal = 1;
                continue;
            }

            std::chrono::high_resolution_clock::time_point begin = std::chrono::high_resolution_clock::now();
            img::EasyImage image = generate_image(conf);
            std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
            std::cout << "Image generation took " << duration << "µs to complete." << std::endl;

            if (image.get_height() > 0 && image.get_width() > 0) {
                std::string fileName(argv[i]);
                std::string::size_type pos = fileName.rfind('.');
                if (pos == std::string::npos) {
                    //filename does not contain a '.' --> append a '.bmp' suffix
                    fileName += ".bmp";
                } else {
                    fileName = fileName.substr(0, pos) + ".bmp";
                }
                try {
                    std::ofstream f_out(fileName.c_str(), std::ios::trunc | std::ios::out | std::ios::binary);
                    f_out << image;

                }
                catch (std::exception &ex) {
                    std::cerr << "Failed to write image to file: " << ex.what() << std::endl;
                    retVal = 1;
                }
            } else {
                std::cout << "Could not generate image for " << argv[i] << std::endl;
            }
        }
    }
    catch (const std::bad_alloc &exception) {
        //When you run out of memory this exception is thrown. When this happens the return value of the program MUST be '100'.
        //Basically this return value tells our automated test scripts to run your engine on a pc with more memory.
        //(Unless of course you are already consuming the maximum allowed amount of memory)
        //If your engine does NOT adhere to this requirement you risk losing points because then our scripts will
        //mark the test as failed while in reality it just needed a bit more memory
        std::cerr << "Error: insufficient memory" << std::endl;
        retVal = 100;
    }
    return retVal;
}
