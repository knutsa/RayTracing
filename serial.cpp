#include <iostream>
#include <fstream>
#include  "rendering.h"

using namespace std;

int main(int argc, char** argv) {
    Camera camera;

    vector<Sphere> scene = {
        Sphere(0.0, -10000.0, 0.0, 9999.0), //Earth
        Sphere({-1.1, 0.0, 2.0}, 1.0, {0.5, 0.5, 0.5}, make_shared<BlankMaterial>(0.5)),
        Sphere({1.1, 0.0, 2.0}, 1.0, {0.5, 0.5, 0.5}, make_shared<BlankMaterial>(0.5)),
        Sphere({0.0, -0.7, 1.3}, 0.3, {0.5, 0.5, 0.5}, make_shared<BlankMaterial>(0.5)),
    };

    auto img = camera.render(scene);

    ofstream output;
    output.open("img.bin");
    int shape[3] = {camera.H, camera.W, 3};
    output.write(reinterpret_cast<char*>(shape), sizeof(int)*3);
    output.write(reinterpret_cast<char*>(img.data()), sizeof(double)*img.size());
    output.close();

    std::cout << "Done!" << endl;
}