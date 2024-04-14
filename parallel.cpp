#include <mpi.h>
#include <fstream>
#include "rendering.h"

using namespace std;

int main(int argc, char **argv) {

    int samples_per_pixel = 5;
    if(argc > 1)
        samples_per_pixel = stoi(argv[1]);
    MPI_Init(&argc, &argv);

    int P, rank;
    MPI_Status* status;
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    Camera camera;
    vector<Sphere> scene = {
        Sphere(0.0, -10000.0, 0.0, 9999.0), //Earth
        Sphere({-1.1, 0.0, 2.0}, 1.0, {0.5, 0.5, 0.5}, make_shared<BlankMaterial>(0.5)),
        Sphere({1.1, 0.0, 2.0}, 1.0, {0.5, 0.5, 0.5}, make_shared<BlankMaterial>(0.5)),
        Sphere({0.0, -0.7, 1.3}, 0.3, {0.5, 0.5, 0.5}, make_shared<BlankMaterial>(0.5)),
    };

    auto img = camera.render(scene, samples_per_pixel);
    vector<double> merged;
    if(rank == 0)
        merged.resize(img.size());
    
    MPI_Reduce(img.data(), merged.data(), img.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if(rank == 0){
        for(int i = 0;i<merged.size();i++) merged[i] /= P;
        ofstream output;
        output.open("img.bin");
        int shape[3] = {camera.H, camera.W, 3};
        output.write(reinterpret_cast<char*>(shape), sizeof(int)*3);
        output.write(reinterpret_cast<char*>(merged.data()), merged.size()*sizeof(double));
        std::cout << "Done! ( in total " << P*samples_per_pixel << " rays cast from " << P << " processors)" <<endl;
    }

    MPI_Finalize();
}