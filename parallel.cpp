#include <mpi.h>
#include <fstream>
#include "rendering.h"

using namespace std;

int main(int argc, char **argv) {


    int samples_per_pixel = 5;
    if(argc > 1)
        samples_per_pixel = stoi(argv[1]);

    int n_small_balls = 100;
    double small_r = 0.1;

    MPI_Init(&argc, &argv);

    int P, rank;
    MPI_Status* status;
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    Camera camera;
    camera.move({1,0,0}, {1,1,0});

    vector<Sphere> scene;
    SceneParams scene_params;
    if(rank == 0){
        scene = generate_random_scene(n_small_balls, small_r);
        scene_params = save_scene_params(scene);
    } else {
        scene_params.second.resize(n_small_balls+3);
        scene_params.first.resize((n_small_balls+3)*max_material_param_size());
    }
    MPI_Bcast(scene_params.first.data(), scene_params.first.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(scene_params.second.data(), scene_params.second.size(), MPI_INT, 0, MPI_COMM_WORLD);
    if(rank != 0)
        scene = parse_scene_params(scene_params);    

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
        std::cout << "Done! In total " << P*samples_per_pixel << " rays per pixel cast from " << P << " processors." <<endl;
    }

    MPI_Finalize();
}