#include <mpi.h>
#include <fstream>
#include "rendering.h"
#include <chrono>
#include <deque>


#define WORK 1
#define DIE 2



using namespace std::chrono;

using namespace std;








vector<double> master(int num_processors, int num_ranges, int num_pixels) {
    
    // Algorithm only works with fewer processors than pixels
    assert(num_ranges <= num_pixels);
    assert(num_ranges >= num_processors);

    
    
    
    
    deque<pair<int, int>> ranges;

    int range_size = 1 + (num_pixels - 1) / (num_ranges);

    // cout << "Seeding work..." << endl;
    for(int work = 0; work < num_ranges; work++) {
        ranges.emplace_back(work*range_size, min((work+1)*range_size, num_pixels));
    }
    vector<double> img(num_pixels*3, 0.0);
    

    // Seed work
    for(int processor = 1; processor < min(num_processors, num_ranges+1); processor++) {
        auto [start, end] = ranges.front();
        ranges.pop_front();
        int work[2] = {start, end};
        // cout << "Sending work " << work[0] << "-" << work[1] << endl;
        MPI_Send(&work, 2, MPI_INT, processor, WORK, MPI_COMM_WORLD);
        
    }
    // cout << "Handling new work..." << endl;
    // Recieve and send new pieces of work while there is still work to be done
    vector<double> result(range_size*3 + 2, 0);
    MPI_Status status;
    while  (!ranges.empty()) {
        
        // Recieve a completed piece of work
        MPI_Recv(static_cast<void*>(result.data()), range_size*3+2, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        int range_start = round(result[0]);
        int range_end = round(result[1]);
        int true_range_size = range_end - range_start;

        for(int i = 0; i < true_range_size*3; i++) {
            img[i + range_start*3] = result[2 + i];
        }

        // Send a new piece of work
        auto [start, end] = ranges.front();
        ranges.pop_front();
        int work[2] = {start, end};
        MPI_Send(&work, 2, MPI_INT, status.MPI_SOURCE, WORK, MPI_COMM_WORLD);

    }

    // cout << "Recieving Last Process Work..." << endl;
    // Receive the last pieces of work
    for(int processor = 1; processor < min(num_processors, num_ranges+1); processor++) {
        // Recieve a completed piece of work
        // cout << "Processor " << processor << endl;
        MPI_Recv(static_cast<void*>(result.data()), range_size*3+2, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        // cout << "Recieved.." << endl;
        int range_start = round(result[0]);
        int range_end = round(result[1]);
        int true_range_size = range_end - range_start;
        
        for(int i = 0; i < true_range_size*3; i++) {
            
            img[i + range_start*3] = result[2 + i];
        }
        MPI_Send(0, 0, MPI_INT, status.MPI_SOURCE, DIE, MPI_COMM_WORLD);
    }


    return img;
}

void slave(Camera& camera, vector<Sphere>& scene, int rank, double tolerance) {
    MPI_Status status;

    // cout << "Initiated Slave Processor " << rank << endl;
    while (true) {
        int work[2];
        MPI_Recv(&work, 2, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        cout << "Slave Processor " << rank << " Recieved work " << work[0] << "-" << work[1] << endl;
        if (status.MPI_TAG == DIE) {
            return;
        }
        int msg_len = 3*(work[1] - work[0]) + 2;

        vector<double> result;
        result.push_back((double)work[0]);
        result.push_back((double)work[1]);
        auto rendered = camera.render_range_adaptive(scene, work[0], work[1], tolerance);
        result.insert(result.end(), rendered.begin(), rendered.end());
        MPI_Send(static_cast<void*>(result.data()), msg_len, MPI_DOUBLE, 0, WORK, MPI_COMM_WORLD);

    }
}

vector<Sphere> generate_common_scene(int n_small_balls = 100, double small_r = 0.1) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    vector<Sphere> scene;
    SceneParams scene_params;
    
    // Distribute the scene parameters to all processes
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
    return scene;
}



double timing_seconds(double tolerance, int groups_per_processor, int num_runs=10) {


    int P, rank;
    MPI_Status* status;
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    Camera camera;
    camera.move({1,0,0}, {1,1,0});
    double sum_time = 0;

    for(int run = 1; run <= num_runs; run++) {
        vector<Sphere> scene = generate_common_scene();
        
        int num_pixels = camera.H * camera.W;

        MPI_Barrier(MPI_COMM_WORLD);

        if (rank == 0) {
            auto start = high_resolution_clock::now();
            auto img = master(P, P*groups_per_processor, num_pixels);
            auto end = high_resolution_clock::now();
            sum_time += (double)duration_cast<microseconds>(end - start).count() / 1000000.0;
            cout << "Average time: " << sum_time / (double)run << "s" << endl;
        } else{
            slave(camera, scene, rank, tolerance);
        }
    }
    return sum_time / (double)num_runs;

}


void generate_image(int groups_per_processor, double tolerance, string fname = "image") {
    int P, rank;
    MPI_Status* status;
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    Camera camera;
    camera.move({1,0,0}, {1,1,0});

    vector<Sphere> scene = generate_common_scene();
    
    int num_pixels = camera.H * camera.W;
    if (rank == 0) {
        auto img = master(P, P*groups_per_processor, num_pixels);
        std::cout << "Done!";
        
        ofstream output;
        output.open(fname+ ".bin");
        int shape[3] = {camera.H, camera.W, 3};
        output.write(reinterpret_cast<char*>(shape), sizeof(int)*3);
        output.write(reinterpret_cast<char*>(img.data()), img.size()*sizeof(double));

    } else{
        slave(camera, scene, rank, tolerance);
    }
}


int main(int argc, char **argv) {


    int groups_per_processor = 4;
    double tolerance = 0.05;
    if(argc > 1)
        tolerance = stod(argv[1]);
    if(argc > 2)
        groups_per_processor = stoi(argv[2]);
    

    int n_small_balls = 100;
    double small_r = 0.1;

    MPI_Init(&argc, &argv);

    generate_image(groups_per_processor, tolerance);


    // auto img = camera.render(scene, samples_per_pixel);
    // vector<double> merged;
    // if(rank == 0)
    //     merged.resize(img.size());
    
    // MPI_Reduce(img.data(), merged.data(), img.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    // if(rank == 0){
    //     for(int i = 0;i<merged.size();i++) merged[i] /= P;
    //     ofstream output;
    //     output.open("img.bin");
    //     int shape[3] = {camera.H, camera.W, 3};
    //     output.write(reinterpret_cast<char*>(shape), sizeof(int)*3);
    //     output.write(reinterpret_cast<char*>(merged.data()), merged.size()*sizeof(double));
    //     std::cout << "Done! In total " << P*samples_per_pixel << " rays per pixel cast from " << P << " processors." <<endl;
    // }

    MPI_Finalize();
}