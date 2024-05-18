#include <mpi.h>
#include <fstream>
#include "rendering.h"
#include <chrono>
#include <deque>

#define VERBOSE 1

#define WORK 1
#define DIE 2

using namespace std::chrono;

using namespace std;


vector<double> master(int num_processors, int num_ranges,const Camera& camera,const vector<Sphere>& scene, double tolerance, int max_samples) {
    
    int num_pixels = camera.H*camera.W;
    assert(num_ranges <= num_pixels); // Algorithm only works with fewer processors than pixels
    assert(num_ranges >= num_processors);

    int range_size = num_pixels/num_ranges;
    int num_rendered = 0;
    vector<double> img(num_pixels*4, 0.0);
    int active_processors = 1;

    for(int processor = 1; processor < num_processors; processor++) { // Seed work
        if(num_rendered < num_pixels){
            int work[2] = {num_rendered, min(num_rendered+range_size, num_pixels)};
            num_rendered = work[1];
            MPI_Send(&work, 2, MPI_INT, processor, WORK, MPI_COMM_WORLD);
            active_processors++;
        } else
            MPI_Send(0, 0, MPI_INT, processor, DIE, MPI_COMM_WORLD);
    }
    vector<double> result(range_size*4 + 2, 0); // Recieve and send new pieces of work while there is still work to be done
    MPI_Status status;
    MPI_Request request;
    int flag;
    while  (num_rendered < num_pixels) {
        
        MPI_Irecv(result.data(), result.size(), MPI_DOUBLE, MPI_ANY_SOURCE, WORK, MPI_COMM_WORLD, &request); // Recieve a completed piece of work async
        MPI_Test(&request, &flag, &status);
        while(!flag){ //Waiting for completed range from workers
            if(num_rendered<num_pixels){
                auto rendered = camera.render_range_adaptive(scene, num_rendered, num_rendered+1, tolerance, max_samples);
                for(int i = 0;i<rendered.size();i++) img[num_rendered*4 + i] = rendered[i];
                num_rendered++;
            }
            MPI_Test(&request, &flag, &status);
        }

        int range_start = round(result[0]), range_end = round(result[1]);
        for(int i = 0;i<(range_end-range_start)*4;i++)
            img[i + range_start*4] = result[2 + i];

        if(num_rendered < num_pixels){ // Send a new piece of work
            int work[2] = {num_rendered, min(num_rendered+range_size, num_pixels)};
            num_rendered = work[1];
            MPI_Send(&work, 2, MPI_INT, status.MPI_SOURCE, WORK, MPI_COMM_WORLD);
        } else{
            MPI_Send(0, 0, MPI_INT, status.MPI_SOURCE, DIE, MPI_COMM_WORLD);
            active_processors--;
        }
    }

    // Receive the last pieces of work
    while(active_processors > 1){
        MPI_Recv(result.data(), result.size(), MPI_DOUBLE, MPI_ANY_SOURCE, WORK, MPI_COMM_WORLD, &status);
        int range_start = round(result[0]), range_end = round(result[1]);
        
        for(int i = 0; i < (range_end-range_start)*4; i++)
            img[i + range_start*4] = result[2 + i];
        MPI_Send(0, 0, MPI_INT, status.MPI_SOURCE, DIE, MPI_COMM_WORLD);
        active_processors--;
    }

    return img;
}

void slave(const Camera& camera,const vector<Sphere>& scene, int rank, double tolerance, int max_samples) {
    MPI_Status status;

    while (true) {
        int work[2];
        MPI_Recv(&work, 2, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if (status.MPI_TAG == DIE)
            return;
#if VERBOSE >= 2
        std::cout << "Slave Processor " << rank << " Recieved work " << work[0] << "-" << work[1] << endl;
#endif

        vector<double> result;
        result.push_back((double)work[0]);
        result.push_back((double)work[1]);
        auto rendered = camera.render_range_adaptive(scene, work[0], work[1], tolerance, max_samples);
        result.insert(result.end(), rendered.begin(), rendered.end());
        MPI_Send(result.data(), result.size(), MPI_DOUBLE, 0, WORK, MPI_COMM_WORLD);
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

void save_img(vector<double>& img, vector<int> shape, string fname){
    ofstream output;
    output.open(fname+ ".bin");
    output.write(reinterpret_cast<char*>(shape.data()), sizeof(int)*3);
    output.write(reinterpret_cast<char*>(img.data()), img.size()*sizeof(double));
}

void generate_image_adaptive(const Camera& camera, const vector<Sphere>& scene, int num_ranges, double tolerance, int max_samples, string fname = "image_adaptive") {
    int P, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    if (rank == 0) {
        auto start = high_resolution_clock::now();
        auto img = master(P, num_ranges, camera, scene, tolerance, max_samples);
        auto elapsed = chrono::duration_cast<chrono::nanoseconds>(high_resolution_clock::now()-start).count() / 1e9;
#if VERBOSE >= 1
        std::cout << "Done!" << "Time: " << elapsed << "s" << " total pixels: " << camera.H*camera.W << endl;
#endif
        save_img(img, {camera.H, camera.W, 4}, fname);
    } else{
        slave(camera, scene, rank, tolerance, max_samples);
    }
}
void generate_image_fixed(const Camera& camera, const vector<Sphere>& scene, int samples_per_processor, string fname = "image_fixed"){
    int P, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    auto start = high_resolution_clock::now();
    auto img = camera.render(scene, samples_per_processor);
    vector<double> merged;
    if(rank == 0)
        merged.resize(img.size());
    
    MPI_Reduce(img.data(), merged.data(), img.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if(rank == 0){
        for(int i = 0;i<merged.size();i++) merged[i] /= P;
        auto elapsed = duration_cast<nanoseconds>(high_resolution_clock::now() - start).count() / 1e9;
        save_img(img, {camera.H, camera.W, 3}, fname);
#if VERBOSE >= 1
        std::cout << "Image generated. Time: " << elapsed << "s" << endl;
#endif
    }
}
double timing_seconds_adaptive(const Camera& camera, const vector<Sphere>& scene, double tolerance, int max_samples, int num_ranges, int num_runs=10) {
    int P, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double sum_time = 0;

    for(int run = 1; run <= num_runs; run++) {
        
        int num_pixels = camera.H * camera.W;

        MPI_Barrier(MPI_COMM_WORLD);

        if (rank == 0) {
            auto start = high_resolution_clock::now();
            auto img = master(P, num_ranges, camera, scene, tolerance, max_samples);
            auto end = high_resolution_clock::now();
            sum_time += (double)duration_cast<microseconds>(end - start).count() / 1000000.0;
#if VERBOSE >= 1
            std::cout << "Average time: " << sum_time / (double)run << "s" << endl;
#endif
        } else{
            slave(camera, scene, rank, tolerance, max_samples);
        }
    }
    return sum_time / (double)num_runs;
}
double timing_seconds_fixed(const Camera& camera, const vector<Sphere>& scene, int samples_per_processor, int num_runs = 10){
    int P, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    double sum_time = 0.0;

    for(int run = 1;run<=num_runs;run++){
        MPI_Barrier(MPI_COMM_WORLD);
        auto start = high_resolution_clock::now();
        auto img = camera.render(scene, samples_per_processor);
        vector<double> merged;
        if(rank == 0)
            merged.resize(img.size());
        
        MPI_Reduce(img.data(), merged.data(), img.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if(rank == 0){
            for(int i = 0;i<merged.size();i++) merged[i] /= P;
            sum_time += duration_cast<nanoseconds>(high_resolution_clock::now() - start).count() / 1e9;
#if VERBOSE >= 1
            std::cout << "Average time: " << sum_time / (double)run << "s" << endl;
#endif
        }
    }
    return sum_time / (double)num_runs;
}

int main(int argc, char **argv) {

    int num_runs = 5;
    int samples_per_pixel = 40;
    double tolerance = 0.01;
    int n_small_balls = 100;
    Camera camera;
    camera.move({1,0,0}, {1,1,0});

    int num_ranges; 
    int P, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int samples_per_processor = (int) ceil(samples_per_pixel/ (double) P);
    int max_samples = samples_per_processor*P;
    auto scene = generate_common_scene(n_small_balls);
    vector<double> times;
    vector<int> ranges_vec;

    for(num_ranges = P;num_ranges<camera.H*camera.W; num_ranges *= 10){
        times.push_back(timing_seconds_adaptive(camera, scene, tolerance, max_samples, num_ranges, num_runs));
        ranges_vec.push_back(num_ranges);
    }
    generate_image_adaptive(camera, scene, 800, tolerance, max_samples, "adaptive_"+to_string(P)+"_"+to_string(samples_per_processor)+"_"+to_string(tolerance));
    double time_fixed = timing_seconds_fixed(camera, scene, samples_per_processor, num_runs);
    generate_image_fixed(camera, scene, samples_per_processor, "fixed_"+to_string(P)+"_"+to_string(samples_per_processor));

    if(rank == 0){
        ofstream output;
        output.open("data"+to_string(P)+"_"+to_string(samples_per_processor)+"_"+to_string(tolerance)+".txt");
        output << "Num runs used:"<<to_string(num_runs)<<"\n";
        output<<"Adaptive avg times:\n";
        for(int i = 0;i<ranges_vec.size();i++){
            output << "\t" << "num_ranges="<<ranges_vec[i] << ", " << "avg_time="<<times[i]<<"s\n";
        }
        output << "Fixed avg time:" << time_fixed << "s";
    }

    MPI_Finalize();
}