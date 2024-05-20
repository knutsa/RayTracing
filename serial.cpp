#include <iostream>
#include <fstream>
#include  "rendering.h"
#include <chrono>

using namespace std::chrono;
using namespace std;

double time_seconds_fixed(const Camera& camera, const vector<Sphere>& scene, int samples_per_processor, int num_runs = 10){
    double sum_time = 0.0;

    for(int run = 1;run<=num_runs;run++){
        auto start = high_resolution_clock::now();
        auto img = camera.render(scene, samples_per_processor);
        sum_time += duration_cast<nanoseconds>(high_resolution_clock::now() - start).count() / 1e9;

        std::cout << "Average time: " << sum_time / (double)run << "s" << endl;
    }
    return sum_time / (double)num_runs;
}
double time_seconds_adaptive(const Camera& camera, const vector<Sphere>& scene, double tolerance, int max_samples, int num_runs = 10){
    double sum_time = 0.0;

    for(int run = 1;run<=num_runs;run++){
        auto start = high_resolution_clock::now();
        auto img = camera.render_range_adaptive(scene, 0, camera.H*camera.W, tolerance, max_samples);
        sum_time += duration_cast<nanoseconds>(high_resolution_clock::now() - start).count() / 1e9;

        std::cout << "Average time: " << sum_time / (double)run << "s" << endl;
    }
    return sum_time / (double)num_runs;
}
void save_img(vector<double>& img, vector<int> shape, string fname){
    ofstream output;
    output.open(fname+ ".bin");
    output.write(reinterpret_cast<char*>(shape.data()), sizeof(int)*3);
    output.write(reinterpret_cast<char*>(img.data()), img.size()*sizeof(double));
}

int main(int argc, char** argv) {
    int samples_per_pixel = 40;
    double tolerance = 0.01;
    int n_small_balls = 100;
    Camera camera;
    camera.move({1,0,0}, {1,1,0});
    auto scene = generate_random_scene(n_small_balls);

    if(argc > 1){ //Generate Output image
        string image_name = argv[1];
        auto img = camera.render_range_adaptive(scene, 0, camera.H*camera.W, tolerance, samples_per_pixel);
        save_img(img, {camera.H, camera.W, 4}, image_name+"adaptive_P=1");
        img = camera.render(scene, samples_per_pixel);
        save_img(img, {camera.H, camera.W, 3},image_name+"fixed_P=1");
    } else { //Benchmarks
        int num_runs = 5;
        auto time_adaptive = time_seconds_adaptive(camera, scene, tolerance, samples_per_pixel, num_runs);
        auto time_fixed = time_seconds_fixed(camera, scene, samples_per_pixel, num_runs);
        ofstream output;
        output.open("data_serial_"+to_string(samples_per_pixel)+"_"+to_string(tolerance)+".txt");
        output << "Num runs used:"<<to_string(num_runs)<<"\n";
        output<<"Adaptive avg time: " << time_adaptive << "s\n";
        output << "Fixed avg time: " << time_fixed << "s";
    }


}