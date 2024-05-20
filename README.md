# RayTracing with MPI

This repo contains an implementation in `C++` using `MPI`  of a parallelized ray tracing algorithm, including implementation of both fixed and adaptive sampling schemes.

## How to run

To run the code compile one of the files `serial.cpp` (for sequential code) or `parallel.cpp` (for the parallelized version) and then run it. If the file is run with a command line argument it generates an ouput image stored into a file named according to the argument, otherwise it performs some benchmarking and outputs the results into a .txt file.

Example for how to run sequential vesrion:
```bash
g++ serial.cpp -o my_serial.out
./my_serial.out
```
This benchmarks the fixed and adaptive rendering versions and stores output into a .txt file.
```bash
./my_serial my_image_name
```
creates two binary files "my_image_name_adaptive.bin" and "my_image_name_fixed.bin". These can be visualized using the vis.py file.

Example for how to run parallelized version:
```bash
mpic++ parallel.cpp -o my_parallel.out
mpirun -np 8 ./my_parallel.out
```
or
```bash
mpirun -np 8 ./my_parallel.out my_image_name
```
This does the same as the sequential code (only faster).