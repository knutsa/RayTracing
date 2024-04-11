#include <iostream>
#include <vector>
#include <random>
#include <fstream>

using namespace std;


const int MAX_DEPTH = 5;
const int H = 900, W = 1600;
const double fx = 100.0, fy = 100.0, cx = 800, cy = 450;

default_random_engine RNG{1337};
uniform_real_distribution U;

class Vec3D {
    public:
        double arr[3];
        Vec3D(double x, double y, double z) : arr{x,y,z} {}
        Vec3D() : arr{0,0,0} {}

        double operator[](int i) const {return arr[i];}
        Vec3D operator+(const Vec3D& other) const {return Vec3D(arr[0]+other.arr[0], arr[1]+other.arr[1], arr[2]+other.arr[2]);}
        Vec3D operator+=(const Vec3D& other){
            arr[0] += other.arr[0];
            arr[1] += other.arr[1];
            arr[2] += other.arr[2];
            return *this;
        }
        Vec3D operator+(double a) const {return Vec3D(arr[0]+a, arr[1]+a, arr[2]+a);}
        Vec3D operator-(double a) const {return Vec3D(arr[0]-a, arr[1]-a, arr[2]-a);}
        Vec3D operator-(const Vec3D& other) const { return (*this)+(other*-1); }
        Vec3D operator-=(const Vec3D& other) {return (*this)+=(other*-1); }

        Vec3D operator*(double a) const { return Vec3D(arr[0]*a, arr[1]*a, arr[2]*a); }
        Vec3D operator*=(double a){
            arr[0] *= a;
            arr[1] *= a;
            arr[2] *= a;
            return *this;
        }
        Vec3D operator/(double a) const { return (*this)*(1/a); }
        Vec3D operator/=(double a){ return (*this)*=(1/a);}

        double norm2() const {return arr[0]*arr[0]+arr[1]*arr[1]+arr[2]*arr[2];}
        double norm() const {return sqrt(this->norm2());}
};
double dot(const Vec3D& x, const Vec3D& y){ return x[0]*y[0]+x[1]*y[1]+x[2]*y[2]; }
void print_vec3d(const Vec3D& x, string caption){
    std::cout << caption << "(" <<x[0] << ","<<x[1]<<","<<x[2]<<")"<<endl;
}
const Vec3D SKY_COLOR(1.0, 1.0, 1.0), BACKGROUND_COLOR(0, 0.0, 0.5);

struct Sphere {
    Vec3D center;
    double radius;
    //Maybe add optical properties ... absorbtion, diffusive properties ...
    Sphere(double cx, double cy, double cz, double R) : center(cx, cy, cz), radius(R) {}

    double first_intersection(const Vec3D& source, const Vec3D& dir, double t_min, double t_max){
        //pq-formula
        double p_half = dot(dir, source-center)/dir.norm2(),
            q = ((source-center).norm2()-radius*radius)/dir.norm2();
        double discriminant = (p_half*p_half - q);
        if(discriminant < 0)
            return INFINITY;
        double sqrt_disc = sqrt(discriminant);
        double root = -p_half - sqrt_disc;
        if(root >= t_min && root < t_max)
            return root;
        root = -p_half + sqrt_disc;
        if(root >= t_min && root < t_max)
            return root;
        
        return INFINITY;
    }
    Vec3D normal_at(const Vec3D& x){
        auto res = x - center;
        res /= res.norm();
        return res;
    }
};

pair<double, double> inds2uv(int i, int j){
    return pair<double, double>{j, H-1-i};
}
Vec3D uv2ray(double u, double v) {
    //u, v horizontal, resp vertical coordinates in image
    return Vec3D((u-cx)/fx, (v-cy)/fy, 1.0);
}
Vec3D color_free_ray(const Vec3D& source,const Vec3D& dir){
    // double p = max(dir[1]/dir.norm(), 0.0); //y-coordinate points up
    double p = 0.5*(dir[1]/dir.norm() + 1.0); //y-coordinate points up
    return SKY_COLOR*p + BACKGROUND_COLOR*(1-p);
}

Vec3D unit_hemisphere_sample(const Vec3D& n){
    double z = 2*U(RNG)-1, phi = U(RNG)*2*M_PI;
    double xy_radius = sqrt(1 - z*z);
    auto res = Vec3D(xy_radius*cos(phi), xy_radius*sin(phi), z);

    if(dot(res,n) < 0) res *= -1;
    return res;
}

Vec3D color_ray(const Vec3D& source,const Vec3D& dir, const vector<Sphere>& scene, int depth = 0, double t_min = 0.0, double t_max = INFINITY){
    if(depth > MAX_DEPTH)
        return BACKGROUND_COLOR;

    double closest_hit = t_max;
    Vec3D intersection, normal;

    for(auto sphere : scene){
        double t = sphere.first_intersection(source, dir, t_min, t_max);
        if(t < closest_hit){
            closest_hit = t;
            intersection = source+dir*t;
            normal = sphere.normal_at(intersection);
        }    
    }

    if(closest_hit < t_max){ //Something was hit!
        // auto new_dir = unit_hemisphere_sample(normal);
        // return color_ray(intersection, new_dir, scene, depth+1, t_min, t_max)*0.5;
        return (normal+1)/2;
    }

    //Nothing Hit
    return color_free_ray(source, dir);
}


int main(){
    Vec3D zero;
    print_vec3d(zero, "ZERO: ");

    print_vec3d(SKY_COLOR, "SKY COLOR: ");
    print_vec3d(BACKGROUND_COLOR, "BACKGROUND_COLOR: ");

    vector<double> res(H*W*3);

    vector<Sphere> scene = {
        //Sphere(0.0, 0.0, 0.0, 1000.0), //Earth
        Sphere(0.0, 0.0, 5.5, 5.0)
    };

    for(int i = 0;i<H;i++){
        for(int j = 0;j<W;j++){
            auto uv = inds2uv(i, j);
            auto ray = uv2ray(uv.first, uv.second);
            // std::cout << ray[1]/ray.norm() << endl;
            // auto color = color_free_ray(zero, ray);
            auto color = color_ray(zero, ray, scene);
            // std::cout << color[0] << " " << color[1] << " " << color[2] << endl;

            for(int k=0;k<3;k++) res[(i*W+j)*3+k] = color[k];
        }
    }

    ofstream output;
    output.open("img.bin");
    int shape[3] = {H, W, 3};
    output.write(reinterpret_cast<char*>(&shape), sizeof(int)*3);
    output.write(reinterpret_cast<char*>(res.data()), sizeof(double)*H*W*3);

    std::cout << "Done"<<endl;
}
