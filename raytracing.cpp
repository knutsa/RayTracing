#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <memory>
#include <cassert>

using namespace std;


const int MAX_DEPTH = 5, SAMPLES_PER_PIXEL = 10;
const int H = 900, W = 1600;
const double fx = 500.0, fy = 500.0, cx = 800, cy = 450;

random_device rd;
default_random_engine RNG{rd()};
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
        Vec3D operator*(const Vec3D& other) const { return Vec3D(arr[0]*other[0], arr[1]*other[1], arr[2]*other[2]); }
        Vec3D operator/(double a) const { return (*this)*(1/a); }
        Vec3D operator/=(double a){ return (*this)*=(1/a);}

        double norm2() const {return arr[0]*arr[0]+arr[1]*arr[1]+arr[2]*arr[2];}
        double norm() const {return sqrt(this->norm2());}
};
double dot(const Vec3D& x, const Vec3D& y){ return x[0]*y[0]+x[1]*y[1]+x[2]*y[2]; }
void print_vec3d(const Vec3D& x, string caption){
    std::cout << caption << "(" <<x[0] << ","<<x[1]<<","<<x[2]<<")"<<endl;
}
const Vec3D SKY_COLOR(1.0, 1.0, 1.0), BACKGROUND_COLOR(0, 0.0, 0.0);

Vec3D unit_sphere_sample(){
    double z = 2*U(RNG)-1, phi = U(RNG)*2*M_PI;
    double xy_radius = sqrt(1 - z*z);
    return Vec3D(xy_radius*cos(phi), xy_radius*sin(phi), z);
}
Vec3D unit_hemisphere_sample(const Vec3D& n){
    auto res = unit_sphere_sample();

    if(dot(res,n) < 0) res *= -1;
    return res;
}
Vec3D lambertian_sample(const Vec3D& n){
    auto res = n + unit_sphere_sample();
    res /= res.norm();
    return res;
}

struct MaterialBase { //Abstract base struct for optic material
    virtual Vec3D scatter(const Vec3D& normal, const Vec3D& ray_dir) = 0;
    MaterialBase() {}
};
struct UniDiffusionMaterial : MaterialBase {
    Vec3D scatter(const Vec3D& normal, const Vec3D& ray_dir) override {
        return unit_hemisphere_sample(normal);
    }
    UniDiffusionMaterial() {}
};
struct LambertianMaterial : MaterialBase {
    Vec3D scatter(const Vec3D& normal, const Vec3D& ray_dir) override {
        return lambertian_sample(normal);
    }
    LambertianMaterial() {}
};
struct BlankMaterial : MaterialBase {
    double fuzziness;

    BlankMaterial() : fuzziness(0.0) {}
    BlankMaterial(double fuzziness) : fuzziness(fuzziness) {}

    Vec3D scatter(const Vec3D& normal, const Vec3D& ray_dir) override {

        auto res = ray_dir - normal*(2*dot(ray_dir, normal));
        if(fuzziness == 0) return res;
        res += unit_sphere_sample()*fuzziness;
        if(dot(res, normal) < 0) return Vec3D();
        return res;
    }
};
//...Refractive Material...

struct Sphere {
    Vec3D center;
    double radius;
    Vec3D albedo; //rgb albedo
    shared_ptr<MaterialBase> material; //Material properties for diffusion/reflection/refraction

    Sphere(double cx, double cy, double cz, double R) :
        center(cx, cy, cz), radius(R), albedo(0.5, 0.5, 0.5), material(make_shared<LambertianMaterial>()) { }
    Sphere(vector<double> center, double R, vector<double> albedo, shared_ptr<MaterialBase> mat) :
        center(center[0], center[1], center[2]), radius(R), albedo(albedo[0], albedo[1], albedo[2]), material(mat) {}


    double first_intersection(const Vec3D& source, const Vec3D& dir, double t_min, double t_max) const {
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
    Vec3D normal_at(const Vec3D& x) const {
        auto res = x - center;
        res /= res.norm();
        return res;
    }
};

//uv coordinates to lower left corner of pixel
pair<double, double> inds2uv(int i, int j){
    return pair<double, double>{j, H-1-i};
}
//u, v horizontal, resp vertical coordinates in image
Vec3D uv2ray(double u, double v) {
    return Vec3D((u-cx)/fx, (v-cy)/fy, 1.0);
}
Vec3D color_free_ray(const Vec3D& source,const Vec3D& dir){
    // double p = max(dir[1]/dir.norm(), 0.0); //y-coordinate points up
    double p = 0.5*(dir[1]/dir.norm() + 1.0); //y-coordinate points up
    return SKY_COLOR*p + BACKGROUND_COLOR*(1-p);
}

Vec3D color_ray(const Vec3D& source,const Vec3D& dir, const vector<Sphere>& scene, int depth = 0, double t_min = 0.0, double t_max = INFINITY){
    if(depth > MAX_DEPTH)
        return BACKGROUND_COLOR;

    double closest_hit = t_max;
    const Sphere * hit;

    for(int i = 0;i<scene.size();i++){
        double t = scene[i].first_intersection(source, dir, t_min, t_max);
        if(t < closest_hit){
            closest_hit = t;
            hit = &scene[i];
        }    
    }

    if(closest_hit < t_max){ //Something was hit!
        auto intersection = source+dir*closest_hit;
        auto normal = hit->normal_at(intersection);
        auto new_dir = hit->material->scatter(normal, dir);
        if(new_dir.norm2() == 0) return BACKGROUND_COLOR; //Absorbed by medium

        return color_ray(intersection, new_dir, scene, depth+1, t_min, t_max)*hit->albedo;
        //auto new_dir = lambertian_sample(normal);
        // return color_ray(intersection, new_dir, scene, depth+1, t_min, t_max)*0.5;
        // return (normal+1)/2;
    }

    //Nothing Hit
    return color_free_ray(source, dir);
}


int main(int argc, char** argv){

    Vec3D zero;
    print_vec3d(zero, "ZERO: ");
    print_vec3d(SKY_COLOR, "SKY COLOR: ");
    print_vec3d(BACKGROUND_COLOR, "BACKGROUND_COLOR: ");

    vector<double> img(H*W*3);
    LambertianMaterial my_mat;

    vector<Sphere> scene = {
        Sphere(0.0, -10000.0, 0.0, 9999.0), //Earth
        Sphere({-1.1, 0.0, 2.0}, 1.0, {0.5, 0.5, 0.5}, make_shared<BlankMaterial>(0.5)),
        Sphere({1.1, 0.0, 2.0}, 1.0, {0.5, 0.5, 0.5}, make_shared<BlankMaterial>(0.5)),
        Sphere({0.0, -0.7, 1.3}, 0.3, {0.5, 0.5, 0.5}, make_shared<BlankMaterial>(0.5)),
    };

    for(int i = 0;i<H;i++){
        for(int j = 0;j<W;j++){
            auto uv = inds2uv(i, j);
            Vec3D avg_color;
            for(int k = 0;k<SAMPLES_PER_PIXEL;k++){
                auto ray = uv2ray(uv.first + U(RNG), uv.second + U(RNG));
                avg_color += color_ray(zero, ray, scene, 0, 0.0001)/SAMPLES_PER_PIXEL;

            }
            for(int k=0;k<3;k++) img[(i*W+j)*3+k] = avg_color[k];
        }
    }

    ofstream output;
    output.open("img.bin");
    int shape[3] = {H, W, 3};
    output.write(reinterpret_cast<char*>(&shape), sizeof(int)*3);
    output.write(reinterpret_cast<char*>(img.data()), sizeof(double)*img.size());
    std::cout << "Done"<<endl;
}
