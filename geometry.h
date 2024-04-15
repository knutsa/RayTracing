#ifndef geometry
#define geometry

#include <random>
#include <memory>
#include "vec3d.h"

#include <cassert>

//Sampling
double uniform_double(){
    static std::random_device rd;
    static std::default_random_engine RNG{rd()};
    static std::uniform_real_distribution<double> U;
    return U(RNG);
}
double uniform_double(double a, double b){
    return a + (b-a)*uniform_double();
}
Vec3D unit_sphere_sample(){
    double z = 2*uniform_double()-1, phi = uniform_double()*2*M_PI;
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


//Objects in scene
enum class MaterialTypes {UniformDiffusion, Lambertian, Blank};
int encode_mat_type(MaterialTypes mt){
    if(mt == MaterialTypes::UniformDiffusion) return 0;
    if(mt == MaterialTypes::Lambertian) return 1;
    if(mt == MaterialTypes::Blank) return 2;
    return -1;
}
struct MaterialBase { //Abstract base struct for optic material
    virtual Vec3D scatter(const Vec3D& normal, const Vec3D& ray_dir) = 0;
    virtual MaterialTypes MAT_TYPE() = 0;
    virtual std::vector<double> parameter_vec() { return {}; } //used to save/communicate a scene
    MaterialBase() {}
};
struct UniDiffusionMaterial : MaterialBase {
    Vec3D scatter(const Vec3D& normal, const Vec3D& ray_dir) override {
        return unit_hemisphere_sample(normal);
    }
    UniDiffusionMaterial() {}
    MaterialTypes MAT_TYPE() override {return MaterialTypes::UniformDiffusion;}
};
struct LambertianMaterial : MaterialBase {
    Vec3D scatter(const Vec3D& normal, const Vec3D& ray_dir) override {
        return lambertian_sample(normal);
    }
    LambertianMaterial() {}
    MaterialTypes MAT_TYPE() override {return MaterialTypes::Lambertian; }
};
struct BlankMaterial : MaterialBase {
    double fuzziness;

    BlankMaterial() : fuzziness(0.0) {}
    BlankMaterial(double fuzziness) : fuzziness(fuzziness) {}

    MaterialTypes MAT_TYPE() override {return MaterialTypes::Blank;}
    Vec3D scatter(const Vec3D& normal, const Vec3D& ray_dir) override {

        auto res = ray_dir - normal*(2*dot(ray_dir, normal));
        if(fuzziness == 0) return res;
        res += unit_sphere_sample()*fuzziness;
        if(dot(res, normal) < 0) return Vec3D(); //zero ray for absorbtion
        return res;
    }
    // std::vector<double> parameter_vec() override { return {fuzziness}; }
    std::vector<double> parameter_vec() override {return {fuzziness};}
};

struct Sphere {
    Vec3D center;
    double radius;
    Vec3D albedo; //rgb albedo
    std::shared_ptr<MaterialBase> material; //Material properties for diffusion/reflection/refraction

    Sphere(double cx, double cy, double cz, double R) :
        center(cx, cy, cz), radius(R), albedo(0.5, 0.5, 0.5), material(std::make_shared<LambertianMaterial>()) { }
    Sphere(Vec3D center, double R, Vec3D albedo, std::shared_ptr<MaterialBase> mat) :
        center(center), radius(R), albedo(albedo), material(mat) {}


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

std::vector<Sphere> generate_random_scene(int n_small_spheres, double small_r = 0.1) {
    std::vector<Sphere> scene;
    //Big Spheres
    scene.push_back(Sphere(0.0, -100000.0, 0.0, 99999.0)); //Earth
    scene.push_back(Sphere({-1.1, 0.0, 3.0}, 1.0, {0.5, 0.5, 0.5}, std::make_shared<BlankMaterial>(0.0)));
    scene.push_back(Sphere({1.1, 0.0, 3.0}, 1.0, {0.5, 0.5, 0.5}, std::make_shared<BlankMaterial>(0.0)));
    int k = 0, attempts = 0;
    while(k < n_small_spheres && attempts<n_small_spheres*1000){

        Vec3D pos(uniform_double(-3.0, 3.0), -1.0+small_r, uniform_double(0.5, 6.0));
        bool collides = false;
        for(int i = 1;i<scene.size();i++){
            if((pos-scene[i].center).norm() <= small_r*1.05 + scene[i].radius)
                collides = true;
        }

        if(!collides){
            Vec3D albedo(uniform_double(), uniform_double(), uniform_double());
            if(uniform_double() < 0.8){ //Diffuse
                scene.push_back(Sphere(pos, small_r, albedo, std::make_shared<LambertianMaterial>()));
            }
            else{
                scene.push_back(Sphere(pos, small_r, albedo, std::make_shared<BlankMaterial>(uniform_double()*0.5)));
            }
            k++;
        }
        attempts++;
    }
    return scene;
}
using SceneParams = std::pair<std::vector<double>, std::vector<int> >;
SceneParams save_scene_params(const std::vector<Sphere>& scene){
    std::vector<double> data;
    std::vector<int> material_types;

    for(int i = 0;i<scene.size();i++){
        data.push_back(scene[i].center[0]);
        data.push_back(scene[i].center[1]);
        data.push_back(scene[i].center[2]);
        data.push_back(scene[i].radius);
        data.push_back(scene[i].albedo[0]);
        data.push_back(scene[i].albedo[1]);
        data.push_back(scene[i].albedo[2]);
        if(scene[i].material->MAT_TYPE() == MaterialTypes::Blank){
            assert(scene[i].material->parameter_vec().size() == 1);
            data.push_back(scene[i].material->parameter_vec()[0]);
        }
        material_types.push_back(encode_mat_type(scene[i].material->MAT_TYPE()));
    }

    return std::pair<std::vector<double>, std::vector<int> >(data, material_types);
}
std::vector<Sphere> parse_scene_params(SceneParams scene_params) {
    auto data = scene_params.first;
    auto material_types = scene_params.second;
    std::vector<Sphere> res;

    int n = material_types.size(), ind = 0;
    for(int i = 0;i<n;i++){
        if(material_types[i] == encode_mat_type(MaterialTypes::UniformDiffusion) ){ // UniDiffuse
            res.push_back(
                Sphere({data[ind], data[ind+1], data[ind+2]}, data[ind+3], {data[ind+4], data[ind+5], data[ind+6]}, std::make_shared<UniDiffusionMaterial>())
            );
            ind += 7;
        }
        if(material_types[i] == encode_mat_type(MaterialTypes::Lambertian) ){ // Lambertian Diffusive
            res.push_back(
                Sphere({data[ind], data[ind+1], data[ind+2]}, data[ind+3], {data[ind+4], data[ind+5], data[ind+6]}, std::make_shared<LambertianMaterial>())
            );
            ind += 7;
        }
        if(material_types[i] == encode_mat_type(MaterialTypes::Blank)){ // Blank
            res.push_back(
                Sphere({data[ind], data[ind+1], data[ind+2]}, data[ind+3], {data[ind+4], data[ind+5], data[ind+6]}, std::make_shared<BlankMaterial>(data[ind+7]))
            );
            ind += 8;
        }
    }

    return res;
}
int max_material_param_size(){
    return 8;
}

#endif