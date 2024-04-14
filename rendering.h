#ifndef rendering
#define rendering

#include <random>
#include <memory>
#include "vec3d.h"

//Sampling
double uniform_double(){
    static std::random_device rd;
    static std::default_random_engine RNG{rd()};
    static std::uniform_real_distribution<double> U;
    return U(RNG);
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
        if(dot(res, normal) < 0) return Vec3D(); //zero ray for absorbtion
        return res;
    }
};

struct Sphere {
    Vec3D center;
    double radius;
    Vec3D albedo; //rgb albedo
    std::shared_ptr<MaterialBase> material; //Material properties for diffusion/reflection/refraction

    Sphere(double cx, double cy, double cz, double R) :
        center(cx, cy, cz), radius(R), albedo(0.5, 0.5, 0.5), material(std::make_shared<LambertianMaterial>()) { }
    Sphere(std::vector<double> center, double R, std::vector<double> albedo, std::shared_ptr<MaterialBase> mat) :
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

class Camera {
    public:
        const int MAX_DEPTH;
        const int H, W;
        const double fx, fy, cx, cy;
        const Vec3D SKY_COLOR, BACKGROUND_COLOR;
        // const int MAX_DEPTH = 5;
        // const int H = 900, W = 1600;
        // const double fx = 500.0, fy = 500.0, cx = 800, cy = 450;
        // const Vec3D SKY_COLOR, BACKGROUND_COLOR;

        Camera() : H(900), W(1600), fx(500.0), fy(500.0), cx(800), cy(450), SKY_COLOR(1,1,1), BACKGROUND_COLOR(0,0,0), MAX_DEPTH(5) {}
        Camera(int h, int w, double fx, double fy, double cx, double cy, Vec3D sky_color, Vec3D background_color, int max_depth) : H(h), W(w), fx(fx), fy(fy), cx(cx), cy(cy), SKY_COLOR(sky_color), BACKGROUND_COLOR(background_color), MAX_DEPTH(max_depth) {}

        //uv coordinates to lower left corner of pixel
        std::pair<double, double> inds2uv(int i, int j){
            return std::pair<double, double>{j, H-1-i};
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

        Vec3D color_ray(const Vec3D& source,const Vec3D& dir, const std::vector<Sphere>& scene, double t_min = 0.0, double t_max = INFINITY, int depth = 0){
            if(depth > MAX_DEPTH)
                return BACKGROUND_COLOR;

            double closest_hit = t_max;
            const Sphere* hit;

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

                return color_ray(intersection, new_dir, scene, t_min, t_max, depth+1)*hit->albedo;
                //auto new_dir = lambertian_sample(normal);
                // return color_ray(intersection, new_dir, scene, t_min, t_max, depth+1)*0.5;
                // return (normal+1)/2;
            }

            //Nothing Hit
            return color_free_ray(source, dir);
        }

        std::vector<double> render(const std::vector<Sphere>& scene, int samples_per_pixel = 5) {
            std::vector<double> img(H*W*3);
            Vec3D zero;

            for(int i = 0;i<H;i++){
                for(int j = 0;j<W;j++){
                    auto uv = inds2uv(i, j);
                    Vec3D avg_color;
                    for(int k = 0;k<samples_per_pixel;k++){
                        auto ray = uv2ray(uv.first + uniform_double(), uv.second + uniform_double());
                        avg_color += color_ray(zero, ray, scene, 0.0001)/samples_per_pixel;
                    }
                    for(int k=0;k<3;k++) img[(i*W+j)*3+k] = avg_color[k];
                }
            }
            return img;
        }

};

#endif