#ifndef rendering
#define rendering

#include "geometry.h"

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