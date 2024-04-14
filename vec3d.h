#ifndef vec3d
#define vec3d

#include <iostream>
#include <math.h>

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
void print_vec3d(const Vec3D& x, std::string caption){
    std::cout << caption << "(" <<x[0] << ","<<x[1]<<","<<x[2]<<")"<<std::endl;
}

#endif