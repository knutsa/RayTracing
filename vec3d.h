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
        double illuminance() const {return (arr[0] * 0.3) + (arr[1] * 0.59) + (arr[2] * 0.11);}
};
double dot(const Vec3D& x, const Vec3D& y){ return x[0]*y[0]+x[1]*y[1]+x[2]*y[2]; }
void print_vec3D(const Vec3D& x, std::string caption){
    std::cout << caption << "(" <<x[0] << ","<<x[1]<<","<<x[2]<<")"<<std::endl;
}

class Mat3D {
    public:
        double arr[9];
        Mat3D(double m[9]) { for(int i = 0;i<9;i++) arr[i]=m[i]; }
        Mat3D(double a11, double a12, double a13, double a21, double a22, double a23, double a31, double a32, double a33) {
            arr[0] = a11; arr[1] = a12; arr[2] = a13;
            arr[3] = a21; arr[4] = a22; arr[5] = a23;
            arr[6] = a31; arr[7] = a32; arr[8] = a33;
        }
        //Defaults to 3x3 Identity/eye
        Mat3D() { for(int i=0;i<9;i++) arr[i]=((i%3)==i/3);}

        double at(int i, int j) const {return arr[i*3+j];}


        Mat3D operator=(const double data[9]) { for(int i=0;i<9;i++) arr[i]=data[i]; return *this; }
        Mat3D operator+(const Mat3D& other) const {
            Mat3D res;
            for(int i = 0;i<3;i++)
                for(int j = 0;j<3;j++)
                    res.arr[i*3+j] = arr[i*3+j] + other.arr[i*3+j];
            return res;
        }
        Mat3D operator*(double a) const {
            Mat3D res;
            for(int i = 0;i<3;i++)
                for(int j = 0;j<3;j++)
                    res.arr[i*3+j] = arr[i*3+j]*a;
            return res;
        }
        Mat3D operator-(const Mat3D& other) const { return (*this)+(other*-1); }

        Vec3D transform(const Vec3D& x) const {
            return Vec3D(
                arr[0]*x[0] + arr[1]*x[1] + arr[2]*x[2],
                arr[3]*x[0] + arr[4]*x[1] + arr[5]*x[2],
                arr[6]*x[0] + arr[7]*x[1] + arr[8]*x[2]
            );
        }
};

void print_mat3D(const Mat3D& A, std::string caption){
    std::cout << caption << std::endl;
    for(int i = 0;i<3;i++){
        for(int j = 0;j<3;j++)
            std::cout << A.at(i, j) << " ";
        std::cout << std::endl;
    }
}

Mat3D matmul3D(const Mat3D& A, const Mat3D& B){
    Mat3D res;
    for(int i = 0;i<3;i++)
        for(int j = 0;j<3;j++)
            for(int k = 0;k<3;k++)
                res.arr[i*3+j] += A.at(i,k)*B.at(k,j);

    return res;
}
Mat3D Rodriguez_formula(const Vec3D& rvec){
    double alpha = rvec.norm();
    if(alpha == 0) return Mat3D();
    auto k = rvec/alpha;
    Mat3D K = {
        0.0, -k[2], k[1],
        k[2], 0.0, -k[0],
        -k[1], k[0], 0.0       
    };

    return Mat3D() + K*sin(alpha) + matmul3D(K, K)*(1-cos(alpha));
}

#endif