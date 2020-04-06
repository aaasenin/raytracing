#include <cassert>
#include <iostream>
#include <cmath>

float color_mul(float a, float b) {
    return std::max(0.f, std::min(255.f, a * b));
}
float color_sum(float a, float b) {
    return std::max(0.f, std::min(255.f, a + b));
}
float color_sub(float a, float b) {
    return std::max(0.f, std::min(255.f, a - b));
}
class Vect3D {
    // one class for 3D vectors and 3D colors
  private:
    float x_, y_, z_;
    bool color_; // if vector meant to be color => 0 <= x_ <= 255
  public:
    Vect3D(float x=0, float y=0, float z=0, bool color=false);
    
    Vect3D(const Vect3D& other);
	
    Vect3D& operator=(const Vect3D& other);
    
    Vect3D operator-() const;
    
    Vect3D reflect(const Vect3D& normal) const;
    
    Vect3D refract(const Vect3D& normal, float n_1, float n_2) const;
    
    float norm() const { return sqrt(x_*x_ + y_*y_ + z_*z_); }
    
    float red() const { return x_; }
    
    float green() const { return y_; }
    
    float blue() const  { return z_; }
    
    bool check() const;
    
    bool is_color() const;
    
    bool is_black() const;
    
    bool is_insideBox(const Vect3D& vmin, const Vect3D& vmax) const;
    
    const Vect3D& set_color(bool color);
    
    Vect3D& normalize();

    float operator[](size_t i) const;
	
    friend Vect3D operator+(const Vect3D& a, const Vect3D& b);
    
    friend Vect3D operator-(const Vect3D& a, const Vect3D& b);
	
    friend float operator*(const Vect3D& a, const Vect3D& b);
	
    friend Vect3D operator*(float a, const Vect3D& b);
	
    friend Vect3D operator*(const Vect3D& a, float b);
	
    friend std::ostream& operator<<(std::ostream& out, const Vect3D& v);
};
Vect3D(float x=0, float y=0, float z=0, bool color=false) {
    x_ = x; y_ = y; z_ = z; color_ = color;
}
Vect3D(const Vect3D& other) {
    x_ = other.x_; y_ = other.y_; z_ = other.z_; color_ = other.color_;
}
Vect3D& operator=(const Vect3D& other) {
    x_ = other.x_; y_ = other.y_; z_ = other.z_; color_ = other.color_;
    return *this;
}
Vect3D operator-() const {
    assert(!color_);
    return Vect3D(-x_, -y_, -z_);
}
Vect3D reflect(const Vect3D& normal) const {
    float bias;
    if ((*this) * normal > 0) bias = 1e-3;
    else bias = -1e-3;
    return ((*this) * normal) * normal * 2 - (*this) + bias * normal;
}
Vect3D refract(const Vect3D& normal, float n_1, float n_2) const {
    Vect3D n = normal; n.normalize();
    Vect3D l = *this; l.normalize();

    float c_1 = - (n * l);
    float r = n_1 / n_2;
    if (c_1 < 0) {
        c_1 = -c_1;
        n = -n;
        r = 1.f / r;
    }
    float c_2 = 1 - r * r * (1 - c_1 * c_1);
    if (c_2 <= 0 || c_1 == 0) 
        return Vect3D(0, 0, 0);
    return r * (*this) + (r * c_1 - sqrt(c_2)) * n;       
} 
bool check() const {
    return (int(x_ * 0.5f + 10000) + int(z_ * 0.5f)) & 1 ? true : false;
}
bool is_color() const { 
    assert((x_ >= 0) && (y_ >= 0) && (z_ >= 0));
    assert((x_ <= 255) && (y_ <= 255) && (z_ <= 255));
    return color_; 
}
bool is_black() const {
    if (x_ == 0.f && y_ == 0.f && z_ == 0.f) return true;
    else return false;
}
bool is_insideBox(const Vect3D& vmin, const Vect3D& vmax) const {
    if (x_ >= vmin.x_ && y_ >= vmin.y_ && z_ >= vmin.z_ && x_ <= vmax.x_ && y_ <= vmax.y_ && z_ <= vmax.z_)
        return true;
    return false;
}
const Vect3D& set_color(bool color) {
    color_ = color;
    return *this;
}
Vect3D& Vect3D::normalize() {
    assert(!(this->is_black()));
	float norm_ = (*this).norm();
    norm_ = norm_ ? (1.f / norm_) : 0.f;
	x_ = x_ * norm_; y_ = y_ * norm_; z_ = z_ * norm_;
	return *this;
}
float Vect3D::operator[](size_t i) const {
	assert(i < 3);
	return i == 0 ? x_ : (i == 1 ? y_ : z_);
}
Vect3D operator+(const Vect3D& a, const Vect3D& b) {
    if ((a.color_) || (b.color_))
	    return Vect3D(color_sum(a.x_, b.x_), color_sum(a.y_, b.y_), color_sum(a.z_, b.z_), true);
    else return Vect3D(a.x_ + b.x_, a.y_ + b.y_, a.z_ + b.z_);
}
Vect3D operator-(const Vect3D& a, const Vect3D& b) {
    if ((a.color_) || (b.color_))
	    return Vect3D(color_sub(a.x_, b.x_), color_sub(a.y_, b.y_), color_sub(a.z_, b.z_), true);
    else return Vect3D(a.x_ - b.x_, a.y_ - b.y_, a.z_ - b.z_);
}
float operator*(const Vect3D& a, const Vect3D& b) {
	return a.x_ * b.x_ + a.y_ * b.y_ + a.z_ * b.z_;
}
Vect3D operator*(float a, const Vect3D& b) {
    if (b.color_) 
        return Vect3D(color_mul(a, b.x_), color_mul(a, b.y_), color_mul(a, b.z_), true);
	return Vect3D(a * b.x_, a * b.y_, a * b.z_);
}
Vect3D operator*(const Vect3D& a, float b) {
    if (a.color_)
        return Vect3D(color_mul(a.x_, b), color_mul(a.y_, b), color_mul(a.z_, b), true);
	return Vect3D(a.x_ * b, a.y_ * b, a.z_ * b);
}
std::ostream& operator<<(std::ostream& out, const Vect3D& v) {
	out << '(' << v.x_ << ", " << v.y_  << ", " << v.z_ << ')';
	return out;
}
