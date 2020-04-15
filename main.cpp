#include <iostream>
#include <limits>
#include <fstream>
#include <cassert>
#include <cmath>
#include <string>
#include <cstring>
#include <vector>
#include "vect.h"
#include "omp.h"

#define BLACK  Vect3D(0,     0,   0, true)
#define WHITE  Vect3D(255, 255, 255, true)
#define GREEN  Vect3D(0,   255,   0, true)
#define YELLOW Vect3D(255, 255,   0, true)
#define BLUE   Vect3D(0,     0, 255, true)
#define NYA    Vect3D(0,   255, 148, true)

//      material       color                                 ambientColor                 specColor                    ambS  difS  glow   speS  rflx   rfcS  rfcC
#define GUNMETALL      Material(Vect3D(95,  102, 103, true), Vect3D(41,  52,  57,  true), Vect3D(192, 192, 192, true), 0.2f, 1.2f, 10.f,  1.0f, 0.30f, 0.0f, 1.5f)
#define GLASS          Material(WHITE,                       Vect3D(192, 192, 192, true), Vect3D(192, 192, 192, true), 0.1f, 0.0f, 60.f,  1.0f, 0.00f, 1.0f, 1.5f)
#define SPACEGRAYCHECK Material(Vect3D(37,  40,  42,  true), Vect3D(37,  40,  42,  true), Vect3D(192, 192, 192, true), 0.0f, 1.0f, 100.f, 0.2f, 0.1f,  0.0f, 0.0f)
#define GOLDCHECK      Material(Vect3D(223, 204, 183, true), Vect3D(223, 204, 183, true), Vect3D(223, 204, 183, true), 0.0f, 1.0f, 50.f,  1.0f, 0.3f,  0.0f, 0.0f)
#define MIRROR         Material(WHITE,                       WHITE,                       Vect3D(192, 192, 192, true), 0.1f, 0.0f, 50.f,  1.0f, 1.f,   0.0f, 0.0f)
#define GOLDY          Material(Vect3D(255, 245, 160, true), Vect3D(249, 176, 81,  true), Vect3D(255, 232, 138, true), 0.3f, 1.0f, 80.f,  1.0f, 0.02f, 0.0f, 0.0f)
#define GREENY         Material(Vect3D(103, 183, 156, true), Vect3D(53,  106, 88,  true), Vect3D(141, 210, 179, true), 0.3f, 1.0f, 80.f,  1.0f, 0.1f,  0.0f, 0.0f)
#define PURPLE         Material(Vect3D(131, 91,  214, true), Vect3D(82,  44,  127, true), Vect3D(116, 77,  158, true), 0.3f, 1.0f, 120.f, 0.2f, 0.01f, 0.0f, 0.0f)
#define ROSY           Material(Vect3D(177, 43,  70,  true), Vect3D(155, 38,  58,  true), Vect3D(202, 90,  115, true), 0.3f, 1.0f, 30.f,  0.6f, 0.01f, 0.0f, 0.0f)

class Light {
  private:
    Vect3D point_;
    float intensity_;
  public:
    Light(const Vect3D& point=Vect3D(0.f, 0.f, 0.f), float intensity=1): point_(point), intensity_(intensity) {}
    
    Light(const Light& other): point_(other.point_), intensity_(other.intensity_) {}
    
    Light& operator=(const Light& other) {
        point_ = other.point_;
        intensity_ = other.intensity_;
        return *this;
    }
    const Vect3D& get_point() const {
        return point_;
    }
    float get_intensity() const {
        return intensity_;
    }
};

class Material {
  private:
    Vect3D color_;
    Vect3D ambientColor_;
    Vect3D specularColor_;
    float ambientSensitivity_; // ambient = ambientColor_ * ambientSensitivity_
    float diffuseSensitivity_; // diffuse = color_ * diffuseSensitivity_ * k(angle)
    float glow_; // exp
    float specularSensitivity_;
    float reflexSensitivity_;
    float refractSensitivity_;
    float refractCoef_;
  public:
    Material(const Vect3D& color=WHITE) {
        color_= color; ambientColor_ = color; 
        specularColor_ = WHITE; ambientSensitivity_ = 0.25f;  
        diffuseSensitivity_ = 1.5f; glow_ = 15.f; 
        specularSensitivity_ = 1.f; reflexSensitivity_ = 1.f;
        refractSensitivity_ = 0.f; refractCoef_ = 1.5f;

        assert(color_.is_color()); 
        assert(ambientColor_.is_color());
    }
    Material(const Vect3D& color, float ambientSensitivity, float glow) { 
        color_= color; ambientColor_ = color; 
        specularColor_ = WHITE; ambientSensitivity_ = ambientSensitivity;  
        diffuseSensitivity_ = 1.f; glow_ = glow; 
        specularSensitivity_ = 1.f; reflexSensitivity_ = 1.f;
        refractSensitivity_ = 0.f; refractCoef_ = 1.5f;
        
        assert(color_.is_color()); 
        assert(ambientColor_.is_color());              
    }
    Material(const Vect3D& color, const Vect3D& ambientColor, const Vect3D& specularColor, float ambientSensitivity, 
             float diffuseSensitivity, float glow, float specularSensitivity, float reflexSensitivity, float refractSensitivity, float refractCoef) 
    {
        color_= color; ambientColor_ = ambientColor; 
        specularColor_ = specularColor; ambientSensitivity_ = ambientSensitivity;  
        diffuseSensitivity_ = diffuseSensitivity; glow_ = glow; 
        specularSensitivity_ = specularSensitivity; reflexSensitivity_ = reflexSensitivity;
        refractSensitivity_ = refractSensitivity; refractCoef_ = refractCoef;
        
        assert(color_.is_color()); 
        assert(ambientColor_.is_color());
        assert(specularColor_.is_color());    
    }
    Material(const Material& other) {
        color_= other.color_; ambientColor_ = other.ambientColor_; 
        specularColor_ = other.specularColor_; ambientSensitivity_ = other.ambientSensitivity_;  
        diffuseSensitivity_ = other.diffuseSensitivity_; glow_ = other.glow_; 
        specularSensitivity_ = other.specularSensitivity_; reflexSensitivity_ = other.reflexSensitivity_;
        refractSensitivity_ = other.refractSensitivity_; refractCoef_ = other.refractCoef_;
    }
    Material& operator=(const Material& other) {
        color_= other.color_; ambientColor_ = other.ambientColor_; 
        specularColor_ = other.specularColor_; ambientSensitivity_ = other.ambientSensitivity_;  
        diffuseSensitivity_ = other.diffuseSensitivity_; glow_ = other.glow_; 
        specularSensitivity_ = other.specularSensitivity_; reflexSensitivity_ = other.reflexSensitivity_;
        refractSensitivity_ = other.refractSensitivity_; refractCoef_ = other.refractCoef_;
        return *this;
    }
    const Vect3D& get_color() const {
        return color_;
    }
    Vect3D get_ambient() const {
        return ambientColor_ * ambientSensitivity_;
    }
    Vect3D get_diffuse() const {
        return color_ * diffuseSensitivity_;
    }
    Vect3D get_specular() const {
        return specularColor_ * specularSensitivity_;
    }
    float get_glow() const {
        return glow_;
    }
    float get_reflex() const {
        return reflexSensitivity_;
    }
    float get_refract() const {
        return refractSensitivity_;
    }
    float get_refractCoef() const {
        return refractCoef_;
    }
};

class Plane {
  private:
    Vect3D normal_;
    float dist_;
  public:
    Plane(const Vect3D& normal, float dist): dist_(dist), normal_(normal) {} 
    
    bool ray_intersect(const Vect3D& point, const Vect3D& dirctn, float& t) const;
    
    const Vect3D& get_normal() const {
        return normal_;
    }
    float get_dist() const {
        return dist_;
    }
};
bool Plane::ray_intersect(const Vect3D& point, const Vect3D& dirctn, float& t) const {
    Vect3D d = dirctn; d.normalize();
    if (d * normal_ == 0) {
        if (fabs(point * normal_ - dist_) < 1e-5) {
            t = 0.f;
            return true;
        }
        else return false;
    }  
    float t_;  
    t_ = (dist_ - (point * normal_)) / (d * normal_);
    
    if (t_ < 0)
        return false;
    else {
        t = t_; return true;
    }
}
    
class Sphere {
  private:
    float radius_;
    Vect3D center_;
    Material material_;
  public:
    Sphere(): radius_(0), center_(), material_() {}
    Sphere(float r, const Vect3D& c, const Material& m): radius_(r), center_(c), material_(m) {}
    Sphere(const Vect3D& c, float r, const Material& m): radius_(r), center_(c), material_(m) {}
    Sphere(const Sphere& other) {
        radius_ = other.radius_;
        center_ = other.center_;
        material_ = other.material_;
    }
    Sphere& operator=(const Sphere& other) {
        radius_ = other.radius_;
        center_ = other.center_;
        material_ = other.material_;
    }
    
    bool ray_intersect(const Vect3D& point, const Vect3D& dirctn, float& t) const;
    
    const Vect3D& get_color() const {
        return material_.get_color();
    }
    Vect3D get_ambient() const {
        return material_.get_ambient();
    }
    Vect3D get_diffuse() const {
        return material_.get_diffuse();
    }
    Vect3D get_specular() const {
        return material_.get_specular();
    }
    float get_glow() const {
        return material_.get_glow();
    }
    float get_reflex() const {
        return material_.get_reflex();
    }
    float get_refract() const {
        return material_.get_refract();
    }
    float get_refractCoef() const {
        return material_.get_refractCoef();
    }
    const Vect3D& get_center() const {
        return center_;
    }
};

bool Sphere::ray_intersect(const Vect3D& point, const Vect3D& dirctn, float& t) const {
    t = std::numeric_limits<float>::max();
    Vect3D d = dirctn;
    d.normalize();
    Vect3D vpc = center_ - point;
    
    if ((vpc * d) < 0.f) {
        if (vpc.norm() > radius_) 
            return false;
        else if (vpc.norm() == radius_) 
            t = radius_;
        else {
            Vect3D pc = (vpc * d) * d + point;
            float pc_c_dist = (pc - center_).norm();
            float dist = sqrt(radius_ * radius_ - pc_c_dist * pc_c_dist);
            t = dist - (pc - point).norm();   
        }    
        return true;
    }
    else {
        Vect3D pc = (vpc * d) * d + point;
        if ((pc - center_).norm() > radius_) 
            return false;
        else {
            float pc_c_dist = (pc - center_).norm();
            float dist = sqrt(radius_ * radius_ - pc_c_dist * pc_c_dist);
            if (vpc.norm() > radius_)
                t = (pc - point).norm() - dist;
            else 
                t = (pc - point).norm() + dist;
            return true;
        }
    }
}

class Scene {
  private:
    Plane plane_;
    Material checkColor_1;
    Material checkColor_2;
    std::vector<Sphere> spheres_;
    std::vector<Light> lights_;
    Vect3D bgColor_;
    size_t sph_count;
    size_t lights_count;
    size_t reflex_depth_;
    bool refraction_;
  public:
    Scene(const Vect3D& bgColor, size_t reflex_depth=4): plane_(Plane(Vect3D(0, 0, 0), 0)), checkColor_1(), checkColor_2(), spheres_(), lights_(), 
                                                         bgColor_(bgColor), sph_count(0), lights_count(0), reflex_depth_(reflex_depth), refraction_(false) { 
        assert(bgColor_.is_color()); 
    }    
    Vect3D ray_throw(const Vect3D& point, const Vect3D& dirctn, size_t depth) const;
    
    Vect3D refract_throw(const Vect3D& point, const Vect3D& dirctn, const Sphere& sph, float refractSensitivity, Vect3D& refracted_hit) const ;
    
    void push_sphere(const Sphere& sph) {
        sph_count++;
        spheres_.push_back(sph);
    }
    void push_light(const Light& light) {
        lights_count++;
        lights_.push_back(light);
    }
    void set_plane(const Plane& plane) {
        plane_ = plane;
    }
    void set_checkColors(const Material& color_1, const Material& color_2) {
        checkColor_1 = color_1;
        checkColor_2 = color_2;
    }
};

Vect3D Scene::ray_throw(const Vect3D& point, const Vect3D& dirctn, size_t depth) const {
    if (dirctn.is_black() || depth > reflex_depth_) return BLACK;
    bool plane_intersection = false; 
    bool sph_intersection = false;
    
    float min_dist = std::numeric_limits<float>::max();
    size_t min_idx = sph_count;
    for (size_t i = 0; i < sph_count; i++) {
        float curr_dist = 0.f;
        if ((spheres_[i].ray_intersect(point, dirctn, curr_dist)) && (curr_dist < min_dist)) {
            sph_intersection = true;
            min_dist = curr_dist;
            min_idx = i;
        }
    }
    assert(!plane_.get_normal().is_black()); //checking if plane is defined
      
    float plane_dist = std::numeric_limits<float>::max();
    if (plane_.ray_intersect(point, dirctn, plane_dist)) plane_intersection = true; 
    if (!sph_intersection && !plane_intersection) { // no plane or spheres intersection
        if (depth > 0) return bgColor_ * 0.5f;
        return bgColor_;
    }
    else {    
        Sphere curr_sph; Material plane_mat;
        Vect3D d = dirctn; d.normalize();
        Vect3D hit, normal_hit, out_normal;
        
        if (plane_intersection && sph_intersection) {
            if (plane_dist < min_dist) sph_intersection = false;
            else plane_intersection = false;
        }
        if (plane_intersection) {
            hit = d * plane_dist + point;
            plane_mat = hit.check() ? checkColor_1 : checkColor_2;
            normal_hit = plane_.get_normal();
            out_normal = -normal_hit;
        }
        if (sph_intersection) {
            curr_sph = spheres_[min_idx];
            hit = d * min_dist + point;
            normal_hit = (curr_sph.get_center() - hit).normalize();
            out_normal = (hit - curr_sph.get_center()).normalize();
        }
        float colorIntensity = 0.f;
        float glowIntensity = 0.f;
        float power = plane_intersection ? plane_mat.get_glow() : curr_sph.get_glow();

        for (size_t i = 0; i < lights_count; i++) {

            Vect3D light_to_hit = hit - lights_[i].get_point();
            float light_to_hit_dist = light_to_hit.norm();
            light_to_hit.normalize();

            bool shadow = false;
            for (size_t j = 0; j < sph_count; j++) {
                if (j == min_idx) continue;
                float curr_dist = 0.f;
                if (spheres_[j].ray_intersect(lights_[i].get_point(), light_to_hit, curr_dist) && curr_dist < light_to_hit_dist) {
                    if (spheres_[j].get_refract() != 1.f)
                        shadow = true; 
		    break;
                }
            }
            if (shadow) {
                shadow = false;
                continue;
            }

            float curr_colorIntensity = light_to_hit * (-out_normal) * lights_[i].get_intensity();
            if (curr_colorIntensity > 0) colorIntensity += curr_colorIntensity;

            float curr_glowIntensity = ((light_to_hit * out_normal) * out_normal * 2 - light_to_hit) * d;
            if (curr_glowIntensity > 0) {
                curr_glowIntensity = std::pow(curr_glowIntensity, power) * lights_[i].get_intensity();
                glowIntensity += curr_glowIntensity;
            }
        }
        float rfctIntensity = !depth ? colorIntensity : 1.f;
        float bias = d * (out_normal) < 0 ? 1e-5 : -1e-5;
        Vect3D ambt = plane_intersection ? plane_mat.get_ambient() : curr_sph.get_ambient();
        Vect3D diff = plane_intersection ? plane_mat.get_diffuse() : curr_sph.get_diffuse();
        Vect3D spec = plane_intersection ? plane_mat.get_specular() : curr_sph.get_specular();
        Vect3D rfct, rflx;
        if (plane_intersection) {
            rfct = plane_mat.get_refract() ? plane_mat.get_refract() * Scene::ray_throw(hit - bias * out_normal, d.refract(out_normal, 1.f, plane_mat.get_refractCoef()), depth + 1) * rfctIntensity: BLACK;
            rflx = plane_mat.get_reflex() ? plane_mat.get_reflex() * Scene::ray_throw(hit + bias * out_normal, (-d).reflect(out_normal), depth + 1) * colorIntensity : BLACK;
        }
        else {
            rfct = curr_sph.get_refract() ? curr_sph.get_refract() * Scene::ray_throw(hit - bias * out_normal, d.refract(out_normal, 1.f, curr_sph.get_refractCoef()), depth + 1) * rfctIntensity: BLACK;
            rflx = curr_sph.get_reflex() ? curr_sph.get_reflex() * Scene::ray_throw(hit + bias * out_normal, (-d).reflect(out_normal), depth + 1) * colorIntensity : BLACK;
        }
        return ambt + diff * colorIntensity + spec * glowIntensity + rflx + rfct;
    }
}
void draw(const char* fileName, const std::vector<Vect3D> buf, size_t width=1280, size_t height=720) {
    const std::string directory = std::string("./") + std::string(fileName) + std::string(".ppm");
    
    std::ofstream ofs;
    ofs.open(directory);
    ofs << "P6\n" << width << " " << height << "\n255\n";
    
    for (size_t i = 0; i < height * width; ++i) {
        assert(buf[i].is_color());
        for (size_t j = 0; j < 3; j++)   
            ofs << (char)buf[i][j];
    }
    ofs.close();
}

void drawBMP(const char* fileName, unsigned char *buf, size_t width=1280, size_t height=720) {
    const int fileHeaderSize = 14;
    const int infoHeaderSize = 40;
    
    unsigned char bmp_file_header[14] = { 'B', 'M', 0, 0, 0, 0, 0, 0, 0, 0, 54, 0, 0, 0, };
	unsigned char bmp_info_header[40] = { 40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 24, 0, };
	unsigned char bmp_pad[3] = { 0, 0, 0, };
    
    int paddingSize = (4 - (width * 3) % 4) % 4;
    int size = fileHeaderSize + infoHeaderSize + (3 * width + paddingSize) * height;
    
	bmp_file_header[2]  = static_cast<unsigned char>(size      );
	bmp_file_header[3]  = static_cast<unsigned char>(size >>  8);
	bmp_file_header[4]  = static_cast<unsigned char>(size >> 16);
	bmp_file_header[5]  = static_cast<unsigned char>(size >> 24);

	bmp_info_header[4]  = static_cast<unsigned char>(width      );
	bmp_info_header[5]  = static_cast<unsigned char>(width >>  8);
	bmp_info_header[6]  = static_cast<unsigned char>(width >> 16);
	bmp_info_header[7]  = static_cast<unsigned char>(width >> 24);

	bmp_info_header[8]  = static_cast<unsigned char>(height      );
	bmp_info_header[9]  = static_cast<unsigned char>(height >>  8);
	bmp_info_header[10] = static_cast<unsigned char>(height >> 16);
	bmp_info_header[11] = static_cast<unsigned char>(height >> 24);
    
    FILE *file = fopen(fileName, "wb");
    
    if (file)
	{
		fwrite(bmp_file_header, 1, 14, file);
		fwrite(bmp_info_header, 1, 40, file);

		for (int i = 0; i < height; i++)
		{
			fwrite(buf + (width * i * 3), 3, width, file);
			fwrite(bmp_pad, 1, ((4 - (width * 3) % 4) % 4), file);
		}
		fclose(file);
	}
}

void renderCamera(const Scene& scene, const char* fileName="out.bmp", bool BMP=true, size_t antialiasing_ppx=4, size_t num_threads=8, float camFOV=80.f, 
                  const Vect3D& camPos=Vect3D(0.f, 0.f, 0.f), const Vect3D& camDirctn=Vect3D(0.f, 0.f, 1.f)) 
{
    // Camera views along OZ axis (+inf direction)
    const int width = 1280;
    const int height = 720;
    
    assert(width % 2 == 0);
    assert(height % 2 == 0);
    assert(antialiasing_ppx == 1 || antialiasing_ppx == 2 || antialiasing_ppx == 4);
    assert(BMP);
    
    camFOV = camFOV * (M_PI / 180.f);
    const float spaceHeight = 2.f * tan(camFOV / 2.f);
    const float spaceWidth = (spaceHeight / height) * width;
    
    const float x_tick = spaceWidth / width;
    const float y_tick = spaceHeight / height;

    unsigned char raveledImg[height][width][3];

    size_t i, j;
    omp_set_num_threads(num_threads);
    #pragma omp parallel private(i, j)
    #pragma omp for collapse(2) schedule(dynamic)
    for(i = 0; i < height; i++){
        for(j = 0; j < width; j++){

	    float x = j * x_tick + x_tick / 2.f - spaceWidth / 2.f;
	    float y = -(-(i * y_tick) - y_tick / 2.f + spaceHeight / 2.f);

	    Vect3D pixel;

	    if (antialiasing_ppx == 1)
	        pixel = scene.ray_throw(camPos, camDirctn + Vect3D(x, y, 0.f), 0);
	    else if (antialiasing_ppx == 2)
	        pixel = (scene.ray_throw(camPos, camDirctn + Vect3D(x - x_tick / 4, y, 0.f), 0).set_color(false) + 
		         scene.ray_throw(camPos, camDirctn + Vect3D(x + x_tick / 4, y, 0.f), 0).set_color(false)).set_color(true) * (0.5f);
	    else if (antialiasing_ppx == 4)
	        pixel = (scene.ray_throw(camPos, camDirctn + Vect3D(x - x_tick / 4, y - y_tick / 4, 0.f), 0).set_color(false) + 
		         scene.ray_throw(camPos, camDirctn + Vect3D(x - x_tick / 4, y + y_tick / 4, 0.f), 0).set_color(false) +
		         scene.ray_throw(camPos, camDirctn + Vect3D(x + x_tick / 4, y - y_tick / 4, 0.f), 0).set_color(false) + 
		         scene.ray_throw(camPos, camDirctn + Vect3D(x + x_tick / 4, y + y_tick / 4, 0.f), 0).set_color(false)).set_color(true) * (0.25f);

	    raveledImg[i][j][2] = (unsigned char) pixel.red(); 
	    raveledImg[i][j][1] = (unsigned char) pixel.green(); 
	    raveledImg[i][j][0] = (unsigned char) pixel.blue(); 
        }
    }
    drawBMP(fileName, (unsigned char *)raveledImg, width, height);
}

int main(int argc, char** argv) {
    const char* fileName;
    size_t num_threads;
    for (size_t i = 0; i < argc; i++) {
        if (!strcmp(argv[i], "-out"))
            fileName = argv[i + 1];
        else if (!strcmp(argv[i], "-scene")) {
            if ((argv[i + 1][0] - '1')) return 0;
        }
        else if (!strcmp(argv[i], "-threads"))
            num_threads = atoi(argv[i + 1]);
    }
    Scene my_scene(BLACK);
    
    my_scene.push_sphere(Sphere(Vect3D(-1, 0, 15), 4.f, SPACEGRAYCHECK));
    my_scene.push_sphere(Sphere(Vect3D(1.5, -1.5, 8), 2.f, GLASS));
    my_scene.push_sphere(Sphere(Vect3D(-10, 5, 17), 6.f, MIRROR));
    my_scene.push_sphere(Sphere(Vect3D(5, 5, 20), 3.f, GOLDY));
    my_scene.push_sphere(Sphere(Vect3D(10, 10, 25), 2.5f, PURPLE));
    my_scene.push_sphere(Sphere(Vect3D(-3.5, -2.5, 10), 1.5f, ROSY));
    
    my_scene.set_plane(Plane(Vect3D(0, -1, 0), 4.f));
    my_scene.set_checkColors(GOLDCHECK, SPACEGRAYCHECK);
    
    my_scene.push_light(Light(Vect3D(5, 4, -5), 1.f));
    
    renderCamera(my_scene, fileName, true, 4, num_threads);
    
    return 0;
}
