#include "trace.H"
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <getopt.h>
#ifdef __APPLE__
#define MAX std::numeric_limits<double>::max()
#else
#include <values.h>
#define MAX DBL_MAX
#endif

// return the determinant of the matrix with columns a, b, c.
double det(const SlVector3 &a, const SlVector3 &b, const SlVector3 &c) {
    return a[0] * (b[1] * c[2] - c[1] * b[2]) +
        b[0] * (c[1] * a[2] - a[1] * c[2]) +
            c[0] * (a[1] * b[2] - b[1] * a[2]);
}

inline double sqr(double x) {return x*x;} 

bool checkRange(const double& _inValue, const double& _min, const double& _max) {
    if (_inValue >= _min && _inValue <= _max) 
        return true;
    
    return false;
}

bool Triangle::intersect(const Ray &r, double t0, double t1, HitRecord &hr) const {

    // Step 1 Ray-triangle test
    hr.bHit = false;
    double aTemp, bTemp, cTemp, d, e, f, g, h, i, j, k, l;
    SlVector3 col1, col2, col3, col4;
    bool bCheck = false;
    aTemp = this->a.x() - this->b.x();
    bTemp = this->a.y() - this->b.y();
    cTemp = this->a.z() - this->b.z();
    col1.set(aTemp, bTemp, cTemp);

    d = this->a.x() - this->c.x();
    e = this->a.y() - this->c.y();
    f = this->a.z() - this->c.z();
    col2.set(d, e, f);

    g = r.d.x();
    h = r.d.y();
    i = r.d.z();
    col3.set(g, h, i);

    j = this->a.x() - r.e.x();
    k = this->a.y() - r.e.y();
    l = this->a.z() - r.e.z();
    col4.set(j, k, l);

    double divisor = det(col1, col2, col3);
    double betaD = det(col4, col2, col3);
    double gammaD = det(col1, col4, col3);
    double tD = det(col1, col2, col4);

    if (divisor == 0.f) return false;
    
    hr.t = tD / divisor; 
    bCheck = checkRange(hr.t, t0, t1);
    //std::cout << "hrt: " << bCheck << " t: " << hr.t; 
    if (!bCheck) return false;

    hr.gamma = gammaD / divisor;
    bCheck = checkRange(hr.gamma, 0, 1);
    //std::cout << "hr gamma: " << bCheck << hr.gamma; 
    if (!bCheck) return false;

    hr.beta = betaD / divisor;
    bCheck = checkRange(hr.beta, 0, (1 - hr.gamma));
    //std::cout << "hrbeta: " << bCheck << hr.beta; 
    if (!bCheck) return false;

    hr.alpha = 1 - hr.beta - hr.gamma;

    // UpdateHitRecord
    hr.bHit = true;
    
    hr.n = cross(b-a, c-a);
    normalize(hr.n);


    return hr.bHit;
}

bool TrianglePatch::intersect(const Ray &r, double t0, double t1, HitRecord &hr) const {
    bool temp = Triangle::intersect(r,t0,t1,hr);
    if (temp) {
        hr.n = hr.alpha * n1 + hr.beta * n2 + hr.gamma * n3;
        normalize(hr.n);
        //std::cout << hr.n.x() << hr.n.y() << hr.n.z();
    }
    return false;
}


bool Sphere::intersect(const Ray &r, double t0, double t1, HitRecord &hr) const {

    // Step 1 Sphere-triangle test
    double a, b, cc, delta;
    hr.bHit = false;
    //bool bHit = false
    a = dot(r.d, r.d);
    b = 2 * dot(r.d, (r.e - this->c));
    cc = dot((r.e - this->c), (r.e - this->c)) - rad * rad;
    delta = b * b - 4 * a * cc;
    if (delta >= 0) {
        hr.t = (-b-sqrt(delta))/(2*a); // Need smallest greater than t.
        if (checkRange(hr.t, t0, t1)) {
            hr.bHit = true;
            hr.p = r.e + hr.t * r.d;
            hr.n = hr.p - this->c;
            normalize(hr.n);   
        }
    }

    
    return hr.bHit;
}




Tracer::Tracer(const std::string &fname) {
    std::ifstream in(fname.c_str(), std::ios_base::in);
    std::string line;
    char ch;
    Fill fill;
    bool coloredlights = false;
    while (in) {
        getline(in, line);
        switch (line[0]) {
            case 'b': {
                std::stringstream ss(line);
                ss>>ch>>bcolor[0]>>bcolor[1]>>bcolor[2];
                break;
            }

            case 'v': {
                getline(in, line);
                std::string junk;
                std::stringstream fromss(line);
                fromss>>junk>>eye[0]>>eye[1]>>eye[2];

                getline(in, line);
                std::stringstream atss(line);
                atss>>junk>>at[0]>>at[1]>>at[2];

                getline(in, line);
                std::stringstream upss(line);
                upss>>junk>>up[0]>>up[1]>>up[2];

                getline(in, line);
                std::stringstream angless(line);
                angless>>junk>>angle;

                getline(in, line);
                std::stringstream hitherss(line);
                hitherss>>junk>>hither;

                getline(in, line);
                std::stringstream resolutionss(line);
                resolutionss>>junk>>res[0]>>res[1];
                break;
            }

            case 'p': {
                bool patch = false;
                std::stringstream ssn(line);
                unsigned int nverts;
                if (line[1] == 'p') {
                    patch = true;
                    ssn>>ch;
                }
                ssn>>ch>>nverts;
                std::vector<SlVector3> vertices;
                std::vector<SlVector3> normals;
                for (unsigned int i=0; i<nverts; i++) {
                    getline(in, line);
                    std::stringstream ss(line);
                    SlVector3 v,n;
                    if (patch) ss>>v[0]>>v[1]>>v[2]>>n[0]>>n[1]>>n[2];
                    else ss>>v[0]>>v[1]>>v[2];
                    vertices.push_back(v);
                    normals.push_back(n);
                }
                bool makeTriangles = false;
                if (vertices.size() == 3) {
                    if (patch) {
                        surfaces.push_back(std::pair<Surface *, Fill>(new TrianglePatch(vertices[0], vertices[1], vertices[2], 
                        normals [0], normals [1], normals [2]), fill));
                    } else {
                        surfaces.push_back(std::pair<Surface *, Fill>(new Triangle(vertices[0], vertices[1], vertices[2]), fill));
                    }
                } else if (vertices.size() == 4) {
                    SlVector3 n0 = cross(vertices[1] - vertices[0], vertices[2] - vertices[0]);
                    SlVector3 n1 = cross(vertices[2] - vertices[1], vertices[3] - vertices[1]);
                    SlVector3 n2 = cross(vertices[3] - vertices[2], vertices[0] - vertices[2]);
                    SlVector3 n3 = cross(vertices[0] - vertices[3], vertices[1] - vertices[3]);
                    if (dot(n0,n1) > 0 && dot(n0,n2) > 0 && dot(n0,n3) > 0) {
                        makeTriangles = true;
                        if (patch) {
                            surfaces.push_back(std::pair<Surface *, Fill>(new TrianglePatch(vertices[0], vertices[1], vertices[2], 
                            normals[0], normals[1], normals[2]), fill));
                            surfaces.push_back(std::pair<Surface *, Fill>(new TrianglePatch(vertices[0], vertices[2], vertices[3], 
                            normals[0], normals[2], normals[3]), fill));
                        } else {
                            surfaces.push_back(std::pair<Surface *, Fill>(new Triangle(vertices[0], vertices[1], vertices[2]), fill));
                            surfaces.push_back(std::pair<Surface *, Fill>(new Triangle(vertices[0], vertices[2], vertices[3]), fill));
                        }
                    }
                    if (!makeTriangles) {
                        std::cerr << "I didn't make triangles.  Poly not flat or more than quad.\n";
                    }
                }
                break;
            }

            case 's' : {
                std::stringstream ss(line);
                SlVector3 c;
                double r;
                ss>>ch>>c[0]>>c[1]>>c[2]>>r;
                surfaces.push_back(std::pair<Surface *, Fill>(new Sphere(c,r), fill));
                break;
            }
	  
            case 'f' : {
                std::stringstream ss(line);
                ss>>ch>>fill.color[0]>>fill.color[1]>>fill.color[2]>>fill.kd>>fill.ks>>fill.shine>>fill.t>>fill.ior;
                break;
            }

            case 'l' : {
                std::stringstream ss(line);
                Light l;
                ss>>ch>>l.p[0]>>l.p[1]>>l.p[2];
                if (!ss.eof()) {
                    ss>>l.c[0]>>l.c[1]>>l.c[2];
                    coloredlights = true;
                }
                lights.push_back(l);
                break;
            }

            default:
            break;
        }
    }
    if (!coloredlights) for (unsigned int i=0; i<lights.size(); i++) lights[i].c = 1.0/sqrt(lights.size());
    im = new SlVector3[res[0]*res[1]];
    shadowbias = 1e-6;
    samples = 1;
    aperture = 0.0;
}

Tracer::~Tracer() {
    if (im) delete [] im;
    for (unsigned int i=0; i<surfaces.size(); i++) delete surfaces[i].first;
}


SlVector3 Tracer::shade(const HitRecord &hr) const {
    if (color) return hr.f.color;

    SlVector3 color(0.0);
    HitRecord dummy;

    for (unsigned int i=0; i<lights.size(); i++) {
        const Light &light = lights[i];
        bool shadow = false;
        bool isShadowOn = true;

        if (isShadowOn) {
            Ray shadowCheck(hr.p, (light.p-hr.p));
            normalize(shadowCheck.d);
            for (int i = 0; i < surfaces.size(); i++){
                shadow = surfaces[i].first->intersect(shadowCheck, shadowbias, mag(light.p - hr.p), dummy);
                if (dummy.bHit) {
                    shadow = true;
                    break;
                }
            }
        }

        if (!shadow) {
            
            // Step 2 do shading here
            SlVector3 l = light.p - hr.p;
            normalize(l);
            SlVector3 r = -l + (2 * dot(l, hr.n)) * hr.n;
            normalize(r);
            SlVector3 ambient = 0.05 * light.c;
            SlVector3 diffuse, specular;
            diffuse = hr.f.kd * light.c * hr.f.color * std::max(dot(l, hr.n), 0.0);
            specular = hr.f.ks * light.c * hr.f.color * std::pow(std::max(dot(r, hr.n), 0.0), hr.f.shine);
            color += (ambient + diffuse + specular);
        }
    }

    return color;
}

SlVector3 Tracer::reflect(const Ray& r, const HitRecord& hr) const {
    SlVector3 blackColor(0.0);
    
    if (hr.raydepth == 0) {
        //std::cout << "TracerReflect Stop";
        return blackColor;
    }

    SlVector3 newDir = r.d - (2 * dot(r.d, hr.n)) * hr.n;
    normalize(newDir);
    Ray ray(hr.p, newDir);
     
    HitRecord hrNew;
    HitRecord dummy;
    hrNew.t = MAX;
    bool hit = false;
    // intersects with each object, find the smallest t
    for (int i = 0; i < surfaces.size(); i++){
        hit = surfaces[i].first->intersect(ray, shadowbias, hrNew.t, dummy);
        dummy.f = surfaces[i].second;
        if (dummy.bHit) {
            hrNew = dummy;
        }
    }  

    // for each object, process light
    if (hrNew.bHit) {
    hrNew.p = ray.e + hrNew.t * ray.d;
    hrNew.v = ray.e - hrNew.p;
    normalize(hrNew.v);
    // add all the lights

    //std::cout << "raydepth: " << hrNew.raydepth;
    blackColor = min(shade(hrNew), SlVector3(1.0));
    hrNew.raydepth = hr.raydepth - 1;
    return blackColor + hrNew.f.ks * reflect(ray, hrNew);
    }

    return blackColor;
}

SlVector3 Tracer::trace(const Ray &r, double t0, double t1) const {
    HitRecord hr;
    hr.bHit = false;
    SlVector3 color(bcolor);
  
    bool hit = false;

    // Step 1 See what a ray hits 
    for (int i = 0; i < surfaces.size(); i++){
        HitRecord hrTemp;
        hrTemp.bHit = false;
        hit = surfaces[i].first->intersect(r, t0, t1, hrTemp);
        hrTemp.f = surfaces[i].second;
        // compare hr.t 
        if (hrTemp.bHit) {
            if (hrTemp.t < t1) {
                hr = hrTemp;
                t1 = hrTemp.t;
            }
        }
    }  

    if (hr.bHit) {   
        // set HitRecord
        hr.p = r.e + hr.t * r.d;
        hr.v = r.e - hr.p;
        normalize(hr.v);
        hr.raydepth = maxraydepth;
        //color = hr.n;
        color = shade(hr);
        color += hr.f.ks * reflect(r, hr);
    }
    return color; 
}

void Tracer::traceImage() {
    // set up coordinate system
    SlVector3 w = eye - at;
    w /= mag(w);
    SlVector3 u = cross(up,w);
    normalize(u);
    SlVector3 v = cross(w,u);
    normalize(v);

    double d = mag(eye - at);
    double h = tan((M_PI/180.0) * (angle/2.0)) * d;
    double l = -h;
    double r = h;
    double b = h;
    double t = -h;

    SlVector3 *pixel = im;

    for (unsigned int j=0; j<res[1]; j++) {
        for (unsigned int i=0; i<res[0]; i++, pixel++) {

            SlVector3 result(0.0,0.0,0.0);

            for (int k = 0; k < samples; k++) {

                double rx = 1.1 * rand() / RAND_MAX;
                double ry = 1.1 * rand() / RAND_MAX;

                double x = l + (r-l)*(i+rx)/res[0];
                double y = b + (t-b)*(j+ry)/res[1];
                SlVector3 dir = -d * w + x * u + y * v;
	
                Ray r(eye, dir);
                normalize(r.d);
                result += trace(r, hither, MAX);

            }
            (*pixel) = result / samples;
        }
    }
}

void Tracer::writeImage(const std::string &fname) {
#ifdef __APPLE__
    std::ofstream out(fname, std::ios::out | std::ios::binary);
#else
    std::ofstream out(fname.c_str(), std::ios_base::binary);
#endif
    out<<"P6"<<"\n"<<res[0]<<" "<<res[1]<<"\n"<<255<<"\n";
    SlVector3 *pixel = im;
    char val;
    for (unsigned int i=0; i<res[0]*res[1]; i++, pixel++) {
        val = (unsigned char)(std::min(1.0, std::max(0.0, (*pixel)[0])) * 255.0);
        out.write (&val, sizeof(unsigned char));
        val = (unsigned char)(std::min(1.0, std::max(0.0, (*pixel)[1])) * 255.0);
        out.write (&val, sizeof(unsigned char));
        val = (unsigned char)(std::min(1.0, std::max(0.0, (*pixel)[2])) * 255.0);
        out.write (&val, sizeof(unsigned char));
    }
    out.close();
}


int main(int argc, char *argv[]) {
    int c;
    double aperture = 0.0;
    int samples = 1;
    int maxraydepth = 5;
    bool color = false;
    while ((c = getopt(argc, argv, "a:s:d:c")) != -1) {
        switch(c) {
            case 'a':
            aperture = atof(optarg);
            break;
            case 's':
            samples = atoi(optarg);
            break;
            case 'c':
            color = true;
            break;
            case 'd':
            maxraydepth = atoi(optarg);
            break;
            default:
            abort();
        }
    }

    if (argc-optind != 2) {
        std::cout<<"usage: trace [opts] input.nff output.ppm"<<std::endl;
        for (unsigned int i=0; i<argc; i++) std::cout<<argv[i]<<std::endl;
        exit(0);
    }	

    Tracer tracer(argv[optind++]);
    tracer.aperture = aperture;
    tracer.samples = samples;
    tracer.color = color;
    tracer.maxraydepth = maxraydepth;
    tracer.traceImage();
    tracer.writeImage(argv[optind++]);
};
