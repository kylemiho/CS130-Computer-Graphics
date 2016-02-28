/**
 * ray_tracer.cpp
 * CS230
 * -------------------------------
 * Implement ray tracer here.
 */

#define SET_RED(P, C)   (P = (((P) & 0x00ffffff) | ((C) << 24)))
#define SET_GREEN(P, C)  (P = (((P) & 0xff00ffff) | ((C) << 16)))
#define SET_BLUE(P, C) (P = (((P) & 0xffff00ff) | ((C) << 8)))

#include "ray_tracer.h"
#include <cmath>

using namespace std;

const double Object::small_t=1e-6;
//--------------------------------------------------------------------------------
// utility functions
//--------------------------------------------------------------------------------
double sqr(const double x)
{
    return x*x;
}

Pixel Pixel_Color(const Vector_3D<double>& color)
{
    Pixel pixel=0;
    SET_RED(pixel,(unsigned char)(min(color.x,1.0)*255));
    SET_GREEN(pixel,(unsigned char)(min(color.y,1.0)*255));
    SET_BLUE(pixel,(unsigned char)(min(color.z,1.0)*255));
    return pixel;
}
//--------------------------------------------------------------------------------
// Shader
//--------------------------------------------------------------------------------
Vector_3D<double> Phong_Shader::
Shade_Surface(const Ray& ray,const Object& intersection_object,const Vector_3D<double>& intersection_point,const Vector_3D<double>& same_side_normal) const
{
    Vector_3D<double> color;
    Vector_3D<double> ambientCoefficient = color_ambient;
    Vector_3D<double> diffuseReflectance = color_diffuse;
    // TODO: determine the color
    //compute each R G B value seperately
    
    //ambient
    ambientCoefficient.x *= world.lights.at(0)->Emitted_Light(ray).x;
    ambientCoefficient.y *= world.lights.at(0)->Emitted_Light(ray).y;
    ambientCoefficient.z *= world.lights.at(0)->Emitted_Light(ray).z;
    
    color += ambientCoefficient;

    for (unsigned i = 0; i < world.lights.size(); ++i)
    {
        //test for shadows
        Vector_3D<double> lightVector = (((world.lights.at(i))->position)-intersection_point);
        //tempray
        Ray tempRay;
        tempRay.endpoint = intersection_point;
        tempRay.direction = lightVector;
        tempRay.direction.Normalize();
        tempRay.t_max = 999999;
        tempRay.semi_infinite = true;
        bool shadowDetected = false;
        for (unsigned j = 0; j < world.objects.size(); ++j)
        {
            ////calculate if object intersects
            if (world.objects.at(j)->Intersection(tempRay))
            {
                if (!(tempRay.current_object == &intersection_object))
                {
                    shadowDetected = true;
                }
            }
        }

        if (!(world.enable_shadows && shadowDetected))
        {
            Vector_3D<double> normal = intersection_object.Normal(intersection_point);
            //diffuse shading
            Vector_3D<double> diffuseIntensity = world.lights.at(i)->Emitted_Light(ray);
            
            lightVector.Normalize();
            double dotProduct = lightVector.x*normal.x + lightVector.y*normal.y + lightVector.z*normal.z;

            if (dotProduct < 0)
            {
                dotProduct = 0;
            }
           // cout << "RED DIFFUSE INTENSITY: " << diffuseReflectance.x*(diffuseIntensity) << endl;
            color.x += dotProduct*diffuseReflectance.x*(diffuseIntensity.x);
            color.y += dotProduct*diffuseReflectance.y*(diffuseIntensity.y);
            color.z += dotProduct*diffuseReflectance.z*(diffuseIntensity.z);
            //specular shading
         Vector_3D<double> eyeVector = ray.endpoint - intersection_point;
            
            eyeVector.Normalize();
            
            Vector_3D<double> tempH = lightVector+eyeVector;
            tempH.Normalize();
 
            double max = (tempH.x*same_side_normal.x+tempH.y*same_side_normal.y+tempH.z*same_side_normal.z);
            if (max < 0)
            {
                max = 0;
            }
             
            color.x += color_specular.x*(world.lights.at(i)->Emitted_Light(ray).x)*pow(max,specular_power);
            color.y += color_specular.y*(world.lights.at(i)->Emitted_Light(ray).y)*pow(max,specular_power);
            color.z += color_specular.z*(world.lights.at(i)->Emitted_Light(ray).z)*pow(max,specular_power);
        }
    }   

    return color;
}

Vector_3D<double> Reflective_Shader::
Shade_Surface(const Ray& ray,const Object& intersection_object,const Vector_3D<double>& intersection_point,const Vector_3D<double>& same_side_normal) const
{

    Vector_3D<double> color;
    Vector_3D<double> color2;
    Vector_3D<double> ambientCoefficient = color_ambient;
    Vector_3D<double> diffuseReflectance = color_diffuse;
    // TODO: determine the color
    //compute each R G B value seperately
    
    //ambient
    if (ray.recursion_depth == 0)
    {
        ambientCoefficient.x *= world.lights.at(0)->Emitted_Light(ray).x;
        ambientCoefficient.y *= world.lights.at(0)->Emitted_Light(ray).y;
        ambientCoefficient.z *= world.lights.at(0)->Emitted_Light(ray).z;
        color += ambientCoefficient;
    }
    


    for (unsigned i = 0; i < world.lights.size(); ++i)
    {
        //test for shadows
        Vector_3D<double> lightVector = (((world.lights.at(i))->position)-intersection_point);
        //tempray
        Ray tempRay;
        tempRay.endpoint = intersection_point;
        tempRay.direction = lightVector;
        tempRay.direction.Normalize();
        tempRay.t_max = 999999;
        tempRay.semi_infinite = true;
        bool shadowDetected = false;
        for (unsigned j = 0; j < world.objects.size(); ++j)
        {
            ////calculate if object intersects
            if (world.objects.at(j)->Intersection(tempRay))
            {
                if (!(tempRay.current_object == &intersection_object))
                {
                    shadowDetected = true;
                }
            }
        }
            Vector_3D<double> normal = intersection_object.Normal(intersection_point);
        if (!(world.enable_shadows && shadowDetected))
        {

            //diffuse shading
            Vector_3D<double> diffuseIntensity = world.lights.at(i)->Emitted_Light(ray);
            
            lightVector.Normalize();
            double dotProduct = lightVector.x*normal.x + lightVector.y*normal.y + lightVector.z*normal.z;

            if (dotProduct < 0)
            {
                dotProduct = 0;
            }
           // cout << "RED DIFFUSE INTENSITY: " << diffuseReflectance.x*(diffuseIntensity) << endl;
            color.x += dotProduct*diffuseReflectance.x*(diffuseIntensity.x);
            color.y += dotProduct*diffuseReflectance.y*(diffuseIntensity.y);
            color.z += dotProduct*diffuseReflectance.z*(diffuseIntensity.z);
            //specular shading
         Vector_3D<double> eyeVector = ray.endpoint - intersection_point;
            
            eyeVector.Normalize();
            
            Vector_3D<double> tempH = lightVector+eyeVector;
            tempH.Normalize();
 
            double max = (tempH.x*same_side_normal.x+tempH.y*same_side_normal.y+tempH.z*same_side_normal.z);
            if (max < 0)
            {
                max = 0;
            }
             
            color.x += color_specular.x*(world.lights.at(i)->Emitted_Light(ray).x)*pow(max,specular_power);
            color.y += color_specular.y*(world.lights.at(i)->Emitted_Light(ray).y)*pow(max,specular_power);
            color.z += color_specular.z*(world.lights.at(i)->Emitted_Light(ray).z)*pow(max,specular_power);
            Vector_3D<double> reflection = ray.direction - normal * (2*(lightVector.x*normal.x + lightVector.y*normal.y + lightVector.z*normal.z));
            //(normal*(2*(lightVector.x*normal.x + lightVector.y*normal.y + lightVector.z*normal.z))) - lightVector;   
            reflection.Normalize();
            
            Ray reflectionRay;
            reflectionRay.endpoint = intersection_point;
            reflectionRay.direction = reflection;
        
            reflectionRay.t_max = 999999;
            reflectionRay.semi_infinite = true;
            reflectionRay.recursion_depth = ray.recursion_depth + 1;
            if (reflectionRay.recursion_depth < 5)
            {
                //cout << reflectionRay.recursion_depth << endl;

                color2 = world.Cast_Ray(reflectionRay,ray)*color_specular;//*(world.lights.at(i)->Emitted_Light(ray))*pow(max,specular_power);

                if ((color2.x > 0.2) || (color2.y > 0.2) || (color.z > 0.2))
                {
                   // color.y = 1;

                    cout << color.x << " " << color.y << " " << color.z << endl;
                }
                if (color.x > 1)
                    color.x = 1;
                if (color.y > 1)
                    color.y = 1;
                if (color.z > 1)
                    color.z = 1;
                         //           color.y = 1;

            }
            
        }
    }
        color += color2;
        return color;
}

Vector_3D<double> Flat_Shader::
Shade_Surface(const Ray& ray,const Object& intersection_object,const Vector_3D<double>& intersection_point,const Vector_3D<double>& same_side_normal) const
{
    return color;
}

//--------------------------------------------------------------------------------
// Objects
//--------------------------------------------------------------------------------
// determine if the ray intersects with the sphere
// if there is an intersection, set t_max, current_object, and semi_infinite as appropriate and return true

//http://www.lighthouse3d.com/tutorials/maths/ray-sphere-intersection/
//http://kylehalladay.com/blog/tutorial/math/2013/12/24/Ray-Sphere-Intersection.html
bool Sphere::
Intersection(Ray& ray) const
{
   // cout << "Sphere: " << this << endl;
    //cout << "Center- x: " << center.x << " y: " << center.y << " z: " << center.z << endl;
    //return false;
    //cout << "Ray Direction- x: " << ray.direction.x << " y: " << ray.direction.y << " z: " << ray.direction.z << endl;
    //cout << "Ray Endpoint- x: " << ray.endpoint.x << " y: " << ray.endpoint.y << " z: " << ray.endpoint.z << endl;
    //                                             vector d = starting point of ray -> end point of ray
    Vector_3D<double> v = center - ray.endpoint; //vector v =  starting point of ray --> center of sphere
    //v.Normalize();
                                              
    //calculate v dot u (dot product) to determine if sphere is in front or behind of ray
    double temp = v.x* ray.direction.x + v.y * ray.direction.y +v.z *ray.direction.z;
   // cout << "ray->center dot ray direction: " << temp << endl;
    
    //if (temp > 0 )  //sphere is along ray -- http://www.lighthouse3d.com/tutorials/maths/line-and-rays/
    //{
        //calculate proj, or closest point of ray to sphere
        //calculate distance1, or the distance from proj to sphere center
        //if distance1 < radius, 2 intersections, if distance 1 == radius, 1 intersection
        //if distance1 > radius, 0 intersections
        //distance1 = abs(spherecenter - proj)
        //q is sphere center
        
        Vector_3D<double> q = center;
        // u dot v * u / magnitude u
        //proj = p + puv |or| ray.endpoint + ray.endpoint*v*ray.direction
     //   cout << "V magnitude: " << v.Length_Squared() << endl;
     //  temp = temp / ray.direction.Length_Squared();
        
        Vector_3D<double> qstar =  ray.direction* temp;
        qstar = qstar + ray.endpoint;
        //can find distance between q and qstar
        Vector_3D<double> distance1 = q - qstar;
        double distance2 = distance1.Length();
        //calculate t
        double b = sqrt(pow((center.x - qstar.x),2) + pow((center.y - qstar.y),2) + pow((center.z - qstar.z),2));
        double t =  (radius*radius) - (b*b);
       // | pc - p | - c
        double t2 = sqrt(pow((ray.endpoint.x - qstar.x),2) + pow((ray.endpoint.y - qstar.y),2) + pow((ray.endpoint.z - qstar.z),2));
        t2 = t2 - t;
        if ((distance2 <= radius) && (t2 < ray.t_max)) //if true, then intersection exists also check if it is closest intersection
        {
            ray.semi_infinite = false;
            ray.current_object = this;
            //tmax is length from ray-endpoint to intersection
            ray.t_max = t2;
       //    cout << "Distance: " << distance2 << "  Radius: " << radius << endl;
       //          cout << "Q- x: " << q.x << " y: " << q.y << " z: " << q.z << endl;
       //          cout << "Q*- x: " << qstar.x << " y: " << qstar.y << " z: " << qstar.z << endl;
       //    cout << "INTERSECTED" << endl << endl;
            return true;
        }
       //          cout << "Q- x: " << q.x << " y: " << q.y << " z: " << q.z << endl;
        //         cout << "Q*- x: " << qstar.x << " y: " << qstar.y << " z: " << qstar.z << endl;
       //          cout << "Distance: " << distance2 << " Radius: " << radius << endl << endl;
        
   // } else          //sphere is opposite side of ray
  //  {
        //DO LATER
        //http://www.lighthouse3d.com/tutorials/maths/ray-sphere-intersection/
  //  }
  

    return false;
}

Vector_3D<double> Sphere::
Normal(const Vector_3D<double>& location) const
{
    Vector_3D<double> normal;
    normal = (location - center);
    normal.Normalize();
    
    return normal;
}

// determine if the ray intersects with the sphere
// if there is an intersection, set t_max, current_object, and semi_infinite as appropriate and return true
bool Plane::
Intersection(Ray& ray) const
{

    ////if the dot product of the ray direction and the normal is <90, then the ray is guaranteed to intersect the plane
    Vector_3D<double> v = x1 - ray.endpoint; //v = vector from ray endpoint -> "center of plane"
    
    double topProduct = normal.x * v.x + normal.y * v.y + normal.z * v.z;
    double bottomProduct = normal.x * ray.direction.x + normal.y * ray.direction.y + normal.z * ray.direction.z;
    
    double t = topProduct/bottomProduct;
    if ((t > 0) && (t < ray.t_max)) //check if it intersects and if it is closest point
    {
    
    //t = n_dot_(x1-e)/(n_dot_d)
    ray.t_max = t;
    ray.semi_infinite = false;
    ray.current_object = this;
    return true;
    } else
    {
        return false;
    }
}

Vector_3D<double> Plane::
Normal(const Vector_3D<double>& location) const
{
    return normal;
}
//--------------------------------------------------------------------------------
// Camera
//--------------------------------------------------------------------------------
// Find the world position of the input pixel
Vector_3D<double> Camera::
World_Position(const Vector_2D<int>& pixel_index)
{
 // cout << "Index- x: " << pixel_index.x << " y: " << pixel_index.y << endl;
 //   cout << "PIXELGRIDCENTER- x:" << focal_point.x << " y: " << focal_point.y <<  "z: " << focal_point.z << endl;
    Vector_3D<double> result;
    int minPixelSize = 0;
    if(film.pixel_grid.m < film.pixel_grid.n)
    {
        minPixelSize = film.pixel_grid.m;
    } else
    {
        minPixelSize = film.pixel_grid.n;
    }
    result = (horizontal_vector*((pixel_index.x-((film.pixel_grid.m-1)*0.5))/minPixelSize) + vertical_vector*((pixel_index.y-((film.pixel_grid.n-1)*0.5))/minPixelSize))+focal_point;
 // cout << "GridPoint- x: " << result.x << " y: " << result.y <<  "z: " << result.z << endl;
    return result;
}
//--------------------------------------------------------------------------------
// Render_World
//--------------------------------------------------------------------------------
// Find the closest object of intersection and return a pointer to it
//   if the ray intersects with an object, then ray.t_max, ray.current_object, and ray.semi_infinite will be set appropriately
//   if there is no intersection do not modify the ray and return 0
const Object* Render_World::
Closest_Intersection(Ray& ray)
{
    Object* tempObject = 0;
    for (unsigned i = 0; i < objects.size(); ++i)
    {
        //calculate if object intersects
        if (objects.at(i)->Intersection(ray))
        {
            tempObject = objects.at(i);
        }
    }
    return tempObject;
}

// set up the initial view ray and call 
void Render_World::
Render_Pixel(const Vector_2D<int>& pixel_index)
{
    // TODO
    Ray ray; // TODO: set up the initial view ray here
    
    ray.endpoint = camera.position; //beginning of ray is start of camera
    Vector_3D<double> pixelGrid = camera.World_Position(pixel_index); //
    ray.direction = pixelGrid - ray.endpoint;
    ray.direction.Normalize();
    ray.t_max = 999999;
    ray.recursion_depth = 0;
    
    Ray dummy_root;
    Vector_3D<double> color=Cast_Ray(ray,dummy_root);
    camera.film.Set_Pixel(pixel_index,Pixel_Color(color));
}

// cast ray and return the color of the closest intersected surface point, 
// or the background color if there is no object intersection
Vector_3D<double> Render_World::
Cast_Ray(Ray& ray,const Ray& parent_ray)
{
    // TODO
    Vector_3D<double> color;

    if (Closest_Intersection(ray)) //intersection occurs
    {
        Vector_3D<double> intersect = ray.endpoint + ray.direction * ray.t_max;
        Vector_3D<double> normal = ray.current_object->Normal(intersect);
        //normal.Normalize();

        color = (ray.current_object)->material_shader->Shade_Surface(ray,*(ray.current_object),intersect,normal);
       // cout << "intersection object:" << ray.current_object << endl;

    }
    

    return color;
}
