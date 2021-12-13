#version 450

layout (local_size_x = 1, local_size_x = 1) in;
layout (binding=0, rgba8) uniform image2D output_texture;
layout (binding=1, rgba32f) uniform image2D debug_texture;
float camera_pos_z = 5.0;

struct Ray
{
  vec3 origin;
  vec3 dir;
};

struct Light
{
  vec3 position;
};

struct Sphere
{
  vec3 center;
  float radius;
};

struct Collision
{
  float t;              // distance along ray
  vec3 p;               // world position of collision
  vec3 n;               // surface normal at collision
  bool inside;          // ray inside the object or not
  int object_index;     // index of the object
};


// Adapted from scratchapixel's sphere intersection 
Collision collision(Sphere sphere, Ray ray)
{
  Collision c;
  c.inside = false;

  vec3 oc = sphere.center - ray.origin;
  float t_ca = dot(oc, ray.dir);
  float d2 = dot(oc, oc) - t_ca * t_ca;

  // Ray goes outside of sphere
  if (d2 > sphere.radius * sphere.radius) 
  {
    c.t = -1;
    return c;
  }

  float t_hc = sqrt(sphere.radius*sphere.radius - d2);
  float t0 = t_ca - t_hc;
  float t1 = t_ca + t_hc;
  float t_near = min(t0, t1);
  float t_far = max(t0, t1);
  c.t = t_near;

  // Sphere behind the ray, no intersection
  if (t_far < 0.0)
  {
    c.t = -1.0;
    return c;
  }

  // Ray inside the sphere
  if (t_near < 0.0)
  {
    c.t = t_far;
    c.inside = true;
  }

  c.p = ray.origin + ray.dir * t0; 
  c.n = normalize(c.p - sphere.center);

  // Flip normal if ray inside
  if (c.inside) c.n *= -1.0;

  return c;
}

vec3 raytrace(Ray ray)
{
  vec3 ambient = vec3(0.05, 0.05, 0.05);
  Sphere sphere;
  sphere.center = vec3(0, 0, -3);
  sphere.radius = 5.0;

  Light light;
  light.position = vec3(0, 1, 3);

  Collision c = collision(sphere, ray);
  if (c.t > 0.0)
  {
    vec3 light_dir = normalize(light.position - c.p);
    vec3 light_refl = normalize(reflect(light_dir, c.n));
    float cos_theta = dot(light_dir, c.n);
    float cos_phi = dot(normalize(-ray.dir), light_refl);

    float diffuse = 0.9 * max(cos_theta, 0.0);
    float specular = 0.9 * pow(max(cos_phi, 0.0), 30);
    
    return ambient + vec3(1,0,0) * (diffuse);
  }
  
  return ambient;
}


void main()
{
  int width = int(gl_NumWorkGroups.x); // one workgroup = one invocation = one pixel
  int height = int(gl_NumWorkGroups.y);
  ivec2 pixel = ivec2(gl_GlobalInvocationID.xy);
  
  // Convert this pixel's screen space location to world space
  float x_pixel = 2.0 * pixel.x / width - 1.0;
  float y_pixel = 2.0 * pixel.y / height - 1.0;
  vec3 pixel_world = vec3(x_pixel, y_pixel, camera_pos_z - 1.0);
  
  // Get this pixel's world-space ray
  Ray ray;
  ray.origin = vec3(0.0, 0.0, camera_pos_z);
  ray.dir = normalize(pixel_world - ray.origin);
  
  // Cast the ray out into the world and intersect the ray with objects
  vec3 color = raytrace(ray);
  //vec3 color = vec3(float(pixel.x), float(pixel.y), 1);
  imageStore(output_texture, pixel, vec4(color, 1.0));
  imageStore(debug_texture, pixel, vec4(color, 1.0));
}