#version 450

layout (local_size_x = 1, local_size_x = 1) in;
layout (binding=0, rgba8) uniform image2D output_texture;
layout (binding=1, rgba32f) uniform image2D debug_texture;
float camera_pos_z = 5.0;

float EPSILON = 1e-8;

float seed = 10;	//seed initialized in main
float rnd() { return fract(sin(seed++)*43758.5453123); }

struct Ray
{
  vec3 origin;
  vec3 dir;
};

struct Light
{
  vec3 position;
};

struct Material
{
  vec3 diffusivity;
  vec3 emmision;
};


struct Sphere
{
  vec3 center;
  float radius;
  Material mat;
};

struct Triangle
{
  vec3 p0;
  vec3 p1;
  vec3 p2;
  Material mat;
};

struct Collision
{
  float t;              // distance along ray
  vec3 p;               // world position of collision
  vec3 n;               // surface normal at collision
  bool inside;          // ray inside the object or not
  Material mat;
  Ray ray;
  int object_index;     // index of the object
};

struct Sample
{
  vec3 value;
  float pdf;
};

Light light = Light(vec3(0, 1, 3));

Material emissive = Material(
  vec3(1, 1, 1),
  vec3(1, 1, 1)
);

Material red = Material(
  vec3(1, 0, 0),
  vec3(0, 0, 0)
);

Material blue = Material(
  vec3(0, 0, 1),
  vec3(0, 0, 0)
);

Sphere sphere = Sphere(vec3(0, -1.5, -2.5), 1.5, blue);

Triangle triangle0 = Triangle(
  vec3(5, -3, -5),
  vec3(-5, -3, -5),
  vec3(5, -3, 0),
  red
);

Triangle triangle1 = Triangle(
  vec3(-5, -3, 0),
  vec3(5, -3, 0),
  vec3(-5, -3, -5),
  red
);

Triangle triangle2 = Triangle(
  vec3(5, 2, -5),
  vec3(5, 2, 0),
  vec3(-5, 2, -5),
  emissive
);

Triangle triangle3 = Triangle(
  vec3(-5, 2, 0),
  vec3(-5, 2, -5),
  vec3(5, 2, 0),
  emissive
);

Collision collisions[5];

vec3 ambient = vec3(0.05, 0.05, 0.05);

// Adapted from scratchapixel's sphere collision 
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

  // Sphere behind the ray, no collision
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

  c.p = ray.origin + ray.dir * c.t; 
  c.n = normalize(c.p - sphere.center);
  c.mat = sphere.mat;
  c.ray = ray;

  // Flip normal if ray inside
  if (c.inside) c.n *= -1.0;

  return c;
}

// Möller-Trumbore algorithm : 
//   - express problem in barycentric coordinate P = wA + uB + vC
//   - also : P = O + tD
//   - reorganise equation of collision (unknowns should be t, u, v)
//   - solve system using Cramer's rule
Collision collision(Triangle triangle, Ray ray)
{
  vec3 edge1 = triangle.p1 - triangle.p0;
  vec3 edge2 = triangle.p2 - triangle.p0;

  Collision c;
  c.inside = false;

  vec3 h = cross(ray.dir, edge2);
  float a = dot(edge1, h);

  // Parallel ray and triangle
  if (a > -EPSILON && a < EPSILON)
  {
    c.t = -1;
    return c;
  }

  float f = 1.0 / a;

  vec3 s = ray.origin - triangle.p0;
  float u = dot(s, h) * f;
  if (u < 0.0 || u > 1.0)
  {
    c.t = -1;
    return c;
  }

  vec3 q = cross(s, edge1);
  float v = dot(ray.dir, q) * f;
  if (v < 0.0 || u + v > 1.0)
  {
    c.t = -1;
    return c;
  }

  float t = dot(edge2, q) * f;
  if (t < EPSILON)
  {
    c.t = -1;
    return c;
  }

  c.t = t;
  c.p = ray.origin + ray.dir * c.t; 
  c.n = normalize(cross(edge1, edge2));
  c.mat = triangle.mat;
  c.ray = ray;

  return c;
}

Collision find_nearest(Ray ray)
{
  collisions[0] = collision(sphere, ray);
  collisions[1] = collision(triangle0, ray);
  collisions[2] = collision(triangle1, ray);
  collisions[3] = collision(triangle2, ray);
  collisions[4] = collision(triangle3, ray);

  float infinity = 1.0 / 0.0;
  float min_val = infinity;
  int min_i = -1;
  for (int i = 0; i < 5; i++)
  {
    if (collisions[i].t > 0 && min_val > collisions[i].t)
    {
      min_val = collisions[i].t;
      min_i = i;
    }
  }
  min_i = (min_i == -1) ? 0 : min_i;
  return collisions[min_i];
}

vec2 uniform_sample_triangle(vec2 u)
{
  float su0 = sqrt(u.x);
  return vec2(1 - su0, u.y * su0);
}

Sample area_sample(Triangle t, vec3 origin, vec3 n)
{
  vec2 b = uniform_sample_triangle(vec2(rnd(), rnd()));
  vec3 sample_pt =  t.p0 * b.x + t.p1 * b.y + t.p2 * (1 - b.x - b.y); 

  vec3 dir = normalize(sample_pt - origin);

  float area_pdf = dot(n, n) / 2;

  return Sample(dir, area_pdf);
}

// MIS could be added here
vec3 uniform_sample_one_light(Collision c)
{
  float pdf = 0.0;

  // TODO: pick random light
  Triangle light = triangle2;

  // Sample a ray direction from light to collision point
  Sample light_sample = area_sample(light, c.p, c.n);   // FIXME weird to put n as arg
  vec3 wi = light_sample.value;
  
  // FIXME adapt way by comparing normal
  Ray ray_in = Ray(c.p, wi);

  Collision light_col = collision(light, ray_in);
  if (light_col.t <= 0)
  {
    return vec3(0, 0, 0);
  }

  vec3 wo = -c.ray.dir;
  //vec3 f = evaluate_bsdf(wo, wi);
  vec3 f = vec3(1);

  // check for FIX pbr - sampling lights 
  float cos1 = dot(light_col.n, wi);
  float cos2 = dot(c.n, -wi);
  float dist = dot(light_col.p - c.p, light_col.p - c.p);

  return (c.mat.emmision * f * cos1 * cos2) / (dist * dist * light_sample.pdf);
}


vec3 pathtrace(Ray ray)
{
  vec3 L = vec3(0);                    // Total radiance estimate
  vec3 throughput = vec3(1);           // Current path throughput

  int max_bounces = 3;
  bool specular_bounce = false;

  for (int bounces = 0; ; bounces++)
  {
    // Intersect ray with scene
    Collision c = find_nearest(ray);

    // Stop if no collision or no more bounce
    if (c.t <= 0 || bounces > max_bounces) break;

    // Account for the emission if :
    //  - it is the initial collision
    //  - previous was specular BSDF so no direct illumination estimate (Dirac distribution) 
    if (bounces == 0 || specular_bounce)
    {
      if (c.t > 0)
      {
        L += throughput * c.mat.emmision;
      }
      else
      {
        // TODO: infinite area light sources
      }
    }

    // TODO: Compute scattering functions and skip over medium boundaries

    // Direct lighting estimation (end of the current path)
    L += throughput * uniform_sample_one_light(c);

  }

  return vec3(0,0,0);
}

vec3 phong(Collision c, Ray ray)
{
  vec3 light_dir = normalize(light.position - c.p);
  vec3 light_refl = normalize(reflect(light_dir, c.n));
  float cos_theta = dot(light_dir, c.n);
  float cos_phi = dot(normalize(-ray.dir), light_refl);

  float diffuse = 0.9 * max(cos_theta, 0.0);
  float specular = 0.9 * pow(max(cos_phi, 0.0), 30);
  
  return ambient + c.mat.diffusivity * (diffuse);
}

vec3 raytrace(Ray ray)
{
  Collision c = find_nearest(ray);
  if (c.t > 0.0) return phong(c, ray);
  
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