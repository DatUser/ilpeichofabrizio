#version 450

layout (local_size_x = 1, local_size_x = 1) in;
layout (binding=0, rgba8) uniform image2D output_texture;
float camera_pos_z = 5.0;

struct Ray
{
  vec3 origin;
  vec3 dir;
};


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
  //vec3 color = raytrace(world_ray);
  vec3 color = vec3(float(pixel.x) / 64, float(pixel.y) / 64, 1);
  imageStore(output_texture, pixel, vec4(color, 1.0));
}