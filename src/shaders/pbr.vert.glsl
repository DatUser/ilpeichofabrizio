#version 450

#define M_PI 3.1415926535897932384626433832795 

precision highp float;

layout (location = 0) in vec3 in_position;
layout (location = 1) in vec3 in_normal;
layout (location = 2) in vec2 in_uv;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;
uniform float time;

out vec3 worldPos;//position of pixel in world
out vec3 vNormalWS;
out vec2 uv;

//Here sphere is considered to be at (0,0,0)
vec2 to_spherical(vec3 p)
{
  float sphere_x = tan(p.x) / M_PI;
  float val = sqrt(p.x * p.x + p.y * p.y);
  float sphere_y = atan(val / p.z) ;// (M_PI);

  return vec2(sphere_x, sphere_y);
}

void
main()
{
  vec3 cameraPos = vec3(0.0, 0.0, 4.0);
  vec3 modelPos = vec3(0.0, 0.0, 0.0);
  vec4 p = model * vec4(in_position, 1.0);

  gl_Position = projection * view * p;

  worldPos = p.xyz;
  vNormalWS = normalize(in_normal);
  uv = in_uv;
}
