#version 450

precision highp float;

layout (location = 0) in vec3 in_position;
layout (location = 1) in vec3 in_normal;
layout (location = 2) in vec2 in_uv;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;
uniform float time;
//out vec3 worldPos;
//out vec3 vNormalWS;
//out vec3 vviewDir;

void main()
{
  vec3 cameraPos = vec3(0.0, 0.0, 4.0);
  vec3 modelPos = vec3(0.0, 0.0, 0.0);

  vec4 p = model * vec4(in_position, 1.0);

  gl_Position = projection * view * p;;

  //worldPos = in_position.xyz;
  //vNormalWS = normalize(in_normal);
  //vviewDir = normalize(cameraPos - in_position.xyz);
}