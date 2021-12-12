#version 450

precision highp float;

layout (location = 0) in vec3 in_position;
layout (location = 1) in vec3 in_normal;

uniform mat4 localToProjection;

out vec3 worldPos;
out vec3 vNormalWS;
out vec3 vviewDir;

void main()
{
  vec3 cameraPos = vec3(0.0, 0.0, 4.0);
  vec3 modelPos = vec3(0.0, 0.0, 0.0);

  gl_Position = localToProjection * vec4(in_position, 1.0);

  worldPos = in_position.xyz;
  vNormalWS = normalize(in_normal);
  vviewDir = normalize(cameraPos - in_position.xyz);
}