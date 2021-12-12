#pragma once

#include <vector>
#include "program.hh"
#include "vector3.hh"

class Object
{
public:
  Object(Vector3 position);
  ~Object();

  static Object* load_object();

  static void init_object();
  static void attach_shader(program* program);
  static void render();

private:
  Vector3 position;
  std::vector<float> vertices;
  std::vector<float> normals;
  std::vector<float> texCoords;

  std::vector<program*> shaders;
};