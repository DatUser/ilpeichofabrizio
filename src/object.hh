#pragma once

#include <vector>
#include <glm/vec3.hpp>

#include "program.hh"

class Object
{
public:
  Object(glm::vec3 position);
  ~Object();

  static Object* load_object();

  static void init_object();
  static void attach_shader(program* program);
  static void render();

private:
  glm::vec3 position;
  std::vector<float> vertices;
  std::vector<float> normals;
  std::vector<float> texCoords;

  std::vector<program*> shaders;
};