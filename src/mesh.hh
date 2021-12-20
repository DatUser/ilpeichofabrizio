#pragma once

#include <vector>
#include <iostream>
#include "program.hh"

struct Vertex
{
  glm::vec3 pos;
  glm::vec3 normal;
  glm::vec2 uvs;
};

class Mesh
{
public:
  Mesh(std::vector<Vertex>& vertices, std::vector<unsigned int>& indices);
  ~Mesh();

  void attach_shader(program* instance);
  //Shader needs to be attached before mesh init
  void init_mesh();
  void render() const;

private:
  std::vector<Vertex> vertices_;
  std::vector<unsigned int> indices_;

  //VAO
  GLuint VAO_;
  //VBOs w/ number of elements
  GLuint VBO_; 
  //EBO
  GLuint EBO_;

  //Shader
  program* instance_;

  template <typename T>
  void bind_buffer(std::vector<T>& buffer, GLuint index, GLenum type);
  void attrib_array(std::size_t idx, std::size_t size);
  
};
