#include "mesh.hh"

Mesh::Mesh(std::vector<Vertex>& vertices, std::vector<unsigned int>& indices)
: vertices_(vertices),
  indices_(indices),
  VAO_(0),
  VBO_(0),
  EBO_(0),
  instance_(nullptr)
{
}

Mesh::~Mesh()
{
}

void Mesh::attach_shader(program* instance)
{
  instance_ = instance;
}

void Mesh::init_mesh()
{
  if (!instance_)
  {
    std::cerr << "A shader should be attached to mesh before setup" << std::endl; 
    return;
  }
  glGenVertexArrays(1, &VAO_); TEST_OPENGL_ERROR();
  glGenBuffers(1, &VBO_); TEST_OPENGL_ERROR();
  glGenBuffers(1, &EBO_); TEST_OPENGL_ERROR();

  glBindVertexArray(VAO_);
  bind_buffer(vertices_, VBO_, GL_ARRAY_BUFFER);
  bind_buffer(indices_, EBO_, GL_ELEMENT_ARRAY_BUFFER);

  //We'll optimize it later
  std::size_t sizes[] = {3, 3, 2};
  for (std::size_t i = 0; i < 3; ++i)
    attrib_array(i, sizes[i]);

  glBindVertexArray(0);
}

void Mesh::render() const
{
  instance_->use();

  glBindVertexArray(VAO_); TEST_OPENGL_ERROR();
  glDrawElements(GL_TRIANGLES, indices_.size(), GL_UNSIGNED_INT, 0); TEST_OPENGL_ERROR();
  glBindVertexArray(0); TEST_OPENGL_ERROR();

}

//nb is the number of elements per coordinate
template <typename T>
void Mesh::bind_buffer(std::vector<T>& buffer, GLuint index, GLenum type)
{
  //Generate buffer
  glBindBuffer(type, index); TEST_OPENGL_ERROR();
  glBufferData(type, buffer.size() * sizeof(T), &buffer[0], GL_STATIC_DRAW); TEST_OPENGL_ERROR();
}

void Mesh::attrib_array(std::size_t idx, std::size_t size)
{
  glEnableVertexAttribArray(idx); TEST_OPENGL_ERROR();
  glVertexAttribPointer(idx, size, GL_FLOAT, GL_FALSE, sizeof(Vertex),
    (void*) (sizeof(float) * idx)); TEST_OPENGL_ERROR();
}