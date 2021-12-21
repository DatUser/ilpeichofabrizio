#include "objloader.hh"

#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"

std::vector<Mesh*> load_obj(std::string& path)
{
  /*
  tinyobj::ObjReader reader;
  tinyobj::ObjReaderConfig reader_config;
  reader_config.mtl_search_path = "../resources/";

  if (!reader.ParseFromFile(path, reader_config))
  {
    if (!reader.Error().empty()) std::cerr << "TinyObjReader: " << reader.Error();

    return nullptr;
  }

  if (!reader.Warning().empty()) std::cout << "TinyObjReader: " << reader.Warning();
  */

  tinyobj::attrib_t attrib;
  std::vector<tinyobj::shape_t> shapes;
  std::vector<tinyobj::material_t> materials;

  //std::string warn;
  std::string err;

  bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &err, path.c_str());

  std::vector<Mesh*> meshes{};

  if (!err.empty()) {
    std::cerr << err << std::endl;
  }

  if (!ret) std::cerr << path << ": Could not be parsed !" << std::endl;
  else std::cout << path << ": Successfully loaded !" << std::endl;

  //Here we load all meshes
  for (auto shape : shapes)
  {
    size_t idx_offset = 0;
    std::vector<Vertex> vertices{};
    std::vector<unsigned int> indices{};
    //iterate through triangles
    unsigned int indice = 0;
    for (size_t face_v : shape.mesh.num_face_vertices)
    {
      //iterate inside of triangles
      for (size_t v = 0; v < face_v; ++v)
      {
        unsigned int idx = shape.mesh.indices[idx_offset + v].vertex_index;

        glm::vec3 pos = glm::vec3(
                  attrib.vertices[3 * idx],
                  attrib.vertices[3 * idx + 1],
                  attrib.vertices[3 * idx + 2]
                  );
        //glm::vec3* norm = attrib.normals;
        //glm::vec2* uvs = attrib.texcoords;

        Vertex vertex = { pos, glm::vec3(0.0), glm::vec2(0.0) };
        vertices.push_back(vertex);
        indices.push_back(indice);
      }
      indice += 1;
    }

    Mesh* mesh = new Mesh(vertices, indices);
    meshes.push_back(mesh);
  }

  std::cout << "Shape size: " << shapes.size() << std::endl;
  std::cout << "Face size: " << shapes[0].mesh.num_face_vertices.size() << std::endl;
  std::cout << "Pos size: " << attrib.vertices.size() << std::endl;
  std::cout << "Norm size: " << attrib.normals.size() << std::endl;
  std::cout << "Tex size: " << attrib.texcoords.size() << std::endl;

  return meshes;
}