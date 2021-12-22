#include "scene.hh"

#include <fstream>
#include <iostream>

#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"

Scene::Scene(const std::string& path)
{
  std::ifstream file(path);
  if (!file) 
  {
    std::cerr << "Error while opening file: " << path << std::endl;
    exit(1);
  }

  json j = json::parse(file);
  parse_scene(j);

  build_bvh();
}

void Scene::parse_scene(json j)
{
  for (auto const& j_mod : j.at("models"))
  {
    auto obj_file = j_mod.at("objFile").get<std::string>();
    auto mtl_basedir = j_mod.at("mtlBasedir").get<std::string>();
    auto scale = j_mod.at("scale").get<float>();
    auto t = j_mod.at("translation").get<std::array<float, 3>>();
    glm::vec3 translation = glm::vec3(t[0], t[1], t[2]);

    // Read obj file and add triangles
    load_model(obj_file, mtl_basedir, scale, translation);
  }
}

void Scene::load_model(const std::string& obj_file, const std::string& mtl_basedir,
                       float scale, const glm::vec3& translation)
{
  tinyobj::attrib_t attrib;
  std::vector<tinyobj::shape_t> shapes;
  std::vector<tinyobj::material_t> mtls;

  std::string err;
  
  const char* mtl_basedir_cstr = (mtl_basedir.empty()) 
    ? NULL
    : mtl_basedir.c_str();
  
  bool ret = tinyobj::LoadObj(&attrib, &shapes, &mtls, &err,
                              obj_file.c_str(), mtl_basedir_cstr);

  if (!err.empty()) std::cerr << err << std::endl;

  if (!ret) std::cerr << "Model file could not be parsed !" << std::endl;
  else std::cout << "Model successfully loaded !" << std::endl;

  // Materials
  for (const auto& mtl : mtls) {
    materials_.push_back({
      glm::vec4(mtl.diffuse[0], mtl.diffuse[1], mtl.diffuse[2], 1.0),    // kd
      glm::vec4(mtl.emission[0], mtl.emission[1], mtl.emission[2], 1.0)  // ke
    });
  }

  std::cout << "Material parsed" << std::endl;

  // Vertices
  for (size_t i = 0; i < attrib.vertices.size(); i++)
  {
    Vertex vertex(
      attrib.vertices[3 * i],
      attrib.vertices[3 * i + 1],
      attrib.vertices[3 * i + 2],
      0.f
    );
    vertices_.push_back(vertex);
  }

  std::cout << "Vertices parsed" << std::endl;

  // Triangles
  for (const auto& shape : shapes)
  {
    // Iterate through triangles
    size_t idx_offset = 0;
    size_t i_face = 0;
    for (size_t face_v : shape.mesh.num_face_vertices)
    {
      Triangle triangle;
      triangle.mat_id = shape.mesh.material_ids[i_face];

      // Iterate inside of triangles
      for (size_t v = 0; v < face_v; ++v)
      {
        size_t idx = shape.mesh.indices[idx_offset + v].vertex_index;
        triangle.vertices_index[v] = idx; 
      }

      if (triangle.mat_id != -1 && materials_[triangle.mat_id].ke != glm::vec4(0.f))
      {
        lights_.push_back(triangle);
      }

      triangles_.push_back(triangle);
      idx_offset += face_v;
      i_face++;
    }
  }

  std::cout << "Triangles parsed" << std::endl;

  materials_.push_back({glm::vec4(0.5, 0.5, 0.5, 1), glm::vec4(0, 0, 0, 0)});
  materials_.push_back({glm::vec4(1, 1, 1, 1), 20.f*glm::vec4(1, 1, 1, 1)});

  std::cout << "Loaded " << vertices_.size()<< " vertices" << std::endl;

  std::cout << "Shape size: " << shapes.size() << std::endl;
  std::cout << "Face size: " << shapes[0].mesh.num_face_vertices.size() << std::endl;
  std::cout << "Indices size: " << shapes[0].mesh.indices.size() << std::endl;
  std::cout << "Pos size: " << attrib.vertices.size() << std::endl;
  std::cout << "Norm size: " << attrib.normals.size() << std::endl;
  std::cout << "Tex size: " << attrib.texcoords.size() << std::endl;
}



void Scene::build_bvh()
{

}
