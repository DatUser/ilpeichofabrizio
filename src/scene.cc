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
  // Materials
  std::vector<Material> materials;
  for (auto const& j_mat : j.at("materials"))
  {
    auto kd = j_mat.at("Kd").get<std::array<float, 3>>();
    auto ke = j_mat.at("Ke").get<std::array<float, 3>>();
    Material mat = {
      glm::vec3(kd[0], kd[1], kd[2]),
      glm::vec3(ke[0], ke[1], ke[2])
    };
    materials.push_back(mat);
  }                

  std::cout << "Materials: " << materials.size() << std::endl;

  // Triangles
  std::vector<Triangle> triangles;
  for (auto const& j_tri : j.at("triangles"))
  {
    auto p1 = j_tri.at("p1").get<std::array<float, 3>>();
    auto p2 = j_tri.at("p2").get<std::array<float, 3>>();
    auto p3 = j_tri.at("p3").get<std::array<float, 3>>();
    auto mat_id = j_tri.at("matId").get<uint8_t>();
    Triangle t = {
      glm::vec3(p1[0], p1[1], p1[2]),
      glm::vec3(p2[0], p2[1], p2[2]),
      glm::vec3(p3[0], p3[1], p3[2]),
      mat_id
    };
    triangles.push_back(t);
  }                

  std::cout << "Triangles: " << triangles.size() << std::endl;

  // Obj files
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
  std::vector<tinyobj::material_t> materials;

  std::string err;
  
  const char* mtl_basedir_cstr = (mtl_basedir.empty()) 
    ? NULL
    : mtl_basedir.c_str();
  
  bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &err,
                              obj_file.c_str(), mtl_basedir_cstr);

  if (!err.empty()) std::cerr << err << std::endl;

  if (!ret) std::cerr << "Model file could not be parsed !" << std::endl;
  else std::cout << "Model successfully loaded !" << std::endl;

  std::vector<glm::vec3> vertices;

  // Here we load all meshes
  for (const auto& shape : shapes)
  {
    size_t idx_offset = 0;

    // Iterate through triangles
    for (size_t face_v : shape.mesh.num_face_vertices)
    {
      // Iterate inside of triangles
      for (size_t v = 0; v < face_v; ++v)
      {
        size_t idx = shape.mesh.indices[idx_offset + v].vertex_index;

        glm::vec3 vertex = glm::vec3(
                  attrib.vertices[3 * idx],
                  attrib.vertices[3 * idx + 1],
                  attrib.vertices[3 * idx + 2]
        );
        vertices.push_back(vertex);
      }
      idx_offset += face_v;
    }
  }

  std::cout << "Loaded " << vertices.size()<< " vertices" << std::endl;

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
