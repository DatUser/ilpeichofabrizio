#pragma once

#include <string>
#include <vector>
#include <glm/vec3.hpp>
#include <glm/vec4.hpp>
#include <nlohmann/json.hpp>

using json = nlohmann::json;
using Vertex = glm::vec4;

struct Material
{
  glm::vec4 kd;
  glm::vec4 ke;
  glm::vec4 ks;
  glm::vec4 kt;
};

struct Triangle
{
  glm::vec3 vertices_index;
  int mat_id;
};

struct BVHNode
{
  glm::vec3 box_min;
  uint8_t left_child;
  glm::vec3 box_max;
  uint8_t count;
};

class Scene
{
  public:
    Scene(const std::string& path);

    const std::vector<Vertex>& get_vertices() { return vertices_; };
    const std::vector<Material>& get_materials() { return materials_; };
    const std::vector<Triangle>& get_triangles() { return triangles_; };
    const std::vector<Triangle>& get_lights() { return lights_ ; }
    const std::vector<BVHNode>& get_bvh() { return bvh_; };

  private:
    void parse_scene(json j);
    void load_model(const std::string& obj_file, const std::string& mtl_basedir,
                    float scale, const glm::vec3& translation);
    void build_bvh();

    std::vector<Vertex> vertices_;
    std::vector<Material> materials_;  
    std::vector<Triangle> triangles_;
    std::vector<Triangle> lights_;
    std::vector<BVHNode> bvh_;
};
