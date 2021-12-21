#pragma once

#include <string>
#include <vector>
#include <glm/vec3.hpp>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

struct Triangle
{
  glm::vec3 p1;
  glm::vec3 p2;
  glm::vec3 p3;
  uint8_t mat_id;
};

struct Material
{
  glm::vec3 kd;
  glm::vec3 ke;
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

    const std::vector<Triangle>& get_triangles() { return triangles_; };
    const std::vector<Material>& get_materials() { return materials_; };
    const std::vector<BVHNode>& get_bvh() { return bvh_; };

  private:
    void parse_scene(json j);
    void load_model(const std::string& obj_file, const std::string& mtl_basedir,
                    float scale, const glm::vec3& translation);
    void build_bvh();

    std::vector<Triangle> triangles_;
    std::vector<Material> materials_;  
    std::vector<BVHNode> bvh_;  
};
