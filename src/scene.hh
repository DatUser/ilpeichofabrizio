#pragma once
#define GLM_FORCE_SWIZZLE

#include <string>
#include <vector>
#include <glm/glm.hpp>
#include <glm/vec3.hpp>
#include <glm/vec4.hpp>
#include <nlohmann/json.hpp>


using json = nlohmann::json;
using Vertex = glm::vec4;
using Vertex3 = glm::vec3;

#pragma align(32)
struct BVHNode {
  Vertex3 box_min; // 16 byte aligned -> 12 = (3 * 4)
  union {
    int firstChildNodeID; // for inner nodes
    int firstTriangleID; // for leaf nodes
  }; //-> 16 = (12 + 4)
  //--- 16 ---
  Vertex3 box_max; // 16 byte aligned -> 12 = (3 * 4)
  int num_triangles; // 0 flags inner node -> 16 = (12 + 4)
};

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

struct Box
{
  Vertex3 bmin;
  Vertex3 bmax;
};

struct BVHTriangle
{
  BVHTriangle(Vertex3 verts[], Box& bounds, Triangle triangle)
  {
    this->verts[0] = verts[0];
    this->verts[1] = verts[1];
    this->verts[2] = verts[2];
    this->bounds = bounds;
    this->triangle = triangle;
  }

  Vertex3 verts[3];
  Box bounds;
  Triangle triangle;
  //Other stuff
};

struct
{
  bool operator()(BVHTriangle& a, BVHTriangle& b) const
  { return a.bounds.bmin[0] < b.bounds.bmin[0]; }
} vertexX; 

struct
{
  bool operator()(BVHTriangle& a, BVHTriangle& b) const
  { return a.bounds.bmin[1] < b.bounds.bmin[1]; }
} vertexY;

struct
{
  bool operator()(BVHTriangle& a, BVHTriangle& b) const
  { return a.bounds.bmin[2] < b.bounds.bmin[2]; }
} vertexZ;

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
    std::pair<std::vector<BVHTriangle>, Box> parse_scene(json j);
    std::pair<std::vector<BVHTriangle>, Box>
    load_model(const std::string& obj_file, const std::string& mtl_basedir,
                    float scale, const glm::vec3& translation);
    void build_bvh(std::vector<BVHTriangle>& tris, Box& bounds, int begin, int end,
                  std::vector<BVHNode>& tree, int idx);

    std::vector<Vertex> vertices_;
    std::vector<Material> materials_;  
    std::vector<Triangle> triangles_;
    std::vector<Triangle> lights_;
    std::vector<BVHNode> bvh_;
};
