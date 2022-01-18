#include "scene.hh"

#include <fstream>
#include <iostream>

#include "tiny_obj_loader.h"

std::vector<Triangle> BtriToTri(std::vector<BVHTriangle>& btris)
{
  std::vector<Triangle> tris;

  for (auto btri : btris)
    tris.push_back(btri.triangle);

  return tris;
}

Scene::Scene(const std::string& path)
{
  std::ifstream file(path);
  if (!file) 
  {
    std::cerr << "Error while opening file: " << path << std::endl;
    exit(1);
  }

  json j = json::parse(file);
  auto btris_bounds = parse_scene(j);
  auto btris = btris_bounds.first;
  auto bounds = btris_bounds.second;
  std::vector<BVHNode> tree;

  build_bvh(btris, bounds, 0, btris.size(), tree, 0);
  bvh_ = tree;
  triangles_ = BtriToTri(btris);
}

Vertex3 minimize(const Vertex3& a, const Vertex3& b)
{
  return Vertex3(std::min(a[0], b[0]),
                    std::min(a[1], b[1]),
                    std::min(a[2], b[2]));
}

Vertex3 maximize(const Vertex3& a, const Vertex3& b)
{
  return Vertex3(std::max(a[0], b[0]),
                    std::max(a[1], b[1]),
                    std::max(a[2], b[2]));
}

std::pair<std::vector<BVHTriangle>, Box> Scene::parse_scene(json j)
{
  auto c = j.at("camera").get<std::array<float, 3>>();
  cam_pos_ = Vertex3(c[0], c[1], c[2]);

  Vertex3 bmin = Vertex3(std::numeric_limits<float>::max());
  Vertex3 bmax = Vertex3(std::numeric_limits<float>::min());

  std::vector<BVHTriangle> btris;

  for (auto const& j_mod : j.at("models"))
  {
    auto obj_file = j_mod.at("objFile").get<std::string>();
    auto mtl_basedir = j_mod.at("mtlBasedir").get<std::string>();
    auto scale = j_mod.at("scale").get<float>();
    auto t = j_mod.at("translation").get<std::array<float, 3>>();
    Vertex3 translation = Vertex3(t[0], t[1], t[2]);

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
    size_t mat_offset = materials_.size();
    for (const auto& mtl : mtls) {
      materials_.push_back({
        glm::vec4(mtl.diffuse[0], mtl.diffuse[1], mtl.diffuse[2], 0.0),                         // kd
        glm::vec4(mtl.emission[0], mtl.emission[1], mtl.emission[2], 0.0),                      // ke
        glm::vec4(mtl.specular[0], mtl.specular[1], mtl.specular[2], 0.0),                      // ks
        glm::vec4(mtl.transmittance[0], mtl.transmittance[1], mtl.transmittance[2], mtl.ior)    // kt + ior
      });
      //std::cout << mtl.specular[0] << " " <<  mtl.specular[1] << " " << mtl.specular[2] << std::endl;    // kt + ior
      //std::cout << mtl.transmittance[0] << " " <<  mtl.transmittance[1] << " " << mtl.transmittance[2] << std::endl;    // kt + ior
    }

    std::cout << "Material parsed" << std::endl;

    // Vertices
    size_t vertex_offset = vertices_.size();
    for (size_t i = 0; i < attrib.vertices.size(); i+=3)
    {
      Vertex vertex(
        attrib.vertices[i] * scale + translation.x,
        attrib.vertices[i + 1] * scale + translation.y,
        attrib.vertices[i + 2] * scale + translation.z,
        0.f
      );
      vertices_.push_back(vertex);

      Vertex3 v3 = vertex.xyz();
      bmin = minimize(v3, bmin);
      bmax = maximize(v3, bmax);
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
        triangle.mat_id = shape.mesh.material_ids[i_face] + mat_offset;

        //Actual vertices for BVH triangle
        Vertex3 verts[3];
        // Iterate inside of triangles
        for (size_t v = 0; v < face_v; ++v)
        {
          size_t idx = shape.mesh.indices[idx_offset + v].vertex_index + vertex_offset;
          triangle.vertices_index[v] = idx;
          verts[v] = vertices_[idx].xyz();
        }

        Vertex3 bmin_tri = minimize(minimize(verts[0], verts[1]), verts[2]);
        Vertex3 bmax_tri = maximize(maximize(verts[0], verts[1]), verts[2]);
        Box bounds = { bmin_tri, bmax_tri };
        BVHTriangle btri(verts, bounds, triangle);
        btris.push_back(btri);

        if (triangle.mat_id != -1 && materials_[triangle.mat_id].ke != glm::vec4(0.f))
        {
          lights_.push_back(triangle);
        }


        //triangles_.push_back(triangle);
        idx_offset += face_v;
        i_face++;
      }
    }

    std::cout << "Lights " << lights_.size() << std::endl;

    std::cout << "Triangles parsed" << std::endl;

    std::cout << "Loaded " << vertices_.size()<< " vertices" << std::endl;

    std::cout << "Shape size: " << shapes.size() << std::endl;
    std::cout << "Face size: " << shapes[0].mesh.num_face_vertices.size() << std::endl;
    std::cout << "Indices size: " << shapes[0].mesh.indices.size() << std::endl;
    std::cout << "Pos size: " << attrib.vertices.size() << std::endl;
    std::cout << "Norm size: " << attrib.normals.size() << std::endl;
    std::cout << "Tex size: " << attrib.texcoords.size() << std::endl;

    //std::cout << "Bounds are bmin: (" << bmin[0] << ", " << bmin[1] << ", " << bmin[2] << ")\n";
    //std::cout << "Bounds are bmax: (" << bmax[0] << ", " << bmax[1] << ", " << bmax[2] << ")\n";
  }
  
  return { btris , { bmin, bmax }};
}

Box group(const Box& b1, const Box& b2)
{
  return { minimize(b1.bmin, b2.bmin), maximize(b1.bmax, b2.bmax) };
}

//We consider the axis with largest gap the "best axis"
//sensible to edge cases
int get_best_axis(Box& box)
{
  Vertex3 diff = box.bmax - box.bmin;
  int max = std::max(std::max(diff[0], diff[1]), diff[2]);

  if (max == diff[0])
    return 0;
  if (max == diff[1])
    return 1;

  return 2;
}

float area(Box& box)
{
  Vertex3 diff = box.bmax - box.bmin;
  return diff[0] * diff[1] * diff[2];
}

//taabb time for intersection between ray and box
//ttri time for intersection between ray and tri
float SAH(float taabb, float ttri, Box& bounds, Box& boundsL, int ntrisL,
          Box& boundsR, int ntrisR)
{
  float aBounds = area(bounds);
  return 2 * taabb + area(boundsL) / aBounds * ntrisL * ttri
                  + area(boundsR) /  aBounds * ntrisR * ttri;
}

//Box BtriToBox(BVHTriangle& btri)
//{
//  return {
//    .bmin=(minimize(minimize(btri.verts[0], btri.verts[1]), btri.verts[2])),
//    .bmax=(maximize(maximize(btri.verts[0], btri.verts[1]), btri.verts[2]))
//  };
//}

Vertex3 compute_bmaxL(Box& bounds, Box& boundsR, int axis)
{
  if (axis == 0)
    return Vertex3(boundsR.bmin[0], bounds.bmax[1], bounds.bmax[2]);
  if (axis == 1)
    return Vertex3(bounds.bmax[0], boundsR.bmin[1], bounds.bmax[2]);

  return Vertex3(bounds.bmax[0], bounds.bmax[1], boundsR.bmin[2]);
}

Vertex3 compute_vmidL(Box& bounds, int axis)
{
  if (axis == 0)
    return Vertex3(bounds.bmin[0] + (bounds.bmax[0] - bounds.bmin[0]) / 2 , bounds.bmax[1], bounds.bmax[2]);
  if (axis == 1)
    return Vertex3(bounds.bmax[0], bounds.bmin[1] + (bounds.bmax[1] - bounds.bmin[1]) / 2 ,  bounds.bmax[2]);

  return Vertex3(bounds.bmax[0], bounds.bmax[1], bounds.bmin[2] + (bounds.bmax[2] - bounds.bmin[2]) / 2);
}

Vertex3 compute_vmidR(Box& bounds, int axis)
{
  if (axis == 0)
    return Vertex3(bounds.bmin[0] + (bounds.bmax[0] - bounds.bmin[0]) / 2 , bounds.bmin[1], bounds.bmin[2]);
  if (axis == 1)
    return Vertex3(bounds.bmin[0], bounds.bmin[1] + (bounds.bmax[1] - bounds.bmin[1]) / 2 ,  bounds.bmin[2]);

  return Vertex3(bounds.bmin[0], bounds.bmin[1], bounds.bmin[2] + (bounds.bmax[2] - bounds.bmin[2]) / 2);
}

//We consider time to build a triangle is 1
//We will make further experiences later
float taabb = 1;
float ttri = 1.2;
void Scene::build_bvh(std::vector<BVHTriangle>& tris, Box& bounds, int begin, int end,
                      std::vector<BVHNode>& tree, int idx)
{
  //Purpose of this algorithm is build an optimized BVH
  //based on time to calcutate intersection with triangle and plane
  float best_cost = ttri * tris.size();
  int best_axis = -1;
  int best_event = -1;

  //Find best axis
  int axis = get_best_axis(bounds);
  axis = std::rand() % 3;
  std::cout << "Sorted by axis: " << axis << std::endl;
  
  //Sort using best axis
  if (axis == 0)// X 
    //Could be replaced by qsort
    std::sort(tris.begin() + begin, tris.begin() + end, vertexX);
  else if (axis == 1)// Y 
    //Could be replaced by qsort
    std::sort(tris.begin() + begin, tris.begin() + end, vertexY);
  else// Z 
    //Could be replaced by qsort
    std::sort(tris.begin() + begin, tris.begin() + end, vertexZ);

  std::vector<float> left_area;
  std::vector<Box> left_boxes;
  std::vector<float> right_area;

  //left_area.push_back(std::numeric_limits<float>::max())
  /*Box left_box = tris[begin].bounds;
  for (int i = begin; i < end; ++i)
  {
    left_box = group(left_box, tris[i].bounds);
    left_area.push_back(area(left_box));
    left_boxes.push_back(left_box);
  }

  Box right_box = tris[end - 1].bounds;
  float cost = best_cost;

  Box best_boxr = right_box;
  Box best_boxl = left_box;//left_area[end - 2];
  for (int i = end - 1; i > begin; --i)
  {
    right_box = group(right_box, tris[i].bounds);
    //Vertex3 bmaxL = compute_bmaxL(bounds, right_box, axis);
    left_box = left_boxes[i - 1];//{ bounds.bmin, bmaxL };
    cost = SAH(taabb, ttri, bounds, left_box, i, right_box, end - i);

    if (cost < best_cost)
    {
      best_cost = cost;
      best_event = i;
      best_boxr = right_box;
      best_boxl = left_box;
    }
  }*/

  //std::cout << "Best event is: " << best_event << std::endl;

  //std::cout << "Area of computed left: " << area(best_boxl) <<std::endl;
  //std::cout << "Actual area of left: " << left_area[best_event - 1] << std::endl;

  if ((int)tree.size() <= idx)
    tree.resize(idx + 1);
  if (end - begin < 30)
  //if (best_event == -1) //then {found no partition better than leaf}
  {
    BVHNode node = { bounds.bmin, begin, bounds.bmax, end - begin };
    tree[idx] = node;
    //tree.push_back(node);
  }
  else
  {
    BVHNode node = { bounds.bmin, 2 * idx + 1, bounds.bmax, 0 };
    tree[idx] = node;
    //tree.push_back(node);

    int mid = begin + (end - begin) / 2;
    //Vertex3 vmidL = compute_vmidL(bounds, axis);//bounds.bmin + (bounds.bmax - bounds.bmin) / 2;
    //Vertex3 vmidR = compute_vmidR(bounds, axis);

    

    Box best_boxl = { bounds.bmin, bounds.bmin };//vmidL };
    for (int i = begin; i < mid; ++i)
      best_boxl = group(best_boxl, tris[i].bounds);
    Box best_boxr = { /*vmidR*/bounds.bmax, bounds.bmax };
    for (int i = end - 1; i >= mid; --i)
      best_boxr = group(best_boxr, tris[i].bounds);
    //int mid = best_event;
    build_bvh(tris, best_boxl, begin, mid, tree, 2 * idx + 1);
    build_bvh(tris, best_boxr, mid, end, tree, 2 * idx + 2);
  }
}
