#pragma once

#include <vector>
#include <glm/vec3.hpp>
#include "mesh.hh"

class Sphere
{
public:
    Sphere(glm::vec3 center, float radius);
    
    //Returns vertices and indices
    std::pair<std::vector<Vertex>, std::vector<unsigned int>>
    generate_vertices(unsigned int stacks, unsigned int sectors);

    std::vector<float>& get_normals() { return normals; }
    std::vector<float>& get_tex() { return tex; }

    ~Sphere() = default;
private:
    glm::vec3 center;
    float radius;

    std::vector<float> normals;
    std::vector<float> tex;
};