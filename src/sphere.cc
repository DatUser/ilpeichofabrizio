#include "sphere.hh"

#include <iostream>
#include <cmath>


Sphere::Sphere(glm::vec3 center, float radius)
: center(center),
  radius(radius),
  normals(std::vector<float>{}),
  tex(std::vector<float>{})
{    
}

/*
** x: left (->)
** y: up
** z: away from you
*/
std::pair<std::vector<Vertex>, std::vector<unsigned int>>
Sphere::generate_vertices(unsigned int lines, unsigned int cols)
{
    //std::vector<float> verts{};//(lines * cols * 3);
    std::vector<Vertex> verts{};
    std::vector<unsigned int> indices{};//((lines - 1) * (cols - 1) * 4);

    if (lines < 2 || cols < 3)
    {
        std::cerr << "Lines: " << lines << " | Cols: " << cols << std::endl
                    << "Lines must bigger than 1 and Cols bigger than 2"; 
        return {verts, indices};
    }

    verts.push_back({
        glm::vec3(0.0, radius, 0.0),
        glm::vec3(0.0, radius, 0.0),
        glm::vec2(0.0, 0.0)});

    //Top vertex
    //verts.push_back(0);
    //verts.push_back(radius);
    //verts.push_back(0);

    //Top normal
    //verts.push_back(0);
    //verts.push_back(1);
    //verts.push_back(0);

    for (unsigned int i = 0; i < lines - 1; ++i)
    {
        float line = i;
        double r = M_PI * (i + 1) / lines;
        for (unsigned int j = 0; j < cols; ++j)
        {
            double s = 2 *  M_PI * j / cols;

            double x = std::sin(r) * std::cos(s) * radius;
            double y = std::cos(r) * radius;
            double z = std::sin(r) * std::sin(s) * radius;

            float col = j;
            float u = (line + 1) / lines;
            float v = col / cols;

            verts.push_back({
                glm::vec3(x, y, z),
                glm::vec3(x, y, z),
                glm::vec2(u, v)});

            //verts.push_back(x * radius);
            //verts.push_back(y * radius);
            //verts.push_back(z * radius);

            //verts.push_back(x);
            //verts.push_back(y);
            //verts.push_back(z);
        }
        
    }
    verts.push_back({
        glm::vec3(0.0, -radius, 0.0),
        glm::vec3(0.0, -radius, 0.0),
        glm::vec2(1.0, 1.0)});

    //Bottom vertex
    //verts.push_back(0);
    //verts.push_back(-radius);
    //verts.push_back(0);
    
    //Bottom normal
    //verts.push_back(0);
    //verts.push_back(-1);
    //verts.push_back(0);

    //Top and Bot
    for (unsigned int i = 0; i < cols; ++i)
    {
        unsigned int i0 = i + 1;
        unsigned int i1 = (i + 1) % cols + 1;
        indices.push_back(0);
        indices.push_back(i1);
        indices.push_back(i0);

        i0 = i + cols * (lines - 2) + 1;
        i1 = (i + 1) % cols + cols * (lines - 2) + 1;
        indices.push_back(verts.size() / 6 - 1);
        indices.push_back(i0);
        indices.push_back(i1);
    }

    tex.push_back(0);
    tex.push_back(0);
    for (unsigned int j = 0; j < lines - 2; j++)
    {
        unsigned int j0 = j * cols + 1;
        unsigned int j1 = (j + 1) * cols + 1;
        for (unsigned int i = 0; i < cols; i++)
        {
            unsigned int i0 = j0 + i;
            unsigned int i1 = j0 + (i + 1) % cols;
            unsigned int i2 = j1 + (i + 1) % cols;
            unsigned int i3 = j1 + i;
            

            indices.push_back(i0);
            indices.push_back(i1);
            indices.push_back(i3);

            indices.push_back(i3);
            indices.push_back(i1);
            indices.push_back(i2);

            //We skip the first line
            //tex.push_back(u);            
            //tex.push_back(v);
        }
    }
    //tex.push_back(1);
    //tex.push_back(1);

    return {verts, indices};
}
