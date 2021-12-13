#include <glm/ext/matrix_clip_space.hpp>
#include <glm/ext/matrix_float4x4.hpp>
#include <glm/ext/matrix_transform.hpp>
#include <glm/ext/quaternion_common.hpp>
#include <glm/ext/vector_float3.hpp>
#include <glm/fwd.hpp>
#include <glm/matrix.hpp>
#include <glm/trigonometric.hpp>
#include <iostream>
#include <GL/glew.h>
#include <GL/freeglut.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/vec3.hpp>
#include <vector>
#include <stdio.h>

#include <imgui.h>
#include "../bindings/imgui_impl_opengl3.h"
#include "../bindings/imgui_impl_glut.h"

#include "program.hh"
#include "sphere.hh"

#define TEST_OPENGL_ERROR()                                                             \
  do {                                                                                  \
    GLenum err = glGetError();                                                          \
    if (err != GL_NO_ERROR) std::cerr << "OpenGL ERROR: "                               \
                                      << gluErrorString(err)                            \
                                      << " file " << __FILE__                           \
                                      << " line " << __LINE__ << std::endl;             \
  } while(0)

#define numVAOs 1
#define numVBOs 2

GLuint vao[numVAOs];
GLuint vbo[numVBOs];

int width = 64;
int height = 64;

int workGroupsX = width;
int workGroupsY = height;
int workGroupsZ = 1;

GLuint screenTextureID; // The texture ID of the full screen texture
unsigned char *screenTexture; // The screen texture RGBA8888 color data

program* pathtrace_shader;
program* display_shader;


void display()
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);TEST_OPENGL_ERROR();

  // PHASE 1 : Path tracing compute shader
  pathtrace_shader->use();
  std::cout << "pathtracing" << std::endl;

  glBindImageTexture(0, screenTextureID, 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RGBA8);

  glDispatchCompute(workGroupsX, workGroupsY, workGroupsZ);
  glMemoryBarrier(GL_ALL_BARRIER_BITS);

  std::cout << +screenTexture[0] << std::endl;

  // PHASE 2 : Display the resulting texture
  display_shader->use();

  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, screenTextureID);
  glBindBuffer(GL_ARRAY_BUFFER, vbo[0]);
  glVertexAttribPointer(0, 3, GL_FLOAT, false, 0, 0);
  glEnableVertexAttribArray(0);
  glBindBuffer(GL_ARRAY_BUFFER, vbo[1]);
  glVertexAttribPointer(1, 2, GL_FLOAT, false, 0, 0);
  glEnableVertexAttribArray(1);
  glDrawArrays(GL_TRIANGLES, 0, 6);

  glutSwapBuffers(); TEST_OPENGL_ERROR();
}

void resize(int width, int height)
{
  glViewport(0, 0, width, height);TEST_OPENGL_ERROR();
}

bool initGlut(int &argc, char* argv[])
{
  glutInit(&argc, argv);TEST_OPENGL_ERROR();
  glutInitContextVersion(4, 5);TEST_OPENGL_ERROR();
  glutInitContextProfile(GLUT_CORE_PROFILE|GLUT_DEBUG);TEST_OPENGL_ERROR();
  glutInitDisplayMode(GLUT_RGBA|GLUT_DOUBLE|GLUT_DEPTH);TEST_OPENGL_ERROR();
  glutInitWindowSize(1024, 1024);TEST_OPENGL_ERROR();
  glutInitWindowPosition(10, 10);TEST_OPENGL_ERROR();
  glutCreateWindow ("Test OpenGL-POGL");TEST_OPENGL_ERROR();
  glutDisplayFunc(display);TEST_OPENGL_ERROR();
  glutReshapeFunc(resize);TEST_OPENGL_ERROR();
  return true;
}

void init()
{
  glEnable(GL_DEPTH_TEST);TEST_OPENGL_ERROR();
  glDepthFunc(GL_LESS);TEST_OPENGL_ERROR();
}

void init_uniform(program* instance)
{
  glm::vec3 camPos = glm::vec3(0.0f, 0.0f, 4.0f);
  glm::vec3 camUp = glm::vec3(0.0f, 1.0f, 0.0f);
  glm::vec3 camFront = glm::vec3(0.0f, 0.0f, -1.0f);

  glm::mat4 proj = glm::perspective(glm::radians(45.0f), 1.0f, 0.1f, 100.0f); 
  glm::mat4 model = glm::mat4(1.0f);
  glm::mat4 view = glm::lookAt(camPos, camPos + camFront, camUp);

  glm::mat4 locToProj = proj * view * model;

  GLint projection_location =
    glGetUniformLocation(instance->get_id(), "localToProjection");TEST_OPENGL_ERROR();
  glUniformMatrix4fv(projection_location, 1, GL_FALSE, &locToProj[0][0]);
  TEST_OPENGL_ERROR();
}



int main(int argc, char** argv)
{

  initGlut(argc, argv);

  GLenum err = glewInit();
  if (err != GLEW_OK)
    std::cout << glewGetErrorString(err) << std::endl;

  init();

  screenTexture = (unsigned char*)malloc(sizeof(unsigned char) * 4 * width * height);
  memset(screenTexture, 0, sizeof(char) * 4 * width * height);
  for (int i = 0; i < height; i++) {
    for (int j = 0; j < width; j++) {
      screenTexture[i * width * 4 + j * 4 + 0] = 42;
      screenTexture[i * width * 4 + j * 4 + 1] = 128;
      screenTexture[i * width * 4 + j * 4 + 2] = 255;
      screenTexture[i * width * 4 + j * 4 + 3] = 255;
    } 
  }

  // Create the OpenGL Texture on which to ray cast the scene
  glGenTextures(1, &screenTextureID);
  glBindTexture(GL_TEXTURE_2D, screenTextureID);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, (const void *)screenTexture);

  glGenVertexArrays(1, vao); TEST_OPENGL_ERROR();
  glGenBuffers(numVBOs, vbo); TEST_OPENGL_ERROR();

  // Create quad vertices and texture coordinates for rendering the finished texture to the window
  const float windowQuadVerts[ ] = {
    -1.0f, 1.0f, 0.3f, -1.0f,-1.0f, 0.3f, 1.0f, -1.0f, 0.3f,
    1.0f, -1.0f, 0.3f, 1.0f, 1.0f, 0.3f, -1.0f, 1.0f, 0.3f
  };
  const float windowQuadUVs[ ] = {
    0.0f, 1.0f, 0.0f, 0.0f, 1.0f, 0.0f,
    1.0f, 0.0f, 1.0f, 1.0f, 0.0f, 1.0f
  };
  
  glBindVertexArray(vao[0]);
  glBindBuffer(GL_ARRAY_BUFFER, vbo[0]); // vertex positions
  glBufferData(GL_ARRAY_BUFFER, sizeof(windowQuadVerts), windowQuadVerts, GL_STATIC_DRAW);
  glBindBuffer(GL_ARRAY_BUFFER, vbo[1]); // texture coordinates
  glBufferData(GL_ARRAY_BUFFER, sizeof(windowQuadUVs), windowQuadUVs, GL_STATIC_DRAW);

  std::string pathtrace_c("../src/shaders/pathtrace.comp.glsl");
  auto pathtrace_files = std::vector<std::pair<GLenum, std::string>>({
    { GL_COMPUTE_SHADER, pathtrace_c},
  });

  std::cout << "Compiling pathtrace shader" << std::endl;
  pathtrace_shader = program::make_program(pathtrace_files);
  std::cout << "Pathtrace shader compiled" << std::endl;

  std::string display_v("../src/shaders/display.vert.glsl");
  std::string display_f("../src/shaders/display.frag.glsl");
  
  auto display_files = std::vector<std::pair<GLenum, std::string>>({
    { GL_VERTEX_SHADER, display_v},
    { GL_FRAGMENT_SHADER, display_f}
  });
  display_shader = program::make_program(display_files);
  std::cout << "Display shader compiled" << std::endl;

  glutMainLoop();

  delete pathtrace_shader;
  delete display_shader;

  return 0;
}
