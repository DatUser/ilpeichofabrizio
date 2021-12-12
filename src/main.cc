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
#include <vector>

#include <imgui.h>
#include "../bindings/imgui_impl_opengl3.h"
#include "../bindings/imgui_impl_glut.h"

#include "program.hh"
#include "vector3.hh"
#include "sphere.hh"

#define TEST_OPENGL_ERROR()                                                             \
  do {                                                                                  \
    GLenum err = glGetError();                                                          \
    if (err != GL_NO_ERROR) std::cerr << "OpenGL ERROR: "                               \
                                      << gluErrorString(err)                            \
                                      << " file " << __FILE__                           \
                                      << " line " << __LINE__ << std::endl;             \
  } while(0)


std::vector<float> vertices;
std::vector<unsigned int> indices;

GLuint VBO;
GLuint VAO;
GLuint EBO;

void display()
{
  //glClearColor(1.0, 0.0, 0.0, 1.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);TEST_OPENGL_ERROR();
  //glMatrixMode(GL_MODELVIEW);
  //glLoadIdentity();
  //
  ImGui_ImplOpenGL3_NewFrame();
  ImGui_ImplGLUT_NewFrame();
  ImGui::NewFrame(); 

  static glm::vec4 color(1.0, 1.0, 1.0, 1.0);
  ImGui::ColorEdit3("color", &color[0]);


  glBindVertexArray(VAO);TEST_OPENGL_ERROR();
  glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0);TEST_OPENGL_ERROR();
  glBindVertexArray(0);TEST_OPENGL_ERROR();


  ImGui::Render();
  ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

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
  glutCreateWindow ("Test OpenGL−POGL");TEST_OPENGL_ERROR();
  glutDisplayFunc(display);TEST_OPENGL_ERROR();
  glutReshapeFunc(resize);TEST_OPENGL_ERROR();
  return true;
}

/*bool initGlew()
{
  return (glewInit() ==GLEW_OK);
}*/

void init()
{
  glEnable(GL_DEPTH_TEST);TEST_OPENGL_ERROR();
  glDepthFunc(GL_LESS);TEST_OPENGL_ERROR();
}

void init_vbo(program* instance)
{if (!instance) return;
  Point3 center(0, 0, 0);
  Sphere s(center, 1);

  auto data = s.generate_vertices(500, 500);
  vertices = data.first;
  indices = data.second;

  glGenVertexArrays(1, &VAO); TEST_OPENGL_ERROR();
  glGenBuffers(1, &VBO); TEST_OPENGL_ERROR();
  glGenBuffers(1, &EBO); TEST_OPENGL_ERROR();

  glBindVertexArray(VAO); TEST_OPENGL_ERROR();
  glBindBuffer(GL_ARRAY_BUFFER, VBO); TEST_OPENGL_ERROR();
  glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), &vertices[0], GL_STATIC_DRAW);   TEST_OPENGL_ERROR();

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO); TEST_OPENGL_ERROR();
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), 
               &indices[0], GL_STATIC_DRAW); TEST_OPENGL_ERROR();

  // Vertex pos
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE,
      6 * sizeof(float), (void*)0); TEST_OPENGL_ERROR();
  glEnableVertexAttribArray(0); TEST_OPENGL_ERROR();

  // Vertex normal
  glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), 
      (void*)(3 * sizeof(float)));TEST_OPENGL_ERROR();
  glEnableVertexAttribArray(1); TEST_OPENGL_ERROR();

  glBindVertexArray(0);
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

  std::string file_v("../src/shaders/pbr.vert.glsl");
  std::string file_f("../src/shaders/pbr.frag.glsl");
  
  auto shaders_src = std::vector<std::pair<GLenum, std::string>>({
                            { GL_VERTEX_SHADER, file_v},
                            { GL_FRAGMENT_SHADER, file_f}
                          });
  program* instance = program::make_program(shaders_src);

  instance->use();

  init_vbo(instance);
  init_uniform(instance);

  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  ImGuiIO& io = ImGui::GetIO();
  (void)io;

	ImGui::StyleColorsDark();

	ImGui_ImplGLUT_Init();
	ImGui_ImplGLUT_InstallFuncs();
	ImGui_ImplOpenGL3_Init();

  glutMainLoop();
  
  ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGLUT_Shutdown();
	ImGui::DestroyContext();

  delete instance;

  return 0;
}
