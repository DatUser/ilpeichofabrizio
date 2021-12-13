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


std::vector<float> vertices;
std::vector<unsigned int> indices;

GLuint VBO;
GLuint VAO;
GLuint EBO;

glm::vec2 start_pos;
glm::vec2 end_pos;

program* instance; 

void display()
{
  ImGui_ImplOpenGL3_NewFrame();
  ImGui_ImplGLUT_NewFrame();

  glClearColor(0.45, 0.6, 0.99, 1.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);TEST_OPENGL_ERROR();
  //glMatrixMode(GL_MODELVIEW);TEST_OPENGL_ERROR();
  //glLoadIdentity();TEST_OPENGL_ERROR();
  //
  //ImGui::NewFrame();

  ImGui::Begin("Parameters");

  static glm::vec3 color(1.0, 1.0, 1.0);
  ImGui::ColorEdit3("color", &color[0]);
  instance->set_vec3("albedo", color);TEST_OPENGL_ERROR();


  glBindVertexArray(VAO);TEST_OPENGL_ERROR();
  glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0);TEST_OPENGL_ERROR();
  glBindVertexArray(0);TEST_OPENGL_ERROR();

  ImGui::End();

  ImGui::Render();
	//ImGuiIO& io = ImGui::GetIO();
	//glViewport(0, 0, (GLsizei)io.DisplaySize.x, (GLsizei)io.DisplaySize.y);
  ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

  glutSwapBuffers(); TEST_OPENGL_ERROR();
  glutPostRedisplay(); TEST_OPENGL_ERROR();
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
  glm::vec3 center(0, 0, 0);
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

void mouseFunc(int glut_button, int state, int x, int y)
{
  ImGuiIO& io = ImGui::GetIO();
  io.MousePos = ImVec2((float) x, (float) y);

  int button;
  switch (glut_button)
  {
  case GLUT_LEFT_BUTTON:
    button = 0;
    break;
  case GLUT_RIGHT_BUTTON:
    button = 1;
    break;
  case GLUT_MIDDLE_BUTTON:
    button = 2;
    break;
  default:
    button = -1;
    break;
  }

  if (button >= 0)
  {
    switch (state)
    {
    case GLUT_DOWN:
      io.MouseDown[button] = true;
      start_pos.x = (!button) ? x : start_pos.x;
      start_pos.x = (!button) ? y : start_pos.y;
      std::cout << "Click: " << x << " - " << y << std::endl;
      break;
    case GLUT_UP:
      io.MouseDown[button] = false;
      end_pos.x = (!button) ? x : end_pos.x;
      end_pos.x = (!button) ? y : end_pos.y;
      std::cout << "Drop: " << x << " - " << y << std::endl;
      break;
    
    default:
      break;
    }
  }
}

void motionFunc(int x, int y)
{
  end_pos.x = x;
  end_pos.y = y;
  ImGuiIO& io = ImGui::GetIO();
  io.MousePos = ImVec2((float) x, (float) y);
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
  instance = program::make_program(shaders_src);

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
  glutMotionFunc(motionFunc);
  glutMouseFunc(mouseFunc);
	ImGui_ImplOpenGL3_Init();

  glutMainLoop();
  
  ImGui_ImplOpenGL3_Shutdown();
	ImGui_ImplGLUT_Shutdown();
	ImGui::DestroyContext();

  delete instance;

  return 0;
}
