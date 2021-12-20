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

#include "../SOIL2/src/SOIL2/SOIL2.h"

#include "program.hh"
#include "sphere.hh"

//std::vector<float> vertices;
//std::vector<unsigned int> indices;
//std::vector<float> tex;
Mesh* mesh;

GLuint VBOs[3];
GLuint VBO;
GLuint VAO;
GLuint EBO;
GLuint TBO;
GLuint aoTex;
GLuint colorTex;
GLuint metalTex;
GLuint roughTex;

glm::vec2 start_pos;
glm::vec2 end_pos;

program* instance;
float anim_time = 0.0;

void set_ui()
{
  static glm::vec3 color(1.0, 1.0, 1.0);
  ImGui::ColorEdit3("color", &color[0]);
  instance->set_vec3("albedo", color);TEST_OPENGL_ERROR();

  static glm::vec3 f0(0.04);//default water value
  ImGui::ColorEdit3("f0", &f0[0]);
  instance->set_vec3("f0", f0);TEST_OPENGL_ERROR();

  static float roughness = 0.0;
  ImGui::SliderFloat("roughness", &roughness, 0.00, 1.0);
  instance->set_float("roughness", roughness);

  static float metalness = 0.0;
  ImGui::SliderFloat("metalness", &metalness, 0.00, 1.0);
  instance->set_float("metalness", metalness);

  static float ao = 1.0;
  ImGui::SliderFloat("Ambient Occlusion", &ao, 0.00, 1.0);
  instance->set_float("a_occlusion", ao);
}

void display()
{
  ImGui_ImplOpenGL3_NewFrame();
  ImGui_ImplGLUT_NewFrame();

  glClearColor(0.45, 0.6, 0.99, 1.0);TEST_OPENGL_ERROR();
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);TEST_OPENGL_ERROR();
  //glMatrixMode(GL_MODELVIEW);TEST_OPENGL_ERROR();
  //glLoadIdentity();TEST_OPENGL_ERROR();
  //
  //ImGui::NewFrame();
  /*if (anim_time > 2.0)
  {
    glm::vec3 camPos = glm::vec3(4.0f, 0.0f, 0.0f);
    glm::vec3 camUp = glm::vec3(0.0f, 1.0f, 0.0f);
    glm::vec3 camFront = glm::vec3(-1.0f, 0.0f, 0.0f);
    instance->set_vec3("cameraPos", camPos);
    glm::mat4 view = glm::lookAt(camPos, camPos + camFront, camUp);
    instance->set_matrix4("view", view);
  }
  if (anim_time > 4.0)
  {
    glm::vec3 camPos = glm::vec3(0.0f, 0.0f, -4.0f);
    glm::vec3 camUp = glm::vec3(0.0f, 1.0f, 0.0f);
    glm::vec3 camFront = glm::vec3(0.0f, 0.0f, 1.0f);
    instance->set_vec3("cameraPos", camPos);
    glm::mat4 view = glm::lookAt(camPos, camPos + camFront, camUp);
    instance->set_matrix4("view", view);
  }
  if (anim_time > 6.0)
  {
    glm::vec3 camPos = glm::vec3(-4.0f, 0.0f, 0.0f);
    glm::vec3 camUp = glm::vec3(0.0f, 1.0f, 0.0f);
    glm::vec3 camFront = glm::vec3(1.0f, 0.0f, 0.0f);
    instance->set_vec3("cameraPos", camPos);
    glm::mat4 view = glm::lookAt(camPos, camPos + camFront, camUp);
    instance->set_matrix4("view", view);
  }*/

  ImGui::Begin("Parameters");

  set_ui();

  //glActiveTexture(GL_TEXTURE0);
  //glBindTexture(GL_TEXTURE_2D, colorTex);

  mesh->render();
  /*glBindVertexArray(VAO);TEST_OPENGL_ERROR();
  glBindBuffer(GL_ARRAY_BUFFER, VBOs[1]);TEST_OPENGL_ERROR();
  glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0);TEST_OPENGL_ERROR();
  glBindVertexArray(0);TEST_OPENGL_ERROR();*/

  ImGui::End();

  ImGui::Render();
  ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

  //glLoadIdentity();

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
{
  if (!instance) return;

  glm::vec3 center(0, 0, 0);
  Sphere s(center, 1);

  auto data = s.generate_vertices(500, 500);
  auto vertices = data.first;
  auto indices = data.second;
  //tex = s.get_tex();

  mesh = new Mesh(vertices, indices);
  mesh->attach_shader(instance);
  mesh->init_mesh();
  /*
  glGenVertexArrays(1, &VAO); TEST_OPENGL_ERROR();
  glGenBuffers(3, VBOs); TEST_OPENGL_ERROR();
  //glGenBuffers(1, &EBO); TEST_OPENGL_ERROR();
  //glGenBuffers(1, &TBO); TEST_OPENGL_ERROR();

  glBindVertexArray(VAO); TEST_OPENGL_ERROR();

  // Vertex pos
  glBindBuffer(GL_ARRAY_BUFFER, VBOs[0]); TEST_OPENGL_ERROR();
  glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), &vertices[0], GL_STATIC_DRAW);   TEST_OPENGL_ERROR();
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0); TEST_OPENGL_ERROR();
  glEnableVertexAttribArray(0); TEST_OPENGL_ERROR();

  // Vertex normal
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, VBOs[1]); TEST_OPENGL_ERROR();
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), 
               &indices[0], GL_STATIC_DRAW); TEST_OPENGL_ERROR();
  glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, 0);TEST_OPENGL_ERROR();
  glEnableVertexAttribArray(1); TEST_OPENGL_ERROR();

  // Vertex tex coords
  glBindBuffer(GL_ARRAY_BUFFER, VBOs[1]);TEST_OPENGL_ERROR();
  glBufferData(GL_ARRAY_BUFFER, tex.size() * sizeof(float), &tex[0], GL_STATIC_DRAW);
  TEST_OPENGL_ERROR();
  glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 0, 0); TEST_OPENGL_ERROR();
  glEnableVertexAttribArray(2);

  glBindVertexArray(0);*/
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
      break;
    case GLUT_UP:
      io.MouseDown[button] = false;
      end_pos.x = (!button) ? x : end_pos.x;
      end_pos.x = (!button) ? y : end_pos.y;
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
  instance->set_vec3("cameraPos", camPos);


  glm::mat4 proj = glm::perspective(glm::radians(45.0f), 1.0f, 0.1f, 100.0f); 
  glm::mat4 model = glm::mat4(1.0f);
  glm::mat4 view = glm::lookAt(camPos, camPos + camFront, camUp);

  //glm::mat4 modelView = model * view;
  //glm::mat4 projection = proj;

  instance->set_matrix4("model", model);
  instance->set_matrix4("view", view);
  instance->set_matrix4("projection", proj);
}

void load_texture(GLuint* id, const char* filename)
{
  *id = SOIL_load_OGL_texture(filename, SOIL_LOAD_AUTO, SOIL_CREATE_NEW_ID,
                                SOIL_FLAG_MIPMAPS);
  
  if (!*id)
    std::cout << "Could not load " << filename << std::endl;
  
  glBindTexture(GL_TEXTURE_2D, *id);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  /*colorTex = SOIL_load_OGL_texture("../resources/metal/color.jpg", SOIL_LOAD_RGB,
                                SOIL_CREATE_NEW_ID, SOIL_FLAG_MIPMAPS);
  metalTex = SOIL_load_OGL_texture("../resources/metal/metal.jpg", SOIL_LOAD_RGB,
                                SOIL_CREATE_NEW_ID, SOIL_FLAG_MIPMAPS);
  roughTex = SOIL_load_OGL_texture("../resources/metal/rough.jpg", SOIL_LOAD_RGB,
                                SOIL_CREATE_NEW_ID, SOIL_FLAG_MIPMAPS);*/
}

void anim() {
  GLint anim_time_location;
  anim_time_location = glGetUniformLocation(instance->get_id(), "time");TEST_OPENGL_ERROR();
  glUniform1f(anim_time_location, anim_time);TEST_OPENGL_ERROR();
  anim_time += 0.1;
  glutPostRedisplay();
}

void timer(int value) {
  (void) value;
  anim();
  glutTimerFunc(33, timer, 0);
}

void init_anim() {
  glutTimerFunc(33, timer, 0);
}

void display_compute()
{
  int wk_grp_cnt[3];
  int wk_grp_siz[3];
  int wk_grp_inv;

  glGetIntegeri_v(GL_MAX_COMPUTE_WORK_GROUP_COUNT, 0, &wk_grp_cnt[0]);
  glGetIntegeri_v(GL_MAX_COMPUTE_WORK_GROUP_COUNT, 1, &wk_grp_cnt[1]);
  glGetIntegeri_v(GL_MAX_COMPUTE_WORK_GROUP_COUNT, 2, &wk_grp_cnt[2]);

  glGetIntegeri_v(GL_MAX_COMPUTE_WORK_GROUP_SIZE, 0, &wk_grp_siz[0]);
  glGetIntegeri_v(GL_MAX_COMPUTE_WORK_GROUP_SIZE, 1, &wk_grp_siz[1]);
  glGetIntegeri_v(GL_MAX_COMPUTE_WORK_GROUP_SIZE, 2, &wk_grp_siz[2]);

  glGetIntegerv(GL_MAX_COMPUTE_WORK_GROUP_INVOCATIONS, &wk_grp_inv);

  printf("max num of wk gps: %i %i %i\n", wk_grp_cnt[0], wk_grp_cnt[1], wk_grp_cnt[2]);
  printf("max siz of wk gps: %i %i %i\n", wk_grp_siz[0], wk_grp_siz[1], wk_grp_siz[2]);
  printf("max local wk gps inv: %i\n", wk_grp_inv);
}

int main(int argc, char** argv)
{

  initGlut(argc, argv);

  GLenum err = glewInit();
  if (err != GLEW_OK)
    std::cout << glewGetErrorString(err) << std::endl;

  init();

  //std::string file_f("../src/shaders/test.frag.glsl");
  std::string file_v("../src/shaders/pbr.vert.glsl");
  std::string file_f("../src/shaders/pbr.frag.glsl");
  //std::string file_v("../src/shaders/glint.vert.glsl");
  //std::string file_f("../src/shaders/glint.frag.glsl");
  
  auto shaders_src = std::vector<std::pair<GLenum, std::string>>({
                            { GL_VERTEX_SHADER, file_v},
                            { GL_FRAGMENT_SHADER, file_f}
                          });
  instance = program::make_program(shaders_src);

  instance->use();

  init_vbo(instance);
  init_uniform(instance);
  load_texture(&colorTex, "../resources/metal/color.jpg");
  init_anim();

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