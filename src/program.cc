#include "program.hh"

// Loads shader from file
std::string load_from_file(std::string &file)
{
  std::ifstream in(file);

  if (in.is_open())
  {
    std::string data{};
    std::string line{};

    while (getline(in, line))
      data += line + '\n';

    in.close();

    return data;
  }

  std::cout << "Unable to open file";
  return std::string{};
}

program *program::make_program(std::vector<std::pair<GLenum, std::string>> &shaders_src)
{
  // Create new program
  program *instance = new program();

  // Load and compile shaders
  for (auto shader : shaders_src)
    instance->shaders.push_back(load_and_compile(shader.first, shader.second));

  // Create program and attach shaders
  instance->create_program();

  return instance;
}

GLuint load_and_compile(GLenum shader_type, std::string &file)
{
  // Create shader
  GLuint shader = glCreateShader(shader_type);
  TEST_OPENGL_ERROR();

  // Load shader from file
  std::string data = load_from_file(file);
  const GLchar *source = data.c_str();

  glShaderSource(shader, 1, &source, 0);
  TEST_OPENGL_ERROR();

  // Compiling shader
  glCompileShader(shader);
  TEST_OPENGL_ERROR();

  GLint isCompiled = 0;
  glGetShaderiv(shader, GL_COMPILE_STATUS, &isCompiled);
  TEST_OPENGL_ERROR();

  if (isCompiled == GL_FALSE)
  {
    GLint maxLength = 0;
    glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &maxLength);

    // MaxLength includes the NULL character
    std::vector<GLchar> errorLog(maxLength);
    glGetShaderInfoLog(shader, maxLength, &maxLength, &errorLog[0]);

    std::cout << "Out: " << std::string(errorLog.begin(), errorLog.end()) << '\n';
    // Exit with error
    glDeleteShader(shader);

    errx(1, "Could not compile shader: %s", std::string(errorLog.begin(), errorLog.end()).c_str());
  }

  // Return compiled shader
  return shader;
}

void program::create_program()
{
  // Create program
  id = glCreateProgram();
  TEST_OPENGL_ERROR();

  // Attach shaders to program
  for (GLuint shader : shaders)
    glAttachShader(id, shader);

  // Link program
  glLinkProgram(id);
  TEST_OPENGL_ERROR();

  GLint isLinked = 0;
  glGetProgramiv(id, GL_LINK_STATUS, (int *)&isLinked);
  TEST_OPENGL_ERROR();

  // Check error
  if (isLinked != GL_TRUE)
  {
    GLint maxLength = 0;
    glGetProgramiv(id, GL_INFO_LOG_LENGTH, &maxLength);

    // The maxLength includes the NULL character
    std::vector<GLchar> infoLog(maxLength);
    glGetProgramInfoLog(id, maxLength, &maxLength, &infoLog[0]);

    // Delete program
    glDeleteProgram(id);

    // Delete shaders
    for (GLuint shader : shaders)
      glDeleteShader(shader);

    errx(1, "Error creating program: %s", std::string(infoLog.begin(), infoLog.end()).c_str());
  }
}

void program::use()
{
  glUseProgram(id);
  TEST_OPENGL_ERROR();
}

void program::set_int(const std::string &name, int value)
{
  GLint64 uniformId = glGetUniformLocation(id, name.c_str());
  glUniform1i(uniformId, value);
}

void program::set_float(const std::string &name, float value)
{
  GLint64 uniformId = glGetUniformLocation(id, name.c_str());
  glUniform1f(uniformId, value);
}

void program::set_vec3(const std::string &name, const glm::vec3 &value)
{
  GLint64 uniformId = glGetUniformLocation(id, name.c_str());
  glUniform3fv(uniformId, 1, &value[0]);
}

void program::set_vec4(const std::string &name, const glm::vec4 &value)
{
  GLint64 uniformId = glGetUniformLocation(id, name.c_str());
  glUniform4fv(uniformId, 1, &value[0]);
}

void program::set_matrix4(const std::string &name, const glm::mat4 &value)
{
  GLint64 uniformId = glGetUniformLocation(id, name.c_str());
  glUniformMatrix4fv(uniformId, 1, GL_FALSE, &value[0][0]);
}