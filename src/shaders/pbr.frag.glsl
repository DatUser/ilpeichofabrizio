#version 450

precision highp float;

layout (binding = 0) uniform sampler2D tex;

out vec4 outFragColor;

//uniform vec3 albedo;
uniform vec3 f0; //Pre-computed (default is water value)
uniform float roughness;
uniform float metalness;// = 0.1;//1 if metallic 0 otherwise
uniform float a_occlusion;
float lightPower = 4.0;



in vec3 worldPos;
in vec3 vNormalWS;
in vec2 uv;

const float PI = 3.1415926535897932384626433832795;
const float E = 0.000001;

float isPos(float nb)
{
  return (nb > 0.0) ? 1.0 : 0.0;
}


// Use dotNH for microdetails
vec3 fresnelSchlick(float dotHV, vec3 albedo)
{
  //vec3 F0 = f0;
  vec3 F0 = mix(f0, albedo, metalness);
  vec3 f90 = vec3(1.0);//Pre-computed (here we use water value)
  return F0 + (f90 - F0) * pow(1.0 - dotHV, 5.0);
}

float distribGGX(float dotNH, float alpha2)
{
  float dotNH2 = pow(dotNH, 2.0);
  float bot = dotNH2 * (alpha2 - 1.0) + 1.0;
  return alpha2 / (PI * bot * bot + E);
}

float geometrySmith(float dotNV, float dotNL, float alpha2)
{
  float kdirect = pow(roughness + 1.0, 2.0) / 8.0;
  float kIBL = alpha2 / 2.0;
  float k = kdirect;
  float Gobstruction = dotNV / (dotNV * (1.0 - k) + k);
  float Gshadowing = dotNL / (dotNL * (1.0 - k) + k);
  return Gshadowing * Gobstruction;
}

vec3 diffuseLambert(vec3 albedo, vec3 kd)
{
  return kd * albedo / PI;
}

void
main()
{
  //normalized cam pos
  vec3 cameraPos = vec3(0.0, 0.0, 4.0);
  vec3 viewDir = normalize(cameraPos - normalize(worldPos));
  vec3 NormalWS = normalize(vNormalWS);

  vec3 kd = vec3(1.0);

  //temporary
  vec3 lightPos = vec3(0.0, 3.0, 3.0);
  vec3 lightColor = vec3(1.0, 1.0, 1.0);

  //incoming direction of light
  vec3 lightDir = normalize(lightPos - worldPos);

  //bissector of v and lightdir
  vec3 h = normalize(lightDir + viewDir);

  //Storing results
  float dotNV = max(abs(dot(NormalWS, viewDir)), 0.0);
  float dotNL = max(dot(NormalWS, lightDir), 0.0);
  float dotVN = max(dot(viewDir, NormalWS), 0.0);
  float dotNH = max(dot(NormalWS, h), 0.0);
  float dotLH = max(dot(lightDir, h), 0.0);
  float dotVH = max(dot(viewDir, h), 0.0);
  float alpha = roughness * roughness;
  float alpha2 = alpha * alpha;

  vec3 L = vec3(0.0);
    // **DO NOT** forget to do all your computation in linear space.
  //vec3 albedo = sRGBToLinear(vec4(albedo, 1.0)).rgb;
  vec3 albedo = texture(tex, uv).rgb;
  //for (int i = 0; i < 50; i++)
  {

  //Calculating Normal Distribution
  float nDistrib = distribGGX(dotNH, alpha2);

  //Calculate Schlick Fresnel approximation
  //Represents ks
  vec3 nFresnel = fresnelSchlick(dotLH, albedo); 

  //Calculate Smith GGX 
  float nGeometric = geometrySmith(dotNV, dotNL, alpha2);

  //Computing Cook-Torrance GGX model
  vec3 specular = (nDistrib * nFresnel * nGeometric) /
    (4.0 * dotNV * dotNL + E);


  //Computing diffuse Lambert
  kd = (kd - nFresnel) * (1.0 - metalness);
  vec3 diffuse = diffuseLambert(albedo, kd);


  float dist = length(lightPos - worldPos);
  float attenuation = 1.0 / (dist * dist);
  vec3 radiance = lightColor * attenuation;
  L += (diffuse + specular) *  dotNL * radiance * lightColor;
  }

  vec3 ambient = vec3(0.03) * albedo * a_occlusion;
  vec3 color = ambient + L;
  color = color / (color + vec3(1.0));
  color = pow(color, vec3(1.0 / 2.2));

  // **DO NOT** forget to apply gamma correction as last step.
  outFragColor.rgba = vec4(/*color*/vec3(uv.x, 0.0, 0.0), 1.0);
  //outFragColor.rgba = LinearTosRGB(vec4(diffuse, 1.0));
  //outFragColor.rgba = vec4(vec3(dotVH), 1.0);
  //outFragColor.rgba = vec4(vec3(dotNV), 1.0);
  //outFragColor.rgba = vec4(vec3(nFresnel), 1.0);
  //outFragColor.rgba = vec4(normalize(vec3(0.0) - worldPos), 1.0);
  //outFragColor.rgba = vec4(abs(NormalWS), 1.0);
  //outFragColor.rgba = vec4((viewDir), 1.0);
  //outFragColor.rgba = vec4((lightDir), 1.0);
  //outFragColor.rgba = vec4(abs(h), 1.0);
  //outFragColor.rgba = vec4(normalize(abs(worldPos)), 1.0);
}