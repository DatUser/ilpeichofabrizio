#version 450

precision highp float;

out vec4 outFragColor;

in vec3 worldPos;
in vec3 vNormalWS;
in vec3 vviewDir;

vec3 albedo = vec3(1.0, 1.0, 1.0);

const float PI = 3.1415926535897932384626433832795;
float metalness = 0.0;//1 if metallic 0 otherwise

float isPos(float nb)
{
  return (nb > 0.0) ? 1.0 : 0.0;
}

void
main()
{
  vec3 ks = vec3(1.0);
  vec3 kd = vec3(1.0);

  vec3 NormalWS = normalize(vNormalWS);
  vec3 viewDir = normalize(vviewDir);


  //temporary
  vec3 lightPos = vec3(2.0, 4.0, 4.0);

  //incoming direction of light
  vec3 lightDir = normalize(lightPos - worldPos);

  //outcoming direction of light
  //vec3 v = reflect(lightDir, vNormalWS);
  //bissector of v and lightdir
  vec3 h = normalize(lightDir + viewDir);

  //Storing results
  float dotNV = max(dot(NormalWS, viewDir), 0.0);
  float dotNL = max(dot(NormalWS, lightDir), 0.0);
  float dotVL = max(dot(viewDir, lightDir), 0.0);
  float dotVH = max(dot(viewDir, h), 0.0);;
  //float dotNH = dot(vNormalWS, h);

  //Calculating Normal Distribution
  //roughness
  float roughness = 0.5;
  float alpha = roughness;//pow(roughness, 2.0);
  float alpha2 = pow(alpha, 2.0);

  //halfway vector between light dir and normal
  //vec3 h_distrib = normalize(lightDir + viewDir);

  float dotNH = max(dot(NormalWS, h), 0.0);
  float dotNH2 = pow(dotNH, 2.0);
  float nDistrib = alpha2 //* isPos(viewLight)
  / (PI * pow(dotNH2 * (alpha2 - 1.0) + 1.0, 2.0));

  //Calculate Fresnel
  vec3 f0 = vec3(0.04);//Pre-computed (here we use water value)
  f0 = mix(f0, albedo, metalness);

  float f90 = 1.0;//Pre-computed (here we use water value)
  vec3 nFresnel = f0 + (f90 - f0) * pow(1.0 - dotNV, 5.0);
  ks = nFresnel;

  //Calculate Schlick Geometric Function
  float kdirect = pow(roughness + 1.0, 2.0) / 8.0;
  float kIBL = alpha2 / 2.0;
  float k = kdirect;//Pre-computed
  float Gobstruction = dotNV / (dotNV * (1.0 - k) + k);
  float Gshadowing = dotNL / (dotNL * (1.0 - k) + k);
  float nGeometric = Gobstruction * Gshadowing;

  //Computing Cook-Torrance GGX model
  vec3 specular = (nDistrib * nFresnel * nGeometric) /
    (4.0 * dotNV * dotNL);

  //Computing diffuse Lambert
  //kd = (kd - nFresnel) * (1.0 - metalness);
  vec3 lightColor = vec3(1.0, 1.0, 1.0);
  float dist = length(lightPos - worldPos);
  //float attenuation = 1.0 / (dist * dist);
  //vec3 radiance = lightColor * attenuation;// * dotNL;

  kd = (1.0 - ks);
  vec3 diffuse = kd * albedo / PI;


  vec3 L = (diffuse + specular) *  dotNL * lightColor;// * radiance;

  //vec3 ambient = vec3(0.03) * albedo;
  //vec3 color = L;// +ambient

  // **DO NOT** forget to apply gamma correction as last step.
  outFragColor.rgba = vec4(L, 1.0);
}
