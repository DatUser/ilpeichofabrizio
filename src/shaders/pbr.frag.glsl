#version 450

precision highp float;

out vec4 outFragColor;

uniform vec3 albedo;
uniform float roughness;
uniform float metalness;// = 0.1;//1 if metallic 0 otherwise

vec4 sRGBToLinear( in vec4 value ) {
	return vec4( mix( pow( value.rgb * 0.9478672986 + vec3( 0.0521327014 ), vec3( 2.4 ) ), value.rgb * 0.0773993808, vec3( lessThanEqual( value.rgb, vec3( 0.04045 ) ) ) ), value.a );
}

vec4 LinearTosRGB( in vec4 value ) {
	return vec4( mix( pow( value.rgb, vec3( 0.41666 ) ) * 1.055 - vec3( 0.055 ), value.rgb * 12.92, vec3( lessThanEqual( value.rgb, vec3( 0.0031308 ) ) ) ), value.a );
}

in vec3 worldPos;
in vec3 vNormalWS;
in vec3 vviewDir;

const float PI = 3.1415926535897932384626433832795;

float isPos(float nb)
{
  return (nb > 0.0) ? 1.0 : 0.0;
}

vec3 RGBMToLinear(vec4 value)
{
  return value.xyz * value.w * 6.0;
}

void
main()
{
  vec3 viewDir = normalize(vviewDir);
  vec3 vNormalWS = normalize(vNormalWS);
  //vec2 pos2D = sphericalTo2D(vNormalWS);
  //vec2 pos2D_spec = sphericalTo2D(reflect(viewDir));

  //float ks = 1.0;
  vec3 kd = vec3(1.0);

  //temporary
  vec3 lightPos = vec3(1.0, 1.0, 3.0);

  //incoming direction of light
  vec3 lightDir = normalize(lightPos - worldPos);

  //outcoming direction of light
  //vec3 v = reflect(lightDir, vNormalWS);
  //bissector of v and lightdir
  vec3 h = normalize(lightDir + viewDir);

  //Storing results
  float dotNV = max(dot(vNormalWS, viewDir), 0.0);
  float dotNL = max(dot(vNormalWS, lightDir), 0.0);
  float dotVN = max(dot(viewDir, vNormalWS), 0.0);
  //float dotHV = dot(h, v);
  //float dotNH = dot(vNormalWS, h);

    // **DO NOT** forget to do all your computation in linear space.
  vec3 albedo = sRGBToLinear(vec4(albedo, 1.0)).rgb;
  //vec3 diffEnv = RGBMToLinear(texture(myTexture, pos2D));
  //vec3 specEnv = RGBMToLinear(texture(texSpec, pos2D));

  //Calculating Normal Distribution
  //roughness
  //float roughness = 0.1;
  float alpha = roughness;//pow(roughness, 2.0);
  float alpha2 = pow(alpha, 2.0);

  //halfway vector between light dir and normal
  //vec3 h_distrib = normalize(lightDir + viewDir);

  float dotNH = max(dot(vNormalWS, h), 0.0);
  float dotNH2 = pow(dotNH, 2.0);
  float bot = dotNH2 * (alpha2 - 1.0) + 1.0;
  float nDistrib = alpha2 //* isPos(viewLight)
  / (PI * pow(dotNH2 * (alpha2 - 1.0) + 1.0, 2.0));

  //Calculate Fresnel
  vec3 f0 = vec3(0.04);//Pre-computed (here we use water value)
  f0 = mix(f0, albedo, metalness);

  float f90 = 1.0;//Pre-computed (here we use water value)
  vec3 nFresnel = f0 + (f90 - f0) * pow(1.0 - dotNH, 5.0);

  //Calculate Schlick GGX 
  float kdirect = pow(roughness + 1.0, 2.0) / 8.0;
  float kIBL = alpha2 / 2.0;
  float k = kIBL;//kdirect;//Pre-computed
  float Gobstruction = dotNV / (dotNV * (1.0 - k) + k);
  float Gshadowing = dotNL / (dotNL * (1.0 - k) + k);
  float nGeometric = Gobstruction * Gshadowing;
  //use 0.0001 as epsilon to avoid invalid operation

  //Computing Cook-Torrance GGX model
  vec3 specular = (nDistrib * nFresnel * nGeometric) /
    (4.0 * dotNV * dotNL);


  //Compute specular IBL
  vec3 R = 2.0 * dotVN * vNormalWS - viewDir;
  //vec4 prefColor = texture(texSpec, vec2(roughness, dotNV));
  //specular = specEnv * viewDir;

  //Computing diffuse Lambert
  kd = (kd - nFresnel) * (1.0 - metalness);
  vec3 lightColor = vec3(1.0, 1.0, 1.0);
  float dist = length(lightPos - worldPos);
  float attenuation = 1.0 / (dist * dist);

  //Calculating Fresnel for env
  vec3 nFresnelEnv = f0 + (f90 - f0) * pow(1.0 - dotNH, 5.0);

  //Calculating kd for env
  vec3 kdEnv = (vec3(1.0) - nFresnelEnv) * (1.0 - metalness);

  vec3 diffuse = kd * albedo / PI;// + diffEnv * kdEnv;


  vec3 L = (diffuse + specular) *  dotNL * lightColor;
  //L += diffEnv * kdEnv;

  // **DO NOT** forget to apply gamma correction as last step.
  outFragColor.rgba = LinearTosRGB(vec4(L, 1.0));
}