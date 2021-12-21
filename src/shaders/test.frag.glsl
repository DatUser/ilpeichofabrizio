#version 450

//layout (local_size_x = 16, local_size_y = 16) in;

// Copyright (c) 2016, NVIDIA CORPORATION.  All rights reserved.
// Created by Tobias Zirr (KIT, NVIDIA) and Anton Kaplanyan (NVIDIA)
// Based on https://www.shadertoy.com/view/Xds3zN created by Inigo Quilez
// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

// For a less restrictive license, find the isolated shading code at:
// http://alphanew.net/releases/glints/shading.txt

// Example code for the paper "Real-time_ Rendering of Procedural Multiscale Materials",
// Tobias Zirr (NVIDIA / KIT), Anton Kaplanyan (NVIDIA),
// in ACM SIGGRAPH Symposium on Interactive 3D Graphics and Games, February 2016.
// More info on https://research.nvidia.com/publication/real-time_-rendering-procedural-multiscale-materials

float hash( float n ) { return fract(sin(mod(n, 3.14))*753.5453123); }
vec2 hash2( float n ) { return vec2(hash(n), hash(1.1 + n)); }

//----------------------------------------------------------------------



//----------------------------------------------------------------------

// heightfield
float mapLandscape(vec2 a) {
    float h = 0.0;
    h += pow(sin(2. + a.x / 37.) * sin(.5 + a.y / 31.), 2.);
    h += .5 * pow(sin(a.x / 11.) * sin(a.y / 13.), 2.);
    h += .2 * sin(1.1 + a.x / 3.1) * sin(.3 + a.y / 3.7);
    h *= min(length(a) / 10., 5.);
    h += .1 * sin(.9 + a.x / 1.7);
    h += .1 * sin(.4 + a.y / 1.4);
    h += .05 * sin(1.1 + a.y / .9);
    h += 15. * (1. - cos(a.x / 51.));
    return h;
}

//----------------------------------------------------------------------

// raycasting
vec2 castRayLandscape( in vec3 ro, in vec3 rd ) {
    float delt = 0.1;
    const float mint = 0.5;
    const float maxt = 90.0;
    float lh = 0.0;
    float ly = 0.0;
    float t = mint;
    for( int i = 0; i < 200; ++i )
    {
        if (t < maxt);
        else break;
        
        vec3  p = ro + rd*t;
        float h = mapLandscape( p.xz );
        if( p.y < h )
        {
            // interpolate the intersection distance
            return vec2(t - delt + delt*(lh-ly)/(p.y-ly-h+lh), 1.99);
        }
        // allow the error to be proportinal to the distance
        delt = max(0.1*t, delt);
        lh = h;
        ly = p.y;
        
        t += delt;
    }
    return vec2(maxt, -1.);
}
vec3 calcNormalLandscape( in vec3 pos ) {
    vec3 eps = vec3( 0.01, 0.0, 0.0 );
    vec3 a = eps.xyz;
    a.y = mapLandscape(pos.xz + a.xz) - mapLandscape(pos.xz - a.xz);
    a.xz *= 2.;
    vec3 b = eps.zyx;
    b.y = mapLandscape(pos.xz + b.xz) - mapLandscape(pos.xz - b.xz);
    b.xz *= 2.;
    return normalize( cross(b, a) );
}

vec2 castRay( in vec3 ro, in vec3 rd ) {
    vec2 resLandscape = castRayLandscape(ro, rd);
    //if (resLandscape.y >= 0. && (m < 0. || resLandscape.x < t))
    // res = resLandscape;
    
    return resLandscape;
}


//----------------------------------------------------------------------

// math
float compMax(vec2 v) { return max(v.x, v.y); }
float maxNrm(vec2 v) { return compMax(abs(v)); }
mat2 inverse2(mat2 m) {
    return mat2(m[1][1], -m[0][1], -m[1][0], m[0][0]) / (m[0][0] * m[1][1] - m[0][1] * m[1][0]);
}

float erfinv(float x) {
    float w, p;
    w = -log((1.0-x)*(1.0+x));
    if(w < 5.000000) {
        w = w - 2.500000;
        p = 2.81022636e-08;
        p = 3.43273939e-07 + p*w;
        p = -3.5233877e-06 + p*w;
        p = -4.39150654e-06 + p*w;
        p = 0.00021858087 + p*w;
        p = -0.00125372503 + p*w;
        p = -0.00417768164 + p*w;
        p = 0.246640727 + p*w;
        p = 1.50140941 + p*w;
    }
    else {
        w = sqrt(w) - 3.000000;
        p = -0.000200214257;
        p = 0.000100950558 + p*w;
        p = 0.00134934322 + p*w;
        p = -0.00367342844 + p*w;
        p = 0.00573950773 + p*w;
        p = -0.0076224613 + p*w;
        p = 0.00943887047 + p*w;
        p = 1.00167406 + p*w;
        p = 2.83297682 + p*w;
    }
    return p*x;
}

// ray differentials
void calcDpDxy( in vec3 ro, in vec3 rd, in vec3 rdx, in vec3 rdy, in float t, in vec3 nor, 
                out vec3 dpdx, out vec3 dpdy ) {
    dpdx = 2.*t*(rdx*dot(rd,nor)/dot(rdx,nor) - rd) * sign(dot(rd, rdx));
    dpdy = 2.*t*(rdy*dot(rd,nor)/dot(rdy,nor) - rd) * sign(dot(rd, rdy));
}

// some microfacet BSDF geometry factors
// (divided by NoL * NoV b/c cancelled out w/ microfacet BSDF)
float geometryFactor(float NoL, float NoV, vec2 roughness) {
    float a2 = roughness.x * roughness.y;
    NoL = abs(NoL);
    NoV = abs(NoV);

    float G_V = NoV + sqrt((NoV - NoV * a2) * NoV + a2);
    float G_L = NoL + sqrt((NoL - NoL * a2) * NoL + a2);
    return 1. / (G_V * G_L);
}

//----------------------------------------------------------------------

// ugly inefficient WebGL implementation of simple bit shifts for
// multilevel coherent grid indices. See comment in multilevelGridIdx.
int multilevelGridIdx1(inout int idx) {
    for (int i = 0; i < 32; ++i) {
        if (idx / 2 == (idx + 1) / 2)
          idx /= 2;
        else
            break;
    }
    return idx;
}
ivec2 multilevelGridIdx(ivec2 idx) {
//  return idx >> findLSB(idx); // findLSB not supported by Shadertoy WebGL version
    return ivec2(multilevelGridIdx1(idx.x), multilevelGridIdx1(idx.y));
}

//----------------------------------------------------------------------

// stable binomial 'random' numbers: interpolate between result for
// two closest binomial distributions where log_{.9}(p_i) integers
float binomial_interp(float u, float N, float p) {
    if(p >= 1.)
        return N;
    else if(p <= 1e-10)
        return 0.;

    // convert to distribution on ints while retaining expected value
    float cN = ceil(N);
    int iN = int(cN);
    p = p * (N / cN);
    N = cN;

    // round p to nearest powers of .9 (more stability)
    float pQ = .9;
    float pQef = log2(p) / log2(pQ);
    float p2 = exp2(floor(pQef) * log2(pQ));
    float p1 = p2 * pQ;
    vec2 ps = vec2(p1, p2);

    // compute the two corresponding binomials in parallel
    vec2 pm = pow(1. - ps, vec2(N));
    vec2 cp = pm;
    vec2 r = vec2(N);

    float i = 0.0;
    // this should actually be < N, no dynamic loops in ShaderToy right now
    for(int ii = 0; ii <= 17; ++ii)
    {
        if(u < cp.x)
            r.x = min(i, r.x);
        if(u < cp.y) {
            r.y = i;
            break;
        }
        // fast path
        if(ii > 16)
        {
            float C = 1. / (1. - pow(p, N - i - 1.));
            vec2 U = (u - cp) / (1. - cp);
            vec2 A = (i + 1. + log2(1. - U / C) / log2(p));
            r = min(A, r);
            break;
        }

        i += 1.;
        pm /= 1. - ps;
        pm *= (N + 1. - i) / i;
        pm *= ps;
        cp += pm;
    }

    // interpolate between the two binomials according to log p (akin to mip interpolation)
    return mix(r.y, r.x, fract(pQef));
}
// resort to gaussian distribution for larger N*p
float approx_binomial(float u, float N, float p) {
    if (p * N > 5.)
    {
        float e = N * p;
        float v = N * p * max(1. - p, 0.0);
        float std = sqrt(v);
        float k = e + erfinv(mix(-.999999, .999999, u)) * std;
        return min(max(k, 0.), N);
    }
    else
        return binomial_interp(u, N, p);
}

//----------------------------------------------------------------------

vec3 glints(vec2 texCO, vec2 duvdx, vec2 duvdy, mat3 ctf
  , vec3 lig, vec3 nor, vec3 view
  , vec2 roughness, vec2 microRoughness, float searchConeAngle, float variation, float dynamicRange, float density)
{
   vec3 col = vec3(0.);

    // Compute pixel footprint in texture space, step size w.r.t. anisotropy of the footprint
    mat2 uvToPx = inverse2(mat2(duvdx, duvdy));
    vec2 uvPP = 1. / vec2(maxNrm(uvToPx[0]), maxNrm(uvToPx[1]));

    // material
    vec2 mesoRoughness = sqrt(max(roughness * roughness - microRoughness * microRoughness, vec2(1.e-12))); // optimizer fail, max 0 removed

    // Anisotropic compression of the grid
    vec2 texAnisotropy = vec2( min(mesoRoughness.x / mesoRoughness.y, 1.)
                             , min(mesoRoughness.y / mesoRoughness.x, 1.) );

    // Compute half vector (w.r.t. dir light)
    vec3 hvW = normalize(lig + view);
    vec3 hv = normalize(hvW * ctf);
    vec2 h = hv.xy / hv.z;
    vec2 h2 = 0.75 * hv.xy / (hv.z + 1.);
    // Anisotropic compression of the slope-domain grid
    h2 *= texAnisotropy;

    // Compute the Gaussian probability of encountering a glint within a given finite cone
    vec2 hppRScaled = h / roughness;
    float pmf = (microRoughness.x * microRoughness.y) / (roughness.x * roughness.y)
        * exp(-dot(hppRScaled, hppRScaled)); // planeplane h
    pmf /= hv.z * hv.z * hv.z * hv.z; // projected h
//  pmf /= dot(lig, nor) * dot(view, nor); // projected area, cancelled out by parts of G, ...
    float pmfToBRDF = 1. / (3.14159 * microRoughness.x * microRoughness.y);
    pmfToBRDF /= 4.; // solid angle o
    pmfToBRDF *= geometryFactor(dot(lig, nor), dot(view, nor), roughness); // ... see "geometryFactor"
    // phenomenological: larger cones flatten distribution
    float searchAreaProj = searchConeAngle * searchConeAngle / (4. * dot(lig, hvW) * hv.z); // * PI
    pmf = mix(pmf, 1., clamp(searchAreaProj, 0.0, 1.0)); // searchAreaProj / PI
    pmf = min(pmf, 1.);
    
    // noise coordinate (decorrelate interleaved grid)
    texCO += vec2(100.);
    // apply anisotropy _after_ footprint estimation
    texCO *= texAnisotropy;

    // Compute AABB of pixel in texture space
    vec2 uvAACB = max(abs(duvdx), abs(duvdy)) * texAnisotropy; // border center box
    vec2 uvb = texCO - 0.5 * uvAACB;
    vec2 uve = texCO + 0.5 * uvAACB;

    vec2 uvLongAxis = uvAACB.x > uvAACB.y ? vec2(1.0, 0.0) : vec2(0.0, 1.0);
    vec2 uvShortAxis = 1.0 - uvLongAxis;

    // Compute skew correction to snap axis-aligned line sampling back to longer anisotropic pixel axis in texture space
    vec2 skewCorr2 = -(uvToPx * uvLongAxis) / (uvToPx * uvShortAxis);
    float skewCorr = abs((uvToPx * uvShortAxis).x) > abs((uvToPx * uvShortAxis).y) ? skewCorr2.x : skewCorr2.y;
    skewCorr *= dot(texAnisotropy, uvShortAxis) / dot(texAnisotropy, uvLongAxis);

    float isoUVPP = dot(uvPP, uvShortAxis);
    // limit anisotropy
    isoUVPP = max(isoUVPP, dot(uvAACB, uvLongAxis) / 16.0);

     // Two virtual grid mips: current and next
    float fracMip = log2(isoUVPP);
    float lowerMip = floor(fracMip);
    float uvPerLowerC = exp2(lowerMip);

    // Current mip level and cell size
    float uvPC = uvPerLowerC;
    float mip = lowerMip;

    int iter = 0;
    int iterThreshold = 60;

    for (int i = 0; i < 2; ++i)
    {
        float mipWeight = 1.0 - abs(mip - fracMip);

        vec2 uvbg = min(uvb + 0.5 * uvPC, texCO);
        vec2 uveg = max(uve - 0.5 * uvPC, texCO);

        // Snapped uvs of the cell centers
        vec2 uvbi = floor(uvbg / uvPC);
        vec2 uvbs = uvbi * uvPC;
        vec2 uveo = uveg + uvPC - uvbs;

        // Resulting compositing values for a current layer
        float weight = 0.0;
        vec3 reflection = vec3(0.0);

        // March along the long axis
        vec2 uvo = vec2(0.0), uv = uvbs, uvio = vec2(0.0), uvi = uvbi;
        for (int iter1 = 0; iter1 < 18; ++iter1) // horrible WebGL-compatible static for loop
        {
            // for cond:
            if (dot(uvo, uvLongAxis) < dot(uveo, uvLongAxis) && iter < iterThreshold);
            else break;

            // Snap samples to long anisotropic pixel axis
            float uvShortCenter = dot(texCO, uvShortAxis) + skewCorr * dot(uv - texCO, uvLongAxis);

            // Snapped uvs of the cell center
            uvi += (floor(uvShortCenter / uvPC) - dot(uvi, uvShortAxis)) * uvShortAxis;
            uv = uvi * uvPC;
            float uvShortEnd = uvShortCenter + uvPC;

            vec2 uvb2 = uvbg * uvLongAxis + uvShortCenter * uvShortAxis;
            vec2 uve2 = uveg * uvLongAxis + uvShortCenter * uvShortAxis;

            // March along the shorter axis
            for (int iter2 = 0; iter2 < 4; ++iter2) // horrible WebGL-compatible static for loop
            {
                // for cond:
                if (dot(uv, uvShortAxis) < uvShortEnd && iter < iterThreshold);
                else break;

                // Compute interleaved cell index
                ivec2 cellIdx = ivec2(uvi + vec2(.5));
                cellIdx = multilevelGridIdx(cellIdx);

                // Randomize a glint based on a texture-space id of current grid cell
                vec2 u2 = hash2(float( (cellIdx.x + 1549 * cellIdx.y) ));
                // Compute index of the cone
                vec2 hg = h2 / (microRoughness + searchConeAngle);
                vec2 hs = floor(hg + u2) + u2 * 533.;    // discrete cone index in paraboloid hv grid
                ivec2 coneIdx = ivec2(hs);

                // Randomize glint sizes within this layer
                float var_u = hash(float( (cellIdx.x + cellIdx.y * 763 + coneIdx.x + coneIdx.y * 577) ));
                float mls = 1. + variation * erfinv(mix(-.999, .999, var_u));
                if (mls <= 0.0) mls = fract(mls) / (1. - mls);
                mls = max(mls, 1.e-12);

                // Bilinear interpolation using coverage made by areas of two rects
                vec2 mino = max(1.0 - max((uvb2 - uv) / uvPC, 0.0), 0.0);
                vec2 maxo = max(1.0 - max((uv - uve2) / uvPC, 0.0), 0.0);
                vec2 multo = mino * maxo;
                float coverageWeight = multo.x * multo.y;

                float cellArea = uvPC * uvPC;
                // Expected number of glints 
                float eN = density * cellArea;
                float sN = max(eN * mls, min(1.0, eN));
                eN = eN * mls;

                // Sample actually found number of glints
                float u = hash(float(coneIdx.x + coneIdx.y * 697));
                float lN = approx_binomial(u, sN, pmf);
#if 0
                // Colored glints
                if (false) {
                    vec3 glintColor = hue_colormap(fract(u + u2.y));
                    glintColor = mix(vec3(1.0f), glintColor, 1. / sqrt(max(lN, 1.0)));
                }
#endif
                // Ratio of glinting vs. expected number of microfacets
                float ratio = lN / eN;
                
                // limit dynamic range (snow more or less unlimited)
                ratio = min(ratio, dynamicRange * pmf);
                
                // convert to reflectance
                ratio *= pmfToBRDF;
#if 0
                // Grid
                reflection += vec3(u);
                weight += coverageWeight;
#else
                // Accumulate results
                reflection += coverageWeight * ratio;
                weight += coverageWeight;
#endif

                // for incr:
                uv += uvPC * uvShortAxis, uvi += uvShortAxis, ++iter;
            }

            // for incr:
              uvo += uvPC * uvLongAxis, uv = uvbs + uvo
            , uvio += uvLongAxis, uvi = uvbi + uvio;
        }

#ifdef DEBUG
        // Normalization
        if (weight < 1.e-15) {
            col = vec3(0.0, 1.0, 1.0);
            break;
        }
#endif

        reflection = reflection / weight;

        // Compositing of two layers
        col += mipWeight * reflection;

        // for incr:
        uvPC *= 2., mip += 1.;
    }

    return col;
}

vec3 render( in vec3 ro, in vec3 rd, in vec3 rdx, in vec3 rdy )
{
    // sun and sky
    vec3 lig = normalize( vec3(0.6, .9, 0.5) );
    vec3 lightPower = vec3(9.);
    vec3 sky = vec3(0.7, 0.9, 1.0) + 1. + rd.y*0.8;
    vec3 col = sky * lightPower;
    // ray cast
    vec2 res = castRay(ro,rd);
    float t = res.x;
    float m = res.y;
    // shade hit
    if( m>-0.5 ) {
        // hit information
        vec3 pos = ro + t*rd;
        //vec3 nor = (m < 2.) ? calcNormalLandscape(pos) : calcNormal( pos );
        vec3 nor = calcNormalLandscape(pos);
        
        mat3 texProjFrame = mat3( vec3(1,0,0), vec3(0,0,1), vec3(0,1,0) );
        if (abs(nor.x) > abs(nor.y) && abs(nor.x) > abs(nor.z)) {
            texProjFrame = mat3( vec3(0,0,1), vec3(0,1,0), vec3(1,0,0) );
        }
        else if (abs(nor.z) > abs(nor.x) && abs(nor.z) > abs(nor.y)) {
            texProjFrame = mat3( vec3(1,0,0), vec3(0,1,0), vec3(0,0,1) );
        }
        
        vec3 bitang = normalize(cross(nor, texProjFrame[0]));
        vec3 tang = cross(bitang, nor);
        mat3 ctf = mat3(tang, bitang, nor);
        
        // texturing
        vec3 dposdx, dposdy;
        calcDpDxy( ro, rd, rdx, rdy, t, nor, dposdx, dposdy );
        // planar projections
        float texScaling = 1.;
        vec2 texCO = texScaling * (pos * texProjFrame).xy;
        vec2 duvdx = texScaling * (dposdx * texProjFrame).xy
           , duvdy = texScaling * (dposdy * texProjFrame).xy;
        // computing these manually from ray differentials to handle edges in ray casting,
        // can simply use standard derivatives in real-world per-object fragment shaders:
//      duvdx = dFdx(texCO), duvdy = dFdy(texCO);
        
        // useful information
        float occ = 0.5;
        //softshadow( pos, lig, 0.02, 25. );
        float amb = clamp( 0.5+0.5*nor.y, 0.0, 1.0 );
        float dif = clamp( dot( nor, lig ), 0.0, 1.0 );
        float fre = 1. - pow(1. - dif, 2.5);
        float dfr = 1. - pow(1. - clamp( dot( nor, -rd ), 0.0, 1.0 ), 2.5);
        dfr *= fre;
        
        // some default material        
        col = 0.45 + 0.3*sin( vec3(0.05,0.08,0.10)*(m-1.0) );
        col *= .4;
        float specularity = fract(m);
        
        // configure multiscale material (snow parameters)
        vec2 roughness = vec2(.6);
        vec2 microRoughness = roughness * .024;
        float searchConeAngle = .01;
        float variation = 100.;
        float dynamicRange = 50000.;
        float density = 5.e8;
    
        col = mix(vec3(.1,.2,.5), vec3(.95,.8,.75), (1. - abs(rd.y)) * dif);
        
        // standard diffuse lighting
        col *= lightPower * mix(.02, 1., occ * dif);
        
        // multiscale specular lighting
        if (specularity > 0.0 && dif > 0.0 && dot(-rd, nor) > 0.0)
			col += specularity * glints(texCO, duvdx, duvdy, ctf, lig, nor, -rd, roughness, microRoughness, searchConeAngle, variation, dynamicRange, density)
				* lightPower * mix(.05, 1., occ);
    }

  return col;
}

//----------------------------------------------------------------------

mat3 setCamera( in vec3 ro, in vec3 ta, float cr )
{
    vec3 cw = normalize(ta-ro);
    vec3 cp = vec3(sin(cr), cos(cr),0.0);
    vec3 cu = normalize( cross(cw,cp) );
    vec3 cv = normalize( cross(cu,cw) );
    return mat3( cu, cv, cw );
}

vec2 iResolution  = vec2(1024, 1024);
//uniform vec4 iMouse;
uniform float time;

out vec4 fragColor;

void main()
{
    vec2 fragCoord = gl_FragCoord.xy;
    vec2 q = fragCoord.xy/iResolution.xy;
    vec2 p = -1.0+2.0*q;
    p.x *= iResolution.x/iResolution.y;
    //vec2 mo = iMouse.xy/iResolution.xy;
     
    float time_ = 15.0 + time;

    // camera
    vec2 mo = vec2(0.0);
    //float time_ = 1505.0;
    float ds = 1.5 + sin(time_ / 2.);
    vec3 ro = vec3( -0.5+ds * 8.5*cos(0.1*time_ + 6.0*mo.x)
                   , 10.0 - 9.5*mo.y
                   , 0.5 + ds * 8.5*sin(0.1*time_ + 6.0*mo.x) );
    ro.y /= 1. + .01 * dot(ro.xz, ro.xz);
    ro.y += mapLandscape(ro.xz);
    vec3 ta = vec3( -0.5, -0.4, 0.5 );
  
    // camera-to-world transformation
    mat3 ca = setCamera(ro, ta, 0.0 );
    
    // ray direction
    vec3 rd = ca * normalize( vec3(p.xy,2.0) );
    vec2 rds = -sign(p + .001);
    vec3 rdx = ca * rds.x * normalize( vec3(p.xy + rds.x * vec2(1./iResolution.y,0),2.0) );
    vec3 rdy = ca * rds.y * normalize( vec3(p.xy + rds.y * vec2(0,1./iResolution.y),2.0) );

    // render 
    vec3 col = render(ro, rd, rdx, rdy );
    // tonemap, gamma
    col *= 1.0 / (max(max(col.r, col.g), col.b) + 1.0);
    col = pow( col, vec3(0.4545) );
    fragColor = vec4( col, 1.0 );
    //fragColor = vec4( vec3(1.0), 1.0 );
}