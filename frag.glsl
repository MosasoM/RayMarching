precision mediump float;
uniform float time;
uniform vec2 resolution;
uniform vec2 mouse;

struct Object{
    vec3 pos;
    vec3 rot;
    int kind;
    float[10] params;
};

struct Camera{
    vec3 pos;
    float fov;
};

struct transform{
    vec3 pos;
    vec3 rot;
};

struct intersection{
    bool hit;
    vec3 p;
    float distance;
    vec3 normal;

};

vec3 calc_ray(Camera cam,float x,float z);
float distance_func(Object obj,vec3 rayhead);
vec3 calc_norm(Object obj,vec3 hitpos);

float atan2(float y, float x);
float df_sphere(vec3 pos, float size);
vec3 sp_norm (vec3 pos,float size);
float df_box(vec3 box,vec3 size);
vec3 box_norm(vec3 pos,vec3 size);
float df_torus(vec3 p, vec2 t);
vec3 torus_norm(vec3 p, vec2 t);
float df_cylinder(vec3 pos, vec3 size);
vec3 sylinder_norm(vec3 pos, vec3 size);
float df_prism(vec3 p,vec2 h);
vec3 prism_norm(vec3 p,vec2 h);
float df_octa(vec3 p,float s);
vec3 octa_norm(vec3 p, float t);

vec3 microBRDF(float a,float dotNH, float dotNV, float dotNL, float dotVH,vec3 speccolor);
float D_GGX(float a, float dotNH);
float G_Smith_Schlick_GGX(float a, float dotNV, float dotNL);
vec3 F_Schlick(vec3 specularColor, float dotVH);
vec3 rept(vec3 pos,float freq);
vec3 rgbnormalize(vec3 rgb);
vec3 rot_by_Rodrigues(vec3 p, vec3 axis, float angle);
vec3 normed_phong(vec3 speccolor,float power,vec3 view, vec3 norm, vec3 lightDir);


const float PI = 3.1415926535;
vec2 oo = vec2(0.0,0.0);
const int marching_max = 128;
// const float eps = 0.001;


void main(){

vec2 st = (gl_FragCoord.xy-resolution.xy)/min(resolution.x,resolution.y);
float st_x = st.x;
float st_z = st.y;

Camera cam;
cam.pos = vec3(0.0,0.0,0.0);
cam.fov = (PI*30.0)/(2.0*180.0);

Object sphere;
sphere.pos = vec3(2.5,15.0,2.5);
sphere.rot = vec3(0.0,0.0,0.0);
sphere.kind = 1;
sphere.params[0] = 1.0;

vec3 ray;
ray = calc_ray(cam,st_x,st_z);


float max_dis = 100.0;
float min_dis = 0.01;
float rlen = 0.0;
bool hitflag = false;
const int max_loop = 100;

for (int i = 0; i < max_loop; ++i){
    float d = distance_func(sphere,cam.pos+ray*rlen);
    if (rlen > max_dis){
        break;
    }
    if (abs(d) < min_dis){
        hitflag = true;
        break;
    }else{
        rlen += d;
    }
}

if (hitflag){
    gl_FragColor = vec4(1.0,1.0,1.0,1.0);
}else{
    gl_FragColor = vec4(0.0, 0.0 ,0.0, 1.0);
}



// float dist = 10000000.0; // length to next object
// float rlen = 0.0; //additional length to ray pos
// vec3 rpos = cpos; //position of the tip of ray


// for (int i = 0; i < marching_max; ++i){
//     for (int j = 0; j < 8; ++ j){
//         float d = df_sphere(rpos-sp_all[j],0.2);
//         dist = min(dist,d);
//     }
//     rlen += dist;
//     rpos = cpos + ray*rlen;
// }

// bool flag = true;
// for (int i = 0; i < 8; ++ i){
//         float d = df_sphere(rpos-sp_all[i],0.2);
//         if (d<eps){
//             flag = false;
//             break;
//         }
// }

// if (flag){     
//     gl_FragColor = vec4(vec3(0.7), 1.0);
// }else{
//     for (int i = 0; i < 8; ++ i){
//         float d = df_sphere(rpos-sp_all[i],0.2);
//         if (d<eps){
//             vec3 temp = vec3(abs(sp_all[i].x),abs(sp_all[i].y),abs(sp_all[i].z));
//             gl_FragColor = vec4((normalize(sp_all[i])+0.2)*sin(time/800.0), 1.0);
//             break;
//         }
//     }
// }

}

vec3 calc_ray(Camera cam, float x, float z){
    return normalize(vec3(sin(cam.fov)*x ,cos(cam.fov),sin(cam.fov)*z));
}

float distance_func(Object obj,vec3 rayhead){
    vec3 p = rayhead-obj.pos;
    if (obj.kind == 1){
        return length(p)-obj.params[0];
    }
}

vec3 calc_norm(Object obj,vec3 hitpos){
    float eps = 0.001;
    return normalize(
        vec3(
            distance_func(obj,hitpos+vec3(eps,0.0,0.0))-distance_func(obj,hitpos+vec3(-eps,0.0,0.0)),
            distance_func(obj,hitpos+vec3(0.0,eps,0.0))-distance_func(obj,hitpos+vec3(0.0,-eps,0.0)),
            distance_func(obj,hitpos+vec3(0.0,0.0,eps))-distance_func(obj,hitpos+vec3(0.0,0.0,-eps))
        )
    );
}



// float atan2(float y, float x){
//     return x == 0.0 ? sign(y)*pi/2.0 : atan(y, x);
// }

// float df_sphere(vec3 pos,float size){
//     return length(pos) - size;
// }

// vec3 sp_norm(vec3 pos,float size){
//     return normalize(
//         vec3(
//             df_sphere(pos+vec3(eps,0.0,0.0),size)-df_sphere(pos+vec3(-eps,0.0,0.0),size),
//             df_sphere(pos+vec3(0.0,eps,0.0),size)-df_sphere(pos+vec3(0.0,-eps,0.0),size),
//             df_sphere(pos+vec3(0.0,0.0,eps),size)-df_sphere(pos+vec3(0.0,0.0,-eps),size)
//         )
//     );
// }

// float df_box(vec3 pos, vec3 size){
//     vec3 d = abs(pos)-size;
//     return length(max(d,0.0))+ min(max(d.x,max(d.y,d.z)),0.0);
// }

// vec3 box_norm(vec3 pos,vec3 size){
//         return normalize(
//         vec3(
//             df_box(pos+vec3(eps,0.0,0.0),size)-df_box(pos+vec3(-eps,0.0,0.0),size),
//             df_box(pos+vec3(0.0,eps,0.0),size)-df_box(pos+vec3(0.0,-eps,0.0),size),
//             df_box(pos+vec3(0.0,0.0,eps),size)-df_box(pos+vec3(0.0,0.0,-eps),size)
//         )
//     );
// }

// float df_torus (vec3 p, vec2 t){
//     vec2 q = vec2(length(p.xz)-t.x,p.y);
//     return length(q)-t.y;
// }

// vec3 torus_norm (vec3 p, vec2 t){
//         return normalize(
//         vec3(
//             df_torus(p+vec3(eps,0.0,0.0),t)-df_torus(p+vec3(-eps,0.0,0.0),t),
//             df_torus(p+vec3(0.0,eps,0.0),t)-df_torus(p+vec3(0.0,-eps,0.0),t),
//             df_torus(p+vec3(0.0,0.0,eps),t)-df_torus(p+vec3(0.0,0.0,-eps),t)
//         )
//     );
// }

// float df_prism(vec3 p,vec2 h){
//     vec3 q = abs(p);
//     return max(q.z-h.y,max(q.x*0.866025+p.y*0.5,-p.y)-h.x*0.5);
// }

// vec3 prism_norm(vec3 p,vec2 h){
//          return normalize(
//         vec3(
//             df_prism(p+vec3(eps,0.0,0.0),h)-df_prism(p+vec3(-eps,0.0,0.0),h),
//             df_prism(p+vec3(0.0,eps,0.0),h)-df_prism(p+vec3(0.0,-eps,0.0),h),
//             df_prism(p+vec3(0.0,0.0,eps),h)-df_prism(p+vec3(0.0,0.0,-eps),h)
//         )
//     );
// }

// float df_octa(vec3 p,float s){
//     p = abs(p);
//     float m = p.x+p.y+p.z-s;
//     vec3 q;
//          if( 3.0*p.x < m ) q = p.xyz;
//     else if( 3.0*p.y < m ) q = p.yzx;
//     else if( 3.0*p.z < m ) q = p.zxy;
//     else return m*0.57735027;
    
//     float k = clamp(0.5*(q.z-q.y+s),0.0,s); 
//     return length(vec3(q.x,q.y-s+k,q.z-k)); 
// }

// vec3 octa_norm(vec3 p, float t){
//         return normalize(
//         vec3(
//             df_octa(p+vec3(eps,0.0,0.0),t)-df_octa(p+vec3(-eps,0.0,0.0),t),
//             df_octa(p+vec3(0.0,eps,0.0),t)-df_octa(p+vec3(0.0,-eps,0.0),t),
//             df_octa(p+vec3(0.0,0.0,eps),t)-df_octa(p+vec3(0.0,0.0,-eps),t)
//         )
//     );
// }

// vec3 microBRDF(float a,float dotNH, float dotNV, float dotNL, float dotVH,vec3 speccolor){
//     float D = D_GGX(a,dotNH);
//     float G = G_Smith_Schlick_GGX(a,dotNV,dotNL);
//     vec3 F = F_Schlick(speccolor,dotVH);

//     return ((F*D*G)/(4.0*dotNV*dotNV+0.0001));
    
// }


// float D_GGX(float a, float dotNH) {
//   float a2 = a*a;
//   float dotNH2 = dotNH*dotNH;
//   float d = dotNH2 * (a2 - 1.0) + 1.0;
//   return a2 / (pi * d * d);
// }

// float G_Smith_Schlick_GGX(float a, float dotNV, float dotNL) {
//   float k = a*a*0.5 + 0.0001;
//   float gl = dotNL / (dotNL * (1.0 - k) + k);
//   float gv = dotNV / (dotNV * (1.0 - k) + k);
//   return gl*gv;
// }

// vec3 F_Schlick(vec3 specularColor, float dotVH) {
//   return specularColor + ((1.0 - specularColor) * pow((1.0 - dotVH), 5.0));
// }

// vec3 rept (vec3 pos, float freq){
//     return mod(pos,freq)-freq/2.0;
// }

// vec3 rgbnormalize(vec3 rgb){
//     return vec3(rgb/255.0);
// }

// vec3 rot_by_Rodrigues(vec3 p, vec3 axis, float angle){
//     vec3 a = normalize(axis);
//     float theta = angle*pi/180.0;
//     float c = cos(theta);
//     float s = sin(theta);
//     float r = 1.0-c;

//     mat3 m = mat3(
//         a.x * a.x * r + c,        a.y * a.x * r + a.z * s,  a.z * a.x * r - a.y * s,
//         a.x * a.y * r - a.z * s,  a.y * a.y * r + c,        a.z * a.y * r + a.x * s,
//         a.x * a.z * r + a.y * s,  a.y * a.z * r - a.x * s,  a.z * a.z * r + c
//     );
//     return m*p;
// }

// vec3 normed_phong(vec3 speccolor,float power,vec3 view, vec3 norm, vec3 lightDir){
//      vec3 R = -view + ( 2.0 * dot( norm, view ) * norm );
//      return speccolor * pow( max (dot(lightDir,R), 0.0 ), power ) * ( ( power + 1.0 )/ ( 2.0 * pi ) );
// }
