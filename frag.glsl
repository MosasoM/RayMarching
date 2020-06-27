precision mediump float;
uniform float time;
uniform vec2 resolution;
uniform vec2 mouse;

#define SIZE_OF_OBJS_ARRAY 100

struct Material{
    vec3 albedo;
    vec3 f0;
    float roughness;
    bool glow;
    int kind;
};

struct Light{
    vec3 power;
    vec3 pos;
    vec3 rot;
    int kind;
};

struct Object{
    vec3 pos;
    vec3 rot;
    int kind;
    float[10] params;
    Material material;
};

struct Camera{
    vec3 pos;
    float fov;
};


vec3 calc_ray(in Camera cam,in float x,in float z);
float distance_func(in Object obj,in vec3 rayhead);
vec3 calc_norm(in Object obj,in vec3 hitpos);
vec3 material_color(in Material mat);
vec3 calc_light(in Light light,in vec3 hitpos);

void raymarching(in vec3 origin,in vec3 ray,in Object[SIZE_OF_OBJS_ARRAY] objs,inout bool hitflag,inout int hitnum);
//inonutで渡さないと参照代入したいやつに関数内で代入が行われないと、mainで代入した値じゃなく各型ごとの初期値が勝手に代入されてバグる。



vec3 microBRDF(float a,float dotNH, float dotNV, float dotNL, float dotVH,vec3 speccolor);
float D_GGX(float a, float dotNH);
float G_Smith_Schlick_GGX(float a, float dotNV, float dotNL);
vec3 F_Schlick(vec3 specularColor, float dotVH);
vec3 rept(vec3 pos,float freq);
vec3 rgbnormalize(vec3 rgb);
vec3 rot_by_Rodrigues(vec3 p, vec3 axis, float angle);
vec3 normed_phong(vec3 speccolor,float power,vec3 view, vec3 norm, vec3 lightDir);


const float PI = 3.1415926535;
const int marching_max = 128;
const int obj_num = 2;


void main(){

    vec2 st = (gl_FragCoord.xy-resolution.xy)/min(resolution.x,resolution.y);
    float st_x = st.x;
    float st_z = st.y;

    Material mate1;
    Material mate2;

    mate1.albedo = vec3(1.0,0.0,1.0);
    mate2.albedo = vec3(1.0,0.0,0.0);

    mate1.kind = 1;
    mate2.kind = 1;

    Camera cam;
    cam.pos = vec3(0.0,0.0,0.0);
    cam.fov = (PI*30.0)/(2.0*180.0);

    Object objs[SIZE_OF_OBJS_ARRAY];

    objs[0].pos = vec3(2.5,15.0,2.5);
    objs[0].rot = vec3(0.0,0.0,0.0);
    objs[0].kind = 1;
    objs[0].params[0] = 1.0;
    objs[0].material = mate1;

    objs[1].pos = vec3(-2.5,15.0,-2.5);
    objs[1].rot = vec3(0.0,0.0,0.0);
    objs[1].kind = 1;
    objs[1].params[0] = 1.0;
    objs[1].material = mate2;



    vec3 ray;
    ray = calc_ray(cam,st_x,st_z);

    bool hitflag = false;
    int hitnum = -1;
    raymarching(cam.pos,ray,objs,hitflag,hitnum);



    if (hitflag){
        vec3 col;
        vec3 origin = cam.pos;
        
        for (int i = 0; i < obj_num; ++ i){//なんとobjs[hitnum]は通らない。可読性のためにcolorを分離したいからこうなった。
            if (i == hitnum){
                col = material_color(objs[i].material);
            }
        };
        gl_FragColor = vec4(col,1.0);
    }else{
        gl_FragColor = vec4(0.0, 0.0 ,0.0, 1.0);
    }


}

void raymarching(in vec3 origin,in vec3 ray,in Object[SIZE_OF_OBJS_ARRAY] objs,inout bool hitflag,inout int hitnum){
    float max_dis = 100.0;
    float min_dis = 0.001;
    float rlen = 0.0;
    const int max_loop = 100;
    for (int i = 0; i < max_loop; ++i){
        float shortest = 1e9;    
        for (int j = 0; j < obj_num; ++ j){
            float d = distance_func(objs[j],origin+ray*rlen);
            if (d < shortest){
                shortest = d;
                hitnum = j;
            }
        }
        if (rlen > max_dis){
            break;
        }
        if (abs(shortest) < min_dis){
            hitflag = true;
            break;
        }else{
            rlen += shortest;
        }
    }
}

vec3 calc_ray(in Camera cam,in float x,in float z){
    return normalize(vec3(sin(cam.fov)*x ,cos(cam.fov),sin(cam.fov)*z));
}

float distance_func(in Object obj,in vec3 rayhead){
    vec3 p = rayhead-obj.pos;
    if (obj.kind == 1){
        return length(p)-obj.params[0];
    }
}

vec3 calc_norm(in Object obj,in vec3 hitpos){
    float eps = 0.001;
    return normalize(
        vec3(
            distance_func(obj,hitpos+vec3(eps,0.0,0.0))-distance_func(obj,hitpos+vec3(-eps,0.0,0.0)),
            distance_func(obj,hitpos+vec3(0.0,eps,0.0))-distance_func(obj,hitpos+vec3(0.0,-eps,0.0)),
            distance_func(obj,hitpos+vec3(0.0,0.0,eps))-distance_func(obj,hitpos+vec3(0.0,0.0,-eps))
        )
    );
}

vec3 material_color(in Material mat){
    if(mat.kind==1){
        return mat.albedo/PI;
    }
}
vec3 calc_light(in Light light,in vec3 hitpos){
    if(light.kind==1){
        vec3 ray_direc = vec3(1.0,1.0,-1.0);
        ray_direc = normalize(ray_direc);
        return vec3(1.0,1.0,1.0);
    }
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
