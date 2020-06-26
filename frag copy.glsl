precision mediump float;
uniform float time;
uniform vec2 resolution;
uniform vec2 mouse;

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


const float pi = 3.1415926535;
vec2 oo = vec2(0.0,0.0);
const int marching_max = 128;
const float eps = 0.001;

vec3 sp_all[8];


// vec3 sp_all[2] = vec3[](vec3(-0.2,0.4,0.1),vec3(0.15,0.3,0.5));

void main(){

sp_all[0] = rot_by_Rodrigues(vec3(-0.2,0.4+0.2*sin(time/2000.0+0.1),0.1),vec3(0.0,1.0,0.0),time/100.0);
sp_all[1] = rot_by_Rodrigues(vec3(0.15,0.3+0.3*sin(time/500.0),0.5),vec3(0.0,1.0,0.0),time/50.0);
sp_all[2] = rot_by_Rodrigues(vec3(-0.5,-0.1+0.3*sin(time/4000.0+1.0),-0.4),vec3(0.0,1.0,0.0),time/200.0);
sp_all[3] = rot_by_Rodrigues(vec3(0.4,0.7+0.2*sin(time/3000.0-0.1),0.2),vec3(0.0,1.0,0.0),time/300.0);
sp_all[4] = rot_by_Rodrigues(vec3(-0.7,0.8+0.1*sin(time/1000.0+0.7),-0.7),vec3(0.0,1.0,0.0),time/150.0);
sp_all[5] = rot_by_Rodrigues(vec3(0.3,-0.6+0.3*sin(time/1500.0-0.3),0.6),vec3(0.0,1.0,0.0),time/500.0);
sp_all[6] = rot_by_Rodrigues(vec3(0.9,-0.1+0.2*sin(time/2000.0-0.7),-0.3),vec3(0.0,1.0,0.0),time/120.0);
sp_all[7] = rot_by_Rodrigues(vec3(-0.8,0.17+0.3*sin(time/800.0+0.6),-0.8),vec3(0.0,1.0,0.0),time/250.0);

vec2 st = (gl_FragCoord.xy-resolution.xy)/min(resolution.x,resolution.y);
float angle = 60.0;
float fov = angle * 0.5 * pi /180.0;
vec3 cpos = vec3(0.0,0.2,-2.5); //postion of camera
vec3 ray = normalize(vec3(sin(fov)*st.x , sin(fov)*st.y , cos(fov)));

// transform octa_trans;
// octa_trans.pos = vec3(sin(time/300.0),sin(2.0+time/300.0),1.0);
// octa_trans.rot = vec3(time/20.0,0.0,time/20.0);

float dist = 10000000.0; // length to next object
float rlen = 0.0; //additional length to ray pos
vec3 rpos = cpos; //position of the tip of ray


for (int i = 0; i < marching_max; ++i){
    for (int j = 0; j < 8; ++ j){
        float d = df_sphere(rpos-sp_all[j],0.2);
        dist = min(dist,d);
    }
    // vec3 temp_pos = rpos-octa_trans.pos;
    // temp_pos = rot_by_Rodrigues(temp_pos,vec3(1.0,0.0,0.0),octa_trans.rot.x);
    // temp_pos = rot_by_Rodrigues(temp_pos,vec3(0.0,1.0,0.0),octa_trans.rot.y);
    // temp_pos = rot_by_Rodrigues(temp_pos,vec3(0.0,0.0,1.0),octa_trans.rot.z);
    // dist = df_octa(temp_pos,1.0);
    rlen += dist;
    rpos = cpos + ray*rlen;
}
// vec3 temp_pos = rpos-octa_trans.pos;
// temp_pos = rot_by_Rodrigues(temp_pos,vec3(1.0,0.0,0.0),octa_trans.rot.x);
// temp_pos = rot_by_Rodrigues(temp_pos,vec3(0.0,1.0,0.0),octa_trans.rot.y);
// temp_pos = rot_by_Rodrigues(temp_pos,vec3(0.0,0.0,1.0),octa_trans.rot.z);
// float is_hit = df_octa(temp_pos,1.0);

bool flag = true;
for (int i = 0; i < 8; ++ i){
        float d = df_sphere(rpos-sp_all[i],0.2);
        if (d<eps){
            flag = false;
            break;
        }
}

if (flag){     
    gl_FragColor = vec4(vec3(0.7), 1.0);
}else{
    for (int i = 0; i < 8; ++ i){
        float d = df_sphere(rpos-sp_all[i],0.2);
        if (d<eps){
            vec3 temp = vec3(abs(sp_all[i].x),abs(sp_all[i].y),abs(sp_all[i].z));
            gl_FragColor = vec4((normalize(sp_all[i])+0.2)*sin(time/800.0), 1.0);
            break;
        }
    }
}


// if(abs(is_hit) < eps){
//     /*******geometric************/
//     //rpos = rot_by_Rodrigues(rpos,vec3(0.0,1.0,0.0),10.0);
//     //vec3 norm = torus_norm(rpos-sp_pos,vec2(1.0,0.5));
//     vec3 norm = octa_norm(temp_pos,1.0);
//     vec3 lightDir = normalize(vec3(1.0,1.0,-1.0));
//     lightDir = rot_by_Rodrigues(lightDir,vec3(1.0,0.0,0.0),octa_trans.rot.x);
//     lightDir = rot_by_Rodrigues(lightDir,vec3(0.0,1.0,0.0),octa_trans.rot.y);
//     lightDir = rot_by_Rodrigues(lightDir,vec3(0.0,0.0,1.0),octa_trans.rot.z);
//     float lightpow = 1.0*dot(lightDir,norm)*pi;
//     vec3 envlight = vec3(0.2,0.2,0.2);
//     vec3 viewDir = normalize(cpos-rpos);
//     vec3 h = normalize(lightDir+viewDir);
//     float dotNL = saturate(dot(norm,lightDir));
//     float dotNV = saturate(dot(norm,viewDir));
//     float dotNH = saturate(dot(norm,h));
//     float dotVH = saturate(dot(viewDir,h));
//     float dotLV = saturate(dot(lightDir,viewDir));
//     /****************************/

//     /***********materila**********/
//     float metallic = 0.5;
//     float roughness = 0.5;
//     vec3 rgb_in = vec3(255,255,255);
//     rgb_in = rgb_in/255.0;
//     vec3 albedo = rgb_in;
//     float a = roughness * roughness;
//     vec3 diffcolor = mix(albedo,vec3(0.0),metallic);
//     vec3 speccolor = mix(vec3(0.04),albedo,metallic);
//     /*****************************/

//     vec3 diffuseBRDF = albedo/pi;
//     vec3 test5 = normed_phong(speccolor,30.0,viewDir,norm,lightDir);
//     vec3 easy = metallic*lightpow*test5 + (1.0-metallic)*lightpow*diffuseBRDF+envlight*albedo;
//     gl_FragColor = vec4(vec3(easy), 1.0);
// }else{
//     gl_FragColor = vec4(vec3(0.0), 1.0);
// }
}



float atan2(float y, float x){
    return x == 0.0 ? sign(y)*pi/2.0 : atan(y, x);
}

float df_sphere(vec3 pos,float size){
    return length(pos) - size;
}

vec3 sp_norm(vec3 pos,float size){
    return normalize(
        vec3(
            df_sphere(pos+vec3(eps,0.0,0.0),size)-df_sphere(pos+vec3(-eps,0.0,0.0),size),
            df_sphere(pos+vec3(0.0,eps,0.0),size)-df_sphere(pos+vec3(0.0,-eps,0.0),size),
            df_sphere(pos+vec3(0.0,0.0,eps),size)-df_sphere(pos+vec3(0.0,0.0,-eps),size)
        )
    );
}

float df_box(vec3 pos, vec3 size){
    vec3 d = abs(pos)-size;
    return length(max(d,0.0))+ min(max(d.x,max(d.y,d.z)),0.0);
}

vec3 box_norm(vec3 pos,vec3 size){
        return normalize(
        vec3(
            df_box(pos+vec3(eps,0.0,0.0),size)-df_box(pos+vec3(-eps,0.0,0.0),size),
            df_box(pos+vec3(0.0,eps,0.0),size)-df_box(pos+vec3(0.0,-eps,0.0),size),
            df_box(pos+vec3(0.0,0.0,eps),size)-df_box(pos+vec3(0.0,0.0,-eps),size)
        )
    );
}

float df_torus (vec3 p, vec2 t){
    vec2 q = vec2(length(p.xz)-t.x,p.y);
    return length(q)-t.y;
}

vec3 torus_norm (vec3 p, vec2 t){
        return normalize(
        vec3(
            df_torus(p+vec3(eps,0.0,0.0),t)-df_torus(p+vec3(-eps,0.0,0.0),t),
            df_torus(p+vec3(0.0,eps,0.0),t)-df_torus(p+vec3(0.0,-eps,0.0),t),
            df_torus(p+vec3(0.0,0.0,eps),t)-df_torus(p+vec3(0.0,0.0,-eps),t)
        )
    );
}

float df_prism(vec3 p,vec2 h){
    vec3 q = abs(p);
    return max(q.z-h.y,max(q.x*0.866025+p.y*0.5,-p.y)-h.x*0.5);
}

vec3 prism_norm(vec3 p,vec2 h){
         return normalize(
        vec3(
            df_prism(p+vec3(eps,0.0,0.0),h)-df_prism(p+vec3(-eps,0.0,0.0),h),
            df_prism(p+vec3(0.0,eps,0.0),h)-df_prism(p+vec3(0.0,-eps,0.0),h),
            df_prism(p+vec3(0.0,0.0,eps),h)-df_prism(p+vec3(0.0,0.0,-eps),h)
        )
    );
}

float df_octa(vec3 p,float s){
    p = abs(p);
    float m = p.x+p.y+p.z-s;
    vec3 q;
         if( 3.0*p.x < m ) q = p.xyz;
    else if( 3.0*p.y < m ) q = p.yzx;
    else if( 3.0*p.z < m ) q = p.zxy;
    else return m*0.57735027;
    
    float k = clamp(0.5*(q.z-q.y+s),0.0,s); 
    return length(vec3(q.x,q.y-s+k,q.z-k)); 
}

vec3 octa_norm(vec3 p, float t){
        return normalize(
        vec3(
            df_octa(p+vec3(eps,0.0,0.0),t)-df_octa(p+vec3(-eps,0.0,0.0),t),
            df_octa(p+vec3(0.0,eps,0.0),t)-df_octa(p+vec3(0.0,-eps,0.0),t),
            df_octa(p+vec3(0.0,0.0,eps),t)-df_octa(p+vec3(0.0,0.0,-eps),t)
        )
    );
}

vec3 microBRDF(float a,float dotNH, float dotNV, float dotNL, float dotVH,vec3 speccolor){
    float D = D_GGX(a,dotNH);
    float G = G_Smith_Schlick_GGX(a,dotNV,dotNL);
    vec3 F = F_Schlick(speccolor,dotVH);

    return ((F*D*G)/(4.0*dotNV*dotNV+0.0001));
    
}


float D_GGX(float a, float dotNH) {
  float a2 = a*a;
  float dotNH2 = dotNH*dotNH;
  float d = dotNH2 * (a2 - 1.0) + 1.0;
  return a2 / (pi * d * d);
}

float G_Smith_Schlick_GGX(float a, float dotNV, float dotNL) {
  float k = a*a*0.5 + 0.0001;
  float gl = dotNL / (dotNL * (1.0 - k) + k);
  float gv = dotNV / (dotNV * (1.0 - k) + k);
  return gl*gv;
}

vec3 F_Schlick(vec3 specularColor, float dotVH) {
  return specularColor + ((1.0 - specularColor) * pow((1.0 - dotVH), 5.0));
}

vec3 rept (vec3 pos, float freq){
    return mod(pos,freq)-freq/2.0;
}

vec3 rgbnormalize(vec3 rgb){
    return vec3(rgb/255.0);
}

vec3 rot_by_Rodrigues(vec3 p, vec3 axis, float angle){
    vec3 a = normalize(axis);
    float theta = angle*pi/180.0;
    float c = cos(theta);
    float s = sin(theta);
    float r = 1.0-c;

    mat3 m = mat3(
        a.x * a.x * r + c,        a.y * a.x * r + a.z * s,  a.z * a.x * r - a.y * s,
        a.x * a.y * r - a.z * s,  a.y * a.y * r + c,        a.z * a.y * r + a.x * s,
        a.x * a.z * r + a.y * s,  a.y * a.z * r - a.x * s,  a.z * a.z * r + c
    );
    return m*p;
}

vec3 normed_phong(vec3 speccolor,float power,vec3 view, vec3 norm, vec3 lightDir){
     vec3 R = -view + ( 2.0 * dot( norm, view ) * norm );
     return speccolor * pow( max (dot(lightDir,R), 0.0 ), power ) * ( ( power + 1.0 )/ ( 2.0 * pi ) );
}
