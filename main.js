window.addEventListener('DOMContentLoaded', init);
let mat,mesh;

function geo(){
 let g = new THREE.BufferGeometry();
 let vert = new Float32Array([
     1.0,1.0,0.0,
     1.0,-1.0,0.0,
     -1.0,-1.0,0.0,
     1.0,1.0,0.0,
     -1.0,1.0,0.0,
     -1.0,-1.0,0.0
 ])
 g.addAttribute('position',new THREE.BufferAttribute(vert,3))
 return g;
}

let shade_par={
    uniforms:{
        "time":{value:1.0},
        "resolution": {value: new THREE.Vector2(pzm.con.init_width, pzm.con.init_height)},
        "mouse" : {value: new THREE.Vector2(0.0,0.0)},
    },
    vertexShader:[
        "precision mediump float;",
        "void main(){",
        "gl_Position = vec4( position, 1.0 );",
        "}"

    ].join("\n"),
    fragmentShader:"",
    side: THREE.DoubleSide,
    transparent: true
}

let gp_op = {
  stop:false
}

let stats = new Stats();
stats.showPanel(0); // 0: fps, 1: ms, 2: mb, 3+: custom
let statson = false;

const mouse = new THREE.Vector2(0.0,0.0);
function init(){
  if (statson){
    document.body.appendChild( stats.dom );
  }
  
  init_main_scene(pzm.par,pzm.con)
  const canvas = document.querySelector('#main_canvas');
  canvas.addEventListener('mousemove', handleMouseMove);
  function handleMouseMove(event) {
    const element = event.currentTarget;
    const x = event.clientX - element.offsetLeft;
    const y = event.clientY - element.offsetTop;
    const w = element.offsetWidth;
    const h = element.offsetHeight;
    mouse.x = ( x / w ) * 2 - 1;
    mouse.y = -( y / h ) * 2 + 1;
  }
    var xmlhttp = new XMLHttpRequest();
    xmlhttp.onreadystatechange = function() {
        if (xmlhttp.readyState == 4) {
          if (xmlhttp.status == 200) { 
          } else {
            alert("status = " + xmlhttp.status);
          } 
        }
      }
      xmlhttp.open("GET", "./frag.glsl");
      xmlhttp.send();
      xmlhttp.onload=(function() {
        data = xmlhttp.responseText
        shade_par.fragmentShader = data;
        let geometry = geo()
        mat = new THREE.ShaderMaterial(shade_par)
        mesh = new THREE.Mesh(geometry,mat);
        pzm.par.scene.add(mesh);
        pzm.con.sttime = performance.now();
      tick();
      //rend();
      })


}

function rend(){
  time = performance.now();
    mat.uniforms.time.value = time;
    mat.uniforms.mouse.value = mouse;
    pzm.par.renderer.render(pzm.par.scene, pzm.par.camera);
}

function tick() {
  stats.begin();
  if (!gp_op.stop){
    time = performance.now();
    mat.uniforms.time.value = (time-pzm.con.sttime)/1000;
    mat.uniforms.mouse.value = mouse;
    pzm.par.renderer.render(pzm.par.scene, pzm.par.camera);
    pzm.par.control_cam.update();
    pzm.par.control_light.update();
  }
  stats.end();
  pzm.par.frames += 1;
  if (pzm.par.frames < 180){
    requestAnimationFrame(tick);
  }
  }




function init_main_scene(main,con){
    main.renderer = new THREE.WebGLRenderer({
        canvas: document.querySelector('#main_canvas'),
        preserveDrawingBuffer: true
      });
      main.renderer.setPixelRatio(window.devicePixelRatio);
      main.renderer.setSize(con.init_width, con.init_height);
      main.scene = new THREE.Scene();
      main.camera = new THREE.PerspectiveCamera(45, con.init_width / con.init_height, 1, 50000);
      main.light = lighting(main.scene);
      main.control_cam = new THREE.OrbitControls(main.camera, main.renderer.domElement);
      main.control_light = new THREE.OrbitControls(main.light, main.renderer.domElement);
      main.control_light.enablePan = false;
      main.control_light.enableZoom = false;
      main.control_cam.enableKeys = false;
      main.control_cam.enablePan = false;
      main.control_light.enableKeys = false;
      main.camera.position.set(0, 0, 25);
      main.renderer.setClearAlpha(0);
}

function lighting(scene) {
    var light = new THREE.DirectionalLight(0xFFFFFF);
    light.intensity = 2;
    light.position.set(0, 0, 100);
    scene.add(light);
    return light;
  }