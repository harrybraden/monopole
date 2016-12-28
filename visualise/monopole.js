var UNIT = 60
var KVAL = 30
var INTENSITY = 0.5
var DATA = null

var scene = new THREE.Scene()
  , camera = new THREE.PerspectiveCamera(
        75, window.innerWidth / window.innerHeight, 1, 10000 )
  , renderer = new THREE.WebGLRenderer();

//camera.position.x = 100;
//camera.position.y = 100;
camera.position.z = 200;
camera.lookAt(0,0,0)


var drawDebug = function(){
  scene.background = new THREE.Color( 0x333333 );
  var geometry = new THREE.BoxGeometry( 3, 3, 3 );
  var origin = new THREE.Mesh( 
                      geometry,
                      new THREE.MeshBasicMaterial( { color: 0xff0000 } )
                    )
  scene.add(origin);
  var up = new THREE.Mesh( 
                      geometry,
                      new THREE.MeshBasicMaterial( { color: 0x00ff00 } )
                    )
  up.position.y = 100
  scene.add(up);
}

var drawImagePlane = function(buf, z, k){
  var texture = new THREE.DataTexture(
                      buf,
                      60,
                      60,
                      THREE.RGBAFormat
                    )
  var material = new THREE.MeshBasicMaterial( { 
    map: texture,
    //color: 'blue',
    transparent: true,
    alphaMap: texture,
    alphaTest: INTENSITY
  } )

  texture.needsUpdate = true
  material.side = THREE.DoubleSide

  var planeTR = new THREE.Mesh(
      new THREE.PlaneGeometry(UNIT, UNIT),
      material
  )

  var planeTL = planeTR.clone()
  var planeBL = planeTR.clone()
  var planeBR = planeTR.clone()

  planeTL.rotateY(Math.PI)
  planeTL.translateX(UNIT/2)
  planeTL.translateY(UNIT/2)
  planeTL.translateZ(-z )

  planeBL.rotateY(Math.PI)
  planeBL.rotateX(Math.PI)
  planeBL.translateX(UNIT/2)
  planeBL.translateY(UNIT/2)
  planeBL.translateZ(z )
  
  planeBR.rotateX(Math.PI)
  planeBR.translateX(UNIT/2)
  planeBR.translateY(UNIT/2)
  planeBR.translateZ(-z )

  planeTR.translateX(UNIT/2)
  planeTR.translateY(UNIT/2)
  planeTR.translateZ(z )

  scene.add(planeTL);
  scene.add(planeTR);
  scene.add(planeBL);
  scene.add(planeBR);
}

var _animationState = {}

function animate() {
  var ctx = DATA

  // If same params, skip
  if (_animationState.k == KVAL &&
      _animationState.intensity == INTENSITY){
    renderer.render(scene, camera);
    requestAnimationFrame( animate );
    return
  }
  console.log('.', KVAL, INTENSITY)
  _animationState.k = KVAL
  _animationState.intensity = INTENSITY

  while (scene.children.length) {
      scene.remove(scene.children[0]);
  }

  drawDebug()
  for (var z = 0; z < 60; z++) {
    var data = getDataFor(ctx, z, KVAL)
    drawImagePlane(data, z, KVAL)
    drawImagePlane(data, -z, KVAL)
  }

  renderer.render(scene, camera);
  requestAnimationFrame( animate );
}

var loadData = function(cb){
  var img = new Image()
  img.src = './data.png'
  img.onload = function(){
    var canvas = document.createElement('canvas')
      , context = canvas.getContext('2d')

    canvas.width = img.width;
    canvas.height = img.height;
    context.drawImage(img, 0, 0, img.naturalWidth, img.naturalHeight);
    cb(context)
  }
}

var getDataFor = function(ctx, z, k) {
  var imdata = ctx.getImageData(k * 60, z*60 , 60, 60);
  var data = imdata.data

  var res = new Uint8Array(60*60*4)
  for (var i = 0; i < 60*60; i++){
    res[i*4]     = data[i*4] // r
    res[i*4 + 1] = data[i*4 + 1] // g
    res[i*4 + 2] = data[i*4 + 2] // b
    res[i*4 + 3] = data[i*4 + 3] // a
  }
  return res

  //var array = new Uint8Array(img.data.buffer) // chrome bug

}

function init() {
    var ambientLight = new THREE.AmbientLight(0x555555);
    scene.add(ambientLight);

    controls = new THREE.OrbitControls(camera, renderer.domElement)
    renderer.setSize( window.innerWidth, window.innerHeight );
    document.body.appendChild( renderer.domElement );

    document.getElementById('k').onchange = function(e){
      KVAL=e.target.value
    }
    document.getElementById('intensity').onchange = function(e){
      INTENSITY=e.target.value
    }

}

loadData(function(ctx){
  DATA = ctx
  init()
  animate()
})
