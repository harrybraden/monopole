var TILE_SIZE = 60
var UNIT = TILE_SIZE
var KVAL = 0.3
var INTENSITY = 0.6
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
  scene.background = new THREE.Color( 0xEEEEEE );
  var geometry = new THREE.BoxGeometry( 3, 3, 3 );
  var origin = new THREE.Mesh( 
                      geometry,
                      new THREE.MeshBasicMaterial( { color: 0xff0000 } )
                    )
  scene.add(origin);

  var gridHelperY = new THREE.GridHelper( UNIT, 6, 0x444444, 0x00ff00);
	scene.add( gridHelperY );

  //var gridHelperX = new THREE.GridHelper( UNIT, 10, 0x444444, 0x0000ff);
	//gridHelperX.geometry.rotateX(Math.PI/2)
	//scene.add( gridHelperX );
}

var drawImagePlane = function(buf, z){
  var texture = new THREE.DataTexture(
                      buf,
                      60,
                      60,
                      THREE.RGBAFormat
                    )
  var material = new THREE.MeshBasicMaterial( { 
    map: texture,
    color: 0xff9999,
    transparent: true,
    alphaMap: texture,
    alphaTest: INTENSITY / 1.5
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
  for (var z = 0; z < TILE_SIZE; z++) {
    var data = getDataFor(ctx, z, KVAL)
    drawImagePlane(data, z)
    drawImagePlane(data, -z)
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
  // K is in range 0.01 - 0.99
  var koffset = Math.round((k - 0.01) * 100)
  var imdata = ctx.getImageData(koffset * TILE_SIZE, z* TILE_SIZE , TILE_SIZE, TILE_SIZE);
  var data = imdata.data

  var res = new Uint8Array(TILE_SIZE*TILE_SIZE*4)
  for (var i = 0; i < TILE_SIZE*TILE_SIZE; i++){
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
      KVAL = e.target.value
			document.getElementById('kVal').textContent = KVAL
    }
    document.getElementById('intensity').onchange = function(e){
      INTENSITY=e.target.value
			document.getElementById('intensityVal').textContent = INTENSITY
    }
}

loadData(function(ctx){
  DATA = ctx
  init()
  animate()
})
