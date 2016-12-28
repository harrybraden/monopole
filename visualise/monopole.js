var UNIT = 100

var scene = new THREE.Scene()
  , camera = new THREE.PerspectiveCamera(
        75, window.innerWidth / window.innerHeight, 1, 10000 )
  , renderer = new THREE.WebGLRenderer();

//camera.position.x = 100;
//camera.position.y = 100;
camera.position.z = 200;
camera.lookAt(0,0,0)


var drawDebug = function(){
  scene.background = new THREE.Color( 0xbbbbbb );
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

var drawImagePlane = function(img, z, k){
  var array = new Uint8Array(img.data.buffer) // chrome bug
  var texture = new THREE.DataTexture(
                      array,
                      img.width,
                      img.height,
                      THREE.RGBAFormat//THREE.UnsignedByteType
                    )
  var material = new THREE.MeshBasicMaterial( { 
    //map: texture,
    color: 'blue',
    alphaMap: texture
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
  planeTL.translateZ(-z)

  planeBL.rotateY(Math.PI)
  planeBL.rotateX(Math.PI)
  planeBL.translateX(UNIT/2)
  planeBL.translateY(UNIT/2)
  planeBL.translateZ(z)
  
  planeBR.rotateX(Math.PI)
  planeBR.translateX(UNIT/2)
  planeBR.translateY(UNIT/2)
  planeBR.translateZ(-z)

  planeTR.translateX(UNIT/2)
  planeTR.translateY(UNIT/2)
  planeTR.translateZ(z)

  scene.add(planeTL);
  scene.add(planeTR);
  scene.add(planeBL);
  scene.add(planeBR);
}

function animate() {
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
  var data = ctx.getImageData(k * 60, z*60 , 60, 60);

  var res = new Uint8Array(60*60)
  for (var i = 0; i < 60*60; i++){
    res[i] = data[i*4]
  }
  return res

  //var array = new Uint8Array(img.data.buffer) // chrome bug

}

function init() {
    var ambientLight = new THREE.AmbientLight(0x555555);
    scene.add(ambientLight);

    renderer.setSize( window.innerWidth, window.innerHeight );
    document.body.appendChild( renderer.domElement );
    animate()
}

loadData(function(ctx){
  init()
  var data = ctx.getImageData(0, 0, 60, 60);
  drawDebug()
  drawImagePlane(data, -10, 0)
})
