/* Photo-inspired commuter carriage, built in metres; no travel animation.
 * All meshes remain under #world, sharing the recorded AP stimulus transform.
 */
(() => {
 let train;
 function build(){
  const T=AFRAME.THREE;
  train=document.createElement('a-entity');train.id='train-scene';
  document.querySelector('#world').appendChild(train);
  const group=new T.Group();train.setObject3D('carriage',group);
  // Deterministic fine-grain surfaces: cloth, brushed metal and vinyl.
  function surface(kind,color){
   const c=document.createElement('canvas');c.width=c.height=256;const ctx=c.getContext('2d');
   ctx.fillStyle=color;ctx.fillRect(0,0,256,256);
   let seed=437;function rand(){seed=(1664525*seed+1013904223)>>>0;return seed/4294967296;}
   for(let i=0;i<26000;i++){
    const v=rand()>.5?255:0;ctx.fillStyle=`rgba(${v},${v},${v},${kind==='cloth'?.14:.055})`;
    const x=rand()*256,y=rand()*256;ctx.fillRect(x,y,kind==='metal'?12:1,1);
   }
   if(kind==='cloth'){ctx.fillStyle='rgba(150,177,240,.13)';for(let y=0;y<256;y+=3)ctx.fillRect(0,y,256,1);}
   const map=new T.CanvasTexture(c);map.wrapS=map.wrapT=T.RepeatWrapping;
   map.repeat.set(kind==='cloth'?3:2,2);map.colorSpace=T.SRGBColorSpace;
   return new T.MeshStandardMaterial({map,roughness:kind==='metal'?.35:.88,metalness:kind==='metal'?.65:0});
  }
  const wall=new T.MeshStandardMaterial({color:'#e2e3dc',roughness:.6});
  const steel=surface('metal','#aeb4b5'),floor=surface('vinyl','#969c9c'),cloth=surface('cloth','#183b99');
  const rubber=new T.MeshStandardMaterial({color:'#19211e',roughness:.85});
  const trim=new T.MeshStandardMaterial({color:'#d0d4d1',metalness:.65,roughness:.25});
  const dark=new T.MeshStandardMaterial({color:'#444b4b',roughness:.7});
  const basic=color=>new T.MeshBasicMaterial({color});
  function mesh(geometry,material,x,y,z){const m=new T.Mesh(geometry,material);m.position.set(x,y,z);group.add(m);return m;}
  function box(x,y,z,w,h,d,mat){return mesh(new T.BoxGeometry(w,h,d),mat,x,y,z);}
  function round(w,h,r){const s=new T.Shape(),x=-w/2,y=-h/2;
   s.moveTo(x+r,y);s.lineTo(x+w-r,y);s.quadraticCurveTo(x+w,y,x+w,y+r);
   s.lineTo(x+w,y+h-r);s.quadraticCurveTo(x+w,y+h,x+w-r,y+h);
   s.lineTo(x+r,y+h);s.quadraticCurveTo(x,y+h,x,y+h-r);s.lineTo(x,y+r);s.quadraticCurveTo(x,y,x+r,y);return s;
  }
  function panel(x,y,z,w,h,d,r,mat){return mesh(new T.ExtrudeGeometry(round(w,h,r),{depth:d,bevelEnabled:true,bevelThickness:.012,bevelSize:.012,bevelSegments:3,steps:1,curveSegments:12}),mat,x,y,z);}
  function tube(a,b,r=.018,mat=trim){const av=new T.Vector3(...a),bv=new T.Vector3(...b),v=bv.clone().sub(av);const m=mesh(new T.CylinderGeometry(r,r,v.length(),12),mat,...av.clone().add(bv).multiplyScalar(.5).toArray());m.quaternion.setFromUnitVectors(new T.Vector3(0,1,0),v.normalize());return m;}
  // Standing in the aisle, looking straight at the opposite bench and windows.
  box(0,-.07,0,10,.14,3.3,floor);box(0,2.55,0,10,.12,3.3,wall);
  box(0,.42,-1.62,10,.84,.13,wall);box(0,2.29,-1.62,10,.52,.13,wall);
  box(0,1.28,1.62,10,2.56,.12,wall);
  box(-5,1.28,0,.12,2.56,3.3,wall);box(5,1.28,0,.12,2.56,3.3,wall);
  // Windows alternate with broad ivory pillars and stainless sliding doors.
  for(const x of [-3.9,-1.3,1.3,3.9]){
   box(x-1.25,1.45,-1.62,.28,1.18,.15,wall);
   panel(x,1.47,-1.57,2.2,1.2,.045,.12,steel);
   panel(x,1.47,-1.51,2.1,1.1,.025,.105,rubber);
   // Window is a textured plane behind the frame, not an opaque wall.
   const c=document.createElement('canvas');c.width=1024;c.height=512;const ctx=c.getContext('2d');
   const sky=ctx.createLinearGradient(0,0,0,512);sky.addColorStop(0,'#e4e9de');sky.addColorStop(.42,'#ccd6bb');sky.addColorStop(1,'#506d35');ctx.fillStyle=sky;ctx.fillRect(0,0,1024,512);
   let seed=97;const rand=()=>{seed=(seed*1664525+1013904223)>>>0;return seed/4294967296;};
   // Fixed leafy landscape. Horizontal strokes evoke the reference photograph,
   // but no texture offset is animated.
   for(let i=0;i<4000;i++){const y=160+rand()*352;ctx.fillStyle=`rgba(${40+rand()*60},${65+rand()*65},${24+rand()*42},.24)`;ctx.fillRect(rand()*1024,y,8+rand()*65,1+rand()*5);}
   for(let i=0;i<12;i++){ctx.fillStyle='rgba(218,225,208,.20)';ctx.fillRect(rand()*1024,0,2,512);}
   const tex=new T.CanvasTexture(c);tex.colorSpace=T.SRGBColorSpace;
   panel(x,1.47,-1.465,1.98,.99,.005,.09,new T.MeshBasicMaterial({map:tex}));
   box(x,1.76,-1.435,2.04,.045,.04,steel);
   box(x,1.79,-1.41,2.02,.012,.02,trim);
   for(const dx of [-.72,.72])box(x+dx,1.995,-1.405,.055,.045,.022,trim);
   // Six individually rounded blue cushions, seams and shaded plinth.
   box(x,.24,-1.28,2.32,.43,.62,steel);
   box(x,.21,-.953,2.15,.24,.018,dark);
   for(let j=0;j<23;j++)box(x-1.04+j*.094,.21,-.94,.009,.22,.012,steel);
   for(let j=0;j<4;j++){
    const sx=x-.87+j*.58;
    const seat=panel(sx,.47,-1.3,.555,.59,.09,.065,cloth);seat.rotation.x=-Math.PI/2;
    const back=panel(sx,.82,-1.48,.555,.64,.12,.075,cloth);back.rotation.x=-.08;
    // Cloth welt along the lower edge.
    tube([sx-.23,.477,-.99],[sx+.23,.477,-.99],.006,cloth);
   }
   for(const dx of [-1.2,1.2]){
    panel(x+dx,.67,-1.3,.065,.76,.42,.025,wall);
    tube([x+dx,.95,-1.22],[x+dx,2.35,-1.22],.017);
   }
   // Overhead luggage rack and hanging straps.
   for(const z of [-1.51,-1.03])tube([x-1.18,2.24,z],[x+1.18,2.24,z],.022);
   for(let j=0;j<13;j++)tube([x-1.14+j*.19,2.24,-1.51],[x-1.14+j*.19,2.24,-1.03],.009);
   tube([x-1.2,2.4,-.81],[x+1.2,2.4,-.81],.018);
   for(const dx of [-.8,0,.8]){
    box(x+dx,2.245,-.81,.026,.28,.018,dark);
    const ring=mesh(new T.TorusGeometry(.065,.012,8,24),wall,x+dx,2.06,-.81);
   }
  }
  // Door at the left edge of the main view, like the supplied photo.
  const doorX=-2.6;
  box(doorX,1.23,-1.40,1.02,2.46,.12,steel);
  for(const dx of [-.255,.255]){
   panel(doorX+dx,1.62,-1.315,.38,.67,.012,.07,rubber);
   panel(doorX+dx,1.62,-1.29,.33,.62,.004,.06,basic('#bccfc3'));
  }
  box(doorX,1.23,-1.315,.014,2.43,.015,dark);
  box(doorX,.035,-1.24,1.06,.025,.3,trim);
  for(const dx of [-.55,.55])tube([doorX+dx,.08,-1.27],[doorX+dx,2.42,-1.27],.022);
  // Blank poster holder, ceiling lights, and skirting finish the carriage.
  panel(-2.03,1.85,-1.34,.38,.51,.012,.012,trim);
  panel(-2.03,1.85,-1.315,.35,.48,.004,.005,basic('#e9e8df'));
  for(const z of [-.65,.65])box(0,2.475,z,9,.035,.13,basic('#fff7df'));
  box(0,.035,-1.5,10,.07,.06,steel);
  // Self-contained lighting also works in the existing preview/Quest scene.
  group.add(new T.HemisphereLight('#ffffff','#737679',1.35));
  const light=new T.DirectionalLight('#fff6e8',1.25);light.position.set(-2,4,3);group.add(light);
 }
 window.TrainScene={select(environment){
  if(!train&&environment==='train')build();
  if(train)train.setAttribute('visible',environment==='train');
  document.querySelector('#sweep-room').setAttribute('visible',environment!=='train');
 }};
})();
