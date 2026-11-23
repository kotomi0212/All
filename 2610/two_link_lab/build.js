// Package the exact same physics and UI into a browser-openable single HTML file.
const fs=require('node:fs'),path=require('node:path');
let html=fs.readFileSync(path.join(__dirname,'index.template.html'),'utf8');
const inline=name=>'<script>\n'+fs.readFileSync(path.join(__dirname,name),'utf8').replace(/<\/script/gi,'<\\/script')+'\n</script>';
html=html.replace('<script src="physics.js"></script>',()=>inline('boot.js')+'\n'+inline('physics.js'));
html=html.replace('<script src="app.js"></script>',()=>inline('app.js'));
if(/<script\b[^>]*\bsrc\s*=/i.test(html))throw new Error('Standalone output still has external scripts');
for(const name of ['index.html','two-link-lab.html'])fs.writeFileSync(path.join(__dirname,name),html,'utf8');
console.log('Built index.html and two-link-lab.html (self-contained, no script downloads).');
