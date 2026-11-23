// Display real startup/runtime errors instead of leaving a blank visualization.
(() => {
 const panel=document.getElementById('startup-status');
 let failed=false,ready=false;
 function show(message){failed=true;panel.hidden=false;panel.setAttribute('role','alert');panel.style.color='#f5b6a8';panel.textContent='起動・実行エラー: '+message;}
 globalThis.addEventListener('error',event=>{
   const error=event.error;
   show(error?error.name+': '+error.message:event.message||'必要なコードを読み込めませんでした。');
 });
 globalThis.addEventListener('unhandledrejection',event=>show(String(event.reason&&event.reason.message||event.reason)));
 globalThis.twoLinkBootComplete=()=>{if(ready||failed)return;ready=true;panel.hidden=true;};
})();
