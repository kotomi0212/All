(function(g){'use strict';
const P={g:9.81,m1:28,l1:.90,c1:.45,I1:2.1,m2:42,l2:.65,c2:.31,I2:1.45,
  kPassiveAnkle:320,kPassiveHip:140,dPassiveAnkle:35,dPassiveHip:22,
  kSelfAnkle:660,kSelfHip:280,dSelfAnkle:130,dSelfHip:70,
  visualTau:.18,visualAnkle:28,visualHip:0,crossKAH:110,crossKHA:70,crossDAH:16,crossDHA:10};
const C=(re,im=0)=>({re,im}),add=(a,b)=>C(a.re+b.re,a.im+b.im),sub=(a,b)=>C(a.re-b.re,a.im-b.im),mul=(a,b)=>C(a.re*b.re-a.im*b.im,a.re*b.im+a.im*b.re),div=(a,b)=>{const d=b.re*b.re+b.im*b.im;return C((a.re*b.re+a.im*b.im)/d,(a.im*b.re-a.re*b.im)/d)},scale=(a,s)=>C(a.re*s,a.im*s),abs=a=>Math.hypot(a.re,a.im),phase=a=>Math.atan2(a.im,a.re);
function mm(A,B){return A.map(r=>B[0].map((_,j)=>r.reduce((sum,v,k)=>sum+v*B[k][j],0)))}
/* Joint torque decomposition used throughout the UI:
 * tau = tau_visual - (Kpassive + Kself + Kcross) phi
 *                  - (Dpassive + Dself + Dcross) phi_dot
 * Passive terms approximate intrinsic muscle/tendon/joint mechanics near upright.
 * Self terms are active local stabilization. Cross terms are active and Model-B only.
 */
function matrices(opt={}){const p={...P,...opt.params},passive=opt.passiveScale??1,self=opt.selfScale??1,cross=opt.control?(opt.crossScale??1):0;const m12=opt.inertia===false?0:p.m2*p.l1*p.c2;const M=[[p.I1+p.m1*p.c1*p.c1+p.m2*p.l1*p.l1,m12],[m12,p.I2+p.m2*p.c2*p.c2]],G=[[p.g*(p.m1*p.c1+p.m2*p.l1),0],[0,p.g*p.m2*p.c2]],B=[[1,-1],[0,1]],Cq=[[1,0],[-1,1]];const Kpassive=[[p.kPassiveAnkle*passive,0],[0,p.kPassiveHip*passive]],Dpassive=[[p.dPassiveAnkle*passive,0],[0,p.dPassiveHip*passive]],Kself=[[p.kSelfAnkle*self,0],[0,p.kSelfHip*self]],Dself=[[p.dSelfAnkle*self,0],[0,p.dSelfHip*self]],Kcross=[[0,p.crossKAH*cross],[p.crossKHA*cross,0]],Dcross=[[0,p.crossDAH*cross],[p.crossDHA*cross,0]],Kjoint=Kpassive.map((r,i)=>r.map((v,j)=>v+Kself[i][j]+Kcross[i][j])),Djoint=Dpassive.map((r,i)=>r.map((v,j)=>v+Dself[i][j]+Dcross[i][j]));return{p,M,G,K:mm(mm(B,Kjoint),Cq),D:mm(mm(B,Djoint),Cq),Kjoint,Djoint,Kpassive,Dpassive,Kself,Dself,Kcross,Dcross,B,Cq}}
function inv2(A){const det=sub(mul(A[0][0],A[1][1]),mul(A[0][1],A[1][0]));return[[div(A[1][1],det),div(scale(A[0][1],-1),det)],[div(scale(A[1][0],-1),det),div(A[0][0],det)]]}function mv(A,v){return A.map(r=>add(mul(r[0],v[0]),mul(r[1],v[1])))}
function response(f,opt={}){const {p,M,G,K,D,B}=matrices(opt),w=2*Math.PI*f,s=C(0,w),vf=div(C(1),add(C(1),scale(s,p.visualTau))),tau=[scale(vf,p.visualAnkle),scale(vf,p.visualHip)],rhs=[sub(tau[0],tau[1]),tau[1]];const Z=[[C(K[0][0]-G[0][0]-w*w*M[0][0],w*D[0][0]),C(K[0][1]-w*w*M[0][1],w*D[0][1])],[C(K[1][0]-w*w*M[1][0],w*D[1][0]),C(K[1][1]-G[1][1]-w*w*M[1][1],w*D[1][1])]];const q=mv(inv2(Z),rhs);return{f,w,q,hip:sub(q[1],q[0]),visual:vf}}
function accel(state,visual,opt={}){const {p,M,G,K,D}=matrices(opt),q=[state[0],state[1]],v=[state[2],state[3]],tau=[p.visualAnkle*state[4],p.visualHip*state[4]],rhs=[tau[0]-tau[1],tau[1]];for(let i=0;i<2;i++)for(let j=0;j<2;j++)rhs[i]+=G[i][j]*q[j]-K[i][j]*q[j]-D[i][j]*v[j];const det=M[0][0]*M[1][1]-M[0][1]*M[1][0],a=[(rhs[0]*M[1][1]-M[0][1]*rhs[1])/det,(M[0][0]*rhs[1]-rhs[0]*M[1][0])/det];return[v[0],v[1],a[0],a[1],(visual-state[4])/p.visualTau]}
function step(x,dt,visual,opt){const f=(y)=>accel(y,visual,opt),k1=f(x),y2=x.map((v,i)=>v+k1[i]*dt/2),k2=f(y2),y3=x.map((v,i)=>v+k2[i]*dt/2),k3=f(y3),y4=x.map((v,i)=>v+k3[i]*dt),k4=f(y4);return x.map((v,i)=>v+dt*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6)}
g.PostureModel={P,matrices,response,step,complex:{abs,phase}};
})(typeof window==='undefined'?globalThis:window);
