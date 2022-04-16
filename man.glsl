vec4 tpos = vec4(0.0, 0.0, 0.0, 0.0);
float angle_tpos;
bool done = false;
float dist = 0.0;

mat2 Rot(float a){float s = sin(a);float c = cos(a);return mat2(c, -s, s, c);}
mat3 RotX(float a){float s = sin(a);float c = cos(a);return mat3(1., .0, .0, .0, c, -s, .0, s, c);}
mat3 RotY(float a){float s = sin(a);float c = cos(a);return mat3(c, .0, s, .0, 1., .0, -s, .0, c);}
mat3 RotZ(float a){float s = sin(a);float c = cos(a);return mat3(c, -s, .0, s, c, .0, .0, .0, 1.);}


float opExtrusion(in vec3 p, in float d, in float h) {
    // d is the distance to the 2D shape using the x and y components of p
    vec2 w = vec2(d, abs(p.z) - h);
    return min(max(w.x, w.y), 0.0) + length(max(w, 0.0));
}

/*
vec4 opElongate( in vec3 p, in vec3 h )
{
    //return vec4( p-clamp(p,-h,h), 0.0 ); // faster, but produces zero in the interior elongated box
    
    vec3 q = abs(p)-h;
    return vec4( max(q,0.0), min(max(q.x,max(q.y,q.z)),0.0) );
}*/


vec3 opTwist(in vec3 p, in float t) {
    float c = cos(t * p.y);
    float s = sin(t * p.y);
    mat2 m = mat2(c, -s, s, c);
    return vec3(m * p.xz, p.y);
}

// from iq
float Noise( in vec3 x )
{
    vec3 p = floor(x);
    vec3 f = fract(x);
  	f = f*f*(3.0-2.0*f);
  	vec2 uv = (p.xy+vec2(37.0,17.0)*p.z) + f.xy;
  	vec2 rg = textureLod( iChannel0, (uv+0.5)/256.0, 0.0).yx;
  	return mix( rg.x, rg.y, f.z );
}

// ref https://www.shadertoy.com/view/Xs33Df
float Noise3D(in vec3 p){
    const vec3 s = vec3(7, 157, 113);
	vec3 ip = floor(p); // Unique unit cell ID.
    vec4 h = vec4(0., s.yz, s.y + s.z) + dot(ip, s);
	p -= ip; // Cell's fractional component.
    p = p*p*(3. - 2.*p);
    h = mix(fract(sin(h)*43758.5453), fract(sin(h + s.x)*43758.5453), p.x);
    h.xy = mix(h.xz, h.yw, p.y);
    return mix(h.x, h.y, p.z); // Range: [0, 1].
	
}

vec2 hash( vec2 p )
{
	p = vec2( dot(p,vec2(127.1,311.7)),
			 dot(p,vec2(269.5,183.3)) );
	return -1.0 + 2.0*fract(sin(p)*43758.5453123);
}

float noise( in vec2 p )
{
	const float K1 = 0.366025404; // (sqrt(3)-1)/2;
	const float K2 = 0.211324865; // (3-sqrt(3))/6;
	
	vec2 i = floor( p + (p.x+p.y)*K1 );
	
	vec2 a = p - i + (i.x+i.y)*K2;
	vec2 o = (a.x>a.y) ? vec2(1.0,0.0) : vec2(0.0,1.0);
	vec2 b = a - o + K2;
	vec2 c = a - 1.0 + 2.0*K2;
	
	vec3 h = max( 0.5-vec3(dot(a,a), dot(b,b), dot(c,c) ), 0.0 );
	
	vec3 n = h*h*h*h*vec3( dot(a,hash(i+0.0)), dot(b,hash(i+o)), dot(c,hash(i+1.0)));
	
	return dot( n, vec3(70.0) );
}
/*
float sdEllipsoid( vec3 p, vec3 r )
{
  float k0 = length(p/r);
  float k1 = length(p/(r*r));
  return k0*(k0-1.0)/k1;
}*/

float sdCylinder(vec3 p, vec3 a, vec3 b, float r) {
	vec3 ab = b-a;
    vec3 ap = p-a;
    
    float t = dot(ab, ap) / dot(ab, ab);
    //t = clamp(t, 0., 1.);
    
    vec3 c = a + t*ab;
    
    float x = length(p-c)-r;
    float y = (abs(t-.5)-.5)*length(ab);
    float e = length(max(vec2(x, y), 0.));
    float i = min(max(x, y), 0.)-0.1;
    
    return e+i;
}

float sdRoundBox( vec3 p, vec3 b, float r )
{
  vec3 q = abs(p) - b;
  return -(length(max(q,0.0)) + min(max(q.x,max(q.y,q.z)),0.0) - r);
}

float Box(vec3 p, vec3 sc, float r){
    return length(max(abs(p)-sc, 0.))-r;
}

float Box2(vec3 p, vec3 sc, float r, float an){
    an = 80.0;
    vec3 cs = vec3(cos(an*3.1415/180.0), sin(an*3.1415/180.0),cos(an*3.1415/180.0));
    float c2 = length(cs);
    
    float a = length(max( abs(p)-(sc)  , 0.0)) - r;
    float b = length(max(abs(p)-vec3(1.0, 2.0, 1.0), 0.));
    return a;
}

float Capsule(vec3 p, vec3 sc, float r){
    vec3 cs = vec3(cos(90.0*3.1415/180.0), sin(90.0*3.1415/180.0),cos(90.0*3.1415/180.0));  
    float a = length(max( abs(p)- (sc*cs)  , 0.0)) - r;
    return a;
}

float opRep( in vec3 p, in vec3 c, vec3 sc, float r )
{
    vec3 q = mod(p+0.5*c,c)-0.5*c;
    return Box( q, sc, r );
}

float opRepS( in vec3 p, in vec3 c, vec4 pos )
{
    vec3 q = mod(p+0.5*c,c)-0.5*c;
    return length(q) - pos.w;
}

float Plane(vec3 p, vec3 n, float h){

    return dot(p, n) + h;
}

vec3 add_tpos(vec3 tp){
    return vec3(tp.x+tpos.x, tp.y+tpos.y, tp.z+tpos.z);
    
}

vec3 rot_y(vec3 tp, float a){
    tp.xz *= Rot(a);
    return tp;

}

float noi(vec3 p){ 
  vec3 f=floor(p),s=vec3(7,157,113);
  p-=f; vec4 h=vec4(0,s.yz,s.y+s.z)+dot(f,s);;
  p=p*p*(3.-2.*p);
  h=mix(fract(sin(h)*43758.5),fract(sin(h+s.x)*43758.5),p.x);
  h.xy=mix(h.xz,h.yw,p.y);
  return mix(h.x,h.y,p.z);  
}

float fbm(vec2 uv)
{
	float f;
	mat2 m = mat2( 1.6,  1.2, -1.2,  1.6 );
	f  = 0.5000*noise( uv ); uv = m*uv;
	f += 0.2500*noise( uv ); uv = m*uv;
	f += 0.1250*noise( uv ); uv = m*uv;
	f += 0.0625*noise( uv ); uv = m*uv;
	f = 0.5 + 0.5*f;
	return f;
}

vec4 eyel = vec4(-0.7,6.5, 3.75,0.3);
vec4 eyer = vec4(0.7 ,6.5, 3.75, 0.3);
vec2 sourcl = vec2(-0.7, 7.5);
vec2 sourcr = vec2(0.7, 7.5);

vec2 GetDist(vec3 p, vec2 uv){


    vec4 s = eyel;
    vec4 s2 = eyer;
    
    vec4 head = vec4(0.0, 6.0, 7.0, 3.);
    
    
    vec2 d=vec2(0.0);;
    float pd = dot(p, vec3(0.0, 1.0, 0.0));
    
     //#e2b1a8
     
    pd = p.y;
   
    //s.z += iTime;
    //s.y += mod(iTime, 3.0);
    vec3 sp = (p)-s.xyz;
    //sp.xz *= Rot(iTime);
    float sd = length((sp))-s.w;
       
       
    
    vec3 sp2 = (p)-s2.xyz;
    //sp2.y += sin(iTime)*0.1;
    float sd2 = length((sp2))-s2.w;
    
    
    vec3 sph = (p)-head.xyz;
    //sp2.y += sin(iTime)*0.1;
    float sdh = length((sph))-head.w;
    
    vec3 bhead = p - vec3(0.0, 4.5, 5.0);
    float sbhd = sdEllipsoid(bhead, vec3(2.5, 2.0, 3.0));
    
    
    
    vec4 tel = vec4(-0.8, 6.75, 3.35, 0.5);
    vec3 tepl = p - tel.xyz;
    float stepld = length(tepl)-tel.w;
    
    vec4 ter = vec4(0.8, 6.75, 3.35, 0.5);
    vec3 tepr = p - ter.xyz;
    float steprd = length(tepr)-ter.w;
    
    vec3 pn = p - vec3(0.0, 6.3-mod(iTime*0.1 ,0.1), 3.5);
    pn *= RotX(-50.0*3.1415/180.0);
    float snose = sdCone(pn, vec2(cos(77.*3.1415/180.0), sin(77.*3.1415/180.0)), 1.2, 0.2);
    
    vec4 stn1 = vec4(-.15, 5.65,2.3, 0.1);
    sp2 = p-stn1.xyz;
    float sdtn1 = length((sp2))-stn1.w;
    
    stn1 = vec4(.15, 5.65,2.3, 0.1);
    sp2 = p-stn1.xyz;
    float sdtn2 = length((sp2))-stn1.w;
    
    vec3 bouche = p - vec3(0.0, 4.9, 2.0);
    float sbchd = sdEllipsoid(bouche, vec3(0.5, clamp(mod(iTime*0.3, 0.2), 0.1, 0.2), 0.5));
    
    pn = p - vec3(0.0, 5.9, 3.15);
    pn *= RotX(-35.0*3.1415/180.0);
    float sbec = sdCone(pn, vec2(cos(75.*3.1415/180.0), sin(75.*3.1415/180.0)), 1.2, 0.15);
    
    vec3 snp = p - vec3(0.0, 3.0, 3.0);
    vec3 a = vec3(0.0, 3.0, 3.0);
    vec3 b = vec3(0.0, 3.0, 1.0);
    
    vec2 se = sdSegment(p, vec3(0,2.5,3.00),  vec3(0.0,0.,3.) );
    float sneckd = se.x-s.y*0.2;
    //d = smax(d,-d2,0.04);
    
    vec3 cap = p - vec3(0.0, 4.0, 4.5);
    cap.y += sin(cap.y+mod(iTime*0.2, 0.2));
    float scap = sdEllipsoid(cap, vec3(3.5, 4.5, 3.5));
    
    vec3 spbc = p - vec3(0.0, 5.0, 5.5);
    float sdbct = Box(spbc, vec3(5.0,5.0, 5.0), 0.0);
    
     cap = p - vec3(0.0, 4.0, 4.0);
     cap.y += sin(cap.y+mod(iTime*0.2, 0.2));
    float scap2 = sdEllipsoid(cap, vec3(3.2, 4.2, 3.2));
    
    float cut = smax(scap, -sdbct, 0.1);
    cut = smax(cut, -scap2, 0.1);
    
    //
    cap = p - vec3(-0.7, 6.5, 3.75);
    float scapel = sdEllipsoid(cap, vec3(0.4, 0.4, 0.4));
    
    spbc = p - vec3(-0.7, 5.15, 3.75);
    spbc.y -= mod(iTime, 0.5);
    float sdbct2 = Box(spbc, vec3(1.0, 1.0, 1.0), 0.0);
    
    float cute = smax(scapel, -sdbct2, 0.4);
    
    cap = p - vec3(0.7, 6.5, 3.75);
    scapel = sdEllipsoid(cap, vec3(0.4, 0.4, 0.4));
    
     spbc = p - vec3(0.7, 5.15, 3.75);
     spbc.y -= mod(iTime, 0.5);
     sdbct2 = Box(spbc, vec3(1.0, 1.0, 1.0), 0.0);
    
    float cute2 = smax(scapel, -sdbct2, 0.4);
    
    
    float dhead = smin(sdh, sbhd, 2.0);
    dhead = smax(dhead, -stepld, 0.5);
    dhead = smax(dhead, -steprd, 0.5);
    dhead = smin(dhead, snose, 0.1);
    dhead = smax(dhead, -sdtn1, 0.1);
    dhead = smax(dhead, -sdtn2, 0.1); 
    dhead = smax(dhead, -sbchd, 0.1); 
    dhead = smin(dhead, sbec, 0.1);
    dhead = smin(dhead, sneckd, 0.3);
             
    float e = mod(iTime, 5.);
    if(e >= 4.5){
        dhead = smin(dhead, cute, 0.1);
        dhead = smin(dhead, cute2, 0.1);
    
    }
    
     
      d.x = sd;
      d.y = -1.0;
      
       
       if(sd2 < d.x)
       {
           d.x = sd2;
          d.y = -2.0;
       }
       
      if(dhead < d.x){
          d.x = dhead;
          d.y = 3.0;
       }
      /* if(scapel < d.x){
          d.x = scapel;
          d.y = 3.0;
       }
       if(sdbct2 < d.x){
          d.x = sdbct2;
          d.y = 1.0;
       }*/
       /*float e = mod(iTime, 5.);
       if(e >= 4.){
           if(cute < d.x){
              d.x = cute;
              d.y = 3.0;
           }
           if(cute2 < d.x){
              d.x = cute2;
              d.y = 3.0;
           }
       }*/
      /* if(scap < d.x){
          d.x = scap;
          d.y = 3.0;
       }*/
       if(cut < d.x){
          d.x = cut;
          d.y = 1.0;
       }
       /*if(sdbct < d.x){
          d.x = sdbct;
          d.y = 3.0;
       }
       
        /*if(stepld < d.x){
          d.x = stepld;
          d.y = 4.0;
       }*/
       
       /*if(sbhd < d.x){
          d.x = sbhd;
          d.y = 4.0;
       }*/
           
     /* if(pd < d.x){
          d.x = pd;
          d.y = 2.0;
       }  */   
       
           
     
    return d;
    
    
}

vec3 RayMarch2(vec3 eye, vec3 viewRayDirection, vec2 uv){
    vec3 t = vec3(0.);
    float max = -100000.0;
    vec2 dd;
    float depth = 0.0, end = 10.0;
    for (int i = 0; i < 100; i++) {
        t.yz = GetDist(eye + t.x * viewRayDirection, uv).xy;
        
                    
        if (t.y < 0.01)break;
                
        t.x += t.y;
        
        if (t.x >= 1000.0)break;
        
    }
    if (t.x >= 1000.0)t.x = -1.0;
    
    return t;


}


float calcAO( in vec3 pos, in vec3 nor )
{
	float occ = 0.0;
    float sca = 1.0;
    for( int i=0; i<5; i++ )
    {
        float hr = 0.01 + 0.12*float(i)/4.0;
        vec3 aopos =  nor * hr + pos;
        float dd = GetDist( aopos , vec2(1.0)).y;
        occ += -(dd-hr)*sca;
        sca *= 0.95;
    }
    return clamp( 1.0 - 3.0*occ, 0.0, 1.0 );    
}

vec3 GetLightM(vec3 p, vec2 uv, vec3 lg, out vec3 n, vec3 ro, vec3 rd){
    vec3 lightpos = lg;
    //lightpos.xz += vec2(sin(iTime), cos(iTime));
    vec3 l = normalize(lightpos-p);
    
    vec2 d = GetDist(p, uv);
    vec2 e = vec2(0.01, 0);
    
    
    
    
    n = d.x - vec3(
        GetDist(p-e.xyy, uv).x,
        GetDist(p-e.yxy, uv).x,
        GetDist(p-e.yyx, uv).x);


    n = normalize(n);
        
    
    float occ = calcAO(p, n);
            
    float dif = clamp(dot(n, l), .0, 1.);
    dif += occ*0.3;
 
    vec3 dd = RayMarch2(p+n*.01, l, uv);
    p = ro + reflect(n, l) * dd.x;
    
        
     
    return vec3(dif);

}
vec3 lightDir = normalize( vec3(0.5,0.6,0.) );
const mat2 m2 = mat2( 0.60, -0.80, 0.80, 0.60 );
//ref: https://www.shadertoy.com/view/Msdfz8
vec3 Cloud(vec3 bgCol,vec3 ro,vec3 rd,vec3 cloudCol,float spd)
{
    vec3 col = bgCol;
    float t = iTime * 0.15* spd;
    vec2 sc = ro.xz + rd.xz*((3.)*40000.0-ro.y)/rd.y;
    vec2 p = 0.00002*sc;
    float f = 0.0;
  	float s = 0.5;
  	float sum =0.;
  	for(int i=0;i<5;i++){
    	p += t;t *=1.5;
    	f += s*textureLod( iChannel0, p/256.0, 0.0).x; p = m2*p*2.02;
    	sum+= s;s*=0.6;
  	}
    float val = f/sum; 
    col = mix( col, cloudCol, smoothstep(0.5,0.8,val) );
    return col;
}
vec3 RayMarchCloud(vec3 ro,vec3 rd){
    vec3 col = vec3(0.0,0.0,0.0);  
    float sundot = clamp(dot(rd,lightDir),0.0,1.0);
    
     // sky      
    col = vec3(0.2,0.5,0.85)*1.1 - rd.y*rd.y*0.5;
    col = mix( col, 0.85*vec3(0.7,0.75,0.85), pow( 1.0-max(rd.y,0.0), 4.0 ) );
    // sun
    col += 0.25*vec3(1.0,0.7,0.4)*pow( sundot,5.0 );
    col += 0.25*vec3(1.0,0.8,0.6)*pow( sundot,64.0 );
    col += 0.4*vec3(1.0,0.8,0.6)*pow( sundot,512.0 );
    // clouds
    col = Cloud(col,ro,rd,vec3(1.0,0.95,1.0),1.);
            // .
    col = mix( col, 1.5*vec3(0.0,0.5,1.0), pow( 1.0-max(rd.y,0.0), 16.0 ) );
    return col;
}


vec3 Eye(vec3 p, vec4 sp, vec2 uv){

    vec3 col = vec3(1.0);
    vec2 pp = p.xy-sp.xy;
    float r = sqrt(dot(pp, pp));
    float a = atan(pp.y, pp.x);

    float ss = 0.5+0.5*sin(5.0*iTime);
    float anim = 1.0 + 0.1*ss*clamp(r, 0.0, 1.0);
    r*=anim;

    if(length(p.xy-sp.xy) < .2){

        col =  vec3(0.0, 0.3, 0.4);

        float f = fbm(5.0*pp);
        col = mix(col, vec3(0.2, 0.5, 0.4), f);

        f = 1. - smoothstep(0.1, 0.15, r);
        col = mix(col, texture(iChannel1, uv).rgb, f);

        /*a += fbm(20.0*pp)*0.05;

        f = smoothstep(0.15, 0.2, fbm(vec2(6.0*r, 20.0*a)) );
        col = mix(col, vec3(1.0), f);

        f = smoothstep(0.4, 0.9, fbm(vec2(10.0*r, 15.0*a)) );
        col *= 1.0-0.5*f;

        f = smoothstep(0.6, 0.8, r);
        col *= 1.0 - 0.5*f;

        f = smoothstep(0.2, 0.25, r);
        col *= f;

        f = 1.0-smoothstep(0.0, 0.5, length(pp - vec2(0.24, 0.2)) );
        col += vec3(1.0, 0.9, 0.8)*f*0.9;

        f = smoothstep(0.7, 0.8, r);
        col = mix(col, vec3(1.0), f);*/

    }
    
    return col;

}


void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Normalized pixel coordinates (from 0 to 1)
    vec2 uv = fragCoord/iResolution.xy;
    uv -= 0.5;
    uv /= vec2(iResolution.y / iResolution.x, 1);
    // Time varying pixel color
    vec3 col = vec3(0.0);
    vec3 lightpos = vec3(0.0, 6., 0.);
    
    //iMouse.xy / iResolution.xy
    vec3 lookat = vec3((iMouse.x/iResolution.x)*100.0, (iMouse.y/iResolution.y)*100.0, 30.0);//vec3(0.0, 1.5, 0.0);
    vec3 ro = vec3(0, 7.0, -4.);
    /*ro.z += iTime;//mod(iTime, 20.0)-4.;
    lightpos.z += iTime;
    lookat.z += iTime;
    /*if(ro.z > 3.0){
        ro.z = 3.0;
        
    }*/
    
    
    float zoom = 1.0;    
    vec3 f = normalize(lookat-ro),
        r = normalize(cross(vec3(0,1,0), f)),
        u = cross(f, r),
        c = ro + f * zoom,
        i = c + uv.x * r + uv.y * u,
        rd = normalize(i-ro);
        
    
        
     vec3 t;
     float dO, dif, dif2, dif3, difT;
     vec3 p;
     
  
     
     t = RayMarch2(ro, rd, uv);
     p = ro + rd * t.x;
     
    float alpha = 1.0;
    vec3 n1, n2, n3, np;
    
  
   vec3 diff2 = GetLightM(p, uv, lightpos, n1, ro, rd);
    
    //vec3 ph = phongIllumination(vec3(1.0), vec3(dif), vec3(1.0, 1.0, 0.0), 50.0, p, ro, uv, np);
    
     if(t.x > 0.0){
         
         
         
     
     
         
         
         if(t.z == 0.0){
             //col = ph * vec3(0.0, 1.0, 0.0);
         }
         else if(t.z == 1.0){
             col = diff2* vec3(1.0, 1.0, 1.0);
         }
         else if(t.z == 2.0){
         
             
             col = diff2*vec3(.9, .9, .9);
             
         
             alpha = 1.0;
             
             
             
             
         }
         else if(t.z == -1.0){
             
             vec4 sp = eyel;
             
             col = Eye(p, sp, uv+0.5)*diff2;
             
         }
         else if(t.z == -2.0){
             
             vec4 sp = eyer;
             
             col = Eye(p, sp, uv+0.5)*diff2;
             
         }
         else if(t.z == 3.0){
            
             col = diff2*vec3(0.88, 0.69, 0.65);
             float f = smoothstep(0.5, 1., fbm(p.xy*25.0));
             col =mix(col,  texture(iChannel0, vec2(fbm(p.xy*25.0)) ).rgb*0.2, f)  ;
             
             //#e2b1a8
             vec2 sl = sourcl;
             sl.y -= 0.4; 
             sl.x -= 0.4;
             float a = 1.0;
             
             vec2 q = p.xy - sl;
             
             for(int i = 0;i < 10;i++){
                 
                 if(length(q) < 0.1){
                     col = vec3(0.0)+noise(q.xy*2.0)*0.5;

                 }
                 float n = sl.x;
                 sl.xy *= Rot(-a*3.1415/180.0);
                 sl.x =  n+0.1;
                 a += 0.1;
                 q = p.xy - sl;
                 
             }
             
             
             sl = sourcr;
             sl.y -= 0.3; 
             sl.x -= 0.4;
             a = 1.0;
             
             q = p.xy - sl;
             
             for(int i = 0;i < 10;i++){
                 
                 if(length(q) < 0.1){
                     col = vec3(0.0)+noise(q.xy*2.0)*0.5;

                 }
                 float n = sl.x;
                 sl.xy *= Rot(-a*3.1415/180.0);
                 sl.x = n+0.1;
                 a += 0.1;
                 q = p.xy - sl;
                 
             }
             
                          
             
         
         
         }
         else if(t.z == 4.0){
             
             col = diff2*vec3(1.0, 1.0, 0.0);
             
         }
         else if(t.z == 5.0){
            
             
             float NdotL = max( 0., dot( n1, lightpos-p ) );
             float SpecularColor = 0.5;
            SpecularColor = SpecularColor + ( 1. - SpecularColor ) * pow( ( 1. - NdotL ),2. );
             
             col = vec3(0.0, 0.8, 1.0) + SpecularColor*0.0001;
             col = mix(col, diff2 * vec3(0.0, 0.8, 1.0), smoothstep(.0, 0.9, diff2));
             
             
             //col =  diff2*vec3(1.0, 0.0, 1.0);
             
         }
         
        /* vec2 j = uv*3.0;
         j.x += 0.0;
         j.y += .1;
         float sparkle = 1./dot(j,j);
                
         col += diff2*(sparkle*sin(mod(iTime*10.0, 3.1415))*0.01) ;
    */
       
         
         
         
     }
     else
     {
         col = RayMarchCloud( ro, rd);
     }
    
   
    // Output to screen
    fragColor = vec4((col),alpha);
}



