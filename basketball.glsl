vec4 tpos = vec4(0.0, 0.0, 0.0, 0.0);
float angle_tpos;
bool done = false;
float dist = 0.0;
vec4 sphball;

mat2 Rot(float a){
    float s = sin(a);
    float c = cos(a);
    return mat2(c, -s, s, c);
}

mat3 RotX(float a){
    float s = sin(a);
    float c = cos(a);
    
    return mat3(1., .0, .0, .0, c, -s, .0, s, c);

}

mat3 RotY(float a){
    float s = sin(a);
    float c = cos(a);
    
    return mat3(c, .0, s, .0, 1., .0, -s, .0, c);

}

mat3 RotZ(float a){
    float s = sin(a);
    float c = cos(a);
    
    return mat3(c, -s, .0, s, c, .0, .0, .0, 1.);

}


float opExtrusion(in vec3 p, in float d, in float h) {
    // d is the distance to the 2D shape using the x and y components of p
    vec2 w = vec2(d, abs(p.z) - h);
    return min(max(w.x, w.y), 0.0) + length(max(w, 0.0));
}


vec4 opElongate( in vec3 p, in vec3 h )
{
    //return vec4( p-clamp(p,-h,h), 0.0 ); // faster, but produces zero in the interior elongated box
    
    vec3 q = abs(p)-h;
    return vec4( max(q,0.0), min(max(q.x,max(q.y,q.z)),0.0) );
}


vec3 opTwist(in vec3 p, in float t) {
    float c = cos(t * p.y);
    float s = sin(t * p.y);
    mat2 m = mat2(c, -s, s, c);
    return vec3(m * p.xz, p.y);
}

float N21(vec2 p){
    p = fract(p*vec2(233.34, 851.73));
    p += dot(p, p+23.45);
    return fract(p.x*p.y);

}

vec2 N22(vec2 p){

    float n = N21(p);
    return vec2(n, N21(p+n));

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

float FBM( in vec3 p )
{
    float n = 0.0;
    n += 0.50000*Noise( p*1.0 );
    n += 0.25000*Noise( p*2.0 );
    n += 0.12500*Noise( p*4.0 );
    n += 0.06250*Noise( p*8.0 );
    n += 0.03125*Noise( p*16.0 );
    return n/0.984375;
}


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

float sdEllipsoid( vec3 p, vec3 r )
{
  float k0 = length(p/r);
  float k1 = length(p/(r*r));
  return k0*(k0-1.0)/k1;
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

vec3 add_tpos(vec3 tp){
    return vec3(tp.x+tpos.x, tp.y+tpos.y, tp.z+tpos.z);
    
}

vec3 rot_y(vec3 tp, float a){
    tp.xz *= Rot(a);
    return tp;

}

float sdfSegment(in vec3 p, in vec3 a, in vec3 b, in float r){

    
    float h = min(1.0, max(0.0, dot(p-a, b-a) / dot(b-a, b-a)));
   // a = a*exp(-4.*h);
    //r =r-sin(9.0*3.1415*h)*0.2;
    //r =r+exp(-4.*h);
    //r = r -4.*h*(1.-h)-0.2*cos(10.*h+4.*iTime);
    return length(p-a-(b-a)*h)-r;
    

}

float opRepSeg( in vec3 p, in vec3 c,  in vec3 a, in vec3 b, in float r)
{
    vec3 q = mod(p+0.5*c,c)-0.5*c;
    return sdfSegment(q, a, b, r);
}



float det( vec2 a, vec2 b ) { return a.x*b.y-b.x*a.y; }
vec4 sdBezier2( vec3 p, vec3 va, vec3 vb, vec3 vc )
{
  vec3 w = normalize( cross( vc-vb, va-vb ) );
  vec3 u = normalize( vc-vb );
  vec3 v =          ( cross( w, u ) );
  //----  
  vec2 m = vec2( dot(va-vb,u), dot(va-vb,v) );
  vec2 n = vec2( dot(vc-vb,u), dot(vc-vb,v) );
  vec3 q = vec3( dot( p-vb,u), dot( p-vb,v), dot(p-vb,w) );
  //----  
  float mn = det(m,n);
  float mq = det(m,q.xy);
  float nq = det(n,q.xy);
  //----  
  vec2  g = (nq+mq+mn)*n + (nq+mq-mn)*m;
  float f = (nq-mq+mn)*(nq-mq+mn) + 4.0*mq*nq;
  vec2  z = 0.5*f*vec2(-g.y,g.x)/dot(g,g);
//float t = clamp(0.5+0.5*(det(z,m+n)+mq+nq)/mn, 0.0 ,1.0 );
  float t = clamp(0.5+0.5*(det(z-q.xy,m+n))/mn, 0.0 ,1.0 );
  vec2 cp = m*(1.0-t)*(1.0-t) + n*t*t - q.xy;
  //----  
  float d2 = dot(cp,cp);
  return vec4(sqrt(d2+q.z*q.z), t, q.z, -sign(f)*sqrt(d2) );
}


vec4 opRepBezier( in vec3 p, in vec3 c, vec3 va, vec3 vb, vec3 vc )
{
    vec3 q = mod(p+0.5*c,c)-0.5*c;
    vec3 qa = mod(va+0.5*c,c)-0.5*c;
    vec3 qb = mod(vb+0.5*c,c)-0.5*c;
    vec3 qc = mod(vc+0.5*c,c)-0.5*c;
    return  sdBezier2(q, va, vb, vc);
}

float _line(vec2 p, vec2 a, vec2 b){

    vec2 pa = p-a;
    vec2 ba = b-a;
    float h = clamp(dot(pa, ba) / dot(ba, ba), 0., 1.);
    return length(pa - ba*h);
}

vec3 Line(vec2 p, vec2 a, vec2 b, float l1, float l2){

    float d = _line(p, a, b);
    float m = smoothstep(l1, l2, d);
    return vec3(m); 

}

// https://iquilezles.org/articles/smin
float smin( float a, float b, float k )
{
    float h = max(k-abs(a-b),0.0);
    return min(a, b) - h*h*0.25/k;
}

// https://iquilezles.org/articles/smin
float smax( float a, float b, float k )
{
    k *= 1.4;
    float h = max(k-abs(a-b),0.0);
    return max(a, b) + h*h*h/(6.0*k*k);
}

float sdCone( in vec3 p, in vec2 c, float h , float r)
{
  // c is the sin/cos of the angle, h is height
  // Alternatively pass q instead of (c,h),
  // which is the point at the base in 2D
  vec2 q = h*vec2(c.x/c.y,-1.0);
    
  vec2 w = vec2( length(p.xz), p.y );
  vec2 a = w - q*clamp( dot(w,q)/dot(q,q), 0.0, 1.0 );
  vec2 b = w - q*vec2( clamp( w.x/q.x, 0.0, 1.0 ), 1.0 );
  float k = sign( q.y );
  float d = min(dot( a, a ),dot(b, b));
  float s = max( k*(w.x*q.y-w.y*q.x),k*(w.y-q.y)  );
  return sqrt(d)*sign(s)-r;
}

vec2 sdSegment(vec3 p, vec3 a, vec3 b)
{
    vec3 pa = p-a, ba = b-a;
	float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
	return vec2( length( pa - ba*h ), h );
}


 
vec2 GetDist(vec3 p, vec2 uv){


    
    vec4 s2 = vec4(3.0 ,3.0, 7.0, 1.0);
    
        
    vec2 d=vec2(0.0);;
    float pd = dot(p, vec3(0.0, 1.0, 0.0));
    
   // pd -=  Noise3D(p+iTime*5.0)*0.2 ;
   // pd +=  0.1*(sin(1.0*p.x)+sin(1.0*p.z));
    
    //float fh = -0.1 + 0.05*(sin(2.0*p.x) + sin(2.0*p.z));
   // pd = p.y - fh;
    pd += -1.0+2.0*smoothstep(-0.5, 0.5, sin(1.0*p.x)+sin(1.0*p.y)+sin(1.0*p.z))*0.1;
    vec4 sph = sphball;
     
    sph.y = sin(sph.y*3.1415+iTime*10.0)*5.0+5.0;
    sph.z -= iTime*2.1;
    
    float x = -1.0 + 2.0*abs(fract(iTime*0.5)-0.5)*2.0;
    sph.x = x*10.5;
    
    float fh = -0.1 - 0.05*(sin(p.x*1.0)+sin(p.z*1.0));
    float gt = fract(iTime+0.1);
    float l = length((p.xyz-sph.xyz).xz);
     fh -= 0.4*sin(gt*10.0+l)*exp(-0.02*l*l )*exp(-0.02*gt ) * smoothstep(0.0, 0.1, gt);
    pd -=  fh;
    
        
    vec3 sp = p - sph.xyz;
    sp *= RotZ(iTime);
    sp *= RotY(iTime);
    sp *= RotX(iTime);
    //float sd = length(sp) - sph.w;
    vec3 dim = vec3(3.0);
    if(sph.y < 2.0){
        dim.x = iTime*0.001+4.0;
    }
    float sd = sdEllipsoid(sp, dim);
    
    
    
    vec3 q = p-sph.xyz;
    q *= RotZ(iTime);
    q *= RotY(iTime);
    q *= RotX(iTime);
    sd += smoothstep(0.0, 1.0, Line(q.xy, vec2(5.0, 0.0), vec2(-5.0, 0.0), 0.1, 0.05).x*0.1);
    sd += smoothstep(0.0, 1.0, Line(q.xy, vec2(0.0, -5.0), vec2(0.0, 5.0), 0.1, 0.05).x*0.1);

    //
    sd += smoothstep(0.0, 1.0, Line(q.xy, vec2(2.5, 2.5), vec2(-2.5, -2.5), 0.1, 0.05).x*0.1);
    sd += smoothstep(0.0, 1.0, Line(q.xy, vec2(-2.5, 2.5), vec2(2.5, -2.5), 0.1, 0.05).x*0.1);

    sd += smoothstep(0.0, 1.0, Line(q.xy, vec2(2.5, 2.5), vec2(-2.5, 2.5), 0.1, 0.05).x*0.1);
    sd += smoothstep(0.0, 1.0, Line(q.xy, vec2(2.5, -2.5), vec2(-2.5, -2.5), 0.1, 0.05).x*0.1);
    
    /*vec3 pm = vec3(x*10.5,sph.y, 14.0);
    pm.z -= iTime*2.10;
    d = sdfMan(p-pm);*/
    
        vec3 bp = p - vec3(25.0, -5.5, 0.0);
    bp.y -= mod(iTime*2.0, 20.0);
    vec3 c = vec3(50.0, 0.0, 25.0);
    bp = mod(bp+0.5*c,c)-0.5*c;;
    vec3 sc = vec3(5.0);
    sc -= mod(iTime*2.0*0.25, 5.0);
   // if(sc.x < 0.0)sc = vec3(5.0);
    float bd = Box(bp, sc, 0.5);
    bd += -1.0+2.0*smoothstep(-0.5, 0.5, sin(1.0*p.x)+sin(1.0*p.y)+sin(1.0*p.z))*0.1;
    
    float gr = smin(pd, bd, 0.5);
      
    //d.x = pd;
    //d.y = 2.0;
       
   // if(sd < d.x){
        d.x = sd;
        d.y = 5.0;
    
    //}
     
    if(gr < d.x){
        d.x = gr;
        d.y = 2.0;
    
    }
      
   // 
    return d;
    
    
}

vec3 RayMarch2(vec3 eye, vec3 viewRayDirection, vec2 uv){
    vec3 t = vec3(0.);
    float max = -100000.0;
    vec2 dd;
    float depth = 0.0, end = 10.0;
    for (int i = 0; i < 100; i++) {
        t.yz = GetDist(eye + t.x * viewRayDirection, uv).xy;
        
                    
        if (abs(t.y) < (t.y*0.01))break;
                
        t.x += t.y;
        
        if (t.x >= 100.0)break;
        
    }
    if (t.x >= 100.0)t.x = -1.0;
    
    return t;


}




float GetLight(vec3 p, vec2 uv, vec3 lg, out vec3 n){
    vec3 lightpos = lg;
    lightpos.xz += vec2(sin(iTime), cos(iTime));
    vec3 l = normalize(lightpos-p);
    
    vec2 d = GetDist(p, uv);
    vec2 e = vec2(0.01, 0);
    
    
   // float dd = d.x;
    n = d.x - vec3(
        GetDist(p-e.xyy, uv).x,
        GetDist(p-e.yxy, uv).x,
        GetDist(p-e.yyx, uv).x);
    
    n = normalize(n);
    
    float dif = clamp(dot(n, l), .0, 1.);
   // float dif = clamp(dot(n, l), 0., 1.);
   // vec3 dd = RayMarch2(p+n*.01, l, uv);
    //if(dd.x < length(lightpos-p))dif *= 0.1;
    return (dif) ;

}


float calcAO( in vec3 pos, in vec3 nor )
{
	float occ = 0.0;
    float sca = 1.0;
    for( int i=0; i<5; i++ )
    {
        float hr = 0.01 + 0.12*float(i)/4.0;
        vec3 aopos =  nor * hr + pos;
        float dd = GetDist( aopos , vec2(1.0)).x;
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
    dif += occ*0.2;
 
    vec3 dd = RayMarch2(p+n*.01, l, uv);
    p = ro + rd * dd.x;
    
    
    
    /*if(dd.x < length(lightpos)){
           if(dd.z == 2.0)
        { 
        
            vec3 col = vec3(length(uv))*texture(iChannel2, vec2(noise(uv+iTime))).rgb;
            
        }
    vec3(0.0, 1.0, 0.0)-p/2.0
    }*/
    
     
    if((dd.x < length(lightpos-p)) && (dd.x > 0.0)){
        if(dd.z == -1.0){
            
             vec3 col = vec3(dif) ;
             return col*0.2;
            

        }
        else if(dd.z == 1.0)
        {
             
             return vec3(dif)*0.2;
        }
        else if(dd.z == 2.0)
        {
             
             return vec3(dif)*0.2;
        }
        else if(dd.z == 3.0)
        {
             //vec3 col = texture(iChannel1, p.xz*0.5).rgb * n.y;
             vec3 col = vec3(dif)*vec3(0.88, 0.69, 0.65); ;
             return col*0.2;
        }
        else if(dd.z == 5.0){
             vec3 col = vec3(1.0, 0.0, 0.0) * vec3(dif);
             return col;
             
         }
        else
            return vec3(dif);
            
     }
     else
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
    	f += s* textureLod( iChannel0, p/256.0, 0.0).x; p = m2*p*2.02;
    	sum+= s;s*=0.8;
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

vec3 Eye(vec3 p, vec4 sp){

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
        col = mix(col, vec3(0.0, 0.0, 0.0), f);

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
    sphball = vec4(0.0,3.0, 15.0,3.0);
    
    //iMouse.xy / iResolution.xy
    vec3 lookat = vec3((iMouse.x/iResolution.x)*10.0, (iMouse.y/iResolution.y)*30.0, 30.0);//vec3(0.0, 1.5, 0.0);
    vec3 ro = vec3(-5, 7.0, 1.0);
    vec4 sph = sphball;
    float speed = 2.1;
    ro.z -=iTime*speed;
    lookat.z -= iTime*speed;
    sph.z -= iTime*speed;
    float an = (iMouse.x/iResolution.x)*10.0;
    ro -= sph.xyz;
    lookat -= sph.xyz;
    ro *= RotY(an);
    lookat *= RotY(an);
    ro += sph.xyz;
    lookat += sph.xyz;
    
    vec3 lightpos = vec3(0.0, 30., -30.);
    lightpos.z -= iTime*speed;
    lightpos -= sph.xyz;
    lightpos *= RotY(an);
    lightpos += sph.xyz;
    //vec3 lightpos2 = vec3(0.0, 5., 16.);
    //lightpos2.z -= iTime;
        
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
   //vec3 diff3 = GetLightM(p, uv, lightpos2, n1, ro, rd);
    
        
     if(t.x > 0.0){
         
         
         
    
     
         
         
         if(t.z == 0.0){
            // col = ph * vec3(0.0, 1.0, 0.0);
         }
         else if(t.z == 1.0){
             col = diff2* vec3(1.0, 0.0, 0.0);
         }
         else if(t.z == 2.0){
         
                                    
            // col = diff2 * vec3(1.0, 0.5, 0.5);
            // col = 1.0-mix(diff2 * vec3(1.0, 0.5, 0.5), 1.0 - texture(iChannel0, p.xz).rgb*n1.y, 2.8);
             //col = mix(col, texture(iChannel1, p.xz).rgb*n1.y, 0.8);
             col = vec3(0.05, 0.09, 0.02);
             float f = -1.0+2.0*smoothstep(-0.5, 0.5, sin(11.0*p.x)+sin(11.0*p.z))*0.5;
             //float f =-1.0+2.0*smoothstep(-0.5, 0.5, sin(1.0*p.x)+sin(1.0*p.y)+sin(1.0*p.z))*0.1;
             col += f*0.05;
             col +=  diff2 * vec3(0.2, 1., 0.2);
             
         }
         else if(t.z == -1.0){
             
                          
         }
         else if(t.z == -2.0){
                          
             
         }
         else if(t.z == 3.0){
            
              col = diff2*vec3(0.88, 0.69, 0.65);
              float f = smoothstep(0.5, 1., fbm(p.xy*25.0));
              col =mix(col,  texture(iChannel0, vec2(fbm(p.xy*25.0)) ).rgb*0.2, f)  ;
         
         }
         else if(t.z == 5.0){
             col = vec3(1.0, 0.0, 0.0) * diff2;
             vec4 sph = sphball;
             sph.y = sin(sph.y*3.1415+iTime*10.0)*5.0+5.0;
            sph.z -= iTime*2.1;

            float x = -1.0 + 2.0*abs(fract(iTime*0.5)-0.5)/0.5;
            sph.x = x*10.5;
            
            vec3 q = p-sph.xyz;
            
            q *= RotZ(iTime);
            q *= RotY(iTime);
            q *= RotX(iTime);
            
            vec3 m = Line(q.xy, vec2(5.0, 0.0), vec2(-5.0, 0.0), 0.1, 0.05);
             m += Line(q.xy, vec2(0.0, -5.0), vec2(0.0, 5.0), 0.1, 0.05);
             
             //
             m += Line(q.xy, vec2(2.5, 2.5), vec2(-2.5, -2.5), 0.1, 0.05);
             m += Line(q.xy, vec2(-2.5, 2.5), vec2(2.5, -2.5), 0.1, 0.05);
             
             m += Line(q.xy, vec2(2.5, 2.5), vec2(-2.5, 2.5), 0.1, 0.05);
             m += Line(q.xy, vec2(2.5, -2.5), vec2(-2.5, -2.5), 0.1, 0.05);
             
             
             col += m;
             
         }
         
       
         
         
         
     }
     else
     {
         col = RayMarchCloud( ro, rd);
         
         
     }
    
   
    // Output to screen
    fragColor = vec4((col),alpha);
}



