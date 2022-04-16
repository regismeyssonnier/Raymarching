vec4 tpos = vec4(0.0, 0.0, 0.0, 0.0);
float angle_tpos;
bool done = false;
const int NUM_STEPS = 8;
const float PI	 	= 3.141592;
const float EPSILON	= 1e-3;
#define EPSILON_NRM (0.1 / iResolution.x)
#define AA

// sea
const int ITER_GEOMETRY = 3;
const int ITER_FRAGMENT = 5;
const float SEA_HEIGHT = 0.6;
const float SEA_CHOPPY = 4.0;
const float SEA_SPEED = 2.0;
const float SEA_FREQ = 0.17;
const vec3 SEA_BASE = vec3(0.2);//vec3(0.0,0.09,0.18);
const vec3 SEA_WATER_COLOR = vec3(0.2);//vec3(0.8,0.9,0.6)*0.6;
#define SEA_TIME (1.0 + iTime * SEA_SPEED)
const mat2 octave_m = mat2(1.6,1.2,-1.2,1.6);

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

mat2 Rot(float a){
    float s = sin(a);
    float c = cos(a);
    return mat2(c, -s, s, c);
}

float hash( vec2 p ) {
	float h = dot(p,vec2(127.1,311.7));	
    return fract(sin(h)*43758.5453123);
}
float noise( in vec2 p ) {
    vec2 i = floor( p );
    vec2 f = fract( p );	
	vec2 u = f*f*(3.0-2.0*f);
    return -1.0+2.0*mix( mix( hash( i + vec2(0.0,0.0) ), 
                     hash( i + vec2(1.0,0.0) ), u.x),
                mix( hash( i + vec2(0.0,1.0) ), 
                     hash( i + vec2(1.0,1.0) ), u.x), u.y);
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

// sea
float sea_octave(vec2 uv, float choppy) {
    uv += noise(uv);        
    vec2 wv = abs(sin(uv));
    vec2 swv = abs(cos(uv));    
    wv = mix(wv,swv,wv);
    return length(wv);//pow(1.0-pow(wv.x * wv.y,0.65),choppy);
}

float map_detailed(vec3 p) {
    float freq = SEA_FREQ;
    float amp = SEA_HEIGHT;
    float choppy = SEA_CHOPPY;
    vec2 uv = p.xz; uv.x *= 0.75;
    
    float d, h = 0.0;    
    for(int i = 0; i < ITER_FRAGMENT; i++) {        
    	d = sea_octave((uv+SEA_TIME)*freq,choppy);
    	d += sea_octave((uv-SEA_TIME)*freq,choppy);
        h += d * amp;        
    	uv *= octave_m; freq *= 1.88; amp *= 0.22;
        //choppy = mix(choppy,1.0,0.2);
    }
    return p.y - h;
}

float map(vec3 p) {
    float freq = SEA_FREQ;
    float amp = SEA_HEIGHT;
    float choppy = SEA_CHOPPY;
    vec2 uv = p.xz; uv.x *= 0.75;
    
    float d, h = 0.0;    
    for(int i = 0; i < ITER_GEOMETRY; i++) {        
    	d = sea_octave((uv+SEA_TIME)*freq,choppy);
    	d += sea_octave((uv-SEA_TIME)*freq,choppy);
        h += d * amp;        
    	uv *= octave_m; freq *= 1.88; amp *= 0.22;
        choppy = mix(choppy,1.0,0.2);
    }
    return p.y - h;
}

vec2 GetDist(vec3 p, vec2 uv){


    
    
    vec4 s = vec4(1.0,0.0, 2.0,1.0);
    vec4 s2 = vec4(2.0 ,0.0, 6.0, 0.2);
    vec4 s3 = vec4(3.0,7.0, 5.0,0.5);
    
    
    
    vec2 d=vec2(0.0);;
    float pd = dot(p, vec3(0.0, 1.0, 0.0));
    
    
    pd = map(p);
       
    d.x = pd;
    d.y = 2.0;
      
    return d;
    
    
}

vec3 RayMarch2(vec3 eye, vec3 viewRayDirection, vec2 uv){
    vec3 t = vec3(0.0, 0.0, 0.0);
   
    vec2 dd;
    float depth = 0.0, end = 10.0;
    for (int i = 0; i < 100; i++) {
        t.yz = GetDist(eye + t.x * viewRayDirection, uv).xy;
        
                    
        if (t.y < 0.001)break;
                
        t.x += t.y;
        
        if (t.x >= 50.0)break;
        
    }
    if (t.x >= 50.0)t.x = -1.0;
    
    return t;


}


// lighting
float diffuse(vec3 n,vec3 l,float p) {
    return pow(dot(n,l) * 0.4 + 0.6,p);
}
float specular(vec3 n,vec3 l,vec3 e,float s) {    
    float nrm = (s + 8.0) / (PI * 8.0);
    return pow(max(dot(reflect(e,n),l),0.0),s) * nrm;
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

// sky
vec3 getSkyColor(vec3 e, vec3 ori, vec3 dir) {
    //return RayMarchCloud(ori, dir);
    e.y = (max(e.y,0.0)*0.8+0.2)*0.8;
    return vec3(pow(1.0-e.y,2.0), 1.0-e.y, 0.6+(1.0-e.y)*0.4) * 1.1;
}

vec3 getSeaColor(vec3 p, vec3 n, vec3 l, vec3 eye, vec3 dist) {  
    float fresnel = clamp(1.0 - dot(n,-eye), 0.0, 1.0);
    fresnel = pow(fresnel,3.0) * 0.5;
        
        float dif = clamp(dot(n, l), .0, 1.);
    vec3 reflected = getSkyColor(reflect(eye,n), eye, dist);    
    vec3 refracted = SEA_BASE + vec3(dif)/*diffuse(n,l,80.0)*/ * SEA_WATER_COLOR * 0.12; 
   //vec3 refracted = SEA_BASE + diffuse(n,l,80.0) * SEA_WATER_COLOR * 0.12; 
    
    
    //vec3 color = vec3(dif) * vec3(0.0, 1.0, 1.0);
    
    vec3 color = mix(refracted,reflected,fresnel);
    
    float atten = max(1.0 - dot(dist,dist) * 0.001, 0.0);
    color += vec3(0.0,1.0, 1.0) * (p.y - SEA_HEIGHT) * 0.18 * atten;
    
    color += vec3(specular(n,l,eye,60.0));
    /**/
    
    
    return color;
}

// tracing
vec3 getNormal(vec3 p, float eps) {
    vec3 n;
    n.y = map_detailed(p);    
    n.x = map_detailed(vec3(p.x+eps,p.y,p.z)) - n.y;
    n.z = map_detailed(vec3(p.x,p.y,p.z+eps)) - n.y;
    n.y = eps;
    return normalize(n);
}




vec3  GetLightSea(vec3 p, vec2 uv, vec3 lg, vec3 ro){
    vec3 lightpos = lg;
    //lightpos.xz += vec2(sin(iTime), cos(iTime));
    vec3 l = normalize(lightpos-p);
    
    //vec2 d = GetDist(p, uv);
    //vec2 e = vec2(0.01, 0);
    
    
   vec3 ori = ro;//vec3(0.0,3.5,5.0);
    vec3 dir = normalize(vec3(uv.xy,1.0));// dir.z += length(uv) * 0.14;
   // dir = normalize(dir) * fromEuler(ang);
    
     //float h = Height(ori,dir,p);
     vec3 t = RayMarch2(ro, dir, vec2(1.0));
     if(t.x > 0.0){
     p = ro + dir * t.x;
    vec3 n = getNormal(p, dot(p-ori,p-ori) * EPSILON_NRM);
    
    return mix(
        getSkyColor(dir, ori, dir),    
       getSeaColor(p, n, l, dir, p-ori)+vec3(0.0, 0.1, 0.10),
    	smoothstep(0.0,-0.02,dir.y));
    
    float dif = dot(refract(l,n, 0.9), l);
    }
    else
    {
        return RayMarchCloud(ro, dir);
    }

}



/**
 * Lighting contribution of a single point light source via Phong illumination.
 * 
 * The vec3 returned is the RGB color of the light's contribution.
 *
 * k_a: Ambient color
 * k_d: Diffuse color
 * k_s: Specular color
 * alpha: Shininess coefficient
 * p: position of point being lit
 * eye: the position of the camera
 * lightPos: the position of the light
 * lightIntensity: color/intensity of the light
 *
 * See https://en.wikipedia.org/wiki/Phong_reflection_model#Description
 */
vec3 phongContribForLight(vec3 k_d, vec3 k_s, float alpha, vec3 p, vec3 eye,
                          vec3 lightPos, vec3 lightIntensity, vec2 uv) {
    
    vec2 d = GetDist(p, uv);
    vec2 e = vec2(0.01, 0);
    
    
    vec3 N = d.x - vec3(
        GetDist(p-e.xyy, uv).x,
        GetDist(p-e.yxy, uv).x,
        GetDist(p-e.yyx, uv).x);
    
    N = normalize(N);
    //vec3 N = estimateNormal(p);
    vec3 L = normalize(lightPos - p);
    vec3 V = normalize(eye - p);
    vec3 R = normalize(reflect(-L, N));
    
    float dotLN = dot(L, N);
    float dotRV = dot(R, V);
    
    if (dotLN < 0.0) {
        // Light not visible from this point on the surface
        return vec3(0.0, 0.0, 0.0);
    } 
    
    if (dotRV < 0.0) {
        // Light reflection in opposite direction as viewer, apply only diffuse
        // component
        return lightIntensity * (k_d * dotLN);
    }
    return lightIntensity * (k_d * dotLN + k_s * pow(dotRV, alpha));
}

/**
 * Lighting via Phong illumination.
 * 
 * The vec3 returned is the RGB color of that point after lighting is applied.
 * k_a: Ambient color
 * k_d: Diffuse color
 * k_s: Specular color
 * alpha: Shininess coefficient
 * p: position of point being lit
 * eye: the position of the camera
 *
 * See https://en.wikipedia.org/wiki/Phong_reflection_model#Description
 */
vec3 phongIllumination(vec3 k_a, vec3 k_d, vec3 k_s, float alpha, vec3 p, vec3 eye, vec2 uv) {
    const vec3 ambientLight = 0.5 * vec3(1.0, 1.0, 1.0);
    vec3 color = ambientLight * k_a;
    
    vec3 light1Pos = vec3(4.0 * sin(iTime),
                          2.0,
                          4.0 * cos(iTime));
    vec3 light1Intensity = vec3(0.4, 0.4, 0.4);
    
    color += phongContribForLight(k_d, k_s, alpha, p, eye,
                                  light1Pos,
                                  light1Intensity,
                                  uv);
    
    vec3 light2Pos = vec3(2.0 * sin(0.37 * iTime),
                          2.0 * cos(0.37 * iTime),
                          2.0);
    vec3 light2Intensity = vec3(0.4, 0.4, 0.4);
    
    color += phongContribForLight(k_d, k_s, alpha, p, eye,
                                  light2Pos,
                                  light2Intensity,
                                  uv);    
    return color;
}


void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Normalized pixel coordinates (from 0 to 1)
    vec2 uv = fragCoord/iResolution.xy;
    uv -= 0.5;
    uv /= vec2(iResolution.y / iResolution.x, 1);
    // Time varying pixel color
    vec3 col = vec3(0.0);
    
    //iMouse.xy / iResolution.xy
    
    vec3 ro = vec3(0, 2.0, -4.);
    //ro.xz *= Rot(iTime);
   //ro.y = sin(iTime)*0.5+0.5;
    vec3 lookat = vec3(0., 0.0, -1.0);//vec3(0.0, 1.5, 0.0);
   ///lookat *= RotY(iTime);
    float zoom = 1.0;    
    vec3 f = normalize(lookat-ro),
        r = normalize(cross(vec3(0,1,0), f)),
        u = cross(f, r),
        c = ro + f * zoom,
        i = c + uv.x * r + uv.y * u,
        rd = normalize(i-ro);
        
    
    vec3 p;
     
    vec3 color = GetLightSea(p, uv, vec3(10.0, 5, -6), ro) ;
    col += color;
        
    
         
     
     
   
    // Output to screen
    fragColor = vec4((col),1.0);
}



