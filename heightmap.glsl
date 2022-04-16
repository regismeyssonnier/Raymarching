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
}
*/

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

float WaterMap( vec3 pos ) {
    return FBM( vec3( pos.xz, iTime*0.3 )) * 1.;
}
vec3 WaterNormal(vec3 pos,float rz){
    float EPSILON =rz*rz* 0.002;
    vec3 dx = vec3( EPSILON, 0.,0. );
    vec3 dz = vec3( 0.,0., EPSILON );
      
    vec3  normal = vec3( 0., 1., 0. );
    float bumpfactor = 0.3 * pow(1.-clamp((rz)/1000.,0.,1.),6.);//
    
    normal.x = -bumpfactor * (WaterMap(pos + dx) - WaterMap(pos-dx) ) / (2. * EPSILON);
    normal.z = -bumpfactor * (WaterMap(pos + dz) - WaterMap(pos-dz) ) / (2. * EPSILON);
    return normalize( normal ); 
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

float noi(vec3 p){ 
  vec3 f=floor(p),s=vec3(7,157,113);
  p-=f; vec4 h=vec4(0,s.yz,s.y+s.z)+dot(f,s);;
  p=p*p*(3.-2.*p);
  h=mix(fract(sin(h)*43758.5),fract(sin(h+s.x)*43758.5),p.x);
  h.xy=mix(h.xz,h.yw,p.y);
  return mix(h.x,h.y,p.z);  
}
vec2 add = vec2(1.0, 0.0);
#define HASHSCALE1 .1031
#define HASHSCALE3 vec3(.1031, .1030, .0973)
#define HASHSCALE4 vec4(1031, .1030, .0973, .1099)

float Hash12(vec2 p)
{
	vec3 p3  = fract(vec3(p.xyx) * HASHSCALE1);
    p3 += dot(p3, p3.yzx + 19.19);
    return fract((p3.x + p3.y) * p3.z);
}
vec2 Hash22(vec2 p)
{
	vec3 p3 = fract(vec3(p.xyx) * HASHSCALE3);
    p3 += dot(p3, p3.yzx+19.19);
    return fract((p3.xx+p3.yz)*p3.zy);

}

float Noise( in vec2 x )
{
    vec2 p = floor(x);
    vec2 f = fract(x);
    f = f*f*(3.0-2.0*f);
    
    float res = mix(mix( Hash12(p),          Hash12(p + add.xy),f.x),
                    mix( Hash12(p + add.yx), Hash12(p + add.xx),f.x),f.y);
    return res;
}

vec2 Noise2( in vec2 x )
{
    vec2 p = floor(x);
    vec2 f = fract(x);
    f = f*f*(3.0-2.0*f);
    float n = p.x + p.y * 57.0;
   vec2 res = mix(mix( Hash22(p),          Hash22(p + add.xy),f.x),
                  mix( Hash22(p + add.yx), Hash22(p + add.xx),f.x),f.y);
    return res;
}
const mat2 rotate2D = mat2(1.3623, 1.7531, -1.7131, 1.4623);
float Terrain( in vec2 p)
{
	vec2 pos = p*0.05;
	float w = (Noise(pos*.25)*0.75+.15);
	w = 66.0 * w * w;
	vec2 dxy = vec2(0.0, 0.0);
	float f = .0;
	for (int i = 0; i < 5; i++)
	{
		f += w * Noise(pos);
		w = -w * 0.4;	//...Flip negative and positive for variation
		pos = rotate2D * pos;
	}
	float ff = Noise(pos*.002);
	
	f += pow(abs(ff), 5.0)*275.-5.0;
	return f;
}

float TerrainRM(in vec2 p){

    vec2 q = p * 0.05;
    float w = Noise(q);
    w = 50.0*w;
    float f = 0.0;
    for(int i = 0;i < 5;i++){
        f += w * Noise(q*0.85);
        w = -w * 0.4;
        q = rotate2D * q;
    
    }
    
    float ff = Noise(q * 0.002);
    f += ff*250.-80.0;
   
    return f;


}

float TerrainRM2(in vec2 p){

    vec2 q = p * 0.05;
    float w = Noise(q);
    w = 50.0*w;
    float f = 0.0;
    for(int i = 0;i < 5;i++){
        f += w * Noise(q*0.85);
        w = -w * 0.4;
        q = rotate2D * q;
    
    }
    
    float ff = Noise(q * 0.002);
   f += ff*250.-80.0;
   
    
    
    for (int i = 0; i < 6; i++)
	{
		f += w * Noise(q*0.85);
		w =  - w * 0.4;
		q = rotate2D * q;
	}

    return f;

}
// High def version only used for grabbing normal information.
float Terrain2( in vec2 p)
{
	// There's some real magic numbers in here! 
	// The Noise calls add large mountain ranges for more variation over distances...
	vec2 pos = p*0.05;
	float w = (Noise(pos*.25)*0.75+.15);
	w = 66.0 * w * w;
	vec2 dxy = vec2(0.0, 0.0);
	float f = .0;
	for (int i = 0; i < 5; i++)
	{
		f += w * Noise(pos);
		w =  - w * 0.4;	//...Flip negative and positive for varition	   
		pos = rotate2D * pos;
	}
	float ff = Noise(pos*.002);
	f += pow(abs(ff), 5.0)*275.-5.0;
	
/*
	treeCol = Trees(p);
	f += treeCol;
	if (treeCol > 0.0) return f;

	*/
	// That's the last of the low resolution, now go down further for the Normal data...
	for (int i = 0; i < 6; i++)
	{
		f += w * Noise(pos);
		w =  - w * 0.4;
		pos = rotate2D * pos;
	}
	
	
	return f;
}

vec2 GetDist2(vec3 p, vec2 uv){


    vec4 s = vec4(-3.0,3.0, 9.0,1.0);
    vec4 s2 = vec4(3.0 ,3.0, 7.0, 1.0);
    
    
    
    
    vec2 d=vec2(0.0);;
    float pd = dot(p, vec3(0.0, 1.0, 0.0));
    
    pd = p.y - TerrainRM2(p.xz);
    
    if(p.y <= 0.0){
        pd = p.y - Noise3D(p*0.7)*0.5 - sin(p.x+iTime*2.0)*0.2-cos(p.z+iTime*2.0)*0.1+sin(p.y*3.0+iTime)*0.1;
    }
    
     
  /*  pd = pd - noi(p*0.9)*1.2 ;
    pd -= noi(p*3.0)*0.01;
   pd -= noi(p*10.0)*0.1;*/
   
    
    vec3 sp = (p)-s.xyz;
    float sd = length((sp))-s.w;
       
       
    
    vec3 sp2 = (p)-s2.xyz;
    //sp2.y += sin(iTime)*0.1;
    float sd2 = length((sp2))-s2.w;
   
    
   
      
      if(sd < sd2){
          d.x = sd;
          d.y = -1.0;
       }
       else
       {
           d.x = sd2;
          d.y = -1.0;
       }
       
       if(pd < d.x){
          d.x = pd;
          d.y = 2.0;
       }
       
     
       
       
      
   // 
    return d;
    
    
}

vec2 GetDist(vec3 p, vec2 uv){


    vec4 s = vec4(-3.0,3.0, 9.0,1.0);
    vec4 s2 = vec4(3.0 ,3.0, 7.0, 1.0);
    
    
    
    
    vec2 d=vec2(0.0);;
    float pd = dot(p, vec3(0.0, 1.0, 0.0));
    
    pd = p.y - TerrainRM(p.xz);
    
    if(p.y <= 0.0){
        pd = p.y - Noise3D(p*0.7)*0.5 - sin(p.x+iTime*2.0)*0.2-cos(p.z+iTime*2.0)*0.1+sin(p.y*3.0+iTime)*0.1;
    }
    
     
  /*  pd = pd - noi(p*0.9)*1.2 ;
    pd -= noi(p*3.0)*0.01;
   pd -= noi(p*10.0)*0.1;*/
   
    
    vec3 sp = (p)-s.xyz;
    float sd = length((sp))-s.w;
       
       
    
    vec3 sp2 = (p)-s2.xyz;
    //sp2.y += sin(iTime)*0.1;
    float sd2 = length((sp2))-s2.w;
   
    
   
      
      if(sd < sd2){
          d.x = sd;
          d.y = -1.0;
       }
       else
       {
           d.x = sd2;
          d.y = -1.0;
       }
       
       if(pd < d.x){
          d.x = pd;
          d.y = 2.0;
       }
       
     
       
       
      
   // 
    return d;
    
    
}

vec3 RayMarch2(vec3 eye, vec3 viewRayDirection, vec2 uv){
    vec3 t = vec3(0.);
    float max = -100000.0;
    vec2 dd;
    float depth = 0.0, end = 10.0, delta=0.0, v;
    for (int i = 0; i < 150; i++) {
        t.yz = GetDist(eye + t.x * viewRayDirection, uv).xy;
       
        
        
            
        if (t.y < 0.01) {
           
           break;
        }
        //t.x += t.y;
        
        

        if (t.x >= 10000.0) {
            break;
        }
       
        if(0.01>= 0.3*t.y)v= 0.01;
        else v = 0.3*t.y;
       delta = v + (t.x*0.0065);
		 t.x += delta;
    }
    if (t.x >= 10000.0)t.x = -1.0;
        
       
    
    return t;


}
 
vec3 RayMarch3(vec3 eye, vec3 viewRayDirection, vec2 uv){
    vec3 t = vec3(0.);
    float max = -100000.0;
    vec2 dd;
    float depth = 0.0, end = 10.0, delta=0.0, v;
    for (int i = 0; i < 150; i++) {
        t.yz = GetDist2(eye + t.x * viewRayDirection, uv).xy;
       
        
        
            
        if (t.y < 0.01) {
           
           break;
        }
        //t.x += t.y;
        
        

        if (t.x >= 10000.0) {
            break;
        }
       
        if(0.01>= 0.3*t.y)v= 0.01;
        else v = 0.3*t.y;
       delta = v + (t.x*0.0065);
		 t.x += delta;
    }
    if (t.x >= 10000.0)t.x = -1.0;
        
       
    
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


vec3 Mat2(vec3 p, vec3 n1, vec3 diff2, vec3 lightpos){
    vec3 col = mix(vec3(1.0)*diff2, vec3(Noise(p.xz)*0.1), 0.5);
    
   
   
    if(p.y <= 0.0){
                 
             col = vec3(0.0, 0.8, 1.0) ;
             col += mix(vec3(Noise(p.xz)*0.1) , diff2 * vec3(0.0, 0.8, 1.0), 0.6);
    }
    else if(p.y < 1.45){
         col =  vec3(1.0, 1.0,0.0)*diff2;
    }
    else if( (p.y < 10.35) && (n1.y > 0.9)){
        col = vec3(0.0,1.0, 0.0)*diff2;
       
    }
   
    else if ((p.y > 80.0)&& (n1.y > 0.5))
    {
        //col = mix(vec3(Noise(p.xz)*0.1)*vec3(1.0), vec3(1.0)*diff2+0.3, 0.9);
        //col += 0.3;
        col =  vec3(1.0)*diff2;
       
    }
    else
    {
       col += vec3(150. ,75., 21.)/255.*diff2;
    
    }
    
    
    
    

   
    
    
       
    
    return clamp(col, 0., 1.);
            
}

float calcSoftshadow( in vec3 ro, in vec3 rd, in float mint, in float tmax )
{
	float res = 1.;
    float t = mint;
    for( int i=0; i<16; i++ )
    {
		float h = GetDist2( ro + rd*t , vec2(1.0)).y;
        res = min( res, 10.0*h/t );
        t += clamp( h, 0.01, 0.10 );
        if( h<0.0001 || t>tmax ) break;
    }
    return clamp( res, 0.0, 1.0 );
}

float calcAO( in vec3 pos, in vec3 nor )
{
	float occ = 0.0;
    float sca = 1.0;
    for( int i=0; i<5; i++ )
    {
        float hr = 0.01 + 0.12*float(i)/4.0;
        vec3 aopos =  nor * hr + pos;
        float dd = GetDist2( aopos , vec2(1.0)).y;
        occ += -(dd-hr)*sca;
        sca *= 0.95;
    }
    return clamp( 1.0 - 3.0*occ, 0.0, 1.0 );    
}

vec3 GetLightM(vec3 p, vec2 uv, vec3 lg, out vec3 n, vec3 ro, vec3 rd){
    vec3 lightpos = lg;
    //lightpos.xz += vec2(sin(iTime), cos(iTime));
    //p = p*2.0;
    vec3 l = normalize(lightpos-p);
    
    vec2 d = GetDist2(p, uv);
    
    vec2 e = vec2(0.01, 0);
    
    n = d.x - vec3(
        GetDist2(p-e.xyy, uv).x,
        GetDist2(p-e.yxy, uv).x,
        GetDist2(p-e.yyx, uv).x);
    
    
    n = normalize(n);
    
    float occ = calcAO(p, n);
    
    
            
    float dif = clamp(dot(n, l), .0, 1.);
    //dif *= calcSoftshadow(p, lightpos, 0.01, 10000. );
    dif += occ*0.15;
    
    //
 
    vec3 dd = RayMarch2(p+n*.01, l, uv);
    p = ro + reflect(n, l) * dd.x;
    
    
    
         
    if(dd.y < length(lightpos-p) && dd.x > 0.0){
        if(dd.z == -1.0){
            vec3 colXZ = texture(iChannel0, p.xz*0.1).rgb;
             vec3 colYZ = texture(iChannel0, p.yz*0.1).rgb;
             vec3 colXY = texture(iChannel0, p.xy*0.1).rgb;

             n = abs(n);

             n *= pow(n, vec3(20));
             n /= n.x+n.y+n.z;

             vec3 col = colYZ * n.x + colXZ * n.y + colXY*n.z;

             uv = vec2(atan(p.x, p.z)/6.2832+.5, p.y/3.+0.5);
             vec4 st = texture(iChannel0, uv);

             col = vec3(dif) * mix(col, st.rgb, st.a);
             return col;

        }
        else if(dd.z == 2.0)
        {
            vec3 col;
            
                        
            col = Mat2(p, n, vec3(dif), lightpos);
            
            
            return col;
            
            
            
        }
        else if(dd.z == 3.0)
        {
            vec3 colXZ = texture(iChannel1, p.xz*0.1).rgb;
             vec3 colYZ = texture(iChannel1, p.yz*0.1).rgb;
             vec3 colXY = texture(iChannel1, p.xy*0.1).rgb;
             
             vec3 col = colYZ * n.x + colXZ * n.y + colXY*n.z;  
             col *= vec3(1.0, 0.0, 0.0);
             return col;
        }
        else if(dd.z == 4.0){
             vec3 colXZ = texture(iChannel2, p.xz*0.1).rgb;
             vec3 colYZ = texture(iChannel2, p.yz*0.1).rgb;
             vec3 colXY = texture(iChannel2, p.xy*0.1).rgb;
             
             
             vec3 col = colYZ * n.x + colXZ * n.y + colXY*n.z;
             
             col -= vec3(.3);
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
                          vec3 lightPos, vec3 lightIntensity, vec2 uv, out vec3 N) {
    
    vec2 d = GetDist(p, uv);
    vec2 e = vec2(0.01, 0);
    
    
    N = d.x - vec3(
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
vec3 phongIllumination(vec3 k_a, vec3 k_d, vec3 k_s, float alpha, vec3 p, vec3 eye, vec2 uv, out vec3 n) {
    const vec3 ambientLight = 0.5 * vec3(1.0, 1.0, 1.0);
    vec3 color = ambientLight * k_a;
    
    vec3 light1Pos = vec3(4.0 * sin(iTime),
                          2.0,
                          4.0 * cos(iTime));
    vec3 light1Intensity = vec3(0.4, 0.4, 0.4);
    
    color += phongContribForLight(k_d, k_s, alpha, p, eye,
                                  light1Pos,
                                  light1Intensity,
                                  uv,
                                  n);
    
    vec3 light2Pos = vec3(2.0 * sin(0.37 * iTime),
                          2.0 * cos(0.37 * iTime),
                          2.0);
    vec3 light2Intensity = vec3(0.4, 0.4, 0.4);
    
    color += phongContribForLight(k_d, k_s, alpha, p, eye,
                                  light2Pos,
                                  light2Intensity,
                                  uv,
                                  n);    
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
    vec3 lightpos = vec3(0.0, 100., -50.);
    
    //iMouse.xy / iResolution.xy
    vec3 lookat = vec3((iMouse.x/iResolution.x)*10.0, (iMouse.y/iResolution.y)*30.0, 30.0);//vec3(0.0, 1.5, 0.0);
    vec3 ro = vec3(0, 20.0, -50.);
    float speed = 50.0;
    
    //lightpos.z+= iTime*speed;
    
    float an = (iMouse.x/iResolution.x)*10.0;
    ro -= lookat;
    lookat -= lookat;
    ro *= RotY(an);
    lookat *= RotY(an);
    ro += lookat;
    lookat += lookat;
    
    //vec3 lightpos = vec3(0.0, 30., -30.);
   
    lightpos -= lookat;
    lightpos *= RotY(an);
    lightpos += lookat;
    
     lightpos.z += iTime*speed;
    ro.z += iTime*speed;
    lookat.z += iTime*speed;
    
    vec3 p;
    float h = 0.0;
    for (int i = 0; i < 4; i++)
	{
        
        h += Terrain(ro.xz) ;
        //
        
    }
    h+= 50.0;
    ro.y += h;
    lightpos.y += h;
    //lightpos.y += h;
    if(h <= 0.0){
        lookat.y = ro.y +10.0 ;
    
    }
    else
    {
        lookat.y = ro.y -10.0 ;
    
    }
    
    if(ro.y < 0.0){
        ro.y = 5.0;
    
    }
    
    float zoom = 2.0;    
    vec3 f = normalize(lookat-ro),
        r = normalize(cross(vec3(0,1,0), f)),
        u = cross(f, r),
        c = ro + f * zoom,
        i = c + uv.x * r + uv.y * u,
        rd = normalize(i-ro);
        
    
        
     vec3 t;
     float dO, dif, dif2, dif3, difT;
     //vec3 p;
     
  
     
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
             col = diff2* vec3(1.0, 0.0, 0.0);
         }
         else if(t.z == 2.0){
         
             
                
            col = Mat2(p, n1, diff2, lightpos);
            
            
         
             alpha = 0.1;
             
            
             
             
         }
         else if(t.z == -1.0){
             vec3 colXZ = texture(iChannel0, p.xz*0.1).rgb;
             vec3 colYZ = texture(iChannel0, p.yz*0.1).rgb;
             vec3 colXY = texture(iChannel0, p.xy*0.1).rgb;
             
             
             n1 = abs(n1);
             
             n1 *= pow(n1, vec3(20));
             n1 /= n1.x+n1.y+n1.z;
             
             col = colYZ * n1.x + colXZ * n1.y + colXY*n1.z;
             
             uv = vec2(atan(p.x, p.z)/6.2832+.5, p.y/3.+0.5);
             vec4 st = texture(iChannel0, uv);
            
             
             col = diff2 * mix(col, st.rgb, st.a);
         }
         else if(t.z == 3.0){
            
             
             float NdotL = max( 0., dot( n1, vec3(0.0, 6., 5.)-p ) );
             float SpecularColor;
            SpecularColor = SpecularColor + ( 1. - SpecularColor ) * pow( ( 1. - NdotL ), 5. );
             
             col = diff2*vec3(.003, .005, .005);
             col = mix(diff2 * vec3(0.5, 0.5, 0.5)+SpecularColor*0.001, 
                       vec3(0., 0.5, 5.0)*diff2+vec3(0.2,0.5, 0.5)*SpecularColor*0.0001, 0.7 );
             col += normalize(vec3(0.1, 0.1, 0.1));
         
         
         }
         else if(t.z == 4.0){
             vec3 colXZ = texture(iChannel2, p.xz*0.1).rgb;
             vec3 colYZ = texture(iChannel2, p.yz*0.1).rgb;
             vec3 colXY = texture(iChannel2, p.xy*0.1).rgb;
             
             
             col = colYZ * n1.x + colXZ * n1.y + colXY*n1.z;
             col += diff2 * vec3(1.0, 0.5, 0.5);
         }
         else if(t.z == 5.0){
             vec3 colXZ = texture(iChannel0, p.xz*0.1).rgb;
             vec3 colYZ = texture(iChannel0, p.yz*0.1).rgb;
             vec3 colXY = texture(iChannel0, p.xy*0.1).rgb;
             
             col = colYZ * n1.x + colXZ * n1.y + colXY*n1.z ;
             col += diff2 * vec3(1.0, 0.5, 0.5);
         
             
         }
         
         /*vec2 j = uv*3.0;
         j.x += 0.0;
         j.y += .1;
         float sparkle = 1./dot(j,j);
                
         col += diff2*(sparkle*sin(mod(iTime*10.0, 3.1415))*0.01) ;
    
       */
         
         
         
     }
     else
     {
         col = RayMarchCloud( ro, rd);
         //col = mix(vec3(1.0, 1.0, 1.0), col, 0.3);
         
     }
     
     if(t.x > 1000.1)col = mix(vec3(1.0, 1.0, 1.0), col, exp(-t.x*0.0005));
     //col = (1.0 - exp(-col * 6.0)) * 1.0024;
	
   
    // Output to screen
    fragColor = vec4((col),alpha);
}



