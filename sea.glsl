vec4 tpos = vec4(0.0, 0.0, 0.0, 0.0);
float angle_tpos;
bool done = false;

mat2 Rot(float a){
    float s = sin(a);
    float c = cos(a);
    return mat2(c, -s, s, c);
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



vec2 GetDist(vec3 p, vec2 uv){


    
    
    vec4 s = vec4(1.0,0.0, 2.0,1.0);
    vec4 s2 = vec4(2.0 ,0.0, 6.0, 0.2);
    vec4 s3 = vec4(3.0,7.0, 5.0,0.5);
    
    
    
    vec2 d=vec2(0.0);;
    float pd = dot(p, vec3(0.0, 1.0, 0.0));
    
    
    pd = pd - noi(p*0.91+iTime)*0.5;// -sin(p.x+iTime*2.0)*0.2-cos(p.z+iTime*2.0)*0.1+sin(p.y*3.0+iTime)*0.1;//- length(texture(iChannel0, (uv*.03) ).rgb)*.8;
    
    vec3 sp = (p)-s.xyz;
     sp.y += sin(iTime)*0.5;
    float sd = length((sp))-s.w;
    
    
    
    //sd -= noi(sp) - noi(sp*15.)*0.1;
    
    
    vec3 sp2 = (p)-s2.xyz;
    sp2.y += sin(iTime)*0.1;
    float sd2 = length((sp2))-s2.w;
    
    
    vec3 bx = vec3(7.0, 0.0, 4.);
    vec3 bp = p-bx;
    vec3 sc = vec3(1.0, 1.0, 0.5);
    
    bp.y += sin(iTime*0.5);
    //bp -= noi(bp) - noi(bp*15.)*0.1;
           
    float dbv = Box(bp, sc, 0.5);
        
    vec3 b1 = vec3(-7.0, 1.0, 15.0);
    vec3 b2 = vec3(-7.0, 1.0, -2);
    
    float balle = sdCylinder(p, b1, b2, 4.0);
    balle -= noi(p) - noi(p*5.)*0.1;
    
    b1 = vec3(20.0, 1.0, 15.0);
    b2 = vec3(20.0, 1.0, -2);
    
    float balle2 = sdCylinder(p, b1, b2, 4.0);
    balle2 -= noi(p/2.5) ;//- noi(p*2.)*0.1;
    
    vec3 sp3 = (p)-s3.xyz;
    //sp3.y += sin(iTime)*0.1;
    float sd3 = length((sp3))-s3.w;
    
      
      if(sd < sd2){
          d.x = sd;
          d.y = -1.0;
       }
       else
       {
           d.x = sd2;
          d.y = 0.0;
       }
       if(dbv < d.x){
          d.x = dbv;
          d.y = 1.0;
       }
            
       if(pd < d.x){
          d.x = pd;
          d.y = 2.0;
       }
       
       if(balle < d.x){
          d.x = balle;
          d.y = 3.0;
       }
       
       if(balle2 < d.x){
          d.x = balle2;
          d.y = 5.0;
       }
       
       if(sd3 < d.x){
          d.x = sd3;
          d.y = 4.0;
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
        dd = GetDist(eye + t.x * viewRayDirection, uv);
        t.y = dd.x;
        
        
            
        if (t.y < 0.01) {
           t.z = dd.y;
           return t;
        }
        t.x += t.y;
        
        

        if (t.x >= 70.0) {
            t.z = dd.y;
            return t;
        }
    }
    
    
    
    t.z = dd.y;
    
    
    return t;


}




float GetLight(vec3 p, vec2 uv, vec3 lg, out vec3 n){
    vec3 lightpos = lg;
    //lightpos.xz += vec2(sin(iTime), cos(iTime));
    vec3 l = normalize(lightpos-p);
    
    vec2 d = GetDist(p, uv);
    vec2 e = vec2(0.01, 0);
    
    
    float dd = d.x;
    n = d.x - vec3(
        GetDist(p-e.xyy, uv).x,
        GetDist(p-e.yxy, uv).x,
        GetDist(p-e.yyx, uv).x);
    
    n = normalize(n);
    
    float dif = dot(refract(l,n, 0.9), l);
    //float dif = clamp(dot(n, l), 0., 1.);
    //vec2 dd = RayMarch2(p+n*.01, l, uv);
   //if(dd.x < length(lightpos-p))dif *= 0.1;
    return dif ;

}

vec3 GetLightM(vec3 p, vec2 uv, vec3 lg, out vec3 n, vec3 ro, vec3 rd){
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
  
    vec3 dd = RayMarch2(p+n*.01, l, uv);
    p = ro + reflect(n, l) * dd.x;
   
    //lightpos-p
    if(dd.x < length(vec3(1.0, 1.0, 1.0)-p/2.0)){
        if(dd.z == -1.0){
            return vec3(dif)*vec3(1.0, 0.1, 1.0);
        }
        else if(dd.z == 1.0){
            return vec3(dif)*vec3(1.0, 0.0, 0.0);
        }
        if(dd.z == 0.0){
            return vec3(dif)*vec3(0.0, 1.0, 0.0);
        }
        else if(dd.z == 3.0)
        {
             vec3 colXZ = texture(iChannel1, p.xz*0.1).rgb;
             vec3 colYZ = texture(iChannel1, p.yz*0.1).rgb;
             vec3 colXY = texture(iChannel1, p.xy*0.1).rgb;
             
             vec3 col = (colYZ * n.x + colXZ * n.y + colXY*n.z)*vec3(dif);
             col += vec3(0.0, 0.8, 0.8);
            return col;
             
        }
        else if(dd.z == 5.0){
        
            vec3 colXZ = texture(iChannel1, p.xz*0.1).rgb;
             vec3 colYZ = texture(iChannel1, p.yz*0.1).rgb;
             vec3 colXY = texture(iChannel1, p.xy*0.1).rgb;
             
             vec3 col = (colYZ * n.x + colXZ * n.y + colXY*n.z)*vec3(dif);
              col -= vec3(0.0, 0.4, 0.4);
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
    
    vec3 ro = vec3(0, 3.0, -4.);
    //ro.xz *= Rot(iTime);
   //ro.y = sin(iTime)*0.5+0.5;
    vec3 lookat = vec3((iMouse.x/iResolution.x)*10.0, (iMouse.y/iResolution.y)*30.0, 10.0);//vec3(0.0, 1.5, 0.0);
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
     
    vec3 n1, n2, n3, nl1;
    float alpha = 1.0;
    vec3 diff = GetLightM(p, uv, vec3(-10.0, 5, -6), n1, ro, rd);
    vec3 diff2 = GetLightM(p, uv, vec3(20.0, 5., 3.0), n2, ro, rd);
    dif = GetLight(p, uv, vec3(-10.0, 5, -6), nl1);
    //dif2 = GetLight(p, uv, vec3(5.0, 5, 0), n2);
    dif3 = GetLight(p, uv, vec3(10.0, 1.0, 3.0), n3);
    //phongIllumination(vec3 k_a, vec3 k_d, vec3 k_s, float alpha, vec3 p, vec3 eye)
    vec3 ph = phongIllumination(vec3(1.0), vec3(dif), vec3(1.0, 1.0, 0.0), 50.0, p, ro, uv);
    
     if(t.y >= 0.1){
         //col = vec3(0.0, 0.8, 1.0);
         col = RayMarchCloud(ro, rd);
         
     }
     else{
     
                  
         if(t.z == 0.0){
             col = diff * vec3(0.0, 1.0, 1.0);//vec3(dif)* vec3(0.0, 1.0, 0.0)+ vec3(dif2)* vec3(0.0, 1.0, 0.0);
         }
         else if(t.z == 1.0){
             col = diff* vec3(1.0, 0.0, 0.0);//vec3(dif)* vec3(1.0, 0.0, 0.0)+ vec3(dif2)* vec3(1.0, 0.0, 0.0);
         }
         else if(t.z == 2.0){
             alpha = 1.0;
             // col += mix(vec3(0.0, 0.8, 1.0), vec3(dif)* vec3(0.0, 1.0, 1.0), pow(smoothstep(0., 1.0, dif),2.0) );
            col = mix(diff * vec3(1.0, 0.5, 0.5), diff2 * vec3(1.0, 0.5, 0.5), exp(vec3(0.4)) );
             //col += mix(vec3(0.0, 0.1, 0.1), vec3(dif)* vec3(0.0, 1.0, 1.0), vec3(0.8) );
             //pow(smoothstep(0., 1.0, dif),2.0)
            col += normalize(vec3(0.1, 0.5, 0.5));
         }
         else if(t.z == -1.0){
             col = vec3(dif)*vec3(1.0, 0.1, 1.0);//vec3(dif)* vec3(1.0, 0.1, 1.0)+ vec3(dif2)* vec3(1.0, 0.1, 1.0);
         }
         else if(t.z == 3.0){
             vec3 colXZ = texture(iChannel1, p.xz*0.1).rgb;
             vec3 colYZ = texture(iChannel1, p.yz*0.1).rgb;
             vec3 colXY = texture(iChannel1, p.xy*0.1).rgb;
             
             col = (colYZ * nl1.x + colXZ * nl1.y + colXY*nl1.z)*vec3(dif);
             
             //col = vec3(dif)* vec3(1.0, 0.5, 0.5) + vec3(dif2)* vec3(1.0, 0.5, 0.5);
         }
         else if(t.z == 4.0){
             
             col = vec3(dif3)* vec3(1.0, 0.5, 0.5) + vec3(dif3)*texture(iChannel1, uv+iTime*0.01).rgb;
         }
         else if(t.z == 5.0){
             vec3 colXZ = texture(iChannel1, p.xz*0.1).rgb;
             vec3 colYZ = texture(iChannel1, p.yz*0.1).rgb;
             vec3 colXY = texture(iChannel1, p.xy*0.1).rgb;
             
             col = (colYZ * n3.x + colXZ * n3.y + colXY*n3.z)*vec3(dif3);
             
             //col = vec3(dif)* vec3(1.0, 0.5, 0.5) + vec3(dif2)* vec3(1.0, 0.5, 0.5);
         }
         
         
         vec2 j = uv*3.0;
         j.x += 0.0;
         j.y += .1;
         float sparkle = 1./dot(j,j);
                
         col += vec3(dif)*(sparkle*sin(mod(iTime*10.0, 3.1415))*0.01) ;
    
         
         
         
     }
     
   
    // Output to screen
    fragColor = vec4((col),alpha);
}



