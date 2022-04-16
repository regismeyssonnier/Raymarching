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

mat3 quat(vec3 p, float theta){

    vec3 z = normalize(p);
    
    float t1 =  cos(theta);
    float t2 =  1. - t1;
    float t3 =  z.x*z.x;
    float t6 =  t2*z.x;
    float t7 =  t6*z.y;
    float t8 =  sin(theta);
    float t9 =  t8*z.z;
    float t11 = t6*z.z;
    float t12 = t8*z.y;
    float t15 = z.y*z.y;
    float t19 = t2*z.y*z.z;
    float t20 = t8*z.x;
    float t24 = z.z*z.z;
    return mat3( t1 + t2*t3, t7 - t9, t11 + t12, t7 + t9, t1 + t2*t15, t19 - t20, t11 - t12, t19 + t20, t1 + t2*t24);

}

mat3 quaternion(vec3 v, float an){

    vec3 z = normalize(v);
    
    float a = cos(an/2.0);
    float s = sin(an/2.0);
    float b = s*z.x;
    float c = s*z.y;
    float d = s*z.z;
    float a2 = a*a;
    float b2 = b*b;
    float c2 = c*c;
    float d2 = d*d;
    float _2ab = 2.*a*b;
    float _2ac = 2.*a*c;
    float _2ad = 2.*a*d;
    float _2bc = 2.*b*c;
    float _2bd = 2.*b*d;
    float _2cd = 2.*c*d;
    
    return mat3(a2+b2-c2-d2  , _2bc-_2ad, _2ac+_2bd, 
                _2ad+_2bc, a2-b2+c2-d2  , _2cd-_2ab,
                _2bd-_2ac, _2ab+_2cd, a2-b2-c2+d2);


}

float Box(vec3 p, vec3 sc, float r){
    return length(max(abs(p)-sc, 0.))-r;
}

float Box2d(vec2 p, vec2 sc, float r){
    return length(max(abs(p)-sc, 0.))-r;
}

float N21(vec2 p){
    p = fract(p*vec2(233.34, 851.73));
    p += dot(p, p+23.45);
    return fract(p.x*p.y);

}

vec2 map(vec3 p){

    vec2 d = vec2(0.0);
    
    float r1 = 3.0, r2 = 0.5;
    vec3 pp = p - vec3(0.0, 4.0, 0.0);
    float x = length(pp.xz) - r1;
    float y = pp.y;
    vec2 cp = vec2(x, y);
    float a = atan(pp.x, pp.z);
    cp *= Rot(a*2.5+iTime);
    cp.y = abs(cp.y)-0.7;
    //float tp = length(cp) - r2;
    
    float tp = Box2d(cp, vec2(.1, .1*sin(a)*.5+0.5), 0.05)*.5;
   
    d.x = tp;
    d.y = 1.0;
    
    
    
  
    
    return d;

}

vec3 RM(vec3 ro, vec3 rd, float _d){
    vec3 d = vec3(_d, 0.0, 0.0);
    for(int i = 0;i < 100;i++){
        d.yz = map(ro + d.x * rd).xy;
        
        if(abs(d.y) <(0.001))
            break;
            
        d.x += d.y;
        
        if(d.x > 30.0)break;
    
    }
    if(d.x > 30.0)d.x = -1.0;
    
    return d;
    
}


float calcAO( in vec3 pos, in vec3 nor )
{
	float occ = 0.0;
    float sca = 1.0;
    for( int i=0; i<5; i++ )
    {
        float hr = 0.01 + 0.12*float(i)/4.0;
        vec3 aopos =  nor * hr + pos;
        float dd = map( aopos ).x;
        occ += -(dd-hr)*sca;
        sca *= 0.95;
    }
    return clamp( 1.0 - 3.0*occ, 0.0, 1.0 );    
}




vec3 getLight(vec3 p, vec3 ro, vec3 rd, vec3 lightpos, out vec3 n, vec2 uv){
    vec2 e = vec2(0.01, 0.0);
    vec2 nd = map(p);

    n = nd.x - vec3(map(p - e.xyy).x ,
                        map(p- e.yxy).x ,
                        map(p- e.yyx).x );

    n = normalize(n);
    
        
    
    vec3 l = normalize(lightpos-p);
    
    float occ = calcAO(p, n); 
    float dif = clamp(dot(n, l), 0.0, 1.0);
    dif += occ;
    
   
    
    
    vec3 sh = RM(p+n*0.01, l, 0.0);
    p = ro + rd * sh.x;
    
   
    
    float spec = pow(max( dot( reflect(-l, n), -rd ), 0.), 8.);
    
    vec3 col = vec3(0.0);
       
    if((sh.x > 0.0) && (sh.x < length(lightpos-p))){

        if(sh.z == 1.0){

            col = vec3(dif);
            return col;

            
        }
        else
             return vec3(dif);
     
            
    }
    else
         return vec3(dif);
   
   
}

vec3 Bg(vec3 rd){

    float k = rd.y*.5+.5;
    
    
    vec3 col = mix(vec3(.2, .1, .1), vec3(.2, .5, 1.), k);
    return col;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Normalized pixel coordinates (from 0 to 1)
    vec2 uv = fragCoord/iResolution.xy;
    uv -= 0.5;
    uv /= vec2(iResolution.y/iResolution.x, 1.0);
    
    // Time varying pixel color
    vec3 col = vec3(0.0);
    
       
    vec3 p;
    vec3 ro = vec3(-2.0, 3.0, -7.0); 
    vec3 lookat = vec3((iMouse.x/iResolution.x)*10.0, (iMouse.y/iResolution.y)*30.0, 30.0);
    vec3 lightpos = vec3(0.0, 5.0, -1.0);
        
    vec3 dirc = vec3(0.0, 2.0, 0.0);
        
    float an = (iMouse.x/iResolution.x)*10.0;
    float anx = (iMouse.y/iResolution.y)*10.0;
    ro -= dirc;
    lookat -= dirc;
    ro *= quaternion(vec3(0.0, 0.0, 1.0), anx);
    lookat *= quaternion(vec3(0.0, 0.0, 1.0), anx);
    ro *= quaternion(vec3(1.0, 0.0, 0.0), anx);
    lookat *= quaternion(vec3(1.0, 0.0, 0.0), anx);
    ro *= quaternion(vec3(0.0, 1.0, 0.0), an);
    lookat *= quaternion(vec3(0.0, 1.0, 0.0), an);
    
    ro += dirc;
    lookat += dirc;
      
    vec3 n = cross(vec3(0.0, 2.0, 0.0), vec3(iMouse.xy, 0.0));
    float a = acos(dot(vec3(0.0, 2.0, 0.0), vec3(iMouse.xy, 0.0)));
    
    
    
    lightpos -= dirc;
    lightpos *= quaternion(vec3(0.0, 0.0, 1.0), anx);
    lightpos *= quaternion(vec3(1.0, 0.0, 0.0), anx);
    lightpos *= quaternion(vec3(0.0, 1.0, 0.0), an);
    lightpos += dirc;
    
       
    
    float zoom = 1.0;    
    vec3 f = normalize(lookat-ro),
        r = normalize(cross(vec3(0,1,0), f)),
        u = cross(f, r),
        c = ro + f * zoom,
        i = c + uv.x * r + uv.y * u,
        rd = normalize(i-ro);
        
    col += Bg(rd);
    
    float dan = 0.0;
    vec3 d = RM(ro, rd, dan);
    
    vec4 rc;
    if(d.x > 0.0){
        p = ro + d.x * rd;
        
        vec3 n;
        vec3 dif = getLight(p, ro, rd, lightpos, n, uv);
    
        if(d.z == 1.0){
            
            vec3 r = reflect(rd, n);
            float spec = pow(max(0.0, r.y), 20.);
            col += mix(Bg(r), dif, .5)+spec ;
                     
           
        
        }
        else if(d.z == 2.0){
            col += dif * vec3(1.0, 1.0, 0.1);
        
        
        }

    }                 
                        
    
    

    // Output to screen
    fragColor = vec4(col,1.0);
}