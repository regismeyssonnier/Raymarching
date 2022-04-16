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

vec2 map2(vec3 p, out vec4 resColor){

    vec2 d = vec2(0.0);
    
    //float pl = p.y;
   
   vec3 w = p;
   
    float m = dot(w,w);

    vec4 trap = vec4(abs(w),m);
	float dz = 1.0;
    
    
	for( int i=0; i<4; i++ )
    {

        dz = 8.0*pow(sqrt(m),7.0)*dz + 1.0;
		//dz = 8.0*pow(m,3.5)*dz + 1.0;
        
        float r = length(w);
        float b = 8.0*acos( w.y/r);
        float a = 8.0*atan( w.x, w.z );
        w = p + pow(r,8.0) * vec3( sin(b)*sin(a), cos(b), sin(b)*cos(a) );
       
  
        trap = min( trap, vec4(abs(w),m) );

        m = dot(w,w);
		if( m > 256.0 )
            break;
    }

    resColor = vec4(m,trap.yzw);

    d.x = 0.25*log(m)*sqrt(m)/dz;
    d.y = 1.0;
    
    
    return d;

}

vec2 map(vec3 p){

    vec2 d = vec2(0.0);
    
    
   
   vec3 w = p;
  
    float m = dot(w,w);

    vec4 trap = vec4(abs(w),m);
	float dz = 1.0;
    
    
	for( int i=0; i<4; i++ )
    {

        dz = 8.0*pow(sqrt(m),7.0)*dz + 1.0;
		//dz = 8.0*pow(m,3.5)*dz + 1.0;
        
        float r = length(w);
        float b = 8.0*acos( w.y/r);
        float a = 8.0*atan( w.x, w.z );
        w = p + pow(r,8.0) * vec3( sin(b)*sin(a), cos(b), sin(b)*cos(a) );

        
        trap = min( trap, vec4(abs(w),m) );

        m = dot(w,w);
		if( m > 256.0 )
            break;
    }

    //resColor = vec4(m,trap.yzw);

    d.x = 0.25*log(m)*sqrt(m)/dz;
    d.y = 1.0;
   
    
    return d;

}

vec3 RM(vec3 ro, vec3 rd, float _d){
    vec3 d = vec3(_d, 0.0, 0.0);
    for(int i = 0;i < 200;i++){
        d.yz = map(ro + d.x * rd).xy;
        
        if((d.y) <(0.005))
            break;
            
        d.x += d.y;
        
        if(d.x > 1000.0)break;
    
    }
    if(d.x > 1000.0)d.x = -1.0;
    
    return d;
    
}

vec3 calcNormal( in vec3 pos, in float t, in float px )
{
    vec4 tmp;
    vec2 e = vec2(1.0,-1.0)*0.5773*0.25*px;
    return normalize( e.xyy*map( pos + e.xyy ).x + 
					  e.yyx*map( pos + e.yyx ).x + 
					  e.yxy*map( pos + e.yxy ).x + 
					  e.xxx*map( pos + e.xxx ).x );
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


vec3 getLight(vec3 p, vec3 ro, vec3 rd, vec3 lightpos, out vec3 n){
    vec2 e = vec2(0.01, 0.0);
    vec2 nd = map(p);

    n = nd.x - vec3(map(p - e.xyy).x ,
                        map(p- e.yxy).x ,
                        map(p- e.yyx).x );

    n = normalize(n);
    
    const float fle = 1.5;

    float px = 2.0/(iResolution.y*fle);
   // n = calcNormal(p, 0.0, px);

    //vec3 lightpos = vec3(0.0, 20.0, -20.0);
    vec3 l = normalize(lightpos-p);
    
    float occ = calcAO(p, n); 
    float dif = clamp(dot(n, l), 0.0, 1.0);
    dif += occ;

    /*vec3 sh = RM(p+n*0.01, l);
    p = ro + rd * sh.x;

    if((sh.x > 0.0) && (length(lightpos-p) > sh.x)){

       dif = dif * 0.2;
    }*/
    
    return vec3(dif);

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
    vec3 ro = vec3(-2.0, -0.5, -7.0); 
    vec3 lookat = vec3((iMouse.x/iResolution.x)*10.0, (iMouse.y/iResolution.y)*30.0, 30.0);
    //vec3 rd = normalize(vec3(0.0, 1.0, 0.0));
    vec3 lightpos = vec3(0.0, 50.0, -5.0);
        
    vec3 dirc = vec3(0.0, 2.0, 0.0);
        
    float an = (iMouse.x/iResolution.x)*10.0;
    float anx = (iMouse.y/iResolution.y)*10.0;
    ro -= dirc;
    lookat -= dirc;
    ro *= RotZ(anx);
    lookat *= RotZ(anx);
    ro *= RotX(anx);
    lookat *= RotX(anx);
    ro *= RotY(an);
    lookat *= RotY(an);
    ro += dirc;
    lookat += dirc;
    
    lightpos -= dirc;
    lightpos *= RotZ(anx);
    lightpos *= RotX(anx);
    lightpos *= RotY(an);
    lightpos += dirc;
    
    
    float zoom = abs(fract(iTime/10.0)-0.5)*10.0+2.0;    
    vec3 f = normalize(lookat-ro),
        r = normalize(cross(vec3(0,1,0), f)),
        u = cross(f, r),
        c = ro + f * zoom,
        i = c + uv.x * r + uv.y * u,
        rd = normalize(i-ro);
    
    float dan ;
    if(anx >2.80)dan=4.0;
    else dan = 5.0;
    vec3 d = RM(ro, rd, dan);
    
    vec4 rc;
    if(d.x > 0.0){
        p = ro + d.x * rd;
        
    vec2  sp = (2.0*p.xy-iResolution.xy) / iResolution.y;
    vec2 m2 = map2(p, rc);
    vec3 n;
    vec3 dif = getLight(p, ro, rd, lightpos, n);
    
    const vec3 light1 = vec3(  0.577, 0.577, -0.577 );
     vec3 hal = normalize( (light1)-rd);
     float spe1 = pow( clamp(dot(n,hal),0.0,1.0), 32.0 )*dif.x*(0.04+0.96*pow(clamp(1.0-dot(hal,lightpos),0.0,1.0),5.0));
       
      col = vec3(0.01);
       //col = dif * vec3(1.0, 0.5, 1.0);//rc.rgb;
		col = mix( col, vec3(0.10,0.20,0.50), clamp(rc.y,0.0,1.0) );
	 	col = mix( col, vec3(0.02,0.10,0.50), clamp(rc.z*rc.z,0.0,1.0) );
        col = mix( col, vec3(0.50,0.10,0.02), clamp(pow(rc.w,6.0),0.0,1.0) );
       // col *= 0.5;
        vec3 lin = vec3(0.0); 
		     lin += 7.0*vec3(1.50,1.10,0.70)*dif.x;
        col *= lin;
      //col = dif * rc.rgb;//vec3(1.0, 0.5, 0.5);
      col += spe1*15.0;
       col = sqrt(col);
       col *= 1.0 - 0.05*length(sp);

    }                 
                        
    
    

    // Output to screen
    fragColor = vec4(col,1.0);
}