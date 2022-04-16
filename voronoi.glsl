vec2 N22(vec2 p){
    vec3 a = fract(p.xyx*vec3(123.34, 234.34, 345.65));
    a+=dot(a, a+34.45);
    return fract(vec2(a.x*a.y, a.y*a.z));
}
vec2 hash22(vec2 p) { 

   
    float n = sin(dot(p, vec2(41, 289)));
    //return fract(vec2(262144, 32768)*n); 
    
    // Animated.
    p = fract(vec2(262144, 32768)*n); 
    
    return sin( p*6.2831853 + iTime )*.45 + .5; 
    
}


float Voronoi(in vec2 p){
    
	vec2 g = floor(p), o; p -= g;
	
	vec3 d = vec3(1); // 1.4, etc. "d.z" holds the distance comparison value.
    
	for(int y = -1; y <= 1; y++){
		for(int x = -1; x <= 1; x++){
            
			o = vec2(x, y);
            o += hash22(g + o) - p;
            
			d.z = dot(o, o); 
                     
            d.y = max(d.x, min(d.y, d.z));
            d.x = min(d.x, d.z); 
                       
		}
	}
	
    return max(d.y/1.2 - d.x*1., 0.)/1.2;
    
    
}

vec2 map(vec3 p){

    vec2 d = vec2(0.0);
    bool mt3 = false;
    float c = Voronoi(p.xy/3.0);
    if (c<.07) {c = smoothstep(0.7, 1., 1.-c)*.2;mt3 = true; }
    float pl =  length(1.0- p.z - c) - 2.0 ;//- (1.0- p.z - c);
   
    d.x = pl;
    if(mt3 == false)
        d.y = 2.0;
    else
        d.y = 3.0;
   
    
    return d;

}

vec3 RM(vec3 ro, vec3 rd, float _d){
    vec3 d = vec3(_d, 0.0, 0.0);
    for(int i = 0;i < 32;i++){
        d.yz = map(ro + d.x * rd).xy;
        
        if(abs(d.y) <(0.001))
            break;
            
        d.x += d.y;
        
        if(d.x > 100.0)break;
    
    }
    if(d.x > 100.0)d.x = -1.0;
    
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
    //float edge = 0.;
    // n = nr(p, edge);  
    
    vec3 l = normalize(lightpos-p);
    
    float occ = calcAO(p, n); 
    float dif = clamp(dot(n, l), 0.0, 1.0);
    dif += occ;
    
   
    
   
    vec3 sh = RM(p+n*0.01, l, 0.0);
    p = ro + rd * sh.x;
    
       
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
    vec3 ro = vec3(-2.0, 2.0, -8.0); 
    vec3 lookat = vec3(0.0, 2.0, 30.0);
    //vec3 rd = normalize(vec3(0.0, 1.0, 0.0));
    vec3 lightpos = vec3(0.0, 5.0, -1.0);
        
    vec3 dirc = vec3(0.0, 2.0, 0.0);
        
    float an = (iMouse.x/iResolution.x)*10.0;
    float anx = (iMouse.y/iResolution.y)*2.0;
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
        
    
    
    float dan = 0.0;
    vec3 d = RM(ro, rd, dan);
    
    vec4 rc;
    if(d.x > 0.0){
        p = ro + d.x * rd;
                
        
        vec3 n;
        vec3 dif = getLight(p, ro, rd, lightpos, n, uv);
    
        if(d.z == 2.0){
            vec3 r = reflect(rd, n);
            float spec = pow(max(0.0, r.y), 20.);
            
            float h = Voronoi(p.xy/3.0);
            col = h*vec3(1.0)+spec+texture(iChannel0, reflect(rd, n)).rgb;
           
        
        }
        else if(d.z == 3.0){
            col = dif * vec3(1.0, 0.0, 0.0)*texture(iChannel0, reflect(rd, -n)).rgb;
            
        }
      
        
          

    }                 
    else
    {
        col = vec3(0.952, 0.243, 0.278);
    
    }
    
    

    // Output to screen
    fragColor = vec4(col,1.0);
}