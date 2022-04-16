
vec2 sdfShip(vec3 p){

    vec2 d = vec2(0.0);

    //left
    vec3 q = p - vec3(-4.0, 2.0, 0.0);
    q *= RotY(1.570796327);
    float tp = sdHexPrism(q-vec3(0.0, 2.0, 0.0), vec2(2.0, 0.1))-0.1;
    //d.x = tp;
    //d.y = 1.0;
    
    q = p-vec3(-3.4, 4.0, 0.0);
    q *= RotY(3.1415);
    float tp2 = sdBox3(q, vec3(0.2, 0.2, 0.2), vec3(0.5, 0.5, 0.5), 0.05, 0.6);
    
    q = p-vec3(-2.5, 4.0, 0.0);
    q *= RotZ(1.570796327);
    float tpc = sdCyl(q, 0.1, 1.0);
    
    q = p-vec3(-1.5, 4.0, 0.0);
    q *= RotZ(1.570796327);
    float tp3 = sdCyl2(q, 0.1, 0.3, 0.2);
    
    
    //right
    q = p - vec3(4.0, 2.0, 0.0);
    q *= RotY(1.570796327);
    float tpr = sdHexPrism(q-vec3(0.0, 2.0, 0.0), vec2(2.0, 0.1))-0.1;
    //d.x = tp;
    //d.y = 1.0;
    
    q = p-vec3(3.4, 4.0, 0.0);
    //q *= RotY(3.1415);
    float tpr2 = sdBox3(q, vec3(0.2, 0.2, 0.2), vec3(0.5, 0.5, 0.5), 0.05, 0.6);
    
    q = p-vec3(2.5, 4.0, 0.0);
    q *= RotZ(1.570796327);
    float tprc = sdCyl(q, 0.1, 1.0);
    
    q = p-vec3(1.5, 4.0, 0.0);
    q *= RotZ(-1.570796327);
    float tpr3 = sdCyl2(q, 0.1, 0.3, 0.2);
    
    
    //center
    q = p-vec3(0.0, 4.0, 0.0);
    float cen = length(q)-1.45;
    
    q = p - vec3(0.0, 4.0, -2.5);
    float bcen = Box(q, vec3(1.5), 0.1);
    
    q = p-vec3(0.0, 4.0, 0.0);
    float bcen2 = length(q)-1.25;
    
    
    
    
    float TP = smin(tp, tp2, 0.1);
    TP = smin(TP, tpc, 0.1);
    TP = smin(TP, tp3, 0.1);
    //right
    TP = smin(TP, tpr, 0.1);
    TP = smin(TP, tpr2, 0.1);
    TP = smin(TP, tprc, 0.1);
    TP = smin(TP, tpr3, 0.1);
    ///center
    TP = smin(TP, cen, 0.1);
    TP = smax(TP, -bcen, 0.1);
    //TP = smin(TP, bcen2, 0.1);
    
    d.x = TP;
    d.y = 1.0;
    
    if(bcen2 <= d.x){
        d.x = bcen2;
        d.y = 3.0;
    
    }
    
    if(tp <= d.x){
        d.x = tp;
        d.y = 4.0;
    
    }
    
    if(tpr <= d.x){
        d.x = tpr;
        d.y = 4.0;
    
    }

    return d;

}

vec2 map(vec3 p){

    vec2 d = vec2(0.0);
    
    float terr = TerrainRM(p.xz);
    float pl = p.y-terr;
   
    d.x = pl;
    d.y = 2.0;
    
   
 
    vec3 q = p - vec3(0.0, 2.0+terr, 0.0);
    q.z += iTime*5.0;
    vec2 tp = sdfShip(q*4.0);
    tp.x /= 4.0;
    if(tp.x < d.x){
        d.x = tp.x;
        d.y = tp.y;
    
    }
    
    q = p - vec3(-3.5, 2.0+terr, 2.0);
    q.z += iTime*5.0;
    vec2 tp2 = sdfShip(q*4.0);
    tp2.x /= 4.0;
    if(tp2.x < d.x){
        d.x = tp2.x;
        d.y = tp2.y;
    
    }
    
    q = p - vec3(3.5, 2.0+terr, 2.0);
    q.z += iTime*5.0;
    vec2 tp3 = sdfShip(q*4.0);
    tp3.x /= 4.0;
    if(tp3.x < d.x){
        d.x = tp3.x;
        d.y = tp3.y;
    
    }
    
    
  
    
    return d;

}

vec3 RM(vec3 ro, vec3 rd, float _d){
    vec3 d = vec3(_d, 0.0, 0.0);
    for(int i = 0;i < 100;i++){
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

vec3 lightDir = normalize( vec3(0.5,0.6,0.) );
const mat2 m2 = mat2( 0.60, -0.80, 0.80, 0.60 );
vec3 Cloud(vec3 bgCol,vec3 ro,vec3 rd,vec3 cloudCol,float spd)
{
    vec3 col = bgCol;
    float t = iTime * 0.15* spd;
    //float r = N21(vec2(t)) * 20.5;
    vec2 sc = ro.xz + rd.xz*2.5/rd.y;//*(12.0-ro.y)/rd.y;
    vec2 p = 0.2*sc;
    float f = 0.0;
  	float s = 0.5;
  	float sum =0.;
  	for(int i=0;i<5;i++){
    	p += t;t *=1.5;
    	f += s* Noise(p) /*texture (iChannel0, p/256.0).x/*textureLod( iChannel0, p/256.0, 0.0).x*/; p = m2*p*2.02;
    	sum+= s;s*=0.6;
  	}
    float val = f/sum; 
    col = mix( col, cloudCol, smoothstep(0.5,0.8,val) );
    return col;
}
vec3 RayMarchCloud(vec3 ro,vec3 rd){
    vec3 col = vec3(0.0,0.0,0.0);  
    /*float sundot = clamp(dot(rd,lightDir),0.0,1.0);
    */
     // sky      
    //col = vec3(0.2,0.5,0.85)*1.1 - rd.y*rd.y*0.5;
   // col = mix( col, 0.85*vec3(0.7,0.75,0.85), pow( 1.0-max(rd.y,0.0), 4.0 ) );
    // sun
    /*
    col += 0.25*vec3(1.0,0.7,0.4)*pow( sundot,5.0 );
    col += 0.25*vec3(1.0,0.8,0.6)*pow( sundot,64.0 );
    col += 0.4*vec3(1.0,0.8,0.6)*pow( sundot,512.0 );*/
     col = vec3(0.0, 0.7, 1.0);
    // clouds
    col = Cloud(col,ro,rd,vec3(1.0,0.95,1.0),1.);
            // .
    col = mix( col, 1.5*vec3(0.0,0.5,1.0), pow( 1.0-max(rd.y,0.0), 16.0 ) );
    return col;
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
    
    
    lightpos.z -= iTime*5.0;
    dirc.z -= iTime*5.0;
    ro.z -= iTime*5.0;
    lookat.z -= iTime*5.0;
    
    
    float terr = TerrainRM(ro.xz);
        lightpos.y += terr;
        dirc.y += terr;
        ro.y += terr;
        
        lookat.y += terr;
        
    if(terr <= 0.0){
        lookat.y = ro.y +10.0 ;
    
    }
    else
    {
        lookat.y = ro.y -10.0 ;
    
    }
    
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
    
        if(d.z == 1.0){
            
            vec3 r = reflect(rd, n);
            float spec = pow(max(0.0, r.y), 20.);
            col += mix(Bg(r), dif, .5)+ spec ;
           
           
        
        }
        else if(d.z == 2.0){
            vec3 r = reflect(rd, n);
            float spec = pow(max(0.0, r.y), 5.);
            col = mix(vec3(1.0, 0.5, 0.2)*dif, vec3(Noise(p.xz)*0.1), 0.5)+spec;
            
        
        
        }
        else if(d.z == 3.0){
            col = dif * vec3(1.0, 0.3, 0.2)*texture(iChannel0, reflect(rd, -n)).rgb;
            
        }
        else if(d.z == 4.0){
            
            col = dif * vec3(0.2, 0.2, 0.2) ;
            
        }
        
          

    }                 
    else
    {
        col = RayMarchCloud(ro, rd);
    
    }
    
    

    // Output to screen
    fragColor = vec4(col,1.0);
}