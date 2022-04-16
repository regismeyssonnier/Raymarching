vec2 sdfBat(vec3 p, vec3 q){

    vec2 d = vec2(0.0);
    
    float height = 60.0;
    
    vec3 pdb = p - vec3(q.x, q.y*height, q.z);
    float bat = Box(pdb, vec3(15.0, q.y*height, 15.0), 0.2);
       
    
    d.x = bat;
    d.y = 5.0;
    
    if(q.y*height > 30.0){
        float bloc = sdCylinder(p- vec3(q.x, q.y*height, q.z), vec3(q.x, q.y*height, q.z), vec3(q.x, q.y*height+10.0, q.z), 10.0);

        if(bloc < d.x){
            d.x = bloc;
            d.y = 6.0;
        }
        
        float bloc2 = sdEllipsoid(p - vec3(q.x, q.y*height*2.2, q.z), vec3(5.0, 15.0, 5.0));
        if(bloc2 < d.x){
            d.x = bloc2;
            d.y = 6.0;
        }
        
        float pointe = sdCylinder(p- vec3(q.x, q.y*height, q.z), vec3(q.x, q.y*height, q.z), vec3(q.x, q.y*height+50.0, q.z), 0.4);
        if(pointe < d.x){
            d.x = pointe;
            d.y = 6.0;
        }
        
    }
    
    
    return d;

}

/**********************************************************/

vec2 map(vec3 p, vec2 uv){


   
    vec2 d=vec2(0.0);
    float pl = dot(p, vec3(0.0, 1.0, 0.0));
    
    vec3 c = vec3(70.0, 0.0, 70.0);
    vec2 u = p.xz;
    vec2 id = floor((u+0.5*c.xz)/c.xz-0.5*c.xz);
    
    float m = 0.0;
    
    
    vec3 q = mod(p+0.5*c,c)-0.5*c;
    
    for(int y = -2;y <=2;y++){
        for(int x =-2;x<=2;x++){
        
            m += N21(id-vec2(x, y));
        }
        
    }
    m = m/25.0;
    vec3 s = vec3(0.0, m, 0.0);
    vec2 bat = sdfBat(q, s);
      
      
       d.x = pl;
       d.y = 2.0;
       
    if(bat.x < d.x){
        d.x = bat.x;
        d.y = bat.y;
    
    }
       
      
    return d;
    
    
}

vec3 RayMarch(vec3 eye, vec3 viewRayDirection, vec2 uv){
    vec3 t = vec3(0.);
    float v, delta;
    for (int i = 0; i < 100; i++) {
        t.yz = map(eye + t.x * viewRayDirection, uv).xy;
        
                    
            if (abs(t.y) < (t.y*0.01))break;
                
        t.x += t.y;
        
        if (t.x >= 10000.0)break;
        
        
        
    }
    if (t.x >= 10000.0)t.x = -1.0;
    
    return t;


}

float calcSoftshadow( in vec3 ro, in vec3 rd, in float mint, in float tmax )
{
	float res = 1.;
    float t = mint;
    for( int i=0; i<16; i++ )
    {
		float h = map( ro + rd*t , vec2(1.0)).x;
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
        float dd = map( aopos , vec2(1.0)).x;
        occ += -(dd-hr)*sca;
        sca *= 0.95;
    }
    return clamp( 1.0 - 3.0*occ, 0.0, 1.0 );    
}



vec3 GetLightM(vec3 p, vec2 uv, vec3 lg, out vec3 n, vec3 ro, vec3 rd){
    vec3 lightpos = lg;
    //lightpos.xz += vec2(sin(iTime), cos(iTime));
    vec3 l = normalize(lightpos-p);
    
    vec2 d = map(p, uv);
    vec2 e = vec2(0.01, 0);
    
    n = d.x - vec3(
        map(p-e.xyy, uv).x,
        map(p-e.yxy, uv).x,
        map(p-e.yyx, uv).x);
    
    
    n = normalize(n);
   
    float occ = calcAO(p, n);       
    float dif = clamp(dot(n, l), .0, 1.);
    dif += occ*0.2;
    //dif *= calcSoftshadow(p, lightpos, 0.01, 100. );
 
    vec3 dd = RayMarch(p+n*.01, l, uv);
    p = ro + reflect(n, l) * dd.x;
    
     
    
         
    if((dd.x < length(lightpos-p)) && (dd.x > 0.0)){
        if(dd.z == 2.0)
        {
             
             return vec3(dif);
        }
        else if(dd.z == 3.0)
        {
             return vec3(dif) * vec3(1.0, 0.0, 0.0)*0.4;
        }
        else if(dd.z == 4.0)
        {
             return vec3(dif)* vec3(1.0, 0.3, 0.0)*0.4;;
        }
        else if(dd.z == 6.0){
             vec3 col = vec3(dif)*vec3(0.63, 0.60, 0.58)+noise(vec2(N21(p.xz)))*0.2;
             return col;

         }
        else
            return vec3(dif)*0.2;
            
     }
     else
         return vec3(dif);

}



float zero(vec2 uv){
    
    float s = smoothstep(.2, .19, length(uv-vec2(clamp(uv.x, -0.2, 0.2), clamp(uv.y, -0.2, 0.2))));
    s *= smoothstep(.05, 0.06, length(uv-vec2(clamp(uv.x, -0.2, 0.2), clamp(uv.y, -0.2, 0.2))));
    
    return s;

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
    vec3 lookat = vec3((iMouse.x/iResolution.x)*10.0, (iMouse.y/iResolution.y)*30.0, 30.0);
    vec3 ro = vec3(-25.5, 200.0, -40.0);
    vec3 lightpos = vec3(0.0, 200., -40.);
    
           
    // camera move
    vec3 dirc = vec3(-10.5, 15.0, 0.0);
   
    
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
    
    if(ro.y < 0.0){
        ro.y = 2.0;
        lightpos.y = 2.0;
    }
    
    
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
     
  
     
     t = RayMarch(ro, rd, uv);
     p = ro + rd * t.x;
     
    float alpha = 1.0;
    vec3 n1, n2, n3, np;
  
   
    
    
        
     if(t.x > 0.0){
                 
         vec3 diff2 = GetLightM(p, uv, lightpos, n1, ro, rd);
         
         if(t.z == 1.0){
             col = diff2* vec3(1.0, 0.0, 0.0);
         }
         else if(t.z == 2.0){
             vec3 m;
             vec3 c = vec3(40.0, 0.0, 30.0);
             vec3 q = mod(p+0.5*c,c)-0.5*c;
             float z = -900.0;
             for(int i = 0;i < 50;i++){
                 m+= Line3(p, vec3(35.0, 0.0, z), vec3(35.0, 0.0, z-15.0), 1.0, 0.5);
                 z += 30.0;
                 
             }
             
            col =  m ;
             
             
         }
         else if(t.z == -1.0){
             
         }
         else if(t.z == 3.0){
            
             
             col = diff2 * vec3(1.0, 0.3, 0.2);
         
         }
         else if(t.z == 4.0){
                          
             col = diff2 * vec3(0.1);
             
             
         }
         else if(t.z == 5.0){
             vec3 c = vec3(50.0, 0.0, 50.0);
             vec3 f = fract(p);
             
             vec3 id = floor(p/1.0);
             
             vec3 window = floor(p / 5.0);
             
             if(n1.y != 1.0){//#a39a95
                 bool iin = false;
                 if((mod(id.z, 2.0) == 0.0)&&(mod(id.x, 2.0) == 1.0)){
                     col = diff2*vec3(0.63, 0.60, 0.58)+noise(vec2(N21(p.xz)))*0.2;
                     iin = true;
                 }
                 if((mod(id.x, 2.0) == 0.0)&&(mod(id.z, 2.0) == 1.0)){
                     col = diff2*vec3(0.63, 0.60, 0.58)+noise(vec2(N21(p.xz)))*0.2;
                     iin = true;
                 }
                 if((mod(id.y, 2.0) == 0.0)){
                     col = diff2*vec3(0.63, 0.60, 0.58)+noise(vec2(N21(p.xz)))*0.2;
                     iin = true;
                 }
                 if(iin == false){
                     float v = (N21(id.xy)+N21(id.yz)+N21(id.xz))/3.0;
                     if(v > 0.5){
                         col = diff2*vec3(1.0, 1.0, 0.0)*(v+abs(fract(iTime)-0.5)*0.5);
                     }
                     else
                     {
                         col = diff2 * texture(iChannel0, reflect(rd, -n1)).rgb;
                     }
                 
                 }
             }
             else
             {
                  col = diff2*vec3(0.63, 0.60, 0.58)+noise(vec2(N21(p.xz)))*0.2;
             }
             
         }
         else if(t.z == 6.0){
             col = diff2*vec3(0.63, 0.60, 0.58)+noise(vec2(N21(p.xz)))*0.2;
         
             
         }
         else if(t.z == 7.0){//C
             col = vec3(0.8);
             float f = -1.0+2.0*smoothstep(-0.1, 0.1, sin(5.0*p.z)+sin(3.0*p.y));
             col += f * diff2 ;
             //col = vec3(0.5) * diff2;
         
             
         }
         else if(t.z == 8.0){
             col = vec3(0.0, 1.0, 0.0) * diff2;
         
             
         }
         else if(t.z == 9.0){
             col = vec3(1.0, 1.0, 1.0) * diff2;
         
             
         }
         
         
         
     }
     else
     {
         //col = RayMarchCloud( ro, rd);
        //col = vec3(0.0, 0.0, 0.0);
        float t = iTime*.05;
        for(float i = 0.;i <1.;i+=1./NUM_LAYERS){
            float depth = fract(i+t);
            float scale = mix(20., .5, depth);
            float fade = depth*smoothstep(1., .9, depth);
            float n = Hash21(vec2(fade*120.25, fade*12345.2));
            col += StarLayer(uv*scale+i*453.2, 1.0, iTime)*fade;
        }
     }
     
      //col = mix(vec3(0.0, 1.0, 1.0), col, 0.6);//exp(-t.x*0.05));
   
    // Output to screen
    fragColor = vec4((col),alpha);
}



