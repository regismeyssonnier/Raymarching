
/**********************************************************/

/**********************************************************/

vec2 GetDist(vec3 p, vec2 uv){


   
    vec2 d=vec2(0.0);;
    float pd = dot(p, vec3(0.0, 1.0, 0.0));
    
 
    
    
    
    vec3 pb = p-vec3(0.0, 10.0, 10.0);
    pb.y += sin(pb.y+iTime*10.0)*0.2 + cos(pb.x+iTime*10.0)*0.2 + cos(pb.z+iTime*10.0)*0.2;
    float bulle = length(pb) - 5.0;
    
    pb = p-vec3(12.0, 10.0, 10.0);
    pb.x += abs(fract(iTime/4.0)-0.5)*15.0;
    pb.y += sin(pb.y+iTime*10.0)*0.2 + cos(pb.x+iTime*10.0)*0.2 + cos(pb.z+iTime*10.0)*0.2;
    float bulle2 = length(pb) - 5.0;
    
    float sbulle = smin(bulle, bulle2, 0.9);
    
    
         
       
       d.x = sbulle;
          d.y = 3.0;
          
       
      
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
        
        if (t.x >= 500.0)break;
        
    }
    if (t.x >= 500.0)t.x = -1.0;
    
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
        float dd = GetDist( aopos , vec2(1.0)).x;
        occ += -(dd-hr)*sca;
        sca *= 0.95;
    }
    return clamp( 1.0 - 3.0*occ, 0.0, 1.0 );    
}


vec3 GetLightM(vec3 p, vec2 uv, vec3 lg, out vec3 n, vec3 ro, vec3 rd){
    vec3 lightpos = lg;
    lightpos.xz += vec2(sin(iTime), cos(iTime));
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
    dif += occ;
    
 
    return vec3(dif);

}

float isKeyPressed(float key)
{
	return texture( iChannel1, vec2(key, 0.0) ).x;
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
    vec3 lookat = vec3((iMouse.x/iResolution.x)*10.0, (iMouse.y/iResolution.y)*30.0, 30.0);//vec3(0.0, 1.5, 0.0);
    vec3 ro = vec3(0, 10.0, -15.0);
    vec3 lightpos = vec3(0.0, 50., -60.);
    
   
    float time = mod(iTime, 10.0);
    
        
    // camera move
    vec3 dirc = vec3(0.0, 10.0, 10.0);
    //dirc.z -= time*25.0;
    float an = (iMouse.x/iResolution.x)*10.0;
    float anx = (iMouse.y/iResolution.y)*3.0;
    ro -= dirc;
    lookat -= dirc;
    /*ro *= RotZ(anx);
    lookat *= RotZ(anx);
    ro *= RotX(anx);
    lookat *= RotX(anx);*/
    ro *= RotY(an);
    lookat *= RotY(an);
    ro += dirc;
    lookat += dirc;
    
    
    lightpos -= dirc;
    /*lightpos *= RotZ(anx);
    lightpos *= RotX(anx);*/
    lightpos *= RotY(an);
    lightpos += dirc;
    
    
    
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
    
    
        
     if(t.x > 0.0){
                 
         //Material num
         if(t.z == 0.0){
            // col = ph * vec3(0.0, 1.0, 0.0);
         }
         else if(t.z == 1.0){
             col = diff2* vec3(1.0, 0.0, 0.0);
         }
         else if(t.z == 2.0){
         
                                    
            col = diff2*vec3(1.0, 0.5, 0.5);;
             
             
         }
         else if(t.z == -1.0){
             
         }
         else if(t.z == 3.0){
              //reflect(rd, -n1)
             vec3 col1 = diff2 * texture(iChannel0, reflect(rd, -n1)).rgb;
             vec3 col2 = diff2 * texture(iChannel0, refract(normalize(lightpos-p), n1, 0.9)).rgb;
             float ch = abs(fract(iTime/10.0)-0.5)*2.0;
             if(isKeyPressed(82.5/256.0)>0.0){//r
                 ch = 1.0;
             }
             else if(isKeyPressed(69.5/256.0)>0.0){//e
                 ch = 0.0;
             }
                          
             col = mix(col1, col2, ch);
                          
         
         }
         else if(t.z == 4.0){
                          
             col = diff2 * vec3(1.0, 0.3, 0.0);
             
             
         }
         else if(t.z == 4.1){
             vec2 gv = fract(p.yz*2.0)-0.5;
             vec2 id = floor(p.yz);
             
             float n = N21(id);
             
             if(n <.5)gv.x *= -1.;
             float width = .1;
             float d = abs(abs(gv.x+gv.y)-0.5);
             d = length(gv-sign(gv.x+gv.y+0.001)*.5)-.5;
             float mask = smoothstep(.01, -.01, abs(d)-width);
             
             col = diff2 * vec3(1.0, 0.3, 0.0) + mask;//( (fract(p.y))+(fract(p.z)) )*0.4;
             
             
         }
         else if(t.z == 5.0){
             col = vec3(1.0) * diff2;
         
             
         }
         else if(t.z == 6.0){
             col = vec3(0.021) * diff2;
         
             
         }
         else if(t.z == 7.0){//C
             col = vec3(0.5) * diff2;
         
             
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
         col = texture(iChannel0, rd).rgb;
     }
     
     //col = mix(vec3(0.0, 1.0, 1.0), col, 0.6);//exp(-t.x*0.05));
   
    // Output to screen
    fragColor = vec4((col),alpha);
}



