vec2 sdF1(vec3 p){

    vec2 d=vec2(0.0);
    
    vec3 pdb = p-vec3(0.0, 5.0, 0.0);
    float db = Box(pdb, vec3(9.0, 2.0, 2.0), 0.5);
    if(db > 50.0){
        d.x = db;
        d.y = 3.0;
        return d;
    }
    
    
    pdb = p-vec3(6.0, 8.0, 0.0);
    pdb *= RotZ(0.43633231305555553);
    float db2 = Box(pdb, vec3(5.0, 2.0, 3.0), 0.5);
    
        
    pdb = p - vec3(-17.0, 5.0, 0.0);
    vec4 els = opElongate(pdb, vec3(7.0, 0.0, 9.0));
    float side = els.w + sdEllipsoid(els.xyz, vec3(2.0, 3.0, 2.0));
    
    pdb = p-vec3(-28.5, 5.0, -12.0);
    pdb *= RotY(0.9424777962000002);
    float ct1 = Box(pdb, vec3(8.0, 4.0, 10.0), 0.5);
    
    pdb = p-vec3(-28.5, 5.0, 12.0);
    pdb *= RotY(2.1991148578);
    float ct2 = Box(pdb, vec3(8.0, 4.0, 10.0), 0.5);
    
    pdb = p-vec3(-30.5, 5.0, 0.0);
    float back = Box(pdb, vec3(6.0, 2.5, 3.0), 0.5);
    
    
    pdb = p - vec3(-27.0, 9.0, 0.0);
    els = opElongate(pdb, vec3(7.0, 1.0, 1.5));
    float hgh = els.w + sdEllipsoid(els.xyz, vec3(2.0, 3.0, 2.0));
    
    pdb = p-vec3(-27.0, 14.0, 0.0);
    pdb *= RotZ(-0.43633231305555553);
    float cthgh = Box(pdb, vec3(12.0, 2.0, 5.0), 0.5);
    
    hgh = smax(hgh, -cthgh, 0.3);
    
    pdb = p-vec3(-14.5, 7.0, 0.0);
    float hbc = Box(pdb, vec3(3.0, 2.0, 3.0), 0.5);
    
    pdb = p-vec3(-8.5, 6.0, 0.0);
    float hbc2 = Box(pdb, vec3(4.0, 0.7, 1.5), 0.5);
    
    pdb = p-vec3(7.0, 3.0, 0.0);
    float afr = Box(pdb, vec3(2.0, 0.05, 7.0), 0.1);
    
    /************************/
     vec3 psp = p-vec3(6.5, 4.5, -3.5 );
        
    //bend
    float k = 0.1;
    float c = cos(k*psp.z);
    float s = sin(k*psp.z);
    mat2  m = mat2(c,-s,s,c);
    vec2 qq = m*psp.xz;
    vec3 q  = vec3( qq.x,psp.y,qq.y );
    //vec3  q = vec3(m*psp.xy,psp.z);
   
    
    float sp = Box(q, vec3(0.3, 0.025, 3.5), 0.05);
    //
    psp = p-vec3(7.0, 4.25, -3.5 );
        
    //bend
    k = 0.05; 
    c = cos(k*psp.x);
    s = sin(k*psp.x);
    m = mat2(c,-s,s,c);
    q = vec3(m*psp.xy,psp.z);
   
    
    float sp2 = Box(q, vec3(0.3, 0.025, 3.5), 0.05);
    
     //
    psp = p-vec3(8.0, 4.0, -3.5 );
        
    //bend
    k = 0.2; 
    c = cos(k*psp.x);
    s = sin(k*psp.x);
    m = mat2(c,-s,s,c);
    q = vec3(m*psp.xy,psp.z);
   
    
    float sp3 = Box(q, vec3(0.3, 0.025, 3.5), 0.05);
    
    //
    /*
    psp = p-vec3(6.0, 4.5, 3.5 );
        
    //bend
    k = 0.05; 
    c = cos(k*psp.x);
    s = sin(k*psp.x);
    m = mat2(c,-s,s,c);
    q = vec3(m*psp.xy,psp.z);
   
    
    float sp4 = Box(q, vec3(0.3, 0.025, 3.5), 0.05);*/
     psp = p-vec3(6.5, 4.5, 3.5 );
        
    //bend
    k = 0.1;
    c = cos(k*psp.z);
    s = sin(k*psp.z);
    m = mat2(c,-s,s,c);
    qq = m*psp.xz;
    q  = vec3( qq.x,psp.y,qq.y );
    //vec3  q = vec3(m*psp.xy,psp.z);
   
    
    float sp4 = Box(q, vec3(0.3, 0.025, 3.5), 0.05);
    
    //
    psp = p-vec3(7.0, 4.25, 3.5 );
        
    //bend
    k = 0.05; 
    c = cos(k*psp.x);
    s = sin(k*psp.x);
    m = mat2(c,-s,s,c);
    q = vec3(m*psp.xy,psp.z);
   
    
    float sp5 = Box(q, vec3(0.3, 0.025, 3.5), 0.05);
    
     //
    psp = p-vec3(8.0, 4.0, 3.5 );
        
    //bend
    k = 0.2; 
    c = cos(k*psp.x);
    s = sin(k*psp.x);
    m = mat2(c,-s,s,c);
    q = vec3(m*psp.xy,psp.z);
   
    
    float sp6 = Box(q, vec3(0.3, 0.025, 3.5), 0.05);
    
    
    
    /*************************/
    
    pdb = p-vec3(7.0, 4.0, -7.0);
    float afr2 = Box(pdb, vec3(2.0, 1.0, 0.05), 0.1);
    
    pdb = p-vec3(7.0, 4.0, 7.0);
    float afr3 = Box(pdb, vec3(2.0, 1.0, 0.05), 0.1);
    
    afr = smin(afr, afr2, 0.1);
    afr = smin(afr, afr3, 0.1);
    
    pdb = p-vec3(-37.0, 7.0, 0.0);
    float abc = Box(pdb, vec3(2.0, 0.05, 7.0), 0.1);
    
    pdb = p-vec3(-37.0, 10.0, 7.0);
    float abc2 = Box(pdb, vec3(2.0, 4.0, 0.05), 0.1);
    abc = smin(abc, abc2, 0.1);
    
    pdb = p-vec3(-37.0, 10.0, -7.0);
    float abc3 = Box(pdb, vec3(2.0, 4.0, 0.05), 0.1);
    abc = smin(abc, abc3, 0.1);
        
    pdb = p-vec3(-37.0, 13.0, 0.0);
    pdb *= RotZ(15.0*3.1415/180.0);
    float abc4 = Box(pdb, vec3(2.0, 0.05, 7.0), 0.1);
    abc = smin(abc, abc4, 0.1);
        
    float car = smax(db, -db2, 0.3);
    car = smin(car, side, 0.3);
    car = smax(car, -ct1, 0.3);
    car = smax(car, -ct2, 0.3);  
    car = smin(car, back, 0.3);  
    car = smin(car, hgh, 0.1); 
    car = smax(car, -hbc, 0.1); 
    car = smax(car, -hbc2, 0.1); 
    ///car = smin(car, afr, 0.1); 
    
    pdb = p - vec3(5.0, 5.0, 0.0);
    float bra = sdfSegment(p, vec3(0.0, 5.0, 8.0), vec3(0.0, 5.0, -8.0), 0.1);
    
    vec3 a = vec3(0.0, 5.0, 10.0);
    vec3 b = vec3(0.0, 5.0, 6.0);
    float wh1 = sdCylinder(p, a, b, 4.0);
    
    a = vec3(0.0, 5.0, -10.0);
    b = vec3(0.0, 5.0, -6.0);
    float wh2 = sdCylinder(p, a, b, 4.0);
    
    float brb = sdfSegment(p, vec3(-30.0, 5.0, 8.0), vec3(-30.0, 5.0, -8.0), 0.1);
    
    a = vec3(-30.0, 5.0, 10.0);
    b = vec3(-30.0, 5.0, 6.0);
    float wh3 = sdCylinder(p, a, b, 4.0);
    
    a = vec3(-30.0, 5.0, -10.0);
    b = vec3(-30.0, 5.0, -6.0);
    float wh4 = sdCylinder(p, a, b, 4.0);
    
   
    float bra2 = sdfSegment(p, vec3(0.0, 3.0, 8.0), vec3(-2.0, 5.0, 2.0), 0.1);
      
    float bra3 = sdfSegment(p, vec3(0.0, 3.0, -8.0), vec3(-2.0, 5.0, -2.0), 0.1);
    
    float vol1 = sdfSegment(p, vec3(-11.2, 9.0, 1.0), vec3(-11.2, 7.0, 1.0), 0.2);
    float vol2 = sdfSegment(p, vec3(-11.2, 9.0, -1.0), vec3(-11.2, 7.0, -1.0), 0.2);
    
    pdb = p-vec3(-11.2, 8.0, -0.0);
    float vol3 = Box(pdb, vec3(0.05, 0.5, 1.0), 0.1);
    vol1 = smin(vol1, vol2, 0.1);
    vol1 = smin(vol1, vol3, 0.1);
    
    d.x = car;
    d.y = 3.0;
    
   /* d.x = db;
    d.y = 3.0;
    */
   if(afr < d.x){
        d.x = afr;
        d.y = 4.0;
    
    }
     if(sp < d.x){
        d.x = sp;
        d.y = 4.0;
    
    }
     if(sp2 < d.x){
        d.x = sp2;
        d.y = 4.0;
    
    }
    if(sp3 < d.x){
        d.x = sp3;
        d.y = 4.0;
    
    }
    if(sp4 < d.x){
        d.x = sp4;
        d.y = 4.0;
    
    }
    if(sp5 < d.x){
        d.x = sp5;
        d.y = 4.0;
    
    }
    if(sp6 < d.x){
        d.x = sp6;
        d.y = 4.0;
    
    }
    if(abc < d.x){
        d.x = abc;
        d.y = 4.0;
    
    }
    
    if(bra < d.x){
        d.x = bra;
        d.y = 5.0;
    
    }
    if(bra2 < d.x){
        d.x = bra2;
        d.y = 5.0;
    
    }
    if(bra3 < d.x){
        d.x = bra3;
        d.y = 5.0;
    
    }
    if(wh1 < d.x){
        d.x = wh1;
        d.y = 5.0;
    
    }
    if(wh2 < d.x){
        d.x = wh2;
        d.y = 5.0;
    
    }
    if(brb < d.x){
        d.x = brb;
        d.y = 5.0;
    
    }
    if(wh3 < d.x){
        d.x = wh3;
        d.y = 5.0;
    
    }
    if(wh4 < d.x){
        d.x = wh4;
        d.y = 5.0;
    
    }
    if(vol1 < d.x){
        d.x = vol1;
        d.y = 5.0;
    
    }
   
    return d;

}


/**********************************************************/

vec2 map(vec3 p, vec2 uv){


   
    vec2 d=vec2(0.0);
    float pl = dot(p, vec3(0.0, 1.0, 0.0));
    
    pl -= Noise3D(p*5.0)*0.1;
    
    vec3 pdep = vec3(0.0, 1.0, 0.0);
    pdep.x = -100.0;
    float X = mod(iTime*10.0, 200.0);
    pdep.x += X;
    if(X > 100.0)
        pdep.x = 0.0;
    vec2 df1 = sdF1(p-pdep);   
    
    float ma = sdfSegment(p, vec3(-40.0, 0.0, 20.0), vec3(-40.0, 30.0, 20.0), 0.3);
    
    vec3 pdb = p-vec3(-40.0, 25.0, 28.0);
    pdb.y += sin(p.z+iTime*10.0)*0.5;
    pdb.x += cos(p.z+iTime*10.0)*0.3;
    float flag1 = Box(pdb, vec3(0.05, 5.0, 7.0), 0.5);
    
    float ma2 = sdfSegment(p, vec3(-40.0, 0.0, -20.0), vec3(-40.0, 30.0, -20.0), 0.3);
    
    pdb = p-vec3(-40.0, 25.0, -28.0);
    pdb.y += sin(p.z+iTime*10.0)*0.5;
    pdb.x += cos(p.z+iTime*10.0)*0.3;
    float flag2 = Box(pdb, vec3(0.05, 5.0, 7.0), 0.5);
    
      
      
       d.x = pl;
       d.y = 2.0;
       
       if(df1.x < d.x){
           d.x = df1.x;
           d.y = df1.y;
       }
       if(ma < d.x){
           d.x = ma;
           d.y = 6.0;
       }
       if(flag1 < d.x){
           d.x = flag1;
           d.y = 7.0;
       }
       
       if(ma2 < d.x){
           d.x = ma2;
           d.y = 6.0;
       }
       if(flag2 < d.x){
           d.x = flag2;
           d.y = 7.0;
       }
       
    return d;
    
    
}

vec3 RayMarch(vec3 eye, vec3 viewRayDirection, vec2 uv){
    vec3 t = vec3(0.);
    for (int i = 0; i < 100; i++) {
        t.yz = map(eye + t.x * viewRayDirection, uv).xy;
        
                    
        if (abs(t.y) < 0.001)break;
                
        t.x += t.y;
        
        if (t.x >= 500.0)break;
        
    }
    if (t.x >= 500.0)t.x = -1.0;
    
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
    lightpos.xz += vec2(sin(iTime), cos(iTime));
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
    dif += occ;
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
             return vec3(dif) * vec3(1.0, 0.0, 0.0)*0.4 ;
        }
        else if(dd.z == 4.0)
        {
             return vec3(dif)* vec3(1.0, 0.3, 0.0)*0.4;;
        }
        else
            return vec3(dif)*0.2;
            
     }
     else
         return vec3(dif);

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
    vec3 ro = vec3(-25.5, 5.0, -25.0);
    vec3 lightpos = vec3(0.0, 50., -60.);
    
           
    // camera move
    vec3 dirc = vec3(-10.5, 15.0, 0.0);
    
    float an = (iMouse.x/iResolution.x)*10.0;
    float anx = (iMouse.y/iResolution.y)*3.0;
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
        
    
        
     vec3 t;
     float dO, dif, dif2, dif3, difT;
     vec3 p;
     
  
     
     t = RayMarch(ro, rd, uv);
     p = ro + rd * t.x;
     
    float alpha = 1.0;
    vec3 n1, n2, n3, np;
  
   
    
    
        
     if(t.x > 0.0){
     
         vec3 diff2 = GetLightM(p, uv, lightpos, n1, ro, rd);
                 
         //Material num
         if(t.z == 0.0){
            // col = ph * vec3(0.0, 1.0, 0.0);
         }
         else if(t.z == 1.0){
             col = diff2* vec3(1.0, 0.0, 0.0);
         }
         else if(t.z == 2.0){
             vec3 m = Line(p.xy, vec2(5.0, 0.0), vec2(-5.0, 0.0), 0.1, 0.05);;
             m+= Line3(p, vec3(0.0, 0.0, -15.0), vec3(-30.0, 0.0, -15.0), 1.0, 0.5);
             m+= Line3(p, vec3(0.0, 0.0, 15.0), vec3(-30.0, 0.0, 15.0), 1.0, 0.5);
            col = diff2*texture(iChannel1, uv).rgb + m ;
             
             
         }
         else if(t.z == -1.0){
             
         }
         else if(t.z == 3.0){
             //float spec = pow(max( dot( reflect(-normalize(lightpos-p), n1), -rd ), 0.), 8.);
             vec3 r = reflect(rd, n1);
             float spec = pow(max(0.0, r.y), 4.);
             //#ec382f
             col = diff2 * vec3(0.92, 0.21, 0.18) + spec;
         
         }
         else if(t.z == 4.0){
             vec3 r = reflect(rd, n1);
             float spec = pow(max(0.0, r.y), 4.);
             col = diff2 * vec3(0.1)+spec;
             
             
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
             vec3 r = reflect(rd, n1);
             float spec = pow(max(0.0, r.y), 4.);
             col = vec3(0.1) * diff2 + spec;
         
             
         }
         else if(t.z == 6.0){
             col = vec3(0.021) * diff2;
         
             
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
         col = RayMarchCloud( ro, rd);
         //col = vec3(0.0, 0.7, 1.0);
     }
     
      //col = mix(vec3(0.0, 1.0, 1.0), col, 0.6);//exp(-t.x*0.05));
   
    // Output to screen
    fragColor = vec4((col),alpha);
}



