#define MAX_STEPS 500
#define MAX_DIST 70.
#define SURF_DIST .01
#define EPSILON 0.01

int what = 0;

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

float GetDistPlan(vec3 p, vec2 uv){
    vec4 s = vec4(0,3.5, 7, 0.1);
    
    float planeDist = dot(p, vec3(0.0, 1.0, 0.0))-sin(p.x+iTime*2.0)*0.2-cos(p.z+iTime)*0.5+sin(p.y*2.0+iTime*2.0)*0.3;
           
    vec3 nc = texture(iChannel0, normalize(uv*-iTime*.03) ).rgb;
    
    float d = planeDist - length(noise(uv))*1.5 - length(nc.rgb)*0.1;
    
    return d;

}

float GetDist(vec3 p, vec2 uv){
    vec4 s = vec4(0,3.5, 7,1.5);
    vec4 s2 = vec4(4.0,3.5, 10.0,1.5);
    
    float planeDist = dot(p, vec3(0.0, 1.0, 0.0))-sin(p.x+iTime)*0.2+cos(p.z+iTime)*0.05;
    float sphereDist = length(p-s.xyz)-s.w;
    float sphereDist2 = length(p-s2.xyz)-s2.w;
    //float sphereDist = CubeApproxSDF(p, 2.0);
    
    vec3 nc = texture(iChannel0, uv*2.-iTime*.03).rgb;
    
    //float d = min(sphereDist-length(nc), planeDist-length(nc)*0.1);
    
     float d = min(sphereDist, sphereDist2);//-length(nc);
    return d;

}

float RayMarch2(vec3 eye, vec3 viewRayDirection,vec2 uv){

    float depth = 0.0, end = 10.0;
    for (int i = 0; i < MAX_STEPS; i++) {
        float dist = GetDist(eye + depth * viewRayDirection, uv);
        if (dist < SURF_DIST) {
            // We're inside the scene surface!
            return depth;
        }
        // Move along the view ray
        depth += dist;

        if (depth >= MAX_DIST) {
            // Gone too far; give up
            return MAX_DIST;
        }
    }
    return end;


}

float RayMarchPlan(vec3 eye, vec3 viewRayDirection,vec2 uv){

    float depth = 0.0, end = 10.0;
    for (int i = 0; i < 100; i++) {
        float dist = GetDistPlan(eye + depth * viewRayDirection, uv);
        if (dist < SURF_DIST) {
            // We're inside the scene surface!
            return depth;
        }
        // Move along the view ray
        depth += dist;

        if (depth >= MAX_DIST) {
            // Gone too far; give up
            return MAX_DIST;
        }
    }
    return end;


}

vec3 GetNormal(vec3 p, vec2 uv){
    float d = GetDist(p, uv);
    vec2 e = vec2(0.01, 0);
    
    vec3 n = d - vec3(
        GetDist(p-e.xyy, uv),
        GetDist(p-e.yxy, uv),
        GetDist(p-e.yyx, uv));
        
    return normalize(n);
}



float GetLight(vec3 p, vec2 uv){
    vec3 lightpos = vec3(2.0, 5, -6);
    lightpos.xz += vec2(sin(iTime), cos(iTime));
    vec3 l = normalize(lightpos-p);
    
    
    vec3 n = GetNormal(p, uv);
    
    float dif = clamp(dot(n, l), 0., 1.);
    float d = RayMarch2(p+n*SURF_DIST, l, uv);
    if(d < length(lightpos-p))dif *= 0.1;
    return dif ;

}

vec3 GetNormalPlan(vec3 p, vec2 uv){
    float d = GetDistPlan(p, uv);
    vec2 e = vec2(.01, 0);
    
    vec3 n = d - vec3(
        GetDistPlan(p-e.xyy, uv),
        GetDistPlan(p-e.yxy, uv),
        GetDistPlan(p-e.yyx, uv)*0.9);
        
     
    vec3 tx = texture(iChannel0, uv+iTime*.23 ).rgb*2.0+1.0;
    n = n*500.0 -tx ;
        
    return normalize(n);
}

float GetLightPlan(vec3 p, vec2 uv){
    vec3 lightpos = vec3(2.0, 5, -6);
    lightpos.xz += vec2(sin(iTime), cos(iTime));
    vec3 l = normalize(lightpos-p);
    
    
    vec3 n = GetNormalPlan(p, uv);
    
    float dif = clamp(dot(n, l), 0., 1.);
    float d = RayMarchPlan(p+n*SURF_DIST, l, uv);
    if(d < length(lightpos-p))dif *= 0.1;
    return dif ;

}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Normalized pixel coordinates (from 0 to 1)
    vec2 uv = fragCoord/iResolution.xy;

    // Time varying pixel color
    vec3 col = 0.5 + 0.5*cos(iTime+uv.xyx+vec3(0,2,4));
    
     uv -= 0.5;
     uv /= vec2(iResolution.y / iResolution.x, 1);
    col = vec3(0.0);
    
    vec3 ro, rd, p;
    float d , dif;
    
         
     float fresnel = clamp(1.0 - dot(vec3(1.0, 1.0, 0.0),-vec3(0.0, 1.0, 0.0)), 0.0, 1.0);
    fresnel = pow(fresnel,3.0) * 0.5;
      
      
       ro = vec3(0, 2.0, -5.0);
     rd = normalize(vec3(uv.x, uv.y, 1));
    
     d = RayMarch2(ro, rd, uv);
     p = ro +rd * d;

     dif =GetLight(p, uv);
        
           col += vec3(dif)* vec3(0.2, 1.0, 0.5)*texture(iChannel1, uv*2.-iTime*.01).rgb*2.0 -
           vec3(dif)*vec3(0.0, 1.0, 1.0) + (vec3(173., 79., 9.)/255.) * vec3(dif);
       
       /*2*/
      ro = vec3(0, 4.0, 3.0);
      rd = normalize(vec3(uv.x, uv.y, 1));
    
      d = RayMarchPlan(ro, rd, uv);
      p = ro +rd * d;
    
      dif =GetLightPlan(p, uv);
     if (( dif <= 0.08) && ( dif > 0.05) ){
          col +=  vec3(0.0, 0.8, 1.0);
      }
     else if (( dif <= 0.05)&& ( dif > 0.01)){
         col +=  vec3(0.0, 0.75, 1.0);
      }
     else if (( dif <= 0.01) && ( dif > 0.005)){
         col +=  vec3(0.0, .9, 1.0);
     }
     else if ( dif <= 0.005){
         col +=  vec3(0.0, 1.0, 1.0);
     
     }else{
        col += mix(vec3(0.0, 0.8, 1.0), vec3(dif)* vec3(0.0, 1.0, 1.0), smoothstep(0., 1.0, (dif)));
       
     }
     
     
  
       

    // Output to screen
    fragColor = vec4(col,1.0);
}