//#define webcam
void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // Normalized pixel coordinates (from 0 to 1)
    
    vec2 m = iMouse.xy/iResolution.xy;
    float zoom = pow(5., -m.x*10.);
    
    vec2 uv = (fragCoord.xy-.5*iResolution.xy) / iResolution.y;
    uv *= 1.0;
    vec3 c = vec3(uv, 1.0)*zoom*1.;
    c+= vec3(.1, .2681, .2681);
                         
    
    vec3 z = vec3(0.);
    float iter = 0.;
    
    
    const float max_iter = 200.;
    
    for(float i = 0.; i < max_iter; i++){
        //z = vec2(z.x*z.x - 4.0*z.y*z.y, 4.*z.x*z.y)+c;
        z = vec3(z.x * z.x - z.y * z.y, 2.*z.x*z.y + 2.*z.y*z.z, 2.0*z.x*z.z + z.z * z.z)+c;
        if(length(z) > 2.)break;
                
        iter++;
    
    }
    
    float f = iter/max_iter;
    
    vec3 col = vec3(f);
    
    #ifdef webcam
    if(f > 0.5){
        uv+=0.5;
        col = texture(iChannel0, uv).rgb;
    }
    #endif

    // Output to screen
    fragColor = vec4(col,1.0);
}