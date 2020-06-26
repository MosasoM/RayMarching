float pi_function(float x, float y){
    y = y*2.0;
    x = x*2.0+1.5;
    float f1 = (1.5*exp(-0.62*(y-0.16)*(y-0.16)))/(1.0+exp(-20.0*(5.0*y-1.0)));
    float f21 = (1.5+0.8*pow((y-0.2),3.0));
    float f22 = (1.0+exp(20.0*(5.0*y-1.0)))*(1.0+exp(-(100.0*(y+1.0)-16.0)));
    float f2 = f21/f22;
    float f3 = (0.2*(exp(-(y+1.0)*(y+1.0))+1.0))/(1.0+exp(100.0*(y+1.0)-16.0));
    float f4 = 0.1/exp(2.0*pow((10.0*y-1.2),4.0));
    return x-(f1+f2+f3+f4);
}