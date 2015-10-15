// simpson function which perfoms numerical integration
// of a time series using simpson's formula
// \int_{x_{i}}^{x_{i+2}}f(x)dx ~ \Delta x /3 (f(x_i)+4f(x_{i+1})+f(x_{i+1}) )

//important! for correct integration length
// must be odd
// parameters: series - initial time seies
//             len - length of the time series
/*
void simpson(double *series, unsigned int *len, double *sum){
    *sum = 0;
    if (*len % 2 == 0){
        //unsigned int l = *len-1;
        for(unsigned int i=0; i< *len-2; i=i+2)
            *sum = *sum+(series[i] + 4*series[i+1] + series[i+2])/3;
        *sum = *sum + (series[*len-1]+series[*len-2])/2;
    }
    else{
        for(unsigned int i=0; i< *len-2; i=i+2)
            *sum = *sum+(series[i] + 4*series[i+1] + series[i+2])/3;
    }
}
*/

void simpson (double *series, double *len, double *sum){
double h = 1/(*len);
//double h = 1;
double s0,s1,s2;

s0 = 0;
s1 = 0;
s2 = 0;

for (unsigned int i = 1; i < (int)*len-1; i = i+2){
    s0 = s0+series[i];
    s1 = s1+series[i-1];
    s2 = s2+series[i+1];
}
*sum = h*(s1+4*s0+s2)/3;
/* If n is even, add the last slice separately */
if (((int)*len) % 2 == 0){*sum = *sum+h*(5*series[(int)*len-1]+8*series[(int)*len-2]-series[(int)*len-3])/12;}
}


