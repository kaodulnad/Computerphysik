#include<math.h>
#include<stdio.h>
#include<stdlib.h>

double potential(double Q, double a, double x, double z  ){
    double pot = Q* exp(-x*x/a*a)/(a*sqrt(3.14159265359)*sqrt(x*x+z*z));
    //double pot = exp(x);
    //printf(" %f\n", pot);
    return pot;
}

double trapez(double a, double b, int n,double (*function)(double, double,double,double),double Q, double alpha, double z ){ //alpha gleich a in Potential
    double sum = 0;
    double h = (b -a)/(n -1);
    

    sum += (h/2*(((*function)(Q ,alpha ,a, z ))+((*function)(Q ,alpha ,b, z )) ));

    for(int i = 1; i < n-1; i++){
        sum += h*((*function)(Q ,alpha ,a+i*h, z ) );
    }
    //printf("\nH: %f, sum %f\n", h, sum);
    return sum;
}

int indices(int i, int k, int mmax){
    return (i-k + k*(mmax+1)-k*(k-1)/2);
}

double romberg(double a, double b, double (*function)(double, double,double,double),double Q, double alpha, double z, double fehler, int mmax, int *n0){
    double result;

    double *h = (double*)malloc(sizeof(double)*(mmax+1));
    double *T = (double*)malloc(sizeof(double)*(mmax+1));
    double *Ttilde = (double*)malloc(sizeof(double)*(mmax+1)*(mmax+2)/2);

    //Nulltes Element berechnen
    int n = *n0;

    *h = (b - a)/(double)(*n0-1);
    *T = trapez(a, b, n, function, Q, alpha, z);
    //printf("            tapez : %f\n", *T);
    *Ttilde = *T;

    int m = 1;
    for(m = 1; m<= mmax; m++){

        n = 2*n;
        *(h+m) = (b-a)/(double)(n-1);
        //printf("h: %f,n %i", *(h+m), n);
        *(T+m) = trapez(a, b, n, function, Q, alpha, z);
        *(Ttilde + indices(m, 0, mmax))= *(T +m);
        //printf(" trapez %f ",*(Ttilde + indices(m, 0, mmax)) );

        for(int k = 1; k <= m ; k++){
            *(Ttilde+ indices(m , k , mmax))= -(*(h + m-k))*(*(h + m-k))/((*(h + m ))*(*(h + m ))-(*(h + m-k))*(*(h + m-k)))*(*(Ttilde + indices(m  , k-1, mmax))) 
                                              +(*(h + m)  )*(*(h + m)  )/((*(h + m ))*(*(h + m ))-(*(h + m-k))*(*(h + m-k)))*(*(Ttilde + indices(m-1, k-1, mmax)));
            //printf(" ttilde %f, k: %i, m: %i", *(Ttilde+ indices(m , k , mmax)), k, m);
        if(fabs(*(Ttilde + indices(m  , k, mmax))- *(Ttilde + indices(m-1  , k-1, mmax))) <= fehler){
            double ret = (*(Ttilde + indices(m  , k, mmax)));

            free(h);     
            free(Ttilde);
            free(T); 
            printf(" ping, %i ", m);          
            return  ret;
        }
        } 

    if(m > mmax) {
    printf("Konvergenz nicht erreicht innerhalb mmax = %d Schritte!\n", mmax);
    }
    *n0 = n;


    printf("  \n");
    result = *(Ttilde +indices((int) fmin(m,mmax),(int) fmin(m,mmax), mmax));
    }
    // free(h);
    // free(T);
    // free(Ttilde);
    printf(" romberg verlassen mit result %f, %i, %i \n", result,(int) fmin(m,mmax), m);
    free(h);     
    free(Ttilde);
    free(T);
    return result;
}

int main( ){
    double sum = 0;
    int n = 100;
    sum = romberg(-1000, 1000, &potential, 1, 4, 1, 0.0000000000005, 100, &n);
    printf(" %f\n", sum);
    
    printf(" %f", potential(1, 1, 1, 1));



    return 0;
}

