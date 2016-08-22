//
//  Pilis.c
//  
//
//  Created by Juan Diego Arango on 3/26/16.
//
//

#include "Pilis.h"
#include <math.h>
#include <stdlib.h>

//parametros
double Kb_T=4.10;
//longitudes
double La0=0.75;
double Lb0=5.7;
double La=La0;
double lp=2.7;
//constante resorte
double ka=2.0;
double dt;



double fuerza(double x, double c)
{
    return (Kb_T/(4.0*lp))*((1.0-(x/c))*-2.0-1.0+4.0*(x/c));
}
double probabilidad_de_transicion (double v, double f,double delta_X)
{
    return v*exp((delta_X*f)/(Kb_T))*dt;
}
double nuevo_La(double f, double Numero)
{
    return (f/ka)+La0*Numero;
}




void  pilis(int *Nb,int *Na,int *posicion_libre,double L,int semilla)
{
    
double kAB=0.05;
double kBA=700.0;
double delta_X_AB=0.5;
double delta_X_BA=-0.5;
//monomeros totales.
int Ntot=Na+Nb;
//arrays
double monomeros[Ntot];
//inicializacion de array.
for(int i=posicion_libre; i<Ntot;i=+1)
{
        if (i<Ntot && posicion_libre<=i)
        {monomeros[i]=1;}
        
        else if (i<posicion_libre && 0<=i)
        {monomeros[i]=2;}
        
}

//solucion cubica
double a=La0*Na;
double b=Kb_T/(lp);
double c=Lb0*Nb;
double cp[4];
cp[0]=4*(b+c*ka);
cp[1]=-(c)*(8*c*ka+4*ka*(L-a)+9*b);
cp[2]=(2*c*2)*(2*c*ka+4*ka*(L-a)+3*b);
cp[3]=4*(c*3)*(ka)*(a-L);

double delta0=cp[1]*cp[1]-3*cp[0]+cp[2];
double delta1=2*cp[1]*cp[1]*cp[1]-9*cp[0]*cp[1]*cp[2]+27*cp[0]*cp[0]*cp[3];
double C=pow(((delta1+sqrt(delta0*delta0-4*delta1*delta1*delta1))*0.5),1.0/3.0);
double Lb=-1.0/(3*cp[0])*(cp[1]+C+delta0/C);

double f= fuerza(Lb,c);
    
//nueva longitud.
    
    
La=nuevo_La(f,Na);
    
//numeros alteatorio:
srand48(semilla);
double p1=drand48();
double p2=drand48();
double pa=probabilidad_de_transicion(kAB,f,delta_X_AB);
double pb=probabilidad_de_transicion(kBA,f,delta_X_BA);


if (pa>p1 && (Na)>0)
{
    monomeros[posicion_libre]=2;
    posicion_libre=posicion_libre+1;
}

Na=0;
Nb=0;
for (int i=0 ; i<Ntot; i=+1){
        if(monomeros[i]==1)
        {
            Na=+1;
        }
        
        
        if(monomeros[i]==2)
        {
            Na=+2;
        }
        
        
    }
    
if ((pb>p2) && (Nb>7) && (posicion_libre>7))
{
    posicion_libre=posicion_libre-1;
    monomeros[posicion_libre]=1;
    }
    
    

}




