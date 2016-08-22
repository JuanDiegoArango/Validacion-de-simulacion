//
//  Catch-Bond.c
//  
//
//  Created by Juan Diego Arango on 3/27/16.
//
//

#include "Catch-Bond.h"
#include <math.h>
#include <stdlib.h>


double Kb_T=4.1;
double k10=6.01;
double k20=0.007;
double k12=0.21;
double k21=0.105;
double x10=1.37*0.1;
double x20=1.76*0.1;
double x12=8.58*0.1;
double x21=-4.2*0.1;
double koff;
double x_off;


double salir_de_estado(int estado,double f, double dt)
{   double respuesta;

    if (estado==1.0)
    {
        respuesta=k10*exp(x10*f/(Kb_T))*dt;
    }
    else if (estado==2.0)
    {
        
        respuesta=k20*exp(x20*f/(Kb_T))*dt;
    }
    return respuesta;
}

double cambiar_de_estado(int estado, double f,double dt)
{   double respuesta;
    if (estado==2.0)
    {
        respuesta=k21*exp(x21*f/(Kb_T))*dt;
    }
    else if (estado==1.0)
    {
        respuesta=k12*exp(x12*f/(Kb_T))*dt;
    }
    return respuesta;
}


double entrar(double f,double dt)

{   double respuesta=koff*exp(x_off*f/(Kb_T))*dt;
    return respuesta;
}

void Catch_Bond(int semilla,double fuerza,int *estado, double dt, double tiempo)
{
    
    double t=0;
    double p;
    double ktot;
    double koff;
    double kon;

    while (t<=dt)
    {
        srand48(semilla);
        p=drand48();
    
        if (estado==1)
        {
            ktot=k12*exp(x12*f/(Kb_T))+k10*exp(x10*f/(Kb_T));
            koff=k10*exp(x10*f/(Kb_T))/ktot;
            kon=k21*exp(x12*f/(Kb_T))/ktot;
            
            if(p<koff)
            {   cambiar_de_estado
            
            }
            
            else{}
            
        }
    
        else
        {
            ktot=k21*exp(x21*f/(Kb_T)+k20*exp(x20*f/(Kb_T);
            koff=k20*exp(x20*f/(Kb_T)/ktot;
            kon=k21*exp(x21*f/(Kb_T)/ktot;
        }
        
        


    
}

    
    
return estado;
    
}



