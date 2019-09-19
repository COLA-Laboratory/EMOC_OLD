#include "../headers/global.h"



void ctp1 (SMRT_individual *ind)
{
    double m;
    double *xreal, *obj;
    double g,h,f1=0,f2=0;
    int i;

    obj = ind->obj;
    xreal = ind->variable;
    for (i = 1;i < g_algorithm_entity.algorithm_para.objective_number;i++)
    {
        f1 += pow(xreal[i],2.0);
        f2 += xreal[i];
    }
    obj[0]=xreal[0];
    g=11+f1-10*cos(2*PI*f2);
    h=exp(-obj[0]/g);
    obj[1] = g*h;//f2(x)

//the constraint numbers is J;
    int j,J;
    J=2;
    double a[J+1],b[J+1],C[J+1];
    double alfa,beta,delta;
    double a1, b1, a2, b2, C1, C2;
    delta=1.0/(J+1.0);
    alfa=delta;
    for (j=0;j<J;j++)
    {
        if(j==0)
        {
            a[j]=1.0;
            b[j]=1.0;
        }
        beta=a[j]*exp(-b[j]*alfa);
        a[j+1]=(a[j]+beta)/2.0;
        b[j+1]=-1.0*log(beta/a[j+1])/alfa;
        alfa=alfa+delta;
    }


    ind->cv=0;
    for (j=0;j<J;j++)
    {
        C[j+1]=obj[1]-a[j+1]*exp(-b[j+1]*obj[0]);
        if(C[j+1]>0||C[j+1]==0)
            m=0;
        else
            m=C[j+1];
        ind->cv += m;
    }

}


void ctp2 (SMRT_individual *ind)
{
    double m;
    double *xreal, *obj;
    double g,h,f=0;
    int i;

    obj = ind->obj;
    xreal = ind->variable;
    for(i = 1; i < g_algorithm_entity.algorithm_para.objective_number; i++)
    {
        f += xreal[i];
    }
    g=1+f;
    obj[0] = xreal[0];
    h=1-obj[0]/g;       //f1(x)
    obj[1] = g*h;       //f2(x)

//the constraint numbers is one;
    double C, theta, a,b,c,d,e;
    theta=-0.2*PI;
    a=0.2;
    b=10;
    c=1;
    d=6;
    e=1;
    C=cos(theta)*(obj[1]-e)-sin(theta)*obj[0]-a*pow(sin(b*PI*pow(sin(theta)*(obj[1]-e)+cos(theta)*obj[0],c)),d);
    if(C<0)
    {
        m=C;
    }
    else
    {
        m=0;
    }
    ind->cv=m;
}

void ctp3(SMRT_individual *ind)
{
    double m;
    double *xreal, *obj;
    double g,h,f=0;
    int i;
    obj = ind->obj;
    xreal = ind->variable;
    for(i=1;i<g_algorithm_entity.algorithm_para.objective_number;i++)
    {
        f += xreal[i];
    }
    g=1+f;
    obj[0] = xreal[0];  //f1(x)
    h=1-obj[0]/g;
    obj[1] = g*h;       //f2(x)

//the constraint numbers is one;
    double C, theta, a,b,c,d,e;
    theta=-0.2*PI;
    a=0.1;
    b=10;
    c=1;
    d=0.5;
    e=1;
    C=cos(theta)*(obj[1]-e)-sin(theta)*obj[0]-a*pow(fabs(sin(b*PI*pow(sin(theta)*(obj[1]-e)+cos(theta)*obj[0],c))),d);
    if(C<0)
    {
        m=C;
    }
    else
    {
        m=0;
    }
    ind->cv=m;
}

void ctp4(SMRT_individual *ind)
{
    double m;
    double *xreal, *obj;
    double g,h,f=0;
    int i;

    obj = ind->obj;
    xreal = ind->variable;

    for(i = 1; i < g_algorithm_entity.algorithm_para.objective_number; i++)
    {
        f += xreal[i];
    }
    g=1+f;
    obj[0] = xreal[0];  //f1(x)
    h=1-obj[0]/g;
    obj[1] = g*h;       //f2(x)

//the constraint numbers is one;
    double C, theta, a,b,c,d,e;
    theta=-0.2*PI;
    a=0.75;
    b=10;
    c=1;
    d=0.5;
    e=1;
    C=cos(theta)*(obj[1]-e)-sin(theta)*obj[0]-a*pow(fabs(sin(b*PI*pow(sin(theta)*(obj[1]-e)+cos(theta)*obj[0],c))),d);
    if(C<0)
    {
        m=C;
    }
    else
    {
        m=0;
    }
    ind->cv=m;
}

void ctp5(SMRT_individual *ind)
{
    double m;
    double *xreal, *obj;
    double g,h,f=0;
    int i;

    obj = ind->obj;
    xreal = ind->variable;
    for(i = 1; i < g_algorithm_entity.algorithm_para.objective_number; i++)
    {
        f += xreal[i];
    }
    g=1+f;
    h=1-obj[0]/g;
    obj[0] = xreal[0];  //f1(x)
    obj[1] = g*h;       //f2(x)

//the constraint numbers is one;
    double C, theta, a,b,c,d,e;
    theta=-0.2*PI;
    a=0.1;
    b=10;
    c=2;
    d=0.5;
    e=1;
    C=cos(theta)*(obj[1]-e)-sin(theta)*obj[0]-a*pow(fabs(sin(b*PI*pow(sin(theta)*(obj[1]-e)+cos(theta)*obj[0],c))),d);
    if(C<0)
    {
        m=C;
    }
    else
    {
        m=0;
    }
    ind->cv=m;
}

void ctp6(SMRT_individual *ind)
{
    double m;
    double *xreal, *obj;
    double g,h,f=0;
    int i;

    obj = ind->obj;
    xreal = ind->variable;
    for(i=1;i<g_algorithm_entity.algorithm_para.objective_number;i++)
    {
        f += xreal[i];
    }
    g=1+f*5;
    obj[0] = xreal[0];  //f1(x)
    h=1-obj[0]/g;
    obj[1] = g*h;       //f2(x)

//the constraint numbers is one;
    double C, theta, a,b,c,d,e;
    theta=0.1*PI;
    a=40;
    b=0.5;
    c=1;
    d=2;
    e=-2;
    C=cos(theta)*(obj[1]-e)-sin(theta)*obj[0]-a*pow(fabs(sin(b*PI*pow(sin(theta)*(obj[1]-e)+cos(theta)*obj[0],c))),d);
    if(C<0)
    {
        m=C;
    }
    else
    {
        m=0;
    }
    ind->cv=m;
}

void ctp7(SMRT_individual *ind)
{
    double m1, m2;
    double *xreal, *obj;
    double g,h,f=0;
    int i;

    obj = ind->obj;
    xreal = ind->variable;
    for(i=1;i<g_algorithm_entity.algorithm_para.objective_number;i++)
    {
        f += xreal[i];
    }
    g=1+f;
    obj[0] = xreal[0];  //f1(x)
    h=1-obj[0]/g;
    obj[1] = g*h;       //f2(x)

//the constraint numbers is one;
    double C, theta, a,b,c,d,e;
    theta=-0.05*PI;
    a=40;
    b=5;
    c=1;
    d=6;
    e=0;
    C=0.9877*(obj[1])+0.1564*obj[0]-a*pow(fabs(sin(b*PI*pow(-0.1564*(obj[1])+0.9877*obj[0],c))),d);
    if(C<0)
    {
        m1=C;
    }
    else
    {
        m1=0;
    }
    m2=obj[1]-1+obj[0];
    if(m2>0)
    {
        m2=0;
    }
    else
    {
        m2=-100000;
    }
    ind->cv=m1+m2;
}


void ctp8(SMRT_individual *ind)
{
    double m;
    double *xreal, *obj;
    double g,h,f=0;
    int i;

    obj = ind->obj;
    xreal = ind->variable;
    for(i=1;i<g_algorithm_entity.algorithm_para.objective_number;i++)
    {
        f += xreal[i];
    }
    g=1+f*5;
    obj[0] = xreal[0];  //f1(x)
    h=1-obj[0]/g;
    obj[1] = g*h;       //f2(x)

//the constraint numbers is one;
    double C1, theta1, a1,b1,c1,d1,e1;
    double C2, theta2, a2,b2,c2,d2,e2;

    theta1=0.1*PI;
    a1=40;
    b1=0.5;
    c1=1;
    d1=2;
    e1=-2;
    C1=cos(theta1)*(obj[1]-e1)-sin(theta1)*obj[0]-a1*pow(sin(b1*PI*pow(sin(theta1)*(obj[1]-e1)+cos(theta1)*obj[0],c1)),d1);

    theta2=-0.05*PI;
    a2=40;
    b2=2;
    c2=1;
    d2=6;
    e2=0;
    C2=cos(theta2)*(obj[1]-e2)-sin(theta2)*obj[0]-a2*pow(sin(b2*PI*pow(sin(theta2)*(obj[1]-e2)+cos(theta2)*obj[0],c2)),d2);

    if(C1>0||C1==0)
        m=0;
    else
        m=C1;

    if(C2>0||C2==0)
        ind->cv=m;
    else
        ind->cv=m+C2;
}

