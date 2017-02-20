#include <QCoreApplication>

#include <iomanip>
#include <ctime>
#include <chrono>
#include <iostream>
#include <fstream>
#include <string>

#include <utils.cpp>

using namespace  std;


polinom p(5);

double f(double x)
{
    return p.f(x);
}

double intergal(double a,double b)
{
   return p.intergal(a,b);
}
double deriv(double x)
{
   return p.deg*p.f(x)/x;
}
double deriv2(double x)
{
   return p.deg*(p.deg-1)*p.f(x)/(x*x);
}
double deriv3(double x)
{
   return p.deg*(p.deg-1)*(p.deg-2)*p.f(x)/(x*x*x);
}
double integLeft(double a,double b,int n)
{
    double d=(b-a)/n;
    double sum=0;
    for(int i=0; i<n; i++)
    {
      sum+=f(a+i*d);
    }
    return sum*d;
}

double integTrap(double a,double b,int n)
{
    double d=(b-a)/n;
    double sum=0;
    for(int i=1; i<n; i++)
    {
      sum+=f(a+i*d);
    }
    return (f(a)+2*sum+f(b))*d/2;
}

double integSimson(double a,double b,int n)
{
    double d=(b-a)/(2*n);
    double sum1=0;
    double sum2=0;
    for(int i=1; i<2*n; i++)
    {
    if(i%2==1) sum1+=f(a+i*d);
    if(i%2==0) sum2+=f(a+i*d);
    }

    return (f(a)+4*sum1+2*sum2+f(b))*d/3;

}

double Rn_1(double a,double b,int n)
{
    return (b-a)*((b-a)/n)*deriv(b)/2;
}

double Rn_2(double a,double b,int n)
{
    return (b-a)*((b-a)/n)*((b-a)/n)*deriv2(b)/12;
}
double Rn_3(double a,double b,int n)
{
    return (b-a)*((b-a)/n)*((b-a)/n)*((b-a)/n)*deriv3(b)/2880;
}
int main(int argc, char *argv[])
{

double a=0;
double b=10;
int n=3;

ofstream fout;
fout.open("/home/skutukov/integrate/res/log.txt",ios::app);
double I=intergal(a,b);
double ILN=integLeft(a,b,n);
double IL2N=integLeft(a,b,2*n);

double R_n=Rn_1(a,b,n);
double R_2n=Rn_1(a,b,2*n);
double R_main=IL2N-ILN;
double Iyt=IL2N+R_main;

fout<<endl;
fout<<"degree      "<<p.deg<< '\t'<<"I_N"<<'\t'<<"I_N-I"<<'\t' <<"rN"<<'\t'<<"I_2N"<<'\t'<<"I_L2N-I"<<'\t'<<"R2N"<<'\t'<<"R_гл"<<'\t'<<"Iyt"<<'\t'<<"I-Iyt"<<endl;
fout<<"Left rectangle"<< '\t'<<ILN<<'\t'<<(ILN-I)<<'\t' <<R_n<<'\t'<<IL2N<<'\t'<<(IL2N-I)<<'\t'<<R_2n<<'\t'<<R_main <<'\t'<<Iyt<<'\t'<<I-Iyt<<endl;

I=intergal(a,b);
ILN=integTrap(a,b,n);
IL2N=integTrap(a,b,2*n);

R_n=Rn_2(a,b,n);
R_2n=Rn_2(a,b,2*n);
R_main=(IL2N-ILN)/3;
Iyt=IL2N+R_main;

fout<<"trapeze       "<< '\t'<<ILN<<'\t'<<(ILN-I)<<'\t' <<R_n<<'\t'<<IL2N<<'\t'<<(IL2N-I)<<'\t'<<R_2n<<'\t'<<R_main <<'\t'<<Iyt<<'\t'<<I-Iyt<<endl;

I=intergal(a,b);

ILN=integSimson(a,b,n);
IL2N=integSimson(a,b,2*n);
n=2*n;
R_n=Rn_3(a,b,n);
R_2n=Rn_3(a,b,2*n);
R_main=(IL2N-ILN)/7;
Iyt=IL2N+R_main;

fout<<"Symson         "<< '\t'<<ILN<<'\t'<<(ILN-I)<<'\t' <<R_n<<'\t'<<IL2N<<'\t'<<(IL2N-I)<<'\t'<<R_2n<<'\t'<<R_main <<'\t'<<Iyt<<'\t'<<I-Iyt<<endl;
fout.close();


    return 0;
}
