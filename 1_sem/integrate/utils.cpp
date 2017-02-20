#include <QCoreApplication>

class polinom
{
public:
    int deg;
    polinom(int deg1){
        deg=deg1;
    }

    double f(double x)
    {
        return pow(x,deg);
    }


    double intergal(double a,double b)
    {
       return pow(b,deg+1)/(deg+1)-pow(a,deg+1)/(deg+1);

    }
};
