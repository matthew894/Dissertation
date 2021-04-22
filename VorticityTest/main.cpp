#include <iostream>
#include <cmath>
#include<random>
#include <fstream>
#include <chrono>
#include <iomanip>
#include<vector>
#include<algorithm>
const int n = 200000;
const int y = 50000;
const double h = double(10.0);
const int LargestCount = n;
const int StepCount = 200;
const double over2pi = double(1)/(double(2)*M_PI);
const double angw = double(1);
const double rnought = double(1.0);
const double wnought = angw/double(n);
double overrnought; //= -double(1)/pow(double(1.1)/sqrt(double(n)),2);
const double half = double(1)/double(2);

std::random_device rd;
//std::mt19937 gen(101); //Original master run
std::mt19937 gen(103);
std::uniform_real_distribution<> dis(-1, 1);



int LargestMethod(const double (&OldPositions)[y][2],double (&Velocity)[y][2], const double (&Vorticity)[n],int x,const double (&asd)[n][2]){
    double maxi =0;
    int element = 0;
    for(int i =0;i<y;i++){

        if(maxi<Velocity[i][0]*Velocity[i][0]+Velocity[i][1]*Velocity[i][1]){
            maxi = Velocity[i][0]*Velocity[i][0]+Velocity[i][1]*Velocity[i][1];
            element = i;
        }



    }

    double a = 0,b=0;
    for(int j =0;j<n;j++){

        if(element!=j) {
            double const xdist = OldPositions[element][0] - asd[j][0];
            double const ydist = -OldPositions[element][1] + asd[j][1];
            double const distancesq = (ydist * ydist + xdist * xdist); //done this way as pow is slow
            double const e = exp(half*overrnought*distancesq);
            double const w =  Vorticity[j]*(double(1) - e) * (double(1)+double(2)*e)/ distancesq ;

            a+= w * ydist;
            b+= w * xdist;
        }
    }
    //std::cout<<1.0-std::abs((Velocity[element][0]*Velocity[element][0]+Velocity[element][1]*Velocity[element][1])/(a*a+b*b))<<"\n";
    if(std::abs(1.0-(Velocity[element][0]*Velocity[element][0]+Velocity[element][1]*Velocity[element][1])/(a*a+b*b))>1.0){
        Velocity[element][0] = a;
        Velocity[element][1] = b;
        return x;
    }
    else{
        return y;
    }


}




int main() {
    //Implementing a fourth order Runge Kutta method
    auto start = std::chrono::steady_clock::now();



    double Vorticity[n];
    double Position[n][2];




    for(int i =0;i<n;i++){
        Position[i][0] = dis(gen);
        Position[i][1] = dis(gen);
        if(Position[i][0]*Position[i][0]+Position[i][1]*Position[i][1]>1.0){
            i--;
        }
        else{
            Vorticity[i] = wnought*exp(-(Position[i][0]*Position[i][0]+Position[i][1]*Position[i][1]));

        }

    }

    double a = 0;

    for(int i =0;i<n;i++){
        double dist = 100000;
        for(int j = 0;j<n;j++){
            if(i!=j){

                if(sqrt(pow(Position[i][0]-Position[j][0],2)+pow(Position[i][1]-Position[j][1],2))<dist){

                    dist = sqrt(pow(Position[i][0]-Position[j][0],2)+pow(Position[i][1]-Position[j][1],2));

                }

            }
        }
        a = a+dist;
    }
    a = a/double(n);
    overrnought = -double(1)/pow(1.1*a,2);




    double x =0;
    double radius = 0.5;

    double velis[y][2];
    double asd[y][2];

    for(int i =0;i<y;i++){
        double w =0;
        double angle = dis(gen)*M_PI*2.0;
        double xrand =  cos(angle)*radius;
        double yrand =  sin(angle)*radius;
        asd[i][0] = xrand;
        asd[i][1] = yrand;
        double a[2];
        velis[i][0] = 0;
        velis[i][1] = 0;
        for(int j=0;j<n;j++){



            double const xdist = xrand - Position[j][0];
            double const ydist = -yrand + Position[j][1];
            double const distancesq = double(1)/(ydist * ydist + xdist * xdist); //done this way as pow is slow
            //double const e = exp(half*overrnought*distancesq);
            //double const w =  Vorticity[j]*(double(1) - e) * (double(1)+double(2)*e)/ distancesq ;
            //double const w = Vorticity[j] * (double(1) - exp(distancesq * overrnought)) / distancesq;
            double const w = Vorticity[j]*distancesq;
            velis[i][0] += w * ydist;
            velis[i][1] += w * xdist;
        }











    }

    for(int i =0;i<y;i++){
        i = LargestMethod(asd,velis,Vorticity,i,Position);

    }


    double mean = 0.19585;
    double t = 0;
    double vels[y];
    for(int i =0;i<y;i++){
        vels[i] = velis[i][0]*velis[i][0]+velis[i][1]*velis[i][1];
        t = t+ vels[i];
    }
    std::cout<<t/double(y)<<"\n";
    double sd =0;
    for(int i =0;i<y;i++){
        sd = sd+pow(vels[i]-mean,2);
    }
    std::cout<<sqrt(sd/double(y))<<"\n";
    std::cout<<sqrt(sd/double(y))/mean*100.0<<"\n";
    /*
    double radius = 0.5;
    for(int i =0;i<1;i++){
        double w =0;
        double angle = dis(gen)*M_PI*2.0;
        double xrand =  cos(angle)*radius;
        double yrand =  sin(angle)*radius;
        if(xrand*xrand+yrand*yrand>1.0){
            i--;
        }
        else{
            double a[2];
            a[0] = 0;
            a[1] = 0;
            for(int j=0;j<n;j++){
                double b = dis(gen);
                double c = dis(gen);
                if(b*b+c*c<1.0) {
                    double const xdist = xrand - b;
                    double const ydist = -yrand + c;
                    double const distancesq =
                            double(1) / (ydist * ydist + xdist * xdist); //done this way as pow is slow
                    double const w = wnought*exp(-(b*b+c*c))* distancesq;

                    a[0] += w * ydist;
                    a[1] += w * xdist;
                }
                else {
                    j--;
                }
            }
            std::cout<<a[0]*a[0]+a[1]*a[1];
        }








    }
    */

    return 0;
}
