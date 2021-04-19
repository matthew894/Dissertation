#include <iostream>
#include <cmath>
#include<random>
#include <fstream>
#include <chrono>
#include <iomanip>

const int n = 20000;
const double h = double(1.0);
const int LargestCount = 0;
const int StepCount = 100;
const double over2pi = double(1)/(double(2)*M_PI);
const double angw = double(1);
const double rnought = double(1.0);
const double wnought = angw*over2pi*double(1)/double(n);
double overrnought= -double(1)/pow(double(1.1)/sqrt(double(n)),2);
const double half = double(1)/double(2);

std::random_device rd;
//std::mt19937 gen(101); //Original master run
std::mt19937 gen(103);
std::uniform_real_distribution<> dis(-1, 1);
std::ofstream myfile;



double hamiltonian(const double(&Positions)[n][2], const double (&Vorticity)[n]){

    double ham = 0;
    static double xdist,ydist;

    for(int i =0;i<n;i++){

        for(int j =0;j<n;j++){

            if(i!=j) {

                 xdist = Positions[i][0] - Positions[j][0];
                 ydist = Positions[i][1] - Positions[j][1];

                ham = ham + Vorticity[j]*Vorticity[i]* log(xdist * xdist + ydist * ydist);

            }
        }
    }

    return -double(1)/double(4)*ham;
}

double AngularImpulse(const double(&Positions)[n][2], const double (&Vorticity)[n]){
    double impulse = 0;

    for(int i =0;i<n;i++){

        impulse = impulse+Vorticity[i]*(Positions[i][0]*Positions[i][0]+Positions[i][1]*Positions[i][1]);

    }

    return impulse;

}

double LinearImpulse(const double (&Positions)[n][2], const double (&Vorticity)[n]){
    double impulse = 0;

    for(int i =0;i<n;i++){

        impulse = impulse+Vorticity[i]*(Positions[i][0]+Positions[i][1]);

    }

    return impulse;
}


void ZeroOrder(const double (&OldPositions)[n][2],double (&Velocity)[n][2], const double (&Vorticity)[n]){

    static double xdist,ydist,distancesq,w;


    for(int i =0;i<n;i++){

        Velocity[i][0]=double(0);
        Velocity[i][1]= double(0);
    }

    for(int i =0;i<n;i++){

        for(int j =i+1;j<n;j++){

            xdist = OldPositions[i][0] -OldPositions[j][0];
            ydist = -OldPositions[i][1] + OldPositions[j][1];
            distancesq = double(1)/(ydist * ydist + xdist * xdist); //done this way as pow is slow

            w = Vorticity[j] *distancesq;
            Velocity[i][0] += w * ydist ;
            Velocity[i][1] +=  w * xdist ;

            w = -Vorticity[i] * distancesq;
            Velocity[j][0] += w * ydist ;
            Velocity[j][1] +=  w * xdist ;

        }
    }
}

void ParalellZeroOrder(const double (&OldPositions)[n][2],double (&Velocity)[n][2], const double (&Vorticity)[n]){

    for(int i =0;i<n;i++){

        Velocity[i][0]=double(0);
        Velocity[i][1]= double(0);

    }

#pragma omp parallel for
    for(int i =0;i<n;i++){

        for(int j =0;j<n;j++){

            if(i!=j) {

                double const xdist = OldPositions[i][0] - OldPositions[j][0];
                double const ydist = -OldPositions[i][1] + OldPositions[j][1];
                double const distancesq = double(1)/(ydist * ydist + xdist * xdist); //done this way as pow is slow
                double const w = Vorticity[j]  * distancesq;

                Velocity[i][0] += w * ydist;
                Velocity[i][1] += w * xdist;

            }

        }
    }
}
/*
void OtherParalellZeroOrder(const double (&OldPositions)[n][2],double (&Velocity)[n][2], const double (&Vorticity)[n]){

     static double x[n][n][2] ;

    for(int i =0;i<n;i++){

        Velocity[i][0]=double(0);
        Velocity[i][1]= double(0);

    }

#pragma omp parallel for
    for(int i =0;i<n;i++){

        for(int j =i+1;j<n;j++){

            double const xdist = OldPositions[i][0] - OldPositions[j][0];
            double const ydist = -OldPositions[i][1] + OldPositions[j][1];
            double const distancesq = double(1)/(ydist * ydist + xdist * xdist); //done this way as pow is slow

            x[j][i][0] = distancesq * ydist;
            x[j][i][1] = distancesq * xdist;

        }
    }

#pragma omp parallel for
    for(int i =0;i<n;i++){

        for(int j =0;j<i;j++){

            Velocity[i][0]-= Vorticity[j]* x[i][j][0];
            Velocity[i][1]-= Vorticity[j]* x[i][j][1];

        }

        for(int j = i+1;j<n;j++){
            Velocity[i][0]+= Vorticity[j]* x[j][i][0];
            Velocity[i][1]+= Vorticity[j]* x[j][i][1];

        }



    }
}
 */
void RegFunc(const double (&OldPositions)[n][2],double (&Velocity)[n][2], const double (&Vorticity)[n]){

    static double xdist,ydist,distancesq,w;

    for(int i =0;i<n;i++){

        Velocity[i][0]=double(0);
        Velocity[i][1]= double(0);

    }

    for(int i =0;i<n;i++){

        for(int j =i+1;j<n;j++){

                xdist = OldPositions[i][0] - OldPositions[j][0];
                ydist = -OldPositions[i][1] + OldPositions[j][1];
                distancesq = (ydist * ydist + xdist * xdist); //done this way as pow is slow
                distancesq = (double(1) - exp(distancesq * overrnought)) / distancesq;

                w = Vorticity[j]* distancesq;
                Velocity[i][0] += w* ydist ;
                Velocity[i][1] +=  w * xdist ;

                w = Vorticity[i]* distancesq;
                Velocity[j][0] -= w * ydist ;
                Velocity[j][1] -= w * xdist ;


        }
    }
}
void ParalellRegFunc(const double (&OldPositions)[n][2],double (&Velocity)[n][2], const double (&Vorticity)[n]){



    for(int i =0;i<n;i++){

        Velocity[i][0]=0;
        Velocity[i][1]= 0;

    }
#pragma omp parallel for
    for(int i =0;i<n;i++){

        for(int j =0;j<n;j++){

            if(i!=j) {
                double const xdist = OldPositions[i][0] - OldPositions[j][0];
                double const ydist = -OldPositions[i][1] + OldPositions[j][1];
                double const distancesq = (ydist * ydist + xdist * xdist); //done this way as pow is slow
                double const w = Vorticity[j] * (double(1) - exp(distancesq * overrnought)) / distancesq;

                Velocity[i][0] += w * ydist;
                Velocity[i][1] += w * xdist;
            }

        }
    }
}
/*
void OtherParalellRegFunc(const double (&OldPositions)[n][2],double (&Velocity)[n][2], const double (&Vorticity)[n]){

    static double x[n][n][2] ;

    for(int i =0;i<n;i++){

        Velocity[i][0]=double(0);
        Velocity[i][1]= double(0);

    }

#pragma omp parallel for
    for(int i =0;i<n;i++){

        for(int j =i+1;j<n;j++){

            double const xdist = OldPositions[i][0] - OldPositions[j][0];
            double const ydist = -OldPositions[i][1] + OldPositions[j][1];
            double distancesq = (ydist * ydist + xdist * xdist);
            distancesq = (double(1) - exp(distancesq * overrnought)) / distancesq;

            x[j][i][0] = distancesq * ydist;
            x[j][i][1] = distancesq * xdist;

        }
    }

#pragma omp parallel for
    for(int i =0;i<n;i++){

        for(int j =0;j<i;j++){

            Velocity[i][0]-= Vorticity[j]* x[i][j][0];
            Velocity[i][1]-= Vorticity[j]* x[i][j][1];

        }

        for(int j = i+1;j<n;j++){
            Velocity[i][0]+= Vorticity[j]* x[j][i][0];
            Velocity[i][1]+= Vorticity[j]* x[j][i][1];

        }



    }
}
 */
void  blobvortex(const double (&OldPositions)[n][2],double (&Velocity)[n][2], const double (&Vorticity)[n]){

    static double xdist,ydist,distancesq,w,e;
    for(int i =0;i<n;i++){

        Velocity[i][0]=0;
        Velocity[i][1]= 0;

    }

    for(int i =0;i<n;i++){

        for(int j =i+1;j<n;j++){

                xdist = OldPositions[i][0] - OldPositions[j][0];
                ydist = -OldPositions[i][1] + OldPositions[j][1];
                distancesq = (ydist * ydist + xdist * xdist); //done this way as pow is slow
                e = exp(half*overrnought*distancesq);
                distancesq =  (double(1) - e) * (double(1)+double(2)*e)/ distancesq ;
                w = Vorticity[j]*distancesq;
                Velocity[i][0] += w * ydist ;
                Velocity[i][1] += w * xdist ;
                w = Vorticity[i]*distancesq;
                Velocity[j][0] -= w * ydist ;
                Velocity[j][1] -= w * xdist ;

        }
    }
}

void  Parallelblobvortex(const double (&OldPositions)[n][2],double (&Velocity)[n][2], const double (&Vorticity)[n]){

    for(int i =0;i<n;i++){

        Velocity[i][0]=0;
        Velocity[i][1]= 0;

    }

#pragma omp parallel for
    for(int i =0;i<n;i++){

        for(int j =0;j<n;j++){

            if(i!=j) {
                double const xdist = OldPositions[i][0] - OldPositions[j][0];
                double const ydist = -OldPositions[i][1] + OldPositions[j][1];
                double const distancesq = (ydist * ydist + xdist * xdist); //done this way as pow is slow
                double const e = exp(half*overrnought*distancesq);
                double const w =  Vorticity[j]*(double(1) - e) * (double(1)+double(2)*e)/ distancesq ;

                Velocity[i][0] += w * ydist;
                Velocity[i][1] += w * xdist;
            }
        }
    }
}
/*
void OtherParalellblobvortex(const double (&OldPositions)[n][2],double (&Velocity)[n][2], const double (&Vorticity)[n]){

    static double x[n][n][2] ;

    for(int i =0;i<n;i++){

        Velocity[i][0]=double(0);
        Velocity[i][1]= double(0);

    }

#pragma omp parallel for
    for(int i =0;i<n;i++){

        for(int j =i+1;j<n;j++){

            double const xdist = OldPositions[i][0] - OldPositions[j][0];
            double const ydist = -OldPositions[i][1] + OldPositions[j][1];
            double distancesq = (ydist * ydist + xdist * xdist);
            double const e = exp(half*overrnought*distancesq);
            distancesq = (double(1) - e) * (double(1)+double(2)*e) / distancesq;

            x[j][i][0] = distancesq * ydist;
            x[j][i][1] = distancesq * xdist;

        }
    }

#pragma omp parallel for
    for(int i =0;i<n;i++){

        for(int j =0;j<i;j++){

            Velocity[i][0]-= Vorticity[j]* x[i][j][0];
            Velocity[i][1]-= Vorticity[j]* x[i][j][1];

        }

        for(int j = i+1;j<n;j++){
            Velocity[i][0]+= Vorticity[j]* x[j][i][0];
            Velocity[i][1]+= Vorticity[j]* x[j][i][1];

        }



    }
}
*/


int LargestMethod(const double (&OldPositions)[n][2],double (&Velocity)[n][2], const double (&Vorticity)[n],int x){
    double maxi =0;
    int element = 0;
    for(int i =0;i<n;i++){

        if(maxi<Velocity[i][0]*Velocity[i][0]+Velocity[i][1]*Velocity[i][1]){
            maxi = Velocity[i][0]*Velocity[i][0]+Velocity[i][1]*Velocity[i][1];
            element = i;
        }



    }

    double a = 0,b=0;
    for(int j =0;j<n;j++){

        if(element!=j) {
            double const xdist = OldPositions[element][0] - OldPositions[j][0   ];
            double const ydist = -OldPositions[element][1] + OldPositions[j][1];
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
        return n;
    }

}



int main() {

    //Implementing a fourth order Runge Kutta method
    std::ofstream vorticityfile;
    auto start = std::chrono::steady_clock::now();
    myfile.open("zxcddata.csv");
    vorticityfile.open("zxcdvorticity.csv");

    for(int i =0;i<n;i++){
        vorticityfile<<i;
        myfile << 2*i;
        myfile << ",";
        myfile << 2*i+1;
        if(i!=n-1){
            myfile << ",";
            vorticityfile <<",";
        }
    }


    myfile<<"\n";
    vorticityfile<<"\n";
    double Velocity[4][n][2]; //Velocity for each stage of runge kutta
    double Vorticity[n];
    double TempPos[n][2];
    double Position[n][2];
    const double halfh = double(0.5)*h;
    const double sixthh = double(1)/double(6)*h;
    double c = 0;

    for(int i =0;i<n/2;i++){
        Position[i][0] = dis(gen);
        Position[i][1] = dis(gen);
        if(Position[i][0]*Position[i][0]+Position[i][1]*Position[i][1]>rnought*rnought){
            i--;
        }
        else{
            Vorticity[i] = wnought*exp(-(Position[i][0]*Position[i][0]+Position[i][1]*Position[i][1])/pow(rnought,2));
            vorticityfile<< Vorticity[i];
            vorticityfile<<",";
            c= c + Vorticity[i];
        }
    }

    for(int i = n/2;i<n;i++){
        Position[i][0] = dis(gen);
        Position[i][1] = dis(gen);
        if(Position[i][0]*Position[i][0]+Position[i][1]*Position[i][1]>rnought*rnought){
            i--;
        }
        else{
            Vorticity[i] = wnought*exp(-(Position[i][0]*Position[i][0]+Position[i][1]*Position[i][1])/pow(rnought,2));
            vorticityfile<<Vorticity[i];
            if(i!=n-1){

                vorticityfile <<",";
            }

            c= c + Vorticity[i];
            Position[i][0] = Position[i][0]+2.0   ;

        }

    }


    /*
    for(int i =0;i<n;i++){
        Position[i][0] = dis(gen);
        Position[i][1] = dis(gen);
        if(Position[i][0]*Position[i][0]+Position[i][1]*Position[i][1]>1.0){
            i--;
        }
        else{
            Vorticity[i] = wnought*exp(-(Position[i][0]*Position[i][0]+Position[i][1]*Position[i][1]));
            c= c + Vorticity[i];
            vorticityfile<<Vorticity[i];
            if(i!=n-1){

                vorticityfile <<",";
            }
        }

    }
    */
    /*

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
    std::cout<<"hi"<<a<<"\n";
    std::cout<<overrnought;
    overrnought = -double(1)/pow(1.1*a,2);
    std::cout<<overrnought;
    */
    /*
    Position[0][0] = -1.0;
    Position[0][1] = -1.0;
    Position[1][0] = 1.0;
    Position[1][1] = 1.0;

    Vorticity[0] = 1.0;
    Vorticity[1] = 1.0;

    vorticityfile<<1.0;
    vorticityfile<<",";
    vorticityfile<<1.0;
    */

    vorticityfile.close();
    std::cout<<"vorticity"<<c<<"\n";
    double initialhamiltonian = hamiltonian(Position,Vorticity);
    double initialAngularImpulse = AngularImpulse(Position,Vorticity);
    double initialLinearImpulse = LinearImpulse(Position,Vorticity);
    double initialpos = pow(Position[0][0]-Position[1][0],2)+pow(Position[0][1]-Position[1][1],2);
    int counter = 0;

    for(int step =1;step<StepCount;step++){
        counter++;

        if(step%(StepCount/100)==0)std::cout<<double(step)/double(StepCount)*100.0<<"%"<<"\n";

        if(1==1) {
            for (int i = 0; i < n; i++) {
                myfile << Position[i][0];
                myfile << ",";
                myfile << Position[i][1];
                if (i != n - 1) {
                    myfile << ",";
                }

            }
            counter = 0;
        }

        myfile << "\n";


        ParalellZeroOrder(Position,Velocity[0],Vorticity);
        for(int i =0;i<LargestCount;i++){
            i = LargestMethod(Position,Velocity[0],Vorticity,i);

        }

        for(int i = 0;i<n;i++){

            TempPos[i][0] = Position[i][0]+halfh*Velocity[0][i][0];
            TempPos[i][1] = Position[i][1]+halfh*Velocity[0][i][1];

        }

        ParalellZeroOrder(TempPos,Velocity[1],Vorticity);
        for(int i =0;i<LargestCount;i++){
            i =LargestMethod(TempPos,Velocity[1],Vorticity,i);
        }

        for(int i = 0;i<-n;i++){

            TempPos[i][0] = Position[i][0]+halfh*Velocity[1][i][0];
            TempPos[i][1] = Position[i][1]+halfh*Velocity[1][i][1];

        }

        ParalellZeroOrder(TempPos,Velocity[2],Vorticity);
        for(int i =0;i<LargestCount;i++){
            i =LargestMethod(TempPos,Velocity[2],Vorticity,i);
        }

        for(int i = 0;i<n;i++){

            TempPos[i][0] = Position[i][0]+h*Velocity[2][i][0];
            TempPos[i][1] = Position[i][1]+h*Velocity[2][i][1];

        }

        ParalellZeroOrder(TempPos,Velocity[3],Vorticity);
        for(int i =0;i<LargestCount;i++){
            i =LargestMethod(TempPos,Velocity[3],Vorticity,i);
        }

        for(int i =0;i<n;i++){
            //Add all at once so data is only pulled to the cpu once instead of repeatedly after each time
            Position[i][0] = Position[i][0] + sixthh*(Velocity[0][i][0]+Velocity[1][i][0]*double(2)+Velocity[2][i][0]*double(2)+Velocity[3][i][0]);
            Position[i][1] = Position[i][1] + sixthh*(Velocity[0][i][1]+Velocity[1][i][1]*double(2)+Velocity[2][i][1]*double(2)+Velocity[3][i][1]);

        }


    }
    double x1 = 0;
    double y1 = 0;
    for(int i =0;i<n/2;i++){
        x1 += Position[i][0];
        y1 +=Position[i][1];
    }
    x1 = x1/(n/2.0);
    y1 = y1/(n/2.0);
    double x2 = 0;
    double y2 = 0;
    for(int i =n/2;i<n;i++){
        x2 += Position[i][0];
        y2 +=Position[i][1];
    }
    x2 = x2/(n/2.0);
    y2 = y2/(n/2.0);
    std::cout<<"radius: "<<sqrt(pow(x1-x2,2)+pow(y1-y2,2))<<"\n";


    if(n==2)std::cout<<"PositionConservation: " <<1 - initialpos/(pow(Position[0][0]-Position[1][0],2)+pow(Position[0][1]-Position[1][1],2))<<"\n";

    initialhamiltonian = 1.0-hamiltonian(Position,Vorticity)/initialhamiltonian;
    initialAngularImpulse = 1.0-AngularImpulse(Position,Vorticity)/initialAngularImpulse;
    initialLinearImpulse = 1.0-LinearImpulse(Position,Vorticity)/initialLinearImpulse;
    std::cout<<"HamiltonianConservation: "<<initialhamiltonian<<"\n";
    std::cout<<"AngularConservation: "<<initialAngularImpulse<<"\n";
    std::cout<<"LinearConservation: "<<initialLinearImpulse<<"\n";
    myfile.close();


    myfile.open("zxcddata2.csv");

    myfile<<h<<","<<initialhamiltonian<<","<<initialAngularImpulse<<","<<initialLinearImpulse<<","<<n;
    myfile.close();


    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
    return 0;

}
/*

 HamiltonianConservation: 7.09988
AngularConservation: 0.271781
LinearConservation: 1

 */