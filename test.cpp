#include "main.cpp"
#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>

int step;
long double steps, endtime;


int main(){
    System sys;
    sys.SetStartTime(0);
	
    sys.timestep = 0.096;
    sys.LoadFile("burn_test.start");
    sys.SetOutToFile("burn_test.txt",10);

    sys.initialize();

    endtime = 30000;
    steps = endtime/sys.timestep;

    struct timeval t1,t2;
    double timeuse;
    gettimeofday(&t1,NULL);

	sys.initialize_threads();

    for(step=0;step<steps;step++){
        sys.Output();
        sys.Step();
    }

    gettimeofday(&t2,NULL);
    timeuse = t2.tv_sec - t1.tv_sec + (t2.tv_usec - t1.tv_usec)/1000000.0;
    std::cout<<timeuse<<" s \n";
    printf("Use Time:%f\ns",timeuse);

}
