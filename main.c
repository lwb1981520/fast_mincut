#include "hi_treem.c"

int main(argc, argv)

    int argc;
char *argv[];

{
    PreprocData *pd = walloc(1,sizeof(PreprocData));
    pd->P = 20; 
    pd->total = 100; 
    pd->C = 10000;

    loadGraphData(pd); 

    initPreprocData(pd); 

struct timespec time_start={0,0},time_end={0,0};
clock_gettime(CLOCK_REALTIME,&time_start);
    preProc(pd); 

    preProc_cut_fingerprint(pd);

    cumsum_cut_fingerprint(pd);

    clock_gettime(CLOCK_REALTIME,&time_end);
	
    double tm = 10e9*time_end.tv_sec +time_end.tv_nsec - 10e9*time_start.tv_sec - time_start.tv_nsec;
    tm = tm/10e9;
	printf("c total proc time %f \n", tm);
	
    calcuRandomPairs(100,pd); 
    
  
    exit(0);

}
