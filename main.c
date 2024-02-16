#include "hi_treem.c"

int main(argc, argv)

    int argc;
char *argv[];

{
    PreprocData *pd = walloc(1,sizeof(PreprocData));
    pd->P = 20; // P% percent of total passes in mode 1, the remaining in mode 2
    pd->total = 100; //number of total passes
    pd->C = 10000;

    loadGraphData(pd); //load graph data from standard input

    initPreprocData(pd); //init data structure

struct timespec time_start={0,0},time_end={0,0};
clock_gettime(CLOCK_REALTIME,&time_start);
    preProc(pd); // preproc by traversing the graph for pd->total times

    preProc_cut_fingerprint(pd);

    cumsum_cut_fingerprint(pd);

    clock_gettime(CLOCK_REALTIME,&time_end);
	
    double tm = 10e9*time_end.tv_sec +time_end.tv_nsec - 10e9*time_start.tv_sec - time_start.tv_nsec;
    tm = tm/10e9;
	printf("c total proc time %f \n", tm);
	
    calcuRandomPairs(100,pd); // randomly choose 100 node pairs and calcu their min-cut and output
    //calcuRandomPairsNeighboring(1000,pd); // randomly choose 100 node pairs and calcu their min-cut and output
  
    exit(0);

}
