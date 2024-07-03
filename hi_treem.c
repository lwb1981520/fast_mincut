

  
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <values.h>
#include <math.h>
#include <time.h>

#include "types_treem.h"  
#include "parser_treem.c" 
#include "timer.c"        

#define MAX_LONG LONG_MAX

#define min(s, t) ((s) < (t) ? (s) : (t))
#define max(s, t) ((s) > (t) ? (s) : (t))











void *walloc(unsigned int num, unsigned int size)
{
  void *ptr = calloc(num, size);
  assert(ptr != NULL);
  return ptr;
}
























RandomData* initrand2(RandomData *rd )
{
  srand((int)(timer() * 1000));

  rd->randNumIdx = 0;

  for (int i = 0; i < rd->maxLen; i++)
  {
    rd->randNums[i] = (cType)rand();
  }
  
  return rd;

}

RandomData* initrand(cType len)
{
  RandomData *rd = walloc(1, sizeof(RandomData));
  rd->maxLen = len;
  rd->randNums = (cType *)walloc(rd->maxLen, sizeof(cType));
  
  initrand2(rd);

  return rd;
}


cType mrand(RandomData* rd)
{
  if(rd->randNumIdx >= rd->maxLen){
		initrand2(rd);
  }	  
  return rd->randNums[rd->randNumIdx++];
}


void HeapAdjustDown(sType *idx, edgeP * edges ,int start,int end)  
{  
    sType tempIdx = idx[start];  

    int i = 2*start+1;      
    
    
    
    while(i<=end)  
    {  
        if(i+1<=end && edges[idx[i+1]].tmp > edges[idx[i]].tmp )    
            i++;  

        if(edges[idx[i]].tmp <= edges[tempIdx].tmp )   
            break;  

        idx[start] = idx[i];

        start = i;  
        i = 2*start+1;  
    }  

    idx[start] = tempIdx;  
}  
  
void HeapSort(sType *idx, edgeP * edges, int len)  
{  

    int i;  
    for(i=(len-1)/2;i>=0;i--){  
        HeapAdjustDown(idx,edges,i,len-1);  
    }

    for(i=len-1;i>0;i--)
    {  
        
        sType temp = idx[i];  
        idx[i] = idx[0];  
        idx[0] = temp;  

        HeapAdjustDown(idx,edges,0,i-1);  
    }  

}  
  

void deOrderEdgeByRandomCap(nodeP *np,PreprocData *pd)
{

  cType cnt = np->nIdx;

  sType *idxs = np->orderedEdges;
  edgeP *pedges = np->edges;

  for (int i = 0; i < cnt; i++)
  {
    pedges[i].tmp = -1*mrand(pd->rd) % pedges[i].cap;
  }

    assert(cnt<4 || idxs[2]!=idxs[3]);
    HeapSort(idxs,pedges,cnt);
    assert(cnt<4 || idxs[2]!=idxs[3]);  
}


void aOrderEdgeByAvgCV(nodeP *np,PreprocData *pd)
{
  if (np->nIdx == 0)
  {
    return;
  }
  int cnt = np->nIdx;

  sType *idxs = np->orderedEdges;
  edgeP *pedges = np->edges;

  for (int i = 0; i < cnt; i++)
  {
    edgeP *pp = pedges +i;
    long acv = pp->avgCV;

    if(acv == 0 ){
      pp->tmp = MAX_LONG;
    }
    else{
      pedges[i].tmp = mrand(pd->rd) % acv; 
    }
  }

    HeapSort(idxs,pedges,cnt);
}


/*
  the function to add more data to a traversal tree to accelerate the searching in the tree
  The idea is to precalcuate minimal cv value of a span of nodes in a traversal tree, e.g., when SPAN_LEN = 100, and a node has depth of 200, then the algorithm will pre-calculate the minimal cv value of the nodes between the node (dep=200) and an acestor(dep=101)
  upid is the id of the last SPAN node, mcv is the min cv among all previous nodes in the recent SPAN
  lastDepMCV is the depth of the node depth that has the minimal cv in the span
  lastJointNodeId is the last ancestor node id that has more than one child nodes
  lastJointMCV is the cv of lastJoineNodeId
*/
void buildAcc(PreprocData *pd, cType curN, cType upid, long mcv, cType lastDepMCV,cType lastJointNodeId, long lastJointMCV)
{
  nodeP* nodes = pd->gd->nodes;
  assert(curN <= pd->gd->N && curN >= 1);
  int cnt = (nodes + curN)->nIdx;
  edgeP *pedges = (nodes + curN)->edges;

  long curCV = pd->gpcv[curN];
  cType curDep = pd->gpdep[curN];

  assert(curDep == 0 || curCV >0);

  mcv = min(mcv, curCV);
  if(mcv == curCV){
    lastDepMCV = curDep;
  }

  pd->gpaccmcv[curN] = mcv;
  pd->gpaccup[curN] = upid;
  pd->gpaccposmcv[curN] = curDep - lastDepMCV;

  if (curDep % pd->SPAN_LEN == 0)
  {
    upid = curN;
    mcv = MAX_LONG;
    lastDepMCV = 0; 
  }


  cType childCnt = pd->gpaccjointid[curN];
  lastJointMCV = min(lastJointMCV,curCV);
  pd->gpaccjointmcv[curN] = lastJointMCV;
  pd->gpaccjointid[curN] = lastJointNodeId;

  assert(pd->gpdep[curN] == 0 || pd->gpaccjointmcv[curN] > 0.1);
  
  
  if(childCnt > 1 ){
    lastJointNodeId = curN;
    lastJointMCV = MAX_LONG;
  }
  else{
    
  }


  while (cnt > 0)
  {
    if (pd->gpfa[pedges->endNode] == curN)
    {
      buildAcc(pd,pedges->endNode, upid, mcv,lastDepMCV,lastJointNodeId,lastJointMCV);
    }
    pedges++;
    cnt--;
  }

  
}


cType gRoot = 0;
cType gRDep = 0;
void markCut(cType curN, PreprocData *pd)
{

  nodeP* nodes = pd->gd->nodes;
  assert(curN <= pd->gd->N && curN >= 1);

  assert((nodes + curN)->nIdx > 0);

  short *curS = pd->gps + curN;

  assert(*curS == 0);

  long *curCV = pd->gpcv + curN;
  cType *curDep = pd->gpdep + curN;

  if(gRDep < *curCV){
    gRDep = *curCV;
    gRoot = curN;
  }
  
  *curS = 1;
  *curCV = 0;
  nodeP *np = nodes + curN;
  edgeP *pedges = np->edges;
  int cnt = np->nIdx;

  if(pd->gpdep[curN] > 0){
    
    
    
    

    assert(pd->gpdep[pd->gpfa[curN]] == pd->gpdep[curN] -1);

  }

  if (pedges == NULL)
  {
    *curS = 2;
    return;
  }

  if (np->orderedEdges == NULL)
  {
    np->orderedEdges = (sType *)walloc(cnt + 1, sizeof(sType));
    for (int i = 0; i < cnt; i++)
    {
      np->orderedEdges[i] = i;
    }
  }

  long cap;
  sType *idxs = np->orderedEdges;

  if (pd->mode == 1)
  {
    deOrderEdgeByRandomCap(np,pd); 
  }
  else if (pd->mode == 2)
  {
    aOrderEdgeByAvgCV(np,pd);
  }
  else if(pd->mode == 11){
    aOrderEdgeByAvgCV(np,pd);
  }

  np->totalCap = 0;
  for (int ni = 0; ni < cnt; ni++)
  {
    

    edgeP *eh = pedges + idxs[ni];
    cType zn = eh->endNode;
    np->totalCap  += eh->cap;

    assert(zn != 0);
    assert(zn != curN);

    short zs = pd->gps[zn];

    assert(!(pd->gpfa[zn] == curN && zs == 2));


    if (zs == 1)
    {
      cap = eh->cap;
      *curCV += cap;
      pd->gpcv[zn] -= cap;
    }
    else if (zs == 0)
    {
      pd->gpfa[zn] = curN;
      pd->gpdep[zn] = *curDep + 1;
      
      markCut(zn,pd);
      assert(pd->gpdep[zn] == pd->gpdep[curN] + 1);
      *curCV += pd->gpcv[zn];
    }
    else
    {
      
    }
  }


  if(pd->mode == 1){

    for (int ni = 0; ni < cnt; ni++)
    {
      

      edgeP *eh = pedges + idxs[ni];
      cType zn = eh->endNode;
      

      assert(zn != 0);
      assert(zn != curN);
      short zs = pd->gps[zn];
    

      
      if (zs == 1 && pd->gpdep[zn] != *curDep - 1)
      {
          cType weight = eh-> w;
          if(eh->avgCV == 0){
            eh->avgCV = MAX_LONG;
          }
          eh->avgCV = min(eh->avgCV, *curCV);
          eh->w = weight+1;
          
          edgeP *reh = eh->rev;

          if(reh->avgCV == 0){
            reh->avgCV = MAX_LONG;
          }

          weight = reh-> w;
          reh->avgCV = min(reh->avgCV, *curCV);
          reh->w = weight+1;

      }
      
    }
  }

    
  assert(pd->gpdep[curN] == 0 || *curCV >0);

  *curS = 2;

}


long solveMaxFlowAccVER4(long minCandi, NodePropArr np, cType s, cType t, int SPAN_LEN)
{
  cType *pDep = np.pdep;
  long *pCV = np.pcv;
  cType *pFa = np.pfa;
  cType *paccup = np.pacc_upid;
  long *paccmcv = np.pacc_upmincv;
  cType *paccposmcv = np.pacc_pos_upmincv;
  cType *jup = np.pacc_jointid;
  long *jmcv = np.pacc_jointmcv;

  assert(s != t);
  if (pDep[s] < pDep[t])
  {
    cType tmp = s;
    s = t;
    t = tmp;
  }

  cType depT = pDep[t];
  assert(pDep[s] >= depT);

  long mcv = MAX_LONG;
  cType ups = paccup[s];

  while (pDep[ups] > depT)
  {
    assert(pDep[ups] % SPAN_LEN == 0);
    mcv = min(mcv, paccmcv[s]);
    s = ups;
    ups = paccup[ups];
  } 

  assert(pDep[ups] <=depT && pDep[s] >= depT);

  cType upt = t;

  if (pDep[t] % SPAN_LEN != 0)
  {
    upt = paccup[t];
  }

  assert(pDep[ups] == pDep[upt]);

  while (ups != upt)
  {
    mcv = min(mcv, min(paccmcv[s], paccmcv[t]));

    s = ups;
    ups = paccup[ups];

    t = upt;
    upt = paccup[upt];

    assert(pDep[s] % SPAN_LEN == 0);
    assert(pDep[t] == pDep[s]);
  }        
  
  assert(ups == upt);
  if(s == t){
    return mcv;
  }
  
  cType min_bound2 = min(paccmcv[s],paccmcv[t]);

  if(min_bound2 >= mcv || min_bound2 >= minCandi){
    
    return mcv;
  }

  if(min_bound2 == paccmcv[s] && pDep[t]+paccposmcv[s] < pDep[s]){
    return mcv;
  }


  
  
  

  if(pDep[s] != pDep[t]){

    if(pDep[s] < pDep[t]){
      cType tmp = s;
      s = t;
      t = tmp;
    }  

    cType jups = jup[s];
    cType depT = pDep[t];
    while(pDep[jups] > depT){
      mcv = min(mcv, jmcv[s]);
      if(mcv == min_bound2){
        return mcv;
      }

      s = jups;
      jups = jup[s];

    }

    if(pDep[jups] == depT){
      if(jups == t){
        return min(mcv,jmcv[s]);
      }
      else{
        
        assert(pDep[jups] == depT && pDep[s] > depT);
        goto STEP_CHECK_IN_TWOLINES;
      }
    }
    else{
      assert(pDep[jups] < depT && pDep[s] > depT);
      while(pDep[s] > depT){
        mcv = min(mcv,pCV[s]);
        if(mcv == min_bound2){
          return mcv;
        }
        s = pFa[s];        
      }

      if(s == t){
        return mcv;
      }

      assert(s!=t && pDep[s] == depT);
      
      goto STEP_CHECK_IN_TWOLINES;

    }
  }
  else{
      
      if(s == t){
        return mcv;
      }    
      goto STEP_CHECK_IN_TWOLINES;
  }

STEP_CHECK_IN_TWOLINES:

  while (s != t)
  { 
    
    
    if(pDep[s] > pDep[t]){
      mcv = min(mcv, jmcv[s]);
      if(mcv == min_bound2){
        return mcv;
      }      
      s = jup[s];
    }
    else if(pDep[s] < pDep[t]){
      mcv = min(mcv, jmcv[t]);
      if(mcv == min_bound2){
        return mcv;
      }      
      t = jup[t];
    }
    else{
      mcv = min(mcv, min(jmcv[s],jmcv[t]));
      if(mcv == min_bound2){
        return mcv;
      }         
      s = jup[s];
      t = jup[t];
    }

  }

  return mcv;

}

long solveMaxFlowAccVER3(long minCandi, NodePropArr np, cType s, cType t)
{
  cType *pDep = np.pdep;
  long *pCV = np.pcv;
  cType *pFa = np.pfa;


  assert(s != t);
  
  if (pDep[s] < pDep[t])
  {
    cType tmp = s;
    s = t;
    t = tmp;
  }

  cType depT = pDep[t];
  assert(pDep[s] >= depT);

  

  long mcv = minCandi;

  while(pDep[s] > depT){
    
    
    
    
    assert(pDep[pFa[s]] == pDep[s] -1);
    mcv = min(mcv,pCV[s]); 
    s = pFa[s];  
    
  }

  

  if(s == t){
    return mcv;
  }
  
  
  
  assert(pDep[s] == depT);

  while (s != t)
  { 

    mcv = min(mcv, min(pCV[s], pCV[t]));
    s = pFa[s];
    t = pFa[t];

  }

  return mcv;

}


void loadGraphData(PreprocData *pd){
  pd->gd = walloc(1,sizeof(GraphData));  

  parse(&(pd->gd->N), &(pd->gd->M), &(pd->gd->nodes));

  printf("c nodes:       %10ld\nc arcs:        %10ld\nc\n", pd->gd->N, pd->gd->M);
}


void initPreprocData(PreprocData *pd){
  pd->rd = NULL;
  pd->SPAN_LEN = (int)(sqrt(pd->gd->N));

  pd->roots = (cType *)walloc(pd->total + 2, sizeof(cType));
  pd->allResults = walloc(pd->total+2, sizeof(NodePropArr));

  NodePropArr * allResults = pd->allResults;
  cType len = pd->gd->N + 2;
  
  pd->gcut_fp = (cType *)walloc(len* pd->C, sizeof(cType));
  pd->gcut_fp_nh = (long *)walloc(len, sizeof(long));
  memset(pd->gcut_fp, 0, len* pd->C * sizeof(cType));
  memset(pd->gcut_fp_nh, 0, len * sizeof(long));

  for(cType curN=1; curN<=pd->gd->N; curN++){
    pd->gcut_fp_nh[curN] = MAX_LONG;
  }  

  pd->garr = (cType *)walloc(100* pd->C, sizeof(cType));
  pd->gh = MAX_LONG;
  memset(pd->garr, 0, 100* pd->C * sizeof(long));    

  for (int i = 0; i < pd->total; i++)
  {
    allResults[i].pfa = (cType *)walloc(len, sizeof(cType));
    allResults[i].pdep = (cType *)walloc(len, sizeof(cType));
    allResults[i].pcv = (long *)walloc(len, sizeof(long));
    allResults[i].ps = (short *)walloc(len, sizeof(short));
    
    
    
    
    

    memset(allResults[i].pfa, 0, len * sizeof(cType));
    memset(allResults[i].pdep, 0, len * sizeof(cType));
    memset(allResults[i].pcv, 0, len * sizeof(long));
    memset(allResults[i].ps, 0, len * sizeof(short));
    
    
    
    
    
  }  
}



void preProc(PreprocData *pd){
  double tm;
  double totalProcTime = 0;

  struct timespec time_start={0,0},time_end={0,0};

  NodePropArr *allResults = pd->allResults;
  
  cType root;
  gRoot = 1;

  for (int ipass = 0; ipass < pd->total; ipass++)
  {
    if(pd->rd != NULL){
      free(pd->rd);
      pd->rd = NULL;
    }
    pd->rd = initrand(pd->gd->M*2+7);    
    
    pd->gpfa = allResults[ipass].pfa;
    pd->gpdep = allResults[ipass].pdep;
    pd->gpcv = allResults[ipass].pcv;
    pd->gps = allResults[ipass].ps;
    
    
    
    
    

    if(pd->total == 110){
      pd-> mode = 11; 
    }
    else{
      pd->mode = ipass < pd->P * pd->total / 100 ? 1 : 2;
    }

    if(ipass % 4 == 0){
      gRoot = 1 +  ((mrand(pd->rd)) % pd->gd->N);
    }

    root = gRoot;
    pd->roots[ipass] = root;
    gRDep = 0;
    pd->gpdep[root] = 0;
    
    printf("pass %d, randidx %ld, root is %ld\n",ipass, 0,root);
    
    clock_gettime(CLOCK_REALTIME,&time_start);
    markCut(root,pd);


    pd->gpcv[root] = MAX_LONG;

    
    clock_gettime(CLOCK_REALTIME,&time_end);
    tm = 10e9*time_end.tv_sec +time_end.tv_nsec - 10e9*time_start.tv_sec - time_start.tv_nsec;
    tm = tm/10e9;
    totalProcTime += tm;
    printf("c proctime for onepass: %10.06f\n",tm);
    if (ipass % 10 == 0)
    {
      printf("c the %d passes\n", ipass);
    }

    
    free(pd->gps);
  }

  printf("c preprocess times %10.6f\n", totalProcTime);

}


void cut_fingerprint(cType curN, PreprocData *pd)
{
  nodeP* nodes = pd->gd->nodes;
  long *curCV = pd->gpcv + curN;
  
  
  nodeP *np = nodes + curN;
  edgeP *pedges = np->edges;
  int cnt = np->nIdx;
  
    
  
  
  cType ov;
  long oh;

  long ulimit = min(np->totalCap, pd->C);

  if(*curCV <= ulimit){
    ov = pd->garr[*curCV]; oh = pd->gh;
    pd->garr[*curCV] = curN; pd->gh = min(pd->gh, *curCV);
  }

  if(pd->gh <= ulimit){ 
    
    
    for( long i=pd->gh; i<=ulimit; i++){
      pd->gcut_fp[curN*pd->C+i] += pd->garr[i];

      
      /*
        叠加计算的含义是：单次序列不是每个割值的个数，而是每个割值及以前割值的总个数
        动机：让某个值不一样时能向前传递
        比较方法：
          认为相同的idx处，之前都没有不同，需要向后搜索
          
        潜在问题1：如果garr叠加序列，idx=3不一样，但到了idx=5又一样了咋办？ 就会出现误判，等于两个地方都不一样，但是出现加和后一样了

      */
      
      
    }


	
  }
  
    pd->gcut_fp_nh[curN] = min (pd->gcut_fp_nh[curN], pd->gh);  
    pd->gcut_fp_nh[curN] = min (pd->gcut_fp_nh[curN], ulimit);  

  for (int ni = 0; ni < cnt; ni++)
  {
    

    edgeP *eh = pedges + ni;
    cType zn = eh->endNode;
    if(pd->gpfa[zn] == curN){
      cut_fingerprint(zn,pd);
    }
  }

  if(*curCV <= ulimit){
    pd->garr[*curCV] = ov; pd->gh = oh;
  }

  assert (pd->gcut_fp_nh[curN] > 0 && pd->gcut_fp_nh[curN] != MAX_LONG);
  
}

int isCumSum = 0;
void cumsum_cut_fingerprint(PreprocData *pd)
{
  isCumSum = 1;
  nodeP* nodes = pd->gd->nodes;
  for(cType curN=1; curN <= pd->gd->N; curN++){

    long ulimit = min((nodes+curN)->totalCap, pd->C);

      
      long sum = 0;
      for( long i=pd->gcut_fp_nh[curN]; i<=ulimit; i++){
        sum+=pd->gcut_fp[curN*pd->C+i];
        pd->gcut_fp[curN*pd->C+i] = sum;

      }

      

  }

}


void preProc_cut_fingerprint(PreprocData *pd){
  double tm;
  double totalProcTime = 0;
  struct timespec time_start={0,0},time_end={0,0};
  NodePropArr *allResults = pd->allResults;
  
  
 
  

  for (int ipass = 0; ipass < pd->total; ipass++)
  {
    if(pd->rd != NULL){
      free(pd->rd);
      pd->rd = NULL;
    }
    pd->rd = initrand(pd->gd->M*2);    
    
    pd->gpfa = allResults[ipass].pfa;
    
    pd->gpcv = allResults[ipass].pcv;
    
    
    
    
    
    

    
    
    
    
    
    
    

    
    
    printf("pass %d, root is %ld\n",ipass, pd->roots[ipass]);
    
    clock_gettime(CLOCK_REALTIME,&time_start);
    pd->gh = MAX_LONG;
    memset(pd->garr, 0, 100* pd->C * sizeof(long));   
    cut_fingerprint(pd->roots[ipass],pd);
    clock_gettime(CLOCK_REALTIME,&time_end);

    free(pd->gpfa);
    
    free(pd->gpcv);
    

    tm = 10e9*time_end.tv_sec +time_end.tv_nsec - 10e9*time_start.tv_sec - time_start.tv_nsec;
    tm = tm/10e9; 
    totalProcTime += tm;
    printf("c proctime for cut fingerprint onepass: %10.06f\n", tm);
    if (ipass % 10 == 0)
    {
      printf("c the %d passes\n", ipass);
    }
  }

  printf("c cut fingerprint preprocess times %10.6f\n", totalProcTime);

}


long solveMaxFlowAccVER4_CF(cType s, cType t, PreprocData *pd)
{
  assert(pd->gcut_fp_nh[s] > 0);
  
  long mv  = MAX_LONG; 
  printf("result2: s %ld t %ld\n",pd->gcut_fp_nh[s] ,pd->gcut_fp_nh[t] );
  if(pd->gcut_fp_nh[s] != pd->gcut_fp_nh[t]){
    mv = min(pd->gcut_fp_nh[s], mv);
    mv = min(pd->gcut_fp_nh[t], mv);
    return mv;
  }  

  for(long i=pd->gcut_fp_nh[s]; i<mv; i++){
    if(pd->gcut_fp[s*pd->C+i] != pd->gcut_fp[t*pd->C+i]){
      printf("@@@@@@ count %d\n",1+i - pd->gcut_fp_nh[s]);
      return i;
    }
  }

  return mv; 
}


long solveMaxFlowAccVER5_CF(cType s, cType t, PreprocData *pd)
{

  int count = 1;  

  assert(pd->gcut_fp_nh[s] > 0);
  
  long mv  = min((pd->gd->nodes+s)->totalCap,(pd->gd->nodes+t)->totalCap);

  

  
  if(pd->gcut_fp_nh[s] != pd->gcut_fp_nh[t] ){
    mv = min(pd->gcut_fp_nh[s], mv);
    mv = min(pd->gcut_fp_nh[t], mv);
    
    return mv;
  }  

  
  
  long lidx = pd->gcut_fp_nh[s];
  long hidx = mv;

  

  
  
  
  

  while(1==1){
	if(lidx > hidx){
		printf("c s is %ld, hidx %ld - lidx %ld =  %ld \n",s, hidx, lidx, hidx - lidx);
		printf("result2: s %d %ld t %d %ld\n",s, pd->gcut_fp_nh[s] , t, pd->gcut_fp_nh[t] );
	}

    assert(lidx < hidx);
	if(lidx == hidx){
		return hidx;
	}
	
    long ix = (lidx + hidx)/2;
    if(pd->gcut_fp[s*pd->C+ix] != pd->gcut_fp[t*pd->C+ix]){
      
      hidx = ix;
    }
    else{
      
      lidx = ix;
    }

    
    

    if(lidx+1 == hidx){
      printf("@count : %ld final\n",count);
      return hidx;
    }

    count++;

  }

}



void calcuRandomPairs(int numOfPairs, PreprocData *pd){
  
  double totalTime = 0;
  long mv = MAX_LONG, mv2 = MAX_LONG;

  double curTime = 0;
  cType ns, nt;

  if(pd->rd != NULL){
    free(pd->rd);
    pd->rd = NULL;
  }
  pd->rd = initrand(pd->gd->M*2);  
	struct timespec time_start={0,0},time_end={0,0};
  long countNotTcap = 0;
  if(isCumSum == 0){
    printf("use ver4 \n");
  }
  else{
    printf("use Ver5 \n");
  }

  for (int ipair = 0; ipair < numOfPairs;)
  {

    
    ns = 1 + ((mrand(pd->rd) * mrand(pd->rd)) % (pd->gd->N));
    nt = 1 + ((mrand(pd->rd) * mrand(pd->rd)) % (pd->gd->N));
    if (ns != nt)
    {
      printf("tcap : %ld %ld\n",(pd->gd->nodes+ns)->totalCap,(pd->gd->nodes+nt)->totalCap);
      mv = min((pd->gd->nodes+ns)->totalCap,(pd->gd->nodes+nt)->totalCap);
      clock_gettime(CLOCK_REALTIME,&time_start);
      if(isCumSum == 0){
        mv = solveMaxFlowAccVER4_CF(ns,nt,pd);
      }
      else{
        mv = solveMaxFlowAccVER5_CF(ns,nt,pd);
        printf("~~result :VER4 %ld, VER5 %ld\n",mv,mv2);
      }
      
      
      
      
        
        
        
        
        
        
        
        
      
      	 
      clock_gettime(CLOCK_REALTIME,&time_end);
	    curTime = 10e9*time_end.tv_sec +time_end.tv_nsec - 10e9*time_start.tv_sec - time_start.tv_nsec;
	    curTime = curTime/10e9;
      totalTime += curTime;
      ipair++;
      printf("c hi_treem_res(n,s,mflow,tm) %lu %lu %12.01f %12.10f\n", ns, nt, 1.0 * mv, curTime);
      printf("mv %ld tcap %ld\n",mv,min ((pd->gd->nodes+ns)->totalCap,(pd->gd->nodes+nt)->totalCap));
      if(mv < min ((pd->gd->nodes+ns)->totalCap,(pd->gd->nodes+nt)->totalCap)){
        printf("      --mv %ld tcap %ld\n",mv,min ((pd->gd->nodes+ns)->totalCap,(pd->gd->nodes+nt)->totalCap));
        countNotTcap ++;
      }

    }
  }



  printf("c run ok! average time %10.6f countNotTcap %ld\n", totalTime / numOfPairs,countNotTcap);


}

long getTotalCap(PreprocData *pd, long curN){
  nodeP* nodes = pd->gd->nodes;
  nodeP *np = nodes + curN;
  edgeP *pedges = np->edges;
  int cnt = np->nIdx;
  long totalCap = 0;
  for(int i=0; i<cnt; i++){
    totalCap += pedges[i].cap;
  }

  return totalCap;
}




















































