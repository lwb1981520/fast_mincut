




typedef unsigned long long int excessType; 


typedef unsigned long cType;
typedef unsigned int sType;

typedef 
   struct edgeProp
{

   cType endNode;
   cType cap;
   cType w; 
   cType avgCV; 
   long tmp;

   struct edgeProp* rev; 

}edgeP;


typedef  
   struct nodeProp
{
   edgeP* edges;
   cType maxEdges;
   cType nIdx;
   cType totalCap;
   
   sType* orderedEdges;

} nodeP;

typedef
struct NodePropExtra_
{
   cType fa;
   cType dep;
   long cv;
   short s;


} NodePropExtra;


typedef
struct NodePropArr_{

   cType * pfa;
   cType * pdep;
   long * pcv;
   short * ps;
   cType * pacc_upid; 
   long* pacc_upmincv; 
   cType * pacc_pos_upmincv; 
   cType * pacc_jointid; 
   long* pacc_jointmcv; 
   
} NodePropArr;






typedef
struct GraphData_{

  long N, M;
  nodeP *nodes;

} GraphData;



typedef
struct RandomData_{
  cType *randNums;
  cType randNumIdx;
  cType len;
  cType maxLen;
} RandomData;



typedef
struct PreprocData_{

  
  GraphData *gd;

  
  NodePropArr* allResults;

  
  RandomData* rd;

  
  cType *gpfa;
  cType *gpdep;
  long *gpcv;
  short *gps;
  cType *gpaccup;
  long *gpaccmcv;
  cType *gpaccposmcv;
  cType *gpaccjointid;
  long *gpaccjointmcv;
  

  cType *roots; 

  int mode; 

  int P; 
  int total; 
  int SPAN_LEN; 

  long C ; 
  cType* gcut_fp; 
  long* gcut_fp_nh; 
  cType* garr; 
  long gh;


} PreprocData;
