/* defs.h */

//#ifdef EXCESS_TYPE_LONG
///typedef unsigned long excessType;
//#else
typedef unsigned long long int excessType; /* change to double if not supported */
//#endif

typedef unsigned long cType;
typedef unsigned int sType;

typedef 
   struct edgeProp
{

   cType endNode;
   cType cap;
   cType w; /*weight*/
   cType avgCV; /*average CV of cut set where this edge belongs to */
   long tmp;

   struct edgeProp* rev; //reverse

}edgeP;


typedef  /* node */
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
   cType * pacc_upid; // up node id
   long* pacc_upmincv; // min of [currnet node, upnode)
   cType * pacc_pos_upmincv; // the diff of depth of nearest min value
   cType * pacc_jointid; //the nearest joint id;
   long* pacc_jointmcv; //the min of [cur, joint node)
   
} NodePropArr;

/////////////////////////////////////////
///////////////////////////////////////// The definition of data structure
/////////////////////////////////////////

//data structure that holds graph data
typedef
struct GraphData_{

  long N, M;
  nodeP *nodes;

} GraphData;


//data structure that holds data for randomrization
typedef
struct RandomData_{
  cType *randNums;
  cType randNumIdx;
  cType len;
  cType maxLen;
} RandomData;


//data structure that holds preprocessing data
typedef
struct PreprocData_{

  //graph data
  GraphData *gd;

  //holds data in all passes
  NodePropArr* allResults;

  //pre-generated data
  RandomData* rd;

  //<BEING> hold the hot data in current pass
  cType *gpfa;
  cType *gpdep;
  long *gpcv;
  short *gps;
  cType *gpaccup;
  long *gpaccmcv;
  cType *gpaccposmcv;
  cType *gpaccjointid;
  long *gpaccjointmcv;
  //<END>

  cType *roots; //records root in each pass

  int mode; //the traversing mode, as explained in the paper

  int P; // P% percent of total passes in mode 1, the remaining in mode 2
  int total; //number of total passes
  int SPAN_LEN; //the length of a ancestor span, used for acceleration in traversal trees

  long C ; //the maximum index for per-node cut fingerprinter
  cType* gcut_fp; //cut fingerprint for all nodes
  long* gcut_fp_nh; //first index of per-node cut fingerprint
  cType* garr; //临时变量，每次要清零
  long gh;


} PreprocData;
