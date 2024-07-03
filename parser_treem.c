
void adjustEdgesIfNecessary(nodeP* np){
    if(np->edges == NULL || np->nIdx >= np->maxEdges){
    np->maxEdges += 5;
    edgeP* newep = (edgeP*) calloc(np->maxEdges,sizeof(edgeP));
    if(np->edges != NULL) { 
      memcpy(newep,np->edges,np->nIdx*sizeof(edgeP));
      for(int i=0; i<np->nIdx; i++){
        newep[i].rev->rev = newep+i;
      }
      free(np->edges); 
    }
    np->edges = newep;
  }
}




int parse( n_ad, m_ad, nodes_ad)


long    *n_ad;                 
long    *m_ad;                 
nodeP    **nodes_ad;            

{

#define MAXLINE       100000000	
#define ARC_FIELDS      3	
#define NODE_FIELDS     2	
#define P_FIELDS        3       
#define PROBLEM_TYPE "max"      































long    no_lines=0,             
        no_plines=0,            
        
        
        no_alines=0;            
        

char*    in_line = (char*) calloc ( MAXLINE, 1);    
char     pr_type[3];             
        

int     err_no;                 


#define EN1   0
#define EN2   1
#define EN3   2
#define EN4   3
#define EN6   4
#define EN10  5
#define EN7   6
#define EN8   7
#define EN9   8
#define EN11  9
#define EN12 10
#define EN13 11
#define EN14 12
#define EN16 13
#define EN15 14
#define EN17 15
#define EN18 16
#define EN21 17
#define EN19 18
#define EN20 19
#define EN22 20

static char *err_message[] = 
  { 
    "more than one problem line.",
    "wrong number of parameters in the problem line.",
    "it is not a Max Flow problem line.",
    "bad value of a parameter in the problem line.",
    "can't obtain enough memory to solve this problem.",
    "more than one line with the problem name.",
    "can't read problem name.",
    "problem description must be before  description.",
    "this parser doesn't support multiply sources and sinks.",
    "wrong number of parameters in the node line.",
    "wrong value of parameters in the node line.",
    " ",
    "source and sink descriptions must be before arc descriptions.",
    "too many arcs in the input.",
    "wrong number of parameters in the arc line.",
    "wrong value of parameters in the arc line.",
    "unknown line type in the input.",
    "reading error.",
    "not enough arcs in the input.",
    "source or sink doesn't have incident arcs.",
    "can't read anything from the input file."
  };


/* The main loop:
        -  reads the line of the input,
        -  analises its type,
        -  checks correctness of parameters,
        -  puts data to the arrays,
        -  does service functions
*/
long n,m,head,tail,cap;
while (fgets(in_line, MAXLINE, stdin) != NULL )
  {
  no_lines ++;
  if(no_lines % 10000 == 0){
    printf("c line %ld\n",no_lines);
  }

  switch (in_line[0])
    {
      case 'c':                  
      case '\n':                 
      case '\0':                 
                break;

      case 'p':                  
                if ( no_plines > 0 )
                   
                   { err_no = EN1 ; goto error; }

                no_plines = 1;
   
                if (
        
                    sscanf ( in_line, "%*c %3s %ld %ld", pr_type, &n, &m )
                != P_FIELDS
                   )
		    
		    { err_no = EN2; goto error; }

                if ( strcmp ( pr_type, PROBLEM_TYPE ) )
		    
		    { err_no = EN3; goto error; }

                if ( n <= 0  || m <= 0 )
		    
		    { err_no = EN4; goto error; }

        
                *n_ad = n;
                *m_ad = m;
                *nodes_ad =  (nodeP*) calloc ( n+2, sizeof(nodeP) );
                memset(*nodes_ad,0,(n+2)*sizeof(nodeP));

                break;

      case 'n':		         
		
    
    

    
		
 
		
    
    

		
    
    

		
		
		
		    
		
		
		

		
		
		

		

		
		
		

		
		
		

		
		
    
		
		


		break;

      case 'a':                    
		
    
    

		if ( no_alines >= m )
                  
                  { err_no = EN16; goto error; }
		
		if (
                    
                    sscanf ( in_line,"%*c %ld %ld %ld",
                                      &head, &tail, &cap )
                    != ARC_FIELDS
                   ) 
                    
                    { err_no = EN15; goto error; }

		if ( tail < 0  ||  tail > n  ||
                     head < 0  ||  head > n  
		   )
                    
		    { err_no = EN17; goto error; }

               
	
  
  no_alines ++;
  nodeP* np = *nodes_ad+head;
  adjustEdgesIfNecessary(np);

  edgeP* newedge = np->edges+np->nIdx;
  newedge->endNode = tail;
  newedge->w = 1;
  newedge->cap = cap;
  np->nIdx ++;

  edgeP* newedge_prev = newedge;
  
  np = *nodes_ad + tail;
  adjustEdgesIfNecessary(np);

  newedge = np->edges + np->nIdx;
  newedge->endNode = head;
  newedge->w = 1;
  newedge->cap = cap;
  np->nIdx ++;
  
  newedge_prev->rev = newedge;
  newedge->rev = newedge_prev;

		break;

	default:
		
		err_no = EN18; goto error;
		break;

    } 
}     

 

if ( feof (stdin) == 0 ) 
  { err_no=EN21; goto error; } 

if ( no_lines == 0 ) 
  { err_no = EN22; goto error; } 

if ( no_alines < m ) 
  { err_no = EN19; goto error; } 











 




























	    








































































  









return (0);


 error:  

printf ( "\nline %ld of input - %s\n", 
         no_lines, err_message[err_no] );

exit (1);

}

