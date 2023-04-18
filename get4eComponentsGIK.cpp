//Runs with ./get4eComponentsGIK <input_graph> <output_4components>

#include <stdio.h>
#include <stdlib.h>
#include <chrono>

void get_adj(long,long,long*,long**,long**);
void read_graph(char*,long*,long**,long**);

void DFS(long,long*,long*,long**,long**,long**);
void get_low(long,long*,long*,long*,long*,long*,long**);
void get_l_and_bcount(long,long*,long*,long*,long*,long*,long**,long**);
void get_lowChildren(long,long*,long*,long*,long*,long**,long**);
void get_M(long,long*,long*,long*,long*,long*,long*,long**,long**);

long get_4e_components_high(long,long*,long*,long**);
long get_4e_components_connected_high(long,long*,long*,long**);
long get_4e_components_2econnected_high(long,long*,long*,long**);
long get_4e_components_3econnected_high(long,long*,long*,long**,long*,long**);
void get_high(long,long*,long*,long*,long*,long*,long**,long**);
void get_allM(long,long*,long*,long*,long*,long*,long*,long*,long**,long**,long**);
void get_2low(long,long*,long*,long*,long*,long*,long**,long**,long**,long**);
long find(long*,long);
void unite(long*,long*,long,long);

using namespace std::chrono;

long numberOfComponents=0;
long numberOf2eComponents=0;
long numberOf3eComponents=0;
long numberOf4eComponents=0;
long numberOf1Cuts=0;
long numberOf2Cuts=0;
long numberOf3Cuts=0;
double timeForComponents=0;
double timeFor2eComponents=0;
double timeFor3eComponents=0;
double timeFor4eComponents=0;

int main(int n_args, char** args)
{
   long n; long* adj; long* firstOut;
   read_graph(args[1],&n,&adj,&firstOut);
   long* C;
high_resolution_clock::time_point t1 = high_resolution_clock::now();
   long k=get_4e_components_high(n,adj,firstOut,&C);
high_resolution_clock::time_point t2 = high_resolution_clock::now();
   FILE* fp = fopen(args[2],"w");
   fprintf(fp,"%ld\n",n);
   for(long i=0;i<n;i++){fprintf(fp,"%ld\n",C[i]);}
duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
   fprintf(fp,"%f\n",time_span.count());
   fprintf(fp,"%f\n",timeFor4eComponents);
   fclose(fp);
   free(adj); free(firstOut); free(C);

printf("number of components: %ld\n",numberOfComponents);
printf("time to compute components: %f\n\n",timeForComponents);

printf("number of 1-cuts: %ld\n",numberOf1Cuts);
printf("number of 2e-components: %ld\n",numberOf2eComponents);
printf("time to compute 2e-components: %f\n\n",timeFor2eComponents);

printf("number of 2-cuts: %ld\n",numberOf2Cuts);
printf("number of 3e-components: %ld\n",numberOf3eComponents);
printf("time to compute 3e-components: %f\n\n",timeFor3eComponents);

printf("number of 3-cuts: %ld\n",numberOf3Cuts);
printf("number of 4e-components: %ld\n",numberOf4eComponents);
printf("time to compute 4e-components: %f\n\n",timeFor4eComponents);

printf("total time: %f\n",time_span.count());

   return 0;
}

long get_4e_components_high(long n, long* adj, long* firstOut, long** C)
{
high_resolution_clock::time_point t1;
high_resolution_clock::time_point t2;
duration<double> time_span;
double temp_time=0;

t1 = high_resolution_clock::now();
   (*C) = (long*)malloc(sizeof(long)*n);
   for(long i=0;i<n;i++){(*C)[i]=-1;}
   long* Q = (long*)malloc(sizeof(long)*n);
   long* edges = (long*)malloc(sizeof(long)*firstOut[n]);
   char* found = (char*)malloc(sizeof(char)*n);
   long* map = (long*)malloc(sizeof(long)*n);
   long* imap = (long*)malloc(sizeof(long)*n);
   for(long i=0;i<n;i++){found[i]=0;}
   long k=0;
   for(long r=0;r<n;r++)
   {
      if(found[r]){continue;}
      long Nr=0;
      long first=0; long last=0;
      Q[last++]=r; found[r]=1; map[r]=Nr; imap[Nr++]=r;
      while(first!=last)
      {
         long x=Q[first++];
         for(long i=firstOut[x];i<firstOut[x+1];i++)
         {
            long y=adj[i];
            if(!found[y])
            {
               Q[last++]=y; found[y]=1; map[y]=Nr; imap[Nr++]=y;
            }
         } 
      }
      long nC = last;
      if(nC>1)
      {
high_resolution_clock::time_point temp_t1 = high_resolution_clock::now();
         long edgeIndx=0;
         for(long i=0;i<nC;i++)
         {
            long x=Q[i];
            for(long j=firstOut[x];j<firstOut[x+1];j++)
            {
               long y=adj[j];
               if(x<y){edges[2*edgeIndx]=map[x];edges[2*edgeIndx+1]=map[y];edgeIndx++;}
            }
         }
         long* adjC; long* firstOutC; 
         get_adj(nC,edgeIndx,edges,&adjC,&firstOutC);
         long* C2;
         long k2 = get_4e_components_connected_high(nC,adjC,firstOutC,&C2);
         for(long i=0;i<nC;i++){(*C)[imap[i]]=C2[i]+k;}
         k+=k2;
         free(adjC); free(firstOutC); free(C2);
high_resolution_clock::time_point temp_t2 = high_resolution_clock::now();
duration<double> temp_time_span =  duration_cast<duration<double>>(temp_t2 - temp_t1);
temp_time += temp_time_span.count();
      }
      else
      {
         (*C)[r]=k; 
         k++;
numberOf2eComponents++;
numberOf3eComponents++;
numberOf4eComponents++;
      }    
numberOfComponents++; 
   }
   free(Q); free(edges); free(found);  
   free(map); free(imap);
t2 = high_resolution_clock::now();

time_span = duration_cast<duration<double>>(t2 - t1);
timeForComponents = time_span.count()-temp_time;
   return k;
}

long get_4e_components_connected_high(long n, long* adj, long* firstOut, long** C)
{
high_resolution_clock::time_point t1;
high_resolution_clock::time_point t2;
duration<double> time_span;
double temp_time=0;

t1 = high_resolution_clock::now();
   (*C) = (long*)malloc(sizeof(long)*n);
   long* dfs; long* idfs; long* p;
   DFS(n,adj,firstOut,&dfs,&idfs,&p);
   long* low;
   get_low(n,adj,firstOut,dfs,idfs,p,&low);
   long* Q = (long*)malloc(sizeof(long)*n);
   long* map = (long*)malloc(sizeof(long)*n);
   long* imap = (long*)malloc(sizeof(long)*n);   
   char* found = (char*)malloc(sizeof(char)*n);
   long* edges = (long*)malloc(sizeof(long)*firstOut[n]);
   for(long i=0;i<n;i++){found[i]=0;}
   long k=0;
for(long v=1;v<n;v++){numberOf1Cuts+=low[v]==v;}
   for(long r=0;r<n;r++)
   {
      if(found[r]){continue;}
      long Nr=0;
      long first=0; long last=0;
      Q[last++]=r; found[r]=1; map[r]=Nr; imap[Nr++]=r;
      while(first!=last)
      {
         long x=Q[first++];
         for(long i=firstOut[x];i<firstOut[x+1];i++)
         {
            long y=adj[i];
            if(found[y]){continue;}
            if((x==p[y]&&low[y]==y)||(y==p[x]&&low[x]==x)){continue;}
            Q[last++]=y; found[y]=1; map[y]=Nr; imap[Nr++]=y;
         }
      }
      long nC = last;
      if(nC>1)
      {
high_resolution_clock::time_point temp_t1 = high_resolution_clock::now();
         long edgeIndx=0;
         for(long i=0;i<last;i++)
         {
            long x=Q[i];
            for(long j=firstOut[x];j<firstOut[x+1];j++)
            {
               long y=adj[j];
               if(y<x){continue;}
               if((x==p[y]&&low[y]==y)||(y==p[x]&&low[x]==x)){continue;}
               edges[2*edgeIndx]=map[x]; edges[2*edgeIndx+1]=map[y]; edgeIndx++;
            }
         }
         long* adjC; long* firstOutC;
         get_adj(nC,edgeIndx,edges,&adjC,&firstOutC);
         long* C3;
         long k3=get_4e_components_2econnected_high(nC,adjC,firstOutC,&C3);
         for(long i=0;i<nC;i++){(*C)[imap[i]]=C3[i]+k;}
         k+=k3;
         free(adjC); free(firstOutC); free(C3);
high_resolution_clock::time_point temp_t2 = high_resolution_clock::now();
duration<double> temp_time_span =  duration_cast<duration<double>>(temp_t2 - temp_t1);
temp_time += temp_time_span.count();
      }
      else
      { 
         (*C)[r]=k;
         k++;
numberOf3eComponents++;
numberOf4eComponents++;
      }
numberOf2eComponents++;
   }
   free(dfs); free(idfs); free(p); free(low);
   free(Q); free(found); free(map); free(imap); free(edges);
t2 = high_resolution_clock::now();

time_span = duration_cast<duration<double>>(t2 - t1);
timeFor2eComponents += time_span.count()-temp_time;
   return k;
}

long get_4e_components_2econnected_high(long n, long* adj, long* firstOut, long** C3)
{
high_resolution_clock::time_point t1;
high_resolution_clock::time_point t2;
duration<double> time_span;
double temp_time=0;

t1 = high_resolution_clock::now();
   long* dfs; long* idfs; long* p;
   DFS(n,adj,firstOut,&dfs,&idfs,&p);
   long* low;
   get_low(n,adj,firstOut,dfs,idfs,p,&low);
   long* l; long* bcount;
   get_l_and_bcount(n,adj,firstOut,dfs,idfs,p,&l,&bcount);
   long* low1C; long* low2C;
   get_lowChildren(n,dfs,idfs,p,low,&low1C,&low2C);
   long* M; long* nextM;
   get_M(n,dfs,idfs,l,low,low1C,low2C,&M,&nextM);
   long* prevM = (long*)malloc(sizeof(long)*n);
   for(long i=0;i<n;i++){prevM[i]=-1;}
   for(long i=1;i<n;i++){if(nextM[i]!=-1){prevM[nextM[i]]=i;}}

   long* vEdgeStack = (long*)malloc(sizeof(long)*4*n);
   long* vEdgeFirst = (long*)malloc(sizeof(long)*n);
   long* vEdgeNext = (long*)malloc(sizeof(long)*4*n);
   long* firstVertex = (long*)malloc(sizeof(long)*4*n);
   for(long i=0;i<n;i++){vEdgeFirst[i]=-1;}
   long SP=0;

   char* isCutEdge = (char*)malloc(sizeof(char)*n);
   char* isCutEdgeLow = (char*)malloc(sizeof(char)*n);
   for(long i=0;i<n;i++){isCutEdge[i]=0;isCutEdgeLow[i]=0;}

   long* cutEdgeStack = (long*)malloc(sizeof(long)*4*n);
   long* cutEdgeFirst = (long*)malloc(sizeof(long)*n);
   long* cutEdgeNext = (long*)malloc(sizeof(long)*2*n);
   for(long i=0;i<n;i++){cutEdgeFirst[i]=-1;}
   long cutEdgeSP=0;
long* cutEdgeCount = (long*)malloc(sizeof(long)*n);
for(long i=0;i<n;i++){cutEdgeCount[i]=0;}
   for(long m=1;m<n;m++)
   {
      if(M[m]!=m){continue;}
      long u=m;
      while(u!=-1)
      {
         long z=nextM[u];
         if(z!=-1 && bcount[z]==bcount[u])
         {
            long last;
            while(z!=-1 && bcount[z]==bcount[u]){last=z; z=nextM[z];}
            if(bcount[u]==1)
            {
               isCutEdgeLow[m]=1;
               if(u!=m)
               {
                  vEdgeNext[SP]=vEdgeFirst[u]; vEdgeFirst[u]=SP; vEdgeStack[SP]=m; firstVertex[SP++]=u;
                  vEdgeNext[SP]=vEdgeFirst[m]; vEdgeFirst[m]=SP; vEdgeStack[SP]=u; firstVertex[SP++]=u;
               }
               if(p[last]!=low[m])
               {
                  vEdgeNext[SP]=vEdgeFirst[p[last]]; vEdgeFirst[p[last]]=SP; vEdgeStack[SP]=low[m]; firstVertex[SP++]=u;
                  vEdgeNext[SP]=vEdgeFirst[low[m]]; vEdgeFirst[low[m]]=SP; vEdgeStack[SP]=p[last]; firstVertex[SP++]=u;
               }
               cutEdgeNext[cutEdgeSP]=cutEdgeFirst[u]; cutEdgeFirst[u]=cutEdgeSP;
               cutEdgeStack[2*cutEdgeSP]=m; cutEdgeStack[2*cutEdgeSP+1]=low[m]; cutEdgeSP++;
               cutEdgeCount[u]++;
            }
            else
            {
               vEdgeNext[SP]=vEdgeFirst[u]; vEdgeFirst[u]=SP; vEdgeStack[SP]=p[last]; firstVertex[SP++]=u;
               vEdgeNext[SP]=vEdgeFirst[p[last]]; vEdgeFirst[p[last]]=SP; vEdgeStack[SP]=u; firstVertex[SP++]=u;
            }
            long v=u;
            while(v!=nextM[last])
            {
               isCutEdge[v]=1;
               if(v!=last&&p[v]!=nextM[v])
               {
                  vEdgeNext[SP]=vEdgeFirst[p[v]]; vEdgeFirst[p[v]]=SP; vEdgeStack[SP]=nextM[v]; firstVertex[SP++]=u;
                  vEdgeNext[SP]=vEdgeFirst[nextM[v]]; vEdgeFirst[nextM[v]]=SP; vEdgeStack[SP]=p[v]; firstVertex[SP++]=u;
               }
               cutEdgeNext[cutEdgeSP]=cutEdgeFirst[u]; cutEdgeFirst[u]=cutEdgeSP;
               cutEdgeStack[2*cutEdgeSP]=v; cutEdgeStack[2*cutEdgeSP+1]=p[v]; cutEdgeSP++;
               cutEdgeCount[u]++;
               v=nextM[v];
            }
         }
         else if(bcount[u]==1)
         {
            isCutEdge[u]=1; isCutEdgeLow[m]=1;
            if(u!=m)
            {
               vEdgeNext[SP]=vEdgeFirst[u]; vEdgeFirst[u]=SP; vEdgeStack[SP]=m; firstVertex[SP++]=u;
               vEdgeNext[SP]=vEdgeFirst[m]; vEdgeFirst[m]=SP; vEdgeStack[SP]=u; firstVertex[SP++]=u;
            }
            if(p[u]!=low[m])
            {
               vEdgeNext[SP]=vEdgeFirst[p[u]]; vEdgeFirst[p[u]]=SP; vEdgeStack[SP]=low[m]; firstVertex[SP++]=u;
               vEdgeNext[SP]=vEdgeFirst[low[m]]; vEdgeFirst[low[m]]=SP; vEdgeStack[SP]=p[u]; firstVertex[SP++]=u;
            }
            cutEdgeNext[cutEdgeSP]=cutEdgeFirst[u]; cutEdgeFirst[u]=cutEdgeSP;
            cutEdgeStack[2*cutEdgeSP]=u; cutEdgeStack[2*cutEdgeSP+1]=p[u]; cutEdgeSP++;
            cutEdgeNext[cutEdgeSP]=cutEdgeFirst[u]; cutEdgeFirst[u]=cutEdgeSP;
            cutEdgeStack[2*cutEdgeSP]=m; cutEdgeStack[2*cutEdgeSP+1]=low[m]; cutEdgeSP++;
            cutEdgeCount[u]=2;
         }
         u=z;
      } 
   } 

   long* C = (long*)malloc(sizeof(long)*n);
   long* Q = (long*)malloc(sizeof(long)*n);
   long* map = (long*)malloc(sizeof(long)*n);
   long* imapStack = (long*)malloc(sizeof(long)*n);
   long* imapFirst = (long*)malloc(sizeof(long)*n);
   long* imapNext = (long*)malloc(sizeof(long)*n);
   for(long i=0;i<n;i++){C[i]=-1;imapFirst[i]=-1;}
   SP=0;
   long k=0;
   for(long r=0;r<n;r++)
   {
      if(C[r]!=-1){continue;}
      long first=0; long last=0;
      C[r]=k; Q[last++]=r;
      long Nr=0;
      map[r]=Nr++;
      imapNext[SP]=imapFirst[k]; imapFirst[k]=SP; imapStack[SP++]=r;
      while(first!=last)
      {
         long v=Q[first++];
         for(long i=firstOut[v];i<firstOut[v+1];i++)
         {
            long u=adj[i];
            if((u==p[v]&&isCutEdge[v])||(v==p[u]&&isCutEdge[u])){continue;}
            if((M[u]==u&&v==low[u]&&isCutEdgeLow[u])||(M[v]==v&&u==low[v]&&isCutEdgeLow[v])){continue;}
            if(C[u]==-1)
            {
               C[u]=k; Q[last++]=u; map[u]=Nr++;
               imapNext[SP]=imapFirst[k]; imapFirst[k]=SP; imapStack[SP++]=u;
            }
         }
         for(long i=vEdgeFirst[v];i!=-1;i=vEdgeNext[i])
         {
            long u=vEdgeStack[i];
            if(C[u]==-1)
            {
               C[u]=k; Q[last++]=u; map[u]=Nr++;
               imapNext[SP]=imapFirst[k]; imapFirst[k]=SP; imapStack[SP++]=u;
            }
         }
      }
      k++;
   }

numberOf3eComponents+=k;
for(long v=1;v<n;v++){numberOf2Cuts+=bcount[v]==1;}
for(long x=1;x<n;x++)
{
   if(M[x]!=x){continue;}
   long u=x;
   while(u!=-1)
   {
      long v=u;
      long temp_num=1;
      while(nextM[v]!=-1 && bcount[nextM[v]]==bcount[v])
      {
         temp_num++;
         v=nextM[v];
      }
      numberOf2Cuts+=(temp_num*(temp_num-1))/2;
      u=nextM[v];
   }
}
t2 = high_resolution_clock::now();

time_span = duration_cast<duration<double>>(t2 - t1);
timeFor3eComponents += time_span.count();

t1 = high_resolution_clock::now();
   long* nC = (long*)malloc(sizeof(long)*k);
   long** adjC = (long**)malloc(sizeof(long*)*k);
   long** firstOutC = (long**)malloc(sizeof(long*)*k);
   long* edges = (long*)malloc(sizeof(long)*(firstOut[n]+n));
   for(long c=0;c<k;c++)
   {
      nC[c]=0;
      long eIndx=0;
      for(long t=imapFirst[c];t!=-1;t=imapNext[t])
      {
         nC[c]++;
         long x=imapStack[t];
         for(long i=firstOut[x];i<firstOut[x+1];i++)
         {
            long y=adj[i];
            if(dfs[x]<dfs[y]){continue;}
            if((y==p[x]&&isCutEdge[x])||(x==M[x]&&y==low[x]&&isCutEdgeLow[x])){continue;}
            edges[2*eIndx]=map[x]; edges[2*eIndx+1]=map[y]; eIndx++;
         }
         for(long i=vEdgeFirst[x];i!=-1;i=vEdgeNext[i])
         {
            long y=vEdgeStack[i];
            if(y>x){continue;}
            edges[2*eIndx]=map[x]; edges[2*eIndx+1]=map[y]; eIndx++; 
         }
      }
      if(nC[c]==1){continue;}
      get_adj(nC[c],eIndx,edges,adjC+c,firstOutC+c);
   }

   long* nCcuts = (long*)malloc(sizeof(long)*k);
   long** Ccuts = (long**)malloc(sizeof(long*)*k);

   (*C3) = (long*)malloc(sizeof(long)*n);
   long c_num=0;
   for(long c=0;c<k;c++)
   {
      if(nC[c]>1)
      {
         long* C4;
         long k4=get_4e_components_3econnected_high(nC[c],adjC[c],firstOutC[c],&C4,nCcuts+c,Ccuts+c);
         for(long t=imapFirst[c];t!=-1;t=imapNext[t])
         {
            long x=imapStack[t];
            (*C3)[x]=C4[map[x]]+c_num;
         }
         c_num+=k4;
         free(C4);
numberOf4eComponents+=k4;
      }
      else
      {
         long x=imapStack[imapFirst[c]];
         (*C3)[x]=c_num;
         c_num++;    
numberOf4eComponents++;
      }
   }
t2 = high_resolution_clock::now();

time_span = duration_cast<duration<double>>(t2 - t1);
timeFor4eComponents += time_span.count();

   /*get the number of 3-edge cuts*/
   long* imap = (long*)malloc(sizeof(long)*n);
   for(long c=0;c<k;c++)
   {
      if(nC[c]==1){continue;}
      SP=1; long temp;
      for(long i=imapFirst[c];i!=-1;i=imapNext[i])
      {
         imap[nC[c]-SP]=imapStack[i]; SP++;
      }
      for(long i=0;i<nCcuts[c];i++)
      {
         for(long t=0;t<6;t++){Ccuts[c][6*i+t]=imap[Ccuts[c][6*i+t]];}
      }
   }

   long* cut3Stack_ = (long*)malloc(sizeof(long)*12*n);
   long* cut3First_ = (long*)malloc(sizeof(long)*n);
   long* cut3Next_ = (long*)malloc(sizeof(long)*12*n);
   long* CIndx_ = (long*)malloc(sizeof(long)*12*n);
   long* cutIndx_ = (long*)malloc(sizeof(long)*12*n);
   long* edgeIndx_ = (long*)malloc(sizeof(long)*12*n);
   for(long i=0;i<n;i++){cut3First_[i]=-1;}
   SP=0;
   for(long c=0;c<k;c++)
   {
      if(nC[c]==1){continue;}
      for(long i=0;i<nCcuts[c];i++)
      {
         for(long t=0;t<3;t++)
         {
            long x=Ccuts[c][6*i+2*t]; long y=Ccuts[c][6*i+2*t+1];
            cut3Next_[SP]=cut3First_[x]; cut3First_[x]=SP; cut3Stack_[SP]=y; CIndx_[SP]=c; cutIndx_[SP]=i; edgeIndx_[SP++]=t;
            cut3Next_[SP]=cut3First_[y]; cut3First_[y]=SP; cut3Stack_[SP]=x; CIndx_[SP]=c; cutIndx_[SP]=i; edgeIndx_[SP++]=t;
         }
      }
   }
   long* cut3Stack = (long*)malloc(sizeof(long)*12*n);
   long* cut3First = (long*)malloc(sizeof(long)*n);
   long* cut3Next = (long*)malloc(sizeof(long)*12*n);
   long* CIndx = (long*)malloc(sizeof(long)*12*n);
   long* cutIndx = (long*)malloc(sizeof(long)*12*n);
   long* edgeIndx = (long*)malloc(sizeof(long)*12*n);
   for(long i=0;i<n;i++){cut3First[i]=-1;}
   SP=0;
   for(long x=0;x<n;x++)
   {
      for(long i=cut3First_[x];i!=-1;i=cut3Next_[i])
      {
         long y=cut3Stack_[i];
         cut3Next[SP]=cut3First[y]; cut3First[y]=SP; cut3Stack[SP]=x; 
         CIndx[SP]=CIndx_[i]; cutIndx[SP]=cutIndx_[i]; edgeIndx[SP++]=edgeIndx_[i];
      }
   } 
   long* vEdgesStackS = (long*)malloc(sizeof(long)*4*n);
   long* vEdgesFirstS = (long*)malloc(sizeof(long)*n);
   long* vEdgesNextS = (long*)malloc(sizeof(long)*4*n);
   long* firstVertexS = (long*)malloc(sizeof(long)*4*n);
   for(long i=0;i<n;i++){vEdgesFirstS[i]=-1;}
   SP=0;
   for(long x=0;x<n;x++)
   {
      for(long i=vEdgeFirst[x];i!=-1;i=vEdgeNext[i])
      {
         long y=vEdgeStack[i];
         vEdgesNextS[SP]=vEdgesFirstS[y]; vEdgesFirstS[y]=SP; vEdgesStackS[SP]=x; firstVertexS[SP++]=firstVertex[i];
      }   
   }

   long** corVertex = (long**)malloc(sizeof(long*)*k);
   for(long c=0;c<k;c++){if(nC[c]!=1){corVertex[c]=(long*)malloc(sizeof(long)*nCcuts[c]*3);}}
   for(long c=0;c<k;c++){if(nC[c]!=1){for(long i=0;i<3*nCcuts[c];i++){corVertex[c][i]=-1;}}}

   for(long x=0;x<n;x++)
   {
      for(long i=cut3First[x];i!=-1;i=cut3Next[i])
      {
         long y=cut3Stack[i];
         if(y<x){continue;}
         long j=vEdgesFirstS[x];
         while(j!=-1)
         {
            long z=vEdgesStackS[j];
            if(z>y){j=vEdgesNextS[j]; vEdgesFirstS[x]=j; continue;}
            if(z<y){break;}
            long c=CIndx[i]; long t=cutIndx[i]; long e=edgeIndx[i];
            corVertex[c][3*t+e]=firstVertexS[j];
            vEdgesFirstS[x]=j; 

            long i1=cut3Next[i];
            if(i1==-1){break;}
            long y1=cut3Stack[i1];
            if(y1!=y){break;}
            if(cutIndx[i1]!=t){break;}
            i=i1;
            j=vEdgesNextS[j]; vEdgesFirstS[x]=j;
            if(j==-1){break;}
            z=vEdgesStackS[j];
            if(z!=y){break;}
            e=edgeIndx[i];
            corVertex[c][3*t+e]=firstVertexS[j];

            i1=cut3Next[i];
            if(i1==-1){break;}
            y1=cut3Stack[i1];
            if(y1!=y){break;}
            if(cutIndx[i1]!=t){break;}
            i=i1;
            j=vEdgesNextS[j]; vEdgesFirstS[x]=j;
            if(j==-1){break;}
            z=vEdgesStackS[j];
            if(z!=y){break;}
            e=edgeIndx[i];
            corVertex[c][3*t+e]=firstVertexS[j];
            break;
         } 
      }
   }

   long n3Cuts=0;
   for(long c=0;c<k;c++)
   {
      if(nC[c]==1){continue;}
      for(long i=0;i<nCcuts[c];i++)
      {
         long num=1;
         for(long t=0;t<3;t++)
         {
            if(corVertex[c][3*i+t]!=-1){num*=cutEdgeCount[corVertex[c][3*i+t]];}
         }
         n3Cuts+=num;
      } 
   }
numberOf3Cuts+=n3Cuts;
   /*end: get the number of 3-edge cuts*/

   free(vEdgeStack); free(vEdgeFirst); free(vEdgeNext); free(firstVertex);
   free(vEdgesStackS); free(vEdgesFirstS); free(vEdgesNextS); free(firstVertexS);

   free(cutEdgeStack); free(cutEdgeFirst); free(cutEdgeNext); free(cutEdgeCount);

   free(cut3Stack_); free(cut3First_); free(cut3Next_); free(CIndx_); free(cutIndx_); free(edgeIndx_);
   free(cut3Stack); free(cut3First); free(cut3Next); free(CIndx); free(cutIndx); free(edgeIndx);

   free(imap);
   for(long i=0;i<k;i++){if(nC[i]>1){free(Ccuts[i]); free(corVertex[i]);}}
   free(nCcuts); free(Ccuts); free(corVertex);

   for(long i=0;i<k;i++){if(nC[i]>1){free(adjC[i]);free(firstOutC[i]);}}
   free(adjC); free(firstOutC); free(edges); free(nC);

   free(dfs); free(idfs); free(p); free(low);
   free(l); free(bcount); free(low1C); free(low2C);
   free(M); free(nextM); free(prevM);
   free(isCutEdge); free(isCutEdgeLow);
   free(C); free(Q); free(map); free(imapStack); free(imapFirst); free(imapNext);

   return c_num;
}

long get_4e_components_3econnected_high(long n, long* adj, long* firstOut, long** C, long* n3cuts, long** cuts3)
{
   long* cuts = (long*)malloc(sizeof(long)*n*12);
   long* dfs; long* idfs; long* p;
   DFS(n,adj,firstOut,&dfs,&idfs,&p);
   long* low; long* lowD; long* low2; long* low2D;
   get_2low(n,adj,firstOut,dfs,idfs,p,&low,&lowD,&low2,&low2D);
   long* low1C; long* low2C;
   get_lowChildren(n,dfs,idfs,p,low,&low1C,&low2C);
   long* l; long* bcount;
   get_l_and_bcount(n,adj,firstOut,dfs,idfs,p,&l,&bcount);
   long* M; long* nextM;
   get_M(n,dfs,idfs,l,low,low1C,low2C,&M,&nextM);
   long* high; long* highD;
   get_high(n,adj,firstOut,dfs,idfs,p,&high,&highD);
   long* Ml; long* Mlow1; long* Mlow2;
   get_allM(n,dfs,idfs,l,low,M,low1C,low2C,&Ml,&Mlow1,&Mlow2);
   long* currentVertex = (long*)malloc(sizeof(long)*n);

   long* ND = (long*)malloc(sizeof(long)*n);
   for(long i=0;i<n;i++){ND[i]=1;}
   for(long i=n-1;i>0;i--){long v=idfs[i];ND[p[v]]+=ND[v];}

   long* size = (long*)malloc(sizeof(long)*2*n);
   char* type = (char*)malloc(sizeof(char)*2*n);
   long k=0;

   //one tree-edge
   for(long v=1;v<n;v++)
   {
      if(bcount[v]==2)
      {
         cuts[6*k+0]=v; cuts[6*k+1]=p[v];
         cuts[6*k+2]=low[v]; cuts[6*k+3]=lowD[v];
         cuts[6*k+4]=low2[v]; cuts[6*k+5]=low2D[v];
         type[k]=0; size[k]=ND[v];
         k++;      
      }
   }

   //two tree-edges
   //upper case
   for(long i=0;i<n;i++){currentVertex[i]=i;}
   for(long i=n-1;i>0;i--)
   {
      long v=idfs[i];
      long m=Ml[v];
      if(m==-1){continue;} 
      long u=currentVertex[m];
      while(nextM[u]!=-1 && dfs[nextM[u]]>dfs[v]){u=nextM[u];}
      currentVertex[m]=u;
      if(dfs[high[u]]<dfs[v] && bcount[v]==bcount[u]+1)
      {
         cuts[6*k+0]=u; cuts[6*k+1]=p[u];
         cuts[6*k+2]=v; cuts[6*k+3]=p[v];
         cuts[6*k+4]=M[v]; cuts[6*k+5]=l[M[v]];
         type[k]=1; size[k]=ND[v]-ND[u];
         k++;
      }
   }
   for(long i=0;i<n;i++){currentVertex[i]=i;}
   for(long i=n-1;i>0;i--)
   {
      long v=idfs[i];
      long m=Mlow1[v];
      if(m==-1){continue;} 
      long u=currentVertex[m];
      while(nextM[u]!=-1 && dfs[nextM[u]]>dfs[v]){u=nextM[u];}
      currentVertex[m]=u;
      if(dfs[high[u]]<dfs[v] && bcount[v]==bcount[u]+1)
      {
         cuts[6*k+0]=u; cuts[6*k+1]=p[u];
         cuts[6*k+2]=v; cuts[6*k+3]=p[v];
         cuts[6*k+4]=Mlow2[v]; cuts[6*k+5]=l[Mlow2[v]];
         type[k]=1; size[k]=ND[v]-ND[u];
         k++;
      }
   }
   for(long i=0;i<n;i++){currentVertex[i]=i;}
   for(long i=n-1;i>0;i--)
   {
      long v=idfs[i];
      long m=Mlow2[v];
      if(m==-1){continue;} 
      long u=currentVertex[m];
      while(nextM[u]!=-1 && dfs[nextM[u]]>dfs[v]){u=nextM[u];}
      currentVertex[m]=u;
      if(dfs[high[u]]<dfs[v] && bcount[v]==bcount[u]+1)
      {
         cuts[6*k+0]=u; cuts[6*k+1]=p[u];
         cuts[6*k+2]=v; cuts[6*k+3]=p[v];
         cuts[6*k+4]=Mlow1[v]; cuts[6*k+5]=l[Mlow1[v]];
         type[k]=1; size[k]=ND[v]-ND[u];
         k++;
      }
   }

   //lower case
   for(long u=1;u<n;u++)
   {
      if(nextM[u]!=-1 && bcount[nextM[u]]==bcount[u]-1)
      {
         cuts[6*k+0]=u; cuts[6*k+1]=p[u];
         cuts[6*k+2]=nextM[u]; cuts[6*k+3]=p[nextM[u]];
         cuts[6*k+4]=highD[u]; cuts[6*k+5]=high[u];
         type[k]=1; size[k]=ND[nextM[u]]-ND[u];
         k++;
      }
   }
   for(long i=0;i<n;i++){currentVertex[i]=i;}
   for(long i=n-1;i>0;i--)
   {
      long u=idfs[i];
      long m=Ml[u];
      if(m==-1){continue;} 
      long v=currentVertex[m];
      while(v!=-1 && dfs[v]>=dfs[u]){v=nextM[v];}
      currentVertex[m]=v;
      if(v==-1){continue;}
      if(bcount[u]==bcount[v]+1 && M[u]!=M[v])
      {
         cuts[6*k+0]=u; cuts[6*k+1]=p[u];
         cuts[6*k+2]=v; cuts[6*k+3]=p[v];
         cuts[6*k+4]=highD[u]; cuts[6*k+5]=high[u];
         type[k]=1; size[k]=ND[v]-ND[u];
         k++;
      }
   }
   for(long i=0;i<n;i++){currentVertex[i]=i;}
   for(long i=n-1;i>0;i--)
   {
      long u=idfs[i];
      long m=Mlow1[u];
      if(m==-1){continue;} 
      long v=currentVertex[m];
      while(v!=-1 && dfs[v]>=dfs[u]){v=nextM[v];}
      currentVertex[m]=v;
      if(v==-1){continue;}
      if(bcount[u]==bcount[v]+1)
      {
         cuts[6*k+0]=u; cuts[6*k+1]=p[u];
         cuts[6*k+2]=v; cuts[6*k+3]=p[v];
         cuts[6*k+4]=highD[u]; cuts[6*k+5]=high[u];
         type[k]=1; size[k]=ND[v]-ND[u];
         k++;
      }
   }


   //three tree-edges
   //u and v are not related as ancestor and descendant
   for(long i=0;i<n;i++){currentVertex[i]=i;}
   for(long i=n-1;i>0;i--)
   {
      long w=idfs[i];
      long m1=Mlow1[w]; long m2=Mlow2[w];
      if(m1==-1||m2==-1){continue;}
      long u=currentVertex[m1];
      while(nextM[u]!=-1 && dfs[nextM[u]]>dfs[w]){u=nextM[u];}
      currentVertex[m1]=u;
      long v=currentVertex[m2];
      while(nextM[v]!=-1 && dfs[nextM[v]]>dfs[w]){v=nextM[v];}
      currentVertex[m2]=v;
      if(bcount[w]==bcount[u]+bcount[v] && dfs[high[u]]<dfs[w] && dfs[high[v]]<dfs[w])
      {
         cuts[6*k+0]=u; cuts[6*k+1]=p[u];
         cuts[6*k+2]=v; cuts[6*k+3]=p[v];
         cuts[6*k+4]=w; cuts[6*k+5]=p[w];
         type[k]=2; size[k]=ND[w]-ND[u]-ND[v];
         k++;
      }
   }

   long* lastM = (long*)malloc(sizeof(long)*n);
   for(long v=1;v<n;v++){if(nextM[v]==-1){lastM[M[v]]=v;}}   

   //u,v,w M[v]!=M[w]
   for(long i=0;i<n;i++){currentVertex[i]=i;}
   for(long i=n-1;i>0;i--)
   {
      long v=idfs[i];
      long m1=Mlow1[v]; long m2=Mlow2[v];
      if(m1==-1||m2==-1){continue;}
      long w=currentVertex[m1];
      while(w!=-1 && dfs[w]>=dfs[v]){w=nextM[w];}
      currentVertex[m1]=w;
      long u=lastM[m2];
      if(w!=-1 && dfs[high[u]]<dfs[v] && bcount[v]==bcount[u]+bcount[w])
      {
         cuts[6*k+0]=u; cuts[6*k+1]=p[u];
         cuts[6*k+2]=v; cuts[6*k+3]=p[v];
         cuts[6*k+4]=w; cuts[6*k+5]=p[w];
         type[k]=2; size[k]=ND[w]-ND[v]+ND[u];
         k++;
      }
   }

   //w=nextM[v]
   long* A = (long*)malloc(sizeof(long)*firstOut[n]);
   for(long i=0;i<firstOut[n];i++){A[i]=-1;}
   long* Hstack = (long*)malloc(sizeof(long)*n);
   long* firstH = (long*)malloc(sizeof(long)*n);
   long* nextH = (long*)malloc(sizeof(long)*n); 
   for(long i=0;i<n;i++){firstH[i]=-1; nextH[i]=-1;}
   long SP=0;
   for(long i=1;i<n;i++)
   {
      long v=idfs[i];
      long h=high[v];
      nextH[SP]=firstH[h];
      firstH[h]=SP;
      Hstack[SP++]=v;
   }  
   long* stack = (long*)malloc(sizeof(long)*n);
   for(long h=0;h<n;h++)
   {
      long uIndx = firstH[h];
      SP=0;
      while(uIndx!=-1)
      {
         long zIndx = nextH[uIndx];
         if(zIndx==-1){break;}
         long u=Hstack[uIndx]; long z=Hstack[zIndx];
         if(nextM[u]==-1)
         {
            stack[SP++]=u;
            A[bcount[u]]=u; 
         }
         if(!(dfs[z]<=dfs[u] && dfs[u]<dfs[z]+ND[z]))
         {
            while(SP!=0)
            {
               long uTemp = stack[SP-1];
               A[bcount[uTemp]]=-1;
               SP--;
            } 
         }
         if(nextM[z]!=-1)
         {
            long v=z; long w=nextM[v];
            if(A[bcount[v]-bcount[w]]!=-1)
            {
               u=A[bcount[v]-bcount[w]];
               if(dfs[low[u]]>=dfs[w])
               {
                  cuts[6*k+0]=u; cuts[6*k+1]=p[u];
                  cuts[6*k+2]=v; cuts[6*k+3]=p[v];
                  cuts[6*k+4]=w; cuts[6*k+5]=p[w];
                  type[k]=2; size[k]=ND[w]-ND[v]+ND[u];
                  k++;              
               }
            }
         }
         uIndx=zIndx;
      }
      while(SP!=0)
      {
         long uTemp = stack[SP-1];
         A[bcount[uTemp]]=-1; 
         SP--;
      } 
   }

   //w!=nextM[v]
   long* Ustack = (long*)malloc(sizeof(long)*n);
   long* firstU = (long*)malloc(sizeof(long)*n);
   long* nextU = (long*)malloc(sizeof(long)*n);
   for(long i=0;i<n;i++){firstU[i]=-1; nextU[i]=-1;}
   long USP=0;
   for(long h=0;h<n;h++)
   {
      long uIndx=firstH[h];
      SP=0;
      while(uIndx!=-1)
      {
         long zIndx=nextH[uIndx];
         if(zIndx==-1){break;}
         long u=Hstack[uIndx]; long z=Hstack[zIndx];
         if(nextM[u]==-1){stack[SP++]=u;}
         if(!(dfs[z]<=dfs[u] && dfs[u]<dfs[z]+ND[z])){SP=0;}
         if(nextM[z]!=-1)
         {
            long v=z;
            while(SP!=0 && dfs[low[stack[SP-1]]]<dfs[lastM[M[v]]]){SP--;}
            while(SP!=0 && dfs[low[stack[SP-1]]]<dfs[nextM[v]])
            {
               SP--;
               long u=stack[SP];
               nextU[USP]=firstU[v];
               firstU[v]=USP;
               Ustack[USP++]=u;  
            }  
         }
         uIndx=zIndx;
      }
   }

   long* lowestW = (long*)malloc(sizeof(long)*n);
   for(long i=0;i<n;i++){lowestW[i]=nextM[i];}
   for(long i=1;i<n;i++)
   {
      long v=idfs[i];
      long indx=firstU[v];
      while(indx!=-1)
      {
         long u=Ustack[indx];
         long w=lowestW[v];
         while(dfs[w]>dfs[low[u]]){w=lowestW[w];}
         lowestW[v]=w;
         if(bcount[v]==bcount[u]+bcount[w])
         {
            cuts[6*k+0]=u; cuts[6*k+1]=p[u];
            cuts[6*k+2]=v; cuts[6*k+3]=p[v];
            cuts[6*k+4]=w; cuts[6*k+5]=p[w];
            type[k]=2; size[k]=ND[w]-ND[v]+ND[u];
            k++;  
         }
         indx=nextU[indx];
      }
   }

   long* sizeStack = (long*)malloc(sizeof(long)*k);
   long* sizeFirst = (long*)malloc(sizeof(long)*n);
   long* sizeNext = (long*)malloc(sizeof(long)*k);
   for(long i=0;i<n;i++){sizeFirst[i]=-1;}
   SP=0;
   for(long i=0;i<k;i++)
   {
      long s=size[i];
      sizeNext[SP]=sizeFirst[s]; sizeFirst[s]=SP; sizeStack[SP++]=i;
   } 
   long* cutStack = (long*)malloc(sizeof(long)*3*k);
   long* cutFirst = (long*)malloc(sizeof(long)*n);
   long* cutNext = (long*)malloc(sizeof(long)*3*k);
   for(long i=0;i<n;i++){cutFirst[i]=-1;}
   SP=0;
   for(long s=0;s<n;s++)
   {
      for(long t=sizeFirst[s];t!=-1;t=sizeNext[t])
      {
         long i=sizeStack[t];
         long u=-1; long v=-1; long w=-1;
         u=cuts[6*i+0];
         if(type[i]>0){v=cuts[6*i+2];}
         if(type[i]>1){w=cuts[6*i+4];}
         cutNext[SP]=cutFirst[u]; cutFirst[u]=SP; cutStack[SP++]=i;
         if(v!=-1)
         {
            cutNext[SP]=cutFirst[v]; cutFirst[v]=SP; cutStack[SP++]=i;
         }
         if(w!=-1)
         {
            cutNext[SP]=cutFirst[w]; cutFirst[w]=SP; cutStack[SP++]=i;
         }
      }
   }

   long* cactusParent = (long*)malloc(sizeof(long)*3*n);
   long* parentCut = (long*)malloc(sizeof(long)*3*n);
   char* marked = (char*)malloc(sizeof(char)*k);
   for(long i=0;i<k;i++){marked[i]=0;}
   long* phi = (long*)malloc(sizeof(long)*k);
   for(long i=0;i<k;i++){phi[i]=-1;}

   long* Cstack = (long*)malloc(sizeof(long)*3*n);
   long* cactusNode = (long*)malloc(sizeof(long)*n);
   (*C) = (long*)malloc(sizeof(long)*n);
   (*C)[0]=0; cactusNode[0]=0;
   long n_comp=1; long Nr=1;
   for(long i=1;i<n;i++)
   {
      long v=idfs[i];
      long x=cactusNode[p[v]];
      SP=0;
      while(x!=0)
      {
         long c=parentCut[x];
         char found=0;
         if(type[c]==0&&v==cuts[6*c+0])
         {
            marked[c]=1; Cstack[SP++]=c; found=1;
         }
         else if(type[c]==1&&(v==cuts[6*c+0]||v==cuts[6*c+2]))
         {
            marked[c]=1; Cstack[SP++]=c; found=1;
         }
         else if(type[c]==2&&(v==cuts[6*c+0]||v==cuts[6*c+2]||v==cuts[6*c+4]))
         {
            marked[c]=1; Cstack[SP++]=c; found=1;
         }
         if(!found){break;}
         x=cactusParent[x];
      }
      char new_comp=0;
      for(long t=cutFirst[v];t!=-1;t=cutNext[t])
      {
         long c=cutStack[t];
         if(marked[c]){continue;}
         if(phi[c]==-1)
         {
            phi[c]=Nr; parentCut[Nr]=c; cactusParent[Nr++]=x; new_comp=1;
         }
         x=phi[c];
      } 
      cactusNode[v]=x;
      while(SP!=0){marked[Cstack[SP-1]]=0;SP--;}
      if(new_comp){x=n_comp++;}
      (*C)[v]=x;
   }

   free(sizeStack); free(sizeFirst); free(sizeNext);
   free(cutStack); free(cutFirst); free(cutNext); free(Cstack);
   free(cactusParent); free(parentCut); free(marked); free(phi); free(cactusNode);

   free(dfs); free(idfs); free(p);
   free(low); free(lowD); free(low2); free(low2D);
   free(low1C); free(low2C); free(l); free(bcount);
   free(M); free(nextM); free(high); free(highD);
   free(Ml); free(Mlow1); free(Mlow2); free(currentVertex); 
   free(lastM); free(A); free(ND); 
   free(Hstack); free(firstH); free(nextH); free(stack);
   free(Ustack); free(firstU); free(nextU);
   free(lowestW);
   free(size); free(type);
   //free(cuts);
   *n3cuts=k; *cuts3=cuts;
   return n_comp;
}

void get_2low(long n, long* adj, long* firstOut, long* dfs, long* idfs, long* p, long** low1, long** low1D, long** low2, long** low2D)
{
   (*low1) = (long*)malloc(sizeof(long)*n);
   (*low1D) = (long*)malloc(sizeof(long)*n);
   (*low2) = (long*)malloc(sizeof(long)*n);
   (*low2D) = (long*)malloc(sizeof(long)*n);
   for(long i=0;i<n;i++){(*low1)[i]=i; (*low2)[i]=i;}
   char* foundP = (char*)malloc(sizeof(char)*n);
   char* foundC = (char*)malloc(sizeof(char)*n);
   for(long i=0;i<n;i++){foundP[i]=0;foundC[i]=0;}
   for(long i=n-1;i>0;i--)
   {
      long v=idfs[i];
      for(long j=firstOut[v];j<firstOut[v+1];j++)
      {
         long u=adj[j];
         if(u==p[v]&&!foundP[v]){foundP[v]=1;continue;}
         if(dfs[u]<dfs[v])
         {
            if(dfs[u]<=dfs[(*low1)[v]])
            { 
               (*low2)[v]=(*low1)[v]; (*low2D)[v]=(*low1D)[v];
               (*low1)[v]=u; (*low1D)[v]=v;
            }
            else if(dfs[u]<dfs[(*low2)[v]])
            {
               (*low2)[v]=u; (*low2D)[v]=v;
            }  
         }
         else if(v==p[u]&&!foundC[u])
         {
            if(dfs[(*low2)[u]]<=dfs[(*low1)[v]])
            {
               (*low1)[v]=(*low1)[u]; (*low1D)[v]=(*low1D)[u];    
               (*low2)[v]=(*low2)[u]; (*low2D)[v]=(*low2D)[u];    
            }
            else if(dfs[(*low1)[u]]<=dfs[(*low1)[v]])
            {
               if(dfs[(*low1)[v]]<dfs[(*low2)[v]])
               {
                  (*low2)[v]=(*low1)[v]; (*low2D)[v]=(*low1D)[v];
               }
               else if(dfs[(*low2)[u]]<dfs[(*low2)[v]])
               {
                  (*low2)[v]=(*low2)[u]; (*low2D)[u]=(*low2D)[u];
               }
               (*low1)[v]=(*low1)[u]; (*low1D)[v]=(*low1D)[u];
            }
            else if(dfs[(*low1)[u]]<dfs[(*low2)[v]])
            {
               (*low2)[v]=(*low1)[u]; (*low2D)[v]=(*low1D)[u];
            }
            foundC[u]=1;
         }
      }
   }
   free(foundP); free(foundC);
}

void get_high(long n, long* adj, long* firstOut, long* dfs, long* idfs, long* p, long** high, long** highD)
{
   (*high) = (long*)malloc(sizeof(long)*n);
   (*highD) = (long*)malloc(sizeof(long)*n);
   long* P = (long*)malloc(sizeof(long)*n);
   long* size = (long*)malloc(sizeof(long)*n);
   long* repr = (long*)malloc(sizeof(long)*n);
   char* found_p = (char*)malloc(sizeof(char)*n);
   for(long i=0;i<n;i++){(*high)[i]=-1; P[i]=i; size[i]=1; repr[i]=i; found_p[i]=0;}
   for(long i=n-1;i>=0;i--)
   {
      long y=idfs[i];
      for(long j=firstOut[y];j<firstOut[y+1];j++)
      {
         long x=adj[j];
         if(p[x]==y&&!found_p[x]){found_p[x]=1;continue;}
         if(dfs[x]<dfs[y]){continue;}
         long u=repr[find(P,x)];
         while(u!=y)
         {
            (*high)[u]=y; (*highD)[u]=x;
            long next=repr[find(P,p[u])];
            unite(P,size,u,p[u]);
            repr[find(P,u)]=y;
            u=next;
         } 
      }
   }
   free(P); free(size); free(repr); free(found_p);
}

long find(long* p, long x)
{
   long r=x;
   while(p[r]!=r){r=p[r];}
   while(x!=r){long next=p[x];p[x]=r;x=next;}
   return r;
}
void unite(long* p, long* size, long x, long y)
{
   long r1=find(p,x);
   long r2=find(p,y);
   if(size[r1]<size[r2]){p[r1]=r2;size[r2]+=size[r1];}
   else{p[r2]=r1;size[r1]+=size[r2];}
}

void get_allM(long n, long* dfs, long* idfs, long* l, long* low, long* M, long* low1C, long* low2C, long** Ml, long** Mlow1, long** Mlow2)
{
   (*Ml) = (long*)malloc(sizeof(long)*n);
   (*Mlow1) = (long*)malloc(sizeof(long)*n);
   (*Mlow2) = (long*)malloc(sizeof(long)*n);
   for(long i=0;i<n;i++)
   {
      (*Ml)[i]=-1; (*Mlow1)[i]=-1; (*Mlow2)[i]=-1;
   }

   long* currentM = (long*)malloc(sizeof(long)*n);
   for(long i=0;i<n;i++){currentM[i]=i;}
   for(long i=n-1;i>0;i--)
   {
      long v=idfs[i];
      long m=M[v];
      if(dfs[l[m]]>=dfs[v]){continue;}
      if(low1C[m]==-1 || dfs[low[low1C[m]]]>=dfs[v]){continue;}
      if(low2C[m]!=-1 && dfs[low[low2C[m]]]<dfs[v])
      {
         (*Ml)[v]=m; 
         continue;
      }
      long tempM=low1C[m];
      m=currentM[low1C[m]];
      while(1)
      {
         if(dfs[l[m]]<dfs[v]){break;}
         if(low2C[m]!=-1 && dfs[low[low2C[m]]]<dfs[v]){break;}
         m=currentM[low1C[m]];
      }
      (*Ml)[v]=m; 
      currentM[tempM]=m;
   }

   for(long i=0;i<n;i++){currentM[i]=i;}
   for(long i=n-1;i>0;i--)
   {
      long v=idfs[i];
      long m=M[v];
      if(dfs[l[m]]<dfs[v]){continue;}
      if(low1C[m]==-1 || dfs[low[low1C[m]]]>=dfs[v]){continue;}
      long tempM=low1C[m];
      m=currentM[low1C[m]];
      while(1)
      {
         if(dfs[l[m]]<dfs[v]){break;}
         if(low2C[m]!=-1 && dfs[low[low2C[m]]]<dfs[v]){break;}
         m=currentM[low1C[m]];
      }
      (*Mlow1)[v]=m; 
      currentM[tempM]=m;
   }

   for(long i=0;i<n;i++){currentM[i]=i;}
   for(long i=n-1;i>0;i--)
   {
      long v=idfs[i];
      long m=M[v];
      if(dfs[l[m]]<dfs[v]){continue;}
      if(low2C[m]==-1 || dfs[low[low2C[m]]]>=dfs[v]){continue;}
      long tempM=low2C[m];
      m=currentM[low2C[m]];
      while(1)
      {
         if(dfs[l[m]]<dfs[v]){break;}
         if(low2C[m]!=-1 && dfs[low[low2C[m]]]<dfs[v]){break;}
         m=currentM[low1C[m]];
      }
      (*Mlow2)[v]=m; 
      currentM[tempM]=m;
   }

   free(currentM);
}

void get_M(long n, long* dfs, long* idfs, long* l, long* low, long* low1C, long* low2C, long** M, long** nextM)
{
   (*M) = (long*)malloc(sizeof(long)*n);
   (*nextM) = (long*)malloc(sizeof(long)*n);
   for(long i=0;i<n;i++){(*M)[i]=-1; (*nextM)[i]=-1;}
   for(long i=n-1;i>0;i--)
   {
      long v=idfs[i];
      long c=v; long m=v;
      while(1)
      {
         if(dfs[l[m]]<i){(*M)[v]=m;break;}
         if(low2C[m]!=-1 && dfs[low[low2C[m]]]<i){(*M)[v]=m;break;}
         c=low1C[m]; m=(*M)[c]; 
      }
      if(c!=v)
      {
         (*nextM)[c]=v;
      }
   } 
}

void get_lowChildren(long n, long* dfs, long* idfs, long* p, long* low, long** low1C, long** low2C)
{
   (*low1C) = (long*)malloc(sizeof(long)*n);
   (*low2C) = (long*)malloc(sizeof(long)*n);
   for(long i=0;i<n;i++){(*low1C)[i]=-1; (*low2C)[i]=-1;}
   for(long i=1;i<n;i++)
   {
      long x=idfs[i];
      long y=p[x];
      if((*low1C)[y]==-1){(*low1C)[y]=x;}
      else if(dfs[low[x]]<dfs[low[(*low1C)[y]]]){(*low1C)[y]=x;}
   }
   for(long i=1;i<n;i++)
   {
      long x=idfs[i];
      long y=p[x];
      if(x!=(*low1C)[y])
      {
        if((*low2C)[y]==-1){(*low2C)[y]=x;}
        else if(dfs[low[x]]<dfs[low[(*low2C)[y]]]){(*low2C)[y]=x;}
      }
   }
}

void get_l_and_bcount(long n, long* adj, long* firstOut, long* dfs, long* idfs, long* p, long** l, long** bcount)
{
   (*l) = (long*)malloc(sizeof(long)*n);
   (*bcount) = (long*)malloc(sizeof(long)*n);
   for(long i=0;i<n;i++){(*l)[i]=i; (*bcount)[i]=0;}
   char* found_p = (char*)malloc(sizeof(char)*n);
   for(long i=0;i<n;i++){found_p[i]=0;}
   char* found_c = (char*)malloc(sizeof(char)*n);
   for(long i=0;i<n;i++){found_c[i]=0;}
   for(long i=n-1;i>0;i--)
   {
      long x=idfs[i];
      for(long t=firstOut[x];t<firstOut[x+1];t++)
      {
         long y=adj[t];
         if(dfs[y]<dfs[x])
         {
            if(y==p[x]&&!found_p[x]) 
            {
               found_p[x]=1;
            }  
            else
            {
               if(dfs[y]<dfs[(*l)[x]]){(*l)[x]=y;}
               (*bcount)[x]++; (*bcount)[y]--;
            }
         }
         else if(x==p[y]&&!found_c[y])
         {
            (*bcount)[x]+=(*bcount)[y];
            found_c[y]=1;     
         }
      }
   }
   free(found_p); free(found_c);
}

void get_low(long n, long* adj, long* firstOut, long* dfs, long* idfs, long* p, long** low)
{
   *low = (long*)malloc(sizeof(long)*n);
   for(long i=0;i<n;i++){(*low)[i]=i;}
   char* found_p = (char*)malloc(sizeof(char)*n);
   for(long i=0;i<n;i++){found_p[i]=0;}
   for(long i=n-1;i>0;i--)
   {
      long v=idfs[i];
      for(long t=firstOut[v];t<firstOut[v+1];t++)
      {
         long u=adj[t];
         if(dfs[u]<dfs[v])
         {
            if(u==p[v] && !found_p[v])
            {
               found_p[v]=1;
            } 
            else
            {
               if(dfs[u]<dfs[(*low)[v]]){(*low)[v]=u;}
            }
         } 
         else if(v==p[u])
         {
            if(dfs[(*low)[u]]<dfs[(*low)[v]]){(*low)[v]=(*low)[u];}
         }
      }
   }
   free(found_p);
}

void DFS(long n, long* adj, long* firstOut, long** dfs, long** idfs, long** p)
{
   *dfs = (long*)malloc(sizeof(long)*n);
   *idfs = (long*)malloc(sizeof(long)*n);
   *p = (long*)malloc(sizeof(long)*n);
   for(long i=0;i<n;i++){(*dfs)[i]=-1; (*p)[i]=-1;}
   long* temp_vertex = (long*)malloc(sizeof(long)*n);
   long* temp_out = (long*)malloc(sizeof(long)*n);
   
   long nr=0;
   (*dfs)[0]=nr; (*idfs)[nr++]=0;
   temp_vertex[0]=0; temp_out[0]=firstOut[0];
   long SP=0;
   while(SP!=-1)
   {
      long v=temp_vertex[SP];
      char descend=0;
      for(long i=temp_out[SP];i<firstOut[v+1];i++)
      {
         long u=adj[i];
         if((*dfs)[u]==-1)
         {
            (*dfs)[u]=nr; (*idfs)[nr++]=u; (*p)[u]=v;
            temp_vertex[SP+1]=u; temp_out[SP+1]=firstOut[u]; temp_out[SP]=i+1;
            descend=1; break;
         }
      }
      if(descend){SP++;continue;}
      SP--;
   }

   free(temp_vertex); free(temp_out);
}

void get_adj(long n, long m, long* edges, long** adj, long** firstOut)
{
   *adj = (long*)malloc(sizeof(long)*2*m);
   *firstOut = (long*)malloc(sizeof(long)*(n+1));
   for(long i=0;i<=n;i++){(*firstOut)[i]=0;}
   for(long i=0;i<m;i++){(*firstOut)[edges[2*i]+1]++; (*firstOut)[edges[2*i+1]+1]++;}
   for(long i=1;i<=n;i++){(*firstOut)[i]+=(*firstOut)[i-1];}
   long* nextOut = (long*)malloc(sizeof(long)*(n+1));
   for(long i=0;i<=n;i++){nextOut[i]=(*firstOut)[i];}
   for(long i=0;i<m;i++)
   {
      long x=edges[2*i]; long y=edges[2*i+1];
      (*adj)[nextOut[x]++]=y; (*adj)[nextOut[y]++]=x;
   } 
   free(nextOut);
}

void read_graph(char* filename, long* n, long** adj, long** firstOut)
{
   long m;
   FILE* fp = fopen(filename,"r");
   fscanf(fp,"%ld %ld",n,&m);
   long* edges = (long*)malloc(sizeof(long)*2*m);
   for(long i=0;i<m;i++){fscanf(fp,"%ld %ld",edges+2*i,edges+2*i+1);}
   fclose(fp);
   get_adj(*n,m,edges,adj,firstOut);
   free(edges);
}
