//Runs with ./get4eComponentsLinearNew <input_graph> <output_4components>

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

long get_4e_components_linear(long,long*,long*,long**);
long get_4e_components_connected_linear(long,long*,long*,long**);
long get_4e_components_2econnected_linear(long,long*,long*,long**);
long get_4e_components_3econnected_linear(long,long*,long*,long**,long*,long**);
long get_3cuts_contr(long,long*,long*,long**);
long get_3cuts_2tree(long,long*,long*,long**);
void get_l1l2_and_bcount(long,long*,long*,long*,long*,long*,long**,long**,long**);
void get_allM(long,long*,long*,long*,long*,long*,long*,long*,long**,long**,long**);
void get_2low(long,long*,long*,long*,long*,long*,long**,long**,long**,long**);

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
   long k=get_4e_components_linear(n,adj,firstOut,&C);
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

long get_4e_components_linear(long n, long* adj, long* firstOut, long** C)
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
         long k2 = get_4e_components_connected_linear(nC,adjC,firstOutC,&C2);
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

long get_4e_components_connected_linear(long n, long* adj, long* firstOut, long** C)
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
         long k3=get_4e_components_2econnected_linear(nC,adjC,firstOutC,&C3);
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

long get_4e_components_2econnected_linear(long n, long* adj, long* firstOut, long** C3)
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
         long k4=get_4e_components_3econnected_linear(nC[c],adjC[c],firstOutC[c],&C4,nCcuts+c,Ccuts+c);
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


long get_4e_components_3econnected_linear(long n, long* adj, long* firstOut, long** C, long* n3cuts, long** cuts3)
{
   long* cuts;
   long k = get_3cuts_contr(n,adj,firstOut,&cuts);
(*cuts3) = (long*)malloc(sizeof(long)*6*k);
*n3cuts=k; for(long i=0;i<6*k;i++){(*cuts3)[i]=cuts[i];}
   long* dfs; long* idfs; long* p;
   DFS(n,adj,firstOut,&dfs,&idfs,&p);
   long* ND = (long*)malloc(sizeof(long)*n);
   for(long i=0;i<n;i++){ND[i]=1;}
   for(long i=n-1;i>0;i--){long v=idfs[i];ND[p[v]]+=ND[v];}
   char* type = (char*)malloc(sizeof(char)*k);
   long* size = (long*)malloc(sizeof(long)*k);
   for(long i=0;i<k;i++)
   {
      long temp;
      if(dfs[cuts[6*i+0]]<dfs[cuts[6*i+1]]){temp=cuts[6*i+0];cuts[6*i+0]=cuts[6*i+1];cuts[6*i+1]=temp;}
      if(dfs[cuts[6*i+2]]<dfs[cuts[6*i+3]]){temp=cuts[6*i+2];cuts[6*i+2]=cuts[6*i+3];cuts[6*i+3]=temp;}
      if(dfs[cuts[6*i+4]]<dfs[cuts[6*i+5]]){temp=cuts[6*i+4];cuts[6*i+4]=cuts[6*i+5];cuts[6*i+5]=temp;}
      
      for(long t=0;t<3;t++)
      {
         for(long s=t+1;s<3;s++)
         {
            if(dfs[cuts[6*i+2*t]]<dfs[cuts[6*i+2*s]])
            {
               temp=cuts[6*i+2*t]; cuts[6*i+2*t]=cuts[6*i+2*s]; cuts[6*i+2*s]=temp;
               temp=cuts[6*i+2*t+1]; cuts[6*i+2*t+1]=cuts[6*i+2*s+1]; cuts[6*i+2*s+1]=temp;
            }  
         }
      }
   }
   for(long i=0;i<k;i++)
   {
      long x1=cuts[6*i+0]; long y1=cuts[6*i+1];  
      long x2=cuts[6*i+2]; long y2=cuts[6*i+3];  
      long x3=cuts[6*i+4]; long y3=cuts[6*i+5];
      type[i]=0;
      if(y1==p[x1]){type[i]++;}
      if(y2==p[x2]){type[i]++;}
      if(y3==p[x3]){type[i]++;}
      if(type[i]==3) 
      {
         if(x1==x2&&x1==x3){type[i]=0;}
         else if(x1==x2){type[i]=1;x2=x3;cuts[6*i+2]=cuts[6*i+4];}
         else if(x2==x3){type[i]=1;}
         else{type[i]=2;}
      }
      else if(type[i]==2)
      {
         if(y1==p[x1]&&y2==p[x2])
         {
            if(x1==x2){type[i]=0;}
            else{type[i]=1;}
         }
         else if(y2==p[x2]&&y3==p[x3])
         {
            if(x2==x3){type[i]=0;x1=x2;cuts[6*i]=cuts[6*i+2];}
            else{type[i]=1;x1=x2;x2=x3;cuts[6*i]=cuts[6*i+2];cuts[6*i+2]=cuts[6*i+4];}
         }
         else if(y1==p[x1]&&y3==p[x3])
         {
            type[i]=1; x2=x3; cuts[6*i+2]=cuts[6*i+4];  
         } 
      }
      else if(type[i]==1)
      {
         type[i]=0;
         if(y2==p[x2]){x1=x2;cuts[6*i]=cuts[6*i+2];}
         else if(y3==p[x3]){x1=x3;cuts[6*i]=cuts[6*i+4];}
      } 

      if(type[i]==0){size[i]=ND[x1];}
      else if(type[i]==1){size[i]=ND[x2]-ND[x1];}
      else if(type[i]==2)
      {
         if(dfs[x2]<dfs[x1] && dfs[x1]<dfs[x2]+ND[x2])
         {
            size[i]=ND[x3]-ND[x2]+ND[x1]; 
         }
         else
         {
            size[i]=ND[x3]-ND[x2]-ND[x1];
         }
      }
   }

   long* sizeStack = (long*)malloc(sizeof(long)*k);
   long* sizeFirst = (long*)malloc(sizeof(long)*n);
   long* sizeNext = (long*)malloc(sizeof(long)*k);
   for(long i=0;i<n;i++){sizeFirst[i]=-1;}
   long SP=0;
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

   long* Cstack = (long*)malloc(sizeof(long)*3*n);
   long* cactusParent = (long*)malloc(sizeof(long)*3*n);
   long* parentCut = (long*)malloc(sizeof(long)*3*n);
   char* marked = (char*)malloc(sizeof(char)*k);
   for(long i=0;i<k;i++){marked[i]=0;}
   long* phi = (long*)malloc(sizeof(long)*k);
   for(long i=0;i<k;i++){phi[i]=-1;}

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

   free(dfs); free(idfs); free(p); free(ND);
   free(type); free(size);
   free(cuts);
   return n_comp;
}

long get_3cuts_contr(long n, long* adj, long* firstOut, long** cuts)
{
   (*cuts) = (long*)malloc(sizeof(long)*n*12);
   long k=0; long k2;
   long* cuts_2tree;

   k2 = get_3cuts_2tree(n,adj,firstOut,&cuts_2tree);
   for(long i=0;i<k2;i++)
   {
      for(long t=0;t<6;t++){(*cuts)[6*k+t]=cuts_2tree[6*i+t];}k++;
   }
   free(cuts_2tree);
   
   long* dfs; long* idfs; long* p;
   DFS(n,adj,firstOut,&dfs,&idfs,&p);
   char* foundP = (char*)malloc(sizeof(char)*n);
   char* foundC = (char*)malloc(sizeof(char)*n);
   for(long i=0;i<n;i++){foundP[i]=0;foundC[i]=0;}
   long* C = (long*)malloc(sizeof(long)*n);
   for(long i=0;i<n;i++){C[i]=-1;}
   long* Q = (long*)malloc(sizeof(long)*n);
   long cIndx=0;
   for(long r=0;r<n;r++)
   {
      if(C[r]!=-1){continue;}
      long first=0; long last=0;
      Q[last++]=r; C[r]=cIndx;
      while(first!=last)
      {
         long v=Q[first++];
         for(long i=firstOut[v];i<firstOut[v+1];i++)
         {
            long u=adj[i];
            if(v==p[u]&&!foundC[u]){foundC[u]=1;continue;}
            if(u==p[v]&&!foundP[v]){foundP[v]=1;continue;}
            if(C[u]==-1)
            {
               Q[last++]=u; C[u]=cIndx;
            }
         }
      }
      cIndx++;
   }

   if(cIndx==1){free(dfs); free(idfs); free(p); free(foundP); free(C); free(Q); return k;}

   long* stackT_ = (long*)malloc(sizeof(long)*2*n);
   long* firstT_ = (long*)malloc(sizeof(long)*n);
   long* nextT_ = (long*)malloc(sizeof(long)*2*n);
   long* tree_edge_ = (long*)malloc(sizeof(long)*4*n);
   long SP=0;
   for(long i=0;i<n;i++){firstT_[i]=-1;}
   for(long v=1;v<n;v++)
   {
      long u=p[v];
      if(C[u]==C[v]){continue;}
      nextT_[SP]=firstT_[C[u]];
      firstT_[C[u]]=SP;
      stackT_[SP]=C[v];
      tree_edge_[2*SP]=u; tree_edge_[2*SP+1]=v; SP++;
      nextT_[SP]=firstT_[C[v]];
      firstT_[C[v]]=SP;
      stackT_[SP]=C[u];
      tree_edge_[2*SP]=u; tree_edge_[2*SP+1]=v; SP++;
   }

   long* stackT = (long*)malloc(sizeof(long)*2*n);
   long* firstT = (long*)malloc(sizeof(long)*n);
   long* nextT = (long*)malloc(sizeof(long)*2*n);
   long* tree_edge = (long*)malloc(sizeof(long)*4*n);
   SP=0;
   for(long i=0;i<n;i++){firstT[i]=-1;}
   for(long c=cIndx-1;c>=0;c--)
   {
      for(long s=firstT_[c];s!=-1;s=nextT_[s])
      {
         long d=stackT_[s];
         nextT[SP]=firstT[d];
         firstT[d]=SP;
         stackT[SP]=c;
         tree_edge[2*SP]=tree_edge_[2*s]; tree_edge[2*SP+1]=tree_edge_[2*s+1]; SP++;
      }
   }  

   long* adj_new = (long*)malloc(sizeof(long)*2*n);
   long* firstOut_new = (long*)malloc(sizeof(long)*(cIndx+1));
   for(long i=0;i<=cIndx;i++){firstOut_new[i]=0;}
   for(long v=1;v<n;v++)
   {
      long u=p[v];
      if(C[u]==C[v]){continue;}
      firstOut_new[C[u]+1]++; firstOut_new[C[v]+1]++;
   }
   long* currentOut = (long*)malloc(sizeof(long)*(cIndx+1));
   currentOut[0]=0;
   for(long i=1;i<=cIndx;i++){firstOut_new[i]+=firstOut_new[i-1]; currentOut[i]=firstOut_new[i];}
   for(long v=1;v<n;v++)
   {
      long u=p[v];
      if(C[u]==C[v]){continue;}
      adj_new[currentOut[C[v]]++]=C[u]; adj_new[currentOut[C[u]]++]=C[v];
   }

   long* cuts_new;
   k2 = get_3cuts_contr(cIndx,adj_new,firstOut_new,&cuts_new);

   long* stackC_ = (long*)malloc(sizeof(long)*cIndx*12);
   long* firstC_ = (long*)malloc(sizeof(long)*cIndx);
   long* nextC_ = (long*)malloc(sizeof(long)*cIndx*12);
   for(long i=0;i<cIndx;i++){firstC_[i]=-1;}
   long* cutIndx_ = (long*)malloc(sizeof(long)*cIndx*12);
   SP=0;
   for(long i=0;i<k2;i++)
   {
      for(long t=0;t<3;t++)
      {
         long x=cuts_new[6*i+2*t]; long y=cuts_new[6*i+2*t+1];
    
         nextC_[SP]=firstC_[x];
         firstC_[x]=SP;
         stackC_[SP]=y;
         cutIndx_[SP++]=i; 

         nextC_[SP]=firstC_[y];
         firstC_[y]=SP;
         stackC_[SP]=x;
         cutIndx_[SP++]=i; 
      }
   }

   long* stackC = (long*)malloc(sizeof(long)*cIndx*12);
   long* firstC = (long*)malloc(sizeof(long)*cIndx);
   long* nextC = (long*)malloc(sizeof(long)*cIndx*12);
   for(long i=0;i<cIndx;i++){firstC[i]=-1;}
   long* cutIndx = (long*)malloc(sizeof(long)*cIndx*12);
   SP=0;   
   for(long x=cIndx-1;x>=0;x--)
   {
      for(long indx=firstC_[x];indx!=-1;indx=nextC_[indx])
      {
         long y=stackC_[indx];
         nextC[SP]=firstC[y];
         firstC[y]=SP;
         stackC[SP]=x;
         cutIndx[SP++]=cutIndx_[indx];
      } 
   }

   long* currentEdge = (long*)malloc(sizeof(long)*k2);
   for(long i=0;i<k2;i++){currentEdge[i]=0;}
   for(long x=0;x<cIndx;x++)
   {
      for(long indx=firstC[x];indx!=-1;indx=nextC[indx])
      {
         long y=stackC[indx];
         if(y<x){continue;}
         long t=firstT[x];
         while(stackT[t]!=y){t=nextT[t];}
         firstT[x]=t;
         if(nextT[t]!=-1 && stackT[nextT[t]]==y){firstT[x]=nextT[t];}
         long cut_indx=cutIndx[indx];
         (*cuts)[6*k+6*cut_indx+2*currentEdge[cut_indx]]=tree_edge[2*t];
         (*cuts)[6*k+6*cut_indx+2*currentEdge[cut_indx]+1]=tree_edge[2*t+1];
         currentEdge[cut_indx]++;
      }
   }

   k+=k2;

   free(dfs); free(idfs); free(p); free(foundP);
   free(C); free(Q);
   free(stackT_); free(firstT_); free(nextT_); free(tree_edge_);
   free(stackT); free(firstT); free(nextT); free(tree_edge);
   free(adj_new); free(firstOut_new); free(currentOut);
   free(stackC_); free(firstC_); free(nextC_); free(cutIndx_);
   free(stackC); free(firstC); free(nextC); free(cutIndx);
   free(currentEdge);
   free(cuts_new);
   return k;
}

long get_3cuts_2tree(long n, long* adj, long* firstOut, long** cuts)
{
   (*cuts) = (long*)malloc(sizeof(long)*n*12);
   long* dfs; long* idfs; long* p;
   DFS(n,adj,firstOut,&dfs,&idfs,&p);
   long* low1; long* low1D; long* low2; long* low2D;
   get_2low(n,adj,firstOut,dfs,idfs,p,&low1,&low1D,&low2,&low2D);
   long* l1; long* l2; long* bcount;
   get_l1l2_and_bcount(n,adj,firstOut,dfs,idfs,p,&l1,&l2,&bcount);
   long* low1C; long* low2C;
   get_lowChildren(n,dfs,idfs,p,low1,&low1C,&low2C);
   long* M; long* nextM; 
   get_M(n,dfs,idfs,l1,low1,low1C,low2C,&M,&nextM);
   long* Ml; long* Mlow1; long* Mlow2;
   get_allM(n,dfs,idfs,l1,low1,M,low1C,low2C,&Ml,&Mlow1,&Mlow2);
  
   long* ND = (long*)malloc(sizeof(long)*n);
   for(long i=0;i<n;i++){ND[i]=1;}
   for(long i=n-1;i>0;i--){long v=idfs[i];ND[p[v]]+=ND[v];}
   long* prevM = (long*)malloc(sizeof(long)*n);
   for(long i=0;i<n;i++){prevM[i]=-1;}
   for(long i=0;i<n;i++){if(nextM[i]!=-1){prevM[nextM[i]]=i;}}

   long* low3C = (long*)malloc(sizeof(long)*n);
   for(long i=0;i<n;i++){low3C[i]=-1;}
   for(long c=1;c<n;c++)
   {
      long v=p[c];
      if(low1C[v]==c||low2C[v]==c){continue;}
      if(low3C[v]==-1||dfs[low1[c]]<dfs[low1[low3C[v]]]){low3C[v]=c;}
   }

   long k=0;
   long* currentVertex = (long*)malloc(sizeof(long)*n);
   for(long i=0;i<n;i++){currentVertex[i]=i;}
   for(long i=n-1;i>0;i--)
   {
      long v=idfs[i];
      long m=Ml[v];
      if(m==-1){continue;}
      long u=currentVertex[m];
      while(nextM[u]!=-1 && dfs[nextM[u]]>dfs[v]){u=nextM[u];}
      currentVertex[m]=u;
      if(bcount[v]==bcount[u]+1 && dfs[l2[M[v]]]>=dfs[v] && (low2C[M[v]]==-1 || dfs[low1[low2C[M[v]]]]>=dfs[v]))
      {
         (*cuts)[6*k+0]=u; (*cuts)[6*k+1]=p[u];
         (*cuts)[6*k+2]=v; (*cuts)[6*k+3]=p[v];
         (*cuts)[6*k+4]=M[v]; (*cuts)[6*k+5]=l1[M[v]];
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
      if(bcount[v]==bcount[u]+1 && dfs[low2[Mlow2[v]]]>=dfs[v] && (low3C[M[v]]==-1 || dfs[low1[low3C[M[v]]]]>=dfs[v]))
      {
         (*cuts)[6*k+0]=u; (*cuts)[6*k+1]=p[u];
         (*cuts)[6*k+2]=v; (*cuts)[6*k+3]=p[v];
         (*cuts)[6*k+4]=Mlow2[v]; (*cuts)[6*k+5]=l1[Mlow2[v]];
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
      if(bcount[v]==bcount[u]+1 && dfs[low2[Mlow1[v]]]>=dfs[v] && (low3C[M[v]]==-1 || dfs[low1[low3C[M[v]]]]>=dfs[v]))
      {
         (*cuts)[6*k+0]=u; (*cuts)[6*k+1]=p[u];
         (*cuts)[6*k+2]=v; (*cuts)[6*k+3]=p[v];
         (*cuts)[6*k+4]=Mlow1[v]; (*cuts)[6*k+5]=l1[Mlow1[v]];
         k++;
      }
   }

   //calculate the new parameters
   long* SDesc = (long*)malloc(sizeof(long)*n);
   long* SAnc = (long*)malloc(sizeof(long)*n);
   char* foundParent = (char*)malloc(sizeof(char)*n);
   for(long i=0;i<n;i++){SDesc[i]=0; SAnc[i]=0; foundParent[i]=0;}
   for(long d=n-1;d>0;d--)
   {
      long v=idfs[d];
      for(long i=firstOut[v];i<firstOut[v+1];i++)
      {
         long u=adj[i];
         if(u==p[v] && !foundParent[v]){foundParent[v]=1;continue;}
         if(dfs[u]>dfs[v]){continue;}
         SDesc[v]^=v; SAnc[v]^=u; 
         SDesc[u]^=v; SAnc[u]^=u;
      }
      SDesc[p[v]]^=SDesc[v]; SAnc[p[v]]^=SAnc[v];
   }
     
   for(long u=1;u<n;u++)
   {
      if(nextM[u]!=-1 && bcount[u]==bcount[nextM[u]]+1)
      {
         (*cuts)[6*k+0]=u; (*cuts)[6*k+1]=p[u];
         (*cuts)[6*k+2]=nextM[u]; (*cuts)[6*k+3]=p[nextM[u]];
         (*cuts)[6*k+4]=SDesc[u]^SDesc[nextM[u]]; (*cuts)[6*k+5]=SAnc[u]^SAnc[nextM[u]];
         k++;
      }
   }
   for(long i=0;i<n;i++){currentVertex[i]=i;}
   for(long i=n-1;i>0;i--)
   {
      long u=idfs[i];
      long m=Ml[u];
      if(m==-1 || Ml[u]==M[u]){continue;}
      long v=currentVertex[m];
      while(v!=-1 && dfs[v]>=dfs[u]){v=nextM[v];}
      currentVertex[m]=v;
      if(v!=-1 && bcount[u]==bcount[v]+1)
      {
         (*cuts)[6*k+0]=u; (*cuts)[6*k+1]=p[u];
         (*cuts)[6*k+2]=v; (*cuts)[6*k+3]=p[v];
         (*cuts)[6*k+4]=M[u]; (*cuts)[6*k+5]=l1[M[u]];
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
      if(v!=-1 && bcount[u]==bcount[v]+1)
      {
         (*cuts)[6*k+0]=u; (*cuts)[6*k+1]=p[u];
         (*cuts)[6*k+2]=v; (*cuts)[6*k+3]=p[v];
         (*cuts)[6*k+4]=Mlow2[u]; (*cuts)[6*k+5]=l1[Mlow2[u]];
         k++;
      }
   }
  
   for(long v=1;v<n;v++)
   {
      if(bcount[v]==2)
      {
         (*cuts)[6*k+0]=v; (*cuts)[6*k+1]=p[v];
         (*cuts)[6*k+2]=low1D[v]; (*cuts)[6*k+3]=low1[v];
         (*cuts)[6*k+4]=low2D[v]; (*cuts)[6*k+5]=low2[v];
         k++;   
      }
   }
 
   free(dfs); free(idfs); free(p); free(ND);
   free(low1); free(low1D); free(low2); free(low2D);  
   free(l1); free(l2); free(bcount);
   free(low1C); free(low2C); free(low3C);
   free(M); free(nextM); free(prevM);
   free(Ml); free(Mlow1); free(Mlow2); free(currentVertex);
   free(SDesc); free(SAnc); free(foundParent);
   return k;
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



void get_l1l2_and_bcount(long n, long* adj, long* firstOut, long* dfs, long* idfs, long* p, long** l1, long** l2, long** bcount)
{
   (*l1) = (long*)malloc(sizeof(long)*n);
   (*l2) = (long*)malloc(sizeof(long)*n);
   (*bcount) = (long*)malloc(sizeof(long)*n);
   for(long i=0;i<n;i++){(*l1)[i]=i; (*l2)[i]=i; (*bcount)[i]=0;}
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
               if(dfs[y]<=dfs[(*l1)[x]]){(*l2)[x]=(*l1)[x]; (*l1)[x]=y;}
               else if(dfs[y]<dfs[(*l2)[x]]){(*l2)[x]=y;}
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
