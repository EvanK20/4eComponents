//Runs with ./sparsify <input_graph> <output_graph>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <chrono>

long find(long*,long);
void unite(long*,long*,long,long);

using namespace std::chrono;

int main(int nArgs, char** args)
{
   srand(time(NULL));
   //source filename, target filename
   long n; long m;
   FILE* fp = fopen(args[1],"r");
   fscanf(fp,"%ld %ld",&n,&m);
   long* graph_edges = (long*)malloc(sizeof(long)*2*m);
   for(long i=0;i<m;i++){fscanf(fp,"%ld %ld",graph_edges+2*i,graph_edges+2*i+1);}
   fclose(fp);
high_resolution_clock::time_point t1 = high_resolution_clock::now();
   long* edges = (long*)malloc(sizeof(long)*8*n);
   char* is_available = (char*)malloc(sizeof(char)*m);
   for(long i=0;i<m;i++){is_available[i]=1;}
   long* ufparent = (long*)malloc(sizeof(long)*n);
   long* ufsize = (long*)malloc(sizeof(long)*n);
   long edgeIndx=0;
   for(long t=0;t<4;t++)
   {
      for(long x=0;x<n;x++){ufparent[x]=x; ufsize[x]=1;}
      for(long i=0;i<m;i++)
      {
         if(!is_available[i]){continue;}
         long x=graph_edges[2*i]; long y=graph_edges[2*i+1];
         if(find(ufparent,x)!=find(ufparent,y))
         {
            unite(ufparent,ufsize,x,y);
            edges[2*edgeIndx]=x; edges[2*edgeIndx+1]=y; edgeIndx++;
            is_available[i]=0;
         } 
      }
   }
high_resolution_clock::time_point t2 = high_resolution_clock::now();
duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
   fp = fopen(args[2],"w");
   fprintf(fp,"%ld %ld\n",n,edgeIndx);
   for(long i=0;i<edgeIndx;i++){fprintf(fp,"%ld %ld\n",edges[2*i],edges[2*i+1]);}
   fclose(fp);

   free(edges); free(graph_edges); free(is_available);
   free(ufparent); free(ufsize);
   printf("time: %f\n",time_span.count());
   return 0;
}

void unite(long* p, long* size, long x, long y)
{
   long r1=find(p,x);
   long r2=find(p,y);
   if(r1==r2){return;}
   if(size[r1]<size[r2]){p[r1]=r2;size[r2]+=size[r1];}
   else{p[r2]=r1;size[r1]+=size[r2];}
}

long find(long* p, long x)
{
   long r=x;
   while(p[r]!=r){r=p[r];}
   while(x!=r){long next=p[x];p[x]=r;x=next;}
   return r;
}
