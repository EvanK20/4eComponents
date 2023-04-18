//Runs with ./get3cutsNRSSsort <input_graph> <output_3cuts>

#include <stdio.h>
#include <stdlib.h>
#include <chrono>

void get_adj(long,long,long*,long**,long**);
void read_graph(char*,long*,long**,long**);

void DFS(long,long*,long*,long**,long**,long**);
void get_2low(long,long*,long*,long*,long*,long*,long**,long**,long**,long**);
void get_bcount(long,long*,long*,long*,long*,long*,long**);
void get_LandR(long,long*,long*,long*,long*,long*,long**,long**,long**,long**,long**,long**);
void get_high_sort(long,long*,long*,long*,long*,long*,long**,long**,long**);
long find(long*,long);
void unite(long*,long*,long,long);

long get_3cuts_rand(long,long*,long*,long**);
long get_3cuts_2tree_rand(long,long*,long*,long**);

using namespace std::chrono;

int main(int n_args, char** args)
{
   long n; long* adj; long* firstOut;
   read_graph(args[1],&n,&adj,&firstOut);
high_resolution_clock::time_point t1 = high_resolution_clock::now();
   long* cuts;
   long k=get_3cuts_rand(n,adj,firstOut,&cuts);
high_resolution_clock::time_point t2 = high_resolution_clock::now();
   FILE* fp = fopen(args[2],"w");
   fprintf(fp,"%ld\n",k);
   for(long i=0;i<k;i++){fprintf(fp,"%ld %ld %ld %ld %ld %ld\n",cuts[6*i+0],cuts[6*i+1],cuts[6*i+2],cuts[6*i+3],cuts[6*i+4],cuts[6*i+5]);}
duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
   fprintf(fp,"%f\n",time_span.count());
   fclose(fp);
   printf("%f\n",time_span.count());;
   free(adj); free(firstOut); free(cuts);
   return 0;
}

long get_3cuts_rand(long n, long* adj, long* firstOut, long** cuts)
{
   (*cuts) = (long*)malloc(sizeof(long)*n*12);
   long k=0; long k2;
   long* cuts_2tree;

   k2 = get_3cuts_2tree_rand(n,adj,firstOut,&cuts_2tree);
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
   k2 = get_3cuts_rand(cIndx,adj_new,firstOut_new,&cuts_new);

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


long get_3cuts_2tree_rand(long n, long* adj, long* firstOut, long** cuts)
{
   (*cuts) = (long*)malloc(sizeof(long)*12*n);

   long* dfs; long* idfs; long* p;
   DFS(n,adj,firstOut,&dfs,&idfs,&p);

   long* high; long* highD; long* highIndx;

   get_high_sort(n,adj,firstOut,dfs,idfs,p,&high,&highD,&highIndx); 
   long* L; long* R; long* Lhead; long* Rhead; long* Lindx; long* Rindx;
   get_LandR(n,adj,firstOut,dfs,idfs,p,&L,&R,&Lhead,&Rhead,&Lindx,&Rindx); 

   long log2n=1; long tempN=n; while(tempN!=1){log2n++;tempN/=2;}
   long n_bits=3*log2n;
   long w=1;
   while(8*sizeof(long)*w<n_bits){w++;}

   long* hash = (long*)malloc(sizeof(long)*w*firstOut[n]);
   char* foundP = (char*)malloc(sizeof(char)*n); 
   for(long i=0;i<n;i++){foundP[i]=0;}
   for(long d=n-1;d>0;d--)
   {
      long x=idfs[d];
      for(long i=firstOut[x];i<firstOut[x+1];i++)
      {
         long y=adj[i];
         if(dfs[y]>dfs[x]){continue;}
         if(y==p[x]&&!foundP[x]){foundP[x]=1;continue;}
         for(long t=0;t<w;t++){hash[i*w+t]=rand(); hash[i*w+t]+=((long)rand())<<(sizeof(int)*8);}
      }
   }
   long* H = (long*)malloc(sizeof(long)*w*n);
   for(long i=0;i<w*n;i++){H[i]=0;}
   char* foundC = (char*)malloc(sizeof(char)*n);
   for(long i=0;i<n;i++){foundP[i]=0; foundC[i]=0;}
   for(long d=n-1;d>0;d--)
   {
      long x=idfs[d];
      for(long i=firstOut[x];i<firstOut[x+1];i++)
      {
         long y=adj[i];
         if(y==p[x]&&!foundP[x]){foundP[x]=1;continue;}
         if(dfs[y]<dfs[x])
         {
            for(long t=0;t<w;t++){H[x*w+t]^=hash[i*w+t]; H[y*w+t]^=hash[i*w+t];}
         }
         else if(x==p[y]&&!foundC[y])
         {
            foundC[y]=1;
            for(long t=0;t<w;t++){H[x*w+t]^=H[y*w+t];}
         } 
      }
   }

   long N=1; for(long i=1;i<=log2n;i++){N*=2;}
   long** table = (long**)malloc(sizeof(long*)*N);
   long* capacity = (long*)malloc(sizeof(long)*N);
   long* taken = (long*)malloc(sizeof(long)*N);
   for(long i=0;i<N;i++){capacity[i]=1;}
   for(long i=0;i<N;i++){table[i]=(long*)malloc(sizeof(long)*w*capacity[i]);}
   long** corresponds = (long**)malloc(sizeof(long*)*N);
   for(long i=0;i<N;i++){corresponds[i]=(long*)malloc(sizeof(long)*capacity[i]);}
    
   for(long i=0;i<N;i++){taken[i]=0;}
   for(long u=1;u<n;u++)
   {
      long h=(ulong)(H[u*w+0])%N;
      if(taken[h]==capacity[h])
      {
         capacity[h]*=2;
         table[h]=(long*)realloc(table[h],sizeof(long)*w*capacity[h]);
         corresponds[h]=(long*)realloc(corresponds[h],sizeof(long)*capacity[h]);
      }
      for(long t=0;t<w;t++){table[h][taken[h]*w+t]=H[u*w+t];}
      corresponds[h][taken[h]++]=u;
   }

   long cutIndx=0;
   long* low1; long* low1D; long* low2; long* low2D;
   get_2low(n,adj,firstOut,dfs,idfs,p,&low1,&low1D,&low2,&low2D);
   long* bcount;
   get_bcount(n,adj,firstOut,dfs,idfs,p,&bcount);
   for(long v=1;v<n;v++)
   {
      if(bcount[v]==2)
      {
         (*cuts)[6*cutIndx+0]=v; (*cuts)[6*cutIndx+1]=p[v];
         (*cuts)[6*cutIndx+2]=low1D[v]; (*cuts)[6*cutIndx+3]=low1[v];
         (*cuts)[6*cutIndx+4]=low2D[v]; (*cuts)[6*cutIndx+5]=low2[v];
         cutIndx++;
      }
   }

   free(bcount); free(low1); free(low1D); free(low2); free(low2D);

   long* temp_hash = (long*)malloc(sizeof(long)*w);
   //L[v]
   for(long v=1;v<n;v++)
   {
      for(long t=0;t<w;t++){temp_hash[t]=H[v*w+t]^hash[Lindx[v]*w+t];}
      long h=(ulong)(temp_hash[0])%N;
      long indx=-1;
      for(long i=0;i<taken[h];i++)
      {
         char found=1;
         for(long t=0;t<w;t++){if(table[h][i*w+t]!=temp_hash[t]){found=0;break;}}
         if(found){indx=i;break;}
      }
      if(indx==-1){continue;}
      long u=corresponds[h][indx];
      if(dfs[v]<dfs[u])
      {
         (*cuts)[6*cutIndx+0]=v; (*cuts)[6*cutIndx+1]=p[v];
         (*cuts)[6*cutIndx+2]=L[v]; (*cuts)[6*cutIndx+3]=Lhead[v];
         (*cuts)[6*cutIndx+4]=u; (*cuts)[6*cutIndx+5]=p[u]; cutIndx++;
      }      
   }
   //R[v]
   for(long v=1;v<n;v++)
   {
      for(long t=0;t<w;t++){temp_hash[t]=H[v*w+t]^hash[Rindx[v]*w+t];}
      long h=(ulong)(temp_hash[0])%N;
      long indx=-1;
      for(long i=0;i<taken[h];i++)
      {
         char found=1;
         for(long t=0;t<w;t++){if(table[h][i*w+t]!=temp_hash[t]){found=0;break;}}
         if(found){indx=i;break;}
      }
      if(indx==-1){continue;}
      long u=corresponds[h][indx];
      if(dfs[v]<dfs[u])
      {
         (*cuts)[6*cutIndx+0]=v; (*cuts)[6*cutIndx+1]=p[v];
         (*cuts)[6*cutIndx+2]=R[v]; (*cuts)[6*cutIndx+3]=Rhead[v];
         (*cuts)[6*cutIndx+4]=u; (*cuts)[6*cutIndx+5]=p[u]; cutIndx++;
      }      
   }
   //highD[v]
   for(long v=1;v<n;v++)
   {
      for(long t=0;t<w;t++){temp_hash[t]=H[v*w+t]^hash[highIndx[v]*w+t];}
      long h=(ulong)(temp_hash[0])%N;
      long indx=-1;
      for(long i=0;i<taken[h];i++)
      {
         char found=1;
         for(long t=0;t<w;t++){if(table[h][i*w+t]!=temp_hash[t]){found=0;break;}}
         if(found){indx=i;break;}
      }
      if(indx==-1){continue;}
      long u=corresponds[h][indx];
      if(dfs[v]>dfs[u])
      {
         (*cuts)[6*cutIndx+0]=v; (*cuts)[6*cutIndx+1]=p[v];
         (*cuts)[6*cutIndx+2]=highD[v]; (*cuts)[6*cutIndx+3]=high[v];
         (*cuts)[6*cutIndx+4]=u; (*cuts)[6*cutIndx+5]=p[u]; cutIndx++;
      }      
   }

   for(long i=0;i<N;i++){free(table[i]); free(corresponds[i]);}
   free(table); free(capacity); free(taken); free(corresponds);
   free(hash); free(H); free(temp_hash);

   free(dfs); free(idfs); free(p);
   free(high); free(highD); free(highIndx);
   free(L); free(R); free(Lhead); free(Rhead); free(Lindx); free(Rindx);
   free(foundP); free(foundC);
   return cutIndx;
}

void get_high_sort(long n, long* adj, long* firstOut, long* dfs, long* idfs, long* p, long** high, long** highD, long** highIndx)
{
   long* backEdges = (long*)malloc(sizeof(long)*(firstOut[n]-2*n+2));
   long* index = (long*)malloc(sizeof(long)*(firstOut[n]/2-n+1));
   long bp=0;  
   char* foundP = (char*)malloc(sizeof(char)*n);
   for(long i=0;i<n;i++){foundP[i]=0;}
   for(long x=0;x<n;x++)
   {
      for(long i=firstOut[x];i<firstOut[x+1];i++)
      {
         long y=adj[i];
         if(y==p[x]&&!foundP[x]){foundP[x]=1;continue;}
         if(dfs[y]>dfs[x]){continue;}
         backEdges[2*bp]=x; backEdges[2*bp+1]=y; index[bp++]=i;
      }
   }
   
   //sort the back-edges decreasingly w.r.t. their lower end
   long* backEdgesStack = (long*)malloc(sizeof(long)*2*bp);
   long* backEdgesFirst = (long*)malloc(sizeof(long)*n);
   long* backEdgesNext = (long*)malloc(sizeof(long)*2*bp);
   for(long i=0;i<n;i++){backEdgesFirst[i]=-1;}
   long SP=0;
   for(long i=0;i<bp;i++)
   {
      long d=dfs[backEdges[2*i+1]];
      backEdgesNext[SP]=backEdgesFirst[d]; backEdgesFirst[d]=SP; backEdgesStack[SP++]=i;
   }

   long* ufparent = (long*)malloc(sizeof(long)*n);
   long* ufsize = (long*)malloc(sizeof(long)*n);
   long* lowest = (long*)malloc(sizeof(long)*n);
   for(long i=0;i<n;i++){ufparent[i]=i; ufsize[i]=1; lowest[i]=i;}
   (*highIndx) = (long*)malloc(sizeof(long)*n);
   (*high) = (long*)malloc(sizeof(long)*n);
   (*highD) = (long*)malloc(sizeof(long)*n);
   for(long d=n-1;d>=0;d--)
   {
      for(long i=backEdgesFirst[d];i!=-1;i=backEdgesNext[i])
      {
         long t=backEdgesStack[i];
         long x=backEdges[2*t]; long y=backEdges[2*t+1];
         long u=lowest[find(ufparent,x)];
         while(u!=y)
         {
            (*high)[u]=y; (*highD)[u]=x; (*highIndx)[u]=index[t];
            long next=lowest[find(ufparent,p[u])];
            unite(ufparent,ufsize,u,p[u]);
            lowest[find(ufparent,p[u])]=y;
            u=next;
         }
      }
   }

   free(backEdges); free(index); free(foundP);
   free(backEdgesStack); free(backEdgesFirst); free(backEdgesNext);
   free(ufparent); free(ufsize); free(lowest);
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

void get_LandR(long n, long* adj, long* firstOut, long* dfs, long* idfs, long* p, long** L, long** R, long** Lhead, long** Rhead, long** Lindx, long** Rindx)
{
   (*L) = (long*)malloc(sizeof(long)*n);
   (*R) = (long*)malloc(sizeof(long)*n);
   (*Lhead) = (long*)malloc(sizeof(long)*n);
   (*Rhead) = (long*)malloc(sizeof(long)*n);
   (*Lindx) = (long*)malloc(sizeof(long)*firstOut[n]);
   (*Rindx) = (long*)malloc(sizeof(long)*firstOut[n]);

   long* next = (long*)malloc(sizeof(long)*n);
   long* prev = (long*)malloc(sizeof(long)*n);
   for(long i=0;i<n;i++){next[i]=i+1; prev[i]=i-1;}
   long* lastChild = (long*)malloc(sizeof(long)*n);
   for(long i=0;i<n;i++){lastChild[i]=-1;}
   for(long i=n-1;i>0;i--)
   {
      long c=idfs[i];
      if(lastChild[p[c]]==-1){lastChild[p[c]]=c;}
   }
   char* foundP = (char*)malloc(sizeof(char)*n);
   for(long i=0;i<n;i++){foundP[i]=0;}
   long* currentOut = (long*)malloc(sizeof(long)*n);
   for(long i=0;i<n;i++){currentOut[i]=firstOut[i];}
   for(long i=n-1;i>0;i--)
   {
      long v=idfs[i];
      long x=v;
      char found=0;
      while(!found)
      {
         for(long t=currentOut[x];t<firstOut[x+1];t++)
         {
            long y=adj[t];
            if(dfs[y]>dfs[x]){currentOut[x]=t+1;continue;}
            if(y==p[x]&&!foundP[x]){foundP[x]=1;currentOut[x]=t+1;continue;}
            if(dfs[y]<dfs[v])
            {
               (*L)[v]=x; (*Lhead)[v]=y; (*Lindx)[v]=t;
               found=1; break;
            }
            currentOut[x]=t+1;
         }
         if(!found)
         {
            long d=dfs[x];
            long next_vertex=idfs[next[d]];
            if(prev[d]!=-1){next[prev[d]]=next[d];}
            if(next[d]!=n){prev[next[d]]=prev[d];}
            x=next_vertex;
         }
      }
      if(lastChild[v]==-1){x=v;}
      else{x=(*R)[lastChild[v]];}
      found=0;
      while(!found)
      {
         for(long t=currentOut[x];t<firstOut[x+1];t++)
         {
            long y=adj[t];
            if(dfs[y]>dfs[x]){currentOut[x]=t+1;continue;}
            if(y==p[x]&&!foundP[x]){foundP[x]=1;currentOut[x]=t+1;continue;}
            if(dfs[y]<dfs[v])
            {
               (*R)[v]=x; (*Rhead)[v]=y; (*Rindx)[v]=t;
               found=1; break;
            }
            currentOut[x]=t+1;
         }
         if(!found)
         {
            long d=dfs[x];
            long next_vertex=idfs[prev[d]];
            if(prev[d]!=-1){next[prev[d]]=next[d];}
            if(next[d]!=n){prev[next[d]]=prev[d];}
            x=next_vertex;
         }
      }
   }
   free(next); free(prev); free(lastChild); free(foundP); free(currentOut);
}

void get_bcount(long n, long* adj, long* firstOut, long* dfs, long* idfs, long* p, long** bcount)
{
   (*bcount) = (long*)malloc(sizeof(long)*n);
   for(long i=0;i<n;i++){(*bcount)[i]=0;}
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
