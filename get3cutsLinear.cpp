//Runs with ./get3cutsLinear <input_graph> <output_3cuts>

#include <stdio.h>
#include <stdlib.h>
#include <chrono>

void get_adj(long,long,long*,long**,long**);
void read_graph(char*,long*,long**,long**);

void DFS(long,long*,long*,long**,long**,long**);
void get_2low(long,long*,long*,long*,long*,long*,long**,long**,long**,long**);
void get_l1l2_and_bcount(long,long*,long*,long*,long*,long*,long**,long**,long**);
void get_lowChildren(long,long*,long*,long*,long*,long**,long**);
void get_M(long,long*,long*,long*,long*,long*,long*,long**,long**);
void get_allM(long,long*,long*,long*,long*,long*,long*,long*,long**,long**,long**);
void get_lowM(long,long*,long*,long*,long*,long*,long*,long*,long*,long*,long**,long**);
void sortAdjInc(long,long*,long*,long*);

long get_3cuts_linear(long,long*,long*,long**);
long get_3cuts_2tree(long,long*,long*,long**);

using namespace std::chrono;

int main(int n_args, char** args)
{
   long n; long* adj; long* firstOut;
   read_graph(args[1],&n,&adj,&firstOut);
high_resolution_clock::time_point t1 = high_resolution_clock::now();
   long* cuts;
   long k=get_3cuts_linear(n,adj,firstOut,&cuts);
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

long get_3cuts_linear(long n, long* adj, long* firstOut, long** cuts)
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
   k2 = get_3cuts_linear(cIndx,adj_new,firstOut_new,&cuts_new);

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

   long* lowM; long* lowMD;
   get_lowM(n,adj,firstOut,dfs,idfs,p,ND,low1C,M,prevM,&lowM,&lowMD);

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
  
   for(long u=1;u<n;u++)
   {
      if(nextM[u]!=-1 && bcount[u]==bcount[nextM[u]]+1)
      {
         (*cuts)[6*k+0]=u; (*cuts)[6*k+1]=p[u];
         (*cuts)[6*k+2]=nextM[u]; (*cuts)[6*k+3]=p[nextM[u]];
         (*cuts)[6*k+4]=lowMD[u]; (*cuts)[6*k+5]=lowM[u];
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
   free(low1); free(low1D); free(low2); free(low2D);  free(lowM); free(lowMD);
   free(l1); free(l2); free(bcount);
   free(low1C); free(low2C); free(low3C);
   free(M); free(nextM); free(prevM);
   free(Ml); free(Mlow1); free(Mlow2); free(currentVertex);
   return k;
}

void get_lowM(long n, long* adj, long* firstOut, long* dfs, long* idfs, long* p, long* ND, long* low1C, long* M, long* prevM, long** lowM, long** lowMD)
{
   sortAdjInc(n,adj,firstOut,idfs);
   (*lowMD) = (long*)malloc(sizeof(long)*n);
   (*lowM) = (long*)malloc(sizeof(long)*n);
   for(long i=0;i<n;i++){(*lowMD)[i]=-1; (*lowM)[i]=-1;}
   long* currentOut = (long*)malloc(sizeof(long)*n);
   char* foundC = (char*)malloc(sizeof(char)*n);
   for(long i=0;i<n;i++){foundC[i]=0;}
   for(long i=0;i<n;i++){currentOut[i]=firstOut[i];}
   for(long i=n-1;i>0;i--)
   {
      long v=idfs[i];
      if(prevM[v]==-1){continue;}
      long u=prevM[v];
      long y=v;
      while((*lowM)[u]==-1)
      {
         while(currentOut[y]!=firstOut[y+1])
         {
            long x=adj[currentOut[y]];
            if(dfs[x]<dfs[y]){currentOut[y]++;continue;}
            if(y==p[x]&&foundC[x]==0){foundC[x]=1;currentOut[y]++;continue;}
            if(dfs[x]<dfs[M[u]]){currentOut[y]++;continue;}
            if(dfs[x]<dfs[M[u]]+ND[M[u]])
            {
               (*lowMD)[u]=x; (*lowM)[u]=y;
            }
            break;
         } 
         if((*lowM)[u]==-1)
         {
            if(prevM[low1C[y]]==-1){y=low1C[y];}
            else{y=(*lowM)[prevM[low1C[y]]];}
         }
      }
   }
   free(currentOut); free(foundC);
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

void sortAdjInc(long n, long* adj, long* firstOut, long* idfs)
{
   long* adj_copy = (long*)malloc(sizeof(long)*firstOut[n]);
   long* currentOut = (long*)malloc(sizeof(long)*(n+1));
   for(long i=0;i<firstOut[n];i++){adj_copy[i]=adj[i];}
   for(long i=0;i<=n;i++){currentOut[i]=firstOut[i];}
   for(long i=0;i<n;i++)
   {
      long x=idfs[i];
      for(long t=firstOut[x];t<firstOut[x+1];t++)
      {
         long y=adj_copy[t];
         adj[currentOut[y]++]=x;
      }
   }
   free(adj_copy); free(currentOut);
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
