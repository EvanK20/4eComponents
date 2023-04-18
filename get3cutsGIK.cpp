//Runs with ./get3cutsGIK <input_graph> <output_3cuts>

#include <stdio.h>
#include <stdlib.h>
#include <chrono>

void get_adj(long,long,long*,long**,long**);
void read_graph(char*,long*,long**,long**);

void DFS(long,long*,long*,long**,long**,long**);
void get_2low(long,long*,long*,long*,long*,long*,long**,long**,long**,long**);
void get_l_and_bcount(long,long*,long*,long*,long*,long*,long**,long**);
void get_lowChildren(long,long*,long*,long*,long*,long**,long**);
void get_M(long,long*,long*,long*,long*,long*,long*,long**,long**);
void get_allM(long,long*,long*,long*,long*,long*,long*,long*,long**,long**,long**);
void get_high(long,long*,long*,long*,long*,long*,long**,long**);
long find(long*,long);
void unite(long*,long*,long,long);

long get_3cuts_high(long,long*,long*,long**);

using namespace std::chrono;

int main(int n_args, char** args)
{
   long n; long* adj; long* firstOut;
   read_graph(args[1],&n,&adj,&firstOut);
high_resolution_clock::time_point t1 = high_resolution_clock::now();
   long* cuts;
   long k=get_3cuts_high(n,adj,firstOut,&cuts);
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

long get_3cuts_high(long n, long* adj, long* firstOut, long** cuts)
{
   (*cuts) = (long*)malloc(sizeof(long)*n*12);
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

   long k=0;

   //one tree-edge
   for(long v=1;v<n;v++)
   {
      if(bcount[v]==2)
      {
         (*cuts)[6*k+0]=v; (*cuts)[6*k+1]=p[v];
         (*cuts)[6*k+2]=low[v]; (*cuts)[6*k+3]=lowD[v];
         (*cuts)[6*k+4]=low2[v]; (*cuts)[6*k+5]=low2D[v];
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
         (*cuts)[6*k+0]=u; (*cuts)[6*k+1]=p[u];
         (*cuts)[6*k+2]=v; (*cuts)[6*k+3]=p[v];
         (*cuts)[6*k+4]=M[v]; (*cuts)[6*k+5]=l[M[v]];
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
         (*cuts)[6*k+0]=u; (*cuts)[6*k+1]=p[u];
         (*cuts)[6*k+2]=v; (*cuts)[6*k+3]=p[v];
         (*cuts)[6*k+4]=Mlow2[v]; (*cuts)[6*k+5]=l[Mlow2[v]];
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
         (*cuts)[6*k+0]=u; (*cuts)[6*k+1]=p[u];
         (*cuts)[6*k+2]=v; (*cuts)[6*k+3]=p[v];
         (*cuts)[6*k+4]=Mlow1[v]; (*cuts)[6*k+5]=l[Mlow1[v]];
         k++;
      }
   }

   //lower case
   for(long u=1;u<n;u++)
   {
      if(nextM[u]!=-1 && bcount[nextM[u]]==bcount[u]-1)
      {
         (*cuts)[6*k+0]=u; (*cuts)[6*k+1]=p[u];
         (*cuts)[6*k+2]=nextM[u]; (*cuts)[6*k+3]=p[nextM[u]];
         (*cuts)[6*k+4]=highD[u]; (*cuts)[6*k+5]=high[u];
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
         (*cuts)[6*k+0]=u; (*cuts)[6*k+1]=p[u];
         (*cuts)[6*k+2]=v; (*cuts)[6*k+3]=p[v];
         (*cuts)[6*k+4]=highD[u]; (*cuts)[6*k+5]=high[u];
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
         (*cuts)[6*k+0]=u; (*cuts)[6*k+1]=p[u];
         (*cuts)[6*k+2]=v; (*cuts)[6*k+3]=p[v];
         (*cuts)[6*k+4]=highD[u]; (*cuts)[6*k+5]=high[u];
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
         (*cuts)[6*k+0]=u; (*cuts)[6*k+1]=p[u];
         (*cuts)[6*k+2]=v; (*cuts)[6*k+3]=p[v];
         (*cuts)[6*k+4]=w; (*cuts)[6*k+5]=p[w];
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
         (*cuts)[6*k+0]=u; (*cuts)[6*k+1]=p[u];
         (*cuts)[6*k+2]=v; (*cuts)[6*k+3]=p[v];
         (*cuts)[6*k+4]=w; (*cuts)[6*k+5]=p[w];
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
                  (*cuts)[6*k+0]=u; (*cuts)[6*k+1]=p[u];
                  (*cuts)[6*k+2]=v; (*cuts)[6*k+3]=p[v];
                  (*cuts)[6*k+4]=w; (*cuts)[6*k+5]=p[w];
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
            (*cuts)[6*k+0]=u; (*cuts)[6*k+1]=p[u];
            (*cuts)[6*k+2]=v; (*cuts)[6*k+3]=p[v];
            (*cuts)[6*k+4]=w; (*cuts)[6*k+5]=p[w];
            k++;  
         }
         indx=nextU[indx];
      }
   }

   free(dfs); free(idfs); free(p);
   free(low); free(lowD); free(low2); free(low2D);
   free(low1C); free(low2C); free(l); free(bcount);
   free(M); free(nextM); free(high); free(highD);
   free(Ml); free(Mlow1); free(Mlow2); free(currentVertex); 
   free(lastM); free(A); free(ND); 
   free(Hstack); free(firstH); free(nextH); free(stack);
   free(Ustack); free(firstU); free(nextU);
   free(lowestW);
   return k;
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
