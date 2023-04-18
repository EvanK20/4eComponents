//Runs with ./make3eConnected <input_graph> <output_graph>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void get_adj(long,long,long*,long**,long**);

void make_2e_connected(long,long*,long*,long**,long**);
long make_3e_connected(long,long*,long*,long**);

void DFS(long,long*,long*,long**,long**,long**);
void get_low(long,long*,long*,long*,long*,long*,long**);
void get_l_and_bcount(long,long*,long*,long*,long*,long*,long**,long**);
void get_lowChildren(long,long*,long*,long*,long*,long**,long**);
void get_M(long,long*,long*,long*,long*,long*,long*,long**,long**);

int main(int n_args, char** args)
{
   srand(time(NULL));
   FILE* fp = fopen(args[1],"r");
   long n; long m;
   fscanf(fp,"%ld %ld",&n,&m);
   long* initial_edges = (long*)malloc(sizeof(long)*m*2);
   for(long i=0;i<m;i++){fscanf(fp,"%ld %ld",initial_edges+2*i,initial_edges+2*i+1);}
   fclose(fp);

   long* adj; long* firstOut;
   get_adj(n,m,initial_edges,&adj,&firstOut);
   free(initial_edges);

   long* adj2; long* firstOut2;
   make_2e_connected(n,adj,firstOut,&adj2,&firstOut2);
   free(adj); free(firstOut);

   long* edges;
   m = make_3e_connected(n,adj2,firstOut2,&edges);

   fp = fopen(args[2],"w");
   fprintf(fp,"%ld %ld\n",n,m);
   for(long i=0;i<m;i++){fprintf(fp,"%ld %ld\n",edges[2*i],edges[2*i+1]);}
   fclose(fp);

   free(adj2); free(firstOut2);
   free(edges);
   return 0;
}

long make_3e_connected(long n, long* adj, long* firstOut, long** edges)
{
   long* dfs; long* idfs; long* p;
   DFS(n,adj,firstOut,&dfs,&idfs,&p);
   long* low; long* l; long* bcount;
   get_low(n,adj,firstOut,dfs,idfs,p,&low);
   get_l_and_bcount(n,adj,firstOut,dfs,idfs,p,&l,&bcount);
   long* low1C; long* low2C;
   get_lowChildren(n,dfs,idfs,p,low,&low1C,&low2C);
   long* M; long* nextM;
   get_M(n,dfs,idfs,l,low,low1C,low2C,&M,&nextM);

   long* vEdgeStack = (long*)malloc(sizeof(long)*2*n);
   long* vEdgeFirst = (long*)malloc(sizeof(long)*n);
   long* vEdgeNext = (long*)malloc(sizeof(long)*2*n);
   for(long i=0;i<n;i++){vEdgeFirst[i]=-1;}
   char* isCutEdge = (char*)malloc(sizeof(char)*n);
   char* isCutEdge2 = (char*)malloc(sizeof(char)*n);
   for(long i=0;i<n;i++){isCutEdge[i]=0; isCutEdge2[i]=0;}
   long SP=0;
   long* corresponding = (long*)malloc(sizeof(long)*n);
   long* corresponding2 = (long*)malloc(sizeof(long)*n);

   for(long m=1;m<n;m++)
   {
      if(M[m]!=m){continue;}
      long v=m;
      while(v!=-1)
      {
         long z=v;
         if(bcount[v]==1){isCutEdge2[m]=1; isCutEdge[v]=1; corresponding2[m]=v; corresponding[v]=v;}
         while(nextM[z]!=-1 && bcount[nextM[z]]==bcount[v])
         {
            z=nextM[z];  
         }
         if(z==v){v=nextM[v];continue;}
         for(long u=v;u!=nextM[z];u=nextM[u])
         {
            corresponding[u]=v;
            isCutEdge[u]=1;
         }
         if(bcount[v]!=1)
         {
            vEdgeNext[SP]=vEdgeFirst[v]; vEdgeFirst[v]=SP; vEdgeStack[SP++]=p[z];
            vEdgeNext[SP]=vEdgeFirst[p[z]]; vEdgeFirst[p[z]]=SP; vEdgeStack[SP++]=v;
         }
         v=nextM[z];
      } 
   }

   long* C = (long*)malloc(sizeof(long)*n);
   for(long i=0;i<n;i++){C[i]=-1;}
   long k=0;
   long* Q = (long*)malloc(sizeof(long)*n);
   for(long r=0;r<n;r++)
   {
      if(C[r]!=-1){continue;}
      long first=0; long last=0;
      Q[last++]=r; C[r]=k;
      while(first!=last)
      {
         long x=Q[first++];
         for(long i=firstOut[x];i<firstOut[x+1];i++)
         {
            long y=adj[i];
            if(C[y]!=-1){continue;}
            if((x==p[y]&&isCutEdge[y])||(y==p[x]&&isCutEdge[x])){continue;}
            if((M[x]==x&&y==l[x]&&isCutEdge2[x])||(M[y]==y&&x==l[y]&&isCutEdge2[y])){continue;}
            Q[last++]=y; C[y]=k;
         }
         for(long i=vEdgeFirst[x];i!=-1;i=vEdgeNext[i])
         {
            long y=vEdgeStack[i];
            if(C[y]!=-1){continue;}
            Q[last++]=y; C[y]=k;
         }
      }
      k++;
   }

   long* compStack = (long*)malloc(sizeof(long)*n);
   long* compFirst = (long*)malloc(sizeof(long)*k);
   long* compNext = (long*)malloc(sizeof(long)*n);
   long* compSize = (long*)malloc(sizeof(long)*k);
   for(long i=0;i<k;i++){compFirst[i]=-1; compSize[i]=0;}
   SP=0;
   for(long v=0;v<n;v++)
   {
      long c=C[v];
      compNext[SP]=compFirst[c]; compFirst[c]=SP; compStack[SP++]=v; compSize[c]++;
   }

   long** adjC = (long**)malloc(sizeof(long*)*k);
   long** cname = (long**)malloc(sizeof(long*)*k);
   long* capacity = (long*)malloc(sizeof(long)*k);
   long* size = (long*)malloc(sizeof(long)*k);
   for(long i=0;i<k;i++){capacity[i]=1; size[i]=0;}
   for(long i=0;i<k;i++){adjC[i]=(long*)malloc(sizeof(long)*capacity[i]);}
   for(long i=0;i<k;i++){cname[i]=(long*)malloc(sizeof(long)*capacity[i]);}
   for(long x=0;x<n;x++)
   {
      for(long i=firstOut[x];i<firstOut[x+1];i++)
      {
         long y=adj[i];
         if(C[x]==C[y]){continue;}
         long c=C[x];
         if(size[c]==capacity[c])
         {
            capacity[c]*=2;
            adjC[c]=(long*)realloc(adjC[c],sizeof(long)*capacity[c]);
            cname[c]=(long*)realloc(cname[c],sizeof(long)*capacity[c]);
         }          
         long corr;
         if(dfs[x]>dfs[y]&&y==p[x]){corr=corresponding[x];}
         if(dfs[x]>dfs[y]&&y!=p[x]){corr=corresponding2[x];}
         if(dfs[y]>dfs[x]&&x==p[y]){corr=corresponding[y];}
         if(dfs[y]>dfs[x]&&x!=p[y]){corr=corresponding2[y];}
         adjC[c][size[c]]=C[y]; cname[c][size[c]]=corr; size[c]++;
      }
   }

   long* dfsC = (long*)malloc(sizeof(long)*k);
   long* idfsC = (long*)malloc(sizeof(long)*k);
   for(long i=0;i<k;i++){dfsC[i]=-1;}
   long* temp_node = (long*)malloc(sizeof(long)*k);
   long* temp_indx = (long*)malloc(sizeof(long)*k);
   long Nr=0;
   dfsC[0]=Nr; idfsC[Nr++]=0;
   temp_node[0]=0; temp_indx[0]=0;
   SP=0;
   while(SP!=-1)
   {
      long v=temp_node[SP];
      char descend=0;
      for(long i=temp_indx[SP];i<size[v];i++)
      {
         long u=adjC[v][i];
         if(dfsC[u]==-1)
         {
            dfsC[u]=Nr; idfsC[Nr++]=u;
            temp_node[SP+1]=u; temp_indx[SP+1]=0; temp_indx[SP]=i+1;
            for(long t=0;t<size[u];t++)
            {
               long z=adjC[u][t];
               if(z!=v && cname[u][t]==cname[v][i])
               {
                  long temp;
                  temp=adjC[u][t]; adjC[u][t]=adjC[u][size[u]-1]; adjC[u][size[u]-1]=temp;
                  temp=cname[u][t]; cname[u][t]=cname[u][size[u]-1]; cname[u][size[u]-1]=temp;
                  break; 
               }
            }
            descend=1; break;
         }
      }
      if(descend){SP++;continue;}
      SP--;
   }

   long* leaf_stack = (long*)malloc(sizeof(long)*k);
   char* is_available = (char*)malloc(sizeof(char)*k);
   for(long i=0;i<k;i++){is_available[i]=1;}
   SP=0;
   for(long i=0;i<k;i++)
   {
      long c=idfsC[i];
      if(size[c]==2 && is_available[c])
      {
         leaf_stack[SP++]=c; is_available[c]=0;
      }
   }

   (*edges) = (long*)malloc(sizeof(long)*(firstOut[n]+SP+2*(SP%2)));
   long edgeIndx=0;

   for(long i=0;i<SP/2;i++)
   {
      long c=leaf_stack[i];
      long d=leaf_stack[i+SP/2];
      long i1=rand()%compSize[c];
      long j1=rand()%compSize[d];
      long indx1=compFirst[c];
      while(i1-->0){indx1=compNext[indx1];}
      long x=compStack[indx1];
      long indx2=compFirst[d];
      while(j1-->0){indx2=compNext[indx2];}
      long y=compStack[indx2];
      (*edges)[2*edgeIndx]=x; (*edges)[2*edgeIndx+1]=y; edgeIndx++; //prlongf("added (%d,%d)\n",x,y);
   }
   if(SP%2==1)
   {
      long c=leaf_stack[SP-1];
      long d=leaf_stack[rand()%(SP-1)];
      long i1=rand()%compSize[c];
      long j1=rand()%compSize[d];
      long indx1=compFirst[c];
      while(i1-->0){indx1=compNext[indx1];}
      long x=compStack[indx1];
      long indx2=compFirst[d];
      while(j1-->0){indx2=compNext[indx2];}
      long y=compStack[indx2];
      (*edges)[2*edgeIndx]=x; (*edges)[2*edgeIndx+1]=y; edgeIndx++; //prlongf("added (%d,%d)\n",x,y);
   }

   for(long x=0;x<n;x++)
   {
      for(long i=firstOut[x];i<firstOut[x+1];i++)
      {
         long y=adj[i];
         if(x>y){continue;}
         (*edges)[2*edgeIndx]=x; (*edges)[2*edgeIndx+1]=y; edgeIndx++;
      }
   }

   free(dfs); free(idfs); free(p); free(low); free(l); free(bcount);
   free(low1C); free(low2C); free(M); free(nextM);
   free(vEdgeStack); free(vEdgeFirst); free(vEdgeNext);
   free(isCutEdge); free(isCutEdge2); free(corresponding); free(corresponding2);
   free(C); free(Q);
   free(compStack); free(compFirst); free(compNext); free(compSize);
   for(long i=0;i<k;i++){free(adjC[i]); free(cname[i]);}
   free(adjC); free(cname); free(capacity); free(size);
   free(dfsC); free(idfsC); free(temp_node); free(temp_indx);;
   free(leaf_stack); free(is_available);

   return edgeIndx;
}

void make_2e_connected(long n, long* adj, long* firstOut, long** adj_new, long** firstOut_new)
{
   long* edgeStack = (long*)malloc(sizeof(long)*4*n);
   long* edgeFirst = (long*)malloc(sizeof(long)*n);
   long* edgeNext = (long*)malloc(sizeof(long)*4*n);
   for(long i=0;i<n;i++){edgeFirst[i]=-1;}
   long edgeIndx=0;
   
   char* found = (char*)malloc(sizeof(char)*n);
   for(long i=0;i<n;i++){found[i]=0;}
   long* Q = (long*)malloc(sizeof(long)*n);
   long prevRoot=-1;
   for(long r=0;r<n;r++)
   {
      if(found[r]){continue;}
      long first=0; long last=0;
      Q[last++]=r; found[r]=1;
      while(first!=last)
      {
         long x=Q[first++];
         for(long i=firstOut[x];i<firstOut[x+1];i++)
         {
            long y=adj[i];
            if(!found[y]){Q[last++]=y; found[y]=1;}
         }
      }
      if(prevRoot!=-1)
      {
         edgeNext[edgeIndx]=edgeFirst[prevRoot]; edgeFirst[prevRoot]=edgeIndx; edgeStack[edgeIndx++]=r;
         edgeNext[edgeIndx]=edgeFirst[r]; edgeFirst[r]=edgeIndx; edgeStack[edgeIndx++]=prevRoot;
      }
      prevRoot=r;
   }

   long* dfs = (long*)malloc(sizeof(long)*n);
   long* idfs = (long*)malloc(sizeof(long)*n);
   long* p = (long*)malloc(sizeof(long)*n);
   long* low = (long*)malloc(sizeof(long)*n); 
   for(long i=0;i<n;i++){dfs[i]=-1;}
   long* temp_vertex = (long*)malloc(sizeof(long)*n);
   long* temp_out = (long*)malloc(sizeof(long)*n);
   long* temp_next = (long*)malloc(sizeof(long)*n);
   char* foundP = (char*)malloc(sizeof(char)*n);
   for(long i=0;i<n;i++){foundP[i]=0;low[i]=i;}
   long Nr=0;
   dfs[0]=Nr; idfs[Nr++]=0; p[0]=-1;
   temp_vertex[0]=0; temp_out[0]=firstOut[0]; temp_next[0]=edgeFirst[0];
   long SP=0;
   while(SP!=-1)
   {
      long v=temp_vertex[SP];
      char descend=0; 
      for(long i=temp_out[SP];i<firstOut[v+1];i++)
      {
         long u=adj[i];
         if(dfs[u]==-1)
         {
            dfs[u]=Nr; idfs[Nr++]=u; p[u]=v;
            temp_vertex[SP+1]=u; temp_out[SP+1]=firstOut[u]; temp_out[SP]=i; temp_next[SP+1]=edgeFirst[u];
            descend=1; break; 
         }
         if(u==p[v]&&!foundP[v]){foundP[v]=1;continue;}
         if(dfs[u]<dfs[v])
         {
            if(dfs[u]<dfs[low[v]]){low[v]=u;}
         }
         else if(v==p[u])
         {
            if(dfs[low[u]]<dfs[low[v]]){low[v]=low[u];}
         } 
      }
      if(descend){SP++;continue;}
      for(long i=temp_next[SP];i!=-1;i=edgeNext[i])
      {
         long u=edgeStack[i]; 
         if(dfs[u]==-1)
         {
            dfs[u]=Nr; idfs[Nr++]=u; p[u]=v; low[u]=u;
            temp_vertex[SP+1]=u; temp_out[SP+1]=firstOut[u]; temp_next[SP+1]=edgeFirst[u]; temp_next[SP]=edgeNext[i];
            descend=1; break;
         }
      }
      if(descend){SP++;continue;}
      SP--;
   }

   long* nBridgesToVertex = (long*)malloc(sizeof(long)*n);
   for(long i=0;i<n;i++){nBridgesToVertex[i]=0;}
   for(long v=1;v<n;v++)
   {
      if(low[v]==v)
      {
         nBridgesToVertex[v]++; nBridgesToVertex[p[v]]++;;
      }
   }
   long* nBridgesToComp = (long*)malloc(sizeof(long)*n);
   for(long i=0;i<n;i++){nBridgesToComp[i]=0;}
   long* C = (long*)malloc(sizeof(long)*n);
   C[0]=0; nBridgesToComp[0]=nBridgesToVertex[0]; long k=1;
   for(long i=1;i<n;i++)
   {
      long v=idfs[i];
      if(low[v]!=v){C[v]=C[p[v]];}
      else{C[v]=k++;}
      nBridgesToComp[C[v]]+=nBridgesToVertex[v];
   }

   long* compStack = (long*)malloc(sizeof(long)*n);
   long* compFirst = (long*)malloc(sizeof(long)*k);
   long* compNext = (long*)malloc(sizeof(long)*n);
   long* compSize = (long*)malloc(sizeof(long)*k);
   for(long i=0;i<k;i++){compFirst[i]=-1; compSize[i]=0;}
   SP=0;
   for(long v=0;v<n;v++)
   {
      long c=C[v];
      compNext[SP]=compFirst[c]; compFirst[c]=SP; compStack[SP++]=v; compSize[c]++;
   }

   char* available = (char*)malloc(sizeof(char)*k);
   for(long i=0;i<k;i++){available[i]=1;}
   long* comp_stack = (long*)malloc(sizeof(long)*k);
   SP=0;
   for(long i=0;i<n;i++)
   {
      long c=C[idfs[i]];
      if(nBridgesToComp[c]==1&&available[c])
      {
         available[c]=0; comp_stack[SP++]=c;
      }
   }

   for(long t=0;t<SP/2;t++)
   {
      long i=comp_stack[t];
      long j=comp_stack[t+SP/2];
      long i1 = rand()%compSize[i];
      long j1 = rand()%compSize[j];
      long indx1=compFirst[i];
      for(long t=0;t<i1;t++){indx1=compNext[indx1];}
      long indx2=compFirst[j];
      for(long t=0;t<j1;t++){indx2=compNext[indx2];}
      long x=compStack[indx1]; long y=compStack[indx2];
      edgeNext[edgeIndx]=edgeFirst[x]; edgeFirst[x]=edgeIndx; edgeStack[edgeIndx++]=y;
      edgeNext[edgeIndx]=edgeFirst[y]; edgeFirst[y]=edgeIndx; edgeStack[edgeIndx++]=x;
   }
   if(SP%2==1)
   {
      long i=comp_stack[SP-1];
      long j=comp_stack[rand()%(SP-1)];
      long i1 = rand()%compSize[i];
      long j1 = rand()%compSize[j];
      long indx1=compFirst[i];
      for(long t=0;t<i1;t++){indx1=compNext[indx1];}
      long indx2=compFirst[j];
      for(long t=0;t<j1;t++){indx2=compNext[indx2];}
      long x=compStack[indx1]; long y=compStack[indx2];
      edgeNext[edgeIndx]=edgeFirst[x]; edgeFirst[x]=edgeIndx; edgeStack[edgeIndx++]=y;
      edgeNext[edgeIndx]=edgeFirst[y]; edgeFirst[y]=edgeIndx; edgeStack[edgeIndx++]=x;
   }

   long* edges = (long*)malloc(sizeof(long)*(firstOut[n]+edgeIndx));
   edgeIndx=0;
   for(long x=0;x<n;x++)
   {
      for(long i=firstOut[x];i<firstOut[x+1];i++)
      {
         long y=adj[i];
         if(x>y){continue;}
         edges[2*edgeIndx]=x; edges[2*edgeIndx+1]=y; edgeIndx++;
      }
   }
   for(long x=0;x<n;x++)
   {
      for(long i=edgeFirst[x];i!=-1;i=edgeNext[i])
      {
         long y=edgeStack[i];
         if(x>y){continue;}
         edges[2*edgeIndx]=x; edges[2*edgeIndx+1]=y; edgeIndx++;
      }
   }

   get_adj(n,edgeIndx,edges,adj_new,firstOut_new);

   free(edgeStack); free(edgeFirst); free(edgeNext);
   free(found); free(Q); free(edges);
   free(dfs); free(idfs); free(p); free(low);
   free(temp_vertex); free(temp_out); free(temp_next); free(foundP);
   free(compStack); free(compFirst); free(compNext); free(compSize);
   free(nBridgesToComp); free(nBridgesToVertex); free(C);
   free(available); free(comp_stack);
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
