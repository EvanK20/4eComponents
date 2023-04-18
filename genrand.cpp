//Runs with ./genrand n m seed <output_file>

#include <stdio.h>
#include <stdlib.h>

long tolong(char*);

int main(int n_args, char** args)
{
   long n = tolong(args[1]);
   long m = tolong(args[2]);
   srand(tolong(args[3]));
   long* edges = (long*)malloc(sizeof(long)*2*m);
   for(long i=0;i<m;i++)
   {
      long x=rand()%n; long y=(x+1+rand()%(n-1))%n;
      edges[2*i]=x; edges[2*i+1]=y;
   }
   FILE* fp = fopen(args[4],"w");
   fprintf(fp,"%ld %ld\n",n,m);
   for(long i=0;i<m;i++){fprintf(fp,"%ld %ld\n",edges[2*i],edges[2*i+1]);}
   fclose(fp);
   free(edges);
   return 0;
}

long tolong(char* x)
{
   long l=0;
   while(x[l]!=0){l++;}
   long pow10=1;
   long ret=0;
   l--;
   while(l!=-1){ret+=pow10*(x[l]-'0');l--;pow10*=10;}
   return ret;
}
