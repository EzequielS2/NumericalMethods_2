

#include <stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>

double **LeMatriz(char *nome, int *m, int *n){

  FILE *fp = fopen(nome, "r");
  int i, j, m1, n1;
  double **Matriz;

  fscanf(fp, "%d %d", &(*m), &(*n));
  m1=*m, n1=*n;
  Matriz=(double **)malloc(m1*sizeof(double *));
  for(i=0; i<m1; i++){
    Matriz[i]=(double *)malloc(n1*sizeof(double));
  }
  for(i=0; i<m1; i++){
    for(j=0; j<n1; j++){
      fscanf(fp, "%lf", &Matriz[i][j]);
    }
  }
  return Matriz;
}

void ImprimeVetor(double *b, int m)
{
  int i;
  for(i=0; i<m; i++) printf("[%d]= %g\n", i, b[i]);
  puts("");
}

void ImprimaMatriz(double **MT, int m, int n){
    int i, j;
   
    for (i = 0; i < m; i++){
     for (j = 0; j < n; j++) {
         printf("%g\t", MT[i][j]);  
      }
     puts("");
    }   
    puts("");
}

double NormaVetor(double *M, int m, int p){
  int i;
  double aux;
  if(p==0){
    aux=fabs(M[0]);
    for(i=0; i<m; i++){
      if(fabs(M[i])>aux){
      aux=fabs(M[i]);
      }
    }
    return aux;
  }
  else{
    for(i=0; i<m;i++){
      aux+=fabs(M[i]);
    }
    return pow(aux, 1.0/p);
  }
}

double NormaMatriz(double **M, int m, int n, int p){
  int i, j;
  double aux1=0, aux2;
  
  if(p==0){
    for(i=0; i<m; i++){
      aux1+=fabs(M[i][0]);
    }
    for(i=0; i<n;i++){
      aux2=0;
      for(j=0; j<m; j++){
        aux2+=fabs(M[j][i]);
      }
      if(aux2>aux1) aux1=aux2;
    }
    return aux1;
  }
  else if(p==1){
    for(i=0; i<n; i++){
      aux1+=fabs(M[0][i]);
    }
    for(i=0; i< m; i++){
      aux2=0;
      for(j=0; j<n; j++){
        aux2+=fabs(M[i][j]);
      }
      if(aux2>aux1) aux1=aux2; 
    }
  }
  else if(p==2){
    for(i=0; i<m; i++){
      for(j=0; j<n; j++){
        aux1+=pow(fabs(M[i][j]), 2);
      }
    }
    return sqrt(aux1);
  }
}




double *LeVetor(char *nome, int *m){

  FILE *fp = fopen(nome, "r");
  int i;
  double *V;

  fscanf(fp, "%d", &(*m));
  V=calloc(*m , sizeof(double *));
  for(i=0; i<*m; i++){
      fscanf(fp, "%lf", &V[i]);    
  }
  return V;
}


double *Erro(double *V1, double *V2, int m){

  int i, j;
  double *Er;
  Er=calloc(m , sizeof(double *));
  for(i=0; i<m; i++){
    Er[i]=V1[i]-V2[i];
  }

  return Er;
}

double Jacobi(double **Matriz, int m, int n, double *x0, int p){

  int i, j;
  double soma, *x, *b, *err;
  x=calloc(m , sizeof(double *));
  b=calloc(m , sizeof(double *));
  err=calloc(m , sizeof(double *));
  for(i=0; i<m; i++) (b[i])=Matriz[i][n-1];
  for(i=0; i<m; i++) x[i]=x0[i];
  for(i=0; i<m; i++){
        soma = 0;
        for(j=0; j<n-1; j++){
            if(j!=i){
                soma +=Matriz[i][j]*x0[j]/Matriz[i][i];
            }
        }
        x0[i] = (b[i]/Matriz[i][i]) - soma;
  }
  err=Erro(x0, x, m);
  return NormaVetor(err, m,  p);   
}



double Gauss(double **Matriz, int m, int n, double *x0, int p){

int i, j, it=1, k; 

double soma, norm, *x, *b, *err;

x=calloc(m , sizeof(double *));
b=calloc(m , sizeof(double *));
err=calloc(m , sizeof(double *));
for(i=0; i<m; i++) (b[i])=Matriz[i][n-1];
for(i=0; i<m; i++) x[i]=x0[i];
for(i=0; i<m; i++){
    soma=0; 
    for(j = 0; j < n-1; j++) {
      if(j!=i){
        soma+=Matriz[i][j] *x0[j];
      }
	  }
    x0[i] = (b[i]-soma)/Matriz[i][i];
}
err=Erro(x0, x, m);
return NormaVetor(err, m,  p);
}


double Relaxacao(double **Matriz, int m, int n, double *x0, double omega, int p){

int i, j; 
double soma, norm=0, *x, *b, *err;
x=calloc(m , sizeof(double *));
b=calloc(m , sizeof(double *));
err=calloc(m , sizeof(double *));
for(i=0; i<m; i++) (b[i])=Matriz[i][n-1];
for(i=0; i<m; i++) x[i]=x0[i];
for (i = 0; i < m; i++) {
      soma = 0;
      for (j = 0; j < n-1; j++) {
                if (j != i)
                    soma += Matriz[i][j]*x0[j];
      }
      x0[i] = (1-omega)*x0[i]+(omega/Matriz[i][i])*(b[i]-soma); 
}
err=Erro(x0, x, m);
return NormaVetor(err, m,  p);
}


double *MatrizxVetor(double **M, double *v, int m){
  int i, j;
  double *x;
  x=calloc(m , sizeof(double *));
  for(i=0; i<m; i++){
    for(j=0; j<m;j++){
        x[i]+=M[i][j]*v[j];
    }
  }
return x;
}


double ProInt(double *v1, double *v2, int m){
  int i;
  double prod=0;

  for(i=0;i<m;i++){
    prod+=v1[i]*v2[i];
  }
  return prod;
}


void GradienteConjugado(double **M, double *x0, int m, double tol, double p){

          int i,j, k, it=0;
          double *d1, *d2, *x, *Md, *err, *b, *r;
          double prod1, dMd, dx, alp, prod2;
          d1=calloc(m , sizeof(double *));
          d2=calloc(m , sizeof(double *));
          r=calloc(m , sizeof(double *));
          x=calloc(m , sizeof(double *));
          b=calloc(m , sizeof(double *));
          Md=calloc(m , sizeof(double *));
          for(i=0; i<m; i++) (b[i])=M[i][m];
          d2=MatrizxVetor(M, x0,  m); 
          for(i = 0; i < m; i++)
          {     
               d2[i]=d2[i]-b[i];
               d1[i] = -d2[i];
          }
          prod1=ProInt(d2, d2, m);
          do {  
               it++;  
               memcpy(x, x0, m*sizeof(double *));
               memcpy(r, d2, m*sizeof(double *));
               Md=MatrizxVetor(M, d1,m); 
               prod2=prod1;
               alp=(ProInt(d2, d2, m)/ProInt(d1, MatrizxVetor(M, d1,m), m));
               for(i = 0; i < m; i++)
               {    x0[i]=x0[i]+alp*d1[i];
                    d2[i]=d2[i] + alp*Md[i];
               }
               prod1=ProInt(d2,d2,m);
               for(i = 0; i < m; i++){
                    d1[i]=-d2[i]+(prod1/prod2)*d1[i];
               }
               err=Erro(x, x0, m);
               dx=NormaVetor(err, m,  p);
               printf("\n\nIt=%d Dx=%8.4g r0*r1= %8.4g:\n\n", it,dx, ProInt(r, d2,m));
               // ProInt(r, d2,m))=r0*r1
               for( i=0; i<m; i++) printf("x%d=%11.6g ", i, x0[i]);
               puts("");
          }while(dx>tol);
}


double Gradient(double **M, double *x0, int m, double tol, double p, double *res){

  int i,j, k;
  double *r, *re, *b, *err, *x, lbd;
  r=calloc(m , sizeof(double *));
  re=calloc(m , sizeof(double *));
  b=calloc(m , sizeof(double *));
  err=calloc(m , sizeof(double *));
  x=calloc(m , sizeof(double *));
  for(i=0; i<m; i++) (b[i])=M[i][m];
  memcpy(x, x0, m*sizeof(double *));
  r=MatrizxVetor(M, x0,  m); 
  for(i = 0; i < m; i++)
  {     
    r[i]=-r[i]+b[i];
  }
  memcpy(re, r, m*sizeof(double *));
  lbd=ProInt(r, r,  m)/ProInt(r, MatrizxVetor(M, r, m), m);
  for(i=0; i<m; i++){
    x0[i]=x0[i]+lbd*r[i];
  }
  r=MatrizxVetor(M, x0,  m); 
  for(i = 0; i < m; i++)
  {     
    r[i]=-r[i]+b[i];
  }
  *res=ProInt(r, re,  m);
  err=Erro(x, x0, m);
  return NormaVetor(err, m,  p);
} 



int main(int argc, char **argv){

 double **M, *v, *x, dx, res,tolerance=pow(10,-8);
 int m, n, l, i, it=0, p=0;
 M = LeMatriz(argv[1], &m, &n);
 v = LeVetor(argv[2], &l);
 puts("\nLetras a) e b):\n");
 puts("\nChutes iniciais\n");
 for( i=0; i<m; i++) printf("%11.6g ", v[i]);
 puts("\n");
 puts("\n1) Passos do método -Gradientes Conjuagados- (mostrando que os resíuos são ortogonais na terceira coluna: )\n");
 GradienteConjugado(M, v, m,tolerance, p);
 puts("\n2) Passos do método -Máxima descida- (mostrando que os resíuos são ortogonais na terceira coluna): \n");
 v = LeVetor(argv[2], &l);
 it=0;
 do{
  it++;
  dx = Gradient(M, v, m, tolerance, p, &res);
  printf("\n\nIt=%d Dx=%8.4g r0*r1=%8.4g:\n\n", it,dx, res);
  // res=(r0)*(r1) é o produto interno dos resíduos em cada processo 
  for( i=0; i<m; i++) printf("x%d=%11.6g ", i, v[i]);
  puts("");
} while (dx > tolerance);
puts("\n3) Passos do método -Gauss-Siedel-: \n");
it=0;
v = LeVetor(argv[2], &l);
 do{
  it++;
  dx = Gauss(M, m,n,v,p);
  printf("\n\nIt=%d Dx=%8.4g\n\n ", it,dx);
  for( i=0; i<m; i++) printf("x%d=%11.6g ", i, v[i]);
  puts("");
} while (dx > tolerance); 
}

