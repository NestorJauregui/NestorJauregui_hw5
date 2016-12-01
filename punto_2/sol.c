#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265359
#define G 19.8254558E-30
#define n_data 8
#define n_times 20000

double rand_norm(double mu, double sigma, double rand_1, double rand_2);
void cambia(int *contador);

int main(){

  double *radios = malloc(n_data*sizeof(double));
  double *velos = malloc(n_data*sizeof(double));
  double *y_obs = malloc(n_data*sizeof(double));
  double *x_obs = malloc(n_data*sizeof(double));
  double *m_walk = malloc(n_times*sizeof(double));
  double *b_walk = malloc(n_times*sizeof(double));
  double *alpha_walk = malloc(n_times*sizeof(double));
  double *log10Ms_walk = malloc(n_times*sizeof(double));
  double *l_walk = malloc(n_times*sizeof(double));
  double *y_init = malloc(n_data*sizeof(double));
  double *y_prime = malloc(n_data*sizeof(double));
  double *rands = malloc((2+5*n_times)*sizeof(double));

  radios[0] = sqrt(pow(0.324190175, 2.0)+pow(0.090955208, 2.0)+pow(0.02292051, 2.0));
  radios[1] = sqrt(pow(0.701534590, 2.0)+pow(0.168809218, 2.0)+pow(0.037947785, 2.0));
  radios[2] = sqrt(pow(0.982564148, 2.0)+pow(0.191145980, 2.0)+pow(0.000014724, 2.0));
  radios[3] = sqrt(pow(1.104185888, 2.0)+pow(0.826097003, 2.0)+pow(0.044595990, 2.0));
  radios[4] = sqrt(pow(3.266443877, 2.0)+pow(3.888055863, 2.0)+pow(0.057015321, 2.0));
  radios[5] = sqrt(pow(9.218802228, 2.0)+pow(1.788299816, 2.0)+pow(0.335737817, 2.0));
  radios[6] = sqrt(pow(19.930781147, 2.0)+pow(2.555241579, 2.0)+pow(0.267710968, 2.0));
  radios[7] = sqrt(pow(24.323085642, 2.0)+pow(17.606227355, 2.0)+pow(0.177974999, 2.0));
  velos[0] = sqrt(pow(4.627851589, 2.0)+pow(10.390063716, 2.0)+pow(1.273504997, 2.0));
  velos[1] = sqrt(pow(1.725066954, 2.0)+pow(7.205747212, 2.0)+pow(0.198268558, 2.0));
  velos[2] = sqrt(pow(1.126784520, 2.0)+pow(6.187988860, 2.0)+pow(0.000330572, 2.0));
  velos[3] = sqrt(pow(3.260215854, 2.0)+pow(4.524583075, 2.0)+pow(0.014760239, 2.0));
  velos[4] = sqrt(pow(2.076140727, 2.0)+pow(1.904040630, 2.0)+pow(0.054374153, 2.0));
  velos[5] = sqrt(pow(0.496457364, 2.0)+pow(2.005021061, 2.0)+pow(0.054667082, 2.0));
  velos[6] = sqrt(pow(0.172224285, 2.0)+pow(1.357933443, 2.0)+pow(0.002836325, 2.0));
  velos[7] = sqrt(pow(0.664855006, 2.0)+pow(0.935497207, 2.0)+pow(0.034716967, 2.0));

  int contador = 0;
  int i;
  double m = 134456.0;
  double n = 8121.0;
  double k = 28411.0;
  double idum = 231.0;
  for(i=0; i<(2+5*n_times); i++){
    idum = fmod((idum*n+k),m);
    rands[i] = idum/m;
  }

  m_walk[0] = rands[contador];
  cambia(&contador);
  b_walk[0] = rands[contador];
  cambia(&contador);

  for(i=0; i<n_data; i++){
    y_obs[i] = log(velos[i]);
    x_obs[i] = log(radios[i]);
    y_init[i] = x_obs[i]*m_walk[0]+b_walk[0];
  }

  double suma_1 = 0.0;
  double suma_2 = 0.0;
  for(i=0; i<n_data; i++){
    suma_1 += pow(y_obs[i]-y_init[i], 2.0);
  }
  l_walk[0] = exp(-0.5*suma_1);

  double m_prime;
  double b_prime;
  double l_prime;
  double l_init;
  double alpha;
  double betta;
  int j;
  for(i=0; i<(n_times-1); i++){
    m_prime = rand_norm(m_walk[i], 0.1, rands[contador], rands[contador+1]);
    cambia(&contador);
    cambia(&contador);
    b_prime = rand_norm(b_walk[i], 0.1, rands[contador], rands[contador+1]);
    cambia(&contador);
    cambia(&contador);
    for(j=0; j<n_data; j++){
      y_init[j] = (x_obs[j]*m_walk[i])+b_walk[i];
      y_prime[j] = (x_obs[j]*m_prime)+b_prime;
    }
    suma_1 = 0.0;
    suma_2 = 0.0;
    for(j=0; j<n_data; j++){
      suma_1 += pow(y_obs[j]-y_prime[j], 2.0);
      suma_2 += pow(y_obs[j]-y_init[j], 2.0);
    }
    l_prime = exp(-0.5*suma_1);
    l_init = exp(-0.5*suma_2);
    alpha = l_prime/l_init;
    if(alpha>=1.0){
      m_walk[i+1] = m_prime;
      b_walk[i+1] = b_prime;
      l_walk[i+1] = l_prime;
    }
    else{
      betta = rands[contador];
      cambia(&contador);
      if(betta<=alpha){
	m_walk[i+1] = m_prime;
	b_walk[i+1] = b_prime;
	l_walk[i+1] = l_prime;
      }
      else{
	m_walk[i+1] = m_walk[i];
	b_walk[i+1] = b_walk[i];
	l_walk[i+1] = l_init;
      }
    }
  }

  for(i=0; i<n_times; i++){
    alpha_walk[i] = 1.0 - 2.0*m_walk[i];
    log10Ms_walk[i] = (2.0*b_walk[i]-log(G))/log(10.0);
    printf("%f %f\n", alpha_walk[i], log10Ms_walk[i]);
  }

  return 0;

}

double rand_norm(double mu, double sigma, double rand_1, double rand_2){
  double num = sqrt(-2.0*log(rand_1))*sin(2.0*PI*rand_2);
  return mu + sigma*num;
}

void cambia(int *contador){
  *contador += 1;
}
