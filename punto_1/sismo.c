#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265359
#define n_stations 6
#define n_times 20000
#define v 5.0
#define sigma_t 0.1

double rand_norm(double mu, double sigma, double rand_1, double rand_2);
void cambia(int *contador);

int main(){

  double *x = malloc(n_stations*sizeof(double));
  double *y = malloc(n_stations*sizeof(double));
  double *t_obs = malloc(n_stations*sizeof(double));
  double *x_walk = malloc(n_times*sizeof(double));
  double *y_walk = malloc(n_times*sizeof(double));
  double *t_init = malloc(n_stations*sizeof(double));
  double *l_walk = malloc(n_times*sizeof(double));
  double *t_prime = malloc(n_stations*sizeof(double));
  double *rands = malloc((2+5*n_times)*sizeof(double));

  x[0] = 3.0;
  x[1] = 4.0;
  x[2] = 5.0;
  x[3] = 3.0;
  x[4] = 4.0;
  x[5] = 5.0;
  y[0] = 15.0;
  y[1] = 15.0;
  y[2] = 15.0;
  y[3] = 16.0;
  y[4] = 16.0;
  y[5] = 16.0;
  t_obs[0] = 3.12;
  t_obs[1] = 2.98;
  t_obs[2] = 2.84;
  t_obs[3] = 3.26;
  t_obs[4] = 3.12;
  t_obs[5] = 2.98;

  int contador = 0;
  int i;
  double m = 134456.0;
  double n = 8121.0;
  double k = 28411.0;
  double idum = 12.0;
  for(i=0; i<(2+5*n_times); i++){
    idum = fmod(idum*n+k, m);
    rands[i] = idum/m;
  }

  x_walk[0] = 32.0*rands[contador]-12.0;
  cambia(&contador);
  y_walk[0] = 31.0*rands[contador];
  cambia(&contador);

  for(i=0; i<n_stations; i++){
    t_init[i] = sqrt(pow(x[i]-x_walk[0], 2.0)+pow(y[i]-y_walk[0], 2.0))/v;
  }

  double suma_1 = 0.0;
  for(i=0; i<n_stations; i++){
    suma_1 += pow(t_obs[i]-t_init[i], 2.0);
  }
  l_walk[0] = exp(-(1.0/(2*sigma_t))*suma_1);

  double x_prime;
  double y_prime;
  double l_prime;
  double l_init;
  double alpha;
  double betta;
  double suma_2;
  int j;
  for(i=0; i<(n_times-1); i++){
    x_prime = rand_norm(x_walk[i], 0.1, rands[contador], rands[contador+1]);
    cambia(&contador);
    cambia(&contador);
    y_prime = rand_norm(y_walk[i], 0.1, rands[contador], rands[contador+1]);
    cambia(&contador);
    cambia(&contador);
    for(j=0; j<n_stations; j++){
      t_init[j] = sqrt(pow(x[j]-x_walk[i], 2.0)+pow(y[j]-y_walk[i], 2.0))/v;
      t_prime[j] = sqrt(pow(x[j]-x_prime, 2.0)+pow(y[j]-y_prime, 2.0))/v;
    }
    suma_1 = 0.0;
    suma_2 = 0.0;
    for(j=0; j<n_stations; j++){
      suma_1 += pow(t_obs[j]-t_prime[j], 2.0);
      suma_2 += pow(t_obs[j]-t_init[j], 2.0);
    }
    l_prime = exp(-(1.0/(2.0*sigma_t))*suma_1);
    l_init = exp(-(1.0/(2.0*sigma_t))*suma_2);
    alpha = l_prime/l_init;
    if(alpha>=1.0){
      x_walk[i+1] = x_prime;
      y_walk[i+1] = y_prime;
      l_walk[i+1] = l_prime;
    }
    else{
      betta = rands[contador];
      cambia(&contador);
      if(betta<=alpha){
	x_walk[i+1] = x_prime;
	y_walk[i+1] = y_prime;
	l_walk[i+1] = l_prime;
      }
      else{
	x_walk[i+1] = x_walk[i];
	y_walk[i+1] = y_walk[i];
	l_walk[i+1] = l_walk[i];
      }
    }
  }

  for(i=0; i<n_times; i++){
    printf("%f %f\n", x_walk[i], y_walk[i]);
  }

  return 0;

}

double rand_norm(double mu, double sigma, double rand_1, double rand_2){
  double num = sqrt(-2.0*log(rand_1))*cos(2.0*PI*rand_2);
  return mu + sigma*num;
}

void cambia(int *contador){
  *contador += 1;
}
