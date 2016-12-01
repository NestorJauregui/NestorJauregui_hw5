#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX 80
#define n_data 96
#define PI 3.14159265359
#define n_times 20000
#define t_min 0.006
#define t_max 0.799
#define h (t_max-t_min)/n_data

double deriv_x(double x, double y, double alpha, double betta);
double deriv_y(double x, double y, double gamma, double delta);
double rand_norm(double mu, double sigma, double rand_1, double rand_2);
void cambia(int *contador);

int main(){

  FILE* file;
  char* filename = "lotka_volterra_obs.dat";
  char lectura[MAX];
  file = fopen(filename, "r");
  double *x = malloc(n_data*sizeof(double));
  double *y = malloc(n_data*sizeof(double));
  double *rands = malloc((4+10*n_times)*sizeof(double));
  double *alpha_walk = malloc(n_times*sizeof(double));
  double *betta_walk = malloc(n_times*sizeof(double));
  double *gamma_walk = malloc(n_times*sizeof(double));
  double *delta_walk = malloc(n_times*sizeof(double));
  double *x_init = malloc(n_data*sizeof(double));
  double *y_init = malloc(n_data*sizeof(double));
  double *l_walk_x = malloc(n_times*sizeof(double));
  double *l_walk_y = malloc(n_times*sizeof(double));
  double *x_prime = malloc(n_data*sizeof(double));
  double *y_prime = malloc(n_data*sizeof(double));

  int i;
  int ident = 0;
  int ind = 0;
  for(i=0; i<(5+n_data*3); i++){
    fscanf(file, "%s", lectura);
    if(i>=5){
      if(ident==1){
	x[ind] = atof(lectura);
      }
      else if(ident==2){
	y[ind] = atof(lectura);
      }
      ident += 1;
      if(ident==3){
	ident = 0;
	ind += 1;
      }
    }
  }

  fclose(file);

  int contador = 0;
  double m = 134456.0;
  double n = 8121.0;
  double k = 28411.0;
  double idum = 5.0;
  for(i=0; i<(4+10*n_times); i++){
    idum = fmod(idum*n+k, m);
    rands[i] = idum/m;
  }

  alpha_walk[0] = rands[contador];
  cambia(&contador);
  betta_walk[0] = rands[contador];
  cambia(&contador);
  gamma_walk[0] = rands[contador];
  cambia(&contador);
  delta_walk[0] = rands[contador];
  cambia(&contador);

  x_init[0] = x[0];
  y_init[0] = y[0];
  double k1_x;
  double k2_x;
  double k3_x;
  double k4_x;
  double k1_y;
  double k2_y;
  double k3_y;
  double k4_y;
  double average_k_x;
  double average_k_y;
  for(i=1; i<n_data; i++){
    k1_x = h*deriv_x(x[i-1], y[i-1], alpha_walk[0], betta_walk[0]);
    k1_y = h*deriv_y(x[i-1], y[i-1], gamma_walk[0], delta_walk[0]);
    k2_x = h*deriv_x(x[i-1]+0.5*k1_x, y[i-1]+0.5*k1_y, alpha_walk[0], betta_walk[0]);
    k2_y = h*deriv_y(x[i-1]+0.5*k1_x, y[i-1]+0.5*k1_y, gamma_walk[0], delta_walk[0]);
    k3_x = h*deriv_x(x[i-1]+0.5*k2_x, y[i-1]+0.5*k2_y, alpha_walk[0], betta_walk[0]);
    k3_y = h*deriv_y(x[i-1]+0.5*k2_x, y[i-1]+0.5*k2_y, gamma_walk[0], delta_walk[0]);
    k4_x = h*deriv_x(x[i-1]+k3_x, y[i-1]+k3_y, alpha_walk[0], betta_walk[0]);
    k4_y = h*deriv_y(x[i-1]+k3_x, y[i-1]+k3_y, gamma_walk[0], delta_walk[0]);
    average_k_x = (1.0/6.0)*(k1_x+2.0*k2_x+2.0*k3_x+k4_x);
    average_k_y = (1.0/6.0)*(k1_y+2.0*k2_y+2.0*k3_y+k4_y);
    x_init[i] = x[i-1] + average_k_x;
    y_init[i] = y[i-1] + average_k_y;
  }

  double suma_1 = 0.0;
  double suma_2 = 0.0;
  for(i=0; i<n_data; i++){
    suma_1 += pow(x[i]-x_init[i], 2.0);
    suma_2 += pow(y[i]-y_init[i], 2.0);
  }
  l_walk_x[0] = exp(-0.5*suma_1);
  l_walk_y[0] = exp(-0.5*suma_2);

  double alpha_prime;
  double betta_prime;
  double gamma_prime;
  double delta_prime;
  double l_prime_x;
  double l_prime_y;
  double l_init_x;
  double l_init_y;
  double alpha_x;
  double alpha_y;
  double betta_x;
  double betta_y;
  double suma_3;
  double suma_4;
  int j;
  for(i=0; i<(n_times-1); i++){
    alpha_prime = rand_norm(alpha_walk[i], 0.1, rands[contador], rands[contador+1]);
    cambia(&contador);
    cambia(&contador);
    betta_prime = rand_norm(betta_walk[i], 0.1, rands[contador], rands[contador+1]);
    cambia(&contador);
    cambia(&contador);
    gamma_prime = rand_norm(gamma_walk[i], 0.1, rands[contador], rands[contador+1]);
    cambia(&contador);
    cambia(&contador);
    delta_prime = rand_norm(delta_walk[i], 0.1, rands[contador], rands[contador+1]);
    cambia(&contador);
    cambia(&contador);
    for(j=1; j<n_data; j++){
      k1_x = h*deriv_x(x[j-1], y[j-1], alpha_walk[i], betta_walk[i]);
      k1_y = h*deriv_y(x[j-1], y[j-1], gamma_walk[i], delta_walk[i]);
      k2_x = h*deriv_x(x[j-1]+0.5*k1_x, y[j-1]+0.5*k1_y, alpha_walk[i], betta_walk[i]);
      k2_y = h*deriv_y(x[j-1]+0.5*k1_x, y[j-1]+0.5*k1_y, gamma_walk[i], delta_walk[i]);
      k3_x = h*deriv_x(x[j-1]+0.5*k2_x, y[j-1]+0.5*k2_y, alpha_walk[i], betta_walk[i]);
      k3_y = h*deriv_y(x[j-1]+0.5*k2_x, y[j-1]+0.5*k2_y, gamma_walk[i], delta_walk[i]);
      k4_x = h*deriv_x(x[j-1]+k3_x, y[j-1]+k3_y, alpha_walk[i], betta_walk[i]);
      k4_y = h*deriv_y(x[j-1]+k3_x, y[j-1]+k3_y, gamma_walk[i], delta_walk[i]);
      average_k_x = (1.0/6.0)*(k1_x+2.0*k2_x+2.0*k3_x+k4_x);
      average_k_y = (1.0/6.0)*(k1_y+2.0*k2_y+2.0*k3_y+k4_y);
      x_init[j] = x[j-1] + average_k_x;
      y_init[j] = y[j-1] + average_k_y;
      k1_x = h*deriv_x(x[j-1], y[j-1], alpha_prime, betta_prime);
      k1_y = h*deriv_y(x[j-1], y[j-1], gamma_prime, delta_prime);
      k2_x = h*deriv_x(x[j-1]+0.5*k1_x, y[j-1]+0.5*k1_y, alpha_prime, betta_prime);
      k2_y = h*deriv_y(x[j-1]+0.5*k1_x, y[j-1]+0.5*k1_y, gamma_prime, delta_prime);
      k3_x = h*deriv_x(x[j-1]+0.5*k2_x, y[j-1]+0.5*k2_y, alpha_prime, betta_prime);
      k3_y = h*deriv_y(x[j-1]+0.5*k2_x, y[j-1]+0.5*k2_y, gamma_prime, delta_prime);
      k4_x = h*deriv_x(x[j-1]+k3_x, y[j-1]+k3_y, alpha_prime, betta_prime);
      k4_y = h*deriv_y(x[j-1]+k3_x, y[j-1]+k3_y, gamma_prime, delta_prime);
      average_k_x = (1.0/6.0)*(k1_x+2.0*k2_x+2.0*k3_x+k4_x);
      average_k_y = (1.0/6.0)*(k1_y+2.0*k2_y+2.0*k3_y+k4_y);
      x_prime[j] = x[j-1] + average_k_x;
      y_prime[j] = y[j-1] + average_k_y;
    }
    suma_1 = 0.0;
    suma_2 = 0.0;
    suma_3 = 0.0;
    suma_4 = 0.0;
    for(j=0; j<n_data; j++){
      suma_1 = pow(x[j]-x_prime[j], 2.0);
      suma_2 = pow(x[j]-x_init[j], 2.0);
      suma_3 = pow(y[j]-y_prime[j], 2.0);
      suma_4 = pow(y[j]-y_init[j], 2.0);
    }
    l_prime_x = exp(-0.5*suma_1);
    l_init_x = exp(-0.5*suma_2);
    l_prime_y = exp(-0.5*suma_3);
    l_init_y = exp(-0.5*suma_4);
    alpha_x = l_prime_x/l_init_x;
    alpha_y = l_prime_y/l_init_y;
    if(alpha_x>=1.0){
      alpha_walk[i+1] = alpha_prime;
      betta_walk[i+1] = betta_prime;
      l_walk_x[i+1] = l_prime_x;
    }
    else{
      betta_x = rands[contador];
      cambia(&contador);
      if(betta_x<=alpha_x){
	alpha_walk[i+1] = alpha_prime;
	betta_walk[i+1] = betta_prime;
	l_walk_x[i+1] = l_prime_x;
      }
      else{
	alpha_walk[i+1] = alpha_walk[i];
	betta_walk[i+1] = betta_walk[i];
	l_walk_x[i+1] = l_walk_x[i];
      }
    }
    if(alpha_y>=1.0){
      gamma_walk[i+1] = gamma_prime;
      delta_walk[i+1] = delta_prime;
      l_walk_y[i+1] = l_prime_y;
    }
    else{
      betta_y = rands[contador];
      cambia(&contador);
      if(betta_y<=alpha_y){
	gamma_walk[i+1] = gamma_prime;
	delta_walk[i+1] = delta_prime;
	l_walk_y[i+1] = l_prime_y;
      }
      else{
	gamma_walk[i+1] = gamma_walk[i];
	delta_walk[i+1] = delta_walk[i];
	l_walk_y[i+1] = l_walk_y[i];
      }
    }
  }

  for(i=0; i<n_times; i++){
    printf("%f %f %f %f\n", alpha_walk[i], betta_walk[i], gamma_walk[i], delta_walk[i]);
  }

  return 0;

}

double deriv_x(double x, double y, double alpha, double betta){
  return x*(alpha-betta*y);
}

double deriv_y(double x, double y, double gamma, double delta){
  return -y*(gamma-delta*x);
}

double rand_norm(double mu, double sigma, double rand_1, double rand_2){
  double num = sqrt(-2.0*log(rand_2))*cos(2.0*PI*rand_1);
  return mu + sigma*num;
}

void cambia(int *contador){
  *contador += 1;
}
