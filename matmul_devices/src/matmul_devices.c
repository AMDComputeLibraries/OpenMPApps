#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

int numthreads = 256;
int numteams = 56;

#ifdef LARGE
#define N (1000*8)
#else
#define N  (128*8)
#endif
float matA[N][N];
float matB[N][N];
float matC[N][N];
float matE[N][N];

int main(int argc, char **argv) {
  int DevS = omp_get_num_devices();
  char *Env = getenv("OMP_DEV_LIMIT");
  if (Env)
    DevS = atoi(Env);
  fprintf(stderr, "DevS=%d\n",DevS);

// initialize on host
#pragma omp parallel for
  for (int i=0; i < N; i++) {
    for (int j=0; j < N; j++) {
      matA[i][j] = i+j;
      matB[i][j] = i+j;
      matE[i][j] = 0;
    }
  }
// generate expected result on host
#pragma omp parallel for
  for (int i=0; i < N; i++) {
    for (int j=0; j < N; j++) {
      for (int k = 0; k < N; k++) {
        matE[i][j] += matA[i][k] * matB[k][j];
      }
    }
  }

  struct timespec t0,t1,t2;
  fprintf(stderr, "Starting matmul on %d devices\n", DevS);
  clock_gettime(CLOCK_REALTIME, &t0);
  // time the data copy separate from gflops
  for (int dev=0; dev < DevS; dev++) {
    #pragma omp target enter data device(dev) map(to: matA, matB, matC)
  }
  clock_gettime(CLOCK_REALTIME, &t1);
#pragma omp parallel for num_threads(DevS)
  for (int dev=0; dev < DevS; dev++) {
    float tmp;
    const int lb = dev*N/DevS;
    const int ub = (1+dev)*N/DevS;
    // leaving num_threads and teams off, to allow easier env-var control
    #pragma omp target teams device(dev)//num_teams(numteams) //thread_limit(numthreads)
    #pragma omp distribute parallel for private(tmp) collapse(2)
    for (int i = lb; i < ub; i++) {
      for (int j = 0; j < N; j++) {
        tmp = 0.0;
        for (int k = 0; k < N; k++) {
          tmp += matA[i][k] * matB[k][j];
        }
        matC[i][j] = tmp;
      }
    }
   fprintf(stderr, "Done with %d dev\n",dev);
  }
  for (int dev=0; dev <DevS; dev++) {
  #pragma omp target exit data device(dev) map(from: matC[dev*N/DevS:N/DevS])
  }
  clock_gettime(CLOCK_REALTIME, &t2);
  float ops = (float)N *N *N *2.0;
  double m = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec)/1e9;
  fprintf(stderr, "Time %g GBytes/sec %g for copy\n", m, N*N*3.0*4/m/1e9);
  double t = (t2.tv_sec - t1.tv_sec) + (t2.tv_nsec - t1.tv_nsec)/1e9;
  fprintf(stderr, "Time %g GFlops/sec %g for compute\n", t, ops/t/1e9);
  double c = (t2.tv_sec - t0.tv_sec) + (t2.tv_nsec - t0.tv_nsec)/1e9;
  fprintf(stderr, "Time %g GFlops/sec %g for copy + compute\n", c, ops/c/1e9);

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      if (matE[i][j] != matC[i][j]) {
        fprintf(stderr, "Failed %d %d\n",i,j);
        return 1;
      }
    }
  }
  fprintf(stderr, "Passed\n");
  return 0;
}
