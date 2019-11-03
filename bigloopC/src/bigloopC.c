#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

const int numthreads = 256;
const int numteams = 56;
#ifndef NUMTHREADS
#define NUMTHREADS numthreads
#endif
#ifndef NUMTEAMS
#define NUMTEAMS numteams
#endif

#define N (10000)
float matA[N][N];
float matB[N][N];
float matC[N][N];
float matE[N][N];

int main(int argc, char **argv) {
// initialize on host
#pragma omp parallel for
  for (int i=0; i < N; i++) {
    for (int j=0; j < N; j++) {
      matA[i][j] = i+j;
      matB[i][j] = 10;
      matE[i][j] = 0;
    }
  }
// generate expected result on host
#pragma omp parallel for
  for (int i=0; i < 20; i++) {
    for (int j=0; j < N; j++) {
      for (int k = 0; k < N; k++) {
        matE[j][k] = matA[j][k] * matB[j][k];
      }
    }
  }

  struct timespec t0,t1,t2;
  fprintf(stderr, "Starting matmul\n");
  // time the data copy separate from gflops
  clock_gettime(CLOCK_REALTIME, &t0);
  #pragma omp target data map(to: matA, matB) map(from: matC)
  {
    clock_gettime(CLOCK_REALTIME, &t1);
    float tmp;
    // leaving num_threads and teams off, to allow easier env-var control
    #pragma omp target teams num_teams(NUMTEAMS) thread_limit(NUMTHREADS)
    #pragma omp distribute parallel for private(tmp) collapse(2)
    for (int i = 0; i < 20; i++) {
      for (int j = 0; j < N; j++) {
        for (int k = 0; k < N; k++) {
          matC[j][k] = matA[j][k] * matB[j][k];
        }
      }
    }
    clock_gettime(CLOCK_REALTIME, &t2);
  }
  double m = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec)/1e9;
  fprintf(stderr, "Time %g for copy\n", m);
  double t = (t2.tv_sec - t1.tv_sec) + (t2.tv_nsec - t1.tv_nsec)/1e9;
  fprintf(stderr, "Time %g for compute\n", t);
  double c = (t2.tv_sec - t0.tv_sec) + (t2.tv_nsec - t0.tv_nsec)/1e9;
  fprintf(stderr, "Time %g for copy + compute\n", c);

  for (int j = 0; j < N; j++) {
    for (int k = 0; k < N; k++) {
      if (matE[j][k] != matC[j][k]) {
        fprintf(stderr, "Failed %d %d\n",j,k);
        return 1;
      }
    }
  }
  fprintf(stderr, "Passed\n");
  return 0;
}
