#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

#if defined(__OFFLOAD_ARCH_gfx940__) || defined(__OFFLOAD_ARCH_gfx941__) || defined(__OFFLOAD_ARCH_gfx942__) || defined(__OFFLOAD_ARCH_gfx90a__)
#define IS_USM 1
#else
#define IS_USM 0
#endif
#if IS_USM >=1
#pragma omp requires unified_shared_memory
#endif

int numthreads = 256;
int numteams = 56;

#ifdef LARGE
#define N (10*7*5*4*3*2)  // 8400
#else
#define N  (5*7*5*3*8*2)    // 1050
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
  if (DevS > 8) DevS = 8;
  fprintf(stderr, "DevS=%d\n",DevS);

  if (N % DevS != 0) {
    fprintf(stderr, "Matrix Size %d not multiple of DevS %d\n", N , DevS);
    return 1;
  }
#if USM
  fprintf(stderr, "building for USM\n");
#else
  fprintf(stderr, "building for none-USM\n");
#endif
fprintf(stderr,"Init on host\n");
// initialize on host
#pragma omp parallel for
  for (int i=0; i < N; i++) {
    for (int j=0; j < N; j++) {
      matA[i][j] = i+j;
      matB[i][j] = i+j;
      matE[i][j] = 0;
    }
  }
fprintf(stderr,"Generate expected\n");
// generate expected result on host
#pragma omp parallel for
  for (int i=0; i < N; i++) {
    for (int j=0; j < N; j++) {
      for (int k = 0; k < N; k++) {
        matE[i][j] += matA[i][k] * matB[k][j];
      }
    }
  }

  struct timespec t0,t1,t2,t3;
  fprintf(stderr, "device init\n");
  clock_gettime(CLOCK_REALTIME, &t3);
  int alpha=0;
  #pragma omp target map(tofrom:alpha)
  {
  alpha =1;
  }
  fprintf(stderr, "Alpha=%d\n",alpha);
  clock_gettime(CLOCK_REALTIME, &t0);

  float ops = (float)N *N *N *2.0;
  double d = (t0.tv_sec - t3.tv_sec) + (t0.tv_nsec - t3.tv_nsec)/1e9;
  fprintf(stderr, "Time %g GBytes/sec %g for devinit\n", d, N*N*3.0*4/d/1e9);

  fprintf(stderr, "Starting matmul %dx%d on %d devices\n", N, N, DevS);
  // time the data copy separate from gflops
//#if IS_USM
  for (int dev=0; dev < DevS; dev++) {
    #pragma omp target enter data device(dev) map(to: matA, matB, matC)
  }
//#endif
  clock_gettime(CLOCK_REALTIME, &t1);
  double m = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec)/1e9;
  fprintf(stderr, "Time %g GBytes/sec %g for copy\n", m, N*N*3.0*4/m/1e9);
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
    fprintf(stderr, "Done with device # %d\n", dev);
  }
//#if IS_USM
  for (int dev=0; dev <DevS; dev++) {
  #pragma omp target exit data device(dev) map(from: matC[dev*N/DevS:N/DevS])
  }
//#endif
  clock_gettime(CLOCK_REALTIME, &t2);
  double t = (t2.tv_sec - t1.tv_sec) + (t2.tv_nsec - t1.tv_nsec)/1e9;
  fprintf(stderr, "Time %g GFlops/sec %g for compute\n", t, ops/t/1e9);
  double c = (t2.tv_sec - t0.tv_sec) + (t2.tv_nsec - t0.tv_nsec)/1e9;
  fprintf(stderr, "Time %g GFlops/sec %g for copy + compute\n", c, ops/c/1e9);

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      if (matE[i][j] != matC[i][j]) {
        fprintf(stderr, "Failed %d %d %x %x\n",i,j,*(int*)&matE[i][j],*(int*)&matC[i][j]);
        return 1;
      }
    }
  }
  fprintf(stderr, "Passed\n");
  return 0;
}
