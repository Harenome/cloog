/* Generated from test/polylib/3D-split.scop by CLooG 0.18.4-1e13234 gmp bits in 0.00s. */
/* Useful macros. */
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define min(x,y)    ((x) < (y) ? (x) : (y))

#ifdef TIME 
#define IF_TIME(foo) foo; 
#else
#define IF_TIME(foo)
#endif

/* Scattering iterators. */
int c2, c4, c6;
/* Original iterators. */
int i, j, k;

if ((D >= 0) && (H >= 0) && (W >= 0)) {
  for (c2=0;c2<=min(D+H,W+H);c2++) {
    for (c4=max(0,c2-W);c4<=min(H,c2);c4++) {
      for (c6=max(0,c2-D);c6<=min(H,c2);c6++) {
        S1 (c2, c4, c6);
      }
    }
  }
}
