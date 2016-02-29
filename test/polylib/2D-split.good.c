/* Generated from 2D-split.scop by CLooG 0.18.4-1e13234 gmp bits in 0.00s. */
extern void hash(int);

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

#define S1(i,j) { hash(1); hash(i); hash(j); }

void test(int W, int H)
{
  /* Scattering iterators. */
  int c2, c4;
  /* Original iterators. */
  int i, j;
  if ((H >= 0) && (W >= 0)) {
    for (c2=0;c2<=W+H;c2++) {
      for (c4=max(0,c2-W);c4<=min(H,c2);c4++) {
        S1(c2, c4);
      }
    }
  }
}
