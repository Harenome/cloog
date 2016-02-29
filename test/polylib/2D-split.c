/* Generated from 2D-split.scop by CLooG 0.18.4-16567be gmp bits in 0.01s. */
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
int c2, c4;
/* Original iterators. */
int i, j;

if ((H >= 0) && (W >= 0)) {
  for (c2=0;c2<=min(H-1,W-1);c2++) {
    for (c4=0;c4<=c2;c4++) {
      S1(c2, c4);
    }
  }
  for (c2=W;c2<=H-1;c2++) {
    for (c4=c2-W;c4<=c2;c4++) {
      S1(c2, c4);
    }
  }
  for (c2=H;c2<=W-1;c2++) {
    for (c4=0;c4<=H;c4++) {
      S1(c2, c4);
    }
  }
  for (c2=max(H,W);c2<=W+H;c2++) {
    for (c4=c2-W;c4<=H;c4++) {
      S1(c2, c4);
    }
  }
}
