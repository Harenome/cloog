/* Generated from 2D-tiled-split.scop by CLooG 0.18.4-1e13234 gmp bits in 0.01s. */
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

#define S1(i,j,k,l) { hash(1); hash(i); hash(j); hash(k); hash(l); }

void test(int W, int H)
{
  /* Scattering iterators. */
  int c2, c4, c6, c8;
  /* Original iterators. */
  int i, j, k, l;
  if ((H >= 0) && (W >= 0)) {
    for (c2=0;c2<=floord(W+H,32);c2++) {
      for (c4=max(0,ceild(32*c2-W-31,32));c4<=min(floord(H,32),c2);c4++) {
        for (c6=32*c2;c6<=min(min(W+H,32*c2+31),32*c4+W+31);c6++) {
          for (c8=max(32*c4,c6-W);c8<=min(min(H,c6),32*c4+31);c8++) {
            S1(c2,c4,c6,c8);
          }
        }
      }
    }
  }
}
