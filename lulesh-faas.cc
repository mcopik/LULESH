
#include <cstdio>
#include <cstdint>

typedef int32_t Int4_t ;
typedef int64_t Int8_t ;
typedef Int4_t  Index_t ; // array subscript and loop index
typedef double       real8 ;
typedef real8   Real_t ;  // floating point representation

static void CalcElemShapeFunctionDerivatives_faas( Real_t const x[],
                                       Real_t const y[],
                                       Real_t const z[],
                                       //Real_t b[][8],
                                       Real_t* b,
                                       Real_t* const volume );
static void CalcElemNodeNormals(Real_t pfx[8],
                         Real_t pfy[8],
                         Real_t pfz[8],
                         const Real_t x[8],
                         const Real_t y[8],
                         const Real_t z[8]);

extern "C" uint32_t empty(void* args, uint32_t size, void* res)
{
  int iters_remote = size / 24 / 8;
  //printf("Execute with size of input %d, iters %d, first val %f\n", size, remote_iters, static_cast<double*>(args)[1]);
  double* d_args = reinterpret_cast<double*>(args);
  double* o_args = reinterpret_cast<double*>(res);
  //printf("First vals %f %f %f %f of %d and %d iters\n", d_args[0], d_args[1], d_args[2], d_args[3], size, iters_remote);
  //int* src = static_cast<int*>(args), *dest = static_cast<int*>(res);
  //*dest = *src;
  Real_t* x_local = d_args;
  Real_t* y_local = d_args + 8*iters_remote;
  Real_t* z_local = d_args + 16*iters_remote;
  //// NEW OUTPUT - MEMCPY
  Real_t* send_determ = o_args;
  //// NEW OUTPUT 
  Real_t* B = o_args + iters_remote;
  for( Index_t k=0; k<iters_remote; ++k ) {
    // now invoke
    // Volume calculation involves extra work for numerical consistency
    CalcElemShapeFunctionDerivatives_faas(&x_local[8*k], &y_local[8*k], &z_local[8*k],
                                         &B[24*k], &send_determ[k]);

    CalcElemNodeNormals( &B[24*k] , &B[24*k+8], &B[24*k+16],
                          &x_local[8*k], &y_local[8*k], &z_local[8*k]);
  }
  //printf("First vals %f %f %f %f, write %d\n", o_args[0], o_args[1], o_args[2], o_args[3], size + iters_remote * 8);
  return size + iters_remote * 8;
}

static inline
void CalcElemShapeFunctionDerivatives_faas( Real_t const x[],
                                       Real_t const y[],
                                       Real_t const z[],
                                       //Real_t b[][8],
                                       Real_t* b,
                                       Real_t* const volume )
{
  const Real_t x0 = x[0] ;   const Real_t x1 = x[1] ;
  const Real_t x2 = x[2] ;   const Real_t x3 = x[3] ;
  const Real_t x4 = x[4] ;   const Real_t x5 = x[5] ;
  const Real_t x6 = x[6] ;   const Real_t x7 = x[7] ;

  const Real_t y0 = y[0] ;   const Real_t y1 = y[1] ;
  const Real_t y2 = y[2] ;   const Real_t y3 = y[3] ;
  const Real_t y4 = y[4] ;   const Real_t y5 = y[5] ;
  const Real_t y6 = y[6] ;   const Real_t y7 = y[7] ;

  const Real_t z0 = z[0] ;   const Real_t z1 = z[1] ;
  const Real_t z2 = z[2] ;   const Real_t z3 = z[3] ;
  const Real_t z4 = z[4] ;   const Real_t z5 = z[5] ;
  const Real_t z6 = z[6] ;   const Real_t z7 = z[7] ;

  Real_t fjxxi, fjxet, fjxze;
  Real_t fjyxi, fjyet, fjyze;
  Real_t fjzxi, fjzet, fjzze;
  Real_t cjxxi, cjxet, cjxze;
  Real_t cjyxi, cjyet, cjyze;
  Real_t cjzxi, cjzet, cjzze;

  fjxxi = Real_t(.125) * ( (x6-x0) + (x5-x3) - (x7-x1) - (x4-x2) );
  fjxet = Real_t(.125) * ( (x6-x0) - (x5-x3) + (x7-x1) - (x4-x2) );
  fjxze = Real_t(.125) * ( (x6-x0) + (x5-x3) + (x7-x1) + (x4-x2) );

  fjyxi = Real_t(.125) * ( (y6-y0) + (y5-y3) - (y7-y1) - (y4-y2) );
  fjyet = Real_t(.125) * ( (y6-y0) - (y5-y3) + (y7-y1) - (y4-y2) );
  fjyze = Real_t(.125) * ( (y6-y0) + (y5-y3) + (y7-y1) + (y4-y2) );

  fjzxi = Real_t(.125) * ( (z6-z0) + (z5-z3) - (z7-z1) - (z4-z2) );
  fjzet = Real_t(.125) * ( (z6-z0) - (z5-z3) + (z7-z1) - (z4-z2) );
  fjzze = Real_t(.125) * ( (z6-z0) + (z5-z3) + (z7-z1) + (z4-z2) );

  /* compute cofactors */
  cjxxi =    (fjyet * fjzze) - (fjzet * fjyze);
  cjxet =  - (fjyxi * fjzze) + (fjzxi * fjyze);
  cjxze =    (fjyxi * fjzet) - (fjzxi * fjyet);

  cjyxi =  - (fjxet * fjzze) + (fjzet * fjxze);
  cjyet =    (fjxxi * fjzze) - (fjzxi * fjxze);
  cjyze =  - (fjxxi * fjzet) + (fjzxi * fjxet);

  cjzxi =    (fjxet * fjyze) - (fjyet * fjxze);
  cjzet =  - (fjxxi * fjyze) + (fjyxi * fjxze);
  cjzze =    (fjxxi * fjyet) - (fjyxi * fjxet);

  /* calculate partials :
     this need only be done for l = 0,1,2,3   since , by symmetry ,
     (6,7,4,5) = - (0,1,2,3) .
  */
  //b[0][0] =   -  cjxxi  -  cjxet  -  cjxze;
  //b[0][1] =      cjxxi  -  cjxet  -  cjxze;
  //b[0][2] =      cjxxi  +  cjxet  -  cjxze;
  //b[0][3] =   -  cjxxi  +  cjxet  -  cjxze;
  //b[0][4] = -b[0][2];
  //b[0][5] = -b[0][3];
  //b[0][6] = -b[0][0];
  //b[0][7] = -b[0][1];

  //b[1][0] =   -  cjyxi  -  cjyet  -  cjyze;
  //b[1][1] =      cjyxi  -  cjyet  -  cjyze;
  //b[1][2] =      cjyxi  +  cjyet  -  cjyze;
  //b[1][3] =   -  cjyxi  +  cjyet  -  cjyze;
  //b[1][4] = -b[1][2];
  //b[1][5] = -b[1][3];
  //b[1][6] = -b[1][0];
  //b[1][7] = -b[1][1];

  //b[2][0] =   -  cjzxi  -  cjzet  -  cjzze;
  //b[2][1] =      cjzxi  -  cjzet  -  cjzze;
  //b[2][2] =      cjzxi  +  cjzet  -  cjzze;
  //b[2][3] =   -  cjzxi  +  cjzet  -  cjzze;
  //b[2][4] = -b[2][2];
  //b[2][5] = -b[2][3];
  //b[2][6] = -b[2][0];
  //b[2][7] = -b[2][1];
  b[0] =   -  cjxxi  -  cjxet  -  cjxze;
  b[1] =      cjxxi  -  cjxet  -  cjxze;
  b[2] =      cjxxi  +  cjxet  -  cjxze;
  b[3] =   -  cjxxi  +  cjxet  -  cjxze;
  b[4] = -b[2];
  b[5] = -b[3];
  b[6] = -b[0];
  b[7] = -b[1];

  b[8] =   -  cjyxi  -  cjyet  -  cjyze;
  b[8+1] =      cjyxi  -  cjyet  -  cjyze;
  b[8+2] =      cjyxi  +  cjyet  -  cjyze;
  b[8+3] =   -  cjyxi  +  cjyet  -  cjyze;
  b[8+4] = -b[8+2];
  b[8+5] = -b[8+3];
  b[8+6] = -b[8+0];
  b[8+7] = -b[8+1];

  b[16+0] =   -  cjzxi  -  cjzet  -  cjzze;
  b[16+1] =      cjzxi  -  cjzet  -  cjzze;
  b[16+2] =      cjzxi  +  cjzet  -  cjzze;
  b[16+3] =   -  cjzxi  +  cjzet  -  cjzze;
  b[16+4] = -b[16+2];
  b[16+5] = -b[16+3];
  b[16+6] = -b[16+0];
  b[16+7] = -b[16+1];

  /* calculate jacobian determinant (volume) */
  *volume = Real_t(8.) * ( fjxet * cjxet + fjyet * cjyet + fjzet * cjzet);
}

static inline
void SumElemFaceNormal(Real_t *normalX0, Real_t *normalY0, Real_t *normalZ0,
                       Real_t *normalX1, Real_t *normalY1, Real_t *normalZ1,
                       Real_t *normalX2, Real_t *normalY2, Real_t *normalZ2,
                       Real_t *normalX3, Real_t *normalY3, Real_t *normalZ3,
                       const Real_t x0, const Real_t y0, const Real_t z0,
                       const Real_t x1, const Real_t y1, const Real_t z1,
                       const Real_t x2, const Real_t y2, const Real_t z2,
                       const Real_t x3, const Real_t y3, const Real_t z3)
{
   Real_t bisectX0 = Real_t(0.5) * (x3 + x2 - x1 - x0);
   Real_t bisectY0 = Real_t(0.5) * (y3 + y2 - y1 - y0);
   Real_t bisectZ0 = Real_t(0.5) * (z3 + z2 - z1 - z0);
   Real_t bisectX1 = Real_t(0.5) * (x2 + x1 - x3 - x0);
   Real_t bisectY1 = Real_t(0.5) * (y2 + y1 - y3 - y0);
   Real_t bisectZ1 = Real_t(0.5) * (z2 + z1 - z3 - z0);
   Real_t areaX = Real_t(0.25) * (bisectY0 * bisectZ1 - bisectZ0 * bisectY1);
   Real_t areaY = Real_t(0.25) * (bisectZ0 * bisectX1 - bisectX0 * bisectZ1);
   Real_t areaZ = Real_t(0.25) * (bisectX0 * bisectY1 - bisectY0 * bisectX1);

   *normalX0 += areaX;
   *normalX1 += areaX;
   *normalX2 += areaX;
   *normalX3 += areaX;

   *normalY0 += areaY;
   *normalY1 += areaY;
   *normalY2 += areaY;
   *normalY3 += areaY;

   *normalZ0 += areaZ;
   *normalZ1 += areaZ;
   *normalZ2 += areaZ;
   *normalZ3 += areaZ;
}

static inline
void CalcElemNodeNormals(Real_t pfx[8],
                         Real_t pfy[8],
                         Real_t pfz[8],
                         const Real_t x[8],
                         const Real_t y[8],
                         const Real_t z[8])
{
   for (Index_t i = 0 ; i < 8 ; ++i) {
      pfx[i] = Real_t(0.0);
      pfy[i] = Real_t(0.0);
      pfz[i] = Real_t(0.0);
   }
   /* evaluate face one: nodes 0, 1, 2, 3 */
   SumElemFaceNormal(&pfx[0], &pfy[0], &pfz[0],
                  &pfx[1], &pfy[1], &pfz[1],
                  &pfx[2], &pfy[2], &pfz[2],
                  &pfx[3], &pfy[3], &pfz[3],
                  x[0], y[0], z[0], x[1], y[1], z[1],
                  x[2], y[2], z[2], x[3], y[3], z[3]);
   /* evaluate face two: nodes 0, 4, 5, 1 */
   SumElemFaceNormal(&pfx[0], &pfy[0], &pfz[0],
                  &pfx[4], &pfy[4], &pfz[4],
                  &pfx[5], &pfy[5], &pfz[5],
                  &pfx[1], &pfy[1], &pfz[1],
                  x[0], y[0], z[0], x[4], y[4], z[4],
                  x[5], y[5], z[5], x[1], y[1], z[1]);
   /* evaluate face three: nodes 1, 5, 6, 2 */
   SumElemFaceNormal(&pfx[1], &pfy[1], &pfz[1],
                  &pfx[5], &pfy[5], &pfz[5],
                  &pfx[6], &pfy[6], &pfz[6],
                  &pfx[2], &pfy[2], &pfz[2],
                  x[1], y[1], z[1], x[5], y[5], z[5],
                  x[6], y[6], z[6], x[2], y[2], z[2]);
   /* evaluate face four: nodes 2, 6, 7, 3 */
   SumElemFaceNormal(&pfx[2], &pfy[2], &pfz[2],
                  &pfx[6], &pfy[6], &pfz[6],
                  &pfx[7], &pfy[7], &pfz[7],
                  &pfx[3], &pfy[3], &pfz[3],
                  x[2], y[2], z[2], x[6], y[6], z[6],
                  x[7], y[7], z[7], x[3], y[3], z[3]);
   /* evaluate face five: nodes 3, 7, 4, 0 */
   SumElemFaceNormal(&pfx[3], &pfy[3], &pfz[3],
                  &pfx[7], &pfy[7], &pfz[7],
                  &pfx[4], &pfy[4], &pfz[4],
                  &pfx[0], &pfy[0], &pfz[0],
                  x[3], y[3], z[3], x[7], y[7], z[7],
                  x[4], y[4], z[4], x[0], y[0], z[0]);
   /* evaluate face six: nodes 4, 7, 6, 5 */
   SumElemFaceNormal(&pfx[4], &pfy[4], &pfz[4],
                  &pfx[7], &pfy[7], &pfz[7],
                  &pfx[6], &pfy[6], &pfz[6],
                  &pfx[5], &pfy[5], &pfz[5],
                  x[4], y[4], z[4], x[7], y[7], z[7],
                  x[6], y[6], z[6], x[5], y[5], z[5]);
}
