/* $Id: nlmefit.c,v 1.43 1999/03/05 04:15:38 pinheiro Exp $ 

   Routines for calculation of the log-likelihood or restricted
   log-likelihood with mixed-effects models.

   Copyright 1997, 1999 Douglas M. Bates <bates@stat.wisc.edu>,
                        Jose C. Pinheiro <jcp@research.bell-labs.com>

   This file is part of the nlme library for S and related languages
   and is made available under the terms of the GNU General Public
   License, version 2, or at your option, any later version,
   incorporated herein by reference.

   This program is distributed in the hope that it will be
   useful, but WITHOUT ANY WARRANTY; without even the implied
   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
   PURPOSE.  See the GNU General Public License for more
   details.

   You should have received a copy of the GNU General Public
   License along with this program; if not, write to the Free
   Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
   MA 02111-1307, USA
 
*/

#include "S.h"								     
#ifndef R_S_H
#include "nonlin.h"
#endif /* R_S_H */
#include <stdlib.h>

#define DNULLP (double *) 0
				
#ifndef SPLUS_VERSION		/* F77_CALL and F77_NAME used as in S-PLUS */
#ifdef S_VERSION		/* S VERSION 4 */
#define F77_CALL F77_SUB
#define F77_NAME F77_SUB
#define msmnh    dmnh
#endif /* S_VERSION */
#endif /* SPLUS_VERSION */

#ifdef R_S_H
#include "S_compat.h"
#else
extern void F77_NAME(dqrdca)();
extern void F77_NAME(dtrsl)();
extern void F77_NAME(dqrsl)();
#define longint long int
#endif /* R_S_H */

extern void F77_NAME(chol)();
extern void F77_NAME(rs)();
extern void F77_NAME(msmnh)();

typedef struct dim_struct {
  longint
    N,				/* number of observations in original data */
    ZXrows,			/* number of rows in ZXy  */
    ZXcols,			/* number of columns in ZXy */
    Q,				/* number of levels of random effects */
    Srows,			/* number of rows in decomposition */
    *q,				/* dimensions of the random effects */
    *ngrp,			/* number of groups at each level */
    *DmOff,			/* offsets into the DmHalf array */
    *ncol,			/* no. of columns decomposed at each level */
    *nrot,			/* no. of columns rotated at each level */
    **ZXoff,			/* offsets into ZXy */
    **ZXlen,			/* groups lengths */
    **SToff,			/* offsets into storage */
    **DecOff,			/* offsets into decomposition */
    **DecLen;			/* decomposition group lengths */
} *dimPTR;

static longint **
setOffsets(longint ** base, longint * ngrp, longint Qp2)
{
  longint i, **ptrVec = Calloc((size_t) Qp2, longint *);
  for (i = 0; i < Qp2; i++) {
    ptrVec[i] = *base;
    *base += ngrp[i];
  }
  return ptrVec;
}

static dimPTR
dims(longint *pdims)
{				/* constructor for a dims object */
  dimPTR value = Calloc((size_t) 1, struct dim_struct);
  longint *base, Qp2, *ngrp;

  value->N = (int) pdims[0];
  value->ZXrows = pdims[1];
  value->ZXcols = pdims[2];
  value->Q = pdims[3];
  Qp2 = (value->Q) + 2;
  value->Srows = pdims[4];
  value->q = pdims + 5;
  ngrp = value->ngrp = value->q + Qp2;
  value->DmOff = value->ngrp + Qp2;
  value->ncol = value->DmOff + Qp2;
  value->nrot = value->ncol + Qp2;
  base = value->nrot + Qp2;
  value->ZXoff = setOffsets(&base, ngrp, Qp2);
  value->ZXlen = setOffsets(&base, ngrp, Qp2); 
  value->SToff = setOffsets(&base, ngrp, Qp2);
  value->DecOff = setOffsets(&base, ngrp, Qp2);
  value->DecLen = setOffsets(&base, ngrp, Qp2);
  return value;
}

static void
dimFree(dimPTR this)
{
  Free(this->DecOff);
  Free(this->SToff);
  Free(this->ZXlen);
  Free(this->ZXoff);
  Free(this);
}

static size_t
count_DmHalf_pars( dimPTR dd, longint *pdClass )
{
  int i, result;
  for ( i = 0, result = 0; i < dd->Q; i++ ) {
    switch( pdClass[ i ] ) {
    case 0: result += ( (dd->q)[ i ] * ( (dd->q)[ i ] + 1 ) ) / 2;
      break;
    case 1: result += (dd->q)[ i ];
      break;
    case 2: result += 1;
      break;
    case 3: result += 2;
    }
  }
  return (size_t) result;
}
  
extern void compSymm_pd( double*, longint*, double* );
extern void matrixLog_pd( double*, longint*, double* );
static void Delta2MatrixLog( double*, longint*, double* );

static double *
generate_DmHalf( double *DmHalf, dimPTR dd, longint *pdClass, double *pars )
{				/* Expand parameters to DmHalf arrays */
  int i, j, q, Q = dd->Q; double diag;
  for (i = 0; i < Q; i++) {
    q = (dd->q)[ i ];
    switch (pdClass[i]) {
    case 0:			/* default: unstructured */
      matrixLog_pd( DmHalf + (dd->DmOff)[ i ], dd->q + i, pars );
      pars += (q * (q + 1))/2;
      break;
    case 1:			/* diagonal */
      for (j = 0; j < q; j++) {
	DmHalf[ (dd->DmOff)[i] + j * (q + 1) ] = exp( *pars++ );
      }
      break;
    case 2:			/* multiple of identity */
      diag = exp( *pars );
      for (j = 0; j < q; j++) {
	DmHalf[ (dd->DmOff)[i] + j * (q + 1) ] = diag;
      }
      pars++;
      break;
    case 3:			/* compound symmetry */
      compSymm_pd( DmHalf + (dd->DmOff)[ i ], dd->q + i, pars );
      pars += 2;
      break;
    }
  }
  return DmHalf;
}

static double *
generate_theta( double *theta, dimPTR dd, longint *pdClass, double *DmHalf )
{				/* Expand parameters to DmHalf arrays */
  int i, j, q, Q = dd->Q;
  for (i = 0; i < Q; i++) {
    q = (dd->q)[ i ];
    switch (pdClass[i]) {
    case 0:			/* default: unstructured */
      Delta2MatrixLog( theta, dd->q + i, DmHalf + (dd->DmOff)[ i ] );
      theta += (q * (q + 1))/2;
      break;
    case 1:			/* diagonal */
      for (j = 0; j < q; j++) {
	*theta++ = log( DmHalf[ (dd->DmOff)[i] + j * (q + 1) ] );
      }
      break;
    case 2:			/* multiple of identity */
      *theta++ = log( DmHalf[(dd->DmOff)[i]] );
      break;
    case 3:			/* compound symmetry */
      PROBLEM "Haven't written the compound symmetry case for this yet"
	RECOVER(NULL_ENTRY);
      break;
    }
  }
  return theta;
}

static void
d_axpy(double *y, double a, double *x, longint n)
{				/* y <- a * x + y  */
  while (n-- > 0) { *y++ += a * *x++; }
}

static double
d_sum_sqr( double *x, longint n )
{				/* sum(x * x) */
  double accum = 0.0;
  while (n-- > 0) { accum += *x * *x; x++; }
  return accum;
}

static double
d_dot_prod( double *x, longint incx, double *y, longint incy, longint n )
{				/* sum(x * y) */
  double accum = 0.0;
  while (n-- > 0) { accum += *x * *y; x +=incx; y += incy; }
  return accum;
}

#ifdef DEBUG
static void
print_mat( char *msg, double *x, longint ldx, longint nrow,
	   longint ncol )
{				/* print matrix and message */
  int i, j;
  printf( "%s\n", msg );
  for (i = 0; i < nrow; i++) {
    for (j = 0; j < ncol; j++) {
      printf( " %10.5g", x[i + j * ldx ] );
    }
    printf( "\n" );
  }
  printf( "\n" );
}
#endif /* DEBUG */

static double *
copy_mat(double *y, longint ldy, double *x, longint ldx,
	 longint nrow, longint ncol) 
{				/* y <- x */
  double * ret = y;
  while (ncol-- > 0) { Memcpy(y, x, nrow); y += ldy; x += ldx; }
  return ret;
}

static double *
copy_trans(double *y, longint ldy, double *x, longint ldx,
	   longint nrow, longint ncol) /* y <- t(x) */
{
  double * ret = y;
  longint i, j;
  for (i = 0L; i < nrow; i++) {
    for (j = 0L; j < ncol; j++) { y[j] = x[i + j * ldx]; }
    y += ldy;
  }
  return ret;
}
  
static double *
scale_mat(double *y, longint ldy, double a,
	  double *x, longint ldx, longint nrow, longint ncol)
{				/* y <- a * x */
  int i, j;
  double * ret = y;

  for (j = 0; j < ncol; j++) {
    for (i = 0; i < nrow; i++) { y[i] = a * x[i]; }
    y += ldy; x += ldx;
  }
  return ret;
}

static double *
plus_equals_mat(double *y, longint ldy, double *x, longint ldx,
		longint nrow, longint ncol)
{				/* y <- y + x */
  double * ret = y;
  int i, j;

  for (j = 0; j < ncol; j++) {
    for (i = 0; i < nrow; i++) { y[i] += x[i]; }
    y += ldy; x += ldx;
  }
  return ret;
}
      
static double *
mult_mat(double *z, longint ldz,
	 double *x, longint ldx, longint xrows, longint xcols,
	 double *y, longint ldy, longint ycols) 
{				/* z <- x %*% y */
  double *t, *tmp = Calloc((size_t)(xrows * ycols), double);
  int i, j;			/* use tmp so z can be either x or y */

  t = tmp;
  for (i = 0; i < ycols; i++) {
    for (j = 0; j < xcols; j++) {
      d_axpy(t, y[j], x + j * ldx, xrows);
    }
    t += xrows;
    y += ldy;
  }
  copy_mat(z, ldz, tmp, xrows, xrows, ycols);
  Free(tmp);
  return z;
}

static double *
crossprod_mat(double *y, longint ldy, double *x, longint ldx,
	      longint nrow, longint ncol) /* y <- t(x) %*% x */
{
  longint i, j;

  for( i = 0; i < ncol; i++ ) {
    y[ i * ldy + i ] = d_dot_prod( x + i * ldx, 1L, x + i * ldx, 1L, nrow );
    for( j = 0; j < i; j++) {
      y[ i * ldy + j ] = y[ j * ldy + i ] =
	d_dot_prod( x + i * ldx, 1L, x + j * ldx, 1L, nrow );
    }
  }
  return y;
}
  
static void
zero_mat(double *y, longint ldy, longint nrow, longint ncol)
{				/* y[,] <- 0 */
  while (ncol-- > 0) {
    int i;
    for (i = 0; i < nrow; i++) { y[i] = 0.0; }
    y += ldy;
  }
}
      
static double sqrt_eps = 0., xlower = 0.;

typedef struct QR_struct {
  double *mat, *qraux;
  longint *pivot, rank, ldmat, nrow, ncol;
} *QRptr;

static QRptr
QR(double *mat, longint ldmat, longint nrow, longint ncol)
{				/* Constructor for a QR object */
  QRptr value = Calloc((size_t) 1, struct QR_struct);
  longint j;  double *work;

  if (sqrt_eps == 0.) { sqrt_eps = sqrt(DOUBLE_EPS); }
  value->mat = mat;
  value->ldmat = ldmat;
  value->nrow = nrow;
  value->ncol = ncol;
  value->qraux = Calloc((size_t) ncol, double);
  value->pivot = Calloc((size_t) ncol, longint);
  for (j = 0; j < ncol; j++) { (value->pivot)[j] = j; }
  work = Calloc( 2 * ncol, double );
  F77_CALL(dqrdca) (mat, &ldmat, &nrow, &ncol, value->qraux, value->pivot,
		    work, &(value->rank), &sqrt_eps);
  Free(work);
  return value;
}

static void
QRfree(QRptr this)
{				/* destructor for a QR object*/
  Free(this->pivot); 
  Free(this->qraux); 
  Free(this);
}

static longint
QRqty(QRptr this, double *ymat, longint ldy, longint ycol)
{				/* ymat <- qr.qty(this, ymat) */
  longint j, info, task = 1000L;
  for (j = 0; j < ycol; j++) {
    double *col = ymat + j * ldy;
    F77_CALL(dqrsl) (this->mat, &(this->ldmat), &(this->nrow), &(this->ncol),
		     this->qraux, col, DNULLP, col, DNULLP, DNULLP, DNULLP,
		     &task, &info);
  }
  return info;
}

static longint
QRsolve( QRptr this, double *ymat, longint ldy, longint ycol,
	 double *beta, longint ldbeta )
{				/* beta <- qr.beta(this, ymat) */
  longint j, info, task = 1100L;
  double *qty = Calloc( this->nrow, double ),
    *bb = Calloc( this->ncol, double );

  for (j = 0; j < ycol; j++) {
    Memcpy( qty, ymat, this->nrow );
    F77_CALL(dqrsl) (this->mat, &(this->ldmat), &(this->nrow), &(this->ncol),
		     this->qraux, qty, DNULLP, qty, bb, DNULLP,
		     DNULLP, &task, &info);
    Memcpy( beta, bb, this->ncol );
    ymat += ldy;
    beta += ldbeta;
  }
  Free( qty ); Free( bb );
  return info;
}

static double
QRlogAbsDet(QRptr this)
{				/* log(abs(det(upper triangle))) */
  longint j;
  double accum = 0.0;
  for (j = 0; j < this->rank; j++)
    accum += log(fabs(this->mat[j * (this->ldmat + 1L)]));
  return accum;
}

static void
QRstoreR(QRptr this, double *dest, longint ldDest)
{				/* store the R part into dest */
  int i;
  for (i = 0; i < this->ncol; i++) {
    Memcpy(dest + this->pivot[i] * ldDest, this->mat + i * this->ldmat, 
	   ((i + 1) > this->rank) ? this->rank : i + 1);
  }
}

static longint
QR_and_rotate(double *mat, longint ldmat, longint nrow, longint ncol,
              double *DmHalf, longint qi, longint ndecomp,
              double *logdet, double *store, longint ldstr)
     /* Append DmHalf to the bottom of mat and take a QR decomposition
        of the first ndecomp columns.  Apply the rotations to the other
        columns.  Return the rank and increment log(abs(det(R11))). */
{
  longint rank, arow = nrow + qi,  /* number of rows in augmented matrix */ 
    ndrow = ((arow > ndecomp) ? ndecomp : arow); 
  double *aug = Calloc((size_t) arow * ncol, double);
  QRptr aQR;

  copy_mat(aug, arow, mat, ldmat, nrow, ncol);
  copy_mat(aug + nrow, arow, DmHalf, qi, qi, qi);
  aQR = QR(aug, arow, arow, ndecomp);
  if (logdet != DNULLP) { *logdet += QRlogAbsDet(aQR); }
  QRqty(aQR, aug + ndecomp * arow, arow, ncol - ndecomp);
  if (ldstr > 0) {
    QRstoreR(aQR, store, ldstr);
    copy_mat(store + ndecomp * ldstr, ldstr, aug + ndecomp * arow,
             arow, ndrow, ncol - ndecomp);
  }
  if (qi < ndecomp) { zero_mat(mat, ldmat, nrow, ncol); }
  copy_mat(mat + ndecomp * ldmat, ldmat, aug + ndecomp * (arow + 1L),
           arow, arow - ndrow, ncol - ndecomp);
  rank = aQR->rank;
  QRfree(aQR); Free(aug);
  return rank;
}

static longint			/* backsolve and update */
backsolve(double *mat, longint ldmat, longint nupdate, longint ncol,
	  longint nrot, longint ny)
{
  longint i, j, ONE = 1L, info;
  double *y = mat + (int) ((ncol + nrot - ny) * ldmat);

  mat = mat - (int) nupdate;
  for (i = 0; i < ny; i++) {	/* usually ny = 1 but just in case ... */
    F77_CALL(dtrsl) (mat + (int) nupdate, &ldmat, &ncol, y, &ONE, &info);
    if (info != 0) {
      return info;
    }
    for (j = 0; j < ncol; j++) {
      d_axpy(y - (int) nupdate, - y[j], mat + (int) (j * ldmat), nupdate);
    }
    y += (int) ldmat;
  }
  return info;
}

static longint			/* invert an upper-triangular matrix in place*/
invert_upper(double *mat, longint ldmat, longint ncol)
{
  longint i, j, ONE = 1L, info = 0L;
  double *b = Calloc((size_t) ncol, double);

  for (i = ncol; i > 1L; i--) {
    for (j = 0; j < (i - 1L); j++) { b[j] = 0.0; }
    b[((int) i) - 1] = 1.0;
    F77_CALL(dtrsl) (mat, &ldmat, &i, b, &ONE, &info);
    if (info != 0) { Free(b); return info; }
    Memcpy(mat + (i - 1) * ldmat, b, (int) i); 
  }
  if (*mat == 0.0) { Free(b); return 1L; }
  *mat = 1.0 / (*mat);
  Free(b); return 0L;
}

static longint			/* invert a block in the virtual R array */
invert_block(double *mat, longint ldmat, longint nabove,
	     longint ncol, longint nright)
{
  double * tpblk = mat - (int) nabove;
  longint info = invert_upper(mat, ldmat, ncol);

  if (info != 0L) return info;
  if (nright > 0) {
    double *ntri = Calloc((size_t) (ncol * ncol), double),
      *rtblk = mat + ncol * ldmat;
    scale_mat(ntri, ncol, -1.0, mat, ldmat, ncol, ncol);
    mult_mat(rtblk, ldmat, ntri, ncol, ncol, ncol, rtblk, ldmat, nright);
    Free(ntri);
    if (nabove > 0L) {
      double *tmp = Calloc((size_t)(nabove * nright), double);
      plus_equals_mat(rtblk - (size_t)nabove, ldmat,
		      mult_mat(tmp, nabove, tpblk, ldmat, nabove, ncol,
			       rtblk, ldmat, nright),
		      nabove, nabove, nright);
      Free(tmp);
    }
  }
  if (nabove > 0L) {
    mult_mat(tpblk, ldmat, tpblk, ldmat, nabove, ncol, mat, ldmat, ncol);
  }
  return 0L;
}

void				/* return the decomposition for ZXy */
mixed_decomp(double *ZXy, longint *pdims)
{
  dimPTR dd = dims(pdims);	/* Create a dimensions structure */
  longint i, j, Qp2 = (dd->Q) + 2;
  double *dc = Calloc((size_t) ((dd->Srows) * (dd->ZXcols)), double);

  if ((dd->Srows) < (dd->ZXrows)) { /* decomposition is worthwhile */
    for (i = 0; i < Qp2; i++) {
      for(j = 0; j < (dd->ngrp)[i]; j++) {
	QR_and_rotate(ZXy + (dd->ZXoff)[i][j], dd->ZXrows, (dd->ZXlen)[i][j],
		      (dd->ncol)[i] + (dd->nrot)[i], DNULLP, 0L,
		      (dd->ncol)[i], DNULLP, dc + (dd->SToff)[i][j],
		      dd->Srows);
      }
    }
    Memcpy(ZXy, dc, dd->Srows * dd->ZXcols);
  }
  Free(dc); dimFree(dd);
}

static double			/* evaluate the log-likelihood pieces */
internal_loglik(dimPTR dd, double *ZXy, double *DmHalf, longint *RML,
		double *dc, double *lRSS)	
{				/* if dc is NULL, don't attempt storage */
  longint i, j, Q = dd->Q,  Qp2 = Q + 2, qi,
    ldstr = (dc != DNULLP) ? (dd->Srows) : 0L;
  double accum, *dmHlf, *lglk = Calloc( Qp2, double );
  QRptr dmQR;

  for (i = 0; i < Qp2; i++) {
    qi = (dd->q)[i];
    for (j = 0; j < (dd->ngrp)[i]; j++) {
      if (qi >
	  QR_and_rotate(ZXy + (dd->ZXoff)[i][j], dd->ZXrows,
			(dd->ZXlen)[i][j], (dd->ncol)[i] + (dd->nrot)[i],
			DmHalf + (dd->DmOff)[i], qi, (dd->ncol)[i],
			lglk + i, dc + (dd->SToff)[i][j], ldstr))
	{ PROBLEM "Singular precision matrix in level %ld, block %ld",
		  i + 1L, j + 1L  RECOVER(NULL_ENTRY); }
    }
  }
  for(i = 0, accum = 0; i < Q; i++) {
    qi = (dd->q)[i];
    dmHlf = Calloc( (size_t) qi * qi, double );
    dmQR = QR( copy_mat( dmHlf, qi, DmHalf + (dd->DmOff)[i],
			 qi, qi, qi ), qi, qi, qi);
    accum += (dd->ngrp)[i] * QRlogAbsDet( dmQR ) - lglk[i];
    QRfree( dmQR ); Free( dmHlf );
  }
  accum -= *RML * lglk[ Q ] + (dd->N - *RML * dd->ncol[ Q ]) * lglk[Q + 1];
  if (lRSS != DNULLP) *lRSS = lglk[Q + 1]; /* return log(RSS)/2 */
  Free( lglk );
  return accum;
}

static void
internal_estimate(dimPTR dd, double *dc)
{				/* solve for Beta and b_i estimates */
  longint i, j, Qp1 = (dd->Q) + 1L;

  for (i = (dd->Q); i >= 0; i--) {
    for (j = 0; j < (dd->ngrp)[i]; j++) {
      if (backsolve(dc + (dd->SToff)[i][j], dd->Srows,
		    (dd->SToff)[i][j] - (dd->DecOff)[i][j],
		    (dd->ncol)[i], (dd->nrot)[i], (dd->ncol)[Qp1])
	  != 0)
	{ PROBLEM "Singularity in backsolve at level %ld, block %ld",
		  i + 1L, j + 1L RECOVER(NULL_ENTRY); }
    }
  }
}

static void
internal_R_invert(dimPTR dd, double *dc)
{				/* Invert the virtual R matrix in place */
  int i, j;
  for (i = (dd->Q); i >= 0; i--) {
    for (j = 0; j < (dd->ngrp)[i]; j++) {
      invert_block(dc + (dd->SToff)[i][j], dd->Srows,
		   (dd->SToff)[i][j] - (dd->DecOff)[i][j],
		   (dd->ncol)[i], (dd->nrot)[i] - 1L);
    }
  }
}

static double cube_root_eps = 0.;

static double *
pt_prod( double *prod, double *a, double *b, int len )
{				/* prod <- a * b */
  int i; double *ret = prod;
  for (i = 0; i < len; i++) {
    *prod++ = *a++ * *b++;
  }
  return ret;
}

static void
finite_diff_Hess( double (*func)(double*), double *pars, int npar,
		  double *vals )
{				/* use Koshal design for finite-differences */
  int i, j, nTot = 1 + npar + ( npar * ( npar + 1 ) ) / 2;
  double *incr = Calloc( npar, double), *ppt, *xpt, *dpt,
    *parray = Calloc( nTot * npar, double ), /* array of parameters */
    *div = Calloc( nTot, double ), /* divisors */
    *Xmat = Calloc( nTot * nTot, double ); /* regressor matrix */
  QRptr xQR;

  if (!cube_root_eps) cube_root_eps = exp( log( DOUBLE_EPS ) / 3.);
  div[ 0 ] = 1.0;
  ppt = parray + npar * ( 2 * npar + 1 ); /* location of first cross term */
  xpt = Xmat + nTot * ( 2 * npar + 1 );	/* location of first cross column */
  dpt = div + 2 * npar + 1;
  for (i = 0; i < npar; i++) {
    incr[i] = (pars[ i ] != 0.0) ? cube_root_eps * pars[ i ] : cube_root_eps;
    div[ i + 1 ] = 1.0 / incr[ i ];
    div[ npar + i + 1 ] = 2.0 / ( incr[ i ] * incr[ i ] );
    parray[ npar + i * (npar + 1) ] = 1.0;
    parray[ (npar + i) * (npar + 1) ] = -1.0;
    for (j = i + 1; j < npar; j++) {
      ppt[ i ] = ppt[ j ] = 1;
      ppt += npar;
    }
    for (j = 0; j < nTot; j++) {
      Xmat[ j + (i + 1) * nTot ] = parray[ i + j * npar ];
    }
    pt_prod( Xmat + (npar + i + 1) * nTot, Xmat + (i + 1) * nTot,
		       Xmat + (i + 1) * nTot, nTot );
    for (j = 0; j < i; j++) {
      pt_prod( xpt, Xmat + (i + 1) * nTot, Xmat + (j + 1) * nTot, nTot );
      xpt += nTot;
      *dpt++ = 1.0 / ( incr[ i ] * incr[ j ] );
    }
  }
#ifdef DEBUG
  print_mat( "parray", parray, npar, npar, nTot );
#endif /* DEBUG */
  vals[ 0 ] = (*func)( pars );
  Xmat[ 0 ] = 1.0;
  for (i = 1; i < nTot; i++) {
    Xmat[i] = 1.0;		/* column of 1's for constant */
    Memcpy( parray, pars, npar );
    for (j = 0; j < npar; j++) {
      parray[ j ] += parray[ j + i * npar ] * incr[ j ];
    }
    vals[i] = (*func)( parray );
  }
#ifdef DEBUG
  print_mat( "Xmat", Xmat, nTot, nTot, nTot );
#endif /* DEBUG */
  xQR = QR( Xmat, (longint) nTot, (longint) nTot, (longint) nTot );
  QRsolve( xQR, vals, nTot, 1, vals, nTot );
  pt_prod( vals, vals, div, nTot );
				/* re-arrange the Hessian terms */
  xpt = vals + npar + 1;
  Memcpy( div, vals + npar + 1, nTot - ( npar + 1 ) );
  dpt = div + npar;		/* first off-diagonal */
  for (i = 0; i < npar; i++) {
    xpt[ i * ( npar + 1 ) ] = div[ i ];	/* diagonals */
    for (j = 0; j < i; j++) {
      xpt[ i + j * npar ] = xpt[ j + i * npar ] = *dpt++;
    }
  }
  QRfree( xQR ); Free( incr ); Free( parray ); Free( div ); Free( Xmat );
  return;
}

				/* In gcc we can use nested function
				   definitions but not for other compilers */
static double *zxcopy, *zxcopy2, *Delta, *values;
static dimPTR dd;
static longint *setngs, *pdC;
size_t zxdim;

static double 
logLik_fun( double *pars )
{				/* defined for finite differences */
  Memcpy( zxcopy2, zxcopy, zxdim );
  return internal_loglik(dd, zxcopy2, generate_DmHalf( Delta, dd, pdC, pars ),
			 setngs, DNULLP, DNULLP );
}

static double 
negLogLik_fun( double *pars )
{				/* defined for finite differences */
  Memcpy( zxcopy2, zxcopy, zxdim );
  return - internal_loglik(dd, zxcopy2, generate_DmHalf( Delta, dd, pdC, pars ),
			   setngs, DNULLP, DNULLP );
}

void
mixed_loglik(double *ZXy, longint *pdims, double *pars, longint *settings,
	     double *logLik, double *lRSS)
{				/* evaluate the log-likelihood */
  dd = dims(pdims);
  /* settings gives RML, asDelta, gradHess, and pdClass in that order */
  if (settings[ 1 ]) {		/* gradHess not used and pdClass ignored */
    *logLik = internal_loglik( dd, ZXy, pars, settings, DNULLP, lRSS);
  } else {			/* generate the Delta arrays from pars */
    setngs = settings;
    pdC = setngs + 3;
    Delta = Calloc( (dd->DmOff)[ dd->Q ], double );

    if (settings[ 2 ] == 0) {	/* no gradient or Hessian */
      *logLik =
	internal_loglik( dd, ZXy, generate_DmHalf( Delta, dd, pdC, pars ),
			 settings, DNULLP, lRSS );
    } else {
      int npar = count_DmHalf_pars( dd, pdC );
      zxdim = (dd->ZXrows) * (dd->ZXcols);
      zxcopy = Calloc( zxdim, double );
      zxcopy2 = ZXy;

      Memcpy( zxcopy, ZXy, zxdim );
      finite_diff_Hess( logLik_fun, pars, npar, logLik);
      Free( zxcopy );
    }
    Free( Delta );
  }
  dimFree( dd );
}

void				/* loglikelihood and parameter estimates */
mixed_estimate(double *ZXy, longint *pdims, double *DmHalf, longint *RML,
	       double *logLik, double *dc, longint *invert)
{				/* dc receives the decomposed ZXy array */
  dimPTR dd = dims(pdims);

  *logLik = internal_loglik(dd, ZXy, DmHalf, RML, dc, DNULLP);
  internal_estimate(dd, dc);
  if (*invert != 0) { internal_R_invert( dd, dc ); }
  dimFree(dd);
}

void				/* EM iterations for mixed-effects models */
mixed_EM(double *ZXy, longint *pdims, double *DmHalf, longint *nIter,
	 longint *pdClass, longint *RML, double *logLik, double *Ra,
	 double *lRSS)
{
  dimPTR dd = dims(pdims);
  double sigmainv, *res, *pt, 
    *dc = Calloc((size_t) ((dd->Srows) * (dd->ZXcols)), double),
    *zxcopy = Calloc((size_t) ((dd->ZXrows) * (dd->ZXcols)), double);
  double  sqrtDF = sqrt((double) (dd->N - *RML * (dd->ncol[dd->Q])));
  longint i, j, k, offset, nn = *nIter, zero = 0L;

  copy_mat(zxcopy, dd->ZXrows, ZXy, dd->ZXrows, dd->ZXrows, dd->ZXcols);
  while (nn-- > 0) {
    copy_mat(ZXy, dd->ZXrows, zxcopy, dd->ZXrows, dd->ZXrows, dd->ZXcols);
    *logLik = internal_loglik(dd, ZXy, DmHalf, RML, dc, DNULLP);
    internal_estimate( dd, dc );
    internal_R_invert( dd, dc );
    sigmainv = *(dc + (size_t)((dd->Srows) * (dd->ZXcols)) - 1)/sqrtDF;
    sigmainv = 1.0/((sigmainv < 0.0) ? - sigmainv : sigmainv);
    offset = ((dd->ZXcols) - 1L) * (dd->Srows);
    for (i = 0L; i < (dd->Q); i++) {
      longint ncol = (dd->q)[i], 
	nright = (dd->nrot)[i] - (dd->nrot)[(dd->Q) - ( (*RML) ? 0 : 1 )];
      longint nrow = (ncol + nright + 1L) * (dd->ngrp)[i];
      QRptr qq;
      pt = res = Calloc((size_t) (ncol * nrow), double);
      for (j = 0L; j < (dd->ngrp)[i]; j++) {
	copy_trans(pt, nrow, dc + (dd->SToff)[i][j], dd->Srows,
		   ncol, ncol + nright);
	pt += ncol + nright;
	scale_mat(pt++, nrow, sigmainv, dc + offset + (dd->SToff)[i][j],
		  1L, 1L, ncol);
      }
      offset -= (dd->Srows) * ncol;
      qq = QR(res, nrow, nrow, ncol);
      QRstoreR(qq, Ra + (dd->DmOff)[i], ncol);
      QRfree(qq); 
      scale_mat(res, nrow, sqrt(1.0/((dd->ngrp)[i])),
		Ra + (dd->DmOff)[i], ncol, ncol, ncol);
      switch (pdClass[i]) {
      case 0:			/* default: unstructured */
	invert_upper(res, nrow, ncol);
	copy_trans(DmHalf + (dd->DmOff)[i], ncol, res, nrow, ncol, ncol);
	break;
      case 1:			/* diagonal */
	for (j = 0; j < ncol; j++) {
	  DmHalf[ (dd->DmOff)[i] + j * (ncol + 1)] =
	    1. / sqrt( d_sum_sqr( res + j * nrow, j + 1L ) );
	}
	break;
      case 2:			/* multiple of identity */
	{
	  double aux = 0.0;
	  for(j = 0; j < ncol; j++) {
	    aux += d_sum_sqr( res + j * nrow, j + 1L );
	  }
	  aux = sqrt(ncol / aux);
	  for(j = 0; j < ncol; j++) {
	    DmHalf[(dd->DmOff)[i] + j * (ncol + 1)] = aux;
	  }
	}
      break;
      case 3:			/* compound symmetry */
	{
	  double trA = 0.0, trAJ = 0.0, *auxRes = res;
	  longint l;
	  for(j = 0; j < ncol; j++) {
	    for(k = 0; k <= j; k++) {
	      trA += auxRes[k] * auxRes[k];
	      for(l = j + 1; l < ncol; l++) {
		trAJ += auxRes[k] * auxRes[k + l * nrow];
	      }
	    }
	    auxRes += nrow;
	  }
	  trAJ = 2 * trAJ + trA;
	  trA = (ncol - 1) / (ncol * trA - trAJ);
	  trAJ = 1/trAJ - trA;
	  trA = ncol * trA + trAJ;
	  auxRes = DmHalf + (dd->DmOff[i]);
	  for(j = 0; j < ncol; j++) {
	    auxRes[j * (ncol + 1)] = trA;
	    for(k = (j + 1); k < ncol; k++) {
	      auxRes[j * ncol + k] = auxRes[j + k * ncol] = trAJ;
	    }
	  }
#ifdef R_S_H
	  F77_CALL(chol)(auxRes, &ncol, &ncol, auxRes, &l);
#else
	  F77_CALL(chol)(auxRes, &ncol, res, &zero, &zero, &l);
#endif /* R_S_H */
	}
      break;
      }
      Free(res);
    }
  }
  copy_mat(ZXy, dd->ZXrows, zxcopy, dd->ZXrows, dd->ZXrows, dd->ZXcols);
  *logLik = internal_loglik(dd, ZXy, DmHalf, RML, dc, lRSS);
  dimFree(dd); Free(dc); Free(zxcopy);
}

void				/* to be called by Fortran msmnh */
mixed_calcf(longint *n, double *theta, longint *nf,
	    double *f, longint *uiparm, double *urparm,
	    void (*ufparm)(void))
{
  Memcpy( zxcopy2, zxcopy, zxdim );
  *f = - internal_loglik(dd, zxcopy2, generate_DmHalf( Delta, dd, pdC, theta ),
			 setngs, DNULLP, DNULLP );
}

void				/* to be called by Fortran msmnh */
mixed_calcgh(longint *n, double *theta, longint *nf,
	     double *g, double *h, longint *uiparm,
	     double *urparm, void (*ufparm)(void))
{
  longint i, nn = *n;
  double *hpt = values + nn + 1;

  finite_diff_Hess( negLogLik_fun, theta, nn, values );
  Memcpy( g, values + 1, nn );
  for( i = 1; i <= nn; i++ ) {	/* copy upper triangle of Hessian */
    Memcpy( h, hpt, i );
    h += i;
    hpt += nn;
  }
}

#ifndef R_S_H

void				/* both EM and Newton-Raphson iterations */
mixed_combined(double *ZXy, longint *pdims, double *DmHalf, longint *nIter,
	 longint *pdClass, longint *RML, double *logLik, double *R0,
	 double *lRSS, longint *info)
{
  longint i, j, Qp2, ssq, p;
  int ntheta;
  longint *iv, liv, lv, uiparm[1]; /* for msmnh */
  double *theta, *work, *scale, ufparm[1], *Ra, *dc;

  dd = dims(pdims);		/* Using global dd, pdC, setngs, and Delta */
  pdC = pdClass;
  setngs = RML;
  Delta = DmHalf;

  dc = Calloc((size_t) ((dd->Srows) * (dd->ZXcols)), double);
  p = (dd->ncol)[(dd->Q)];

  ssq = 0;
  for(i = 0; i < (dd->Q); i++) {
    ssq += (dd->q)[i] * (dd->q)[i];
  }
  Ra = Calloc((size_t) ssq, double);
  Qp2 =  (dd->Q) + 2L;
  ntheta = count_DmHalf_pars( dd, pdC );
  theta = Calloc((size_t) ntheta, double);
  if (dd->ZXrows > dd->Srows) { /* Take a decomposition */
    mixed_decomp( ZXy, pdims );
    for (i = 0; i < Qp2; i++) {	/* re-write the offsets and lengths */
      for (j = 0; j < (dd->ngrp)[i]; j++) {
	(dd->ZXoff)[i][j] = (dd->DecOff)[i][j];
	(dd->ZXlen)[i][j] = (dd->DecLen)[i][j];
      }
    }
    pdims[1] = pdims[4];	/* ZXrows = Srows in original locations */
    dd->ZXrows = dd->Srows;	/* and in the copy */
  }
  zxdim = (dd->ZXrows) * (dd->ZXcols); /* also global zxdim, zxcopy, and zxcopy2 */
  zxcopy = Calloc( zxdim, double );
  zxcopy2 = ZXy;
  Memcpy( zxcopy, ZXy, zxdim );	/* keep a copy before we mess it up */

  *uiparm = 0; work = DmHalf;	/* check for non-zero entries in DmHalf */
  for( i = 0; i < dd->Q; i++ ) {
    for( j = 0; j < (dd->q)[i] * (dd->q)[i]; j++ ) {
      if (*work++ != 0.) {
	*uiparm = 1;
      }
    }
  }
  if (*uiparm != 0) {
    *uiparm = 0;
  } else {			/* initialize the DmHalf arrays */
    work = ZXy;
    scale = DmHalf;
    for( i = 0; i < dd->Q; i++ ) {
      for ( j = 0; j < (dd->q)[i]; j++ ) {
	*scale = 0.375 * sqrt( d_dot_prod( work, 1L, work, 1L, dd->ZXrows ) /
			       (dd->ngrp)[i]);
	scale += (dd->q)[i] + 1;
	work += dd->ZXrows;
      }
      scale -= (dd->q)[i];	/* have moved too far - step back */
    }
  }
  mixed_EM(ZXy, pdims, DmHalf, nIter, pdClass, RML, logLik, Ra, lRSS);
  generate_theta( theta, dd, pdClass, DmHalf );

  values = Calloc( (size_t) ntheta * (ntheta + 1) + 1, double ); /* global */
  liv = 60;
  iv = Calloc( (size_t) liv, longint );
  lv = 78 + ntheta * (ntheta + 12);
  work = Calloc( (size_t) lv, double );
  scale = Calloc( (size_t) ntheta, double );
  for( i = 0; i < ntheta; i++ ) { scale[i] = 1.; }
  F77_CALL(msmnh) (&ntheta, scale, theta, mixed_calcf, mixed_calcgh,
		   iv, &liv, &lv, work, uiparm, ufparm, abort);
  *info = iv[0];
  Memcpy( zxcopy2, zxcopy, zxdim );
  *logLik = internal_loglik( dd, zxcopy2, generate_DmHalf( Delta, dd, pdC, theta ),
			     setngs, dc, lRSS );
  copy_mat(R0, p, dc + (dd->SToff)[(dd->Q)][0], (dd->Srows), p, p + 1);
  Free( values ); Free(iv); Free(work); Free(scale); Free( zxcopy );
  dimFree( dd ); Free( theta ); Free(dc); Free(Ra);
}

#endif /* R_S_H */

/* functions for calculating df's for fixed effects tests */

static double
inner_perc(double *x, longint *grp, longint n)
     /* percentage of groups for which x is inner */
{
  /* x - column of X matrix to be assessed
     grp - integer vector with groups
     n - length of x and grp
     data are assumed to be ordered by grp */

  longint currGrp, nn = 0, isInner;
  double nInner = 0., nGrp = 0., currVal;

  while (nn < n) {
    currGrp = grp[nn];
    currVal = x[nn];
    nGrp++;
    isInner = 0;
    do {
      if (isInner == 0 && x[nn] != currVal) {
	nInner++;
	isInner = 1;
      }
      nn++;
    } while (nn < n && currGrp == grp[nn]);
  }
  return(nInner/nGrp);
}

void
inner_perc_table(double *X, longint *grps, longint *p, longint *Q,
		 longint *n, double *pTable) 
     /* constructs an p x Q "inner-percentage" table for a fixed effects
	matrix X and a set of grouping vectors grps */
{
  longint i, j, pp = *p, nn = *n, ipp = 0, inn = 0;
  for(i = 0; i < *Q; i++) {
    for(j = 0; j < pp; j++) {
      pTable[j + ipp] = inner_perc(X + j * nn, grps + inn, nn);
    }
    ipp += pp;
    inn += nn;
  }
}

/* gls functions */
void
gls_loglik(double *Xy, longint *pdims, double *logLik, double *lRSS)
{
  longint i, N = pdims[0], p = pdims[1], RML = pdims[2],
    Np1 = N + 1, Nr = N - RML * p, rnkm1;
  QRptr dmQR;

  dmQR = QR(Xy, N, N, p + 1);
  rnkm1 = dmQR->rank - 1;
  *lRSS = log(fabs(dmQR->mat[rnkm1 * Np1]));
  *logLik -= Nr * (*lRSS);
  if (RML == 1) {
    for(i = 0; i < rnkm1; i++) {
      *logLik -= log(fabs(dmQR->mat[i * Np1]));
    }
  }
  QRfree(dmQR);
}
  
void 
gls_estimate(double *Xy, longint *pdims, double *beta, double *sigma, 
	     double *logLik, double *varBeta, longint *rank, longint *pivot)
{
  longint i, N = pdims[0], p = pdims[1], RML = pdims[2], pp1 = p + 1,
    Nr = N - RML * p, rk, rkm1, rkp1;
  QRptr dmQR;
  double *R = Calloc((size_t) (pp1 * pp1), double);

  dmQR = QR(Xy, N, N, pp1);
  *rank = rk = dmQR->rank;
  rkm1 = rk - 1;
  rkp1 = rk + 1;
  Memcpy(pivot, dmQR->pivot, pp1);
  for(i = 0; i < rk; i++) {
    Memcpy(R + i * rk, dmQR->mat + i * N, i + 1);
  }
  *sigma = fabs(R[rk * rk - 1]);
  *logLik -= Nr * log(*sigma);
  *sigma /= sqrt(((double) Nr));
  if (RML == 1) {
    for(i = 0; i < rkm1; i++) {
      *logLik -= log(fabs(R[i * (rkp1)]));
    }
  }
  copy_mat(varBeta, rkm1, R, rk, rkm1, rkm1);
  invert_upper(varBeta, rkm1, rkm1);
  mult_mat(beta, rkm1, varBeta, rkm1, rkm1, rkm1, R + rkm1 * rk, rk,  1L);
  QRfree(dmQR);
  Free(R);
}

/* Positive definite matrices */

static void 
Chol_pd(double *L, longint *q, double *l)
{
  longint i, qq = *q;
  for(i = 0; i < qq; i++) {
    Memcpy(L + i * qq, l, i + 1);
    l += i + 1;
  }
}

void 
matrixLog_pd(double *L, longint *q, double *l)
{
  longint i, j, qq = *q, one = 1L, info = 0L;
  if ( qq == 1 ) {
    *L = exp( *l );
  } else {
    double *vectors = Calloc((size_t) qq * qq, double), 
      *work1 = Calloc((size_t) qq, double), *work2 = Calloc((size_t) qq, double), 
      *values = Calloc((size_t) qq, double);
    Chol_pd(L, q, l);
    for(i = 0; i < qq - 1; i++) {
      copy_mat(L + (i * (qq + 1) + 1), 1L, L + i * (qq + 1) + qq, qq, 1L,
	       qq - (i + 1));
    }
    F77_CALL(rs) (q, q, L, values, &one, vectors, work1, work2, &info);
    for(i = 0; i < qq; i++) {
      values[i] = exp(values[i]);
      for(j = 0; j < qq; j++) {
	vectors[i * qq + j] *= values[i];
      }
    }
    copy_trans(L, qq, vectors, qq, qq, qq);
    Free(vectors); Free(work1); Free(work2); Free(values);
  }
}


void
natural_pd(double *L, longint *q, double *l) /* natural parametrization  */
{
  longint i, j, qp1 = *q + 1, info, zero = 0L;
  double *std = l, *corr = l + *q, *work = Calloc(*q, double);

  for(i = 0; i < *q; i++) std[i] = exp(std[i]);

  for(i = 0; i < *q; i++) {
    L[i * qp1] = std[i] * std[i];
    for(j = i + 1; j < *q; j++) {
      *corr = exp(*corr);
      *corr = (*corr - 1)/(*corr + 1);
      L[i * (*q) + j] = L[j * (*q) + i] = std[i] * std[j] * (*corr);
      corr++;
    }
  }
#ifdef R_S_H
  F77_CALL(chol) (L, q, q, L, &info);
#else
  F77_CALL(chol) (L, q, work, &zero, &zero, &info);
#endif /* R_S_H */
  Free(work);
}

void
compSymm_pd(double *L, longint *q, double *l) /* compound symmetry */
{
  longint i, j, qp1 = *q + 1;
  double aux = exp(l[0]), aux1 = exp(l[1]), aux2;

  aux1 = (aux1 - 1.0/((double) *q - 1.0))/(aux1 + 1.0);
  aux2 = aux * sqrt(1.0 - aux1);
  aux1 = aux * sqrt((1.0 + (*q - 1.0) * aux1) / ((double) *q));

  for(i = 0; i < *q; i++) {
    L[i * (*q)] = aux1;
  }
  for(i = 1; i < *q; i++) {
    aux = -aux2/sqrt(i * (i + 1));
    for(j = 0; j < i; j++) {
      L[i + (j * (*q))] = aux;
    }
    L[i * qp1] = -aux * i;
  }
}

/* Forming the parameter structure from the Delta matrix */
/*  Not sure if these will ever be called from S. */
/*  Will leave open the possibility. */

static void
Delta2MatrixLog( double *theta, longint *q, double *Delta )
{
  longint i, j, qq = *q, one = 1L, info = 0L;
  if ( qq == 1 ) {
      *theta = log(*Delta * *Delta)/2.;
  } else {
    double *vectors = Calloc((size_t) qq * qq, double),
      *DtransD = Calloc((size_t) qq * qq, double),
      *workmat = Calloc((size_t) qq * qq, double),
      *work2 = Calloc((size_t) qq, double), 
      *values = Calloc((size_t) qq, double), *pt;
    crossprod_mat(DtransD, qq, Delta, qq, qq, qq); /* form t(Delta) %*% Delta */
    F77_CALL(rs) (q, q, DtransD, values, &one, vectors, workmat, work2, &info);
    if (info != 0L) {
      PROBLEM "Unable to form eigenvalue-eigenvector decomposition"
	RECOVER(NULL_ENTRY);
    }
    copy_mat(workmat, qq, vectors, qq, qq, qq);
    for(i = 0; i < qq; i++) {
      values[i] = log(values[i])/2;
      for(j = 0; j < qq; j++) {
	workmat[i * qq + j] *= values[i];
      }
    }
    copy_trans(DtransD, qq, workmat, qq, qq, qq);
    mult_mat(workmat, qq, vectors, qq, qq, qq, DtransD, qq, qq);
    for( i = 0, pt = theta; i < qq; i++ ) {
      for( j = 0; j <= i; j++ ) {
	*pt++ = workmat[ i * qq + j ];
      }
    }
    Free(vectors); Free(DtransD); Free(workmat), Free(work2); Free(values);
  }
}

  
/* Correlation Structures */

/* Factor list and Recalc for general corStruct object */

void 
corStruct_factList(double *mat, longint *pdims, double *FactorL, double *logdet)
{
  longint i, j, M = pdims[1], *len = pdims + 4, job = 11L, info, zero = 0L;
  double *work, *work1;
  
  for(i = 0; i < M; i++) {
    longint li = len[i], lisq = li * li, lip1 = li + 1;
    work = Calloc(li, double);
    work1 = Calloc(lisq, double);
    F77_CALL(chol) (mat, &li, work, &zero, &zero, &info);
    for(j = 0; j < li; j++) {
      work1[j * lip1] = 1;
      F77_CALL(dtrsl) (mat, &li, &li, work1 + j * li, &job, &info);
      *logdet -= log(fabs(mat[j * lip1]));
    }
    Memcpy(FactorL, work1, lisq);
    Free(work); Free(work1);
    FactorL += lisq;
    mat += lisq;
  }
}

void
corStruct_recalc(double *Xy, longint *pdims, longint *ZXcol, double *Factor)
{
  longint N = pdims[0], M = pdims[1], *len = pdims + 4, *start = len + M, i;
  for(i = 0; i < M;  i++) {
    mult_mat(Xy + start[i], N, Factor, len[i], len[i], len[i], 
	     Xy + start[i], N, *ZXcol);
    Factor += (len[i] * len[i]);
  }
}

/* symm class - unstructured correlation - based on spherical
   parametrization */ 

void
symm_fullCorr(double *par, longint *maxC, double *crr)
 /* calculates full correlation structure  */
{
  double *work, aux, aux1, *src = par, *src1, *dest;
  longint i, j, n = *maxC;

  /* first get upper-triangular factor */
  dest = work = Calloc(n * (n + 1) / 2, double);
  for(i = 0; i < n; i++) {
    aux = 1.0;
    for(j = 0; j < i; j++) {
      aux1 = exp(*src);
      aux1 = PI * aux1/(1 + aux1); /* untransforming */
      *dest = aux * cos(aux1);
      aux *= sin(aux1);
      dest++; src++;
    }
    *dest = aux;
    dest++;
  }

  /* getting the correlations */
  for(i = 0, dest = crr, src = work; i < n - 1; i++) {
    longint ip1 = i + 1;
    src += i;
    for(j = ip1, src1 = src; j < n; j++) {
      src1 += j;
      *dest = d_dot_prod(src, 1L, src1, 1L, ip1);
      dest++;
    }
  }
  Free(work);
}

static void 
symm_mat(double *crr, longint *time, longint *n, longint *maxC, double *mat)
{
  longint i, j, k, np1 = *n + 1, n1, n2;
  for(i = 0; i < *n; i++) {
    mat[i * np1] = 1.0;
    for(j = i + 1; j < *n; j++) {
      n1 = (time[i] < time[j]) ? time[i] : time[j];
      n2 = time[i] + time[j] - 2 * n1 - 1;
      k = n1 * *maxC - n1 * (n1 + 1) / 2 + n2;
      *(mat + i + j * (*n)) = *(mat + j + i * (*n)) = crr[k];
    }
  }
}

void 
symm_matList(double *pars, longint *time, longint *maxC,
	     longint *pdims, double *mat)
{
  double *crr = Calloc(*maxC * (*maxC - 1) / 2, double);
  longint i, M = pdims[1], *len = pdims + 4;
  /* parameters assumed in unconstrained form */
  symm_fullCorr(pars, maxC, crr);
  for(i = 0; i < M;  i++) {
    symm_mat(crr, time, &len[i], maxC, mat);
    mat += len[i] * len[i];
    time += len[i];
  }
  Free(crr);
}

static void 
symm_fact(double *crr, longint *time, longint *n, longint *maxC, double *mat, 
	  double *logdet)
{
  longint job = 11L, info, zero = 0L, i, nsq = *n * (*n), np1 = *n + 1;
  double *work = Calloc(*n, double), *work1 = Calloc(nsq, double);

  symm_mat(crr, time, n, maxC, mat);
  F77_CALL(chol) (mat, n, work, &zero, &zero, &info);
  for(i = 0; i < *n; i++) {
    work1[i * np1] = 1;
    F77_CALL(dtrsl) (mat, n, n, work1 + i * (*n), &job, &info);
    *logdet -= log(fabs(mat[i * np1]));
  }
  Memcpy(mat, work1, nsq);
  Free(work); Free(work1);
}

void 
symm_factList(double *pars, longint *time, longint *maxC, longint *pdims, 
	      double *FactorL, double *logdet)
{
  double *crr = Calloc(*maxC * (*maxC - 1L) / 2L, double);
  longint i, M = pdims[1], *len = pdims + 4;
  /* parameters assumed in unconstrained form */
  symm_fullCorr(pars, maxC, crr);
  for(i = 0; i < M;  i++) {
    symm_fact(crr, time, &len[i], maxC, FactorL, logdet);
    FactorL += len[i] * len[i];
    time += len[i];
  }
  Free(crr);
}
  
void
symm_recalc(double *Xy, longint *pdims, longint *ZXcol, double *pars, 
	    longint *time, longint *maxC, double *logdet)
{
  longint N = pdims[0], M = pdims[1], *len = pdims + 4, *start = len + M, i;
  double *crr = Calloc(*maxC * (*maxC - 1) / 2, double);
  /* parameters assumed in unconstrained form */
  symm_fullCorr(pars, maxC, crr);
  for(i = 0; i < M;  i++) {
    double *Factor = Calloc((len[i] * len[i]), double);
    symm_fact(crr, time + start[i], &len[i], maxC, Factor, logdet);
    mult_mat(Xy + start[i], N, Factor, len[i], len[i], len[i], 
	     Xy + start[i], N, *ZXcol);
    Free(Factor); 
  }
  Free(crr);
}

/* nat class - unstructured correlation - natural parametrization */ 

void
nat_fullCorr(double *par, longint *maxC, double *crr)
 /* calculates full correlation structure  */
{
  double aux;
  longint i, npar = *maxC * (*maxC - 1) / 2;

  for(i = 0; i < npar; i++) {
    aux = exp(par[i]);
    crr[i] = (aux - 1)/(aux + 1);
  }
}

void 
nat_matList(double *pars, longint *time, longint *maxC,
	     longint *pdims, double *mat)
{
  double *crr = Calloc(*maxC * (*maxC - 1) / 2, double);
  longint i, M = pdims[1], *len = pdims + 4;
  /* parameters assumed in unconstrained form */
  nat_fullCorr(pars, maxC, crr);
  for(i = 0; i < M;  i++) {
    symm_mat(crr, time, &len[i], maxC, mat);
    mat += len[i] * len[i];
    time += len[i];
  }
  Free(crr);
}

void 
nat_factList(double *pars, longint *time, longint *maxC, longint *pdims, 
	      double *FactorL, double *logdet)
{
  double *crr = Calloc(*maxC * (*maxC - 1L) / 2L, double);
  longint i, M = pdims[1], *len = pdims + 4;
  /* parameters assumed in unconstrained form */
  nat_fullCorr(pars, maxC, crr);
  for(i = 0; i < M;  i++) {
    symm_fact(crr, time, &len[i], maxC, FactorL, logdet);
    FactorL += len[i] * len[i];
    time += len[i];
  }
  Free(crr);
}
  
void
nat_recalc(double *Xy, longint *pdims, longint *ZXcol, double *pars, 
	    longint *time, longint *maxC, double *logdet)
{
  longint N = pdims[0], M = pdims[1], *len = pdims + 4, *start = len + M, i;
  double *crr = Calloc(*maxC * (*maxC - 1) / 2, double);
  /* parameters assumed in unconstrained form */
  nat_fullCorr(pars, maxC, crr);
  for(i = 0; i < M;  i++) {
    double *Factor = Calloc((len[i] * len[i]), double);
    symm_fact(crr, time + start[i], &len[i], maxC, Factor, logdet);
    mult_mat(Xy + start[i], N, Factor, len[i], len[i], len[i], 
	     Xy + start[i], N, *ZXcol);
    Free(Factor); 
  }
  Free(crr);
}

/* AR1 class */

static void
AR1_mat(double *par, longint *n, double *mat)
{
  longint i, j;
  double aux;
  for(i = 0; i < *n; i++) {
    for(j = i; j < *n; j++) {
      aux = pow(*par, j - i);
      *(mat + i + j * (*n)) = *(mat + j + i * (*n)) = aux;
    }
  }
}
  
void 
AR1_matList(double *par, longint *pdims, double *mat)
{
  longint i, M = pdims[1], *len = pdims + 4;
  /* par assumed in unconstrained form */
  double aux = exp(*par);
  *par = (aux - 1.0) / (aux + 1.0);
  for(i = 0; i < M;  i++) {
    AR1_mat(par, &len[i], mat);
    mat += len[i] * len[i];
  }
}

static void 
AR1_fact(double *par, longint *n, double *mat, double *logdet)
{
  longint i, np1 = *n + 1;
  double aux = sqrt(1 - *par * (*par)), aux1 = - (*par)/aux;

  *logdet -= (*n - 1) * log(aux);
  aux = 1/aux;
  mat[0] = 1;
  for(i = 1; i < *n; i++) {
    mat[i * np1] = aux;
    mat[i + *n * (i - 1)] = aux1;
  }
}

void 
AR1_factList(double *par, longint *pdims, double *FactorL, 
	     double *logdet)
{
  longint i, M = pdims[1], *len = pdims + 4;
  /* par assumed in unconstrained form */
  double aux = exp(*par);
  *par = (aux - 1.0) / (aux + 1.0);
  for(i = 0; i < M;  i++) {
    AR1_fact(par, &len[i], FactorL, logdet);
    FactorL += len[i] * len[i];
  }
}

void
AR1_recalc(double *Xy, longint *pdims, longint *ZXcol, double *par, 
	   double *logdet) 
{
  longint N = pdims[0], M = pdims[1], *len = pdims + 4,  *start = len + M, i;
  double *Factor, aux = exp(*par);
  /* par assumed in unconstrained form */
  *par = (aux - 1.0) / (aux + 1.0);
  for(i = 0; i < M;  i++) {
    Factor = Calloc(len[i] * len[i], double);
    AR1_fact(par, &len[i], Factor, logdet);
    mult_mat(Xy + start[i], N, Factor, len[i], len[i], len[i], 
	     Xy + start[i], N, *ZXcol);
    Free(Factor);
  }
}


/* Continuous AR1 class */

static void
CAR1_mat(double *par, double *time, longint *n, double *mat)
{
  longint i, j;
  double aux;
  for(i = 0; i < *n; i++) {
    for(j = i; j < *n; j++) {
      aux = pow(*par, fabs(time[j] - time[i]));
      *(mat + i + j * (*n)) = *(mat + j + i * (*n)) = aux;
    }
  }
}
  
void 
CAR1_matList(double *par, double *time, longint *pdims, double *mat)
{
  longint i, M = pdims[1], *len = pdims + 4;
  double aux = exp(*par);
  /* parameter assumed in unconstrained form */
  *par = aux / (1.0 + aux);
  for(i = 0; i < M;  i++) {
    CAR1_mat(par, time, &len[i], mat);
    mat += len[i] * len[i];
    time += len[i];
  }
}

static void 
CAR1_fact(double *par, double *time, longint *n, double *mat, double *logdet)
{
  longint job = 11L, info, zero = 0L, i, nsq = *n * (*n), np1 = *n + 1;
  double *work = Calloc(*n, double), *work1 = Calloc(nsq, double);
  CAR1_mat(par, time, n, mat);
  F77_CALL(chol) (mat, n, work, &zero, &zero, &info);
  for(i = 0; i < *n; i++) {
    work1[i * np1] = 1;
    F77_CALL(dtrsl) (mat, n, n, work1 + i * (*n), &job, &info);
    *logdet -= log(fabs(mat[i * np1]));
  }
  Memcpy(mat, work1, nsq);
  Free(work); Free(work1);
}

void 
CAR1_factList(double *par, double *time, longint *pdims,  
	     double *FactorL, double *logdet)
{
  longint i, M = pdims[1], *len = pdims + 4;
  double aux = exp(*par);
  /* parameter assumed in unconstrained form */
  *par = aux / (1.0 + aux);
  for(i = 0; i < M;  i++) {
    CAR1_fact(par, time, &len[i], FactorL, logdet);
    FactorL += len[i] * len[i];
    time += len[i];
  }
}

void
CAR1_recalc(double *Xy, longint *pdims, longint *ZXcol, 
	   double *par, double *time, double *logdet)
{
  longint N = pdims[0], M = pdims[1], *len = pdims + 4, *start = len + M, i;
  double aux = exp(*par);
  /* parameter assumed in unconstrained form */
  *par = aux / (1.0 + aux);
  for(i = 0; i < M;  i++) {
    double *Factor = Calloc(len[i] * len[i], double);
    CAR1_fact(par, time + start[i], &len[i], Factor, logdet);
    mult_mat(Xy + start[i], N, Factor, len[i], len[i], len[i], 
	     Xy + start[i], N, *ZXcol);
    Free(Factor);
  }
}

/* ARMA class */

static void 
ARMA_transPar(longint N, double *pars, double sgn)
{
  longint i, j, n, n2;
  double ps, D, aux;
  for(n = N - 1; n > -1; n--) {
    if ((ps = pars[n] * pars[n]) >= 1.0) {
      PROBLEM "All parameters must be less than 1 in absolute value"
	RECOVER(NULL_ENTRY);
    }
    if (n) {
      D = 1 - ps;
      n2 = (n - 1)/2;
      for(i = 0; i <= n2; i++) {
	if ((j = n - i -1) > i) {
	  aux = (pars[i] + sgn * pars[j] * pars[n])/D;
	  pars[j] = (pars[j] + sgn * pars[i] * pars[n])/D;
	  pars[i] = aux;
	} else {
	  pars[i] /= (1 - sgn * pars[n]);
	}
      }
    }
    pars[n] = log((1 + pars[n])/(1 - pars[n]));
  }
}

void 
ARMA_unconstCoef(longint *p, longint *q, double *pars)
{
  ARMA_transPar(*p, pars, 1.0);
  ARMA_transPar(*q, pars + *p, -1.0);
}

static void 
ARMA_untransPar(longint N, double *pars, double sgn)
{
  longint i, j;
  double *aux;
  if (N) {
    aux = Calloc(N, double);
    for(i = 0; i < N; i++) {
      pars[i] = exp(-pars[i]);
      aux[i] = pars[i] = (1 - pars[i])/(1 + pars[i]);
      if (i) {
	for(j = 0; j < i; j++) {
	  pars[j] = aux[j] + sgn * pars[i] * aux[i - j - 1];
	}
	Memcpy(aux, pars, i);
      }
    }
    Free(aux);
  }
}

void 
ARMA_constCoef(longint *p, longint *q, double *pars)
{
  ARMA_untransPar(*p, pars, -1.0);
  ARMA_untransPar(*q, pars + *p, 1.0);
}

static void 
ARMA_cross(longint *p, longint *q, double *pars, double *psi)
{
  longint i, j, M = *q + 1, PM;
  M = (*p > M ? *p : M);
  psi[0] = 1;
  for(i = 1; i < M; i++) {
    psi[i] = ((*q < i) ? 0 : pars[*p + i - 1]);
    PM = (*p < i ? *p : i);
    for(j = 0; j < PM; j++) {
      psi[i] += pars[j] * psi[i - j - 1];
    }
  }
}

static void 
ARMA_corr(longint *p, longint *q, longint *maxlag, double *pars, double *psi, 
	  double *crr) 
{
  longint P = *p + 1, Pp1 = P + 1, i, j, k, minPQ, Mlag, maxPQ,
       *pivot = Calloc(P, longint);
  double *coef = Calloc(P * P, double), *src, *qraux = Calloc(P, double),
         *work = Calloc(P * P, double), *work1;
  
  if (!sqrt_eps) sqrt_eps = sqrt(DOUBLE_EPS);
  if ((maxPQ = ((*p > *q) ? *p : *q))) {
    for(i = 0, src = coef; i < P; i++, src += Pp1) {
      crr[i] = 0;
      *src = 1;
    }
    Mlag = ((*maxlag > *q) ? *maxlag : *q);
    Mlag = ((Mlag > *p) ? Mlag : *p) + 1;
    work1 = Calloc(Mlag, double);
    for(i = P; i < Mlag; i++) {
      crr[i] = 0;
    }
    crr[0] = 1;
    for(i = 1, src = pars + *p; i <= *q; i++, src++) {
      crr[0] += (*src) * psi[i];
    }
    if (*p) {
      if ((minPQ = ((*p < *q) ? *p : *q))) {
	for(i = 1, src = pars + *p - 1; i <= minPQ; i++) {
	  for(j = i; j <= *q; j++) {
	    crr[i] += *(src + j) * psi[j - i];
	  }
	}
      }
      for(i = 0, src = coef; i < P; i++, src++) {
	for(j = 0; j < *p; j++) {
	  k = i - j - 1;
	  k = ((k > 0) ? k : -k);
	  *(src + (k * P)) -= pars[j];
	}
      }
      F77_CALL(dqrdca) (coef, &P, &P, &P, qraux, pivot, work, &i, &sqrt_eps);
      if (i < P) {
	PROBLEM "Coeficient matrix not invertible" RECOVER(NULL_ENTRY);
      }
      i = 100L;
      F77_CALL(dqrsl) (coef, &P, &P, &P, qraux, crr, DNULLP, crr, work1, DNULLP, 
		       DNULLP, &i, &j);
      Memcpy(crr, work1, Mlag);
    }
    for(i = P; i <= *q; i++) {
      for(j = 0; j < *p; j++) {
	crr[i] += pars[j] * crr[i - j - 1];
      }
      for(j = i, src = pars + i - 1; j <= *q; j++, src++) {
	crr[i] += *src * psi[j - i];
      }
    }
    for(i = maxPQ + 1; i < Mlag; i++) {
      for(j = 0; j < *p; j++) {
	crr[i] += pars[j] * crr[i - j - 1];
      }
    }
    for(i = 1; i < Mlag; i++) {
      crr[i] /= crr[0];
    }
    Free(qraux); Free(work); Free(coef); Free(pivot); Free(work1);
  }
  crr[0] = 1;
}

static void 
ARMA_fullCorr(longint *p, longint *q, longint *maxlag, double *pars,
	      double *crr)
{
  longint M = *q + 1;
  double *psi;
  M = ((M < *p) ? *p : M);
  psi = Calloc(M, double);
  ARMA_cross(p, q, pars, psi);
  ARMA_corr(p, q, maxlag, pars, psi, crr);
  Free(psi);
}

static void 
ARMA_mat(double *crr, longint *time, longint *n, double *mat)
{
  longint i, j, k;
  for(i = 0; i < *n; i++) {
    for(j = i; j < *n; j++) {
      k = time[j] - time[i];
      k = ((k < 0) ? -k : k);
      *(mat + i + j * (*n)) = *(mat + j + i * (*n)) = crr[k];
    }
  }
}

void 
ARMA_matList(double *pars, longint *p, longint *q, longint *time,
	     longint *maxlag, longint *pdims, double *mat)
{
  double *crr = Calloc(*maxlag + 1L, double);
  longint i, M = pdims[1], *len = pdims + 4;
  /* parameters assumed in unconstrained form */
  ARMA_constCoef(p, q, pars);
  ARMA_fullCorr(p, q, maxlag, pars, crr);
  for(i = 0; i < M;  i++) {
    ARMA_mat(crr, time, &len[i], mat);
    mat += len[i] * len[i];
    time += len[i];
  }
  Free(crr);
}

static void 
ARMA_fact(double *crr, longint *time, longint *n, double *mat, double *logdet)
{
  longint job = 11L, info, zero = 0L, i, nsq = *n * (*n), np1 = *n + 1;
  double *work = Calloc(*n, double), *work1 = Calloc(nsq, double);
  ARMA_mat(crr, time, n, mat);
  F77_CALL(chol) (mat, n, work, &zero, &zero, &info);
  for(i = 0; i < *n; i++) {
    work1[i * np1] = 1;
    F77_CALL(dtrsl) (mat, n, n, work1 + i * (*n), &job, &info);
    *logdet -= log(fabs(mat[i * np1]));
  }
  Memcpy(mat, work1, nsq);
  Free(work); Free(work1);
}

void 
ARMA_factList(double *pars, longint *p, longint *q, longint *time,
	      longint *maxlag, longint *pdims, double *FactorL,
	      double *logdet)
{
  double *crr = Calloc(*maxlag + 1L, double);
  longint i, M = pdims[1], *len = pdims + 4;
  /* parameters assumed in unconstrained form */
  ARMA_constCoef(p, q, pars);
  ARMA_fullCorr(p, q, maxlag, pars, crr);
  for(i = 0; i < M;  i++) {
    ARMA_fact(crr, time, &len[i], FactorL, logdet);
    FactorL += len[i] * len[i];
    time += len[i];
  }
  Free(crr);
}
  
void
ARMA_recalc(double *Xy, longint *pdims, longint *ZXcol, double *pars, 
	   longint *p, longint *q, longint *time, longint *maxlag,
	    double *logdet)
{
  longint N = pdims[0], M = pdims[1], *len = pdims + 4, *start = len + M, i;
  double *crr = Calloc(*maxlag + 1L, double);
  /* parameters assumed in unconstrained form */
  ARMA_constCoef(p, q, pars);
  ARMA_fullCorr(p, q, maxlag, pars, crr);
  for(i = 0; i < M;  i++) {
    double *Factor = Calloc(len[i] * len[i], double);
    ARMA_fact(crr, time + start[i], &len[i], Factor, logdet);
    mult_mat(Xy + start[i], N, Factor, len[i], len[i], len[i], 
	     Xy + start[i], N, *ZXcol);
    Free(Factor); 
  }
}

/* Compound symmetry */

static void
compSymm_mat(double *par, longint *n, double *mat)
{
  longint i, j;
  for(i = 0; i < *n; i++) {
    mat[(*n + 1) * i] = 1.0;
    for(j = i + 1; j < *n; j++) {
      *(mat + i + j * (*n)) = *(mat + j + i * (*n)) = *par;
    }
  }
}
  
void 
compSymm_matList(double *par, double *inf, longint *pdims, double *mat)
{
  longint i, M = pdims[1], *len = pdims + 4;
  /* parameter assumed in unconstrained form */
  double aux = exp(*par);
  *par = (aux + *inf)/(aux + 1.0);
  for(i = 0; i < M;  i++) {
    compSymm_mat(par, &len[i], mat);
    mat += len[i] * len[i];
  }
}

static void 
compSymm_fact(double *par, longint *n, double *mat, double *logdet)
{
  longint i, j, np1 = *n + 1, nsq = *n * (*n);
  double aux, aux1, *work = Calloc(nsq, double);
  aux = 1 + (*n - 1) * (*par);
  *logdet -= log(aux)/2;
  aux = 1/sqrt(aux * (*n));
  for(i = 0; i < nsq; i += *n) {
    work[i] = aux;
  }
  aux = 1 - (*par);
  *logdet -= (*n - 1) * log(aux)/2;
  for(i = 1; i < *n; i++) {
    aux1 = -1/sqrt(aux * i * (i + 1));
    for(j = 0; j < i; j++) {
      work[i + j * (*n)] = aux1;
    }
    work[i * np1] = -i * aux1;
  }
  Memcpy(mat, work, nsq);
  Free(work);
}

void 
compSymm_factList(double *par, double *inf, longint *pdims, 
		 double *FactorL, double *logdet)
{
  longint i, M = pdims[1], *len = pdims + 4;
  /* parameter assumed in unconstrained form */
  double aux = exp(*par);
  *par = (aux + *inf)/(aux + 1.0);
  for(i = 0; i < M;  i++) {
    compSymm_fact(par, &len[i], FactorL, logdet);
    FactorL += len[i] * len[i];
  }
}

void
compSymm_recalc(double *Xy, longint *pdims, longint *ZXcol, double *par,  
		double *inf, double *logdet)
{
  longint N = pdims[0], M = pdims[1], *len = pdims + 4, *start = len + M, i;
  double aux = exp(*par);
  /* parameter assumed in unconstrained form */
  *par = (aux + *inf)/(aux + 1.0);
  for(i = 0; i < M;  i++) {
    double *Factor = Calloc(len[i] * len[i], double);
    compSymm_fact(par, &len[i], Factor, logdet);
    mult_mat(Xy + start[i], N, Factor, len[i], len[i], len[i], 
	     Xy + start[i], N, *ZXcol);
    Free(Factor);
  }
}

/* Huyn-Feldt class */

static void
HF_mat(double *par, longint *time, longint *n, double *mat)
{
  longint i, j, np1 = *n + 1;
  for(i = 0; i < *n; i++) {
    mat[i * np1] = par[time[i]];
    for(j = i + 1; j < *n; j++) {
      *(mat + i + j * (*n)) = *(mat + j + i * (*n)) = 
	0.5 * (par[time[i]] + par[time[j]]) - 1.0;
    }
  }
}
  
void 
HF_matList(double *par, longint *maxC, longint *time, longint *pdims,
	     double *mat)
{
  longint i, M = pdims[1], *len = pdims + 4;
  /* parameter assumed in unconstrained form */
  double inf = -1.0/(2.0 * ((double) *maxC));
  for(i = 0; i < *maxC; i++) {
    par[i] = 2.0 * (exp(par[i]) + inf) + 1.0;
  }
  for(i = 0; i < M;  i++) {
    HF_mat(par, time, &len[i], mat);
    mat += len[i] * len[i];
    time += len[i];
  }
}

static void 
HF_fact(double *par, longint *time, longint *n, double *mat, double *logdet)
{
  longint job = 11L, info, zero = 0L, i, nsq = *n * (*n), np1 = *n + 1;
  double *work = Calloc(*n, double), *work1 = Calloc(nsq, double);
  HF_mat(par, time, n, mat);
  F77_CALL(chol) (mat, n, work, &zero, &zero, &info);
  for(i = 0; i < *n; i++) {
    work1[i * np1] = 1;
    F77_CALL(dtrsl) (mat, n, n, work1 + i * (*n), &job, &info);
    *logdet -= log(fabs(mat[i * np1]));
  }
  Memcpy(mat, work1, nsq);
  Free(work); Free(work1);
}

void 
HF_factList(double *par, longint *maxC, longint *time, longint *pdims,
	     double *FactorL, double *logdet)
{
  longint i, M = pdims[1], *len = pdims + 4;
  /* parameter assumed in unconstrained form */
  double inf = -1.0/(2.0 * ((double) *maxC));
  for(i = 0; i < *maxC; i++) {
    par[i] = 2.0 * (exp(par[i]) + inf) + 1.0;
  }
  for(i = 0; i < M;  i++) {
    HF_fact(par, time, &len[i], FactorL, logdet);
    FactorL += len[i] * len[i];
    time += len[i];
  }
}

void
HF_recalc(double *Xy, longint *pdims, longint *ZXcol, double *par, 
	 longint *time, longint *maxC, double *logdet)
{
  longint N = pdims[0], M = pdims[1], *len = pdims + 4,  *start = len + M, i;
  double inf = -1.0/(2.0 * ((double) *maxC));
  /* parameter assumed in unconstrained form */
  for(i = 0; i < *maxC; i++) {
    par[i] = 2.0 * (exp(par[i]) + inf) + 1.0;
  }
  for(i = 0; i < M;  i++) {
    double *Factor = Calloc(len[i] * len[i], double);
    HF_fact(par, time + start[i], &len[i], Factor, logdet);
    mult_mat(Xy + start[i], N, Factor, len[i], len[i], len[i], 
	     Xy + start[i], N, *ZXcol);
    Free(Factor);
  }
}

/* Spatial correlation structures */

/* Spherical class */

static double
spher_corr(double val)
{
  if (val < 1) return(1.0 - 1.5 * val + 0.5 * pow(val, 3));
  else return(0.0);
}

/* Exponential class */

static double
exp_corr(double val)
{
  return(exp(-val));
}

/* Gaussian class */

static double
Gaus_corr(double val)
{
  return(exp(-(val * val)));
}

/* Linear class */

static double
lin_corr(double val)
{
  if (val < 1) return(1.0 - val);
  else return(0.0);
}

/* Rational class */

static double
ratio_corr(double val)
{
  double val2 = val * val;
  return(1/(1+val2));
}

/* Dummy class */
static double
dummy_corr(double val)
{
  PROBLEM "Unknown spatial correlation class"
    RECOVER(NULL_ENTRY);
  return(0.0);              /* can't occur but this keeps -Wall option happy */
}


/* methods for the virtual class */

static void
spatial_mat(double *par, double *dist, longint *n, longint *nug, 
	    double (*corr)(double ), double *mat)
{
  longint i, j, np1 = *n + 1;
  double aux, *sdist, ratio = 1.0;
  sdist = dist;
  if (*nug) ratio = par[1];
  for(i = 0; i < *n; i++) {
    mat[i * np1] = 1.0;
    for(j = i + 1; j < *n; j++, sdist++) {
      aux = *sdist / *par;
      *(mat + i + j * (*n)) = *(mat + j + i * (*n)) = ratio * corr(aux);
    }
  }
}
  
void 
spatial_matList(double *par, longint *nug, double *dist, longint *pdims,
		double *minD, double *mat)
{
  longint i, M = pdims[1], spClass = pdims[2], *len = pdims + 4, 
    *start = len + M;
  double aux, (*corr)(double ) = dummy_corr;
  /* parameter assumed in unconstrained form */
  par[0] = exp(par[0]);
  if (*nug == 1) {
    aux = exp(par[1]);
    par[1] = 1 / (1.0 + aux);	/* 1 - nugget */
  }
  switch(spClass) {
  case 1:			/* spherical */
    corr = spher_corr;	
    par[0] += *minD;
    break;
  case 2:			/* exponential */
    corr = exp_corr;	
    break;
  case 3:			/* Gaussian */
    corr = Gaus_corr;	
    break;
  case 4:			/* linear */
    corr = lin_corr;
    par[0] +=  *minD;
    break;
  case 5:			/* rational quadratic */
    corr = ratio_corr;
    break;
  default: {PROBLEM "Unknown spatial correlation class"
	      RECOVER(NULL_ENTRY);}
    break;
  }
  for(i = 0; i < M;  i++) {
    spatial_mat(par, dist + start[i], &len[i], nug, corr, mat);
    mat += len[i] * len[i];
  }
}

static void 
spatial_fact(double *par, double *dist, longint *n, longint *nug, 
	     double (*corr) (double ), double *mat, 
	     double *logdet)
{
  longint job = 11L, info, zero = 0L, i, nsq = *n * (*n), np1 = *n + 1;
  double *work = Calloc(*n, double), *work1 = Calloc(nsq, double);
  spatial_mat(par, dist, n, nug, corr, mat);
  F77_CALL(chol) (mat, n, work, &zero, &zero, &info);
  for(i = 0; i < *n; i++) {
    work1[i * np1] = 1;
    F77_CALL(dtrsl) (mat, n, n, work1 + i * (*n), &job, &info);
    *logdet -= log(fabs(mat[i * np1]));
  }
  Memcpy(mat, work1, nsq);
  Free(work); Free(work1);
}

void 
spatial_factList(double *par, longint *nug, double *dist, longint *pdims,  
		 double *minD, double *FactorL, double *logdet)
{
  longint i, M = pdims[1], spClass = pdims[2], *len = pdims + 4, 
    *start = len + M;
  double aux, (*corr)(double ) = dummy_corr;
  /* parameter assumed in unconstrained form */
  par[0] = exp(par[0]);
  if (*nug == 1) {
    aux = exp(par[1]);
    par[1] = 1 / (1.0 + aux);	/* 1 - nugget */
  }

  switch(spClass) {
  case 1:			/* spherical */
    corr = spher_corr;	
    par[0] += *minD;
    break;
  case 2:			/* exponential */
    corr = exp_corr;	
    break;
  case 3:			/* Gaussian */
    corr = Gaus_corr;	
    break;
  case 4:			/* linear */
    corr = lin_corr;
    par[0] +=  *minD;
    break;
  case 5:			/* rational quadratic */
    corr = ratio_corr;
    break;
  default: {PROBLEM "Unknown spatial correlation class"
	      RECOVER(NULL_ENTRY);}
    break;
  }
  for(i = 0; i < M;  i++) {
    spatial_fact(par, dist + start[i], &len[i], nug, corr, FactorL, logdet);
    FactorL += len[i] * len[i];
  }
}

void
spatial_recalc(double *Xy, longint *pdims, longint *ZXcol, double *par, 
	       double *dist, double *minD, longint *nug, double *logdet)
{
  longint N = pdims[0], M = pdims[1], spClass = pdims[2], 
    *len = pdims + 4, *start = len + M, i;
  double aux, (*corr)(double ) = dummy_corr, *sXy;
  /* parameter assumed in unconstrained form */
  par[0] = exp(par[0]);
  if (*nug == 1) {
    aux = exp(par[1]);
    par[1] = 1 / (1.0 + aux);	/* 1 - nugget */
  }

  switch(spClass) {
  case 1:			/* spherical */
    corr = spher_corr;	
    par[0] += *minD;
    break;
  case 2:			/* exponential */
    corr = exp_corr;	
    break;
  case 3:			/* Gaussian */
    corr = Gaus_corr;	
    break;
  case 4:			/* linear */
    corr = lin_corr;
    par[0] +=  *minD;
    break;
  case 5:			/* rational quadratic */
    corr = ratio_corr;
    break;
  default: {PROBLEM "Unknown spatial correlation class"
	      RECOVER(NULL_ENTRY);}
    break;
  }

  for(i = 0, sXy = Xy; i < M;  i++) {
    double *Factor = Calloc(len[i] * len[i], double);
    spatial_fact(par, dist + start[i], &len[i], nug, corr, Factor, logdet);
    mult_mat(sXy, N, Factor, len[i], len[i], len[i], sXy, N, *ZXcol);
    sXy += len[i];
    Free(Factor);
  }
}

#ifndef R_S_H

/* nlme functions and variables     */

				/* Variables that must be initialized */
				/* in fit_nlme */
static int conv_failure;
static longint corOpt, varOpt;
static double *newtheta, *theta, *incr, *add_ons, new_objective, objective;

static struct {			/* Gauss-Newton nonlinear least squares */
  double *residuals, *gradient, *DmHalf, *corFactor, *varWeights, 
    minFactor, tolerance;
  longint nparTot, nrdof, *sgroups, *corDims, maxIter, *npar;
  dimPTR dd;
} nlme;

static long int *
make_sequential(longint *dest, longint *src, longint n)
{                    
  /*  copy the pattern from src to dest */
  /*  but in sequential values starting */
  /*  from 0 */
  long val = 0, *ret = dest, sval;
  if (n <= 0) return dest;
  sval = *src++; *dest++ = val;
  while (--n) {
    if (*src != sval) {sval = *src; val++;}
    src++;
    *dest++ = val;
  }
  return ret;
}

static void
nlme_init(double *DmHalf, double *corFactor, double *varWeights, 
	longint	*groups, longint *pDims, longint *corDims, double *settings)
{
  longint i, *src;
  nlme.DmHalf = DmHalf;
  nlme.corFactor = corFactor;
  nlme.varWeights = varWeights;
  nlme.corDims = corDims;
  nlme.dd = dims(pDims);
  nlme.npar = Calloc(nlme.dd->Q + 1, longint);
  for(i = 0, nlme.nparTot = 0; i <= nlme.dd->Q; i++) {
    nlme.npar[i] = (nlme.dd->ncol)[i] * (nlme.dd->ngrp)[i];
    nlme.nparTot += nlme.npar[i];
  }
  nlme.nrdof = nlme.dd->N - nlme.nparTot;
  nlme.sgroups = groups;
  for(i = 0, src = nlme.sgroups; i < nlme.dd->Q; i++) {
    make_sequential(src, src, nlme.dd->N);
    src += nlme.dd->N;
  }
  nlme.maxIter = (long) settings[0];
  nlme.minFactor = settings[1];
  nlme.tolerance = settings[2];
}

static double
nlme_objective(void)
{
  double RSS, *srcB;
  longint i, j;
  if(varOpt) {			/* variance function correction */
    for(i = 0; i < nlme.dd->N; i++) {
      for(j = 0; j < nlme.dd->ZXcols; j++) {
	*(nl_results[0] + i + j * nlme.dd->N) *= nlme.varWeights[i];
      }
    }
  }
  if(corOpt) {			/* correlation structure correction */
    corStruct_recalc(nl_results[0],nlme.corDims,&nlme.dd->ZXcols,nlme.corFactor);
  }
  nlme.residuals = nl_results[0] + (nlme.dd->ZXcols - 1) * nlme.dd->N;
  nlme.gradient = nl_results[0];
  RSS = d_sum_sqr(nlme.residuals, nlme.dd->N);
  for(i = 0, srcB = newtheta; i < nlme.dd->Q; i++) { 
    double *work = Calloc(nlme.npar[i], double);
    mult_mat(work, (nlme.dd->ncol)[i], nlme.DmHalf + (nlme.dd->DmOff)[i],
	     (nlme.dd->ncol)[i], (nlme.dd->ncol)[i], (nlme.dd->ncol)[i],
	     srcB, (nlme.dd->ncol)[i], (nlme.dd->ngrp)[i]);
    RSS += d_sum_sqr(work, nlme.npar[i]);
    srcB += nlme.npar[i];
    Free(work);
  }
  return(RSS);
}

static void
nlme_correction(void)
{
  long i, j, nq;
  for(i = 0; i < nlme.dd->N; i++) {
    for(j = 0, nq = 0; j < nlme.dd->Q; j++) {
      nlme.residuals[i] += d_dot_prod(nlme.gradient + (nlme.dd->ZXoff)[j][0] + i, 
	    nlme.dd->N, theta+nq+(nlme.sgroups+j*nlme.dd->N)[i]*nlme.dd->ncol[j], 
				    1L, (nlme.dd->ncol)[j]);
      nq += nlme.npar[j];
    }
    nlme.residuals[i] += 
      d_dot_prod(nlme.gradient + (nlme.dd->ZXoff)[nlme.dd->Q][0] + i,
		 nlme.dd->N, theta + nq, 1L, (nlme.dd->ncol)[nlme.dd->Q]);
  }
}

static double
nlme_RegSS(double *incr, double *mat)
{
  long i, j, nq;
  double regSS = 0, aux = 0, *src;
  for(i = 0; i < nlme.dd->N; i++) {
    for(j = 0, nq = 0, aux = 0; j < nlme.dd->Q; j++) {
       aux += d_dot_prod(mat + (nlme.dd->ZXoff)[j][0] + i, 
	   nlme.dd->N, incr+nq+(nlme.sgroups+j*nlme.dd->N)[i]*nlme.dd->ncol[j], 
		       1L, nlme.dd->ncol[j]);
       nq += nlme.npar[j];
    }
    aux += d_dot_prod(mat + (nlme.dd->ZXoff)[nlme.dd->Q][0] + i,
		      nlme.dd->N, incr + nq, 1L, (nlme.dd->ncol)[nlme.dd->Q]);
    regSS += aux * aux;
  }
  for(i = 0, src = incr; i < nlme.dd->Q; i++) { 
    double *work = Calloc((nlme.npar)[i], double);
    mult_mat(work, (nlme.dd->ncol)[i], nlme.DmHalf + (nlme.dd->DmOff)[i],
	     (nlme.dd->ncol)[i], (nlme.dd->ncol)[i], (nlme.dd->ncol)[i],
	     src, (nlme.dd->ncol)[i], (nlme.dd->ngrp)[i]);
    regSS += d_sum_sqr(work, (nlme.npar)[i]);
    src += (nlme.npar)[i];
    Free(work);
  }
  return(regSS);
}

static double
nlme_increment(double *incr)
{
  double regSS, *dc = Calloc(nlme.dd->Srows * nlme.dd->ZXcols, double),
    ll, *dest, *src, 
    *auxGrad = Calloc(nlme.dd->N * (nlme.dd->ZXcols - 1), double);
  long i, j, start, RML = 0;
  if (!sqrt_eps) sqrt_eps = sqrt(DOUBLE_EPS);
  Memcpy(auxGrad, nlme.gradient, (nlme.dd->ZXcols - 1) * nlme.dd->N);
  nlme_correction();
  ll = internal_loglik(nlme.dd, nl_results[0], nlme.DmHalf, &RML, dc, DNULLP);
  internal_estimate(nlme.dd, dc);
  src = dc +  (nlme.dd->ZXcols - 1) * nlme.dd->Srows;
  dest = incr;
  for(i = 0, start = 0; i <= nlme.dd->Q; i++) {
    for(j = 0; j < (nlme.dd->ngrp)[i]; j++) {
      Memcpy(dest, src + ((nlme.dd->SToff)[i][j] - start), (nlme.dd->ncol)[i]);
      dest += (nlme.dd->ncol)[i];
    }
    start += (nlme.dd->ncol)[i] * nlme.dd->Srows;
  }
  for(i = 0; i < nlme.nparTot; i++) {
    incr[i] -= theta[i];
  }
  regSS = nlme_RegSS(incr, auxGrad);	/* Regression Sum of Squares */
  Free(dc); Free(auxGrad);
  return(sqrt(((double) nlme.nrdof) * regSS / 
	      (((double) nlme.nparTot) * (new_objective - regSS))));
}

static longint
nlme_iterate(void)
{
  double factor, criterion;
  longint iteration;
  Memcpy(newtheta, theta, nlme.nparTot);
  spread(theta, nlme.nparTot); eval_model(TRUE);
  new_objective = objective = nlme_objective();
  conv_failure = 0;
  for (factor = 1.0, iteration = 1; iteration <= nlme.maxIter;
       iteration++) {		/* outer iteration loop */
    criterion = nlme_increment(incr); /* increment and convergence criterion */
    if (conv_failure) return(iteration); /* Unable to make increment  */
    if (criterion < nlme.tolerance) return(iteration); /* successful completion */
    do {			/* inner loop for acceptable step size */
      if (factor < nlme.minFactor) {
	conv_failure = 1;
	return(iteration);
      }
      Memcpy(newtheta, theta, nlme.nparTot);
      d_axpy(newtheta, factor, incr, nlme.nparTot);
      spread(newtheta, nlme.nparTot); 
      eval_model(TRUE); 
      new_objective = nlme_objective();
      if (conv_failure) return(iteration); /* unable to evaluate objective */
      factor /= 2.0;
    } while (new_objective >= objective);
    factor *= 4.0;
    if (factor > 1.0)
      factor = 1.0;
    objective = new_objective;
    Memcpy(theta, newtheta, nlme.nparTot);
  }
  conv_failure = 2;		/* Maximum number of iterations exceeded */
  return(iteration - 1);
}

static void
nlme_wrapup(void)
{
  spread(theta, nlme.nparTot); 
  eval_model(TRUE);
  Memcpy(add_ons, nl_results[0], nlme.dd->N * nlme.dd->ZXcols);
  objective = nlme_objective();
  Free(nlme.npar);
  dimFree(nlme.dd);
}

void
fit_nlme(double *ptheta, double *pDmHalf, longint *pgroups, 
	longint *pdims, double *pcorFactor, double *pvarWeights, 
	longint *pcorDims, double *settings, double *additional, 
	longint *pcorOpt, longint *pvarOpt)
{
  if(!sqrt_eps) sqrt_eps = sqrt(DOUBLE_EPS);
  corOpt = *pcorOpt;
  varOpt = *pvarOpt;
  theta = ptheta; 
  add_ons = additional;
  nlme_init(pDmHalf, pcorFactor, pvarWeights, pgroups, pdims, pcorDims, settings);
  newtheta = Calloc(nlme.nparTot, double);
  incr =  Calloc(nlme.nparTot, double);
  settings[4] = (double) nlme_iterate();
  nlme_wrapup();
  settings[3] = conv_failure;
  settings[5] = objective;
  Free(newtheta);
  Free(incr);
}

/*  Quinidine model  */

void 
nlme_one_comp_first (long int *nrow, double *Resp, double *inmat)
{
  long int i, j, nn = *nrow, mm = 0;
  double v, cl, *tl = Calloc(nn, double), *ds = Calloc(nn, double),
    *Subject, *Time, *Dose, *V, *Cl,
         sl = DOUBLE_EPS;	/* sl is last subject number, usually */
				/* an integer but passed as double. */
				/* It is started at an unlikely value. */
  Subject = inmat;
  Time = inmat + nn;
  Dose = inmat + 2 * nn;
  V = inmat + 3 * nn;
  Cl = inmat + 4 * nn;
  for(i = nn; i--; Resp++, Subject++, Time++, Dose++, V++, Cl++) {
    v = *V; cl = *Cl;
    *Resp = 0;
    if (*Subject != sl) {	/* new Subject */
      if (is_na_DOUBLE(Dose)) {
	PROBLEM
	  "First observation on an individual must have a dose"
	  RECOVER(NULL_ENTRY);
      }
      sl = *Subject;
      mm = 0;
      tl[mm] = *Time;
      ds[mm] = *Dose;
    } else {			/* same Subject */
      if (!is_na_DOUBLE(Dose)) { /* Dose measurement */
	mm++;
	tl[mm] = *Time;
	ds[mm] = *Dose;
      } else {			/* Concentration measurement */
	for(j = 0; j <= mm; j++) {
	  *Resp += ds[j] * exp(-cl * (*Time - tl[j]) / v) / v;
	}
      }
    }
  }
  Free(ds); Free(tl);
}

/* nls functions used in predict.nls */

static double
est_delta(double *x, longint i)
{
  double xx;
  if(!sqrt_eps) sqrt_eps = sqrt(DOUBLE_EPS);
  if(!xlower) xlower = 100.*DOUBLE_XMIN;
                                /* should sometime use the strategy of */
				/* the grd routine in dmnf */
  xx = fabs(x[i]);
  if (xx < xlower) return sqrt_eps;
  else return xx*sqrt_eps;
}


void
nls_diff_gradient(longint *pnpar, longint *pnobs, double *theta, double
		  *base, double *gradient, longint *pneg) 
{
  longint i, j, npar = *pnpar, nobs = *pnobs, neg = *pneg; double xx, *gcol, di;
  for(i=0, gcol = gradient; i<npar; i++, gcol += nobs) {
    xx = theta[i];
    theta[i] = theta[i] + (di = est_delta(theta,i));
    spread(theta, npar); eval_model(FALSE);
    if(neg) di = -di;		/* want negative gradient? */
    for(j=0; j<nobs; j++)
      gcol[j] = (nl_results[0][j] - base[j])/di;
    theta[i] = xx;
  }
}


/* gnls functions and variables     */

				/* Variables that must be initialized */
				/* in do_nlme */

static struct {			/* Gauss-Newton nonlinear least squares */
  double *residuals, *gradient, *corFactor, *varWeights, minFactor, tolerance;
  longint npar, N, nrdof, ncol, *corDims, maxIter;
} gnls;

static void
gnls_init(longint *dims, double *corFactor, double *varWeights, 
	  longint *corDims, double *settings)
{
  gnls.corFactor = corFactor;
  gnls.varWeights = varWeights;
  gnls.corDims = corDims;
  gnls.npar = dims[0];
  gnls.N = dims[1];
  gnls.nrdof = gnls.N - gnls.npar;
  gnls.ncol = gnls.npar + 1;
  gnls.maxIter = (long) settings[0];
  gnls.minFactor = settings[1];
  gnls.tolerance = settings[2];
}

static double
gnls_objective(void)
{
  longint i, j;
  if(varOpt) {			/* variance function correction */
    for(i = 0; i < gnls.N; i++) {
      for(j = 0; j < gnls.ncol; j++) {
	*(nl_results[0] + i + j * gnls.N) *= gnls.varWeights[i];
      }
    }
  }
  if(corOpt) {			/* correlation structure correction */
    corStruct_recalc(nl_results[0], gnls.corDims, &gnls.ncol, gnls.corFactor);
  }
  gnls.residuals = nl_results[0] + gnls.npar * gnls.N;
  gnls.gradient = nl_results[0];
  return(d_sum_sqr(gnls.residuals, gnls.N));
}

static double
gnls_RegSS(double *incr)
{
  long i;
  double regSS = 0, aux;
  for(i = 0; i < gnls.N; i++) {
    aux = d_dot_prod(gnls.gradient + i, gnls.N, incr, 1L, gnls.npar);
    regSS += aux * aux;
  }
  return(regSS);
}

static double
gnls_increment(double *incr)
{
  double regSS = 0, *auxRes;
  QRptr aQR;
  long i, j, start;
  if (!sqrt_eps) sqrt_eps = sqrt(DOUBLE_EPS);
  auxRes = Calloc(gnls.N, double);
  Memcpy(auxRes, gnls.residuals, gnls.N);
  aQR = QR(gnls.gradient, gnls.N, gnls.N, gnls.npar);
  QRsolve(aQR, gnls.residuals, gnls.N, 1L, incr, gnls.npar);
  QRqty(aQR, auxRes, gnls.N, 1L);
  for(i=0; i < gnls.npar; i++) {
    regSS += auxRes[i] * auxRes[i];
  }
  QRfree(aQR);
  Free(auxRes);
  return(sqrt(((double) gnls.nrdof) * regSS / 
	      ((double) gnls.npar) * (new_objective - regSS)));
}

static longint
gnls_iterate(void)
{
  double factor, criterion;
  longint iteration;
  Memcpy(newtheta, theta, gnls.npar);
  spread(theta, gnls.npar); eval_model(TRUE);
  new_objective = objective = gnls_objective();
  conv_failure = 0;
  for (factor = 1.0, iteration = 1; iteration <= gnls.maxIter;
       iteration++) {		/* outer iteration loop */
    criterion = gnls_increment(incr); /* increment and convergence
					 criterion */
    if (conv_failure) return(iteration); /* Unable to make increment */
    if (criterion < gnls.tolerance) return(iteration); /* successful completion */
    do {			/* inner loop for acceptable step size */
      if (factor < gnls.minFactor) {
	conv_failure = 1;
	return(iteration);
      }
      Memcpy(newtheta, theta, gnls.npar);
      d_axpy(newtheta, factor, incr, gnls.npar);
      spread(newtheta, gnls.npar); 
      eval_model(TRUE); 
      new_objective = gnls_objective();
      if (conv_failure) return(iteration); /* unable to evaluate objective */
      factor /= 2.0;
    } while (new_objective >= objective);
    factor *= 4.0;
    if (factor > 1.0)
      factor = 1.0;
    objective = new_objective;
    Memcpy(theta, newtheta, gnls.npar);
  }
  conv_failure = 2;		/* Maximum number of iterations exceeded */
  return(iteration - 1);
}

static void
gnls_wrapup(void)
{
  spread(theta, gnls.npar); 
  eval_model(TRUE);
  Memcpy(add_ons, nl_results[0] + gnls.npar * gnls.N, gnls.N);
  objective = gnls_objective();
}

void
fit_gnls(double *ptheta, longint *pdims, double *pcorFactor, double
	 *pvarWeights, longint *pcorDims, double *settings,
	 double *additional, longint *pcorOpt, longint *pvarOpt)
{
  if(!sqrt_eps) sqrt_eps = sqrt(DOUBLE_EPS);
  corOpt = *pcorOpt;
  varOpt = *pvarOpt;
  theta = ptheta; 
  add_ons = additional;
  gnls_init(pdims, pcorFactor, pvarWeights, pcorDims, settings);
  newtheta = Calloc(gnls.npar, double);
  incr =  Calloc(gnls.npar, double);
  settings[4] = (double) gnls_iterate();
  gnls_wrapup();
  settings[3] = conv_failure;
  settings[5] = objective;
  Free(newtheta);
  Free(incr);
}

#endif /* R_S_H */
