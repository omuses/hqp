/*
 *  meschext_hl.C
 *   -- some Meschach add-ons and extensions
 */

/*
    Copyright (C) 1996--1998  Hartmut Linke

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Library General Public
    License as published by the Free Software Foundation; 
    version 2 of the License.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Library General Public License for more details.

    You should have received a copy of the GNU Library General Public
    License along with this library (file COPYING.LIB);
    if not, write to the Free Software Foundation, Inc.,
    59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include	<math.h>
#include        "Meschach.h"


#define MINROWLEN 	10
#define INSERT_IDX(idx) 	(-(idx + 2))
#define MEM_MOVE(src,dst,len) 	MEM_COPY(src,dst,len)
#define SPROW_IDX(r,k) \
   ((k < r->elt[0].col) ? -2 : \
    ((k > r->elt[r->len-1].col) ? -(r->len+2) : sprow_idx(r,k)))

/* 
   sprow_ins_val (Copyright Ruediger Franke)
   -- insert the idx's entry into sparse row r with col j
   -- derived frome sprow_set_val (Copyright D.E. Steward, Z. Leyk)
*/

static Real sprow_ins_val(SPROW *r, int idx, Real val, int j, int type)
{
  int  new_len;
  
  if (!r)
    error(E_NULL,"sprow_ins_val");
  
  /* shift & insert new value */
  if (idx < 0)
    idx = INSERT_IDX(idx);
  if ( r->len >= r->maxlen ) {
    r->len = r->maxlen;
    new_len = max(2*r->maxlen+1,5);
    if (mem_info_is_on()) {
      mem_bytes(type, r->maxlen*sizeof(row_elt),
		new_len*sizeof(row_elt)); 
    }
    
    r->elt = RENEW(r->elt,new_len,row_elt);
    if ( ! r->elt )        /* can't allocate */
      error(E_MEM,"sprow_ins_val");
    r->maxlen = new_len;
  }
  
  if ( idx < r->len )
    MEM_MOVE((char *)(&(r->elt[idx])),(char *)(&(r->elt[idx+1])),
             (r->len-idx)*sizeof(row_elt));

  r->len++;
  r->elt[idx].col = j;

  return r->elt[idx].val = val;
}

/*
 * derived from sprow_inprod (Copyright Ruediger Franke)
 */
Real sprow_inprod_full(const SPROW *r1, const SPROW *r2)
{
  Real sum;
  int j_idx1, j_idx2, j1, j2;
  int len1, len2;
  row_elt *elt1, *elt2;

  j_idx1=0;
  j_idx2=0;
  len1 = r1->len;
  len2 = r2->len;
  elt1 = r1->elt;
  elt2 = r2->elt;
  j1 = elt1->col;
  j2 = elt2->col;
  sum = 0.0;
  while (j_idx1 < len1 && j_idx2 < len2) {
    if (j1 < j2) {
      j_idx1 ++;
      elt1 ++;
      j1 = elt1->col;
    }
    else if (j1 > j2) {
      j_idx2 ++;
      elt2 ++;
      j2 = elt2->col;
    }
    else {
      sum += elt1->val * elt2->val;

      j_idx1 ++;
      elt1 ++;
      j1 = elt1->col;

      j_idx2 ++;
      elt2 ++;
      j2 = elt2->col;
    }
  }

  return sum;
}

/* sprow_cpy -- copies r1 into r_out
   -- cannot be done in-situ
   -- type must be SPMAT or SPROW depending on
      whether r_out is a row of a SPMAT structure
      or a SPROW variable
   -- returns r_out */

SPROW	*sprow_cpy(const SPROW *r1, SPROW *r_out, int type)
{
   int	len1, len_out;
   
   if ( ! r1 )
     error(E_NULL,"sprow_cpy");
   if ( ! r_out )
     error(E_NULL,"sprow_cpy");
   if ( r1 == r_out )
     error(E_INSITU,"sprow_cpy");
   
   /* Initialise */
   len1 = r1->len;	
   len_out = r_out->maxlen;

   if(len1 > len_out) {

     if (mem_info_is_on()) {
       mem_bytes(type, len_out*sizeof(row_elt),
		 len1*sizeof(row_elt)); 
     }
    
     r_out->elt = RENEW(r_out->elt,len1,row_elt);
     if ( ! r_out->elt )        /* can't allocate */
      error(E_MEM,"sprow_cpy");
     r_out->maxlen = len1;

   };

   MEM_COPY((char *)(r1->elt),(char *)(r_out->elt),
	    len1*sizeof(row_elt));

   r_out->len = len1;

   return r_out;
}


/*
  derived from sprow_smlt (Copyright D.E. Steward, Z. Leyk)
  sprow_smlt_full -- sets r_out <- alpha*r1 
   -- can be in situ
   -- returns r_out
*/

SPROW	*sprow_smlt_full(SPROW *r1, double alpha, SPROW *r_out, int type)
{
   int	    i, len;
   row_elt  *elt, *elt_out;
   
   if ( ! r1 )
     error(E_NULL,"sprow_smlt r1");
   if ( ! r_out )
     error(E_NULL,"sprow_smlt r_out");
   
   /* Initialise */
   len = r1->len;
   elt = r1->elt;

   r_out = sprow_resize(r_out,len,type);  
   elt_out = r_out->elt;

   for ( i = 0 ; i < len; elt++,elt_out++,i++)
     elt_out->val = alpha*elt->val;

   return r_out;
}


/* 
   derived from sprow_mltadd (Copyright D.E. Steward, Z. Leyk)
   -- sets r_out <- r1 + alpha.r2
   -- cannot be in situ
   -- type must be SPMAT or SPROW depending on
      whether r_out is a row of a SPMAT structure
      or a SPROW variable
   -- returns r_out 
*/
SPROW	*sprow_mltadd_full(SPROW *r1, SPROW *r2, double alpha,
			   SPROW *r_out, int type)
{
   int	    idx1, idx2, idx_out, len1, len2, len_out;
   row_elt  *elt1, *elt2, *elt_out;
   
   if ( ! r1 || ! r2 )
     error(E_NULL,"sprow_mltadd");
   if ( r1 == r_out || r2 == r_out )
     error(E_INSITU,"sprow_mltadd");
   if ( ! r_out )
     error(E_NULL,"sprow_mltadd");
   
   /* Initialise */
   len1 = r1->len;
   len2 = r2->len;
   len_out = r_out->maxlen;

   idx1 = 0;
   idx2 = 0;
   idx_out = 0;

   /* idx1 = idx2 = idx_out = 0; */
   elt1    = r1->elt;
   elt2    = r2->elt;
   elt_out = r_out->elt;
   
   while ( idx1 < len1 || idx2 < len2 )
   {
      if ( idx_out >= len_out )
      {   /* r_out is too small */
	 r_out->len = idx_out;
	 r_out = sprow_xpd(r_out,0,type);
	 len_out = r_out->maxlen;
	 elt_out = &(r_out->elt[idx_out]);
      }
      if ( idx2 >= len2 || (idx1 < len1 && elt1->col <= elt2->col) )
      {
	 elt_out->col = elt1->col;
	 elt_out->val = elt1->val;
	 if ( idx2 < len2 && elt1->col == elt2->col )
	 {
	    elt_out->val += alpha*elt2->val;
	    elt2++;		idx2++;
	 }
	 elt1++;	idx1++;
      }
      else
      {
	 elt_out->col = elt2->col;
	 elt_out->val = alpha*elt2->val;
	 elt2++;	idx2++;
      }
      elt_out++;	idx_out++;
   }

   r_out->len = idx_out;
   
   return r_out;
}

/*
  derived from sprow_smlt (Copyright D.E. Steward, Z. Leyk)
  sprow_smlt_full -- sets r1 <- alpha*r1 from idx to row->len
   -- in situ
*/
SPROW	*sprow_smlt_idx(SPROW *r1, int idx, double alpha)
{
   int	    len;
   row_elt  *elt;
   
   if ( ! r1 )
     error(E_NULL,"sprow_smlt");

   len = r1->len;
   elt = r1->elt + idx;
   for ( ; idx < len; elt++,idx++)
     elt->val *= alpha;

   return r1;
}

/* spbkp_mltadd -- sets r_out <- r1 + alpha.r2 (Copyright Ruediger Franke)
 * -- derived frome sprow_mltadd (Copyright D.E. Steward, Z. Leyk)
 * -- specialized for BKP: 
 *     + r1 starts with index r1->diag
 *     + start index for r2 is passed as argument (r2_idx0)
 *     + r_out is filled from index 0; r_out->diag is set
 * -- cannot be in situ
 * -- type must be SPMAT or SPROW depending on
 *    whether r_out is a row of a SPMAT structure
 *    or a SPROW variable
 * -- returns r_out
 */
static SPROW *spbkp_mltadd(const SPROW *r1, const SPROW *r2, int r2_idx0,
			   Real s, SPROW *r_out, int type)
{
  int	idx1, idx2, idx_out, len1, len2, len_out;
  row_elt	*elt1, *elt2, *elt_out;
  
  if (!r1 || !r2)
    error(E_NULL,"spbkp_mltadd");
  if (r1 == r_out || r2 == r_out)
    error(E_INSITU,"spbkp_mltadd");
  if (!r_out)
    /* don't use sprow_get() because of Meschach memory management */
    /* r_out = sprow_get(MINROWLEN); */
    error(E_NULL,"spbkp_mltadd");
  
  /* Initialise */
  len1 = r1->len;
  len2 = r2->len;
  len_out = r_out->maxlen;
  idx1    = r1->diag;
  idx2    = r2_idx0;
  idx_out = 0;
  if (idx1 >= 0 || idx2 >= 0)
    r_out->diag = 0;
  else
    r_out->diag = -2;
  idx1    = (idx1 < 0)? INSERT_IDX(idx1): idx1;
  idx2    = (idx2 < 0)? INSERT_IDX(idx2): idx2;
  elt1    = &(r1->elt[idx1]);
  elt2    = &(r2->elt[idx2]);
  elt_out = &(r_out->elt[idx_out]);

  while (idx1 < len1 || idx2 < len2) {
    if (idx_out >= len_out) {
      /* r_out is too small */
      r_out->len = idx_out;
      r_out = sprow_xpd(r_out,0,type);
      len_out = r_out->maxlen;
      elt_out = &(r_out->elt[idx_out]);
    }
    if (idx2 >= len2 || (idx1 < len1 && elt1->col <= elt2->col)) {
      elt_out->col = elt1->col;
      elt_out->val = elt1->val;
      if (idx2 < len2 && elt1->col == elt2->col) {
	elt_out->val += s*elt2->val;
	elt2++;
	idx2++;
      }
      elt1++;
      idx1++;
    }
    else {
      elt_out->col = elt2->col;
      elt_out->val = s*elt2->val;
      elt2++;
      idx2++;
    }
    elt_out++;
    idx_out++;
  }
  r_out->len = idx_out;
  
  return r_out;
}


/* 
   derived from sprow_ip (Copyright D.E. Steward, Z. Leyk)
   -- finds the (partial) inner product of a pair of sparse rows
   -- uses a "merging" approach & assumes column ordered rows
   -- row indices for inner product are all < lim 
*/

double	sprow_ip(SPROW *row1, SPROW *row2, int lim)
{
	int	idx1, idx2, len1, len2, tmp;
	row_elt	*elts1, *elts2;
	Real	sum;

	elts1 = row1->elt;	elts2 = row2->elt;
	len1 = row1->len;	len2 = row2->len;

	sum = 0.0;

	if ( len1 <= 0 || len2 <= 0 )
		return 0.0;
	if ( elts1->col >= lim || elts2->col >= lim )
		return 0.0;

	/* use sprow_idx() to speed up inner product where one row is
		much longer than the other */
	idx1 = idx2 = 0;
	if ( len1 > 2*len2 )
	{
		idx1 = SPROW_IDX(row1,elts2->col);
		idx1 = (idx1 < 0) ? -(idx1+2) : idx1;
		if ( idx1 < 0 )
			error(E_UNKNOWN,"sprow_ip");
		len1 -= idx1;
	}
	else if ( len2 > 2*len1 )
	{
		idx2 = SPROW_IDX(row2,elts1->col);
		idx2 = (idx2 < 0) ? -(idx2+2) : idx2;
		if ( idx2 < 0 )
			error(E_UNKNOWN,"sprow_ip");
		len2 -= idx2;
	}
	if ( len1 <= 0 || len2 <= 0 )
		return 0.0;

	elts1 = &(elts1[idx1]);		elts2 = &(elts2[idx2]);


	for ( ; ; )	/* forever do... */
	{
		if ( (tmp=elts1->col-elts2->col) < 0 )
		{
		    len1--;		elts1++;
		    if ( ! len1 || elts1->col >= lim )
			break;
		}
		else if ( tmp > 0 )
		{
		    len2--;		elts2++;
		    if ( ! len2 || elts2->col >= lim )
			break;
		}
		else
		{
		    sum += elts1->val * elts2->val;
		    len1--;		elts1++;
		    len2--;		elts2++;
		    if ( ! len1 || ! len2 ||
				elts1->col >= lim || elts2->col >= lim )
			break;
		}
	}

	return sum;
}

double	sprow_sqr_full(SPROW *row)
{
	row_elt	*elts;
	int	idx, len;
	Real	sum, tmp;

	sum = 0.0;
	elts = row->elt;	len = row->len;
	for ( idx = 0; idx < len; idx++, elts++ )
	{
		tmp = elts->val;
		sum += tmp*tmp;
	}

	return sum;
}


/* sp_mmt_mlt -- compute A.A^T where A is a given sparse matrix */
SPMAT *sp_mmt_mlt(const SPMAT *A, SPMAT *AAT)
{

  int     i, i_min, i_max, j, j_min, j_max, m, hlen, vlen,
          rlen_aat, clen_aat;
  Real    tmp;
  SPROW   *hrow, *vrow, *c_aat, *r_aat;
  row_elt *helt, *velt, *celt_aat, *relt_aat;
  

  if ( ! A )
    error(E_NULL,"sp_AAT");
  if ( ! AAT )
    error(E_NULL,"sp_AAT");

  if(AAT->m != A->m || AAT->n != A->m)
    error(E_SIZES,"sp_AAT");

  for (i = 0; i < AAT->m; i++)
    AAT->row[i].len = 0;

  m = A->m;

  for( i = 0; i < m; i++) {

    hrow = A->row + i;
    hlen = hrow->len;
    helt = hrow->elt;
    i_min = helt[0].col;
    i_max = helt[hlen-1].col;

    r_aat = AAT->row + i;
    rlen_aat = r_aat->len;

    if (rlen_aat == r_aat->maxlen)
      r_aat = sprow_xpd(r_aat,r_aat->maxlen+10,TYPE_SPMAT);

    relt_aat = r_aat->elt + rlen_aat;
    relt_aat->val = sprow_sqr_full(hrow);
    relt_aat->col = i;
    rlen_aat++;

    for(j = i+1; j < m; j++) {

	vrow = A->row+j;
	vlen = vrow->len;
	velt = vrow->elt;
	j_min = velt[0].col;
	j_max = velt[vlen-1].col;

        if((i_min <= j_max) || (j_min <= i_max)) {

	  tmp = sprow_inprod_full(hrow,vrow);

	  if(tmp != 0) {

	    if (rlen_aat == r_aat->maxlen)
	      r_aat = sprow_xpd(r_aat,r_aat->maxlen+10,TYPE_SPMAT);

	    relt_aat = r_aat->elt + rlen_aat;
	    relt_aat->val = tmp;
	    relt_aat->col = j;
	    rlen_aat++;

	    c_aat = AAT->row + j;
	    clen_aat = c_aat->len;

	    if (clen_aat == c_aat->maxlen)
	      c_aat = sprow_xpd(c_aat,c_aat->maxlen+10,TYPE_SPMAT);
	    
	    celt_aat = c_aat->elt + clen_aat;
	    celt_aat->val = tmp;
	    celt_aat->col = i;
	    c_aat->len = clen_aat + 1;;

	  };

	};

    };

    r_aat->len = rlen_aat;

  };

  return AAT;

}    

/* 
   sp_mmt2_mlt -- compute A.A^T where A is a given sparse matrix 
*/
SPMAT *sp_mmt2_mlt(const SPMAT *A,const SPMAT *AT, SPMAT *AAT)
{

  int     i, idx, j, len, m;
  Real    s;
  SPMAT   *swap_mat;
  SPROW   *row, *swap, tmp_row;
  row_elt *elt;
  

  if ( ! A || !AT || !AAT)
    error(E_NULL,"sp_AAT");

  if(AAT->m != A->m || AAT->n != A->m)
    error(E_SIZES,"sp_AAT");

  if(AAT->m != AT->n || AAT->n != AT->n)
    error(E_SIZES,"sp_AAT");

  for (i = 0; i < AAT->m; i++)
    AAT->row[i].len = 0;

  AAT->flag_diag = 0;
  AAT->flag_col = 0;

  m = A->m;
  swap_mat = sp_get(1,m,MINROWLEN);
  swap = swap_mat->row;

  for( i = 0; i < m; i++) {

    row = A->row + i;
    elt = row->elt;
    len = row->len;
    idx = 0;

    for(; idx < len; idx++, elt++) {

      j = elt->col;
      s = elt->val;

      sprow_mltadd_full(AAT->row+i,AT->row+j,s,swap,TYPE_SPMAT);

      /* exchange swap with row #j */
      MEM_COPY(AAT->row+i, &tmp_row, sizeof(SPROW));
      MEM_COPY(swap, AAT->row+i, sizeof(SPROW));
      MEM_COPY(&tmp_row, swap, sizeof(SPROW));
      
    };

  };

  sp_free(swap_mat);

  return AAT;

} 
 
/* 
   sp_mmt_mlt_u -- compute A.A^T where A is a given sparse matrix 
   only the upper part of AAT is filled
*/
SPMAT *sp_mmt2_mlt_u(const SPMAT *A,const SPMAT *AT, SPMAT *AAT)
{

  int     i, idx, idx1, j, len, m;
  Real    s;
  SPMAT   *swap_mat;
  SPROW   *row, *row1, *swap, tmp_row;
  row_elt *elt;
  

  if ( ! A || !AT || !AAT)
    error(E_NULL,"sp_AAT");

  if(AAT->m != A->m || AAT->n != A->m)
    error(E_SIZES,"sp_AAT");

  if(AAT->m != AT->n || AAT->n != AT->n)
    error(E_SIZES,"sp_AAT");

  for (i = 0; i < AAT->m; i++) {
    AAT->row[i].len = 0;
    AAT->row[i].diag = 0;
  }
  
  AAT->flag_diag = 0;
  AAT->flag_col = 0;

  m = A->m;
  swap_mat = sp_get(1,m,MINROWLEN);
  swap = swap_mat->row;

  for( i = 0; i < m; i++) {

    row = A->row + i;
    elt = row->elt;
    len = row->len;
    idx = 0;

    for(; idx < len; idx++, elt++) {

      j = elt->col;
      s = elt->val;
      row1 = AT->row+j;
      idx1 = SPROW_IDX(row1,i);
      if(idx1 < 0) idx1 = INSERT_IDX(idx1);
      if(idx1 > 0) idx1--;

      spbkp_mltadd(AAT->row+i,row1,0,s,swap,TYPE_SPMAT);

      /* exchange swap with row #j */
      MEM_COPY(AAT->row+i, &tmp_row, sizeof(SPROW));
      MEM_COPY(swap, AAT->row+i, sizeof(SPROW));
      MEM_COPY(&tmp_row, swap, sizeof(SPROW));
      
    };

  };

  sp_free(swap_mat);

  return AAT;

}    


/* spLMsolve 
   -- solve L.OUT=B where L is a sparse lower triangular matrix,
   -- OUT,B are sparse matrices
   -- returns OUT; operation may not be in-situ
   if L is not blocktriangular, OUT may be full !!!

*/
SPMAT *spLMsolve(SPMAT *L, const SPMAT *B,SPMAT *OUT)
{

  int	   i, j_idx, in_situ, n;
  SPMAT    *tmp_mat;
  SPROW    *row, *tmp1_row, *tmp;
  row_elt  *elt;


  if ( L == SMNULL || B == SMNULL )
    error(E_NULL,"spLsolve");
  if ( L->m != L->n )
    error(E_SQUARE,"spLsolve");
  if ( B->m != L->m )
    error(E_SIZES,"spLsolve");
  //  if ( ! L->flag_diag )
    sp_diag_access(L);
  if( !OUT) 
    OUT = sp_get(B->m,B->n,MINROWLEN);
  if( OUT->m != B->m || OUT->n != B->n)
    sp_resize(OUT,B->m,B->n);
  if(B == OUT)
    in_situ = 1;
  else
    in_situ = 0;
    

  tmp_mat = sp_get(1,B->n,MINROWLEN);
  tmp = tmp_mat->row;
    
  n = L->n;
  
  for ( i = 0; i < n; i++ ) {

    row = &(L->row[i]);
    elt = row->elt;

    if(!in_situ)
      sprow_cpy(B->row+i,OUT->row+i,TYPE_SPMAT);

    for ( j_idx = 0; j_idx < row->len; j_idx++, elt++ ) {
      if ( elt->col >= i )
	break;

      tmp1_row = &(OUT->row[elt->col]);
      sprow_mltadd_full(OUT->row+i,tmp1_row,-elt->val,tmp,TYPE_SPMAT);
      sprow_cpy(tmp,OUT->row+i,TYPE_SPMAT);
    }
    if ( row->diag >= 0 ) {
      sprow_smlt_full(OUT->row+i,1.0/row->elt[row->diag].val,
		      OUT->row+i,TYPE_SPMAT);
    }
    else
      error(E_SING,"spLMsolve");
  };

  sp_free(tmp_mat);

  return OUT;
}


/* spLDLTsolve -- solve L.D.L^T.out=b where L is a sparse matrix,
	-- out, b dense vectors
	-- returns out; operation may be in-situ */
VEC	*spLDLTsolve(SPMAT *L,const VEC *b,VEC *out)
{
  int	   col, i, j, j_idx, max_col, n;
  SPROW	   *row;
  row_elt  *elt;
  Real	   diag_val, sum, *out_ve;

  if ( L == SMNULL || b == VNULL )
    error(E_NULL,"spCHsolve");
  if ( L->m != L->n )
    error(E_SQUARE,"spCHsolve");
  if ((int)b->dim != L->m )
    error(E_SIZES,"spCHsolve");
  
  if ( ! L->flag_diag )
    sp_diag_access(L);
  
  col = 0;

  out = v_copy(b,out);
  out_ve = out->ve;

  /* forward substitution: solve L.x=b for x */
  n = L->n;

  for ( i = 0; i < n; i++ ) {
    sum = out_ve[i];
    row = &(L->row[i]);
    elt = row->elt;

    for ( j_idx = 0; j_idx < row->len; j_idx++, elt++ ) {
      if ( elt->col >= i )
	break;
      sum -= elt->val*out_ve[elt->col];
    }
  
    out_ve[i] = sum;
    if ( row->diag < 0 )  error(E_SING,"spLDLTsolve");

  };


  /* backward substitution: solve L^T.out = D^1.x for out */
  for ( i = n-1; i >= 0; i-- ) {

    row = &(L->row[i]);
    /* Note that row->diag >= 0 by above loop */

    elt = row->elt;
    max_col = elt[row->len-1].col + 1;
    if(max_col > col) col = max_col;

    elt = &(row->elt[row->diag]);
    diag_val = elt->val;
    sum = out_ve[i]/diag_val;

    for(j = i+1; j < col; j++) {

    /* scan down column */
      row = &(L->row[j]);
      j_idx = SPROW_IDX(row,i);
      if(j_idx >= 0) {
	elt = &(row->elt[j_idx]);
	sum -= elt->val*out_ve[j];
      };

    };
    
    out_ve[i] = sum;
    
  };

  return out;
}

SPMAT *spCHOLfac(SPMAT *A)
/*
 * -- factorize A in situ into A = LL'
 */
{
  int	i, j, n;
  int  	idx, len;
  Real	aii;
  Real	s;
  SPMAT *swap_mat;
  SPROW *row, *swap, tmp_row;
  row_elt *elt;

  if (!A )
    error(E_NULL, "spCHOLfac");
  if (A->n != A->m)
    error(E_SIZES, "spCHOLfac");

  n = A->n;

  /* don't use sprow_get because of Meschach memory management */
  /* swap = sprow_get(MINROWLEN); */
  swap_mat = sp_get(1, n, MINROWLEN);
  swap = swap_mat->row;

  if (!A->flag_diag)
    sp_diag_access(A);
  A->flag_col = 0;

  for (i = 0; i < n-1; i++) {

    row = A->row + i;
    idx = row->diag;

    if (idx >= 0) {
      aii = row->elt[idx].val;
      if( aii > 0.0 )
	aii = row->elt[idx].val = sqrt(aii);
      else
	error(E_POSDEF,"spCHfactor");
      idx ++;
    }
    else
      error(E_POSDEF,"spCHfactor");

    sprow_smlt_idx(row,idx,1.0/aii);

    elt = row->elt + idx;
    len = row->len;
    for (; idx < len; idx++, elt++) {
      j = elt->col;
      s = elt->val;
      if (s != 0.0) {
	spbkp_mltadd(A->row+j, row, idx, -s, swap, TYPE_SPMAT);
	/* exchange swap with row #j */
	MEM_COPY(A->row+j, &tmp_row, sizeof(SPROW));
	MEM_COPY(swap, A->row+j, sizeof(SPROW));
	MEM_COPY(&tmp_row, swap, sizeof(SPROW));
      }
    }
  }

  if (n > 0) {
    row = A->row + n-1;
    idx = row->diag;

    if (idx >= 0) {
      aii = row->elt[idx].val;
      if( aii > 0.0 )
	row->elt[idx].val = sqrt(aii);
      else
	error(E_POSDEF,"spCHfactor");
    }
    else
      error(E_POSDEF,"spCHfactor");
  }

  /* sprow_free(swap); */
  sp_free(swap_mat);

  A->flag_diag = 1;

  return A;
}


SPMAT *spMODCHOLfac(SPMAT *A,VEC *b,double eps)
/*
 * -- factorize A in situ into A = LL'
 */
{
  int	i, j, n;
  int  	idx, len;
  Real	aii;
  Real	s;
  SPMAT *swap_mat;
  SPROW *row, *swap, tmp_row;
  row_elt *elt;

  if (!A )
    error(E_NULL, "spMODCHOLfac");
  if (A->n != A->m)
    error(E_SIZES, "spMODCHOLfac");
  if ((int)A->n != (int)b->dim)
    error(E_SIZES, "spMODCHOLfac");

  n = A->n;

  /* don't use sprow_get because of Meschach memory management */
  /* swap = sprow_get(MINROWLEN); */
  swap_mat = sp_get(1, n, MINROWLEN);
  swap = swap_mat->row;

  if (!A->flag_diag)
    sp_diag_access(A);
  A->flag_col = 0;

  for (i = 0; i < n-1; i++) {

    row = A->row + i;
    idx = row->diag;

    if (idx >= 0) {
      aii = row->elt[idx].val;
      if( aii > eps ) {
	aii = row->elt[idx].val = sqrt(aii);
	idx ++;	
	sprow_smlt_idx(row,idx,1.0/aii);

	elt = row->elt + idx;
	len = row->len;
	for (; idx < len; idx++, elt++) {
	  j = elt->col;
	  s = elt->val;
	  if (s != 0.0) {
	    spbkp_mltadd(A->row+j, row, idx, -s, swap, TYPE_SPMAT);
	    /* exchange swap with row #j */
	    MEM_COPY(A->row+j, &tmp_row, sizeof(SPROW));
	    MEM_COPY(swap, A->row+j, sizeof(SPROW));
	    MEM_COPY(&tmp_row, swap, sizeof(SPROW));
	  };
	};
      }
      else {
	row->elt[idx].val = 1.0;
	idx ++;
	sprow_smlt_idx(row,idx,0.0);
        b->ve[i] = 0.0;
      };
    }
    else {
      sprow_ins_val(row,idx,1.0,i,TYPE_SPMAT);
      idx ++;
      sprow_smlt_idx(row,idx,0.0);
      b->ve[i] = 0.0;
    };
  }

  if (n > 0) {
    row = A->row + n-1;
    idx = row->diag;

    if (idx >= 0) {
      aii = row->elt[idx].val;
      if( aii > eps )
	row->elt[idx].val = sqrt(aii);
      else {
	row->elt[idx].val = 1.0;
	b->ve[n-1] = 0.0;
      };
    }
    else {
      sprow_ins_val(row,idx,1.0,n-1,TYPE_SPMAT);
      b->ve[n-1] = 0.0;
    };
  }

  /* sprow_free(swap); */
  sp_free(swap_mat);

  A->flag_diag = 1;

  return A;
}


/* 
   derived from spCHsolve (Copyright D.E. Steward, Z. Leyk)
   -- solve L.L^T.out=b where L is a sparse matrix,
   -- out, b dense vectors
   -- returns out; operation may be in-situ 
*/
VEC	*spCHOLsol(SPMAT *L,VEC *b,VEC *x)
{
  int	   i, idx, len, n;
  SPROW	   *row;
  row_elt  *elt;
  Real	   *x_ve, tmp;

  if ( L == SMNULL || b == VNULL )
    error(E_NULL,"spCHsolve");
  if ( L->m != L->n )
    error(E_SQUARE,"spCHsolve");
  if ((int)b->dim != L->m )
    error(E_SIZES,"spCHsolve");
  if ( ! L->flag_diag )
    sp_diag_access(L);
  
  x = v_copy(b,x);
  x_ve = x->ve;

  /* forward substitution: solve L.x=b for x */
  n = L->n;

  for(i = 0; i < n-1; i++) {

    row = L->row + i;
    idx = row->diag;

    if(idx >= 0 ) {

      len = row->len;
      elt = row->elt + idx;

      if(elt->val > 0.0) {
	x_ve[i] /= elt->val;
	tmp = x_ve[i];
	idx++;
	elt++;
      }
      else error(E_SING,"forward neg spCHsolve");
    }
    else error(E_SING,"forward idx spCHsolve");

    for(; idx < len; idx++, elt++)
      x_ve[elt->col] -= tmp * elt->val;

  };

  if (n > 0) {
    row = L->row + n - 1;
    idx = row->diag;

    if(idx >= 0 ) {

      elt = row->elt + idx;

      if(elt->val > 0.0) {
	x_ve[n-1] /= elt->val*elt->val;
      }
      else error(E_SING,"backward neg spCHsolve");
    }
    else error(E_SING,"backward idx spCHsolve");
  }

  /* backward substitution: solve L^T.out = x for out */
  for ( i = n-2; i >= 0; i-- ) {
    
    row = L->row + i;
    idx = row->diag;
    len = row->len;
    elt = row->elt + idx;
    tmp = elt->val;
    idx++;
    elt++;

    for(;idx < len; idx++, elt++)
      x_ve[i] -= elt->val*x_ve[elt->col];

    x_ve[i] = x_ve[i]/tmp;

  };

  return x;
}

/* spLMsolve 
   -- solve L.OUT=B where L is a sparse upper triangular matrix,
   -- OUT,B are sparse matrices
   -- returns OUT; operation may not be in-situ
   if L is not blocktriangular, OUT may be full !!!

*/
SPMAT *spLTMsolve(SPMAT *L,SPMAT *B,SPMAT *OUT)
{

  int	   i, idx, j, len, n;
  Real     s;
  SPMAT    *swap_mat;
  SPROW    *row, *swap, tmp_row;
  row_elt  *elt;


  if ( L == SMNULL || B == SMNULL )
    error(E_NULL,"spLTMsolve");
  if ( L->m != L->n )
    error(E_SQUARE,"spLTMsolve");
  if ( B->m != L->m )
    error(E_SIZES,"spLTMsolve");
  if ( ! L->flag_diag )
    sp_diag_access(L);
  if( !OUT) 
    OUT = sp_get(B->m,B->n,MINROWLEN);
  if( OUT->m != B->m || OUT->n != B->n)
    sp_resize(OUT,B->m,B->n);
  if(B != OUT)
    OUT = sp_copy3(B,OUT);

  OUT->flag_diag = 0;
  OUT->flag_col = 0;

   n = L->n;
   
  /* don't use sprow_get because of Meschach memory management */
  /* swap = sprow_get(MINROWLEN); */
  swap_mat = sp_get(1, n, MINROWLEN);
  swap = swap_mat->row;

  L->flag_col = 0;
  if(!L->flag_diag)
    sp_diag_access(L);
     
  for ( i = 0; i < n; i++ ) {

    row = L->row + i;
    idx = row->diag;

    if(idx >= 0 ) {

      len = row->len;
      elt = row->elt + idx;
      s = elt->val;

      if(s > 0.0) {
	sprow_smlt_full(OUT->row+i,1.0/s,OUT->row+i,TYPE_SPMAT);
	idx++;
	elt++;
      }
      else error(E_SING,"spLTMsolve");
    }
    else error(E_SING,"spLTMsolve");

    for(; idx < len; idx++, elt++) {

      j = elt->col;
      s = elt->val;

      sprow_mltadd_full(OUT->row+j,OUT->row+i,-s,swap,TYPE_SPMAT);

      /* exchange swap with row j */
      MEM_COPY(OUT->row+j, &tmp_row, sizeof(SPROW));
      MEM_COPY(swap, OUT->row+j, sizeof(SPROW));
      MEM_COPY(&tmp_row, swap, sizeof(SPROW));
    }
  }

  sp_free(swap_mat);

  return OUT;
}

/* spmm_mlt -- sparse matrix-matrix multiplication */
MAT  *spmm_mlt(SPMAT *A, MAT *B, MAT *OUT)
{

  int	   i, j, m, n, idx, len;
  Real	   **B_v, s;
  SPROW    *row;
  row_elt  *elt;

  if ( A==(SPMAT *)NULL || B==(MAT *)NULL )
    error(E_NULL,"spmm_mlt");
  if ((int)A->n != (int)B->m )
    error(E_SIZES,"spmm_mlt");
  if ( B == OUT )
    error(E_INSITU,"spmm_mlt");
  
  m = A->m;
  n = B->n;

  B_v = B->me;
  
  if ( OUT==(MAT *)NULL || (int)OUT->m != (int)A->m || 
       (int)OUT->n != (int)B->n )
    OUT = m_resize(OUT,A->m,B->n);
  
  m_zero(OUT);
  
  for (i = 0; i < m; i++) {

    row = A->row + i;
    len = row->len;
    elt = row->elt;
    idx = 0;

    for (; idx < len; idx++, elt++) {

      j = elt->col;
      s = elt->val;

      if (s != 0.0) __mltadd__(OUT->me[i],B_v[j],s,n);
 
    };

  };
  
  return OUT;

}


/* spm_m_spmtr_mlt -- sparse matrix-matrix multiplication */
MAT  *spm_m_spmtr_mlt(SPMAT *A, MAT *B, MAT *OUT)
{

  int	   i, j, k, m, n, idx, len;
  MAT      *mtmp;
  VEC      *vtmp;
  Real	   **B_me, **mtmp_me,**OUT_me, *vtmp_ve,s;
  SPROW    *row;
  row_elt  *elt;

  if ( A==(SPMAT *)NULL || B==(MAT *)NULL )
    error(E_NULL,"m_mlt");
  if ((int)A->n != (int)B->m )
    error(E_SIZES,"m_mlt");
  if ( B == OUT )
    error(E_INSITU,"m_mlt");
  
  if ( OUT==(MAT *)NULL || (int)OUT->m != (int)A->m || 
       (int)OUT->n != (int)A->m )
    OUT = m_resize(OUT,A->m,A->n);
  m_zero(OUT);
    
  mtmp = m_get(A->n,A->m);
  vtmp = v_get(A->n);

  m = A->m;
  n = B->n;

  B_me = B->me;
  mtmp_me = mtmp_me;
  vtmp_ve = vtmp_ve;

  for (i = 0; i < m; i++) {

    row = A->row + i;
    len = row->len;
    elt = row->elt;
    idx = 0;

    __zero__(vtmp_ve,n);

    for (; idx < len; idx++, elt++) {

      j = elt->col;
      s = elt->val;

      if (s != 0.0) __mltadd__(vtmp_ve,B_me[j],s,n);
 
    };

    for(k = 0; k < n; k++) mtmp_me[k][i] = vtmp_ve[k];

  };

  vtmp = v_resize(vtmp,m);
  vtmp->ve = vtmp->ve;
  OUT_me = OUT->me;

  for (i = 0; i < m; i++) {

    row = A->row + i;
    len = row->len;
    elt = row->elt;
    idx = 0;

    __zero__(vtmp_ve,m);

    for (; idx < len; idx++, elt++) {

      j = elt->col;
      s = elt->val;

      if (s != 0.0) __mltadd__(vtmp_ve,mtmp_me[j],s,n);
 
    };

    for(k = 0; k < m; k++) OUT_me[k][i] = vtmp_ve[k];

  };

  M_FREE(mtmp);
  V_FREE(vtmp);
  
  return OUT;

}
