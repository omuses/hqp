
/**************************************************************************
**
** Copyright (C) 1993 David E. Stewart & Zbigniew Leyk, all rights reserved.
**
**			     Meschach Library
** 
** This Meschach Library is provided "as is" without any express 
** or implied warranty of any kind with respect to this software. 
** In particular the authors shall not be liable for any direct, 
** indirect, special, incidental or consequential damages arising 
** in any way from use of the software.
** 
** Everyone is granted permission to copy, modify and redistribute this
** Meschach Library, provided:
**  1.  All copies contain this copyright notice.
**  2.  All modified copies shall carry a notice stating who
**      made the last modification and the date of such modification.
**  3.  No charge is made for this software or works derived from it.  
**      This clause shall not be construed as constraining other software
**      distributed on the same medium as this software, nor is a
**      distribution fee considered a charge.
**
***************************************************************************/

/*
  This file contains basic routines which are used by the functions
  in meschach.a etc.
  These are the routines that should be modified in order to take
  full advantage of specialised architectures (pipelining, vector
  processors etc).
  */

static	char	*rcsid = "$Id: machine.c,v 1.1.1.1 2001/03/01 17:18:43 rfranke Exp $";

#include	"machine.h"

/* __ip__ -- inner product */
Real	__ip__(dp1,dp2,len)
register Real	*dp1, *dp2;
int	len;
{
#ifdef VUNROLL
    register int	j, k;
#endif
    register int	i;
    register Real     sum;

    sum = 0.0;
#ifdef VUNROLL
    i = len >> 3;
    j = 0;
    k = len & 7;

    while (i--)
      {
	sum += dp1[j]*dp2[j]; j++;
	sum += dp1[j]*dp2[j]; j++;
	sum += dp1[j]*dp2[j]; j++;
	sum += dp1[j]*dp2[j]; j++;
	sum += dp1[j]*dp2[j]; j++;
	sum += dp1[j]*dp2[j]; j++;
	sum += dp1[j]*dp2[j]; j++;
	sum += dp1[j]*dp2[j]; j++;
      }
    switch(k)
      {
      case 7: sum += dp1[j]*dp2[j]; j++;
      case 6: sum += dp1[j]*dp2[j]; j++;
      case 5: sum += dp1[j]*dp2[j]; j++;
      case 4: sum += dp1[j]*dp2[j]; j++;
      case 3: sum += dp1[j]*dp2[j]; j++;
      case 2: sum += dp1[j]*dp2[j]; j++;
      case 1: sum += dp1[j]*dp2[j]; j++;
      }
#else   
    for ( i = 0; i < len; i++ )
	sum  += dp1[i]*dp2[i];
#endif
    return sum;
}

/* __mltadd__ -- scalar multiply and add c.f. v_mltadd() */
void	__mltadd__(dp1,dp2,s,len)
register Real	*dp1, *dp2;
register Real s;
register int	len;
{
    register int	i;
#ifdef VUNROLL
    register int	j, k;
    i = len >> 3;
    j = 0;
    k = len & 7;

    while (i--)
      {
	dp1[j] += s*dp2[j]; j++;
	dp1[j] += s*dp2[j]; j++;
	dp1[j] += s*dp2[j]; j++;
	dp1[j] += s*dp2[j]; j++;
	dp1[j] += s*dp2[j]; j++;
	dp1[j] += s*dp2[j]; j++;
	dp1[j] += s*dp2[j]; j++;
	dp1[j] += s*dp2[j]; j++;
      }
    switch(k)
      {
      case 7: dp1[j] += s*dp2[j]; j++;
      case 6: dp1[j] += s*dp2[j]; j++;
      case 5: dp1[j] += s*dp2[j]; j++;
      case 4: dp1[j] += s*dp2[j]; j++;
      case 3: dp1[j] += s*dp2[j]; j++;
      case 2: dp1[j] += s*dp2[j]; j++;
      case 1: dp1[j] += s*dp2[j]; j++;
      }
#else
    for ( i = 0; i < len; i++ )
	dp1[i] += s*dp2[i];
#endif
}

/* __smlt__ scalar multiply array c.f. sv_mlt() */
void	__smlt__(dp,s,out,len)
register Real	*dp, *out;
register Real s;
register int	len;
{
    register int	i;
#ifdef VUNROLL
    register int	j, k;
    i = len >> 3;
    j = 0;
    k = len & 7;

    while (i--)
      {
	out[j] = s*dp[j]; j++;
	out[j] = s*dp[j]; j++;
	out[j] = s*dp[j]; j++;
	out[j] = s*dp[j]; j++;
	out[j] = s*dp[j]; j++;
	out[j] = s*dp[j]; j++;
	out[j] = s*dp[j]; j++;
	out[j] = s*dp[j]; j++;
      }
    switch(k)
      {
      case 7: out[j] = s*dp[j]; j++;
      case 6: out[j] = s*dp[j]; j++;
      case 5: out[j] = s*dp[j]; j++;
      case 4: out[j] = s*dp[j]; j++;
      case 3: out[j] = s*dp[j]; j++;
      case 2: out[j] = s*dp[j]; j++;
      case 1: out[j] = s*dp[j]; j++;
      }
#else
    for ( i = 0; i < len; i++ )
	out[i] = s*dp[i];
#endif
}

/* __add__ -- add arrays c.f. v_add() */
void	__add__(dp1,dp2,out,len)
register Real	*dp1, *dp2, *out;
register int	len;
{
    register int	i;
#ifdef VUNROLL
    register int	j, k;

    i = len >> 3;
    j = 0;
    k = len & 7;

    while (i--)
      {
	out[j] = dp1[j]+dp2[j]; j++;
	out[j] = dp1[j]+dp2[j]; j++;
	out[j] = dp1[j]+dp2[j]; j++;
	out[j] = dp1[j]+dp2[j]; j++;
	out[j] = dp1[j]+dp2[j]; j++;
	out[j] = dp1[j]+dp2[j]; j++;
	out[j] = dp1[j]+dp2[j]; j++;
	out[j] = dp1[j]+dp2[j]; j++;
      }
    switch(k)
      {
      case 7: out[j] = dp1[j]+dp2[j]; j++;
      case 6: out[j] = dp1[j]+dp2[j]; j++;
      case 5: out[j] = dp1[j]+dp2[j]; j++;
      case 4: out[j] = dp1[j]+dp2[j]; j++;
      case 3: out[j] = dp1[j]+dp2[j]; j++;
      case 2: out[j] = dp1[j]+dp2[j]; j++;
      case 1: out[j] = dp1[j]+dp2[j]; j++;
      }
#else
    for ( i = 0; i < len; i++ )
	out[i] = dp1[i] + dp2[i];
#endif
}

/* __sub__ -- subtract arrays c.f. v_sub() */
void	__sub__(dp1,dp2,out,len)
register Real	*dp1, *dp2, *out;
register int	len;
{
    register int	i;
#ifdef VUNROLL
    register int	j, k;

    i = len >> 3;
    j = 0;
    k = len & 7;

    while (i--)
      {
	out[j] = dp1[j]-dp2[j]; j++;
	out[j] = dp1[j]-dp2[j]; j++;
	out[j] = dp1[j]-dp2[j]; j++;
	out[j] = dp1[j]-dp2[j]; j++;
	out[j] = dp1[j]-dp2[j]; j++;
	out[j] = dp1[j]-dp2[j]; j++;
	out[j] = dp1[j]-dp2[j]; j++;
	out[j] = dp1[j]-dp2[j]; j++;
      }
    switch(k)
      {
      case 7: out[j] = dp1[j]-dp2[j]; j++;
      case 6: out[j] = dp1[j]-dp2[j]; j++;
      case 5: out[j] = dp1[j]-dp2[j]; j++;
      case 4: out[j] = dp1[j]-dp2[j]; j++;
      case 3: out[j] = dp1[j]-dp2[j]; j++;
      case 2: out[j] = dp1[j]-dp2[j]; j++;
      case 1: out[j] = dp1[j]-dp2[j]; j++;
      }
#else 
    for ( i = 0; i < len; i++ )
	out[i] = dp1[i] - dp2[i];
#endif
}

/* __zero__ -- zeros an array of floating point numbers */
void	__zero__(dp,len)
register Real	*dp;
register int	len;
{
#ifdef CHAR0ISDBL0
    /* if a floating point zero is equivalent to a string of nulls */
    MEM_ZERO((char *)dp,len*sizeof(Real));
#else
    /* else, need to zero the array entry by entry */
    register int	i;
#ifdef VUNROLL
    register int	j, k;

    i = len >> 3;
    j = 0;
    k = len & 7;
    while (i--)
      {
	dp[j] = 0.0; j++;
	dp[j] = 0.0; j++;
	dp[j] = 0.0; j++;
	dp[j] = 0.0; j++;
	dp[j] = 0.0; j++;
	dp[j] = 0.0; j++;
	dp[j] = 0.0; j++;
	dp[j] = 0.0; j++;
      }
    switch(k)
      {
      case 7: dp[j] = 0.0; j++;
      case 6: dp[j] = 0.0; j++;
      case 5: dp[j] = 0.0; j++;
      case 4: dp[j] = 0.0; j++;
      case 3: dp[j] = 0.0; j++;
      case 2: dp[j] = 0.0; j++;
      case 1: dp[j] = 0.0; j++;
      }
#else 
    for ( i = 0; i < len; i++ )
	dp[i] = 0.0;
#endif
#endif
}

