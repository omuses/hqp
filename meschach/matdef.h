/*
		Type definitions for general purpose maths package
*/
/* RCS id: $Header: /home/ruediger/hqp/hqp/meschach/matdef.h,v 1.3 2002/05/01 17:50:39 rfranke Exp $ */

#include "machine.h"

#define	TRUE	1
#define	FALSE	0

#ifndef RS6000
/* this is defined in <sys/types.h> on RS6000/AIX systems */
typedef	unsigned int	u_int;
#endif

/* vector definition */
typedef	struct	{
		u_int	dim, max_dim;
		double	*ve;
		} VEC;
/* matrix definition */
typedef	struct	{
		u_int	m, n;
		u_int	max_m, max_n, max_size;
		double	**me,*base;	/* base is base of alloc'd mem */
		} MAT;
/* permutation definition */
typedef	struct	{
		u_int	size, max_size, *pe;
		} PERM;
/* integer vector definition */
typedef struct	{
		u_int	dim, max_dim;
		int	*ive;
	        } IVEC;

#define	VNULL	((VEC *)NULL)
#define	MNULL	((MAT *)NULL)
#define	PNULL	((PERM *)NULL)
#define	IVNULL	((IVEC *)NULL)

#define	MAXDIM	2001

#ifndef MALLOCDECL
#ifndef ANSI_C
extern	char	*malloc(), *calloc(), *realloc();
#else
extern	void	*malloc(size_t),
		*calloc(size_t,size_t),
		*realloc(void *,size_t);
#endif
#endif

#define	NEW(type)		((type *)calloc(1,sizeof(type)))
#define	NEW_A(size,type)	((type *)calloc((size)>1?(size):1,sizeof(type)))
#define	RENEW(var,num,type) \
	((var)=(type *)((var) ? realloc((var),(num)*sizeof(type)) : \
				calloc((num),sizeof(type))))

/* useful things to have around... */
#define	min(a,b)	((a) < (b) ? (a) : (b))
#define	max(a,b)	((a) > (b) ? (a) : (b))
#define	PI	3.141592653589793
#define	E	2.718281828459045
#define	cp_vec(in,out)	_cp_vec(in,out,0)
#define	cp_mat(in,out)	_cp_mat(in,out,0,0)
#define	set_col(mat,col,vec)	_set_col(mat,col,vec,0)
#define	set_row(mat,row,vec)	_set_row(mat,row,vec,0)

/* for input routines */
#define	MAXLINE	81
char	line[MAXLINE];

/* Error recovery */
#include	<setjmp.h>
extern	jmp_buf	restart;
#ifdef ANSI_C
extern int	ev_err(char *, int, int, char *);
#else
/* extern int	ev_err(); */
#endif
#define	m_error(err_num,fn_name) ev_err(__FILE__,err_num,__LINE__,fn_name)
#define	E_UNKNOWN	0
#define	E_SIZES		1
#define	E_BOUNDS	2
#define	E_MEM		3
#define	E_SING		4
#define	E_POSDEF	5
#define	E_FORMAT	6
#define	E_INPUT		7
#define	E_NULL		8
#define	E_SQUARE	9
#define	E_RANGE		10
#define	E_INSITU2	11
#define	E_INSITU	12
#define	E_ITER		13
#define	E_CONV		14
#define	E_START		15
#define	E_SIGNAL	16
#define	E_INTERN	17
#define	E_EOF		18

#define	EF_EXIT		0
#define	EF_ABORT	1
#define	EF_JUMP		2
#define	EF_SILENT	3
#define	m_tracecatch(ok_part,function) \
	{	jmp_buf _save;	int _err_num, _old_flag; \
		_old_flag = set_err_flag(EF_JUMP); \
		mem_copy(restart,_save,sizeof(jmp_buf)); \
		if ( (_err_num=setjmp(restart)) == 0 ) \
		{	ok_part; \
			set_err_flag(_old_flag); \
			mem_copy(_save,restart,sizeof(jmp_buf));	} \
		else \
		{	set_err_flag(_old_flag);  \
			mem_copy(_save,restart,sizeof(jmp_buf)); \
			error(_err_num,function);	} \
	}

#ifdef ANSI_C
extern double	__ip__(double *, double *, int);
extern void	__mltadd__(double *, double *, double, int);
extern void	__smlt__(double *, double, double *, int);
extern void	__add__(double *, double *, double *, int);
extern void	__sub__(double *, double *, double *, int);
extern void	__zero__(double *, int);
extern int	set_err_flag(int);
#else
extern double	__ip__();
extern void	__mltadd__(), __smlt__(), __add__(), __sub__(), __zero__();
extern int	set_err_flag();
#endif
