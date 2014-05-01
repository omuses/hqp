/*
		Type definitions for general purpose maths package
*/
/* RCS id: $Header: /home/ruediger/hqp/hqp/meschach/matdef.h,v 1.5 2002/12/09 10:57:47 e_arnold Exp $ */

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
MESCH__BEGIN_DECLS
extern	jmp_buf	restart;
MESCH_API int	ev_err(char *, int, int, char *);
MESCH__END_DECLS
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
			m_error(_err_num,function);	} \
	}

MESCH_API double	__ip__(double *, double *, int);
MESCH_API void	__mltadd__(double *, double *, double, int);
MESCH_API void	__smlt__(double *, double, double *, int);
MESCH_API void	__add__(double *, double *, double *, int);
MESCH_API void	__sub__(double *, double *, double *, int);
MESCH_API void	__zero__(double *, int);
MESCH_API int	set_err_flag(int);
