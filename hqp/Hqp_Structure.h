/*
 * Hqp_Structure.h --
 *  - express sparsity structure of a matrix for permutations
 *  - reference: Weissinger: Sp"arlich besetzte Gleichungssysteme
 *
 * rf, 5/19/94
 */

/*
    Copyright (C) 1994--1998  Ruediger Franke

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

#ifndef Hqp_Structure_H
#define Hqp_Structure_H

#include "Hqp_impl.h"

class Hqp_Program;

class Hqp_BandSize {

 public:

  int	lb;	// number of bands below diagonal
  int	ub;	// number of bands above diagonal
  int	n;	// number of rows (cols)
};

class Hqp_Structure {

 protected:

  // buffer for returning by methods
  Hqp_BandSize	_bandSize;
  IVEC		_neigh_out;

  // data to represent sparsity structure
  int	_nentries;
  int	_n;		// dimension (of symmetric matrix)
  int	*_neighbours;	// list of neighbours of all nodes
  int	_neigh_size;	// size of _neighbours
  int	*_neigh_start;	// offsets into _neighbors of all nodes
  int	*_degree;	// degrees of all nodes
  PERM	*_order;	// hold ordering

  int	add_row(const SPMAT *, int i, int node, int offs=0);
  int	add_col(const SPMAT *, int j, int node, int offs=0);

  void	neigh_grow();
  void	neigh_sort();

  void	realloc(int n, int neigh_size);

 public:

  Hqp_Structure();
  ~Hqp_Structure();

  // initializing
#ifdef SPARSE_COL_ACCESS
  void init_QAC(const SPMAT *Q, const SPMAT *A, const SPMAT *C);
  void init_QA(const SPMAT *Q, const SPMAT *A);
#endif
  void init_spmat(const SPMAT *);

  // performing actions
  void order_rcm();

  // access to results
  int		nentries() {return _nentries;}
  PERM		*px_get(PERM *px=PNULL);
  const IVEC	*neighbours(int node);
  const Hqp_BandSize *bd_size();	// return read only data
};


#endif
