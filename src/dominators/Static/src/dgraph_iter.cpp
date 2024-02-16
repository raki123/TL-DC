/*******************************************************************************
 * Copyright (C) 2006,
 * Loukas Georgiadis, Robert E. Tarjan, and Renato F. Werneck
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 ******************************************************************************/

#include "dgraph.h"

/*--------------------------------------------------------------
 | iterative dominators algorithm of Cooper, Harvey, and Kennedy
 | - dominators initalized with zero
 *--------------------------------------------------------------*/

void DominatorGraph::idfs (int r, int *idom) {
	int v, i, new_idom, N;
	int bsize = n+1;
	int *buffer = new int [2*bsize];
	int *post2label = &buffer[0]; //post-dfs ids to original label
	int *dom = &buffer[bsize];    //dominators (indexed by post-ids)

	resetcounters();

	int *label2post = idom; //idom will not be used until later
	N = postDFS (r, label2post, post2label); //get post-ids
	bool changed;

	for (v=n; v>=0; v--) dom[v] = 0;
	dom[N] = N;

	/*-----------
	 | main loop
	 *----------*/
	do {
		inci(); //increment number of iterations (operation count)
		changed = false;

		for (i=N-1; i>0; i--) { //reverse post-order
			new_idom = 0; //using dom[i] is not faster

			/*----------------------------------------------------
			 | for each incoming arc (v,w), compute nca between v
			 | and the current candidate dominator of w
			 *---------------------------------------------------*/
			int *p, *stop;
			getInBounds (post2label[i], p, stop);
			for (; p<stop; p++) {
				int v = label2post[*p]; //v is the source of the arc
				incc();
				if (dom[v]) {           //find nca between current dom and v 
					new_idom = (new_idom ? intersect(v,new_idom,dom) : v);
					incc();
				}
			}
					
			/*-------------------------------------------------------
			 | if new dominator found, update dom and mark as changed
			 *------------------------------------------------------*/
			incc();
			if (new_idom > dom[i]) {
				dom[i] = new_idom;
				changed = true;
			}
		}
	} while (changed);

	/*-----------------------------------------------------------
	 | restore idoms: unreachable nodes are already zero because
	 | array is shared with label2post
	 *----------------------------------------------------------*/
	idom[r] = r;
	for (i=N-1; i>0; i--) idom[post2label[i]] = post2label[dom[i]];
	
	delete [] buffer;
}


/*-------------------------------------------------
 | iterative dominators algorithm
 | - dominators initialized with parent in bfs
 | - vertices visited in direct pre-order
 *------------------------------------------------*/

void DominatorGraph::ibfs (int r, int *idom) {
	int bsize = n+1;
	int *buffer = new int [2*bsize];
	int *pre2label = &buffer[0];
	int *dom       = &buffer[bsize];
	int *label2pre = idom;          //indexed by label
	resetcounters();

	//find pre-ids, initialize dom with parents in BFS tree
	int N = preBFSp (r, label2pre, pre2label, dom); 

	bool changed = true;

	while (changed) {
		inci(); //increment iteration counter
		changed = false;

		// process vertices in preorder
		for (int i=2; i<=N; i++) {
			int new_idom = dom[i];

			/*----------------------------------------------------
			 | for each incoming arc (v,w), compute nca between v
			 | and the current candidate dominator of v
			 *---------------------------------------------------*/
			int *p, *stop;
			getInBounds (pre2label[i], p, stop);
			for (; p<stop; p++) {
				int v = label2pre[*p];
				incc();
				if (v) new_idom = preIntersect (v, new_idom, dom);
			}
					
			/*-----------------------------------------------------------
			 | if we new dominator found, update dom and mark as changed
			 *----------------------------------------------------------*/
			incc();
			if (new_idom!=dom[i]) {
				dom[i] = new_idom;
				changed = true;
			}
		}
	}

	//get dominators
	for (int i=N; i>0; i--) idom[pre2label[i]] = pre2label[dom[i]];
	delete [] buffer;
}
