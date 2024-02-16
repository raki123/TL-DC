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

/*--------------------------------------------------------
 | Simple Lenguauer-Tarjan algorithm (SLT)
 | - ancestor and parent share an array
 | - recursive compress
 | - bucket processed at the beginning of each iteration
 | - vertex v not inserted in bucket if semi[v]==parent[v]
 *--------------------------------------------------------*/

void DominatorGraph::slt (int r, int *idom) {
	int bsize = n+1;
	int *buffer    = new int [6*bsize];
	int *pre2label = &buffer[0];
	int *parent    = &buffer[bsize];
	int *semi      = &buffer[2*bsize];
	int *label     = &buffer[3*bsize];
	int *dom       = &buffer[4*bsize];
	int *ubucket   = &buffer[5*bsize];

	int *label2pre = idom;          //indexed by label

	resetcounters();

	int i;
	for (i=n; i>=0; i--) {
		label[i] = semi[i] = i;
		ubucket[i] = 0;
	}

	//pre-dfs
	int N = preDFSp (r, label2pre, pre2label, parent);

	// process the vertices in reverse preorder 
	for (i=N; i>1; i--) {
		/*--------------------- 
		 | process i-th bucket
		 *--------------------*/
		for (int v=ubucket[i]; v; v=ubucket[v]) {
			rcompress (v, parent, semi, label, i);
			int u = label[v];
			incc();
			dom[v] = (semi[u]<semi[v]) ? u : i;
		}
		//no need to empty the bucket

		/*---------------------------------------------
		 | check incoming arcs, update semi-dominators
		 *--------------------------------------------*/
		int *p, *stop;
		getInBounds (pre2label[i], p, stop);
		for (; p<stop; p++) {
			int v = label2pre[*p];
			incc();
			if (v) {
				int u; 
				incc();
				if (v<=i) {u=v;} //v is an ancestor of i
				else {
					rcompress (v, parent, semi, label, i);
					u = label[v];
				}
				incc();
				if (semi[u]<semi[i]) semi[i] = semi[u];
			}
		}

		/*---------------------------
		 | process candidate semidom
		 *--------------------------*/
		int s = semi[i];
		incc();
		if (s!=parent[i]) { //if semidominator n not parent: add i to s's bucket
			ubucket[i] = ubucket[s]; 
			ubucket[s] = i;
		} else {
			dom[i] = s; //semidominator is parent: s is a candidate dominator
		}
	}

	/*------------------
	 | process bucket 1
	 *-----------------*/
	for (int v=ubucket[1]; v; v=ubucket[v]) dom[v]=1;

	/*---------------
	 | recover idoms 
	 *--------------*/
	dom[1] = 1;
	idom[r] = r;
	for (i=2; i<=N; i++) {
		incc();
		if (dom[i]!=semi[i]) dom[i]=dom[dom[i]]; //make relative absolute
		idom[pre2label[i]] = pre2label[dom[i]];
   	}

	delete [] buffer; //cleanup stuff
}
