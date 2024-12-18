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

/*---------------------------------------------------------------
 | SEMI-NCA (snca): two-phase algorithm:
 | 1. computes semidominators as in Simple Lengauer-Tarjan (SLT)
 | 2. builds the dominator tree incrementally
 |
 | Notes:
 | - parent and ancestor share an array
 | - recursive compress
 *--------------------------------------------------------------*/

void DominatorGraph::snca (int r, int *idom) {
	int bsize = n+1;
	int *buffer    = new int [5*bsize];
	int *dom       = &buffer[0*bsize]; //not shared
	int *pre2label = &buffer[1*bsize];
	int *parent    = &buffer[2*bsize]; //shared with ancestor
	int *label     = &buffer[3*bsize];
	int *semi      = &buffer[4*bsize];

	int *label2pre = idom;          //indexed by label

	resetcounters();

	//initialize semi and label
	int i;
	for (i=n; i>=0; i--) label[i] = semi[i] = i;

	int N = preDFSp (r, label2pre, pre2label, parent);
	
	/*----------------
	 | semidominators 
	 *---------------*/
	for (i=N; i>1; i--) {
		int *p, *stop;
		dom[i] = parent[i]; //can't put dom and parent together

		//process each incoming arc
		getInBounds (pre2label[i], p, stop);
		for (; p<stop; p++) {
			int v = label2pre[*p];
			if (v) {
				int u;
				incc();
				if (v<=i) {u=v;} //v is an ancestor of i
				else {
					rcompress (v, parent, label, i);
					u = label[v];
				}
				incc();
				if (semi[u]<semi[i]) semi[i] = semi[u];
			}
		}
		label[i] = semi[i];
	}


	/*-----------------------------------------------------------
	 | compute dominators using idom[w]=NCA(I,parent[w],sdom[w])
	 *----------------------------------------------------------*/
	dom[1] = 1;
	idom[r] = r;
	for (i=2; i<=N; i++) {
		int j = dom[i];
		while (j>semi[i]) {j=dom[j]; incc();}
		incc();
		dom[i] = j;
		idom[pre2label[i]] = pre2label[dom[i]];
	}

	//cleanup stuff
	delete [] buffer;
}
