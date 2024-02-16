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

/*****************************************************************
 * 
 * Lengauer-Tarjan computation of semi-dominators (simple version)
 *
 *****************************************************************/


#include "dgraph.h"

int DominatorGraph::semi_dominators (int r) {
	int bsize = n+1;
	int *buffer = new int [5*bsize];
	int *label2pre = &buffer[0];
	int *pre2label = &buffer[bsize];
	int *parent    = &buffer[2*bsize];
	//int *ancestor  = &buffer[3*bsize];
	int *label     = &buffer[3*bsize];
	int *semi      = &buffer[4*bsize];

	resetcounters();

	int npdom = 0; //number of vertices dominated by their parent

	int i;
	for (i=0; i<=n; i++) {
		label[i] = semi[i] = i;
		//ancestor[i] = 0;
	}

	int N = preDFSp (r, label2pre, pre2label, parent);

	for (i=N; i>=2; i--) {
		int w = pre2label[i];
		int *p, *stop;
		getInBounds(w,p,stop);
		for (; p<stop; p++) {
			int v = label2pre[*p];
			if (v) {
				int u;
				//if (!ancestor[v]) {u=v;}
				//if (!parent[v]) {u=v;}
				if (v<=i) {u=v;} //u is an ancestor of i
				else {
					//rcompress(v,ancestor,semi,label);
					rcompress(v,parent,semi,label,i);
					u = label[v];
				}
				if (semi[u]<semi[i]) semi[i] = semi[u];
			}
		}

		//if (semi[i]==parent[i]) npdom++;
		incs();
		//ancestor[i] = parent[i];
	}

	delete [] buffer;
	return npdom;
}
