/* ----------------------------------------------------------------- 
 * Originally written by Chris Chu, Iowa State University:
 -------------- FLUTE - Version 3.1 -----------------
                       by
                Chris C.-N. Chu
       Dept. of ECpE, Iowa State University
             Copyright (c) - 2005
  Iowa State University Research Foundation, Inc.
 ----------------------------------------------------
 * Just fantastic code. Thanks Chris! (WPS).
 * Modified by William Swartz for public use and to reduce memory usage
 * and to remove MAXD constraint.  Flute is now dynamic and reentrant!
 * Copyright (c) 2019 William Swartz.  UTD.  All rights reserved.
----------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#define FLUTE_INTERNAL_VARS
#include "flute.h"
#include "flute_malloc.h"

namespace Flute {

/* ***************** STATIC FUNCTION DEFINITIONS ******************* */
static Tree dmergetree(Tree t1, Tree t2) ;
static Tree vmergetree(Tree t1, Tree t2) ;
static Tree hmergetree(Tree t1, Tree t2, int *s) ;
static void local_refinement( FluteState * flute_p, int d, Tree *tp, int p) ;

/* ***************** STATIC VARIABLE DEFINITIONS ******************* */


FluteState * flute_init(  const char *wirelength_vector_file, const char *routing_tree_file )
{
    FILE *fpwv, *fprt;
    FLUTE_SOLNPTR p ;
    int d, i, j, k, kk, ns, nn, ne;
    char ch ;
    unsigned char line[2000], *linep, c;
    unsigned char charnum[256] ;
    unsigned char divd[256], modd[256], div16[256], mod16[256], dsq;
    FluteState * flute_p ;			/* top module */

    flute_p = UTDCALLOC( 1, FLUTEBOX ) ;
    /* was static int numgrp[10]={0,0,0,0,6,30,180,1260,10080,90720}; */
    flute_p->numgrp[0] = 0 ;
    flute_p->numgrp[1] = 0 ;
    flute_p->numgrp[2] = 0 ;
    flute_p->numgrp[3] = 0 ;
    flute_p->numgrp[4] = 6 ;
    flute_p->numgrp[5] = 30 ;
    flute_p->numgrp[6] = 180 ;
    flute_p->numgrp[7] = 1260 ;
    flute_p->numgrp[8] = 10080 ;
    flute_p->numgrp[9] = 90720 ;
    flute_p->num_used = 0 ;
    flute_p->num_alloc = FLUTE_SOLUTION_STORAGE ;
    flute_p->storage = UTDMALLOC( flute_p->num_alloc, FLUTE_SOLNPTR ) ;

    for (i=0; i<=255; i++) {
        if ('0'<=i && i<='9')
            charnum[i] = i - '0';
        else if (i>='A')
            charnum[i] = i - 'A' + 10;
        else // if (i=='$' || i=='\n' || ... )
            charnum[i] = 0;
	div16[i] = i / 16;
	mod16[i] = i % 16;
    }

    fpwv=fopen( wirelength_vector_file, "r");
    if (!(fpwv)) {
        printf("Error in opening %s\n", wirelength_vector_file ) ;
        return(NULL) ;
    }

#if FLUTE_ROUTING==1
    fprt=fopen( routing_tree_file, "r");
    if (!(fprt)) {
        printf("Error in opening %s\n", routing_tree_file );
        return(NULL) ;
    }
#endif


    for (d = 4; d <= FLUTE_D; d++) {
      for (i = 0; i <= 255; i++) {
	divd[i] = i / d;
	modd[i] = i % d;
      }
      dsq = d*d;

      fscanf(fpwv, "d=%d\n", &d);
#if FLUTE_ROUTING==1    
      fscanf(fprt, "d=%d\n", &d);
#endif
      for (k = 0; k < flute_p->numgrp[d]; k++) {
	ns = (int) charnum[fgetc(fpwv)];

	if (ns == 0) { // same as some previous group
	  fscanf(fpwv, "%d\n", &kk);
	  flute_p->numsoln[d][k] = flute_p->numsoln[d][kk];
	  flute_p->lut[d][k] = flute_p->lut[d][kk];
	} else {
	  ch = fgetc(fpwv); // '\n'
	  flute_p->numsoln[d][k] = ns;
	  p = UTDMALLOC( ns, FLUTE_SOLN ) ;
	  if( ++flute_p->num_used > flute_p->num_alloc ){
	    flute_p->num_alloc *= 2 ;
	    flute_p->storage = UTDREALLOC( flute_p->storage, flute_p->num_alloc, FLUTE_SOLNPTR ) ;
	  }
	  flute_p->storage[flute_p->num_used-1] = p ;
	  flute_p->lut[d][k] = p;
	  for (i = 1; i <= ns; i++) {
	    linep = (unsigned char*) fgets((char*) line, 2000, fpwv);
	    p->parent = charnum[*(linep++)];
	    j = 0;
	    while ( (p->seg[j++] = charnum[*(linep++)]) != 0) ; 
	    ;
	      j = 10;
	    while ((p->seg[j--] = charnum[*(linep++)]) != 0);
#if FLUTE_ROUTING==1    
	    nn = 2 * d - 2;
	    ne = 2 * d - 3;
	    fread(line, 1, d - 2, fprt);
	    linep = line;
	    for (j = d; j < nn; j++) {
	      c = *(linep++);
	      if (c >= dsq) {
		c -= dsq;
		p->neighbor[divd[c]] = j;
		ne--;
	      }
	      p->row[j - d] = divd[c];
	      p->col[j - d] = modd[c];
	      p->neighbor[j] = j; // initialized
	    }
	    fread(line, 1, ne + 1, fprt);
	    linep = line; // last char = \n
	    for (j = 0; j < ne; j++) {
	      c = *(linep++);
	      p->neighbor[div16[c]] = mod16[c];
	    }
#endif
	    p++;
	  }
	}
      }
    }
    return( flute_p ) ;
} /* end flute_init() */

void flute_free( FluteState * flute_p )
{
    int i ;			/* counter */
    if( flute_p ){
      for( i = 0 ; i < flute_p->num_used ; i++ ){
	UTDFREE( flute_p->storage[i] ) ;
      }
      UTDFREE( flute_p->storage ) ;
      UTDFREE(flute_p) ;
    }
} /* end flute_free() */

void free_tree( FluteState * flute_p, Tree tree )
{
    if( tree.deg > 0 ){
      UTDFREE( tree.branch ) ;
    }
    tree.deg = 0 ;

} /* end flute_free_tree() */

FLUTE_DTYPE flute_wl( FluteState * flute_p, int d, FLUTE_DTYPE *x, FLUTE_DTYPE *y, int acc)
{
    FLUTE_DTYPE minval, l, xu, xl, yu, yl ;
    FLUTE_DTYPE *xs, *ys ;
    int i, j, k, minidx, degree ;
    int *s ;
    FLUTE_POINTPTR tmpp, *ptp ;
    FLUTE_POINTPTR pt ;

    /* allocate the dynamic pieces on the heap rather than the stack */
    degree = d + 1 ;
    xs = UTDMALLOC( degree, FLUTE_DTYPE ) ;
    ys = UTDMALLOC( degree, FLUTE_DTYPE ) ;
    s = UTDMALLOC( degree, int ) ;
    pt = UTDMALLOC( degree+1, FLUTE_POINT ) ;
    ptp = UTDMALLOC( degree+1, FLUTE_POINTPTR ) ;

    if (d==2)
        l = ADIFF(x[0], x[1]) + ADIFF(y[0], y[1]);
    else if (d==3) {
        if (x[0] > x[1]) {
            xu = max(x[0], x[2]);
            xl = min(x[1], x[2]);
        }
        else {
            xu = max(x[1], x[2]);
            xl = min(x[0], x[2]);
        }
        if (y[0] > y[1]) {
            yu = max(y[0], y[2]);
            yl = min(y[1], y[2]);
        }
        else {
            yu = max(y[1], y[2]);
            yl = min(y[0], y[2]);
        }
        l = (xu-xl) + (yu-yl);
    } else {
        for (i=0; i<d; i++) {
            pt[i].x = x[i];
            pt[i].y = y[i];
            ptp[i] = &pt[i];
        }
        
        // sort x
        for (i=0; i<d-1; i++) {
            minval = ptp[i]->x;
            minidx = i;
            for (j=i+1; j<d; j++) {
                if (minval > ptp[j]->x) {
                    minval = ptp[j]->x;
                    minidx = j;
                }
            }
            tmpp = ptp[i];
            ptp[i] = ptp[minidx];
            ptp[minidx] = tmpp;
        }

#if FLUTE_REMOVE_DUPLICATE_PIN==1
        ptp[d] = &pt[d];
        ptp[d]->x = ptp[d]->y = -999999;
        j = 0;
        for (i=0; i<d; i++) {
            for (k=i+1; ptp[k]->x == ptp[i]->x; k++)
                if (ptp[k]->y == ptp[i]->y)  // pins k and i are the same
                    break;
            if (ptp[k]->x != ptp[i]->x)
                ptp[j++] = ptp[i];
        }
        d = j;
#endif
        
        for (i=0; i<d; i++) {
            xs[i] = ptp[i]->x;
            ptp[i]->o = i;
        }

        // sort y to find s[]
        for (i=0; i<d-1; i++) {
            minval = ptp[i]->y;
            minidx = i;
            for (j=i+1; j<d; j++) {
                if (minval > ptp[j]->y) {
                    minval = ptp[j]->y;
                    minidx = j;
                }
            }
            ys[i] = ptp[minidx]->y;
            s[i] = ptp[minidx]->o;
            ptp[minidx] = ptp[i];
        }
        ys[d-1] = ptp[d-1]->y;
        s[d-1] = ptp[d-1]->o;
        
        l = flutes_wl( flute_p, d, xs, ys, s, acc);
    }
    UTDFREE( xs ) ;
    UTDFREE( ys ) ;
    UTDFREE( s ) ;
    UTDFREE( pt ) ;
    UTDFREE( ptp ) ;

    return l;
} /* end flute_wl() */

// xs[] and ys[] are coords in x and y in sorted order
// s[] is a list of nodes in increasing y direction
//   if nodes are indexed in the order of increasing x coord
//   i.e., s[i] = s_i as defined in paper
// The points are (xs[s[i]], ys[i]) for i=0..d-1
//             or (xs[i], ys[si[i]]) for i=0..d-1

FLUTE_DTYPE flutes_wl_RDP( FluteState * flute_p, int d, FLUTE_DTYPE *xs, FLUTE_DTYPE *ys, int *s, int acc)
{
    int i, j, ss;

    for (i=0; i<d-1; i++) {
        if (xs[s[i]]==xs[s[i+1]] && ys[i]==ys[i+1]) {
            if (s[i] < s[i+1])
                ss = s[i+1];
            else {
                ss = s[i];
                s[i] = s[i+1];
            }
            for (j=i+2; j<d; j++) {
                ys[j-1] = ys[j];
                s[j-1] = s[j];
            }
            for (j=ss+1; j<d; j++)
                xs[j-1] = xs[j];
            for (j=0; j<=d-2; j++)
                if (s[j] > ss) s[j]--;
            i--;
            d--;
        }
    }
    return flutes_wl_ALLD( flute_p, d, xs, ys, s, acc);
} /* end flutes_wl_RDP() */

// For low-degree, i.e., 2 <= d <= D
FLUTE_DTYPE flutes_wl_LD( FluteState * flute_p, int d, FLUTE_DTYPE *xs, FLUTE_DTYPE *ys, int *s)
{
    int k, pi, i, j;
    FLUTE_SOLNPTR rlist;
    FLUTE_DTYPE dd[2*FLUTE_D-2];  // 0..D-2 for v, D-1..2*D-3 for h
    FLUTE_DTYPE minl, sum, l[MPOWV+1];
    
    if (d <= 3)
        minl = xs[d-1]-xs[0]+ys[d-1]-ys[0];
    else {
        k = 0;
        if (s[0] < s[2]) k++;
        if (s[1] < s[2]) k++;
        
        for (i=3; i<=d-1; i++) {  // p0=0 always, skip i=1 for symmetry
            pi = s[i];
            for (j=d-1; j>i; j--)
                if (s[j] < s[i])
                    pi--;
            k = pi + (i+1)*k;
        }
        
        if (k < flute_p->numgrp[d])  // no horizontal flip
            for (i=1; i<=d-3; i++) {
                dd[i]=ys[i+1]-ys[i];
                dd[d-1+i]=xs[i+1]-xs[i];
            }
        else {
            k=2*flute_p->numgrp[d]-1-k;
            for (i=1; i<=d-3; i++) {
                dd[i]=ys[i+1]-ys[i];
                dd[d-1+i]=xs[d-1-i]-xs[d-2-i];
            }
        }
        
        minl = l[0] = xs[d-1]-xs[0]+ys[d-1]-ys[0];
        rlist = flute_p->lut[d][k];
        for (i=0; rlist->seg[i]>0; i++)
            minl += dd[rlist->seg[i]];
        
        l[1] = minl;
        j = 2;
        while (j <= flute_p->numsoln[d][k]) {
            rlist++;
            sum = l[rlist->parent];
            for (i=0; rlist->seg[i]>0; i++)
                sum += dd[rlist->seg[i]];
            for (i=10; rlist->seg[i]>0; i--)
                sum -= dd[rlist->seg[i]];
            minl = min(minl, sum);
            l[j++] = sum;
        }
    }
    
    return minl;
} /* end flutes_wl_LD() */

// For medium-degree, i.e., D+1 <= d
FLUTE_DTYPE flutes_wl_MD(FluteState * flute_p,int d, FLUTE_DTYPE *xs, FLUTE_DTYPE *ys, int *s, int acc)
{
    float pnlty, dx, dy;
    float *score, *penalty ;
    FLUTE_DTYPE xydiff;
    FLUTE_DTYPE ll, minl, extral;
    FLUTE_DTYPE *x1, *x2, *y1, *y2 ;
    FLUTE_DTYPE *distx, *disty ;
    int i, r, p, maxbp, nbp, bp, ub, lb, n1, n2, newacc;
    int ms, mins, maxs, minsi, maxsi, degree ;
    int return_val ;
    int *si, *s1, *s2 ;

    degree = d + 1 ;
    score = UTDMALLOC( 2 * degree, float ) ;
    penalty = UTDMALLOC( degree, float ) ;

    x1 = UTDMALLOC( degree, FLUTE_DTYPE ) ;
    y1 = UTDMALLOC( degree, FLUTE_DTYPE ) ;
    x2 = UTDMALLOC( degree, FLUTE_DTYPE ) ;
    y2 = UTDMALLOC( degree, FLUTE_DTYPE ) ;
    distx = UTDMALLOC( degree, FLUTE_DTYPE ) ;
    disty = UTDMALLOC( degree, FLUTE_DTYPE ) ;
    si = UTDMALLOC( degree, int ) ;
    s1 = UTDMALLOC( degree, int ) ;
    s2 = UTDMALLOC( degree, int ) ;

    if (s[0] < s[d-1]) {
        ms = max(s[0], s[1]);
        for (i=2; i<=ms; i++)
            ms = max(ms, s[i]);
        if (ms <= d-3) {
            for (i=0; i<=ms; i++) {
                x1[i] = xs[i];
                y1[i] = ys[i];
                s1[i] = s[i];
            }
            x1[ms+1] = xs[ms]; 
            y1[ms+1] = ys[ms]; 
            s1[ms+1] = ms+1;
            
            s2[0] = 0;
            for (i=1; i<=d-1-ms; i++)
                s2[i] = s[i+ms]-ms;
            
            return_val = flutes_wl_LMD( flute_p, ms+2, x1, y1, s1, acc) + 
	                 flutes_wl_LMD( flute_p, d-ms, xs+ms, ys+ms, s2, acc);
	    UTDFREE( score ) ;
	    UTDFREE( penalty ) ;
	    UTDFREE( x1 ) ;
	    UTDFREE( y1 ) ;
	    UTDFREE( x2 ) ;
	    UTDFREE( y2 ) ;
	    UTDFREE( distx ) ;
	    UTDFREE( disty ) ;
	    UTDFREE( si ) ;
	    UTDFREE( s1 ) ;
	    UTDFREE( s2 ) ;
	    return return_val;
        }
    }    
    else {  // (s[0] > s[d-1])
        ms = min(s[0], s[1]);
        for (i=2; i<=d-1-ms; i++)
            ms = min(ms, s[i]);
        if (ms >= 2) {
            x1[0] = xs[ms];
            y1[0] = ys[0];
            s1[0] = s[0]-ms+1;
            for (i=1; i<=d-1-ms; i++) {
                x1[i] = xs[i+ms-1];
                y1[i] = ys[i];
                s1[i] = s[i]-ms+1;
            }
            x1[d-ms] = xs[d-1];
            y1[d-ms] = ys[d-1-ms];
            s1[d-ms] = 0;
            
            s2[0] = ms;
            for (i=1; i<=ms; i++)
                s2[i] = s[i+d-1-ms];
            
            return_val = flutes_wl_LMD( flute_p, d+1-ms, x1, y1, s1, acc) +
                         flutes_wl_LMD( flute_p, ms+1, xs, ys+d-1-ms, s2, acc);
	    UTDFREE( score ) ;
	    UTDFREE( penalty ) ;
	    UTDFREE( x1 ) ;
	    UTDFREE( y1 ) ;
	    UTDFREE( x2 ) ;
	    UTDFREE( y2 ) ;
	    UTDFREE( distx ) ;
	    UTDFREE( disty ) ;
	    UTDFREE( si ) ;
	    UTDFREE( s1 ) ;
	    UTDFREE( s2 ) ;
	    return return_val;
        }
    }
    
// Find inverse si[] of s[]
    for (r=0; r<d; r++)
        si[s[r]] = r;

// Determine breaking directions and positions dp[]
    lb=(d-2*acc+2)/4;
    if (lb < 2) lb = 2;
    ub=d-1-lb;

// Compute scores    
#define AAWL 0.6
#define BBWL 0.3
    float CCWL = 7.4/((d+10.)*(d-3.));
    float DDWL = 4.8/(d-1);
    
    // Compute penalty[]    
    dx = CCWL*(xs[d-2]-xs[1]);
    dy = CCWL*(ys[d-2]-ys[1]);
    for (r = d/2, pnlty = 0; r>=0; r--, pnlty += dx)
        penalty[r] = pnlty,  penalty[d-1-r] = pnlty;
    for (r = d/2-1, pnlty = dy; r>=0; r--, pnlty += dy)
        penalty[s[r]] += pnlty,  penalty[s[d-1-r]] += pnlty;
//#define CCWL 0.16
//    for (r=0; r<d; r++)
//        penalty[r] = abs(d-1-r-r)*dx + abs(d-1-si[r]-si[r])*dy;

    // Compute distx[], disty[]
    xydiff = (xs[d-1] - xs[0]) - (ys[d-1] - ys[0]);
    if (s[0] < s[1])
        mins = s[0], maxs = s[1];
    else mins = s[1], maxs = s[0];
    if (si[0] < si[1])
        minsi = si[0], maxsi = si[1];
    else minsi = si[1], maxsi = si[0];
    for (r=2; r<=ub; r++) {
        if (s[r] < mins)
            mins = s[r];
        else if (s[r] > maxs)
            maxs = s[r];
        distx[r] = xs[maxs] - xs[mins];
        if (si[r] < minsi)
            minsi = si[r];
        else if (si[r] > maxsi)
            maxsi = si[r];
        disty[r] = ys[maxsi] - ys[minsi] + xydiff;
    }

    if (s[d-2] < s[d-1])
        mins = s[d-2], maxs = s[d-1];
    else mins = s[d-1], maxs = s[d-2];
    if (si[d-2] < si[d-1])
        minsi = si[d-2], maxsi = si[d-1];
    else minsi = si[d-1], maxsi = si[d-2];
    for (r=d-3; r>=lb; r--) {
        if (s[r] < mins)
            mins = s[r];
        else if (s[r] > maxs)
            maxs = s[r];
        distx[r] += xs[maxs] - xs[mins];
        if (si[r] < minsi)
            minsi = si[r];
        else if (si[r] > maxsi)
            maxsi = si[r];
        disty[r] += ys[maxsi] - ys[minsi];
    }
    
    nbp=0;
    for (r=lb; r<=ub; r++) {
        if (si[r]==0 || si[r]==d-1)
            score[nbp] = (xs[r+1] - xs[r-1]) - penalty[r]
                - AAWL*(ys[d-2]-ys[1]) - DDWL*disty[r];
        else score[nbp] = (xs[r+1] - xs[r-1]) - penalty[r]
                 - BBWL*(ys[si[r]+1]-ys[si[r]-1]) - DDWL*disty[r];
        nbp++;

        if (s[r]==0 || s[r]==d-1)
            score[nbp] = (ys[r+1] - ys[r-1]) - penalty[s[r]]
                - AAWL*(xs[d-2]-xs[1]) - DDWL*distx[r];
        else score[nbp] = (ys[r+1] - ys[r-1]) - penalty[s[r]]
                 - BBWL*(xs[s[r]+1]-xs[s[r]-1]) - DDWL*distx[r];
        nbp++;
    }

    if (acc <= 3)
        newacc = 1;
    else {
        newacc = acc/2;
        if (acc >= nbp) acc = nbp-1;
    }
    
    minl = (FLUTE_DTYPE) INT_MAX;
    for (i=0; i<acc; i++) {
        maxbp = 0;
        for (bp=1; bp<nbp; bp++)
            if (score[maxbp] < score[bp]) maxbp = bp;
        score[maxbp] = -9e9;

#define BreakPt(bp) ((bp)/2+lb)
#define BreakInX(bp) ((bp)%2==0)
        p = BreakPt(maxbp);
// Breaking in p
        if (BreakInX(maxbp)) {  // break in x
            n1 = n2 = 0;
            for (r=0; r<d; r++) {
                if (s[r] < p) {
                    s1[n1] = s[r];
                    y1[n1] = ys[r];
                    n1++;
                }
                else if (s[r] > p) {
                    s2[n2] = s[r]-p;
                    y2[n2] = ys[r];
                    n2++;
                }
                else { // if (s[r] == p)  i.e.,  r = si[p]
                    s1[n1] = p;  s2[n2] = 0;
                    if (r == d-1 || r == d-2) {
                        y1[n1] = y2[n2] = ys[r-1];
                        extral = ys[r] - ys[r-1];
                    }
                    if (r == 0 || r == 1) {
                        y1[n1] = y2[n2] = ys[r+1];
                        extral = ys[r+1] - ys[r];
                    }
                    else {
                        y1[n1] = y2[n2] = ys[r];
                        extral = 0;
                    }
                    n1++;  n2++;
                }
            }
            ll = extral + flutes_wl_LMD( flute_p, p+1, xs, y1, s1, newacc)
                + flutes_wl_LMD( flute_p, d-p, xs+p, y2, s2, newacc);
        }
        else {  // if (!BreakInX(maxbp))
            n1 = n2 = 0;
            for (r=0; r<d; r++) {
                if (si[r] < p) {
                    s1[si[r]] = n1;
                    x1[n1] = xs[r];
                    n1++;
                }
                else if (si[r] > p) {
                    s2[si[r]-p] = n2;
                    x2[n2] = xs[r];
                    n2++;
                }
                else { // if (si[r] == p)  i.e.,  r = s[p]
                    s1[p] = n1;  s2[0] = n2;
                    if (r == d-1 || r == d-2) {
                        x1[n1] = x2[n2] = xs[r-1];
                        extral = xs[r] - xs[r-1];
                    }
                    if (r == 0 || r == 1) {
                        x1[n1] = x2[n2] = xs[r+1];
                        extral = xs[r+1] - xs[r];
                    }
                    else {
                        x1[n1] = x2[n2] = xs[r];
                        extral = 0;
                    }
                    n1++;  n2++;
                }
            }
            ll = extral + flutes_wl_LMD( flute_p, p+1, x1, ys, s1, newacc) + 
	                  flutes_wl_LMD( flute_p, d-p, x2, ys+p, s2, newacc);
        }
        if (minl > ll) minl = ll;
    }
    return_val = minl ;

EXIT_FLUTE_WL_MD:
    UTDFREE( score ) ;
    UTDFREE( penalty ) ;
    UTDFREE( x1 ) ;
    UTDFREE( y1 ) ;
    UTDFREE( x2 ) ;
    UTDFREE( y2 ) ;
    UTDFREE( distx ) ;
    UTDFREE( disty ) ;
    UTDFREE( si ) ;
    UTDFREE( s1 ) ;
    UTDFREE( s2 ) ;

    return return_val ;
} /* end flutes_wl_MD() */

static int orderx(const void *a, const void *b)
{
    struct point *pa, *pb;

    pa = *(struct point**)a;
    pb = *(struct point**)b;

    if (pa->x < pb->x) return -1;
    if (pa->x > pb->x) return 1;
    return 0;
}

static int ordery(const void *a, const void *b)
{
    struct point *pa, *pb;

    pa = *(struct point**)a;
    pb = *(struct point**)b;

    if (pa->y < pb->y) return -1;
    if (pa->y > pb->y) return 1;
    return 0;
}

Tree flute( FluteState * flute_p, int d, FLUTE_DTYPE *x, FLUTE_DTYPE *y, int acc)
{
    FLUTE_DTYPE *xs, *ys, minval;
    int *s;
    int i, j, k, minidx;
    struct point *pt, **ptp, *tmpp;
    Tree t;
    
    if (d==2) {
        t.deg = 2;
        t.length = ADIFF(x[0], x[1]) + ADIFF(y[0], y[1]);
        t.branch = UTDMALLOC( 2, Branch ) ;
        t.branch[0].x = x[0];
        t.branch[0].y = y[0];
        t.branch[0].n = 1;
        t.branch[1].x = x[1];
        t.branch[1].y = y[1];
        t.branch[1].n = 1;
    }
    else {
        xs = UTDMALLOC( d, FLUTE_DTYPE ) ;
        ys = UTDMALLOC( d, FLUTE_DTYPE ) ;
        s = UTDMALLOC( d, int ) ;
        pt = UTDMALLOC( d+1, FLUTE_POINT ) ;
        ptp = UTDMALLOC( d+1, FLUTE_POINTPTR ) ;

        for (i=0; i<d; i++) {
            pt[i].x = x[i];
            pt[i].y = y[i];
            ptp[i] = &pt[i];
        }

        // sort x
        if (d<200) {
            for (i=0; i<d-1; i++) {
                minval = ptp[i]->x;
                minidx = i;
                for (j=i+1; j<d; j++) {
                    if (minval > ptp[j]->x) {
                        minval = ptp[j]->x;
                        minidx = j;
                    }
                }
                tmpp = ptp[i];
                ptp[i] = ptp[minidx];
                ptp[minidx] = tmpp;
            }
        } else {
            qsort(ptp, d, sizeof(struct point *), orderx);
        }

#if FLUTE_REMOVE_DUPLICATE_PIN==1
        ptp[d] = &pt[d];
        ptp[d]->x = ptp[d]->y = -999999;
        j = 0;
        for (i=0; i<d; i++) {
            for (k=i+1; ptp[k]->x == ptp[i]->x; k++)
                if (ptp[k]->y == ptp[i]->y)  // pins k and i are the same
                    break;
            if (ptp[k]->x != ptp[i]->x)
                ptp[j++] = ptp[i];
        }
        d = j;
#endif
        
        for (i=0; i<d; i++) {
            xs[i] = ptp[i]->x;
            ptp[i]->o = i;
        }

        // sort y to find s[]
        if (d<200) {
            for (i=0; i<d-1; i++) {
                minval = ptp[i]->y;
                minidx = i;
                for (j=i+1; j<d; j++) {
                    if (minval > ptp[j]->y) {
                        minval = ptp[j]->y;
                        minidx = j;
                    }
                }
                ys[i] = ptp[minidx]->y;
                s[i] = ptp[minidx]->o;
                ptp[minidx] = ptp[i];
            }
            ys[d-1] = ptp[d-1]->y;
            s[d-1] = ptp[d-1]->o;
        } else {
            qsort(ptp, d, sizeof(struct point *), ordery);
            for (i=0; i<d; i++) {
                ys[i] = ptp[i]->y;
                s[i] = ptp[i]->o;
            }
        }
        
        t = flutes( flute_p, d, xs, ys, s, acc);

        UTDFREE(xs);
        UTDFREE(ys);
        UTDFREE(s);
        UTDFREE(pt);
        UTDFREE(ptp);
    }

    return t;
} /* end flute() */

// xs[] and ys[] are coords in x and y in sorted order
// s[] is a list of nodes in increasing y direction
//   if nodes are indexed in the order of increasing x coord
//   i.e., s[i] = s_i as defined in paper
// The points are (xs[s[i]], ys[i]) for i=0..d-1
//             or (xs[i], ys[si[i]]) for i=0..d-1

Tree flutes_RDP( FluteState * flute_p, int d, FLUTE_DTYPE xs[], FLUTE_DTYPE ys[], int s[], int acc)
{
    int i, j, ss;
    
    for (i=0; i<d-1; i++) {
        if (xs[s[i]]==xs[s[i+1]] && ys[i]==ys[i+1]) {
            if (s[i] < s[i+1])
                ss = s[i+1];
            else {
                ss = s[i];
                s[i] = s[i+1];
            }
            for (j=i+2; j<d; j++) {
                ys[j-1] = ys[j];
                s[j-1] = s[j];
            }
            for (j=ss+1; j<d; j++)
                xs[j-1] = xs[j];
            for (j=0; j<=d-2; j++)
                if (s[j] > ss) s[j]--;
            i--;
            d--;
        }
    }
    return flutes_ALLD( flute_p, d, xs, ys, s, acc);
}
    
// For low-degree, i.e., 2 <= d <= D
Tree flutes_LD( FluteState * flute_p, int d, FLUTE_DTYPE *xs, FLUTE_DTYPE *ys, int *s)
{
      int k, pi, i, j, k0, k1 ;
      FLUTE_SOLNPTR rlist, bestrlist;
      FLUTE_DTYPE dd[2 * FLUTE_D - 2]; // 0..D-2 for v, D-1..2*D-3 for h
      FLUTE_DTYPE minl, sum;
      FLUTE_DTYPE l[MPOWV + 1];
      int hflip;
      Tree t;

      t.deg = d;
      t.branch = UTDMALLOC((2 * d - 2), Branch) ;
      if (d == 2) {
	      minl = xs[1] - xs[0] + ys[1] - ys[0];
	      t.branch[0].x = xs[s[0]];
	      t.branch[0].y = ys[0];
	      t.branch[0].n = 1;
	      t.branch[1].x = xs[s[1]];
	      t.branch[1].y = ys[1];
	      t.branch[1].n = 1;
      } else if (d == 3) {
	      minl = xs[2] - xs[0] + ys[2] - ys[0];
	      t.branch[0].x = xs[s[0]];
	      t.branch[0].y = ys[0];
	      t.branch[0].n = 3;
	      t.branch[1].x = xs[s[1]];
	      t.branch[1].y = ys[1];
	      t.branch[1].n = 3;
	      t.branch[2].x = xs[s[2]];
	      t.branch[2].y = ys[2];
	      t.branch[2].n = 3;
	      t.branch[3].x = xs[1];
	      t.branch[3].y = ys[1];
	      t.branch[3].n = 3;
      } else {
	      k = 0;
	      if (s[0] < s[2]) k++;
	      if (s[1] < s[2]) k++;

	      k0 = k ;
	      for (i = 3; i <= d - 1; i++) { // p0=0 always, skip i=1 for symmetry
		      pi = s[i];
		      for (j = d - 1; j > i; j--)
			      if (s[j] < s[i])
				      pi--;
		      k = pi + (i + 1) * k;
	      }

	      k1 = k ;
	      if (k < flute_p->numgrp[d]) { // no horizontal flip
		      hflip = 0;
		      for (i = 1; i <= d - 3; i++) {
			      dd[i] = ys[i + 1] - ys[i];
			      dd[d - 1 + i] = xs[i + 1] - xs[i];
		      }
	      } else {
		      hflip = 1;
		      k = 2 * flute_p->numgrp[d] - 1 - k;
		      for (i = 1; i <= d - 3; i++) {
			      dd[i] = ys[i + 1] - ys[i];
			      dd[d - 1 + i] = xs[d - 1 - i] - xs[d - 2 - i];
		      }
	      }

	      minl = l[0] = xs[d - 1] - xs[0] + ys[d - 1] - ys[0];
	      if( k < 0 ){
		fprintf( stderr, "mistake here.\n" ) ;
		k = 0 ;
	      }
	      rlist = flute_p->lut[d][k];
	      for (i = 0; rlist->seg[i] > 0; i++)
		      minl += dd[rlist->seg[i]];
	      bestrlist = rlist;
	      l[1] = minl;
	      j = 2;
	      while (j <= flute_p->numsoln[d][k]) {
		      rlist++;
		      sum = l[rlist->parent];
		      for (i = 0; rlist->seg[i] > 0; i++)
			      sum += dd[rlist->seg[i]];
		      for (i = 10; rlist->seg[i] > 0; i--)
			      sum -= dd[rlist->seg[i]];
		      if (sum < minl) {
			      minl = sum;
			      bestrlist = rlist;
		      }
		      l[j++] = sum;
	      }

	      t.branch[0].x = xs[s[0]];
	      t.branch[0].y = ys[0];
	      t.branch[1].x = xs[s[1]];
	      t.branch[1].y = ys[1];
	      for (i = 2; i < d - 2; i++) {
		      t.branch[i].x = xs[s[i]];
		      t.branch[i].y = ys[i];
		      t.branch[i].n = bestrlist->neighbor[i];
	      }
	      t.branch[d - 2].x = xs[s[d - 2]];
	      t.branch[d - 2].y = ys[d - 2];
	      t.branch[d - 1].x = xs[s[d - 1]];
	      t.branch[d - 1].y = ys[d - 1];
	      if (hflip) {
		      if (s[1] < s[0]) {
			      t.branch[0].n = bestrlist->neighbor[1];
			      t.branch[1].n = bestrlist->neighbor[0];
		      } else {
			      t.branch[0].n = bestrlist->neighbor[0];
			      t.branch[1].n = bestrlist->neighbor[1];
		      }
		      if (s[d - 1] < s[d - 2]) {
			      t.branch[d - 2].n = bestrlist->neighbor[d - 1];
			      t.branch[d - 1].n = bestrlist->neighbor[d - 2];
		      } else {
			      t.branch[d - 2].n = bestrlist->neighbor[d - 2];
			      t.branch[d - 1].n = bestrlist->neighbor[d - 1];
		      }
		      for (i = d; i < 2 * d - 2; i++) {
			      t.branch[i].x = xs[d - 1 - bestrlist->col[i - d]];
			      t.branch[i].y = ys[bestrlist->row[i - d]];
			      t.branch[i].n = bestrlist->neighbor[i];
		      }
	      } else { // !hflip
		      if (s[0] < s[1]) {
			      t.branch[0].n = bestrlist->neighbor[1];
			      t.branch[1].n = bestrlist->neighbor[0];
		      } else {
			      t.branch[0].n = bestrlist->neighbor[0];
			      t.branch[1].n = bestrlist->neighbor[1];
		      }
		      if (s[d - 2] < s[d - 1]) {
			      t.branch[d - 2].n = bestrlist->neighbor[d - 1];
			      t.branch[d - 1].n = bestrlist->neighbor[d - 2];
		      } else {
			      t.branch[d - 2].n = bestrlist->neighbor[d - 2];
			      t.branch[d - 1].n = bestrlist->neighbor[d - 1];
		      }
		      for (i = d; i < 2 * d - 2; i++) {
			      t.branch[i].x = xs[bestrlist->col[i - d]];
			      t.branch[i].y = ys[bestrlist->row[i - d]];
			      t.branch[i].n = bestrlist->neighbor[i];
		      }
	      }
	}
	t.length = minl;

	return t;

} /* end flutes_LD() */


// For medium-degree, i.e., D+1 <= d
Tree flutes_MD( FluteState * flute_p, int d, FLUTE_DTYPE xs[], FLUTE_DTYPE ys[], int s[], int acc)
{
    float *score, *penalty, pnlty, dx, dy;
    int ms, mins, maxs, minsi, maxsi;
    int i, r, p, maxbp, bestbp, bp, nbp, ub, lb, n1, n2, nn1, nn2, newacc;
    int *si, *s1, *s2, degree ;
    Tree t, t1, t2, bestt1, bestt2;
    FLUTE_DTYPE ll, minl, coord1, coord2;
    FLUTE_DTYPE *distx, *disty, xydiff;
    FLUTE_DTYPE *x1, *x2, *y1, *y2;

    degree = d + 1 ;
    score = UTDMALLOC( 2 * degree, float ) ;
    penalty = UTDMALLOC( degree, float ) ;

    x1 = UTDMALLOC( degree, FLUTE_DTYPE ) ;
    y1 = UTDMALLOC( degree, FLUTE_DTYPE ) ;
    x2 = UTDMALLOC( degree, FLUTE_DTYPE ) ;
    y2 = UTDMALLOC( degree, FLUTE_DTYPE ) ;
    distx = UTDMALLOC( degree, FLUTE_DTYPE ) ;
    disty = UTDMALLOC( degree, FLUTE_DTYPE ) ;
    si = UTDMALLOC( degree, int ) ;
    s1 = UTDMALLOC( degree, int ) ;
    s2 = UTDMALLOC( degree, int ) ;

    if (s[0] < s[d-1]) {
        ms = max(s[0], s[1]);
        for (i=2; i<=ms; i++)
            ms = max(ms, s[i]);
        if (ms <= d-3) {
	  for (i=0; i<=ms; i++) {
	    x1[i] = xs[i];
	    y1[i] = ys[i];
	    s1[i] = s[i];
	  }
	  x1[ms+1] = xs[ms]; 
	  y1[ms+1] = ys[ms]; 
	  s1[ms+1] = ms+1;
            
	  s2[0] = 0;
	  for (i=1; i<=d-1-ms; i++)
	    s2[i] = s[i+ms]-ms;

	  t1 = flutes_LMD( flute_p, ms+2, x1, y1, s1, acc);
	  t2 = flutes_LMD( flute_p, d-ms, xs+ms, ys+ms, s2, acc);
	  t = dmergetree(t1, t2);
	  UTDFREE(t1.branch);
	  UTDFREE(t2.branch);
            
	  UTDFREE( score ) ;
	  UTDFREE( penalty ) ;
	  UTDFREE( x1 ) ;
	  UTDFREE( y1 ) ;
	  UTDFREE( x2 ) ;
	  UTDFREE( y2 ) ;
	  UTDFREE( distx ) ;
	  UTDFREE( disty ) ;
	  UTDFREE( si ) ;
	  UTDFREE( s1 ) ;
	  UTDFREE( s2 ) ;

	  return t;
        }
    }
    else {  // (s[0] > s[d-1])
        ms = min(s[0], s[1]);
        for (i=2; i<=d-1-ms; i++)
            ms = min(ms, s[i]);
        if (ms >= 2) {
	  x1[0] = xs[ms];
	  y1[0] = ys[0];
	  s1[0] = s[0]-ms+1;
	  for (i=1; i<=d-1-ms; i++) {
	    x1[i] = xs[i+ms-1];
	    y1[i] = ys[i];
	    s1[i] = s[i]-ms+1;
	  }
	  x1[d-ms] = xs[d-1];
	  y1[d-ms] = ys[d-1-ms];
	  s1[d-ms] = 0;
            
	  s2[0] = ms;
	  for (i=1; i<=ms; i++)
	    s2[i] = s[i+d-1-ms];
            
	  t1 = flutes_LMD( flute_p, d+1-ms, x1, y1, s1, acc);
	  t2 = flutes_LMD( flute_p, ms+1, xs, ys+d-1-ms, s2, acc);
	  t = dmergetree(t1, t2);
	  UTDFREE(t1.branch);
	  UTDFREE(t2.branch);
            
	  UTDFREE( score ) ;
	  UTDFREE( penalty ) ;
	  UTDFREE( x1 ) ;
	  UTDFREE( y1 ) ;
	  UTDFREE( x2 ) ;
	  UTDFREE( y2 ) ;
	  UTDFREE( distx ) ;
	  UTDFREE( disty ) ;
	  UTDFREE( si ) ;
	  UTDFREE( s1 ) ;
	  UTDFREE( s2 ) ;

	  return t;
        }
    }

// Find inverse si[] of s[]
    for (r=0; r<d; r++)
        si[s[r]] = r;
    
// Determine breaking directions and positions dp[]
    lb=(d-2*acc+2)/4;
    if (lb < 2) lb = 2;
    ub=d-1-lb;

// Compute scores    
#define AA 0.6  // 2.0*BB
#define BB 0.3
    float CC = 7.4/((d+10.)*(d-3.));
    float DD = 4.8/(d-1);

    // Compute penalty[]    
    dx = CC*(xs[d-2]-xs[1]);
    dy = CC*(ys[d-2]-ys[1]);
    for (r = d/2, pnlty = 0; r>=2; r--, pnlty += dx)
        penalty[r] = pnlty,  penalty[d-1-r] = pnlty;
    penalty[1] = pnlty,  penalty[d-2] = pnlty;
    penalty[0] = pnlty,  penalty[d-1] = pnlty; 
    for (r = d/2-1, pnlty = dy; r>=2; r--, pnlty += dy)
        penalty[s[r]] += pnlty,  penalty[s[d-1-r]] += pnlty;
    penalty[s[1]] += pnlty,  penalty[s[d-2]] += pnlty;
    penalty[s[0]] += pnlty,  penalty[s[d-1]] += pnlty;
//#define CC 0.16
//#define v(r) ((r==0||r==1||r==d-2||r==d-1) ? d-3 : abs(d-1-r-r))
//    for (r=0; r<d; r++)
//        penalty[r] = v(r)*dx + v(si[r])*dy;

    // Compute distx[], disty[]
    xydiff = (xs[d-1] - xs[0]) - (ys[d-1] - ys[0]);
    if (s[0] < s[1])
        mins = s[0], maxs = s[1];
    else mins = s[1], maxs = s[0];
    if (si[0] < si[1])
        minsi = si[0], maxsi = si[1];
    else minsi = si[1], maxsi = si[0];
    for (r=2; r<=ub; r++) {
        if (s[r] < mins)
            mins = s[r];
        else if (s[r] > maxs)
            maxs = s[r];
        distx[r] = xs[maxs] - xs[mins];
        if (si[r] < minsi)
            minsi = si[r];
        else if (si[r] > maxsi)
            maxsi = si[r];
        disty[r] = ys[maxsi] - ys[minsi] + xydiff;
    }

    if (s[d-2] < s[d-1])
        mins = s[d-2], maxs = s[d-1];
    else mins = s[d-1], maxs = s[d-2];
    if (si[d-2] < si[d-1])
        minsi = si[d-2], maxsi = si[d-1];
    else minsi = si[d-1], maxsi = si[d-2];
    for (r=d-3; r>=lb; r--) {
        if (s[r] < mins)
            mins = s[r];
        else if (s[r] > maxs)
            maxs = s[r];
        distx[r] += xs[maxs] - xs[mins];
        if (si[r] < minsi)
            minsi = si[r];
        else if (si[r] > maxsi)
            maxsi = si[r];
        disty[r] += ys[maxsi] - ys[minsi];
    }

    nbp=0;
    for (r=lb; r<=ub; r++) {
        if (si[r]<=1)
            score[nbp] = (xs[r+1] - xs[r-1]) - penalty[r]
                - AA*(ys[2]-ys[1]) - DD*disty[r];
        else if (si[r]>=d-2)
            score[nbp] = (xs[r+1] - xs[r-1]) - penalty[r]
                - AA*(ys[d-2]-ys[d-3]) - DD*disty[r];
        else score[nbp] = (xs[r+1] - xs[r-1]) - penalty[r]
                 - BB*(ys[si[r]+1]-ys[si[r]-1]) - DD*disty[r];
        nbp++;
        
        if (s[r]<=1)
            score[nbp] = (ys[r+1] - ys[r-1]) - penalty[s[r]]
                - AA*(xs[2]-xs[1]) - DD*distx[r];
        else if (s[r]>=d-2)
            score[nbp] = (ys[r+1] - ys[r-1]) - penalty[s[r]]
                - AA*(xs[d-2]-xs[d-3]) - DD*distx[r];
        else score[nbp] = (ys[r+1] - ys[r-1]) - penalty[s[r]]
                 - BB*(xs[s[r]+1]-xs[s[r]-1]) - DD*distx[r];
        nbp++;
    }

    if (acc <= 3)
        newacc = 1;
    else {
        newacc = acc/2;
        if (acc >= nbp) acc = nbp-1;
    }
    
    minl = (FLUTE_DTYPE) INT_MAX;
    bestt1.branch = bestt2.branch = NULL;
    for (i=0; i<acc; i++) {
        maxbp = 0;
        for (bp=1; bp<nbp; bp++)
            if (score[maxbp] < score[bp]) maxbp = bp;
        score[maxbp] = -9e9;

#define BreakPt(bp) ((bp)/2+lb)
#define BreakInX(bp) ((bp)%2==0)
        p = BreakPt(maxbp);
// Breaking in p
        if (BreakInX(maxbp)) {  // break in x
            n1 = n2 = 0;
            for (r=0; r<d; r++) {
                if (s[r] < p) {
                    s1[n1] = s[r];
                    y1[n1] = ys[r];
                    n1++;
                }
                else if (s[r] > p) {
                    s2[n2] = s[r]-p;
                    y2[n2] = ys[r];
                    n2++;
                }
                else { // if (s[r] == p)  i.e.,  r = si[p]
                    s1[n1] = p;  s2[n2] = 0;
                    y1[n1] = y2[n2] = ys[r];
                    nn1 = n1;  nn2 = n2;
                    n1++;  n2++;
                }
            }

            t1 = flutes_LMD( flute_p, p+1, xs, y1, s1, newacc);
            t2 = flutes_LMD( flute_p, d-p, xs+p, y2, s2, newacc);
            ll = t1.length + t2.length;
            coord1 = t1.branch[t1.branch[nn1].n].y;
            coord2 = t2.branch[t2.branch[nn2].n].y;
            if (t2.branch[nn2].y > max(coord1, coord2))
                ll -= t2.branch[nn2].y - max(coord1, coord2);
            else if (t2.branch[nn2].y < min(coord1, coord2))
                ll -= min(coord1, coord2) - t2.branch[nn2].y;
        }
        else {  // if (!BreakInX(maxbp))
            n1 = n2 = 0;
            for (r=0; r<d; r++) {
                if (si[r] < p) {
                    s1[si[r]] = n1;
                    x1[n1] = xs[r];
                    n1++;
                }
                else if (si[r] > p) {
                    s2[si[r]-p] = n2;
                    x2[n2] = xs[r];
                    n2++;
                }
                else { // if (si[r] == p)  i.e.,  r = s[p]
                    s1[p] = n1;  s2[0] = n2;
                    x1[n1] = x2[n2] = xs[r];
                    n1++;  n2++;
                }
            }

            t1 = flutes_LMD( flute_p, p+1, x1, ys, s1, newacc);
            t2 = flutes_LMD( flute_p, d-p, x2, ys+p, s2, newacc);
            ll = t1.length + t2.length;
            coord1 = t1.branch[t1.branch[p].n].x;
            coord2 = t2.branch[t2.branch[0].n].x;
            if (t2.branch[0].x > max(coord1, coord2))
                ll -= t2.branch[0].x - max(coord1, coord2);
            else if (t2.branch[0].x < min(coord1, coord2))
                ll -= min(coord1, coord2) - t2.branch[0].x;
        }
        if (minl > ll) {
            minl = ll;
	    if( bestt1.branch ){
	      UTDFREE(bestt1.branch) ;
	    }
	    if( bestt2.branch ){
	      UTDFREE(bestt2.branch) ;
	    }
            bestt1 = t1;
            bestt2 = t2;
            bestbp = maxbp;
        }
        else {
	  if( t1.branch ){
            UTDFREE(t1.branch) ;
	  }
	  if( t2.branch ){
            UTDFREE(t2.branch) ;
	  }
        }
    }

#if FLUTE_LOCAL_REFINEMENT==1
    if (BreakInX(bestbp)) {
        t = hmergetree(bestt1, bestt2, s);
        local_refinement( d, &t, si[BreakPt(bestbp)]);
    } else {
        t = vmergetree(bestt1, bestt2);
        local_refinement( d, &t, BreakPt(bestbp));
    }
#else
    if (BreakInX(bestbp)) {
        t = hmergetree(bestt1, bestt2, s);
    } else {
        t = vmergetree(bestt1, bestt2);
    }
#endif
    
    UTDFREE(bestt1.branch);
    UTDFREE(bestt2.branch);

EXIT_FLUTES_MD:
    UTDFREE( score ) ;
    UTDFREE( penalty ) ;
    UTDFREE( x1 ) ;
    UTDFREE( y1 ) ;
    UTDFREE( x2 ) ;
    UTDFREE( y2 ) ;
    UTDFREE( distx ) ;
    UTDFREE( disty ) ;
    UTDFREE( si ) ;
    UTDFREE( s1 ) ;
    UTDFREE( s2 ) ;

    return t;
} /* end flutes_MD() */

static Tree dmergetree(Tree t1, Tree t2)
{
    int i, d, prev, curr, next, offset1, offset2;
    Tree t;

    t.deg = d = t1.deg + t2.deg - 2;
    t.length = t1.length + t2.length;
    t.branch = UTDMALLOC( (2*d-2), Branch ) ;
    offset1 = t2.deg-2;
    offset2 = 2*t1.deg-4;
    
    for (i=0; i<=t1.deg-2; i++) {
        t.branch[i].x = t1.branch[i].x;
        t.branch[i].y = t1.branch[i].y;
        t.branch[i].n = t1.branch[i].n + offset1;
    }
    for (i=t1.deg-1; i<=d-1; i++) {
        t.branch[i].x = t2.branch[i-t1.deg+2].x;
        t.branch[i].y = t2.branch[i-t1.deg+2].y;
        t.branch[i].n = t2.branch[i-t1.deg+2].n + offset2;
    }
    for (i=d; i<=d+t1.deg-3; i++) {
        t.branch[i].x = t1.branch[i-offset1].x;
        t.branch[i].y = t1.branch[i-offset1].y;
        t.branch[i].n = t1.branch[i-offset1].n + offset1;
    }
    for (i=d+t1.deg-2; i<=2*d-3; i++) {
        t.branch[i].x = t2.branch[i-offset2].x;
        t.branch[i].y = t2.branch[i-offset2].y;
        t.branch[i].n = t2.branch[i-offset2].n + offset2;
    }

    prev = t2.branch[0].n + offset2;
    curr = t1.branch[t1.deg-1].n + offset1;
    next = t.branch[curr].n;
    while (curr != next) {
        t.branch[curr].n = prev;
        prev = curr;
        curr = next;
        next = t.branch[curr].n;
    }
    t.branch[curr].n = prev;

    return t;
}

static Tree hmergetree(Tree t1, Tree t2, int *s)
{
    int i, prev, curr, next, extra, offset1, offset2;
    int p, ii, n1, n2, nn1, nn2;
    FLUTE_DTYPE coord1, coord2;
    Tree t;

    t.deg = t1.deg + t2.deg - 1;
    t.length = t1.length + t2.length;
    t.branch = UTDMALLOC( (2*t.deg-2), Branch ) ;
    offset1 = t2.deg-1;
    offset2 = 2*t1.deg-3;

    p = t1.deg - 1;
    n1 = n2 = 0;
    for (i=0; i<t.deg; i++) {
        if (s[i] < p) {
            t.branch[i].x = t1.branch[n1].x;
            t.branch[i].y = t1.branch[n1].y;
            t.branch[i].n = t1.branch[n1].n + offset1;
            n1++;            
        }
        else if (s[i] > p) {
            t.branch[i].x = t2.branch[n2].x;
            t.branch[i].y = t2.branch[n2].y;
            t.branch[i].n = t2.branch[n2].n + offset2;
            n2++;            
        }
        else {
            t.branch[i].x = t2.branch[n2].x;
            t.branch[i].y = t2.branch[n2].y;
            t.branch[i].n = t2.branch[n2].n + offset2;
            nn1 = n1;  nn2 = n2;  ii = i;
            n1++;  n2++;
        }
    }
    for (i=t.deg; i<=t.deg+t1.deg-3; i++) {
        t.branch[i].x = t1.branch[i-offset1].x;
        t.branch[i].y = t1.branch[i-offset1].y;
        t.branch[i].n = t1.branch[i-offset1].n + offset1;
    }
    for (i=t.deg+t1.deg-2; i<=2*t.deg-4; i++) {
        t.branch[i].x = t2.branch[i-offset2].x;
        t.branch[i].y = t2.branch[i-offset2].y;
        t.branch[i].n = t2.branch[i-offset2].n + offset2;
    }
    extra = 2*t.deg-3;
    coord1 = t1.branch[t1.branch[nn1].n].y;
    coord2 = t2.branch[t2.branch[nn2].n].y;
    if (t2.branch[nn2].y > max(coord1, coord2)) {
        t.branch[extra].y = max(coord1, coord2);
        t.length -= t2.branch[nn2].y - t.branch[extra].y;
    }  
    else if (t2.branch[nn2].y < min(coord1, coord2)) {
        t.branch[extra].y = min(coord1, coord2);
        t.length -= t.branch[extra].y - t2.branch[nn2].y;
    }
    else t.branch[extra].y = t2.branch[nn2].y;
    t.branch[extra].x = t2.branch[nn2].x;
    t.branch[extra].n = t.branch[ii].n;
    t.branch[ii].n = extra;
        
    prev = extra;
    curr = t1.branch[nn1].n + offset1;
    next = t.branch[curr].n;
    while (curr != next) {
        t.branch[curr].n = prev;
        prev = curr;
        curr = next;
        next = t.branch[curr].n;
    }
    t.branch[curr].n = prev;
    
    return t;
}

static Tree vmergetree(Tree t1, Tree t2)
{
    int i, prev, curr, next, extra, offset1, offset2;
    FLUTE_DTYPE coord1, coord2;
    Tree t;

    t.deg = t1.deg + t2.deg - 1;
    t.length = t1.length + t2.length;
    t.branch = UTDMALLOC( (2*t.deg-2), Branch ) ;
    offset1 = t2.deg-1;
    offset2 = 2*t1.deg-3;

    for (i=0; i<=t1.deg-2; i++) {
        t.branch[i].x = t1.branch[i].x;
        t.branch[i].y = t1.branch[i].y;
        t.branch[i].n = t1.branch[i].n + offset1;
    }
    for (i=t1.deg-1; i<=t.deg-1; i++) {
        t.branch[i].x = t2.branch[i-t1.deg+1].x;
        t.branch[i].y = t2.branch[i-t1.deg+1].y;
        t.branch[i].n = t2.branch[i-t1.deg+1].n + offset2;
    }
    for (i=t.deg; i<=t.deg+t1.deg-3; i++) {
        t.branch[i].x = t1.branch[i-offset1].x;
        t.branch[i].y = t1.branch[i-offset1].y;
        t.branch[i].n = t1.branch[i-offset1].n + offset1;
    }
    for (i=t.deg+t1.deg-2; i<=2*t.deg-4; i++) {
        t.branch[i].x = t2.branch[i-offset2].x;
        t.branch[i].y = t2.branch[i-offset2].y;
        t.branch[i].n = t2.branch[i-offset2].n + offset2;
    }
    extra = 2*t.deg-3;
    coord1 = t1.branch[t1.branch[t1.deg-1].n].x;
    coord2 = t2.branch[t2.branch[0].n].x;
    if (t2.branch[0].x > max(coord1, coord2)) {
        t.branch[extra].x = max(coord1, coord2);
        t.length -= t2.branch[0].x - t.branch[extra].x;
    }        
    else if (t2.branch[0].x < min(coord1, coord2)) {
        t.branch[extra].x = min(coord1, coord2);
        t.length -= t.branch[extra].x - t2.branch[0].x;
    }
    else t.branch[extra].x = t2.branch[0].x;
    t.branch[extra].y = t2.branch[0].y;
    t.branch[extra].n = t.branch[t1.deg-1].n;
    t.branch[t1.deg-1].n = extra;

    prev = extra;
    curr = t1.branch[t1.deg-1].n + offset1;
    next = t.branch[curr].n;
    while (curr != next) {
        t.branch[curr].n = prev;
        prev = curr;
        curr = next;
        next = t.branch[curr].n;
    }
    t.branch[curr].n = prev;

    return t;
}

static void local_refinement( FluteState * flute_p, int deg, Tree *tp, int p)
{
    int d, dd, i, ii, j, prev, curr, next, root;
    int *SteinerPin, *index, *ss, degree ;
    FLUTE_DTYPE *x, *xs, *ys ;
    Tree tt;

    degree = deg + 1 ;
    SteinerPin = UTDMALLOC( 2 * degree, int ) ;
    index = UTDMALLOC( 2 * degree, int ) ;
    x = UTDMALLOC( degree, FLUTE_DTYPE ) ;
    xs = UTDMALLOC( degree, FLUTE_DTYPE ) ;
    ys = UTDMALLOC( degree, FLUTE_DTYPE ) ;
    ss = UTDMALLOC( degree, int ) ;
 
    d = tp->deg;
    root = tp->branch[p].n;

// Reverse edges to point to root    
    prev = root;
    curr = tp->branch[prev].n;
    next = tp->branch[curr].n;
    while (curr != next) {
        tp->branch[curr].n = prev;
        prev = curr;
        curr = next;
        next = tp->branch[curr].n;
    }
    tp->branch[curr].n = prev;
    tp->branch[root].n = root;
    
// Find Steiner nodes that are at pins
    for (i=d; i<=2*d-3; i++)
        SteinerPin[i] = -1;
    for (i=0; i<d; i++) {
        next = tp->branch[i].n;
        if (tp->branch[i].x == tp->branch[next].x &&
            tp->branch[i].y == tp->branch[next].y)
            SteinerPin[next] = i; // Steiner 'next' at Pin 'i'
    }
    SteinerPin[root] = p;

// Find pins that are directly connected to root    
    dd = 0;
    for (i=0; i<d; i++) {
        curr = tp->branch[i].n;
        if (SteinerPin[curr] == i)
            curr = tp->branch[curr].n;
        while (SteinerPin[curr] < 0)
            curr = tp->branch[curr].n;
        if (curr == root) {
            x[dd] = tp->branch[i].x;
            if (SteinerPin[tp->branch[i].n] == i && tp->branch[i].n != root)
                index[dd++] = tp->branch[i].n;  // Steiner node
            else index[dd++] = i;  // Pin
        }
    }

    if (4 <= dd && dd <= FLUTE_D) {
// Find Steiner nodes that are directly connected to root    
        ii=dd;
        for (i=0; i<dd; i++) {
            curr = tp->branch[index[i]].n;
            while (SteinerPin[curr] < 0) {
                index[ii++] = curr;
                SteinerPin[curr] = INT_MAX;
                curr = tp->branch[curr].n;
            }
        }
        index[ii] = root;
        
        for (ii=0; ii<dd; ii++) {
            ss[ii] = 0;
            for (j=0; j<ii; j++)
                if (x[j] < x[ii])
                    ss[ii]++;
            for (j=ii+1; j<dd; j++)
                if (x[j] <= x[ii])
                    ss[ii]++;
            xs[ss[ii]] = x[ii];
            ys[ii] = tp->branch[index[ii]].y;
        }

        tt = flutes_LD( flute_p, dd, xs, ys, ss);

// Find new wirelength
        tp->length += tt.length;
        for (ii=0; ii<2*dd-3; ii++) {
            i = index[ii];
            j = tp->branch[i].n;
            tp->length -= ADIFF(tp->branch[i].x, tp->branch[j].x)
                + ADIFF(tp->branch[i].y, tp->branch[j].y);
        }
        
// Copy tt into t
        for (ii=0; ii<dd; ii++) {
            tp->branch[index[ii]].n = index[tt.branch[ii].n];
        }
        for (; ii<=2*dd-3; ii++) {
            tp->branch[index[ii]].x = tt.branch[ii].x;
            tp->branch[index[ii]].y = tt.branch[ii].y;
            tp->branch[index[ii]].n = index[tt.branch[ii].n];
        }
        UTDFREE(tt.branch);
    }
    
    UTDFREE( SteinerPin ) ;
    UTDFREE( index ) ;
    UTDFREE( x ) ;
    UTDFREE( xs ) ;
    UTDFREE( ys ) ;
    UTDFREE( ss ) ;
    
    return;
} /* end local_refinement() */


FLUTE_DTYPE flute_wirelength(Tree t)
{
    int i, j;
    FLUTE_DTYPE l=0;

    for (i=0; i<2*t.deg-2; i++) {
        j = t.branch[i].n;
        l += ADIFF(t.branch[i].x, t.branch[j].x)
            + ADIFF(t.branch[i].y, t.branch[j].y);
    }

    return l;
}

void printtree(Tree t)
{
    int i;

    for (i=0; i<t.deg; i++)
        printf(" %-2d:  x=%4g  y=%4g  e=%d\n",
               i, (float) t.branch[i].x, (float) t.branch[i].y, t.branch[i].n);
    for (i=t.deg; i<2*t.deg-2; i++)
        printf("s%-2d:  x=%4g  y=%4g  e=%d\n",
               i, (float) t.branch[i].x, (float) t.branch[i].y, t.branch[i].n);
    printf("\n");
}

// Output in a format that can be plotted by gnuplot
void plottree(Tree t)
{
    int i;

    for (i=0; i<2*t.deg-2; i++) {
        printf("%d %d\n", t.branch[i].x, t.branch[i].y);
        printf("%d %d\n\n", t.branch[t.branch[i].n].x,
               t.branch[t.branch[i].n].y);
    }
}
} // end namespace