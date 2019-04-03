//////////////////////////////////////////////////////////////////
/*
 NOTES:

 1. Decide on valuable info to output
 2. In the R calling function check for integrity, produce a warning
    if a parent has no record, but do add the record. Just send the pedigree
    to a function the simply add those and produce a list of missing records
    of parents
*/
//////////////////////////////////////////////////////////////////



//
// 5 / 12 / 2005. pedSort - was prepped_2.c ..
//
// 5 / 19 / 2005:
// Super Fast Algorithm ... Putting nonparents at the bottom and sorting
// only the parent portion of the pedigree saved about 75% of the time!!!
//
// 5 / 23 / 2005. Call from R ..
//
// 11 / 11 / 2008:
// Included in package GEP
//
// 07/22/2010:
// Create independent pedigree data management package, PDMP


#include <R.h>
#define MaxSwap 100

/* swap(&x[c], &y[d]) is how to call swap to interchange
    a[c] <==> a[d].
    You can use it with variables too,i.e., swap(&x, &y) */
void swap(int *x, int *y) {
  int temp;
  if(x == y) return;
  temp = *x;
  *x = *y;
  *y = temp;
}

int OrderInIVec(int num, int *a, int length) {
  // Order In Integer Vector: returns the order or a 0 if num is NOT in a[]
  // num must not equal 0, a[] must not contain 0's!

  int i;
  if(!num) return 0;
  for(i=1; i<=length; i++) if(a[i] == num) return i;
  return 0;
}

int OrderInDexVec(int ind, int *id, int at, int N, int *dex) {
  // Return actual order in a sorting index vector
  // Do the search only in [at:N]

  int i;
  if(!ind) return 0;
  for(i=at; i<=N; i++)
    if(ind == id[dex[i]]) return i;
  return 0;
}

void niRunThrough(int r, int *id, int *s, int *d, int N, int *dex) {
  // Slightly slower than RunThrough
  // Takes 104% of the time, based on 6 runs of 9391 records

  int x, j;
  do {
    x = 0;
    while((j= OrderInDexVec(s[dex[r]], id, r+1, N, dex)) > r) {
      x++;
      swap(&dex[r], &dex[j]);
    }
    while((j= OrderInDexVec(d[dex[r]], id, r+1, N, dex)) > r) {
      x++;
      swap(&dex[r], &dex[j]);
    }
  } while( x );

}


void RunThrough(int r, int *id, int *s, int *d, int N, int *dex, int nnp, int *NCalled, int *pedRow) {
    // Start at dex[r] and search down the pedigree ...
    // The assumption here is that if a parent is not down then s/he is
    //   -for sure- up without further verification ..

    int Ord;
    int ncalled = (*NCalled);

    ncalled++;
    if(ncalled > MaxSwap) {
	Rprintf("\n ... Looping Pedigree Problem ... \n");
	Rprintf("\nIndividual: %d  With parents: %d   and   %d\n", pedRow[0], pedRow[1], pedRow[2]);
	return;
    }

    if( (Ord=OrderInDexVec(s[dex[r]], id, r+1, N-nnp, dex)) > r ) {
	swap(&dex[r], &dex[Ord]);
	RunThrough(r, id, s, d, N, dex, nnp, &ncalled, pedRow);
    }

    if( (Ord=OrderInDexVec(d[dex[r]], id, r+1, N-nnp, dex)) > r ) {
	swap(&dex[r], &dex[Ord]);
	RunThrough(r, id, s, d, N, dex, nnp, &ncalled, pedRow);
    }
}


void pedSort_Base(int *sr, int *dm, int *dex, int N, int *BaseNo) {
  int i;
  *BaseNo = 0;
  for(i=1; i<=N; i++)
    if(!sr[dex[i]] && !dm[dex[i]]) {
      (*BaseNo)++;
      swap(&dex[i], &dex[(*BaseNo)]);
    }
}

void pedSort_NonParents(int *ind, int *dex, int N, int BaseNo, int *unp, int nnp, int *BaseOverlap) {
  int i, Ord;
  for(i=(BaseNo+1); i<=N; i++) {
    while((Ord=OrderInIVec(ind[dex[i]], unp, nnp))) {
      swap(&dex[i], &dex[N]);
      N--;

      swap(&unp[nnp], &unp[Ord]); // Don't re-consider unp[Ord]!
      nnp--;
    }
  }
  *BaseOverlap = nnp;
}

void pedSortMain(int *ind, int *sr, int *dm, int *unp, int *dex,
		int *nrow_ped, int *length_np, int *nrow_base) {

  int i, n, nnp, BaseNo, BaseOverlap;
  int pedRow[3], NC=0;

  n = *nrow_ped;
  nnp = *length_np;

  pedSort_Base(sr, dm, dex, n, &BaseNo); // Step-1: Put Base Individuals on Top.
  *nrow_base = BaseNo;

  // Step-2: Put Non-Parents down at the end.
  pedSort_NonParents(ind, dex, n, BaseNo, unp, nnp, &BaseOverlap);
  nnp -= BaseOverlap;

  for(i=BaseNo+1; i<=(n-nnp); i++) { // Loop over parents ONLY --
      pedRow[0] = ind[dex[i]]; pedRow[1] = sr[dex[i]]; pedRow[2] =  dm[dex[i]];
      RunThrough(i, ind, sr, dm, n, dex, nnp, &NC, pedRow);
  }

  // for(i=1; i<=n; i++)
  //  fprintf(fpo,"%d %d %d\n", ind[dex[i]], sr[dex[i]], dm[dex[i]]);


  Rprintf("  pedSortMain: Base=  %d, Non-Parents=  %d, Base-NP Overlap= %d\n",BaseNo, nnp, BaseOverlap);

}
