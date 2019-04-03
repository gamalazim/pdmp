
// by row
// starts at [1][1]
int v2m(int r, int c, int ORD) {
  int dex;
  if(r == 0 || c == 0) dex = 0;
  else dex = (ORD*(r-1)+c);
  return (dex);
}


// Build NRM A
void KinshipMain(int *sr, int *dm, int *RORD, double *A) {
  int r, c, idex, idexL, sdex, ddex, pdex;
  int ORD = (*RORD);
  for(r=1; r<=ORD; r++) {
    idex = v2m(r, r, ORD); 
    if(sr[r] > 0 && dm[r] > 0) {
      pdex = v2m(sr[r], dm[r], ORD);
      A[idex] = 1 + .5*A[pdex];
    }
    else A[idex] = 1;
       
    for(c=(r+1); c<=ORD; c++) {
      idex = v2m(r, c, ORD); idexL = v2m(c, r, ORD);
      sdex = v2m(r, sr[c], ORD); ddex = v2m(r, dm[c], ORD);
      A[idex] = A[idex] + .5*(A[sdex] + A[ddex]);
      A[idexL] = A[idex];
    }
  }      
}

// Build the L factor of A = LDL'
void getLMain(int *sr, int *dm, int *RORD, double *L) {
  int r, c, idex, sdex, ddex;
  int ORD = (*RORD);
  for(r=1; r<=ORD; r++) {
    idex = v2m(r, r, ORD);
    L[idex] = 1;

    for(c=1; c<r; c++) {
      idex = v2m(r, c, ORD);
      sdex = v2m(sr[r], c, ORD); ddex = v2m(dm[r], c, ORD);
      L[idex] = .5*(L[sdex] + L[ddex]);
    }
  }
}

// Build NRM A
//void KininvMain(int *sr, int *dm, int *RORD, double *A) {