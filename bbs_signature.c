/**********************************************************************
 *                                                                    *
 *   Create signatures using pairings.  Use paper "Revisiting BBS     *
 *   Signatures" by Stefano Tessaro and Chenzhi Zhufor.               *
 *   Example is for GF(p^11) found from pairing_sweep_alpha and       *
 *   get_curve.                                                       *
 *                                                                    *
 *********************************************************************/

#include "signature.h"
#include <string.h>

#define EMBED_DEGREE 11
#define NUMESG  16

/* create a point based on value i */

void msgpoint(POINT *H, long i, SIG_SYSTEM sig)
{
  unsigned char *buf;
  long j, *bfptr;

  buf = (unsigned char*)malloc(sizeof(long)*2048);
  for(j=0; j<2048; j++)
  {
    bfptr = (long*)(&buf[4*j]);
    *bfptr = i + 1;
  }
  hash0(H, sig, buf, 8192);
  free(buf);
}

int main(int argc, char *argv[])
{
  FILE *key;
  int i, j, k, m;
  SIG_SYSTEM sig;
  POLY_POINT PK[4];
  mpz_t  msghsh[NUMESG], e, xe, ex, sk[4];
  unsigned short *msg;
  POINT *Cgrp, Csum;
  POINT *Hvec, A;
  POLY_POINT A2, C2, eV, R;
  char filename[256];
  POLY wa, wc;
  
/*  read in system parameters from previously generated file  */

  get_system("curve_11_parameters.bin", &sig);
  minit(sig.prime);
  poly_irrd_set(sig.irrd);
  poly_mulprep(sig.irrd);

/*  read in public and private keys for testing signature aggregation  */
  
  key = fopen("key_data_11.skpk", "r");
  if(!key)
  {
    printf("can't find file key_data_11.skpk\n");
    exit(-2);
  }
  fread(&k, sizeof(int), 1, key);

/* only need one private/public key for BBS test */
  
  for(i=0; i<2; i++)
  {
    poly_point_init(&PK[i]);
    mpz_init(sk[i]);
    mpz_inp_raw(sk[i], key);
    poly_point_read(&PK[i], key);
  }
  fclose(key);

/* read in random data (from random.org) for messages to sign  */

  msg = (unsigned short*)malloc(sizeof(unsigned short)*(128*NUMESG+2));
  for(j=0; j<NUMESG; j++)
  {
    sprintf(filename, "message%02d.dat", j);
    key = fopen(filename, "r");
    if(!key)
    {
      printf("can't find %s\n", filename);
      exit(-4);
    }
    for(i=0; i<128; i++)
    fscanf(key, "%hd", &msg[j*128 + i]);
    fclose(key);
  }

/* create H[j] vector of points for each message j 
   and compute message signature with point. */
  
  Hvec = (POINT*)malloc(sizeof(POINT)*NUMESG);
  Cgrp = (POINT*)malloc(sizeof(POINT)*NUMESG);
  for(j=0; j<NUMESG; j++)
  {
    point_init(&Hvec[j]);
    msgpoint(&Hvec[j], j, sig);   // create H_j point
    point_init(&Cgrp[j]);
    mpz_init(msghsh[j]);
    hash1(msghsh[j], (unsigned char*)(&msg[j*128]), 256, sig.tor);  // get hash of message j
    elptic_mul(&Cgrp[j], Hvec[j], msghsh[j], sig.E); // compute m_j*H_j
  }

/* compute C, the combination of all signature points and G_1 */

  point_init(&Csum);
  point_copy(&Csum, sig.G1);
  for(j=0; j<NUMESG; j++)
    elptic_sum(&Csum, Csum, Cgrp[j], sig.E);

/* compute aggregate signature of all documents.
   access to a hardware random number generator 
   is a better choice than this example. */

  mpz_inits(e, xe, ex, NULL);
  mrand(e);                      // technically wrong modulus - don't do this for real
  mod_add(xe, e, sk[0], sig.tor);
  mpz_invert(ex, xe, sig.tor);
  point_init(&A);
  elptic_mul(&A, Csum, ex, sig.E);

/* Output final signature (A, e)  */
  
  gmp_printf("Signature random value is %Zd\n", e);
  point_printf("Signature point: ", A);
  printf("\n");
  poly_point_printf("public key:\n", PK[0]);

/* compute verification of signature */

  poly_init(&wa);
  poly_init(&wc);            // create check results
  poly_point_init(&A2);
  poly_point_init(&C2);      // space for G2 versions of G1 points
  tog2(&A2, A);
  tog2(&C2, Csum);
  poly_point_init(&eV);      // compute e * G2 + public key
  poly_elptic_mul(&eV, sig.G2, e, sig.Ex);
  poly_elptic_sum(&eV, eV, PK[0], sig.Ex);
  poly_point_init(&R);
  poly_point_rand(&R, sig.Ex);
  weil(&wa, A2, eV, R, sig.tor, sig.Ex);
  weil(&wc, C2, sig.G2, R, sig.tor, sig.Ex);
  if(poly_cmp(wa, wc))
    printf("Signature verifies!\n");
  else
    printf("Signature FAILS!\n");
  
  free(msg);
  for(j=0; j<NUMESG; j++)
  {
    point_clear(&Hvec[j]);
    point_clear(&Cgrp[j]);
    mpz_clear(msghsh[j]);
  }
  poly_point_clear(&PK[0]);
  poly_point_clear(&PK[1]);
  mpz_clear(sk[0]);
  mpz_clear(sk[1]);
  mpz_clears(e, xe, ex, NULL);
  point_clear(&A);
  poly_clear(&wa);
  poly_clear(&wc);
  poly_point_clear(&A2);
  poly_point_clear(&C2);
  poly_point_clear(&R);
  poly_point_clear(&eV);
}
