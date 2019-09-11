#ifndef _RANDOM_H_
#define _RANDOM_H_

#define RAND_MAX 0x7fff




void advance_random ();
void warmup_random (double seed);
void randomize ();
double randomperc ();
int rnd (int low, int high);
double rndreal (double low, double high);
void random_permutation (int *perm, int size);

#endif