#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>

unsigned long extendedAlg(mpz_t a, mpz_t b);

void main(){
  // take the first parameter (n)

  // while n is not 0
    // take 4 integers (a, p, Gx, Gy)

    // create a result var (Rx, Ry) passing the same values as Gx and Gy

    // (n-1) times do
      // if Gx == Rx && Gy == Ry
        // set delta equal to (3Xq^2 + a) / (2Yq) mod p
      // else
        // set delta equal to (Yg – Yq) / (Xg-Xq) mod p

      // calculate the new Rx = ( delta * 2 - Xq - Xg ) mod p
      // calculate the new Yx = ( delta * (Xq – Xr) - Yq ) mod p

    // print the result (Rx, Ry)

    // take the new first parameter (n)
}

unsigned long extendedAlg(int a, int b){
  mpz_t x;
  mpz_t y;
  mpz_t quo;
  mpz_t res;
  mpz_t m_ant;
  mpz_t n_ant;
  mpz_t m;
  mpz_t n;
  mpz_t aux;

  mpz_init(quo);
  mpz_init(res);
  mpz_init(aux);
  mpz_init(n_ant);
  mpz_init(m);

  mpz_init_set_ui(x, a);
  mpz_init_set_ui(y, b);

  mpz_init_set_ui(m_ant, 1);
  mpz_init_set_ui(n, 1);

  while(mpz_cmp_ui(y, 1)) {
    // divide x por y, guardando o quociente e o resto
    mpz_tdiv_qr(quo, res, x, y);

    // aux = m_ant - m * quo
    mpz_mul(aux, m, quo);
    mpz_sub(aux, m_ant, aux);

    // m_ant = m
    // m = aux
    mpz_set(m_ant, m);
    mpz_set(m, aux);

    // n = n_ant - n * quo
    mpz_mul(aux, n, quo);
    mpz_sub(aux, n_ant, aux);

    // n_ant = n
    // n = aux
    mpz_set(n_ant, n);
    mpz_set(n, aux);

    // x = y
    mpz_set(x, y);

    // y = res
    mpz_set(y, res);
  }

  if(mpz_cmp_ui(n, 0) < 0){
    mpz_add(n, n, a);
  }

  mpz_clear(x);
  mpz_clear(y);
  mpz_clear(quo);
  mpz_clear(res);
  mpz_clear(aux);
  mpz_clear(m_ant);
  mpz_clear(n_ant);
  mpz_clear(m);
  mpz_clear(n);

  return mpz_get_ui(n);
}

