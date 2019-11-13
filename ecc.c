#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>

unsigned long extendedAlg(int a, int b);

int main(){
  int n, a, p, gx, gy, rx, ry, delta;
  unsigned long inverse;

  // take the first parameter (n)
  scanf("%d", &n);

  // while n is not 0
  while(n != 0){
    // take 4 integers (a, p, Gx, Gy)
    scanf("%d %d %d %d", &a, &p, &gx, &gy);

    // create a result var (Rx, Ry) passing the same values as Gx and Gy
    rx = gx;
    ry = gy;

    n--;
    // (n-1) times do
    while (n != 0){
      if(gx == rx && gy == ry){
        // set delta equal to (3 * Gx^2 + a) / (2 * Gy) mod p
        inverse = extendedAlg(p, (2 * gy));
        delta = ((3 * (gx * gx) + a) % p) * inverse;
        delta = delta % p;
      } else {
        // set delta equal to (Ry – Gy) / (Rx - Gx) mod p
        inverse = extendedAlg(p, (rx - gx));
        delta = ((ry - gy) % p) * inverse;
        delta = delta % p;
      }

      // calculate the new Rx = ( delta * 2 - Gx - Rx ) mod p
      rx = (delta * 2 - gx - rx) % p;
      // calculate the new Yx = ( delta * (Gx – Rx) - Gy ) mod p
      ry = (delta * (gx - rx) - gy) % p;

      n--;
    }
    // print the result (Rx, Ry)
    printf("%d %d\n", rx, ry);

    // take the new first parameter (n)
    scanf("%d", &n);
  }

  return 0;
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
    mpz_add_ui(n, n, a);
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

