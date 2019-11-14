#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>

unsigned long extendedAlg(int a, int b);
int modulo(int a, int b);
int is_inverse(int qx, int qy, int gx, int gy, int p);

int main(){
  int n, a, p, gx, gy, qx, qy, rx, ry, delta;
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
      qx = rx;
      qy = ry;

      printf("(%d, %d) + (%d, %d)\n", qx, qy, gx, gy);

      if (qx == 0 && qy == 0) {
        rx = gx;
        ry = gy;
      } else if(is_inverse(qx, qy, gx, gy, p)){
        printf("(%d, %d) é inverso de (%d, %d)\n", gx, gy, qx, qy);
        rx = 0;
        ry = 0;
      } else {
        if(gx == qx && gy == qy){ // if P == Q
          // set delta equal to (3 * Gx^2 + a) / (2 * Gy) mod p
          inverse = extendedAlg(p, modulo(2 * gy, p));
          delta = modulo((3 * (gx * gx) + a), p) * inverse;
          delta = modulo(delta, p);
          printf("Delta: (3 * %d^2 + %d) / (2 * %d) mod %d = ", gx, a, gy, p);
          printf("%d / %d mod %d = ", modulo((3 * (gx * gx) + a), p), modulo(2 * gy, p), p);
          printf("%d * %lu mod %d = %d\n", modulo((3 * (gx * gx) + a), p), inverse, p, delta);
        } else {
          // set delta equal to (Gy – Qy) / (Gx - Qx) mod p
          inverse = extendedAlg(p, modulo(gx - qx, p));
          delta = modulo((gy - qy), p) * inverse;
          delta = modulo(delta, p);
          printf("Delta: (%d - %d) / (%d - %d) mod %d = ", gy, qy, gx, qx, p);
          printf("%d / %d mod %d = ", modulo((gy - qy), p), modulo(gx - qx, p), p);
          printf("%d * %lu mod %d = %d\n", modulo((gy - qy), p), inverse, p, delta);
        }

        // calculate the new Rx = ( delta^2 - Qx - Gx ) mod p
        rx = modulo((delta * delta - qx - gx), p);
        printf("Rx: ( %d^2 - %d - %d ) mod %d = %d\n", delta, qx, gx, p, rx);

        // calculate the new Yx = ( delta * (Qx – Rx) - Qy ) mod p
        ry = modulo((delta * (qx - rx) - qy), p);
        printf("Ry: ( %d * (%d – %d) - %d ) mod %d = %d\n\n", delta, qx, rx, qy, p, ry);
      }

      n--;
    }
    // print the result (Rx, Ry)
    printf("%d %d\n", rx, ry);

    // take the new first parameter (n)
    scanf("%d", &n);
  }

  return 0;
}

int is_inverse(int qx, int qy, int gx, int gy, int p){
  if (qx != gx) return 0;

  if (gy == modulo(p-qy, p)) return 1;

  return 0;
}

int modulo(int a, int b){
  int res = a % b;
  while(res < 0){
    res += b;
  }
  return res;
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

