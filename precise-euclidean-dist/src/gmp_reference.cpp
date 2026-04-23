#include "gmp_reference.h"

#include <gmp.h>

namespace euclidean {

double gmp_euclidean(const Point3D& p1, const Point3D& p2) {
    // выбрано 256 бит чтобы ошибка была намного меньше double
    constexpr mp_bitcnt_t PREC = 256;

    mpf_t gx1, gy1, gz1, gx2, gy2, gz2;
    mpf_t dx,  dy,  dz,  dx2, dy2, dz2, sum, result;

    mpf_init2(gx1, PREC);  mpf_init2(gy1, PREC);  mpf_init2(gz1, PREC);
    mpf_init2(gx2, PREC);  mpf_init2(gy2, PREC);  mpf_init2(gz2, PREC);
    mpf_init2(dx,  PREC);  mpf_init2(dy,  PREC);  mpf_init2(dz,  PREC);
    mpf_init2(dx2, PREC);  mpf_init2(dy2, PREC);  mpf_init2(dz2, PREC);
    mpf_init2(sum, PREC);  mpf_init2(result, PREC);

    // конвертация double без потери точности
    mpf_set_d(gx1, p1.x);  mpf_set_d(gy1, p1.y);  mpf_set_d(gz1, p1.z);
    mpf_set_d(gx2, p2.x);  mpf_set_d(gy2, p2.y);  mpf_set_d(gz2, p2.z);

    mpf_sub(dx, gx2, gx1);
    mpf_sub(dy, gy2, gy1);
    mpf_sub(dz, gz2, gz1);

    mpf_mul(dx2, dx, dx);
    mpf_mul(dy2, dy, dy);
    mpf_mul(dz2, dz, dz);

    mpf_add(sum, dx2, dy2);
    mpf_add(sum, sum, dz2);
    mpf_sqrt(result, sum);

    // округление к ближайшему double
    const double answer = mpf_get_d(result);

    mpf_clear(gx1);  mpf_clear(gy1);  mpf_clear(gz1);
    mpf_clear(gx2);  mpf_clear(gy2);  mpf_clear(gz2);
    mpf_clear(dx);   mpf_clear(dy);   mpf_clear(dz);
    mpf_clear(dx2);  mpf_clear(dy2);  mpf_clear(dz2);
    mpf_clear(sum);  mpf_clear(result);

    return answer;
}

}