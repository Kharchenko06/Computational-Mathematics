#include "gmp_reference.h"

#include <cassert>
#include <gmp.h>

namespace euclidean {

double gmp_euclidean(std::span<const double> a, std::span<const double> b) {
    assert(a.size() == b.size());

    // выбрано 256 бит чтобы ошибка была намного меньше double
    constexpr mp_bitcnt_t PREC = 256;

    mpf_t ga, gb, diff, sq, sum, result;
    mpf_init2(ga,     PREC);
    mpf_init2(gb,     PREC);
    mpf_init2(diff,   PREC);
    mpf_init2(sq,     PREC);
    mpf_init2(sum,    PREC);
    mpf_init2(result, PREC);

    mpf_set_d(sum, 0.0);

    for (std::size_t i = 0; i < a.size(); ++i) {
        mpf_set_d(ga, a[i]);
        mpf_set_d(gb, b[i]);
        mpf_sub(diff, gb, ga);
        mpf_mul(sq, diff, diff);
        mpf_add(sum, sum, sq);
    }

    mpf_sqrt(result, sum);

    // округление к ближайшему double
    const double answer = mpf_get_d(result);

    mpf_clear(ga);
    mpf_clear(gb);
    mpf_clear(diff);
    mpf_clear(sq);
    mpf_clear(sum);
    mpf_clear(result);

    return answer;
}

}