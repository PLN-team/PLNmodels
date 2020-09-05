#include "packer.h"

// [[Rcpp::export]]
bool cpp_test_packer() {
    bool success = true;
    auto check = [&success](bool check_value, const char * context) {
        if(!check_value) {
            REprintf("Cpp internals failed: %s", context);
            success = false;
        }
    };
    const double epsilon = 1e-6;

    auto z = arma::vec{}; // 0-sized
    auto a = arma::mat(4, 10, arma::fill::randu);
    auto b = arma::vec(7, arma::fill::randu);

    const auto packer = make_packer(z, a, b, b);
    check(packer.size == 4 * 10 + 7 + 7, "packer size computation");
    check(std::get<0>(packer.elements).offset == 0, "packer offset 0");
    check(std::get<1>(packer.elements).offset == 0, "packer offset 1");
    check(std::get<2>(packer.elements).offset == 4 * 10, "packer offset 2");
    check(std::get<3>(packer.elements).offset == 4 * 10 + 7, "packer offset 3");

    auto packed = arma::vec(packer.size);
    packer.pack<0>(packed, z);
    packer.pack<1>(packed, a);
    packer.pack<2>(packed, b);
    packer.pack<3>(packed, b);

    check(packer.unpack<0>(packed).n_elem == 0, "unpack 0");
    check(arma::approx_equal(a, packer.unpack<1>(packed), "absdiff", epsilon), "unpack 1");
    check(arma::approx_equal(b, packer.unpack<2>(packed), "absdiff", epsilon), "unpack 2");
    check(arma::approx_equal(b, packer.unpack<3>(packed), "absdiff", epsilon), "unpack 3");

    packer.pack_double_or_arma<1>(packed, Rcpp::wrap(0.));
    check(packer.unpack<1>(packed).is_zero(), "pack_double_or_arma double(0.) in mat");
    packer.pack_double_or_arma<1>(packed, Rcpp::wrap(a));
    check(arma::approx_equal(a, packer.unpack<1>(packed), "absdiff", epsilon), "pack_double_or_arma mat");

    packer.pack_double_or_arma<2>(packed, Rcpp::wrap(0.));
    check(packer.unpack<2>(packed).is_zero(), "pack_double_or_arma double(0.) in vec");
    packer.pack_double_or_arma<2>(packed, Rcpp::wrap(b));
    check(arma::approx_equal(b, packer.unpack<2>(packed), "absdiff", epsilon), "pack_double_or_arma vec");

    return success;
}