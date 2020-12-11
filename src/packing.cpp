#include "packing.h"

// [[Rcpp::export]]
bool cpp_test_packing() {
    bool success = true;
    auto check = [&success](bool check_value, const char * context) {
        if(!check_value) {
            REprintf("Cpp internals failed: %s\n", context);
            success = false;
        }
    };
    const double epsilon = 1e-6;

    auto z = arma::vec{}; // 0-sized
    auto a = arma::mat(4, 10, arma::fill::randu);
    auto b = arma::vec(7, arma::fill::randu);

    // Metadata offsets
    const auto metadata = tuple_metadata(z, a, b, b);
    check(metadata.packed_size == 4 * 10 + 7 + 7, "metadata size computation");
    check(std::get<0>(metadata.elements).offset == 0, "metadata offset 0");
    check(std::get<1>(metadata.elements).offset == 0, "metadata offset 1");
    check(std::get<2>(metadata.elements).offset == 4 * 10, "metadata offset 2");
    check(std::get<3>(metadata.elements).offset == 4 * 10 + 7, "metadata offset 3");

    // Pack values
    auto packed = std::vector<double>(metadata.packed_size);
    metadata.map<0>(packed.data()) = z;
    metadata.map<1>(packed.data()) = a;
    metadata.map<2>(packed.data()) = b;
    metadata.map<3>(packed.data()) = b;

    // Check mapped values
    check(metadata.map<0>(packed.data()).n_elem == 0, "map 0");
    check(arma::approx_equal(a, metadata.map<1>(packed.data()), "absdiff", epsilon), "map 1");
    check(arma::approx_equal(b, metadata.map<2>(packed.data()), "absdiff", epsilon), "map 2");
    check(arma::approx_equal(b, metadata.map<3>(packed.data()), "absdiff", epsilon), "map 3");

    // Check non-copy or copy semantics
    check(metadata.map<1>(packed.data()).memptr() == packed.data(), "map 1 no-copy");
    check(metadata.copy<1>(packed.data()).memptr() != packed.data(), "copy 1 did copy");

    // set_from_r_value on temporary mapped packed buffer
    set_from_r_sexp(metadata.map<1>(packed.data()), Rcpp::wrap(0.));
    check(metadata.map<1>(packed.data()).is_zero(), "pack_double_or_arma double(0.) in mat");
    set_from_r_sexp(metadata.map<1>(packed.data()), Rcpp::wrap(a));
    check(arma::approx_equal(a, metadata.map<1>(packed.data()), "absdiff", epsilon), "pack_double_or_arma mat");

    set_from_r_sexp(metadata.map<2>(packed.data()), Rcpp::wrap(0.));
    check(metadata.map<2>(packed.data()).is_zero(), "pack_double_or_arma double(0.) in vec");
    set_from_r_sexp(metadata.map<2>(packed.data()), Rcpp::wrap(b));
    check(arma::approx_equal(b, metadata.map<2>(packed.data()), "absdiff", epsilon), "pack_double_or_arma vec");

    return success;
}