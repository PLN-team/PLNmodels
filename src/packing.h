// Packing / unpacking helper.
// Automatically generate a struct containing offsets, size, types of a list of numeric values.
// Provide functions to access a linearized buffer as its components.
// See tests in packing.cpp for usage.
#pragma once

#include <RcppArmadillo.h>
#include <cstddef> // size_t
#include <tuple>   // packer system
#include <utility> // move, forward

#define ARMA_EXTRA_DEBUG

// Stores type, dimensions and offset for a single T object.
// Must be specialised ; see specialisations for double/arma::vec/arma::mat below.
//
// Required API when specialised:
// Constructor(const T & reference_object, std::size_t & current_offset);
// <?> map(double * buffer) const; -> Returns an object wrapping part of the buffer.
// <?> copy(const double * buffer) const; -> Returns a new object copying data from the buffer.

template <typename T> struct PackingMetadata;

template <> struct PackingMetadata<double> {
    std::size_t offset;

    PackingMetadata(const double &, std::size_t & current_offset) {
        offset = current_offset;
        current_offset += 1;
    }

    double & map(double * buffer) const { return buffer[offset]; }
    const double & map(const double * buffer) const { return buffer[offset]; }
    double copy(const double * buffer) const { return buffer[offset]; }
};

template <> struct PackingMetadata<arma::vec> {
    std::size_t offset;
    arma::uword size;

    PackingMetadata(const arma::vec & v, std::size_t & current_offset) {
        offset = current_offset;
        size = v.n_elem;
        current_offset += size;
    }

    arma::vec map(double * buffer) const { return {&buffer[offset], size, false /*copy*/, true /*strict*/}; }
    arma::vec map(const double * buffer) const {
        // arma::vec has no way to represent a read-only vector.
        // For safety, arma::vec mapped from a const buffer must not be modified ; they should be declared const.
        return map(const_cast<double *>(buffer));
    }
    arma::vec copy(const double * buffer) const { return {&buffer[offset], size}; }
};

template <> struct PackingMetadata<arma::mat> {
    std::size_t offset;
    arma::uword rows;
    arma::uword cols;

    PackingMetadata(const arma::mat & m, std::size_t & current_offset) {
        offset = current_offset;
        rows = m.n_rows;
        cols = m.n_cols;
        current_offset += rows * cols;
    }

    arma::mat map(double * buffer) const { return {&buffer[offset], rows, cols, false /*copy*/, true /*strict*/}; }
    arma::mat map(const double * buffer) const {
        // arma::mat has no way to represent a read-only vector.
        // For safety, arma::mat mapped from a const buffer must not be modified ; they should be declared const.
        return map(const_cast<double *>(buffer));
    }
    arma::mat copy(const double * buffer) const { return {&buffer[offset], rows, cols}; }
};

// Stores packing information for multiple objects of types T0,T1,...,TN.
// Created (with type deduction) using tuple_metadata(T0, ..., TN) below.
//
// In what follows, 'i' is an index (0 <= i <= N) representing T_i.
// Using an enum to name these indexes is recommended :
// auto metadata = tuple_metadata(a, b, c);
// enum { A_ID, B_ID, C_ID }; // Index names, respectively {0, 1, 2}
// metadata.map<A_ID>(storage); // Value representing the A value in packed storage.
// metadata.map<A_ID>(storage) = value; // Sets A value in storage to 'value'.
template <typename... Types> struct TuplePackingMetadata {
    std::tuple<PackingMetadata<Types>...> elements; // Packing info for each element (offset, type, dimensions)
    std::size_t packed_size;                        // Total number of doubles when packed

    // metadata.map<i>(buffer) : extract value of T_i in 'buffer', referencing the buffer
    template <std::size_t I> auto map(double * buffer) const -> decltype(std::get<I>(elements).map(buffer)) {
        return std::get<I>(elements).map(buffer);
    }
    template <std::size_t I> auto map(const double * buffer) const -> decltype(std::get<I>(elements).map(buffer)) {
        return std::get<I>(elements).map(buffer);
    }

    // metadata.copy<i>(buffer) : extract value of T_i in 'buffer', creating a new value
    template <std::size_t I> auto copy(double * buffer) const -> decltype(std::get<I>(elements).copy(buffer)) {
        return std::get<I>(elements).copy(buffer);
    }
};

template <typename... Types> TuplePackingMetadata<Types...> tuple_metadata(const Types &... values) {
    // Initialize TuplePackingMetadata<Types...> using brace init, which guarantees evaluation order (required!).
    // Will call each PackingMetadata<T> constructor in order, increasing offset every time.
    // Then the final offset value will be copied into the size field.
    std::size_t current_offset = 0;
    return {
        std::tuple<PackingMetadata<Types>...>{PackingMetadata<Types>(values, current_offset)...},
        current_offset,
    };
}

// Set value from an R SEXP.
// Uniformly fill a vec/mat if given a float value.
// If given a vec/mat, set individual values, expecting the same dimensions.
inline void set_from_r_sexp(double & v, SEXP sexp) {
    v = Rcpp::as<double>(sexp);
}

inline void set_from_r_sexp(arma::vec & v, SEXP sexp) {
    if(Rcpp::is<double>(sexp)) {
        v.fill(Rcpp::as<double>(sexp));
    } else {
        v = Rcpp::as<arma::vec>(sexp);
    }
}
inline void set_from_r_sexp(arma::vec && v, SEXP sexp) {
    set_from_r_sexp(v, sexp); // Allow set of temporary
}

inline void set_from_r_sexp(arma::mat & m, SEXP sexp) {
    if(Rcpp::is<double>(sexp)) {
        m.fill(Rcpp::as<double>(sexp));
    } else {
        m = Rcpp::as<arma::mat>(sexp);
    }
}
inline void set_from_r_sexp(arma::mat && m, SEXP sexp) {
    set_from_r_sexp(m, sexp); // Allow set of temporary
}