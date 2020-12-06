// Packing / unpacking helper.
// Automatically generate a struct containing offsets, size, types of a list of arma values.
// Provide functions to store or extract the arma values into a linearized arma::vec.
// See tests in packer.cpp for usage.

#include <RcppArmadillo.h>
#include <cstddef> // size_t
#include <tuple>   // packer system
#include <utility> // move, forward

// Stores type, dimensions and offset for a single T object
// Must be specialised ; see specialisations for arma::vec/arma::mat below
//
// Required API when specialised:
// Constructor(const T & reference_object, arma::uword & current_offset);
// T unpack(const arma::vec & packed_storage);
// void pack(arma::vec & packed_storage, "T-like arma expression type" expr);
// void pack_double_or_arma(arma::vec & packed_storage, SEXP r_value); (see OptimizerConfiguration)
template <typename T> struct PackedInfo;

// All following implementation use vec.subvec() to access slices of the packed vector.
// The "prefered" way to give indexces is using span(start, end).
// However end is inclusive and unsigned, which causes underflow for 0-sized span of offset 0.
// Thus I use the less intuitive (but correct) form: subvec(offset, arma::size(size, 1))

template <> struct PackedInfo<arma::vec> {
    arma::uword offset;
    arma::uword size;

    PackedInfo(const arma::vec & v, arma::uword & current_offset) {
        offset = current_offset;
        size = v.n_elem;
        current_offset += size;
    }

    arma::vec unpack(const arma::vec & packed) const { return packed.subvec(offset, arma::size(size, 1)); }

    template <typename Expr> void pack(arma::vec & packed, Expr && expr) const {
        packed.subvec(offset, arma::size(size, 1)) = std::forward<Expr>(expr);
    }

    void pack_double_or_arma(arma::vec & packed, SEXP r_value) const {
        if(Rcpp::is<double>(r_value)) {
            packed.subvec(offset, arma::size(size, 1)).fill(Rcpp::as<double>(r_value));
        } else {
            pack(packed, Rcpp::as<arma::vec>(r_value));
        }
    }
};

template <> struct PackedInfo<arma::mat> {
    arma::uword offset;
    arma::uword rows;
    arma::uword cols;

    PackedInfo(const arma::mat & m, arma::uword & current_offset) {
        offset = current_offset;
        rows = m.n_rows;
        cols = m.n_cols;
        current_offset += rows * cols;
    }

    arma::mat unpack(const arma::vec & packed) const {
        return arma::reshape(packed.subvec(offset, arma::size(rows * cols, 1)), arma::size(rows, cols));
    }

    template <typename Expr> void pack(arma::vec & packed, Expr && expr) const {
        // Handles: mat expressions, vec expressions
        packed.subvec(offset, arma::size(rows * cols, 1)) = arma::vectorise(std::forward<Expr>(expr));
    }

    void pack_double_or_arma(arma::vec & packed, SEXP r_value) const {
        if(Rcpp::is<double>(r_value)) {
            packed.subvec(offset, arma::size(rows * cols, 1)).fill(Rcpp::as<double>(r_value));
        } else {
            pack(packed, Rcpp::as<arma::mat>(r_value));
        }
    }
};

// Packer : stores packing information for multiple objects of types T0,T1,...,TN.
// Created (with type deduction) using make_packer(T0, ..., TN) below.
//
// In what follows, 'i' is an index (0 <= i <= N) representing T_i.
// Using an enum to name these indexes is recommended :
// auto packer = make_packer(a, b, c);
// enum { A_ID, B_ID, C_ID }; // Index names, respectively {0, 1, 2}
// packer.pack<A_ID>(storage, value);
template <typename... Types> struct Packer {
    std::tuple<PackedInfo<Types>...> elements; // Packing info for each element (offset, type, dimensions)
    arma::uword size;                          // Total number of packed elements

    // packer.unpack<i>(storage) : extract value of T_i in 'storage'
    template <std::size_t Index> auto unpack(const arma::vec & packed) const
        -> decltype(std::get<Index>(elements).unpack(packed)) {
        return std::get<Index>(elements).unpack(packed);
    }

    // packer.pack<i>(storage, value) : stores 'value' at T_i's location in 'storage'
    template <std::size_t Index, typename Expr> void pack(arma::vec & packed, Expr && expr) const {
        std::get<Index>(elements).pack(packed, std::forward<Expr>(expr));
    }

    // packer.pack_double_or_arma<i>(storage, r_value) :
    // if 'r_value' is a double, fills T_i location in 'storage' with the double.
    // if 'r_value' is an arma::mat/vec, check dimensions and fills T_i's location with it.
    template <std::size_t Index> void pack_double_or_arma(arma::vec & packed, SEXP r_value) const {
        std::get<Index>(elements).pack_double_or_arma(packed, r_value);
    }
};

template <typename... Types> Packer<Types...> make_packer(const Types &... values) {
    // Initialize Packer<Types...> using brace init, which guarantees evaluation order (required here !).
    // Will call each PackedInfo<T> constructor in order, increasing offset every time.
    // Then the final offset value will be copied into the size field.
    arma::uword current_offset = 0;
    return {
        std::tuple<PackedInfo<Types>...>{PackedInfo<Types>(values, current_offset)...},
        current_offset,
    };
}