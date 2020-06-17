context("test-cpp")

test_that("PLN: cpp internals are sane", {
    expect_true(cpp_internal_tests())
})