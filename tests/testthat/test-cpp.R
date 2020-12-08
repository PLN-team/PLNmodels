context("test-cpp")

test_that("PLN: cpp internals are sane", {
    expect_true(cpp_test_nlopt())
    expect_true(cpp_test_packing())
})