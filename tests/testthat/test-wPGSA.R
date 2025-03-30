res1 <- wPGSA(dummy_network_data, dummy_GSE143371_count)
res2 <- wPGSA(dummy_network_data, dummy_GSE143371_count, GSE143371_group)

expect_equal(length(res1), 7)
expect_equal(length(res2), 7)

expect_true(is.null(res1$group_mean))
expect_true(is.null(res1$group_se))
expect_true(is.null(res1$group_t))
expect_true(is.null(res1$pval))
expect_true(is.null(res1$fdr))

expect_true(!is.null(res2$group_mean))
expect_true(!is.null(res2$group_se))
expect_true(!is.null(res2$group_t))
expect_true(!is.null(res2$pval))
expect_true(!is.null(res2$fdr))
