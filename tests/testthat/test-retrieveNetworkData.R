# First Data Retrieval
out <- retrieveNetworkData()

# Test First Data Retrieval
expect_true(is.data.frame(out))

# Second Data Retrieval
out2 <- retrieveNetworkData()

# Test Secound Data Retrieval
expect_true(identical(out, out2))

# Test Cache File
bfc <- BiocFileCache()
info <- bfcinfo(bfc)
expect_true(nrow(info) > 0)
