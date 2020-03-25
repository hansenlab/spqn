.digestMatrix <- function(mat, digits = 6) {
    content <- sprintf(paste0("%.", digits, "f"), mat)
    # NOTE: Handling signed zero as per IEEE specs
    zero <- paste(c("0.", rep("0", digits)), collapse = "")
    content[content == paste0("-", zero)] <- zero
    digest::digest(c(content, rownames(mat), colnames(mat)))
}
