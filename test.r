'''
mat_list = list(snv = matrix(c(1, 0, 1, 1, 1, 0, 0, 1, 1), nrow = 3),
                indel = matrix(c(1, 0, 0, 0, 1, 0, 1, 0, 0), nrow = 3))
rownames(mat_list$snv) = rownames(mat_list$indel) = c("g1", "g2", "g3")
colnames(mat_list$snv) = colnames(mat_list$indel) = c("s1", "s2", "s3")
mat_list
mat_list$indel = mat_list$indel[1:2, 1:2]
mat_list
mat_list = unify_mat_list(mat_list)
mat_list
oncoPrint(mat_list,
    alter_fun_list = list(
        background = function(x, y, w, h) NULL,
        snv = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "red", col = NA)),
        indel = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = "blue", col = NA))
    ), col = c(snv = "red", indel = "blue"))
oncoPrint(mat_list)
'''
