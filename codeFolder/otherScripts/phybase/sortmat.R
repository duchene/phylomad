'sortmat' <- 
function(mat, columns) 
{ 
        m <- do.call("order", as.data.frame(mat[, columns])) 
        mat[m, ] 
} 
