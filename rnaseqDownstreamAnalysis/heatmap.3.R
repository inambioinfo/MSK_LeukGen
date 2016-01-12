heatmap.3(x = as.matrix(res), dist.FUN = dist2, scale = 'row', Rowv = T,Colv = T,
          color.FUN = col_ylgnbu, margin.for.labCol = .5,margin.for.labRow = 1,
          cex.main = 2, trace = 'none',
          adjCol = c(0.1,0), offsetCol=1,
          key = T, cexRow = 2, cexCol = 1)

