# distance matrix
dist2 <- function(x, ...) {as.dist(1-cor(t(x), method="pearson"))}

# color palette
col_ylgnbu <- colorRampPalette(rev(brewer.pal(n = 9, "YlGnBu")))


# custom colors for legend
# cols <- with(sample.info, data.frame(Source_type = levels(Source_type), color = gg_color_hue(nlevels(Source_type)))) 
# cols <- with(sample.info, data.frame(Source_type = levels(Source_type), color = palette()[1:nlevels(Source_type)])) 
# cols <- with(sample.info, data.frame(Source_type = levels(Source_type), color = color_pallete_function(nlevels(Source_type)))) 
# sample.info <- merge(sample.info,cols)

heatmap.2(x = as.matrix(res), distfun = dist2, scale = 'row', col = col_ylgnbu, 
          cex.main = 2, trace = 'none', srtCol = 45,
          main = paste(title), 
          key = T, cexRow = 1, cexCol = 1, ColSideColors = as.character(sample.info$color),margins = c(15,8))

legend("topright", legend = as.character(cols$Source_type), col = as.character(cols$color), lty= 1, lwd = 10)