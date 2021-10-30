myheatDissMat <- function (dissimObj, main = NULL, xlab = NULL, ylab = NULL)
{
  nSbj <- dissimObj$disSimRunInfoObj$nSubjects
  col.labels <- c("0", "0.5", "1")
  colours <- colorRampPalette(c("white", "black"))(10)
  dissMat <- vec2mat(dissimObj$disSimMat, nrow = nSbj)

  heatmap(1 - dissMat, keep.dendro = FALSE, symm = TRUE, Rowv = NULL,
          labRow = FALSE, labCol = FALSE, margins = c(4.5, 4.5),
          col = colours, main = main, xlab = xlab, ylab = ylab)
  plotrix::color.legend(0.95, 0.7, 1, 1, legend = col.labels, colours,
                        gradient = "y", align = "rb")
}
