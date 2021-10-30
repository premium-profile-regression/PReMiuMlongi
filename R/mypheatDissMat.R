mypheatDissMat <- function (dissimObj, true_clusters=NULL)
{
  nSbj <- dissimObj$disSimRunInfoObj$nSubjects
  col.labels <- c("0", "0.5", "1")
  colours <- colorRampPalette(c("white", "black"))(10)
  dissMat <- vec2mat(dissimObj$disSimMat, nrow = nSbj)
  SymMat <- 1 - dissMat

  if(!missing(true_clusters)){
    rownames(SymMat) <- paste0("r", 1:nSbj)
    colnames(SymMat) <- paste0("c", 1:nSbj)
    myannot <- data.frame("true_clusters"=true_clusters)
    #myannot <- as.data.frame(t(dissimObj$disSimRunInfoObj$xMat))
    myannot_row <- apply(t(myannot), 2, as.factor)
    myannot_row <- as.data.frame(myannot_row)
    rownames(myannot_row) <- paste0("r", 1:nSbj)

    myannot_col <- apply(t(myannot), 2, as.factor)
    myannot_col <- as.data.frame(myannot_col)
    rownames(myannot_col) <- paste0("c", 1:nSbj)
    colnames(myannot_row) <- "true_clusters"
    col_annot <- c("#6B0077",   "#7665A4",  "#8DA3CA",  "#B7D5E4",  "#F1F1F1")
    #col_annot <- c("#D3A362",   "#8BF0AD",  "#9AFFFF",  "#FFE0FF",  "#FFEEB3")

out <- pheatmap::pheatmap(SymMat, annotation_row = myannot_row,
             color = colorRampPalette(c( "white", "black"))(50),
             annotation_names_row=F,
             #clustering_callback = callback(hcl,SymMat),
             #annotation_col = myannot_col,
             show_rownames=F, show_colnames=F, cluster_rows=F,cluster_cols=F,
             annotation_colors = list(
               true_clusters= c("1" = col_annot[1], "2" = col_annot[2], "3" = col_annot[3],
                                "4" = col_annot[4], "5" = col_annot[5])
               #c("1" = "darkblue", "2" = "purple", "3" = "salmon",
               #                "4" = "rosybrown", "5" = "Aquamarine")
               #c("1" = "darkblue", "2" = "purple", "3" = "salmon",
               #  "4" = "lightblue", "5" = "lightgreen")
             ))
    }else{

      pheatmap::pheatmap(SymMat,
               color = colorRampPalette(c( "white", "black"))(50),
               #annotation_col = myannot_col,
               show_rownames=F, show_colnames=F )
    }
}
