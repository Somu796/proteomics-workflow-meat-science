barplot_efs <- function(efs_table, filter_by = 0.5, order = TRUE){
  
  if(order == TRUE){
    efs_table <- efs_table[, colSums(efs_table)>filter_by]
    efs_table <- efs_table[, order(colSums(efs_table))]
  }
  
  paranr = length(efs_table[1,])  
  if(paranr>100){
    b =  colSums(efs_table)
    #b= a[order(a)]
    
    barplot(b,
            ylim=c(0,1),
            main= 'Ensemble Feature Selection',
            xlab = "Features",
            ylab = "Importance values",
            axisnames = FALSE
    )
  }
  if(paranr<35) h=10
  else h=(paranr/5)
  
  names=colnames(efs_table)
  cols = c('goldenrod1','navy',
           'royalblue','indianred3',
           'darkolivegreen1','darkgreen',
           'darkolivegreen3','chartreuse4')
  
  par(mar=c(5, 4, 4, 10), xpd=TRUE)
  barplot= barplot(efs_table,
                   xlim=c(0,1),
                   main= 'Ensemble Feature Selection',
                   horiz=T,
                   las=2,
                   names.arg=abbreviate(names),
                   col=cols)
  legend("topright", inset=c(-0.2,0), legend=row.names(efs_table),col=cols, lty=1, lwd=12 )
  text(colSums(efs_table)+0.065,barplot,
       format(round(colSums(efs_table), 2),T))
  # segments(1, 0, 1, 1.25*paranr, lty = 3, col = "gray40")
}