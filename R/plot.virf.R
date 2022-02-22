#' @import ggplot2
#' @import grid
#' @import gridExtra
#' @export

plot.virf <- function(x, ...){

  V1 <- NULL
  l <- NULL

    trianglePlotGrid <- function(plots){
      #take a list of plots and returns a single plot where the elements in the list arranged in a triangular grid
      #plots should be a list of 1 or 3 or 6... plots to be arranged in a trianglular structure with 1 plot in the top row
      ncols <- (-1 + sqrt(1 + 8*length(plots)))/2

      grobs <- list()

      k <-1

      for (i in 1:ncols) {
        for (j in 1:ncols) {
          if (i <= j) {
            grobs[[length(grobs)+1]] <- plots[[k]]
            k <- k+1
          } else {
            grobs[[length(grobs)+1]] <- nullGrob()
          }
        }
      }
      do.call("grid.arrange", c(grobs, ncol=ncols))
    }





      xxc <- colnames(x$VIRF)

      plist <- vector(mode = "list", length = ncol(x$VIRF))

      for (i in 1:ncol(x$VIRF)) {
        xx1 <- data.frame(x$VIRF[,i])
        colnames(xx1) <- 'V1'
        if (grepl('correlation', xxc[i])) {
          plist[[i]] <- ggplot(xx1, aes(x = 1:nrow(x$VIRF), y = V1)) + geom_line() + theme_bw()+ xlab('') + ylab('') + ylim(-1,1) + geom_hline(yintercept = 0, col = 'red')
        } else {
          plist[[i]] <- ggplot(xx1, aes(x = 1:nrow(x$VIRF), y = V1)) + geom_line() + theme_bw()+ xlab('') + ylab('')
        }
      }



    for (i in 1:ncol(x$VIRF)) {
      plist[[i]] <- plist[[i]] + ggtitle(xxc[i])
    }

    trianglePlotGrid(plist)


}
