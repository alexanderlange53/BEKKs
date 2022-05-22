#' @import ggplot2
#' @import grid
#' @import gridExtra
#' @export

plot.virf <- function(x, ...){

  V1 <- NULL
  l <- NULL
  lower <- NULL
  upper <- NULL

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
        xxci <- data.frame(x$VIRF_lower[,i],x$VIRF_upper[,i])
        colnames(xx1) <- 'V1'
        colnames(xxci) <- c("lower","upper")

        if (grepl('correlation', xxc[i])) {
          plist[[i]] <- ggplot(xx1, aes(x = 1:nrow(x$VIRF), y = V1)) + geom_line() + geom_ribbon(data=xxci,aes(ymin=lower,ymax=upper),alpha=0.3) + theme_bw() + xlab('') + ylab('') + geom_hline(yintercept = 0, col = 'red')
        } else {
          plist[[i]] <- ggplot(xx1, aes(x = 1:nrow(x$VIRF), y = V1)) + geom_line() + geom_ribbon(data=xxci,aes(ymin=lower,ymax=upper),alpha=0.3) + theme_bw()+ xlab('') + ylab('') + geom_hline(yintercept = 0, col = 'red')
        }
      }



    for (i in 1:ncol(x$VIRF)) {
      plist[[i]] <- plist[[i]] + ggtitle(xxc[i])
    }

    trianglePlotGrid(plist)


}
