#' @import ggplot2
#' @method plot bekk
#' @export

plot.bekk <- function(x, diagnostic = FALSE, ...){

  if (diagnostic == FALSE) {
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




    if (inherits(x$sigma_t, "ts")) {
      plist <- vector(mode = "list", length = ncol(x$sigma_t))
      xxc <- colnames(x$sigma_t)
      for (i in 1:ncol(x$sigma_t)) {
        if (grepl('correlation', xxc[i])) {
          plist[[i]] <- suppressMessages(autoplot(x$sigma_t[,i]) + theme_bw() + ylim(-1,1) + geom_hline(yintercept = 0, col = 'red'))
        } else {
          plist[[i]] <- autoplot(x$sigma_t[,i]) + theme_bw()
        }
      }
    } else {
      xxc <- colnames(x$sigma_t)

      plist <- vector(mode = "list", length = ncol(x$sigma_t))

      for (i in 1:ncol(x$sigma_t)) {
        xx1 <- data.frame(x$sigma_t[,i])
        colnames(xx1) <- 'V1'
        if (grepl('correlation', xxc[i])) {
          plist[[i]] <- ggplot(xx1, aes(x = 1:nrow(x$sigma_t), y = V1)) + geom_line() + theme_bw()+ xlab('') + ylab('') + ylim(-1,1) + geom_hline(yintercept = 0, col = 'red')
        } else {
          plist[[i]] <- ggplot(xx1, aes(x = 1:nrow(x$sigma_t), y = V1)) + geom_line() + theme_bw()+ xlab('') + ylab('')
        }
      }
    }


    for (i in 1:ncol(x$sigma_t)) {
      plist[[i]] <- plist[[i]] + ggtitle(xxc[i])
    }

    trianglePlotGrid(plist)
  } else {
    dat <- data.frame(l = x$likelihood_iter)
    ggplot(dat, aes(x = 1:nrow(dat), y = l)) + geom_line() + ggtitle('BHHH-algorithm convergence') +
      xlab('Iteration') + ylab('log-likelihood') + theme_bw()
  }

}
