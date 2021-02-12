#' @import ggplot2
#' @method plot bekk
#' @export

plot.bekk <- function(x, ...){
  xx <- cbind(seq(1, nrow(x$sigma_t)), x$sigma_t)
  colnames(xx)[1] <- 'time'
  xx <- melt(xx, id = 'time')
  ggplot(xx) + geom_line(aes(x = time, y = value, group = variable)) +
    facet_wrap(~variable, scales = 'free_y', ncol = ncol(result0$data)) + theme_bw() + ylab('')
}
