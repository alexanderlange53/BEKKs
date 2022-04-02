#' @import ggplot2
#' @import ggfortify
#' @import reshape2
#' @export

plot.backtest <- function(x, ...) {

  obs <- NULL
  V1 <- NULL
  lower <- NULL
  type <- NULL
  upper <- NULL


    if(is.null(x$portfolio_weights)) {
      if (inherits(x$bekk$data, "ts")) {
        autoplot(x$VaR) + autoplot(x$out_sample_returns) + theme_bw() + ylab('VaR')
      } else {
        x$VaR$obs <- 1:nrow(x$VaR)
        VaR <- melt(x$VaR, id = 'obs')
        x$out_sample_returns$obs <- 1:nrow(x$out_sample_returns)
        out_sample_returns <- melt(x$out_sample_returns, id = 'obs')
        ggplot(VaR) + geom_line(aes(x = obs, y = value)) + geom_point(data = out_sample_returns,
                                                                    mapping = aes(x = 1:nrow(x$VaR), y = V1)) + theme_bw() + xlab('') + ylab('VaR') + facet_wrap(~variable, scales = 'free_y', ncol = 1)
      }
    } else {
      if (inherits(x$data, "ts")) {
        ggplot(x$VaR) + geom_point(x$out_sample_returns) + theme_bw() + ylab('VaR') + ggtitle('Portfolio VaR')
      } else {
        ggplot(x$VaR) + geom_line(aes(x = 1:nrow(x$VaR), y = V1)) + geom_point(data = x$out_sample_returns,
                                                                              mapping = aes(x = 1:nrow(x$VaR), y = V1)) + theme_bw() + xlab('') + ylab('VaR') + ggtitle('Portfolio VaR')
      }
    }



}
