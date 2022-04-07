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
        autoplot(x$VaR) + autoplot(x$out_sample_returns, colour = "blue") + theme_bw() + ylab('VaR/returns')
      } else {
        x$VaR$obs <- 1:nrow(x$VaR)
        VaR <- melt(x$VaR, id = 'obs')
        x$out_sample_returns$obs <- 1:nrow(x$out_sample_returns)


        cb_l <- melt(cb_lower, id = c('obs', 'type'))
        cb_u <- melt(cb_upper, id = c('obs', 'type'))
        out_sample_returns <- melt(x$out_sample_returns, id = 'obs')
        VaR$variable = out_sample_returns$variable
        ggplot(VaR) + geom_line(aes(x = obs, y = value, colour = "Estimated VaR")) + geom_point(data = out_sample_returns, mapping = aes(x = obs, y = value, colour = "Returns"),  show.legend = TRUE) + theme_bw() + xlab('') + ylab('VaR') +  scale_color_manual(values = c('black', 'blue'), "") +
          facet_wrap(~variable, scales = 'free_y', ncol = 1)
    }
      }else {
      if (inherits(x$data, "ts")) {
        ggplot(x$VaR) + geom_point(x$out_sample_returns, colour = "blue") + theme_bw() + ylab('VaR/portfolio returns') + ggtitle('Portfolio Backtest')
      } else {
        ggplot(x$VaR) + geom_line(aes(x = 1:nrow(x$VaR), y = V1, colour="VaR")) + geom_point(data = x$out_sample_returns,
                                                                              mapping = aes(x = 1:nrow(x$VaR), y = V1, colour="Returns")) + theme_bw() + xlab('') + ylab('VaR/returns') + ggtitle('Portfolio Backtest')+ scale_color_manual(values = c('blue', 'black'), "")

      }
    }



}
