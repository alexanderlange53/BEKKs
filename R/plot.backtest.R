#' @import ggplot2
#' @import ggfortify
#' @import reshape2
#' @export
#'
plot.backtest <- function(x, ...) {

  obs <- NULL
  type <- NULL
  time <- NULL
  V1 <- NULL
  Portfolio <- NULL


    if(is.null(x$portfolio_weights)) {
      if (inherits(x$VaR, c("ts"))) {
        colnames(x$VaR) = colnames(x$out_sample_returns)

        t_series <- as.data.frame(x$VaR)
        t_series$time <- time(x$VaR)
        t_series <- melt(t_series, id ="time")

        out_sample_returns <- as.data.frame(x$out_sample_returns)
        out_sample_returns$time <- time(x$out_sample_returns)
        out_sample_returns <- melt(out_sample_returns, id ="time")

        ggplot(t_series) + geom_point(data = out_sample_returns, mapping = aes(x = time, y = value, colour = "Returns"),  show.legend = TRUE) + theme_bw() + xlab('') + ylab('Returns/VaR') +
          geom_line(aes(x = time, y = value, colour = "Estimated VaR")) + scale_color_manual(values = c('black', 'blue'), "") +
          facet_wrap(~variable, scales = 'free_y', ncol = 1)+
          theme(legend.position="bottom", legend.title = element_blank())



      }else if (inherits(x$VaR, c("xts","zoo"))) {
        names(x$VaR) = names(x$out_sample_returns)

        t_series <- as.data.frame(x$VaR)
        t_series$time <- time(x$VaR)
        t_series <- melt(t_series, id ="time")

        out_sample_returns <- as.data.frame(x$out_sample_returns)
        out_sample_returns$time <- time(x$out_sample_returns)
        out_sample_returns <- melt(out_sample_returns, id ="time")


        ggplot(t_series) + geom_point(data = out_sample_returns, mapping = aes(x = time, y = value, colour = "Returns"),  show.legend = TRUE) + theme_bw() + xlab('') + ylab('Returns/VaR') +
          geom_line(aes(x = time, y = value, colour = "Estimated VaR")) + scale_color_manual(values = c('black', 'blue'), "") +
          facet_wrap(~variable, scales = 'free_y', ncol = 1)+
          theme(legend.position="bottom", legend.title = element_blank())

        } else {

        names(x$VaR) = names(x$out_sample_returns)

        x$VaR$obs <- 1:nrow(x$VaR)
        VaR <- melt(x$VaR, id = 'obs')
        x$out_sample_returns$obs <- 1:nrow(x$out_sample_returns)



        out_sample_returns <- melt(x$out_sample_returns, id = 'obs')


        ggplot(VaR) + geom_point(data = out_sample_returns, mapping = aes(x = obs, y = value, colour = "Returns"),  show.legend = TRUE) + theme_bw() + xlab('') + ylab('VaR') +
          geom_line(aes(x = obs, y = value, colour = "Estimated VaR")) + scale_color_manual(values = c('black', 'blue'), "") +
          facet_wrap(~variable, scales = 'free_y', ncol = 1)+
          theme(legend.position="bottom", legend.title = element_blank())
    }
      }else {
      if (inherits(x$VaR, c("xts","zoo"))) {

        ggplot(data=x$VaR) + geom_point(data=x$out_sample_returns, mapping = aes(x = time(x$VaR), y = Portfolio, colour="Returns")) +
          geom_line(x$VaR , mapping = aes(x = time(x$VaR), y = Portfolio, colour="Estimated VaR")) + theme_bw() + xlab("") + ylab('Portfolio returns/VaR') + ggtitle('Portfolio Backtest')+ scale_color_manual(values = c('blue', 'black'), "") +
          theme(legend.position="bottom", legend.title = element_blank())
      } else {
        ggplot(x$VaR) + geom_point(data = x$out_sample_returns,
                                                                              mapping = aes(x = 1:nrow(x$VaR), y = Portfolio, colour="Returns")) + theme_bw() + xlab('') + ylab('Returns/VaR') +
        geom_line(aes(x = 1:nrow(x$VaR), y = Portfolio, colour="Estimated VaR")) + ggtitle('Portfolio Backtest')+ scale_color_manual(values = c('blue', 'black'), "") +
          theme(legend.position="bottom", legend.title = element_blank())

      }
    }



}
