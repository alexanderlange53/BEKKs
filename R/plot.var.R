#' @import ggplot2
#' @import From reshape2 melt
#' @export

plot.var <- function(x, ...) {

  if (any(class(x) == 'bekkFit')) {
    if(is.null(x$portfolio_weights)) {
      if (inherits(x$bekk$bekkfit$data, "ts")) {
        autoplot(x$VaR) + theme_bw() + ylab('VaR')
      } else {
        x$VaR$obs <- 1:nrow(x$VaR)
        VaR <- melt(x$VaR, id = 'obs')
        ggplot(VaR) + geom_line(aes(x = obs, y = value)) + theme_bw() + xlab('') + ylab('VaR') + facet_wrap(~variable, scales = 'free_y', ncol = 1)
      }
    } else {
      if (inherits(x$bekk$sigma_t, "ts")) {
        autoplot(x$VaR) + theme_bw() + ylab('VaR') + ggtitle('Portfolio VaR')
      } else {
        ggplot(x$VaR) + geom_line(aes(x = 1:nrow(x$VaR), y = V1)) + theme_bw() + xlab('') + ylab('VaR') + ggtitle('Portfolio VaR')
      }
    }
  } else if (any(class(x) == 'bekkForecast')) {
    if(is.null(x$portfolio_weights)) {
      sample <- x$VaR[1:(nrow(x$VaR)-x$n.ahead),]
      forc <- x$VaR[(nrow(x$VaR)-x$n.ahead+1):nrow(x$VaR),]

      if (inherits(x$bekk$bekkfit$data, "ts")) {
        autoplot(x$VaR) + theme_bw() + ylab('VaR')
      } else {

        sample$obs <- as.character(1:nrow(sample))
        forc$obs <- as.character((nrow(sample)+1):(nrow(sample)+x$n.ahead))
        sample <- sample[(nrow(sample)-4*x$n.ahead):nrow(sample),]
        sample$type <- as.factor('Sample')
        forc$type <- as.factor('Forecast')

        total <- rbind(sample, forc)

        VaR <- melt(total, id = c('obs', 'type'))
        ggplot(VaR, aes(x = obs, y = value, group = type, color = type, linetype = type, shape = type)) + geom_line() + geom_point() + theme_bw() + xlab('') + ylab('VaR') +
          scale_color_manual(values = c('black', 'blue')) + facet_wrap(~variable, scales = 'free_y', ncol = 1) +
          theme(legend.position="bottom", legend.title = element_blank())
      }
    } else {
      sample <- as.data.frame(x$VaR[1:(nrow(x$VaR)-x$n.ahead),])
      forc <- as.data.frame(x$VaR[(nrow(x$VaR)-x$n.ahead+1):nrow(x$VaR),])

      if (inherits(x$bekk$bekkfit$data, "ts")) {


        autoplot(x$VaR) + theme_bw() + ylab('VaR') + ggtitle('Portfolio VaR')
      } else {
        sample$obs <- as.character(1:nrow(sample))
        forc$obs <- as.character((nrow(sample)+1):(nrow(sample)+x$n.ahead))
        sample <- sample[(nrow(sample)-4*x$n.ahead):nrow(sample),]
        sample$type <- as.factor('Sample')
        forc$type <- as.factor('Forecast')

        colnames(sample)[1] <- colnames(forc)[1] <- 'V1'

        total <- rbind(sample, forc)

        VaR <- melt(total, id = c('obs', 'type'))

        ggplot(VaR, aes(x = obs, y = value, group = type, color = type, linetype = type, shape = type)) + geom_line() + geom_point() + theme_bw() + xlab('') + ylab('VaR') +
          scale_color_manual(values = c('black', 'blue')) +
          theme(legend.position="bottom", legend.title = element_blank()) + ggtitle('Portfolio VaR')
      }
    }
  }


}
