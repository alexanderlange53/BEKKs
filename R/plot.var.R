#' @import ggplot2
#' @import From reshape2 melt
#' @import zoo
#' @export

plot.var <- function(x, ...) {

  if (any(class(x) == 'bekkFit')) {
    if(is.null(x$portfolio_weights)) {
      if (inherits(x$bekk$data, "ts")) {
        autoplot(x$VaR) + theme_bw() + ylab('VaR')
      } else {
        x$VaR$obs <- 1:nrow(x$VaR)
        VaR <- melt(x$VaR, id = 'obs')
        ggplot(VaR) + geom_line(aes(x = obs, y = value)) + theme_bw() + xlab('') + ylab('VaR') + facet_wrap(~variable, scales = 'free_y', ncol = 1)
      }
    } else {
      if (inherits(x$bekk$data, "ts")) {
        autoplot(x$VaR) + theme_bw() + ylab('VaR') + ggtitle('Portfolio VaR')
      } else {
        ggplot(x$VaR) + geom_line(aes(x = 1:nrow(x$VaR), y = V1)) + theme_bw() + xlab('') + ylab('VaR') + ggtitle('Portfolio VaR')
      }
    }
  } else if (any(class(x) == 'bekkForecast')) {
    if(is.null(x$portfolio_weights)) {
      sample <- x$VaR[1:(nrow(x$VaR)-x$n.ahead),]
      forc <- x$VaR[(nrow(x$VaR)-x$n.ahead+1):nrow(x$VaR),]
      cb_lower <- x$VaR_lower[(nrow(x$VaR)-x$n.ahead+1):nrow(x$VaR),]
      cb_upper <- x$VaR_upper[(nrow(x$VaR)-x$n.ahead+1):nrow(x$VaR),]

      sample$obs <- as.character(1:nrow(sample))
      forc$obs <- as.character((nrow(sample)+1):(nrow(sample)+x$n.ahead))
      cb_lower$obs <- as.character((nrow(sample)+1):(nrow(sample)+x$n.ahead))
      cb_upper$obs <- as.character((nrow(sample)+1):(nrow(sample)+x$n.ahead))

      sample <- sample[(nrow(sample)-4*x$n.ahead):nrow(sample),]
      sample$type <- as.factor('Sample')
      forc$type <- as.factor('Forecast')
      cb_lower$type <- as.factor('Forecast')
      cb_upper$type <- as.factor('Forecast')

      cb_l <- melt(cb_lower, id = c('obs', 'type'))
      cb_u <- melt(cb_upper, id = c('obs', 'type'))

      cb <- cbind(cb_l, cb_u$value)
      colnames(cb)[4:5] <- c('lower', 'upper')



      total <- rbind(sample, forc)

      VaR <- melt(total, id = c('obs', 'type'))

      cc <- merge(VaR, cb, all.x = TRUE, all.y = TRUE)

      ggplot(cc, aes(x = obs, y = value)) +
        geom_line(aes(y = lower, group = type, color = type, linetype = type), na.rm = TRUE) +
        geom_line(aes(y = upper, group = type, color = type, linetype = type), na.rm = TRUE) +
        geom_line(aes(group = type, color = type)) +
        geom_point(aes(shape = type)) +
        theme_bw() + xlab('') + ylab('VaR') +
        scale_color_manual(values = c('black', 'blue')) +
        facet_wrap(~variable, scales = 'free_y', ncol = 1) +
        theme(legend.position="bottom", legend.title = element_blank())
    } else {
      sample <- as.data.frame(x$VaR[1:(nrow(x$VaR)-x$n.ahead),])
      forc <- as.data.frame(x$VaR[(nrow(x$VaR)-x$n.ahead+1):nrow(x$VaR),])
      cb_lower <- as.data.frame(x$VaR_lower[(nrow(x$VaR)-x$n.ahead+1):nrow(x$VaR),])
      cb_upper <- as.data.frame(x$VaR_upper[(nrow(x$VaR)-x$n.ahead+1):nrow(x$VaR),])

      sample$obs <- as.character(1:nrow(sample))
      forc$obs <- as.character((nrow(sample)+1):(nrow(sample)+x$n.ahead))
      cb_lower$obs <- as.character((nrow(sample)+1):(nrow(sample)+x$n.ahead))
      cb_upper$obs <- as.character((nrow(sample)+1):(nrow(sample)+x$n.ahead))

      sample <- sample[(nrow(sample)-4*x$n.ahead):nrow(sample),]
      sample$type <- as.factor('Sample')
      forc$type <- as.factor('Forecast')

      colnames(sample)[1] <- colnames(forc)[1] <- colnames(cb_lower)[1] <- colnames(cb_upper)[1] <- 'V1'

      cb_lower$type <- as.factor('Forecast')
      cb_upper$type <- as.factor('Forecast')

      cb_l <- melt(cb_lower, id = c('obs', 'type'))
      cb_u <- melt(cb_upper, id = c('obs', 'type'))

      cb <- cbind(cb_l, cb_u$value)
      colnames(cb)[4:5] <- c('lower', 'upper')

      total <- rbind(sample, forc)

      VaR <- melt(total, id = c('obs', 'type'))

      cc <- merge(VaR, cb, all.x = TRUE, all.y = TRUE)

      ggplot(cc, aes(x = obs, y = value)) +
        geom_line(aes(y = lower, group = type, linetype = type), color = 'red', na.rm = TRUE) +
        geom_line(aes(y = upper, group = type, linetype = type), color = 'red', na.rm = TRUE) +
        geom_line(aes(group = type, color = type)) +
        geom_point(aes(shape = type)) +
        theme_bw() + xlab('') + ylab('VaR') +
        scale_color_manual(values = c('black', 'blue')) +
        theme(legend.position="bottom", legend.title = element_blank()) + ggtitle('Portfolio VaR')

          }
  }


}
