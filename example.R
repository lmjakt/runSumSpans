y <- readRDS('example.rds')

plot(y, type='l')

y.d <- diff(y)
plot(y.d)

## prepare values for find_spans

pos <- which(y.d < 0)
val <- rep(1, length(pos))

dyn.load( 'src/run_sum_spans.so')

spans <- as.data.frame( .Call('rs_spans', as.double(val), as.double(pos), 0.1) )

plot( y, type='l' )
rect( spans$start, 0, spans$end, spans$score )

## and that seems to work reasonably well..

tmp2 <- readRDS("example_2.rds")
tmp2.d <- which( diff(tmp2[[1]]) < 0 )


spans <- .Call('find_spans', as.double(rep(1, length(tmp2.d)), as.double(tmp2.d), 0.1))
