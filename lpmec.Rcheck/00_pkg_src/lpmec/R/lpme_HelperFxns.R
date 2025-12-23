f2n <- function(x){as.numeric( as.character( x ) ) }
f2i <- function(f_){lpmec_env$jnp$array(f_,lpmec_env$jnp$int32)}
f2a <- function(x){lpmec_env$jnp$array(x,lpmec_env$jnp$float32)}
ai <- as.integer
lpmec_env <- new.env( parent = emptyenv() )
