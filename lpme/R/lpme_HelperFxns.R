f2n <- function(x){as.numeric( as.character( x ) ) }
f2i <- function(f_){lpme_env$jnp$array(f_,lpme_env$jnp$int32)}
f2a <- function(x){lpme_env$jnp$array(x,lpme_env$jnp$float32)}
ai <- as.integer
lpme_env <- new.env( parent = emptyenv() )
