#' CAVITOPEN - an R implementation of Herve Cochard et al.'s PLC simulator
#' @description Description of the method to be provided. Based on Excel spreadsheet 'Cavitopen 3.xlsx'.
#' @param wps Water potentials (MPa) at which to evaluate the various PLCs.
#' @param rotor_length Length of the rotor used (cm)
#' @param p50_true True P50 for non-open vessels
#' @param s50_true True slope at P50 for non-open vessels
#' @param pcav_open Threshold water potential for cavitation of open vessels
#' @param vessel_length_mean Mean vessel length, used if \code{vessel_dist} is linear.
#' @param vessel_length_max Maximum vessel length, used if \code{vessel_dist} is exponential.
#' @param d_stepsize Step size along segment length to integrate PLC. Probably don't ever have to change this.
#'
#' @examples
#'
#' # The default settings reproduce the example in 'cavitopen 3.xlsx'.
#' run1 <- cavitopen()
#'
#' # Default plot. Has very few options.
#' plot(run1)
#'
#' # Smoother curves
#' run2 <- cavitopen(wps = seq(0, -12, length=101))
#' plot(run2, type='l')
#'
#' # Set parameters.
#' run3 <- cavitopen(vessel_length_max = 16, s50_true=30)
#' plot(run3)
#' @export
cavitopen <- function(wps = seq(0, -12, by=-1), # water potentials requested.
  rotor_length = 27,  # cm
  p50_true = -7.38,
  s50_true = 51.28,
  pcav_open = -0.5,
  vessel_length_mean = 10,  # used only for linear length distribution
  vessel_length_max = 10,   # used only for exponential length distribution
  vessel_dist = c("exp","linear"),
  d_stepsize = 0.0001  # should not have to change this, unless function appears slow.
){

  vessel_dist <- match.arg(vessel_dist)

  # for exponential vessel length distribution
  # I assume this should be lambda, but we are keeping typos!
  lamda <- -log(100)/vessel_length_max

  rotor_length_m <- rotor_length / 100
  rotor_r <- rotor_length_m / 2  # $C$21

  # Distance along rotor
  dist <- seq(0, rotor_r, by = d_stepsize)
  R <- rotor_r - dist

  # Open vessels
  perc_open_lin <- pmax(100 * ( 1 - (R*100)/vessel_length_mean), 0)
  perc_open_exp <- 100 * exp(lamda * R * 100)

  perc_open <- if(vessel_dist == "linear")perc_open_lin else perc_open_exp

  # Function for calculating PLCs for a given pressure.
  # Each block in the xls file basically is a call to this function.
  int_fun <- function(wp){

    # velocity, rad/s
    wrad <- sqrt(-wp*10*160000/ (2 * rotor_length^2)*100)

    # water potential along radius
    P <- -(0.5 * wrad^2 * rotor_r^2/100 - 0.5*wrad^2 * (rotor_r - R)^2/100) / 10

    # PLC due to cavitation, using 'true' P50 and S50.
    PLC_cavit <- 100 / (1 + exp(s50_true/25*(P - p50_true)))

    # Resistance (rel.) due to cavitation.
    R_cavit <- 100 / (100 - PLC_cavit)

    # PLC due to open vessels
    PLC_open <- perc_open
    ii <- min(which(P > pcav_open))
    PLC_open[ii:length(PLC_open)] <- perc_open[ii]

    # Likewise, rel. resistance.
    R_open <- 100 / (100 - PLC_open)

    # Total relative conductivity
    K_tot <- 1/R_open *(1 - PLC_cavit/100)

    # Total PLC
    PLC_tot <- 100 * (1 - K_tot)

    # Total resistance
    R_tot <- 100 / (100 - PLC_tot)

    k_to_plc <- function(x)100 * (1-x)

    PLC_cavit <- k_to_plc(1/mean(R_cavit))
    PLC_center <- 100 / (1 + exp(s50_true/25 * (wp - p50_true)))

    # Average resistance, then convert to PLC.
  return(c(PLC_cavit = PLC_cavit,
           PLC_open = k_to_plc(1/mean(R_open)),
           PLC_tot = k_to_plc(1/mean(R_tot)),
           PLC_center = PLC_center))
  }

  res <- as.data.frame(do.call(rbind, lapply(wps, int_fun)))
  res <- cbind(WP = wps, res)

  class(res) <- c("cavitopen","data.frame")

return(res)
}

#' @export
plot.cavitopen <- function(x, type='o', pch=19, legend=TRUE,  ...){

  with(x, {
    plot(WP, PLC_cavit, type=type, pch=pch,  col=palette()[1], ylim=c(0,100), ...)
    lines(WP, PLC_open, type=type, pch=pch, col=palette()[2])
    lines(WP, PLC_center, type=type, pch=pch, col=palette()[3])
    lines(WP, PLC_tot, type=type, pch=pch, col=palette()[4])
  })

  if(legend){
    legend("topright", c("PLC cavit segment", "PLC open vessels","PLC at center","PLC total"),
           lty=1, pch=pch, col=palette())
  }

}







