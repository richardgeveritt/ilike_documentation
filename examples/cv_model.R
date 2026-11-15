get_data <- function()
{
  utils::data("constant_velocity_y", package = "ilike", envir = environment())
  constant_velocity_y[1:20]
}
