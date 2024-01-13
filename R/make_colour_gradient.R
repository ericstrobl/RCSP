make_colour_gradient = function(x, brewer_palette = "Purples") {
  min_x = min(x)
  max_x = max(x)
  range_x = max_x - min_x
  x_scaled = (x - min_x) / range_x
  
  # Chopping out first colour as it's too light to work well as a
  #   point colour
  colours = scales::brewer_pal("seq", brewer_palette)(5)[2:5]
  
  colour_vals = scales::colour_ramp(colours)(x_scaled)
  colour_vals
}