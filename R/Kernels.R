
#'
#' Kernels
#'


Contact_Kernel = function(contact_matrix){
  function(beta){
    beta*contact_matrix
  }
}

Spatial_Kernel = function(coord_matrix){
  function(par){
    distance_matrix = spDists(x = random_coords, y = random_coords)
    par[1]*exp(-distance/par[2])
  }
}


