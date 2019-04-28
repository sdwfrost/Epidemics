
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
    par[1]*exp(-distance_matrix/par[2])
  }
}

Household_Kernel = function(households){
  function(par){
    household_matrix = sapply(1:length(households), function(X) households[X] == households)
    par[1] + par[2]*household_matrix
  }
}


