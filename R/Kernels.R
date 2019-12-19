
#'
#' Kernels
#'


# ==== Homogeneous Kernels ====

Homogeneous_Kernel = function(N){
  function(inf_par){
    inf_par
  }
}

#' Contact_Kernel Given a matrix containing information on whether individuals in an
#'                epidemic network are able to make infectious contact with eachother
#'                and returns a function which can be given infectious process parameters
#'                to return infectious pressure between individuals
#' @param contact_matrix A Matrix whose ij-th entry specifies whether an infected indiviual i
#'                       can make infectious contact with susceptible individual j
#' @param SYM            A logical arguement. When set to TRUE, it will assume that contact
#'                       is mutual, giving a symmetrical matrix. If FALSE, it mutual infectious
#'                       contact is not assumed.
Contact_Kernel = function(contact_matrix, SYM = TRUE){
  if(SYM == TRUE){
    contact_matrix[lower.tri(contact_matrix)] = Matrix::t(contact_matrix)[lower.tri(contact_matrix)]
  }
  #contact_matrix = as(contact_matrix, "sparseMatrix")
  function(beta){
    beta*contact_matrix
  }
}

Contact_Kernel2 = function(contact_matrix){
  contact_matrix = as(contact_matrix, "sparseMatrix")
  function(beta){
    beta[1] + beta[2]*contact_matrix
  }
}



Spatial_Kernel = function(coord_matrix){
  distance_matrix = spDists(x = random_coords, y = random_coords)
  distance_matrix = Matrix::as(distance_matrix, "sparseMatrix")
  function(par){
    par[1]*exp(-distance_matrix/par[2])
  }
}


