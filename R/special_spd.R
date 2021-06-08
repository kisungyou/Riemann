## SOME SPECIAL FUNCTIONS ON SPD 
#  (00) spd.geometry : return the geometry type for all functions below



# (00) spd.geometry -------------------------------------------------------
#' @keywords internal
#' @noRd
spd.geometry <- function(geometry){
  name.all = tolower(c("AIRM","LERM","Jeffrey","Stein"))
  name.tgt = tolower(geometry)
  
  return(match.arg(name.tgt, name.all))
}
