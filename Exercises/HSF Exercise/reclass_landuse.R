reclass_landuse <- function(x) {
  fct_collapse(factor(x),
               forest = c("30","40","50","60", "70","80", "90","100"),
               grass = c("120", "140"),
               wet = c("160"),
               other_level = "other")
}

# # Example:
# hab <- get_amt_fisher_covars()
# values(hab$landuse) <- reclass_landuse(values(hab$landuse))
# 
# plot(hab$landuse)
