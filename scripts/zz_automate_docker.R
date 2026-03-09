#---
# Extract all packages automatically for Docker file
#---
deps_df <- renv::dependencies()

# Filter relevant files
deps <- deps_df |> 
  # string in column "Source" should contain one of "/00", "/01" or "/05"
  # using grepl
  dplyr::filter(grepl("/00|/01|/05", Source)) |> 
  dplyr::pull(Package) |> 
  unique()

dep_string <- paste0(deps, collapse = " ")
