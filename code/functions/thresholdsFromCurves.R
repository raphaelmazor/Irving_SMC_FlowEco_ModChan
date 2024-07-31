find_closest_values <- function(values_str) {
  # Split the string and convert to a numeric vector
  # values_list <- as.numeric(unlist(strsplit(values_str, "\\s+")))
  values_list <- values_str
  # Sort the list based on absolute values
  sorted_list <- values_list[order(abs(values_list))]
  
  # Initialize variables to store the closest negative and positive values
  closest_negative <- NA
  closest_positive <- NA
  
  # Iterate through the sorted list to find the closest negative and positive values
  for (value in sorted_list) {
    if (value < 0 && is.na(closest_negative)) {
      closest_negative <- value
    } else if (value > 0 && is.na(closest_positive)) {
      closest_positive <- value
    }
    # Break the loop if both values are found
    if (!is.na(closest_negative) && !is.na(closest_positive)) {
      break
    }
  }
  
  # Handle cases where there are only negative or only positive values
  if (all(values_list < 0)) {
    closest_positive <- NA
  } else if (all(values_list > 0)) {
    closest_negative <- NA
  }
  
  return(list(closest_negative = closest_negative, closest_positive = closest_positive))
}

# Example usage
# values_str <- "-13991.42181  -4363.43338  -3187.51854  -2907.19687   -810.51306   -472.79613   -100.18528 -15.45081   1191.04305"
class(threshNAT)
closest_values <- find_closest_values(threshNAT)

cat("Closest negative value:", closest_values$closest_negative, "\n")
cat("Closest positive value:", closest_values$closest_positive, "\n")

closest_values

## save
save(find_closest_values, file = "code/functions/find_closest_values.RData")
