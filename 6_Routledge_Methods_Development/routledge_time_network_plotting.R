# Load required libraries 

##devtools::install_github("reconhub/epicontacts", ref = "timeline")
pacman::p_load(
  rio,          # File import
  epicontacts,  # Epicontacts
  here,         # File locator
  tidyverse,    # Data management + ggplot2 graphics
  remotes       # Package installation from github
)

# As <- A/rowSums(A) # creating normalised matrix
# As[is.nan(As)] <- 0 # unclear if this is right
# saveRDS(As, file = "tempA_mat.rds")
testAmat <- read_rds("tempA_mat.rds")
maximum <- unlist(apply(testAmat, 1, which.max))

test_linelist <- data.frame(id = paste0("a", seq(1:34)),
                            date_onset = as.Date(linelist$Date),
                            outcome = "bloop")
test_contacts <- data.frame(infector = paste0("a", maximum),
                            case_id = paste0("a", not_imported),
                            prob = 1)
test_contacts <- rbind(test_contacts, c("a19", "a24", 0.5))
test_contacts$prob <- as.numeric(test_contacts$prob)

epic <- make_epicontacts(
  linelist = test_linelist,
  contacts = test_contacts,
  id = "case_id",
  from = "infector",
  to = "case_id",
  directed = TRUE
)

plot(
  epic,
  x_axis = "date_onset",
  arrow_size = 0.5,
  node_size = 13,
  node_color = "R_i",
  edge_width = 'prob',
  label = FALSE, 
  height = 700,
  width = 1200
)
