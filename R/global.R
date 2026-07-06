# Declare data-frame column names used in ggplot aes() so R CMD check
# does not flag them as undefined global variables.
utils::globalVariables(c(
  "total_reads", "pct_mouse_10x", "call_10x",
  "our_classification", "mouse_10x", "human_10x"
))
