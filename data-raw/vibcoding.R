# Build the 12 vibcoding datasets from raw CSV
# Run from package root: source("data-raw/vibcoding.R")

raw <- read.csv(
  "/Users/mohammedsaqr/Library/CloudStorage/GoogleDrive-saqr@saqr.me/.tmp/1097498/Git/Saqrvibcodingtna.csv",
  stringsAsFactors = FALSE
)
raw <- raw[order(raw$session_id, raw$id, raw$timestamp), ]

# Numeric session mapping
session_ids <- unique(raw$session_id)
session_map <- setNames(seq_along(session_ids), session_ids)

# ---- Helper: long -> wide (category level) ----
.to_wide <- function(data, action_col = "category", group_col = "session_id") {
  sessions <- split(data[[action_col]], data[[group_col]])
  max_len <- max(vapply(sessions, length, integer(1L)))
  col_names <- paste0("T", seq_len(max_len))
  mat <- vapply(sessions, function(s) {
    c(s, rep(NA_character_, max_len - length(s)))
  }, character(max_len))
  df <- as.data.frame(t(mat), stringsAsFactors = FALSE)
  colnames(df) <- col_names
  df
}

# ---- Split by actor ----
human <- raw[raw$actor == "Human", ]
ai <- raw[raw$actor == "AI", ]

# ---- 9 long-format datasets ----
human_ai          <- raw
human_ai_cat      <- raw
human_ai_super    <- raw
human_detailed    <- human
human_cat         <- human
human_super       <- human
ai_detailed       <- ai
ai_cat            <- ai
ai_super          <- ai

# ---- 2 wide-format datasets (category level) ----
human_wide <- .to_wide(human, "category")
ai_wide    <- .to_wide(ai, "category")

# ---- Edge list ----
sessions_split <- split(raw, raw$session_id)
human_ai_edges <- do.call(rbind, lapply(sessions_split, function(s) {
  if (nrow(s) < 2L) return(NULL)
  n <- nrow(s)
  data.frame(
    from            = s$code[-n],
    to              = s$code[-1L],
    weight          = 1L,
    session_id      = s$session_id[-n],
    session         = session_map[s$session_id[1L]],
    project         = s$project[-n],
    order           = seq_len(n - 1L),
    timepoint       = s$timestamp[-n],
    from_actor      = s$actor[-n],
    to_actor        = s$actor[-1L],
    from_category   = s$category[-n],
    to_category     = s$category[-1L],
    from_superclass = s$superclass[-n],
    to_superclass   = s$superclass[-1L],
    stringsAsFactors = FALSE
  )
}))
rownames(human_ai_edges) <- NULL

# ---- Save all 12 ----
usethis::use_data(
  human_ai, human_ai_cat, human_ai_super,
  human_detailed, human_cat, human_super,
  ai_detailed, ai_cat, ai_super,
  human_wide, ai_wide,
  human_ai_edges,
  overwrite = TRUE
)
