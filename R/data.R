#' Human-AI Vibe Coding Interaction Data
#'
#' Coded interaction sequences from 429 human-AI pair programming sessions
#' across 34 projects. Three coding granularities: code (32 states),
#' category (17 states), and superclass (6 states).
#'
#' @format Long-format data frames with columns:
#' \describe{
#'   \item{id}{Integer. Turn index within the session.}
#'   \item{project}{Character. Project identifier (Project_1 .. Project_34).}
#'   \item{session_id}{Character. Unique session hash.}
#'   \item{timestamp}{Character. ISO 8601 timestamp.}
#'   \item{session_date}{Character. Date of the session (YYYY-MM-DD).}
#'   \item{actor}{Character. \code{"Human"} or \code{"AI"}.}
#'   \item{code}{Character. Fine-grained action code (32 states).}
#'   \item{category}{Character. Mid-level category (17 states).}
#'   \item{superclass}{Character. High-level superclass (6 states).}
#' }
#'
#' @details
#' Nine long-format datasets are provided, filtered by actor and named by
#' granularity level:
#'
#' \tabular{lll}{
#'   \strong{Dataset} \tab \strong{Actor} \tab \strong{Granularity} \cr
#'   \code{human_ai} \tab Both \tab code (32 states) \cr
#'   \code{human_ai_cat} \tab Both \tab category (17 states) \cr
#'   \code{human_ai_super} \tab Both \tab superclass (6 states) \cr
#'   \code{human_detailed} \tab Human \tab code (32 states) \cr
#'   \code{human_cat} \tab Human \tab category (17 states) \cr
#'   \code{human_super} \tab Human \tab superclass (6 states) \cr
#'   \code{ai_detailed} \tab AI \tab code (32 states) \cr
#'   \code{ai_cat} \tab AI \tab category (17 states) \cr
#'   \code{ai_super} \tab AI \tab superclass (6 states) \cr
#' }
#'
#' Two wide-format datasets at category level (rows = sessions,
#' columns = T1, T2, ...):
#'
#' \tabular{ll}{
#'   \code{human_wide} \tab Human actions in wide sequence format \cr
#'   \code{ai_wide} \tab AI actions in wide sequence format \cr
#' }
#'
#' @source Saqr, M. (2026). Human-AI vibe coding interaction study.
#'
#' @examples
#' \donttest{
#' # Build a transition network from human category sequences
#' net <- build_network(human_wide, method = "relative")
#'
#' # Use the edge list directly
#' head(human_ai_edges)
#' }
#'
#' @name vibcoding-data
#' @aliases human_ai human_ai_cat human_ai_super
#'   human_detailed human_cat human_super
#'   ai_detailed ai_cat ai_super
#'   human_wide ai_wide
NULL

#' @rdname vibcoding-data
"human_ai"

#' @rdname vibcoding-data
"human_ai_cat"

#' @rdname vibcoding-data
"human_ai_super"

#' @rdname vibcoding-data
"human_detailed"

#' @rdname vibcoding-data
"human_cat"

#' @rdname vibcoding-data
"human_super"

#' @rdname vibcoding-data
"ai_detailed"

#' @rdname vibcoding-data
"ai_cat"

#' @rdname vibcoding-data
"ai_super"

#' @rdname vibcoding-data
"human_wide"

#' @rdname vibcoding-data
"ai_wide"

#' Human-AI Vibe Coding Edge List
#'
#' Non-aggregated edge list of all consecutive action transitions across
#' 429 sessions. Each row is one transition from one action to the next
#' within a session.
#'
#' @format A data frame with 18,918 rows and 14 columns:
#' \describe{
#'   \item{from}{Character. Source action (code level).}
#'   \item{to}{Character. Target action (code level).}
#'   \item{weight}{Integer. Always 1 (non-aggregated).}
#'   \item{session_id}{Character. Unique session hash.}
#'   \item{session}{Integer. Numeric session ID (1..429).}
#'   \item{project}{Character. Project identifier.}
#'   \item{order}{Integer. Position of the transition within the session.}
#'   \item{timepoint}{Character. ISO 8601 timestamp of the source action.}
#'   \item{from_actor}{Character. Actor of the source action.}
#'   \item{to_actor}{Character. Actor of the target action.}
#'   \item{from_category}{Character. Category of the source action.}
#'   \item{to_category}{Character. Category of the target action.}
#'   \item{from_superclass}{Character. Superclass of the source action.}
#'   \item{to_superclass}{Character. Superclass of the target action.}
#' }
#'
#' @source Saqr, M. (2026). Human-AI vibe coding interaction study.
#'
#' @examples
#' \donttest{
#' # Filter to Human -> AI transitions only
#' handoffs <- human_ai_edges[
#'   human_ai_edges$from_actor == "Human" &
#'   human_ai_edges$to_actor == "AI", ]
#' }
#'
"human_ai_edges"


#' Group Regulation in Collaborative Learning (Long Format)
#'
#' Students' regulation strategies during collaborative learning, in long
#' format. Contains 27,533 timestamped action records from multiple students
#' working in groups across two courses.
#'
#' @format A data frame with 27,533 rows and 6 columns:
#' \describe{
#'   \item{Actor}{Integer. Student identifier.}
#'   \item{Achiever}{Character. Achievement level: \code{"High"} or \code{"Low"}.}
#'   \item{Group}{Numeric. Collaboration group identifier.}
#'   \item{Course}{Character. Course identifier (\code{"A"}, \code{"B"}, or \code{"C"}).}
#'   \item{Time}{POSIXct. Timestamp of the action.}
#'   \item{Action}{Character. Regulation action (e.g., cohesion, consensus,
#'     discuss, synthesis).}
#' }
#'
#' @source Synthetically generated from the \code{group_regulation} dataset
#'   in the \pkg{tna} package.
#'
#' @seealso \code{\link{learning_activities}}, \code{\link{srl_strategies}}
#'
#' @examples
#' \donttest{
#' # Build a transition network per actor
#' net <- build_network(group_regulation_long,
#'                      method = "relative",
#'                      actor = "Actor", action = "Action", time = "Time")
#' net
#'
#' # Group networks by achievement level
#' nets <- build_network(group_regulation_long,
#'                       method = "relative",
#'                       actor = "Actor", action = "Action", time = "Time",
#'                       groups = "Achiever")
#' nets
#' }
#'
"group_regulation_long"


#' Self-Regulated Learning Strategy Frequencies
#'
#' Simulated frequency counts of 9 self-regulated learning (SRL) strategies
#' for 250 university students. Strategies are grouped into three clusters:
#' metacognitive (Planning, Monitoring, Evaluating), cognitive (Elaboration,
#' Organization, Rehearsal), and resource management (Help_Seeking, Time_Mgmt,
#' Effort_Reg). Within-cluster correlations are moderate (0.3--0.6),
#' cross-cluster correlations are weaker.
#'
#' @format A data frame with 250 rows and 9 columns. Each column is an
#'   integer count of how often the student used that strategy.
#'
#' @examples
#' \donttest{
#' net <- build_network(srl_strategies, method = "glasso",
#'                      params = list(gamma = 0.5))
#' net
#' }
#'
"srl_strategies"


#' Online Learning Activity Indicators
#'
#' Simulated binary time-series data for 200 students across 30 time points.
#' At each time point, one or more learning activities may be active (1) or
#' inactive (0). Activities: Reading, Video, Forum, Quiz, Coding, Review.
#' Includes temporal persistence (activities tend to continue across
#' adjacent time points).
#'
#' @format A data frame with 6,000 rows and 7 columns:
#' \describe{
#'   \item{student}{Integer. Student identifier (1--200).}
#'   \item{Reading}{Integer (0/1). Reading activity indicator.}
#'   \item{Video}{Integer (0/1). Video watching indicator.}
#'   \item{Forum}{Integer (0/1). Discussion forum indicator.}
#'   \item{Quiz}{Integer (0/1). Quiz/assessment indicator.}
#'   \item{Coding}{Integer (0/1). Coding practice indicator.}
#'   \item{Review}{Integer (0/1). Review/revision indicator.}
#' }
#'
#' @examples
#' \donttest{
#' net <- build_network(learning_activities, method = "cna",
#'                      actor = "student", codes = c("Reading", "Video",
#'                      "Forum", "Quiz", "Coding", "Review"), window_size = 3)
#' net
#' }
#'
"learning_activities"
