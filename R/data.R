#' Human-AI Vibe Coding Interaction Data (Long Format)
#'
#' Coded turns from 429 human-AI pair-programming ("vibe coding") sessions
#' across 34 projects, held between 5 February and 8 March 2026, in long
#' format with one row per code. `human_long` holds the codes of the human
#' turns and `ai_long` those of the AI turns (428 of the sessions have AI
#' turns). The two data frames share their columns and their sessions, and
#' interleave on `order_in_session`.
#'
#' @format Two data frames with the same 9 columns:
#' \describe{
#'   \item{message_id}{Integer identifier of the message (turn); the codes of
#'     one message share it, and no identifier occurs in both data frames.}
#'   \item{project}{Character project identifier, `Project_1` to
#'     `Project_34`.}
#'   \item{session_id}{Character session identifier, a 12-character hash.}
#'   \item{timestamp}{Integer Unix time (seconds) of midnight, UTC+2, on
#'     `session_date`. It has day resolution, so the rows of one day share
#'     it; `build_network(time = "timestamp")` keeps tied rows in row order,
#'     which follows `order_in_session`. With the default `time_threshold`
#'     (900 seconds), each change of day also starts a new sequence, so the
#'     84 sessions that run over more than one day are split (526 sequences
#'     instead of 429); to keep every session whole, use
#'     `time = "order_in_session"`.}
#'   \item{session_date}{Character date of the row (`YYYY-MM-DD`); 84
#'     sessions run over more than one day.}
#'   \item{code}{Character interaction code. The two data frames use
#'     disjoint code sets: 9 labels in `human_long` (Command, Correct,
#'     Frustrate, Inquire, Interrupt, Refine, Request, Specify, Verify) and
#'     8 in `ai_long` (Ask, Delegate, Execute, Explain, Investigate, Plan,
#'     Repair, Report), 17 in total.}
#'   \item{cluster}{Character high-level cluster of the code. `human_long`
#'     uses Directive, Evaluative and Metacognitive; `ai_long` uses Action,
#'     Communication and Repair.}
#'   \item{code_order}{Integer order of the code within its message,
#'     starting at 1.}
#'   \item{order_in_session}{Integer position of the row in its session,
#'     counting the rows of both data frames together.}
#' }
#'
#' \tabular{ll}{
#'   \code{human_long} \tab Human turns only, 10,796 rows \cr
#'   \code{ai_long}    \tab AI turns only, 8,551 rows \cr
#' }
#'
#' @source Saqr, M. (2026). Human-AI vibe coding interaction study.
#'   \url{https://saqr.me/blog/2026/human-ai-interaction-cograph/}
#'
#' @examples
#' net <- build_network(human_long, method = "tna",
#'                      action = "code", actor = "session_id",
#'                      time = "timestamp")
#' net
#'
#' @name long-data
#' @aliases ai_long
NULL

#' @rdname long-data
"human_long"

#' @rdname long-data
"ai_long"


#' Group Regulation in Collaborative Learning (Long Format)
#'
#' Students' regulation actions during collaborative learning, in long
#' format: 27,533 timestamped actions by 2,000 students working in 200
#' groups across three courses, one row per action.
#'
#' @format A data frame with 27,533 rows and 6 columns:
#' \describe{
#'   \item{Actor}{Integer student identifier, 1 to 2,000.}
#'   \item{Achiever}{Character achievement level, `"High"` or `"Low"`.}
#'   \item{Group}{Numeric collaboration group identifier, 1 to 200; each
#'     student belongs to one group.}
#'   \item{Course}{Character course identifier, `"A"`, `"B"` or `"C"`.}
#'   \item{Time}{POSIXct timestamp of the action.}
#'   \item{Action}{Character regulation action, one of adapt, cohesion,
#'     consensus, coregulate, discuss, emotion, monitor, plan and
#'     synthesis.}
#' }
#'
#' @source Synthetically generated from the \code{group_regulation} dataset
#'   in the \pkg{tna} package.
#'
#' @seealso \code{\link{learning_activities}}, \code{\link{srl_strategies}}
#'
#' @examples
#' net <- build_network(group_regulation_long,
#'                      method = "relative",
#'                      actor = "Actor", action = "Action", time = "Time")
#' net
#'
"group_regulation_long"


#' Self-Regulated Learning Strategy Frequencies
#'
#' Simulated counts of how often 250 university students used 9
#' self-regulated learning (SRL) strategies, in three clusters:
#' metacognitive (Planning, Monitoring, Evaluating), cognitive (Elaboration,
#' Organization, Rehearsal) and resource management (Help_Seeking,
#' Time_Mgmt, Effort_Reg). Strategies of the same cluster were simulated to
#' correlate moderately and those of different clusters more weakly; in the
#' bundled sample the correlations range from 0.18 to 0.49 within clusters
#' and from -0.13 to 0.31 across them.
#'
#' @format A data frame with 250 rows and 9 columns, one row per student
#'   and one column per strategy. Every column is numeric (double) and
#'   holds a whole-number count of how often that student used the
#'   strategy; the observed values range from 0 to 37, with no missing
#'   values.
#'
#' @examples
#' net <- build_network(srl_strategies, method = "glasso",
#'                      params = list(gamma = 0.5))
#' net
#'
"srl_strategies"


#' Online Learning Activity Indicators
#'
#' Simulated binary (one-hot) indicator time series of six online learning
#' activities (Reading, Video, Forum, Quiz, Coding, Review) for 200
#' students over 30 time points. Each row is one time point of one student
#' and marks every activity active then (1) or not (0); at least one
#' activity is active in every row, and several often are. Activities tend
#' to continue from one time point to the next. The rows of each student
#' are consecutive and in time order, so the data are ready for [wtna()]
#' with `actor = "student"`.
#'
#' @format A data frame with 6,000 rows and 7 columns:
#' \describe{
#'   \item{student}{Integer student identifier, 1 to 200, with 30 rows
#'     each.}
#'   \item{Reading}{Integer reading indicator (0/1).}
#'   \item{Video}{Integer video-watching indicator (0/1).}
#'   \item{Forum}{Integer discussion-forum indicator (0/1).}
#'   \item{Quiz}{Integer quiz or assessment indicator (0/1).}
#'   \item{Coding}{Integer coding-practice indicator (0/1).}
#'   \item{Review}{Integer review or revision indicator (0/1).}
#' }
#'
#' @examples
#' head(learning_activities)
#'
"learning_activities"


#' ChatGPT Self-Regulated Learning Scale Scores
#'
#' Scores on five self-regulated learning (SRL) scales for 1,000 responses
#' that ChatGPT generated to a validated SRL questionnaire. The data are
#' part of a larger set comparing the responses of seven large language
#' models with human norms.
#'
#' @format A data frame with 1,000 rows and 5 columns, one row per
#'   response; each column is the numeric mean of one scale's items (13
#'   items for `CSU`, 4 for `TA` and 9 for each of the others), with
#'   observed values from 2 to 6.5 and no missing values:
#' \describe{
#'   \item{CSU}{Cognitive Strategy Use.}
#'   \item{IV}{Intrinsic Value.}
#'   \item{SE}{Self-Efficacy.}
#'   \item{SR}{Self-Regulation.}
#'   \item{TA}{Test Anxiety.}
#' }
#'
#' @source Vogelsmeier, L.V.D.E., Oliveira, E., Misiejuk, K.,
#'   Lopez-Pernas, S., & Saqr, M. (2025). Delving into the psychology of
#'   Machines: Exploring the structure of self-regulated learning via
#'   LLM-generated survey responses. \emph{Computers in Human Behavior},
#'   173, 108769. \doi{10.1016/j.chb.2025.108769}
#'
#' @examples
#' net <- build_network(chatgpt_srl, method = "glasso",
#'                      params = list(gamma = 0.5))
#' net
#'
"chatgpt_srl"


#' Student Engagement Trajectories
#'
#' Wide-format sequences of student engagement over 15 consecutive
#' observations, one row per student and one column per observation. The
#' states are `"Active"`, `"Average"` and `"Disengaged"`. A missing value
#' occurs only as a run at the end of a row, which marks drop-out; two rows
#' are missing throughout.
#'
#' @format A character matrix with 138 rows and 15 columns named `"1"` to
#'   `"15"`, in time order. Each entry is `"Active"`, `"Average"`,
#'   `"Disengaged"` or `NA`. The trailing `NA` runs are what the
#'   right-censored sequence verbs ([actor_endpoints()],
#'   [mark_terminal_state()]) are built to read.
#'
#' @examples
#' sequence_plot(trajectories, main = "Engagement trajectories")
#' sequence_plot(trajectories, k = 3)
#' sequence_plot(trajectories, type = "distribution")
#'
"trajectories"
