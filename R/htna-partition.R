# HTNA partition preservation -------------------------------------------------

# Nestimate does not depend on htna, but htna objects inherit from netobject
# and carry a small, stable node-partition contract.  Keep that contract when
# clustering splits one network into several fitted networks.
.capture_htna_partition <- function(x) {
  if (!inherits(x, "htna") || is.null(x$node_groups)) return(NULL)

  ng <- x$node_groups
  if (!is.data.frame(ng) || !all(c("node", "group") %in% names(ng))) {
    return(NULL)
  }

  group_chr <- as.character(ng$group)
  actor_levels <- x$actor_levels %||% attr(ng, "actor_levels")
  if (is.null(actor_levels) && is.factor(ng$group)) {
    actor_levels <- levels(ng$group)
  }
  if (is.null(actor_levels)) actor_levels <- unique(group_chr)

  list(
    node_groups = data.frame(
      node = as.character(ng$node),
      group = group_chr,
      stringsAsFactors = FALSE
    ),
    actor_levels = as.character(actor_levels)
  )
}

.restore_htna_partition <- function(net, partition) {
  if (is.null(partition)) return(net)

  ng <- partition$node_groups
  actor_levels <- partition$actor_levels %||% unique(as.character(ng$group))
  node_names <- as.character(net$nodes$label)
  lookup <- stats::setNames(as.character(ng$group), as.character(ng$node))
  node_actor <- unname(lookup[node_names])

  missing <- node_names[is.na(node_actor)]
  if (length(missing)) {
    stop(
      "Cannot preserve the HTNA actor partition; these clustered-network ",
      "nodes have no actor assignment: ", paste(unique(missing), collapse = ", "),
      call. = FALSE
    )
  }

  # Levels that do not cover every observed actor would be turned into NA by
  # factor() without any signal, so reject them the same way missing nodes are.
  unknown <- setdiff(node_actor, actor_levels)
  if (length(unknown)) {
    stop(
      "Cannot preserve the HTNA actor partition; these actors are absent from ",
      "'actor_levels': ", paste(unique(unknown), collapse = ", "),
      call. = FALSE
    )
  }

  # Mirror the shape htna itself builds: $node_groups$group is character,
  # $nodes$groups is a factor carrying the full actor level set.
  net$nodes$groups <- factor(node_actor, levels = actor_levels)
  net$node_groups <- data.frame(
    node = node_names,
    group = node_actor,
    stringsAsFactors = FALSE
  )
  attr(net$node_groups, "actor_levels") <- actor_levels
  net$actor_levels <- actor_levels
  class(net) <- unique(c("htna", class(net)))
  net
}

.restore_htna_group <- function(nets, partition) {
  if (is.null(partition)) return(nets)

  # `nets[] <-` keeps the clustering attributes that lapply() would drop.
  nets[] <- lapply(nets, .restore_htna_partition, partition = partition)
  class(nets) <- unique(c("htna_group", "netobject_group", "list"))
  attr(nets, "actor_levels") <- partition$actor_levels
  nets
}
