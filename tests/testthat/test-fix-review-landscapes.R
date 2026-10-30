# Regression oracles for R08/R09 from the 2026-10-30 review.
test_that("landscapes replay weighting, concatenation and bootstrap units", {
  d <- data.frame(T1=c('A','A','B'),T2=c('A','B','A'),
                  T3=c('A',NA,'B'),T4=c('B',NA,'A'))
  n <- build_network(d,method='relative',params=list(weighted=TRUE))
  cups <- summary(stability_landscape(n,iter=0))
  expect_equal(cups$stay,c(1/3,0))
  expect_equal(cups$stay,unname(diag(n$weights)))
  trend <- summary(landscape_trend(n,window=4,iter=0))
  expect_equal(trend$depth,cups$depth)

  d <- d[c('T1','T2')]
  n <- build_network(d,method='relative',params=list(concat=3))
  land <- stability_landscape(n,iter=10,seed=1)
  cups <- summary(land)
  expect_equal(cups$stay,c(2/3,1/2))
  expect_equal(land$fits$network$n_sequences,1L)
  expect_equal(cups$ci_lower,cups$depth)
  expect_equal(cups$ci_upper,cups$depth)
  expect_equal(summary(landscape_trend(n,window=6,iter=0))$depth,cups$depth)

  combined <- build_network(data.frame(T1=c('A','A','B','A','A'),
    T2=c('A','B','A','A','B')),method='relative',
    params=list(concat=3,weighted=TRUE))
  expect_equal(summary(stability_landscape(combined,iter=0))$stay,
               unname(diag(combined$weights)))
})

test_that("merged covariates must be constant within concatenated units", {
  d <- data.frame(T1=rep(c('A','B'),6),T2=rep(c('B','A'),6),
                  grade=rep(1:6,each=2))
  n <- build_network(d,method='relative',params=list(concat=2))
  tr <- landscape_trend(n,along='grade',iter=0,k=5,n_grid=2)
  expect_equal(sort(unique(summary(tr)$value)),c(1.25,5.75))
  d$grade <- seq_len(nrow(d))
  n <- build_network(d,method='relative',params=list(concat=2))
  expect_error(landscape_trend(n,along='grade',iter=0),
               'constant within each concatenated sequence',
               class='nestimate_bad_input')
})

test_that("landscapes refuse transformed weights rather than changing chains", {
  d <- data.frame(T1=c('A','A','B'),T2=c('A','B','A'))
  expect_error(stability_landscape(build_network(d,method='relative',
    scaling='max'),iter=0),class='nestimate_bad_input')
  expect_error(stability_landscape(build_network(d,method='relative',
    threshold=.6),iter=0),class='nestimate_bad_input')
  n <- build_network(d,method='relative')
  n$weights['A','A'] <- 0
  expect_error(stability_landscape(n,iter=0),
               'cannot be reconstructed',class='nestimate_bad_input')
})

test_that("tied trimmed quantiles retain a nondegenerate surface axis", {
  skip_on_cran()
  d <- as.data.frame(matrix(rep(c('A','B','A','B'),each=21),nrow=21))
  d[21,] <- c('A','A','B','A')
  surf <- suppressWarnings(landscape_surface(d,k=5))
  expect_equal(length(unique(surf$slices$value)),25L)
  expect_lt(min(surf$slices$value),max(surf$slices$value))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(),add=TRUE)
  expect_identical(plot(surf),surf)
  d$grade <- c(rep(10,20),20)
  trend <- landscape_trend(d,along='grade',iter=0)
  expect_equal(length(unique(summary(trend)$value)),20L)
  expect_identical(plot(trend,type='surface'),trend)
})

test_that("weighted surface estimates agree with a hand-pooled slice", {
  d <- data.frame(T1=c('A','A','B','B','A','B'),
    T2=c('A','B','A','B','A','A'),T3=c('B',NA,'B','A','B','A'),
    T4=c('A',NA,NA,'B','A','B'))
  net <- build_network(d,method='relative',params=list(weighted=TRUE))
  surf <- suppressWarnings(landscape_surface(net,measure='longitudinal_entropy',
                                             k=5,n_grid=3))
  m <- sequence_metrics(d)$longitudinal_entropy
  y <- sort(unique(surf$slices$value))[2]
  dist <- abs(m-y)
  h <- max(sort(dist)[5],diff(range(m))/50)*1.01
  w <- ifelse(dist<h,(1-(dist/h)^3)^3,0)
  # Each row's A->A and total outgoing A counts, weighted by observed length.
  aa <- c(1,0,0,0,1,1)
  from_a <- c(2,1,1,1,2,2)
  lengths <- rowSums(!is.na(d))
  stay_a <- sum(w*aa/lengths)/sum(w*from_a/lengths)
  got <- surf$slices$depth[surf$slices$value==y & surf$slices$state=='A']
  expect_equal(got,-log1p(-stay_a))
})


test_that("surface measures include prepared boundary states", {
  d <- data.frame(T1=c('A','B','A','B','A','B'),
    T2=c('B','A','A','B','B','A'),T3=c('A','B','B','A','A','B'))
  net <- build_network(d,method='relative',params=list(begin_state='A'))
  surf <- suppressWarnings(landscape_surface(net,measure='self_loop_tendency',
                                             n_grid=3,k=5))
  expected <- c(1/3,0,2/3,1/3,1/3,0)
  endpoints <- quantile(expected,c(.05,.95),names=FALSE)
  grid <- seq(endpoints[1],endpoints[2],length.out=3)
  expect_equal(sort(unique(surf$slices$value)),grid)
  distance <- abs(expected-grid[2])
  bandwidth <- max(sort(distance)[5],diff(range(expected))/50)*1.01
  w <- ifelse(distance<bandwidth,(1-(distance/bandwidth)^3)^3,0)
  stay <- sum(w*c(1,0,2,0,1,0))/sum(w*c(2,2,3,1,2,2))
  row <- surf$slices[surf$slices$value==grid[2] & surf$slices$state=='A',]
  expect_equal(row$depth,-log1p(-stay))
})
