library(causaloptim)  ## version 0.8.0, available from https://github.com/sachsmc/causaloptim
library(data.table)
library(ggplot2)

rmvunif <- function(n, K, alpha) {
  y <- matrix(rgamma(n * K, alpha), ncol = K)
  y / rowSums(y)
}


b2 <- graph_from_literal(X -+ Y, X -+ M1, X -+ M2, Ul -+ X, M1 -+ Y, M2 -+ Y,
                         Ur -+ M1, Ur -+ M2, Ur -+ Y, M1 -+ M2)
V(b2)$leftside <- c(1, 0, 0, 0, 1, 0)
V(b2)$latent <- c(0, 0, 0, 0, 1, 1)
E(b2)$rlconnect <- rep(0, 10)
E(b2)$edge.monotone <- rep(0, 10)

plot(b2)

nde <- "p{Y(X = 1, M1(X = %s), M2(X = %s, M1(X = %s))) = 1} - p{Y(X = 0, M1(X = %s), M2(X = %s, M1(X = %s))) = 1}"

obj <- analyze_graph(b2, constraints = NULL,
                     effectt = sprintf(nde, 0, 0, 0, 0, 0, 0))

###


bnds.funcs <- c("nde.000", "jnie.1", "ms2.nie.1.11", "nie.2.100")
## probabilities are P(Y, M1, M2 | X)
f.bnds.list <- lapply(bnds.funcs, function(x) {

  readRDS(file.path("bndfuncs", paste0(x, ".rds")))

})
names(f.bnds.list) <- bnds.funcs


width <- excl0 <- upper <- lower <- NULL
for(j in 1:length(obj$variables)){
  sim.qs <- rep(0, length(obj$variables))
  sim.qs[j] <- 1
  sim.qs <- as.list(sim.qs)
  names(sim.qs) <- obj$variables

  sim.ps <- as.list(c(obj$R[-1, ] %*% unlist(sim.qs)))
  names(sim.ps) <- obj$parameters

  compd.bnds <- lapply(f.bnds.list, function(f) {

    do.call(f, sim.ps)

  })


  res <- do.call(rbind, compd.bnds)
  res$effect <- names(compd.bnds)
  #res$iter <- i
  res$width <- res$upper - res$lower
  res$exclude0 <- res$upper <= 0 | res$lower >= 0
  width <- rbind(width, res$width)
  excl0 <- rbind(excl0, res$exclude0)
  upper <- rbind(upper, res$upper)
  lower <- rbind(lower, res$lower)

}

## Table 1

for(i in 1:4) {

  print(names(f.bnds.list)[i])
  print(table(lower[, i], upper[, i]))


}


## Remaining figures

rep_one <- function(i, alp) {

  sim.qs <- as.list(rmvunif(1, length(obj$variables),
                            alp)[1, ])
  names(sim.qs) <- obj$variables

  sim.ps <- as.list(c(obj$R[-1, ] %*% unlist(sim.qs)))

  names(sim.ps) <- obj$parameters

  compd.bnds <- lapply(f.bnds.list, function(f) {

    do.call(f, sim.ps)

  })


  res <- do.call(rbind, compd.bnds)
  res$effect <- names(compd.bnds)
  #res$iter <- i
  res$width <- res$upper - res$lower
  res$exclude0 <- res$upper <= 0 | res$lower >= 0
  res

}

set.seed(1320)
outwidth2 <- NULL
outwidth <- outexcl <- NULL
alphas <- seq(0.00001, .001, length.out = 12)
for(j in alphas) {
  sim.results <- lapply(1:1000, rep_one, j)

  outwidth <- rbind(outwidth, rowMeans(sapply(sim.results,
                                              function(x) x$width <.999)))
  outexcl <- rbind(outexcl, rowMeans(sapply(sim.results,
                                            function(x) x$exclude0)))
  outwidth2 <- rbind(outwidth2, data.frame(width = unlist(lapply(sim.results, function(x) x$width)),
                                           effect = bnds.funcs,
                                           alpha = j))

}



bnds.funcs2 <-  c("'NDE-000'", "JNIE[1]", "MS^2-NIE[1]-11",
                  "NIE[2]-100")
names(bnds.funcs2) <- bnds.funcs
outwidth2$label <- bnds.funcs2[outwidth2$effect]

res1 <- data.frame(perc = c(outwidth), alpha = rep(alphas, 4),
                   bound = rep(bnds.funcs2, each = length(alphas)))

res1$lower <- sapply(res1$perc, function(p) {
  binom.test(p * 1000, 1000)$conf.int[1]
})
res1$upper <- sapply(res1$perc, function(p) {
  binom.test(p * 1000, 1000)$conf.int[2]
})

ggplot(res1, aes(x = alpha*1000, y = perc, ymin = lower, ymax = upper)) +
  geom_ribbon(alpha = .2) +
  geom_line() + geom_point() +
  facet_wrap(~ bound, labeller = "label_parsed") + theme_bw() +
  ylab("Proportion with width less than 1")

#ggsave("sim-fig-width.pdf", width = 6.75, height = 5.5)

oog <- data.table(outwidth2)

gtext <- oog[, .(pone = mean(width >= 1)), by = .(label, alpha)]
gtext$width <- 2.5
gtext$text <- sprintf("%.1f", (100 * (1 - gtext$pone)))

labs <- c(".01", round(sort(unique(1000 * outwidth2$alpha)), 1)[-1])
labs[7]<- ".55"

#cairo_pdf("sim-fig-width2.pdf", width = 6.75, height = 4.5)
ggplot(outwidth2, aes(x = factor(alpha * 1000), y = width)) +
  geom_jitter(alpha = .5, width = .3, color = "grey85") +
  geom_boxplot(outlier.alpha = 0, fill = NA) +
  geom_label(data = gtext, aes(x = factor(alpha * 1000),
                               y = width, label = text),
             size = 2.5, color = "#FD5E0F",
             label.padding = unit(0.075, "lines"),
             label.r = unit(0.0, "lines")) +
  facet_wrap(~ label, labeller = "label_parsed") + theme_bw() + xlab("alpha * 1000") +
  scale_x_discrete(labels = labs) +
  scale_y_continuous(breaks = c(0, .5, 1, 1.5, 2, 2.5),
                     labels = c("0.0", "0.5", "1.0", "1.5", "2.0", "% < 1")) +
  theme(axis.text.y = element_text(color = c(rep("black", 5), "#FD5E0F")),
        panel.grid.major.x = element_blank())

#dev.off()

#ggplot(outwidth2, aes(x = alpha * 1000, y = width)) + geom_bin2d(binwidth = c(1/12, .1))

res2 <- data.frame(perc = c(outexcl), alpha = rep(alphas, 4),
                   bound = rep(bnds.funcs2, each = length(alphas)))

res2$lower <- sapply(res2$perc, function(p) {
  binom.test(p * 1000, 1000)$conf.int[1]
})
res2$upper <- sapply(res2$perc, function(p) {
  binom.test(p * 1000, 1000)$conf.int[2]
})

ggplot(res2, aes(x = alpha*1000, y = perc, ymin = lower, ymax = upper)) +
  geom_ribbon(alpha = .2) +
  geom_line() + geom_point() +
  facet_wrap(~ bound, labeller = "label_parsed") + theme_bw() +
  ylab("Proportion excluding 0")
#ggsave("sim-fig-excl.pdf", width = 6.75, height = 5.5)


ggplot(subset(res2, bound == "'NDE-000'"),
       aes(x = alpha*1000, y = perc, ymin = lower, ymax = upper)) +
  geom_ribbon(alpha = .2) +
  geom_line() + geom_point() +
  facet_wrap(~ bound, labeller = "label_parsed") + theme_bw() +
  ylab("Proportion excluding 0")
#ggsave("sim-fig-excl2.pdf", width = 4.25, height = 1.75)

