library(mediation)
library(ggplot2)
library(patchwork)
library(exact2x2)

data(jobs)
names(jobs)

boundsdata<-data.frame(trt = jobs$treat)
boundsdata$outcome <- 1.0 * (jobs$work1 == "psyemp")
boundsdata$M1 <- 1.0 * (jobs$depress2 > 2)
boundsdata$M2 <- jobs$job_dich

by(boundsdata, boundsdata$trt, summary)

con.probs <- list()
typs <- expand.grid(y = 0:1, m1 = 0:1, m2 = 0:1, x = 0:1)
for(j in 1:nrow(typs)) {
  btmp <- boundsdata[boundsdata$trt == typs[j,"x"],]
  con.probs[[with(typs[j, ], sprintf("p%s%s%s_%s", y, m1, m2, x))]] <-
    mean(btmp$outcome == typs[j,"y"] &
           btmp$M1 == typs[j,"m1"] &
           btmp$M2 == typs[j,"m2"])
}

bnds.funcs <- c("nde.000", "jnie.1", "ms2.nie.1.11", "nie.2.100")
f.bnds.list <- lapply(bnds.funcs, function(x) {

  readRDS(file.path("bndfuncs", paste0(x, ".rds")))

})
names(f.bnds.list) <- bnds.funcs

cdebnds <- readRDS(file.path("bndfuncs", "cde-bounds.rds"))


TE <- mean(subset(boundsdata, trt == 1)$outcome) -
  mean(subset(boundsdata, trt == 0)$outcome)

t.test(outcome ~ I(1 - trt), data = boundsdata)
table(boundsdata$trt, boundsdata$outcome)
binomMeld.test(86, 299, 207, 600, parmtype = "difference")

bees <- do.call(rbind, lapply(f.bnds.list, function(f) {
  do.call(f, con.probs)
}))
bees$point <- NA
bees <- rbind(bees, data.frame(lower = c(NA),
                               upper = c(NA),
                               point = c(TE)))

bees$effect <- c("'NDE-000'", "JNIE[1]", "MS^2-NIE[1]-11",
                 "NIE[2]-100", "TE")
bees$place <- c(4, 3, 2, 1, 4.5)

bees2 <- bees
bees2$place <- c(4, 0, 2, 1, 4.5)

p1 <- ggplot(bees[c(1, 2, 5),],
             aes(y = place,
                 x = point, xmin = lower, xmax = upper)) +
  geom_linerange(size = 1) + geom_point(size = 2) + theme_bw() + xlab("Bounds") +
  scale_y_continuous(breaks = bees[c(1, 2, 5),]$place,
                     labels = str2expression(bees[c(1, 2, 5),]$effect)) +
  ylab("Effect") + geom_vline(xintercept = TE, linetype = 3) +
  geom_polygon(data = data.frame(x = c(TE, bees$upper[1], TE - bees$upper[1], TE),
                                 y = c(4, 4, 3, 3)), aes(x = x, y = y),
               inherit.aes = FALSE, alpha = .2, fill = "grey20") +
  geom_polygon(data = data.frame(x = c(TE, bees$lower[1], TE - bees$lower[1], TE),
                                 y = c(4, 4, 3, 3)), aes(x = x, y = y),
               inherit.aes = FALSE, alpha = .2, fill = "grey80") +
  xlim(c(-1.5, 1))

p2<- ggplot(bees2[c(2,3,4),],
            aes(y = place,
                x = point, xmin = lower, xmax = upper)) +
  geom_linerange(size = 1)  + theme_bw() + xlab("Bounds") +
  scale_y_continuous(breaks = bees2[c(2,3,4),]$place,
                     labels = str2expression(bees2[c(2, 3, 4),]$effect)) +
  ylab("Effect") +
  geom_path(data = data.frame(x = c(bees$lower[3], bees$lower[4], sum(bees$lower[4:3])),
                              y = c(2, 1, 0)), aes(x = x, y = y),
            inherit.aes = FALSE, color = "#FD5E0F", arrow = arrow()) +
  geom_path(data = data.frame(x = c(bees$upper[3], bees$upper[4], sum(bees$upper[4:3])),
                              y = c(2, 1, 0)), aes(x = x, y = y),
            inherit.aes = FALSE, color = "#FD5E0F", arrow = arrow())+
  geom_path(data = data.frame(x = c(bees$lower[2], bees$upper[4], bees$lower[2] - bees$upper[4]),
                              y = c(0, 1, 2)), aes(x = x, y = y),
            inherit.aes = FALSE, color = "#5F3A3F", arrow = arrow()) +
  geom_path(data = data.frame(x = c(bees$upper[2], bees$lower[4], bees$upper[2] - bees$lower[4]),
                              y = c(0, 1, 2)), aes(x = x, y = y),
            inherit.aes = FALSE, color = "#5F3A3F", arrow = arrow()) +
  geom_vline(xintercept = c(bees$lower[2], bees$upper[2]), linetype = 3) +
  xlim(c(-1.5, 1.5))


p1 + p2 + plot_layout(ncol = 1)

#ggsave("jobs-fig.pdf", width = 5.5, height = 4.75)

## cdes

cdbees <- do.call(rbind, lapply(cdebnds, function(f) {
  do.call(f, con.probs)
}))
cdbees$bound <- c("'CDE-00'", "'CDE-01'",
                  "'CDE-10'", "'CDE-11'")


knitr::kable(cdbees, digits = 2)
knitr::kable(bees, digits =2 )

###
