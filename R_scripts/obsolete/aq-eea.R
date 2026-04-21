HPC <- F
setwd(
  "~/Library/Mobile Documents/com~apple~CloudDocs/Lavoro/PhD Bergamo/R/GitHub/GRINS-Spoke0-WP2"
)
period_HPC <- seq.Date(as.Date("2018-12-25"), as.Date("2019-01-10"), by = "days")

d <- period_HPC[9]
p <- "AQ-FRK/v.3.0.6/"
source(paste0(p, "Z_settings/Z_settings.R"), local = TRUE)

p <- "AQ-FRK/v.3.0.1/"
source(paste0(p, "A_input/A_input.R"), local = TRUE)
# as.Date(A1_aq_eea@data$time)
# no B
p <- "AQ-FRK/v.3.0.1/"
source(paste0(p, "C_BAUs/C_BAUs.R"), local = TRUE)

# adding eea t-1
d <- period_HPC[8]
p <- "AQ-FRK/v.3.0.6/"
source(paste0(p, "Z_settings/Z_settings.R"), local = TRUE)
p <- "AQ-FRK/v.3.0.1/"
crs_wgs84 <- CRS(SRS_string = "EPSG:4326")
if (!HPC) {
  p <- "AQ-FRK/v.3.0.0/"
  load(paste0(p, "A_input/data/input/AQ_EEA_v100_ST.rda"))
}
slot(AQ_EEA_v100_ST@sp, "proj4string") <- crs_wgs84
A1_aq_eea_t1 <- #temporal subset
  AQ_EEA_v100_ST[, which(index(AQ_EEA_v100_ST@time) %in% sub_t[2])] #v.3.0.1
#spatial subset
A1_aq_eea_t1 <- A1_aq_eea_t1[A1_aq_eea_t1@coords[, 1] > bbox[1] + .25 &
                               A1_aq_eea_t1@coords[, 1] < bbox[2] - .25 &
                               A1_aq_eea_t1@coords[, 2] > bbox[3] + .25 &
                               A1_aq_eea_t1@coords[, 2] < bbox[4] - .25, ]

if (!HPC) {
  p <- "AQ-FRK/v.3.0.0/"
  load(paste0(p, "A_input/data/input/AQ_CAMS_v100_ST.rda"))
}
A2_aq_cams_t1 <- #temporal subset
  AQ_CAMS_v100_ST[, which(index(AQ_CAMS_v100_ST@time) %in% sub_t[2])] #v.3.0.1
#spatial subset
A2_aq_cams_t1 <- A2_aq_cams_t1[A2_aq_cams_t1@coords[, 1] > bbox[1] &
                                 #maybe change?
                                 A2_aq_cams_t1@coords[, 1] < bbox[2] &
                                 A2_aq_cams_t1@coords[, 2] > bbox[3] &
                                 A2_aq_cams_t1@coords[, 2] < bbox[4], ]
names(A2_aq_cams_t1@data)[4] <-
  "AQ_CAMS_NO2"
#
p <- "AQ-FRK/v.3.0.1/"

setwd(
  "~/Library/Mobile Documents/com~apple~CloudDocs/Lavoro/PhD Bergamo/R/GitHub/Italy_AQ_AmortisedLatticeKrig"
)

save.image(paste0("data/starting_Lattice_", d + 1, ".RData"))
library(fields)
library(LatticeKrig)

Z <- as.matrix(over(A1_aq_eea, C2_BAUs_ext)[, c(5:7, 10, 11, 12, 14)])
Z_a2 <- as.matrix(over(A1_aq_eea, A2_aq_cams[, 4]))
Z <- cbind(Z, Z_a2)
head(Z)
# Z[,-8] <- scale(Z[,-8])
Z <- cbind(Z, A1_aq_eea_t1@data$AQ_EEA_NO2)
colnames(Z)[9] <- "eea_t1"
head(Z)

# training
A1_aq_eea_df <-
  A1_aq_eea@data


s <- A1_aq_eea_df[, c(1, 2)]
y <-  A1_aq_eea_df$AQ_EEA_NO2

set.seed(as.numeric(sub_t[1]))
test <- sample(1:length(y), length(y) / 10)
train <- 1:length(y)
train <- train[-which(train %in% test)]

y_train <- y[train]
test <- test[order(test)]
y_test <- y[test]

s_train <- s[train, ]
s_test <- s[test, ]

Z_train <- Z[train, ]
Z_test <- Z[test, ]

naidx_test <- is.na(Z_test[, "eea_t1"])
y_test <- y_test[!naidx_test]
s_test <- s_test[!naidx_test, ]
Z_test <- Z_test[!naidx_test, ]

y_train_mod <- y_train
naidx <- is.na(y_train_mod)
y_train_mod <- y_train_mod[!naidx]
s_train_mod <- s_train[!naidx, ]
Z_train_mod <- Z_train[!naidx, ]
naidx <- is.na(Z_train_mod[, "eea_t1"])
y_train_mod <- y_train_mod[!naidx]
s_train_mod <- s_train_mod[!naidx, ]
Z_train_mod <- Z_train_mod[!naidx, ]

Lmodel <- LatticeKrig(s_train_mod, y_train_mod, Z = Z_train_mod[, -8])

obj <- Lmodel

ncov <- 8

yhat_test <- predict(object = obj, x = s_test, Z = Z_test[, c(1:7, 9)])


# sqrt(mean((y_test - yhat_test)^2,na.rm=T))

ysehat_test <- predictSE(object = obj, x = s_test, Z = Z_test[, c(1:7, 9)])

y_CI_lower <- yhat_test - 1.96 * ysehat_test
y_CI_upper <- yhat_test + 1.96 * ysehat_test

setwd(
  "~/Library/Mobile Documents/com~apple~CloudDocs/Lavoro/PhD Bergamo/R/GitHub/GRINS-Spoke0-WP2"
)
d <- period_HPC[9]
p <- "AQ-FRK/v.3.0.6/"
source(paste0(p, "Z_settings/Z_settings.R"), local = TRUE)
setwd(
  "~/Library/Mobile Documents/com~apple~CloudDocs/Lavoro/PhD Bergamo/R/GitHub/Italy_AQ_AmortisedLatticeKrig"
)
A1_aq_eea_df_test <- A1_aq_eea_df[test, ]
A1_aq_eea_df_test <- A1_aq_eea_df_test[!naidx_test, ]
F3_cv_STI <- cbind(A1_aq_eea_df_test,
                   as.data.frame(Z_test),
                   yhat_test,
                   y_CI_lower,
                   y_CI_upper)

head(F3_cv_STI)

summary(F3_cv_STI$AQ_EEA_NO2 - F3_cv_STI$yhat_test)

rm(AQ_CAMS_v100_ST)
save.image(paste0("data/after_simple_Lattice_", d, ".RData"))


#ctms
Lmodel_ctms <- LatticeKrig(s_train_mod, y_train_mod, Z = Z_train_mod[, c(1:7, 9)], LKinfo = LKinfo)

objc <- Lmodel_ctms

ncov <- 8

yhat_testc <- predict(object = objc, x = s_test, Z = Z_test[, c(1:7, 9)])


# sqrt(mean((y_test - yhat_test)^2,na.rm=T))

ysehat_testc <- predictSE(object = objc, x = s_test, Z = Z_test[, c(1:7, 9)])

y_CI_lowerc <- yhat_testc - 1.96 * ysehat_testc
y_CI_upperc <- yhat_testc + 1.96 * ysehat_testc

# setwd(
#   "~/Library/Mobile Documents/com~apple~CloudDocs/Lavoro/PhD Bergamo/R/GitHub/GRINS-Spoke0-WP2"
# )
# d <- period_HPC[8]
# p <- "AQ-FRK/v.3.0.6/"
# source(paste0(p, "Z_settings/Z_settings.R"), local = TRUE)
# setwd(
#   "~/Library/Mobile Documents/com~apple~CloudDocs/Lavoro/PhD Bergamo/R/GitHub/Italy_AQ_AmortisedLatticeKrig"
# )
# A1_aq_eea_df_test <- A1_aq_eea_df[test, ]
# A1_aq_eea_df_test <- A1_aq_eea_df_test[!naidx_test, ]
# F3_cv_STI <- cbind(A1_aq_eea_df_test,
#                    as.data.frame(Z_test),
#                    yhat_test,
#                    y_CI_lower,
#                    y_CI_upper)

F3_cv_STI <- cbind(F3_cv_STI, yhat_testc, y_CI_lowerc, y_CI_upperc)

head(F3_cv_STI)

summary(F3_cv_STI$AQ_EEA_NO2 - F3_cv_STI$yhat_test)
summary(F3_cv_STI$AQ_EEA_NO2 - F3_cv_STI$yhat_testc)
dev.off()

plot(F3_cv_STI$AQ_EEA_NO2 ~ F3_cv_STI$yhat_test, col = "red")
points(F3_cv_STI$AQ_EEA_NO2 ~ F3_cv_STI$yhat_testc, col = "blue")

save.image(paste0("data/after_simple_Lattice_", d, ".RData"))

# large-scale explaining
res <- cbind(
  F3_cv_STI$AQ_EEA_NO2,
  cbind(1, Z_test[, -8]) %*% c(Lmodel$d.coef[-c(2:3)]),
  F3_cv_STI$yhat_test,
  cbind(1, Z_test[, -8]) %*% c(Lmodel_ctms$d.coef[-c(2:3)]),
  F3_cv_STI$yhat_testc
)

colnames(res) <- c(
  "measurements",
  "large-scale stationary",
  "full stationary",
  "large-scale non stationary",
  "full non stationary"
)

head(res)

#modified stationary
Lmodel2 <- LatticeKrig(
  s_train_mod,
  y_train_mod,
  Z = Z_train_mod[, c(1:7, 9)],
  nlevel = 1,
  NC = 128,
  normalize = FALSE,
  NC.buffer = 0,
  overlap = 2.5
)

obj <- Lmodel2

ncov <- 8

yhat_test2 <- predict(object = obj, x = s_test, Z = Z_test[, c(1:7, 9)])


# sqrt(mean((y_test - yhat_test)^2,na.rm=T))

ysehat_test2 <- predictSE(object = obj, x = s_test, Z = Z_test[, c(1:7, 9)])

y_CI_lower2 <- yhat_test2 - 1.96 * ysehat_test2
y_CI_upper2 <- yhat_test2 + 1.96 * ysehat_test2

F3_cv_STI <- cbind(F3_cv_STI, ysehat_test2, y_CI_lower2, y_CI_upper2)

res2 <-
  cbind(cbind(1, Z_test[, -8]) %*% c(Lmodel2$d.coef[-c(2:3)]), ysehat_test2)

colnames(res2) <- c("large-scale mod stationary", "full mod stationary")

res <- cbind(res, res2)
head(res)

LatticeKrig(s_train_mod, y_train_mod, Z = Z_train_mod[, c(-9)], LKinfo = LKinfo)
