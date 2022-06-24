library(tidyverse)
library(brms)
library(gridExtra)
library(ggthemr)
ggthemr("fresh")

## we read the data

data_cues <- read_csv("../../data/data-cues.csv")
languages <- read_csv("../../data/langs.csv") ## I made this for the locations

data_cues <- left_join(data_cues, languages
                     , by = c("Language" = "name_2"))

trees <- ape::read.tree("../../data/posterior.tree")
tree <- phangorn::maxCladeCred(trees)

tree_cor <- ape::vcv.phylo(tree, corr = TRUE)

## not needed if we use brms, but just in case we change to STAN directly
tree_cor <- tree_cor[data_cues$Language, data_cues$Language]

rescale <- function(x, e = 0.0001) {
    y = ifelse(x == 1, x - e, x)
    ifelse(y == 0, y + e, y)
}

data_cues_3 <- data_cues %>%
    mutate(
        x1 = rescale(Case_Marking)
      , x2 = rescale(Tight_Semantics)
      , x3 = rescale(Rigid_Order)
      , x4 = rescale(Verb_Middle))

## no phylogenetic term

bm_1 <- brm(bf(mvbind(x1,x2,x3,x4) ~ 1
               + (1|p|gr(Language, cov = phylo)))
          , prior = c(prior(normal(0, 1), class = Intercept)
                    , prior(normal(0, 1), class = sd, resp = x1)
                    , prior(normal(0, 1), class = sd, resp = x2)
                    , prior(normal(0, 1), class = sd, resp = x3)
                    , prior(normal(0, 1), class = sd, resp = x4)
                      )
          , family = Beta
          , data = data_cues_3
          , data2 = list(phylo = tree_cor)
          , cores = 4, chains = 4
          , iter = 4000 # this model needs a lot more iters
          , warmup = 1000, seed = 1234
          , control = list(max_treedepth = 14
                         , adapt_delta = 0.99))

bm_1

## Group-Level Effects: 
## ~Language (Number of levels: 30) 
##                                Estimate Est.Error l-95% CI u-95% CI Rhat
## sd(x1_Intercept)                   1.53      0.37     0.87     2.33 1.00
## sd(x2_Intercept)                   0.25      0.15     0.02     0.56 1.00
## sd(x3_Intercept)                   1.55      0.26     1.03     2.07 1.00
## sd(x4_Intercept)                   2.79      0.43     1.93     3.58 1.00
## cor(x1_Intercept,x2_Intercept)     0.37      0.31    -0.35     0.86 1.00
## cor(x1_Intercept,x3_Intercept)    -0.74      0.13    -0.95    -0.43 1.00
## cor(x2_Intercept,x3_Intercept)    -0.16      0.31    -0.76     0.48 1.00
## cor(x1_Intercept,x4_Intercept)    -0.41      0.16    -0.69    -0.05 1.00
## cor(x2_Intercept,x4_Intercept)    -0.35      0.29    -0.84     0.34 1.00
## cor(x3_Intercept,x4_Intercept)     0.33      0.15    -0.00     0.60 1.00
##                                Bulk_ESS Tail_ESS
## sd(x1_Intercept)                   2650     5439
## sd(x2_Intercept)                   1303     3077
## sd(x3_Intercept)                   1627     2619
## sd(x4_Intercept)                   1809     1734
## cor(x1_Intercept,x2_Intercept)     7037     6895
## cor(x1_Intercept,x3_Intercept)     3228     4306
## cor(x2_Intercept,x3_Intercept)     3984     3849
## cor(x1_Intercept,x4_Intercept)     3898     5873
## cor(x2_Intercept,x4_Intercept)     1309     1573
## cor(x3_Intercept,x4_Intercept)     8290     9166

## Population-Level Effects: 
##              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
## x1_Intercept    -1.17      0.47    -2.07    -0.24 1.00     4789     7028
## x2_Intercept    -1.53      0.12    -1.77    -1.30 1.00     5786     5720
## x3_Intercept     0.34      0.45    -0.56     1.22 1.00     6395     7401
## x4_Intercept    -0.34      0.69    -1.69     1.01 1.00     9115     8483

## Family Specific Parameters: 
##        Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
## phi_x1     9.42      8.19     2.46    28.60 1.00     2330     3828
## phi_x2    73.87     48.26    28.38   207.08 1.00     1537     2723
## phi_x3    74.42     73.04     8.99   281.55 1.01      866     2295
## phi_x4    74.48     67.58     5.57   253.13 1.00     1064     1374

## Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
## and Tail_ESS are effective sample size measures, and Rhat is the potential
## scale reduction factor on split chains (at convergence, Rhat = 1).

write_rds(bm_1, "./models/bm-1.rds")

##

data_cues_2 <- data_cues %>%
    mutate(
        x1 = 1000*(Case_Marking)
      , x2 = 1000*(Tight_Semantics)
      , x3 = 1000*(Rigid_Order)
      , x4 = 1000*(Verb_Middle)
      , n = 1000)

#####
## Binomial model
## with this model the correlations are weaker but clearer

bm_a1 <- brm(bf(mvbind(x1) | trials(n) ~ 1 + (1|p|gr(Language, cov = phylo)))
          , prior = c(prior(normal(0, 1), class = Intercept)
                    , prior(normal(0, 1), class = sd, resp = x1)
                    , prior(normal(0, 1), class = sd, resp = x2)
                    , prior(normal(0, 1), class = sd, resp = x3)
                    , prior(normal(0, 1), class = sd, resp = x4)
                      )
          , family = binomial
          , data = data_cues_2
          , data2 = list(phylo = tree_cor)
          , cores = 4, chains = 4
          , iter = 4000 # this model needs a lot more iters
          , warmup = 1000, seed = 1234
          , control = list(max_treedepth = 14
                         , adapt_delta = 0.99))

bm_a1

## Group-Level Effects: 
## ~Language (Number of levels: 30) 
##                                Estimate Est.Error l-95% CI u-95% CI Rhat
## sd(x1_Intercept)                   3.33      0.38     2.68     4.15 1.00
## sd(x2_Intercept)                   0.56      0.08     0.42     0.75 1.00
## sd(x3_Intercept)                   1.70      0.19     1.37     2.12 1.00
## sd(x4_Intercept)                   3.28      0.34     2.67     4.01 1.00
## cor(x1_Intercept,x2_Intercept)     0.16      0.16    -0.16     0.46 1.00
## cor(x1_Intercept,x3_Intercept)    -0.59      0.11    -0.77    -0.36 1.00
## cor(x2_Intercept,x3_Intercept)    -0.06      0.16    -0.38     0.26 1.00
## cor(x1_Intercept,x4_Intercept)    -0.41      0.12    -0.63    -0.15 1.00
## cor(x2_Intercept,x4_Intercept)    -0.26      0.14    -0.53     0.03 1.00
## cor(x3_Intercept,x4_Intercept)     0.32      0.13     0.03     0.56 1.00
##                                Bulk_ESS Tail_ESS
## sd(x1_Intercept)                   4985     6806
## sd(x2_Intercept)                   5373     7629
## sd(x3_Intercept)                   5625     8110
## sd(x4_Intercept)                   7463     8896
## cor(x1_Intercept,x2_Intercept)     6145     7714
## cor(x1_Intercept,x3_Intercept)     7446     8214
## cor(x2_Intercept,x3_Intercept)     7270     8062
## cor(x1_Intercept,x4_Intercept)     7212     7808
## cor(x2_Intercept,x4_Intercept)     6881     8281
## cor(x3_Intercept,x4_Intercept)     6865     8136

## Population-Level Effects: 
##              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
## x1_Intercept    -1.26      0.73    -2.70     0.17 1.00     4004     6191
## x2_Intercept    -1.56      0.19    -1.94    -1.17 1.00     5142     6528
## x3_Intercept     0.23      0.46    -0.68     1.14 1.00     5722     7064
## x4_Intercept    -0.69      0.74    -2.16     0.78 1.00     9318     7829

## Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
## and Tail_ESS are effective sample size measures, and Rhat is the potential
## scale reduction factor on split chains (at convergence, Rhat = 1).

write_rds(bm_a1, "./models/bm-a1.rds")

#####
## with a GP

bm_gp_1 <- brm(bf(mvbind(x1,x2,x3,x4) ~ 1
                  + (1|p|gr(Language, cov = phylo))
                  + gp(longitude, latitude))
             , prior = c(prior(normal(0, 1), class = Intercept)
                       , prior(normal(0, 1), class = sd, resp = x1)
                       , prior(normal(0, 1), class = sd, resp = x2)
                       , prior(normal(0, 1), class = sd, resp = x3)
                       , prior(normal(0, 1), class = sd, resp = x4)
                         )
             , family = Beta
             , data = data_cues_3
             , data2 = list(phylo = tree_cor)
             , cores = 4, chains = 4
             , iter = 4000 # this model needs a lot more iters
             , warmup = 1000, seed = 1234
             , control = list(max_treedepth = 14
                            , adapt_delta = 0.99))

bm_gp_1

## Gaussian Process Terms: 
##                                Estimate Est.Error l-95% CI u-95% CI Rhat
## sdgp(x1_gplongitudelatitude)       0.80      0.58     0.04     2.07 1.00
## sdgp(x2_gplongitudelatitude)       0.39      0.47     0.02     1.53 1.00
## sdgp(x3_gplongitudelatitude)       0.75      0.54     0.04     2.01 1.00
## sdgp(x4_gplongitudelatitude)       2.60      1.22     0.27     5.19 1.00
## lscale(x1_gplongitudelatitude)     0.13      0.98     0.01     0.75 1.00
## lscale(x2_gplongitudelatitude)     0.21      0.67     0.02     1.41 1.00
## lscale(x3_gplongitudelatitude)     0.16      0.78     0.02     0.66 1.00
## lscale(x4_gplongitudelatitude)     0.13      0.22     0.03     0.37 1.00
##                                Bulk_ESS Tail_ESS
## sdgp(x1_gplongitudelatitude)       1469     4595
## sdgp(x2_gplongitudelatitude)       1548     1239
## sdgp(x3_gplongitudelatitude)       3185     4411
## sdgp(x4_gplongitudelatitude)       2093     1883
## lscale(x1_gplongitudelatitude)     2103     2670
## lscale(x2_gplongitudelatitude)     1393     1203
## lscale(x3_gplongitudelatitude)     2183     2686
## lscale(x4_gplongitudelatitude)     2132     2500

## Group-Level Effects: 
## ~Language (Number of levels: 30) 
##                                Estimate Est.Error l-95% CI u-95% CI Rhat
## sd(x1_Intercept)                   1.56      0.37     0.87     2.33 1.00
## sd(x2_Intercept)                   0.21      0.14     0.01     0.52 1.00
## sd(x3_Intercept)                   1.40      0.27     0.88     1.94 1.00
## sd(x4_Intercept)                   2.11      0.55     0.99     3.15 1.00
## cor(x1_Intercept,x2_Intercept)     0.32      0.35    -0.51     0.87 1.00
## cor(x1_Intercept,x3_Intercept)    -0.78      0.14    -0.97    -0.43 1.00
## cor(x2_Intercept,x3_Intercept)    -0.19      0.36    -0.82     0.58 1.00
## cor(x1_Intercept,x4_Intercept)    -0.48      0.21    -0.83    -0.02 1.00
## cor(x2_Intercept,x4_Intercept)    -0.24      0.37    -0.85     0.57 1.00
## cor(x3_Intercept,x4_Intercept)     0.43      0.20     0.00     0.77 1.00
##                                Bulk_ESS Tail_ESS
## sd(x1_Intercept)                   3846     6787
## sd(x2_Intercept)                   2258     5183
## sd(x3_Intercept)                   3434     4838
## sd(x4_Intercept)                   1277     1296
## cor(x1_Intercept,x2_Intercept)     9965     7985
## cor(x1_Intercept,x3_Intercept)     4149     7367
## cor(x2_Intercept,x3_Intercept)     4919     6706
## cor(x1_Intercept,x4_Intercept)     2395     5311
## cor(x2_Intercept,x4_Intercept)     2551     4531
## cor(x3_Intercept,x4_Intercept)     7921     8275

## Population-Level Effects: 
##              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
## x1_Intercept    -1.20      0.58    -2.28     0.03 1.00     4889     6053
## x2_Intercept    -1.45      0.33    -1.84    -0.47 1.00     2808     1062
## x3_Intercept     0.39      0.52    -0.62     1.41 1.00     8927     7445
## x4_Intercept    -0.08      0.79    -1.61     1.50 1.00    12719     8932

## Family Specific Parameters: 
##        Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
## phi_x1    21.15     27.07     3.31    98.33 1.00      805     2132
## phi_x2   100.85     58.63    35.05   252.89 1.00     2292     5285
## phi_x3    77.17     66.94    11.92   263.64 1.00     1331     2986
## phi_x4    97.50     86.56     6.41   323.47 1.00     1005     1129

## Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
## and Tail_ESS are effective sample size measures, and Rhat is the potential
## scale reduction factor on split chains (at convergence, Rhat = 1).
## Warning message:
## There were 6 divergent transitions after warmup. Increasing adapt_delta above 0.99 may help. See http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup 

write_rds(bm_gp_1, "./models/bm-gp-1.rds")

##

bm_gp_a1 <- brm(bf(mvbind(x1,x2,x3,x4) | trials(n) ~ 1
                  + (1|p|gr(Language, cov = phylo))
                  + gp(longitude, latitude))
             , prior = c(prior(normal(0, 1), class = Intercept)
                       , prior(normal(0, 1), class = sd, resp = x1)
                       , prior(normal(0, 1), class = sd, resp = x2)
                       , prior(normal(0, 1), class = sd, resp = x3)
                       , prior(normal(0, 1), class = sd, resp = x4)
                         )
             , family = binomial
             , data = data_cues_2
             , data2 = list(phylo = tree_cor)
             , cores = 4, chains = 4
             , iter = 4000 # this model needs a lot more iters
             , warmup = 1000, seed = 1234
             , control = list(max_treedepth = 14
                            , adapt_delta = 0.99))

## slightly better, could be fixable Increasing adapt_delta

bm_gp_a1

## Gaussian Process Terms: 
##                                Estimate Est.Error l-95% CI u-95% CI Rhat
## sdgp(x1_gplongitudelatitude)       2.13      0.71     1.07     3.84 1.00
## sdgp(x2_gplongitudelatitude)       0.36      0.41     0.02     1.35 1.00
## sdgp(x3_gplongitudelatitude)       0.72      0.45     0.06     1.70 1.00
## sdgp(x4_gplongitudelatitude)       3.74      1.41     0.80     6.83 1.00
## lscale(x1_gplongitudelatitude)     0.03      0.01     0.01     0.05 1.00
## lscale(x2_gplongitudelatitude)     0.20      0.78     0.01     1.35 1.01
## lscale(x3_gplongitudelatitude)     0.20      5.73     0.02     0.47 1.00
## lscale(x4_gplongitudelatitude)     0.11      0.12     0.04     0.28 1.00
##                                Bulk_ESS Tail_ESS
## sdgp(x1_gplongitudelatitude)       2881     4629
## sdgp(x2_gplongitudelatitude)       2160     1911
## sdgp(x3_gplongitudelatitude)       3899     5955
## sdgp(x4_gplongitudelatitude)       2491     1175
## lscale(x1_gplongitudelatitude)     4663     5806
## lscale(x2_gplongitudelatitude)     1175     1983
## lscale(x3_gplongitudelatitude)     3181     2706
## lscale(x4_gplongitudelatitude)     2000     2050

## Group-Level Effects: 
## ~Language (Number of levels: 30) 
##                                Estimate Est.Error l-95% CI u-95% CI Rhat
## sd(x1_Intercept)                   2.04      0.52     0.99     3.06 1.00
## sd(x2_Intercept)                   0.46      0.12     0.20     0.69 1.00
## sd(x3_Intercept)                   1.54      0.24     1.08     2.03 1.00
## sd(x4_Intercept)                   2.26      0.52     1.25     3.34 1.00
## cor(x1_Intercept,x2_Intercept)     0.36      0.23    -0.13     0.77 1.00
## cor(x1_Intercept,x3_Intercept)    -0.74      0.15    -0.95    -0.38 1.00
## cor(x2_Intercept,x3_Intercept)    -0.05      0.22    -0.48     0.39 1.00
## cor(x1_Intercept,x4_Intercept)    -0.56      0.21    -0.87    -0.05 1.00
## cor(x2_Intercept,x4_Intercept)    -0.20      0.23    -0.64     0.27 1.00
## cor(x3_Intercept,x4_Intercept)     0.42      0.18     0.02     0.73 1.00
##                                Bulk_ESS Tail_ESS
## sd(x1_Intercept)                   4258     3672
## sd(x2_Intercept)                   1634     1796
## sd(x3_Intercept)                   5374     5887
## sd(x4_Intercept)                   1907     2261
## cor(x1_Intercept,x2_Intercept)     4308     5715
## cor(x1_Intercept,x3_Intercept)     3999     6127
## cor(x2_Intercept,x3_Intercept)     7425     4828
## cor(x1_Intercept,x4_Intercept)     3342     4496
## cor(x2_Intercept,x4_Intercept)     7845     7159
## cor(x3_Intercept,x4_Intercept)     9096     8425

## Population-Level Effects: 
##              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
## x1_Intercept    -1.62      0.70    -2.97    -0.21 1.00    10747     8952
## x2_Intercept    -1.44      0.33    -1.87    -0.62 1.00     4679     2103
## x3_Intercept     0.26      0.51    -0.75     1.28 1.00    11621     8509
## x4_Intercept    -0.20      0.86    -1.88     1.50 1.00    12855     9014

## Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
## and Tail_ESS are effective sample size measures, and Rhat is the potential
## scale reduction factor on split chains (at convergence, Rhat = 1).
## Warning message:
## There were 2 divergent transitions after warmup. Increasing adapt_delta above 0.99 may help. See http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup

write_rds(bm_gp_a1, "./models/bm-gp-a1.rds")

########################################
########################################

#####
## comparing correlations Beta

bm_gp_1

## no gp

##                                Estimate Est.Error l-95% CI u-95% CI Rhat
## cor(x1_Intercept,x2_Intercept)     0.37      0.31    -0.35     0.86 1.00
## cor(x1_Intercept,x3_Intercept)    -0.74      0.13    -0.95    -0.43 1.00
## cor(x2_Intercept,x3_Intercept)    -0.16      0.31    -0.76     0.48 1.00
## cor(x1_Intercept,x4_Intercept)    -0.41      0.16    -0.69    -0.05 1.00
## cor(x2_Intercept,x4_Intercept)    -0.35      0.29    -0.84     0.34 1.00
## cor(x3_Intercept,x4_Intercept)     0.33      0.15    -0.00     0.60 1.00

## with gp:

##                                Estimate Est.Error l-95% CI u-95% CI Rhat
## cor(x1_Intercept,x2_Intercept)     0.32      0.35    -0.51     0.87 1.00
## cor(x1_Intercept,x3_Intercept)    -0.78      0.14    -0.97    -0.43 1.00
## cor(x2_Intercept,x3_Intercept)    -0.19      0.36    -0.82     0.58 1.00
## cor(x1_Intercept,x4_Intercept)    -0.48      0.21    -0.83    -0.02 1.00
## cor(x2_Intercept,x4_Intercept)    -0.24      0.37    -0.85     0.57 1.00
## cor(x3_Intercept,x4_Intercept)     0.43      0.20     0.00     0.77 1.00

## the spatial component has very little effect on the correlations

#####
## comparing correlations binomial

## no gp

##                                Estimate Est.Error l-95% CI u-95% CI Rhat
## cor(x1_Intercept,x2_Intercept)     0.16      0.16    -0.16     0.46 1.00
## cor(x1_Intercept,x3_Intercept)    -0.59      0.11    -0.77    -0.36 1.00
## cor(x2_Intercept,x3_Intercept)    -0.06      0.16    -0.38     0.26 1.00
## cor(x1_Intercept,x4_Intercept)    -0.41      0.12    -0.63    -0.15 1.00
## cor(x2_Intercept,x4_Intercept)    -0.26      0.14    -0.53     0.03 1.00
## cor(x3_Intercept,x4_Intercept)     0.32      0.13     0.03     0.56 1.00

## with gp:

##                                Estimate Est.Error l-95% CI u-95% CI Rhat
## cor(x1_Intercept,x2_Intercept)     0.36      0.23    -0.13     0.77 1.00
## cor(x1_Intercept,x3_Intercept)    -0.74      0.15    -0.95    -0.38 1.00
## cor(x2_Intercept,x3_Intercept)    -0.05      0.22    -0.48     0.39 1.00
## cor(x1_Intercept,x4_Intercept)    -0.56      0.21    -0.87    -0.05 1.00
## cor(x2_Intercept,x4_Intercept)    -0.20      0.23    -0.64     0.27 1.00
## cor(x3_Intercept,x4_Intercept)     0.42      0.18     0.02     0.73 1.00

## this is super-weird, the spatial component strenghtens the correlations


########################################
## kfold cv
########################################

bm_a1_k <- kfold(bm_a1, K = 10
              , cores = 4, silent = TRUE
              , open_progress=FALSE)

write_rds(bm_a1_k, "./models/bm-a1-k.rds")

bm_gp_a1_k <- kfold(bm_gp_a1, K = 10
                 , cores = 4, silent = TRUE
                 , open_progress=FALSE)

write_rds(bm_gp_a1_k, "./models/bm-gp-a1-k.rds")

#####
## comparing GP vs no GP

loo_compare(bm_gp_a1_k, bm_a1_k)

##          elpd_diff se_diff
## bm_gp_a1     0.0       0.0
## bm_a1    -1945.9     528.6

loo_compare(bm_gp_1_k, bm_1_k)

##         elpd_diff se_diff
## bm_gp_1  0.0       0.0
## bm_1    -6.2       8.6

## this is a very weird result
## with the beta model, there is no effect from the GP
## while with the binomial model there is a massive effect
## more at the end

######
## comparing fit with Bayesian R2

apply(bayes_R2(bm_gp_a1), 2, round, 2)

##      Estimate Est.Error Q2.5 Q97.5
## R2x1     1.00      0.00 1.00  1.00
## R2x2     0.95      0.01 0.92  0.97
## R2x3     1.00      0.00 1.00  1.00
## R2x4     1.00      0.00 1.00  1.00

apply(bayes_R2(bm_gp_1), 2, round, 2)

##      Estimate Est.Error Q2.5 Q97.5
## R2x1     0.85      0.12 0.54  0.99
## R2x2     0.53      0.21 0.08  0.85
## R2x3     0.95      0.05 0.83  0.99
## R2x4     0.98      0.05 0.88  1.00

########################################
########################################

bm_13_24 <- brm(bf(mvbind(x1,x3) ~ 1 + x2 + x4
               + (1|p|gr(Language, cov = phylo)))
             , prior = c(prior(normal(0, 1), class = Intercept)
                       , prior(normal(0, 1), class = b, resp = x1)
                       , prior(normal(0, 1), class = b, resp = x3)
                       , prior(normal(0, 1), class = sd, resp = x1)
                       , prior(normal(0, 1), class = sd, resp = x3)
                      )
          , family = Beta
          , data = data_cues_3
          , data2 = list(phylo = tree_cor)
          , cores = 4, chains = 4
          , iter = 4000 # this model needs a lot more iters
          , warmup = 1000, seed = 1234
          , control = list(max_treedepth = 14
                         , adapt_delta = 0.99))

bm_13_24

## Group-Level Effects: 
## ~Language (Number of levels: 30) 
##                                Estimate Est.Error l-95% CI u-95% CI Rhat
## sd(x1_Intercept)                   1.37      0.35     0.76     2.13 1.00
## sd(x3_Intercept)                   1.48      0.26     0.97     2.01 1.00
## cor(x1_Intercept,x3_Intercept)    -0.79      0.15    -0.99    -0.45 1.00
##                                Bulk_ESS Tail_ESS
## sd(x1_Intercept)                   2841     5302
## sd(x3_Intercept)                   1420     2650
## cor(x1_Intercept,x3_Intercept)     2273     4149

## Population-Level Effects: 
##              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
## x1_Intercept    -0.70      0.60    -1.86     0.51 1.00     7778     8454
## x3_Intercept    -0.16      0.58    -1.30     0.97 1.00     8470     8343
## x1_x2            0.31      0.97    -1.61     2.20 1.00    22261     8267
## x1_x4           -0.91      0.64    -2.17     0.35 1.00     8387     9368
## x3_x2            0.10      0.97    -1.82     2.00 1.00    22778     8008
## x3_x4            0.83      0.59    -0.36     1.97 1.00     8576     8721

## Family Specific Parameters: 
##        Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
## phi_x1     7.57      5.94     2.30    22.13 1.00     2477     4450
## phi_x3    66.73     64.07     8.61   242.93 1.00      847     2037

## Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
## and Tail_ESS are effective sample size measures, and Rhat is the potential
## scale reduction factor on split chains (at convergence, Rhat = 1).

##

bm_14_23 <- brm(bf(mvbind(x1,x4) ~ 1 + x2 + x3
               + (1|p|gr(Language, cov = phylo)))
             , prior = c(prior(normal(0, 1), class = Intercept)
                       , prior(normal(0, 1), class = b, resp = x1)
                       , prior(normal(0, 1), class = b, resp = x4)
                       , prior(normal(0, 1), class = sd, resp = x1)
                       , prior(normal(0, 1), class = sd, resp = x4)
                      )
          , family = Beta
          , data = data_cues_3
          , data2 = list(phylo = tree_cor)
          , cores = 4, chains = 4
          , iter = 4000 # this model needs a lot more iters
          , warmup = 1000, seed = 1234
          , control = list(max_treedepth = 14
                         , adapt_delta = 0.99))

bm_14_23

## Group-Level Effects: 
## ~Language (Number of levels: 30) 
##                                Estimate Est.Error l-95% CI u-95% CI Rhat
## sd(x1_Intercept)                   0.89      0.44     0.10     1.79 1.00
## sd(x4_Intercept)                   2.70      0.44     1.77     3.52 1.00
## cor(x1_Intercept,x4_Intercept)    -0.40      0.27    -0.88     0.19 1.00
##                                Bulk_ESS Tail_ESS
## sd(x1_Intercept)                   1372     2767
## sd(x4_Intercept)                   1479     1756
## cor(x1_Intercept,x4_Intercept)     1157     1460

## Population-Level Effects: 
##              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
## x1_Intercept     0.07      0.59    -1.13     1.20 1.00     8413     8926
## x4_Intercept    -0.71      0.89    -2.48     0.99 1.00     9096     8876
## x1_x2            0.27      0.98    -1.66     2.16 1.00    20708     9098
## x1_x3           -2.24      0.75    -3.66    -0.71 1.00     9553     9046
## x4_x2           -0.26      1.00    -2.24     1.68 1.00    21611     8248
## x4_x3            0.80      0.88    -0.92     2.48 1.00    12037     9020

## Family Specific Parameters: 
##        Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
## phi_x1     6.28      5.04     1.78    19.33 1.00     1623     3778
## phi_x4    67.98     64.15     4.70   240.22 1.00      882     1497

## Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
## and Tail_ESS are effective sample size measures, and Rhat is the potential
## scale reduction factor on split chains (at convergence, Rhat = 1).

##


bm_24_13 <- brm(bf(mvbind(x2,x4) ~ 1 + x1 + x3
               + (1|p|gr(Language, cov = phylo)))
             , prior = c(prior(normal(0, 1), class = Intercept)
                       , prior(normal(0, 1), class = b, resp = x2)
                       , prior(normal(0, 1), class = b, resp = x4)
                       , prior(normal(0, 1), class = sd, resp = x2)
                       , prior(normal(0, 1), class = sd, resp = x4)
                      )
          , family = Beta
          , data = data_cues_3
          , data2 = list(phylo = tree_cor)
          , cores = 4, chains = 4
          , iter = 4000 # this model needs a lot more iters
          , warmup = 1000, seed = 1234
          , control = list(max_treedepth = 14
                         , adapt_delta = 0.99))

bm_24_13

## Group-Level Effects: 
## ~Language (Number of levels: 30) 
##                                Estimate Est.Error l-95% CI u-95% CI Rhat
## sd(x2_Intercept)                   0.25      0.15     0.01     0.56 1.00
## sd(x4_Intercept)                   2.63      0.44     1.71     3.46 1.00
## cor(x2_Intercept,x4_Intercept)    -0.38      0.34    -0.93     0.52 1.01
##                                Bulk_ESS Tail_ESS
## sd(x2_Intercept)                   1202     2722
## sd(x4_Intercept)                   1560     2100
## cor(x2_Intercept,x4_Intercept)      807      623

## Population-Level Effects: 
##              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
## x2_Intercept    -1.78      0.34    -2.45    -1.12 1.00     8775     8214
## x4_Intercept    -0.49      0.91    -2.34     1.29 1.00     9533     8213
## x2_x1            0.50      0.40    -0.30     1.28 1.00     8548     8022
## x2_x3            0.16      0.38    -0.59     0.92 1.00     9334     8494
## x4_x1           -0.92      0.89    -2.66     0.82 1.00    11410     9225
## x4_x3            0.83      0.89    -0.94     2.57 1.00    12799     9388

## Family Specific Parameters: 
##        Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
## phi_x2    77.79     51.70    29.18   220.41 1.00     1452     2414
## phi_x4    66.10     62.37     4.19   234.97 1.00      975     1766

## Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
## and Tail_ESS are effective sample size measures, and Rhat is the potential
## scale reduction factor on split chains (at convergence, Rhat = 1).

########################################
## binomial
########################################

## these results are more striking
## particularly, these models find coefficients of 0 for the variables
## we're conditioning on
## I think here independence is much more clear

bm_13_24a <- brm(bf(mvbind(x1,x3) | trials(n) ~ 1 + x2 + x4
               + (1|p|gr(Language, cov = phylo)))
             , prior = c(prior(normal(0, 1), class = Intercept)
                       , prior(normal(0, 1), class = b, resp = x1)
                       , prior(normal(0, 1), class = b, resp = x3)
                       , prior(normal(0, 1), class = sd, resp = x1)
                       , prior(normal(0, 1), class = sd, resp = x3)
                      )
          , family = binomial
          , data = data_cues_2
          , data2 = list(phylo = tree_cor)
          , cores = 4, chains = 4
          , iter = 4000 # this model needs a lot more iters
          , warmup = 1000, seed = 1234
          , control = list(max_treedepth = 14
                         , adapt_delta = 0.99))

bm_13_24a

## Group-Level Effects: 
## ~Language (Number of levels: 30) 
##                                Estimate Est.Error l-95% CI u-95% CI Rhat
## sd(x1_Intercept)                   3.11      0.38     2.43     3.93 1.00
## sd(x3_Intercept)                   1.67      0.20     1.32     2.13 1.00
## cor(x1_Intercept,x3_Intercept)    -0.57      0.12    -0.76    -0.31 1.00
##                                Bulk_ESS Tail_ESS
## sd(x1_Intercept)                   6593     8377
## sd(x3_Intercept)                   6051     7806
## cor(x1_Intercept,x3_Intercept)     7038     8432

## Population-Level Effects: 
##              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
## x1_Intercept     1.23      2.14    -2.93     5.59 1.00     6382     7414
## x3_Intercept    -0.92      1.13    -3.16     1.30 1.00     7437     8049
## x1_x2            0.00      0.01    -0.01     0.02 1.00     5650     7194
## x1_x4           -0.00      0.00    -0.01    -0.00 1.00     6493     7254
## x3_x2            0.00      0.00    -0.01     0.01 1.00     6826     7535
## x3_x4            0.00      0.00     0.00     0.00 1.00     7668     8732

## Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
## and Tail_ESS are effective sample size measures, and Rhat is the potential
## scale reduction factor on split chains (at convergence, Rhat = 1).

##

bm_14_23a <- brm(bf(mvbind(x1,x4) | trials(n) ~ 1 + x2 + x3
               + (1|p|gr(Language, cov = phylo)))
             , prior = c(prior(normal(0, 1), class = Intercept)
                       , prior(normal(0, 1), class = b, resp = x1)
                       , prior(normal(0, 1), class = b, resp = x4)
                       , prior(normal(0, 1), class = sd, resp = x1)
                       , prior(normal(0, 1), class = sd, resp = x4)
                      )
          , family = binomial
          , data = data_cues_2
          , data2 = list(phylo = tree_cor)
          , cores = 4, chains = 4
          , iter = 4000 # this model needs a lot more iters
          , warmup = 1000, seed = 1234
          , control = list(max_treedepth = 14
                         , adapt_delta = 0.99))

bm_14_23a

## Group-Level Effects: 
## ~Language (Number of levels: 30) 
##                                Estimate Est.Error l-95% CI u-95% CI Rhat
## sd(x1_Intercept)                   2.65      0.36     2.02     3.44 1.00
## sd(x4_Intercept)                   3.05      0.34     2.44     3.78 1.00
## cor(x1_Intercept,x4_Intercept)    -0.32      0.14    -0.57    -0.03 1.00
##                                Bulk_ESS Tail_ESS
## sd(x1_Intercept)                   7457     9164
## sd(x4_Intercept)                   6914     7939
## cor(x1_Intercept,x4_Intercept)     6521     8451

## Population-Level Effects: 
##              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
## x1_Intercept     2.68      1.83    -0.88     6.35 1.00     5647     6797
## x4_Intercept    -0.29      2.01    -4.15     3.65 1.00     5910     7246
## x1_x2            0.01      0.01    -0.01     0.02 1.00     6120     6992
## x1_x3           -0.01      0.00    -0.01    -0.01 1.00     5346     6787
## x4_x2           -0.02      0.01    -0.03    -0.00 1.00     6037     7208
## x4_x3            0.00      0.00     0.00     0.01 1.00     6000     7036

## Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
## and Tail_ESS are effective sample size measures, and Rhat is the potential
## scale reduction factor on split chains (at convergence, Rhat = 1).

##

bm_24_13a <- brm(bf(mvbind(x2,x4) | trials(n) ~ 1 + x1 + x3
               + (1|p|gr(Language, cov = phylo)))
             , prior = c(prior(normal(0, 1), class = Intercept)
                       , prior(normal(0, 1), class = b, resp = x2)
                       , prior(normal(0, 1), class = b, resp = x4)
                       , prior(normal(0, 1), class = sd, resp = x2)
                       , prior(normal(0, 1), class = sd, resp = x4)
                      )
          , family = binomial
          , data = data_cues_2
          , data2 = list(phylo = tree_cor)
          , cores = 4, chains = 4
          , iter = 4000 # this model needs a lot more iters
          , warmup = 1000, seed = 1234
          , control = list(max_treedepth = 14
                         , adapt_delta = 0.99))

bm_24_13a


## Group-Level Effects: 
## ~Language (Number of levels: 30) 
##                                Estimate Est.Error l-95% CI u-95% CI Rhat
## sd(x2_Intercept)                   0.56      0.08     0.42     0.75 1.00
## sd(x4_Intercept)                   3.19      0.36     2.56     3.96 1.00
## cor(x2_Intercept,x4_Intercept)    -0.25      0.15    -0.53     0.07 1.00
##                                Bulk_ESS Tail_ESS
## sd(x2_Intercept)                   4123     5850
## sd(x4_Intercept)                   5862     7064
## cor(x2_Intercept,x4_Intercept)     4930     5960

## Population-Level Effects: 
##              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
## x2_Intercept    -1.79      0.47    -2.72    -0.85 1.00     7949     7310
## x4_Intercept    -1.08      2.64    -6.20     4.24 1.00     6373     7183
## x2_x1            0.00      0.00    -0.00     0.00 1.00     8162     7273
## x2_x3            0.00      0.00    -0.00     0.00 1.00     8124     7471
## x4_x1           -0.00      0.00    -0.01     0.00 1.00     5900     7134
## x4_x3            0.00      0.00    -0.00     0.01 1.00     6304     6552

## Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
## and Tail_ESS are effective sample size measures, and Rhat is the potential
## scale reduction factor on split chains (at convergence, Rhat = 1).

########################################
## plots
########################################

## recall it seems like the binomial model
## fits much better

apply(bayes_R2(bm_gp_a1), 2, round, 2)

##      Estimate Est.Error Q2.5 Q97.5
## R2x1     1.00      0.00 1.00  1.00
## R2x2     0.95      0.01 0.92  0.97
## R2x3     1.00      0.00 1.00  1.00
## R2x4     1.00      0.00 1.00  1.00

apply(bayes_R2(bm_gp_1), 2, round, 2)

##      Estimate Est.Error Q2.5 Q97.5
## R2x1     0.85      0.12 0.54  0.99
## R2x2     0.53      0.21 0.08  0.85
## R2x3     0.95      0.05 0.83  0.99
## R2x4     0.98      0.05 0.88  1.00


bm_1 <- read_rds("./bm-1.rds")
bm_a1 <- read_rds("./bm-a1.rds")
bm_gp_1 <- read_rds("./bm-gp-1.rds")
bm_gp_a1 <- read_rds("./bm-gp-a1.rds")


pp_1_x1 <- pp_check(bm_1, resp = "x1", ndraws = 100)
pp_1_x2 <- pp_check(bm_1, resp = "x2", ndraws = 100)
pp_1_x3 <- pp_check(bm_1, resp = "x3", ndraws = 100)
pp_1_x4 <- pp_check(bm_1, resp = "x4", ndraws = 100)

pp_a1_x1 <- pp_check(bm_a1, resp = "x1", ndraws = 100)
pp_a1_x2 <- pp_check(bm_a1, resp = "x2", ndraws = 100)
pp_a1_x3 <- pp_check(bm_a1, resp = "x3", ndraws = 100)
pp_a1_x4 <- pp_check(bm_a1, resp = "x4", ndraws = 100)

## the binomial model does look better
## even without GP

x1_plot <- grid.arrange(pp_1_x1+ggtitle("Beta")
                     , pp_a1_x1+ggtitle("Binomial")
                     , ncol = 2)

ggsave("beta-vs-binom-x1.pdf", x1_plot
     , height = 4, width = 6)

x2_plot <- grid.arrange(pp_1_x2+ggtitle("Beta")
                     , pp_a1_x2+ggtitle("Binomial")
                     , ncol = 2)

ggsave("beta-vs-binom-x2.pdf", x2_plot
     , height = 4, width = 6)

x3_plot <- grid.arrange(pp_1_x3+ggtitle("Beta")
                     , pp_a1_x3+ggtitle("Binomial")
                     , ncol = 2)

ggsave("beta-vs-binom-x3.pdf", x3_plot
     , height = 4, width = 6)


x4_plot <- grid.arrange(pp_1_x4+ggtitle("Beta")
                     , pp_a1_x4+ggtitle("Binomial")
                     , ncol = 2)

ggsave("beta-vs-binom-x4.pdf", x4_plot
     , height = 4, width = 6)

## Binomial with and without GP
## in terms of posterior predictive checks
## there is no apparent difference

pp_gpa1_x1 <- pp_check(bm_gp_a1, resp = "x1", ndraws = 100)
pp_gpa1_x2 <- pp_check(bm_gp_a1, resp = "x2", ndraws = 100)
pp_gpa1_x3 <- pp_check(bm_gp_a1, resp = "x3", ndraws = 100)
pp_gpa1_x4 <- pp_check(bm_gp_a1, resp = "x4", ndraws = 100)

x1_plot <- grid.arrange(pp_gpa1_x1+ggtitle("GP")
                     , pp_a1_x1+ggtitle("No GP")
                     , ncol = 2)

ggsave("gp-vs-nogp-x1.pdf", x1_plot
     , height = 4, width = 6)

x2_plot <- grid.arrange(pp_gpa1_x2+ggtitle("GP")
                     , pp_a1_x2+ggtitle("No GP")
                     , ncol = 2)

ggsave("gp-vs-nogp-x2.pdf", x2_plot
     , height = 4, width = 6)

x3_plot <- grid.arrange(pp_gpa1_x3+ggtitle("GP")
                     , pp_a1_x3+ggtitle("No GP")
                     , ncol = 2)

ggsave("gp-vs-nogp-x3.pdf", x3_plot
     , height = 4, width = 6)


x4_plot <- grid.arrange(pp_gpa1_x4+ggtitle("GP")
                     , pp_a1_x4+ggtitle("No GP")
                     , ncol = 2)

ggsave("gp-vs-nogp-x4.pdf", x4_plot
     , height = 4, width = 6)

#####
## let's see the phylogenetic effects

## first thing to look at is the sd of the phylogenetic effects

bm_a1

##                                Estimate Est.Error l-95% CI u-95% CI Rhat
## sd(x1_Intercept)                   3.33      0.38     2.68     4.15 1.00
## sd(x2_Intercept)                   0.56      0.08     0.42     0.75 1.00
## sd(x3_Intercept)                   1.70      0.19     1.37     2.12 1.00
## sd(x4_Intercept)                   3.28      0.34     2.67     4.01 1.00

bm_gp_a1

##                                Estimate Est.Error l-95% CI u-95% CI Rhat
## sd(x1_Intercept)                   2.04      0.52     0.99     3.06 1.00
## sd(x2_Intercept)                   0.46      0.12     0.20     0.69 1.00
## sd(x3_Intercept)                   1.54      0.24     1.08     2.03 1.00
## sd(x4_Intercept)                   2.26      0.52     1.25     3.34 1.00

## adding a GP reduces the sd noticeably, especially for x1 and x4
## this makes sence, since x1 and x4 have some clear spatial autocor

space_plot <-
    data_cues %>% pivot_longer(cols = c(Case_Marking, Tight_Semantics
                                  , Rigid_Order, Verb_Middle)
                         , names_to = "term") %>% 
    ggplot(aes(x = longitude, y = latitude
             , color = value)) +
    geom_point()+
    facet_wrap(~term)

ggsave("./spatial-plot.pdf", height = 5, width = 8)

## we can look at the coefficients themselves

re_a1 <- ranef(bm_a1)[[1]]
re_gp_a1 <- ranef(bm_gp_a1)[[1]]

re_a1_df <- lapply(1:4, function(j) {
    mutate(as.data.frame(re_a1[,,j]) %>%
         rownames_to_column("language"), resp = str_c("x", j))
}) %>% bind_rows()%>%
    mutate(model = "no gp")

re_gp_a1_df <- lapply(1:4, function(j) {
    mutate(as.data.frame(re_gp_a1[,,j]) %>%
         rownames_to_column("language"), resp = str_c("x", j))
}) %>% bind_rows() %>%
    mutate(model = "gp")

re_plot_1 <- rbind(re_a1_df, re_gp_a1_df) %>%
    ggplot(aes(x = Estimate, y = language, color = model))+
    geom_point()+
    facet_wrap( ~ resp)

ggsave("re-plot-gp-nogp.pdf", re_plot_1, height = 8, width = 10)

## now it's very clear that for X1 and X4 there is a lot of spatial cor
## which is captured by the GP
## this is why the phylogenetic estimates are much more extreme
## for the model without a GP

## The Beta model also picks up on this for X4 but not X1
## I am not clear on why

re_1 <- ranef(bm_1)[[1]]
re_gp_1 <- ranef(bm_gp_1)[[1]]

re_1_df <- lapply(1:4, function(j) {
    mutate(as.data.frame(re_1[,,j]) %>%
         rownames_to_column("language"), resp = str_c("x", j))
}) %>% bind_rows()%>%
    mutate(model = "no gp")

re_gp_1_df <- lapply(1:4, function(j) {
    mutate(as.data.frame(re_gp_1[,,j]) %>%
         rownames_to_column("language"), resp = str_c("x", j))
}) %>% bind_rows() %>%
    mutate(model = "gp")

re_plot_1 <- rbind(re_1_df, re_gp_1_df) %>%
    ggplot(aes(x = Estimate, y = language, color = model))+
    geom_point()+
    facet_wrap( ~ resp)

ggsave("re-beta-plot-gp-nogp.pdf", re_plot_1, height = 8, width = 10)
