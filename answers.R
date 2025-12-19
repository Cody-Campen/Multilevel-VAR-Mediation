library(tidyverse)
library(ggplot2)

answers_list = list()
times = c()
for(this_seed in (1:200)){
  file_name = paste0("sim_", this_seed,".RData")
  load(file_name)
  answers_list = rbind(answers_list, answers)
}
column_means = colMeans(answers_list)

# These are the parameters we care about in our mediation:
#   indirect_effect[1,1,3]
#   indirect_effect[1,3,2]
#   indirect_effect[1,6,3]
#   direct_effect[1,2]
#   direct_effect[1,3]

# ---- Bias ----
# This is not so good. 
intercepts_NIE = answers_list[,c("intercept_indirect_effect[1,1,3]_relative_bias", 
                            "X_to_intercept[1,3]_relative_bias",
                            "intercept_to_Y[1,1]_relative_bias")] |>
                  as.data.frame() |> 
                  gather(key = "variable", value = "value")

ggplot(intercepts_NIE, aes(value)) + 
  geom_histogram() +
  facet_wrap(~variable) +
  ggtitle("Intercepts indirect effect relative bias")

# This seem not quite as bad, but still not so good. 
AR_NIE = answers_list[,c("AR_indirect_effect[1,1,2]_relative_bias", 
                           "X_to_AR[1,2]_relative_bias",
                           "AR_to_Y[1,1]_relative_bias")] |>
  as.data.frame() |> 
  gather(key = "variable", value = "value")

ggplot(AR_NIE, aes(value)) + 
  geom_histogram() +
  facet_wrap(~variable) +
  ggtitle("AR indirect effect relative bias")

CR_NIE = answers_list[,c("CR_indirect_effect[1,2,3]_relative_bias", 
                           "X_to_CR[2,3]_relative_bias",
                           "CR_to_Y[1,2]_relative_bias")] |>
  as.data.frame() |> 
  gather(key = "variable", value = "value")

ggplot(CR_NIE, aes(value)) + 
  geom_histogram() +
  facet_wrap(~variable) +
  ggtitle("CR indirect effect relative bias")

DE = answers_list[,c("direct_effect[1,2]_relative_bias",
                       "direct_effect[1,3]_relative_bias")] |>
  as.data.frame() |> 
  gather(key = "variable", value = "value")

ggplot(DE, aes(value)) + 
  geom_histogram() +
  facet_wrap(~variable) +
  ggtitle("Direct effects relative bias")

# ---- Coverage ----
coverage = column_means[c("X_to_intercept[1,3]_coverage",
                         "intercept_to_Y[1,1]_coverage",
                         "X_to_AR[1,2]_coverage",
                         "AR_to_Y[1,1]_coverage",
                         "X_to_CR[2,3]_coverage",
                         "CR_to_Y[1,2]_coverage",
                         "direct_effect[1,2]_coverage",
                         "direct_effect[1,3]_coverage")]

coverage_df = data.frame(parameter = factor(names(coverage), levels = names(coverage)),
                         rate = as.numeric(coverage))

ggplot(coverage_df, aes(x = parameter, y = rate)) + 
  geom_point(shape = 15, size = 3) +
  geom_hline(yintercept = .95, linetype = "dashed", color = "red") +
  ylim(.5, 1) +
  ggtitle("Coverage Rates")

CI_widths = column_means[c("X_to_intercept[1,3]_CI_high",
                           "intercept_to_Y[1,1]_CI_high",
                           "X_to_AR[1,2]_CI_high",
                           "AR_to_Y[1,1]_CI_high",
                           "X_to_CR[2,3]_CI_high",
                           "CR_to_Y[1,2]_CI_high",
                           "direct_effect[1,2]_CI_high",
                           "direct_effect[1,3]_CI_high")] -
            column_means[c("X_to_intercept[1,3]_CI_low",
                           "intercept_to_Y[1,1]_CI_low",
                           "X_to_AR[1,2]_CI_low",
                           "AR_to_Y[1,1]_CI_low",
                           "X_to_CR[2,3]_CI_low",
                           "CR_to_Y[1,2]_CI_low",
                           "direct_effect[1,2]_CI_low",
                           "direct_effect[1,3]_CI_low")]

ESS = answers_list[, c("X_to_intercepts[1,3]_ess",
                        "M_fixed_effect[1,1]_ess",
                        "X_to_ARs[1,2]_ess",
                        "M_fixed_effect[1,3]_ess",
                        "X_to_CRs[2,3]_ess",
                        "M_fixed_effect[1,6]_ess",
                        "direct_effect[1,2]_ess",
                        "direct_effect[1,3]_ess")] |>
  as.data.frame()
plot(ESS)
