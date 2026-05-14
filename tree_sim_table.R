library(dplyr)
results <- readRDS("~/Desktop/mich_tree/results/results.rds")

results %>% 
  group_by(n, L, method) %>% 
  summarize(bias = mean(L - L_est), 
            #hausdorff = mean(hausdorff_1 + hausdorff_2, na.rm = TRUE),
            ARI = mean(ARI, na.rm = TRUE),
            AMI = mean(AMI, na.rm = TRUE),
            mean_mse = mean(mean_mse, na.rm = TRUE),
            #ci_length = sum(L_est * avg_len, na.rm = TRUE) / sum(L_est, na.rm = TRUE),
            coverage = sum(n_covered) / sum(n_detected),
            time = mean(time, na.rm = TRUE))
