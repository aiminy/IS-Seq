install.packages("yardstick")
install.packages("tidymodels")


library(yardstick)
library(dplyr)
library(ggplot2)

data(two_class_example)


pr_curve(two_class_example, truth, Class1) %>%
  ggplot(aes(x = recall, y = precision)) +
  geom_path() +
  coord_equal() +
  theme_bw()

df <- pr_curve(two_class_example, truth, Class1)

temp <- df[1,]

temp$recall <- 1

temp %>% ggplot(aes(x = recall, y = precision)) 

temp %>% ggplot() +
  geom_path(aes(recall, precision)) +
  coord_equal()

library(tidymodels)
theme_set(theme_bw())

true_class score

print(dat)

dat <- as_tibble(data.frame(true_class=c(rep("IS",1),rep("Non_IS",1)),score=c(100,0)))

dat <- dat %>%
  mutate(pred_class = if_else(score > 0, "IS", "Non_IS")) %>%
  mutate_at(vars(true_class, pred_class), list(. %>% factor() %>% forcats::fct_relevel("IS")))

dat %>% sens(true_class, pred_class)

dat %>% spec(true_class, pred_class)

dat %>% precision(true_class, pred_class)

dat %>% recall(true_class, pred_class)

pr_dat <- dat %>%
  pr_curve(true_class, score)

pr_dat %>%
  arrange(.threshold) %>% 
  ggplot() +
  geom_path(aes(recall, precision)) +
  coord_equal()

pr_dat_1 <- pr_dat[2,] 

pr_dat_1 <- pr_dat[2,] 

pr_dat_2 <- pr_dat_1[rep(seq_len(nrow(pr_dat_1)), each = 4), ]

pr_dat_2$Method <- c("Read","Umi","Fragment","INSPIIRED")

pr_dat_2 %>% ggplot(aes(recall, precision,color=Method)) + geom_point()+xlim(0,1)+ylim(0,1)

roc_dat <- dat %>% roc_curve(true_class, score)

roc_dat %>%
  arrange(.threshold) %>% 
  ggplot() +
  geom_path(aes(1 - specificity, sensitivity)) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") + 
  coord_equal()

roc_dat_1 <- roc_dat[3,]

roc_dat_1 %>%
  arrange(.threshold) %>% 
  ggplot() +
  geom_path(aes(1 - specificity, sensitivity)) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") + 
  coord_equal()

roc_dat_1 %>%
  arrange(.threshold) %>% 
  ggplot() +
  geom_path(aes(1 - specificity, sensitivity)) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") + 
  coord_equal()

roc_dat_2 <- roc_dat_1[rep(seq_len(nrow(roc_dat_1)), each = 4), ]

roc_dat_2$Method <- c("Read","Umi","Fragment","INSPIIRED")

roc_dat_2 %>% ggplot(aes(1 - specificity, sensitivity,color=Method)) + geom_point()+geom_abline(intercept = 0, slope = 1, linetype = "dotted")+xlim(0,1)+ylim(0,1)


