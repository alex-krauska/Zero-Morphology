## Zero Morphology in On-line Sentence Processing: Exploratory Data Analysis
# Alex Krauska

# load packages
library(tidyverse)
library(modelr)
library(lme4)
library(lmerTest)
library(MASS)
library(languageR)

# load unprocessed eye tracking data set
zm_data <- read_tsv("data/unprocessed/ET_data.txt")

# load processed acceptability rating data set
zm_rating <- readRDS("data/zm_rating.RDS")

# Function provided by Dr. Shayne Sloggett
zscoreET <- function(data, cutoff=5, reject = T) {
  data$logRT = log(data$value)
  means <- with(data, tapply(logRT, list(region, fixationtype, subj), mean, na.rm = 1))
  sds <- with(data, tapply(logRT, list(region, fixationtype, subj), sd, na.rm = 1))
  data$z = 0
  for(row in 1:nrow(data)){
    measure = data$fixationtype[row]
    region = data$region[row]
    subj = as.character(data$subj[row])
    z = abs(data$logRT[row] - means[region, measure, subj])/sds[region, measure, subj]
    if(measure != 'pr' & !is.na(z)){ data$z[row] = z }
  }
  if(reject){ return(subset(data, z < cutoff)) }
  return(data)
}

# Experiment structure: derivation status x base category
exp_str <- data.frame(
  cond = c(21:24),
  derivation = c('derived', 'underived', 'derived', 'underived'),
  base = c('noun', 'noun', 'verb', 'verb'))

# Tidy data: filter for fixation type and NA, select relevant info, add logRT, 
# apply zscoreET, join with acceptability judgement ratings and experiment 
# structure, save as RDS
zm_data <- zm_data %>% 
  filter(is.na(value)==FALSE, fixationtype != "pr", fixationtype != "rr") %>% 
  dplyr::select("subj", "cond", "item", "value", "region", "fixationtype") %>% 
  mutate(logRT = log(value),
         item = item - 72) %>% 
  zscoreET() %>%
  left_join(zm_rating, by = c("item", "cond" = "condition")) %>% 
  left_join(exp_str, by = "cond")

# save zm_data as an RDS for faster access in the future
zm_data %>% 
  saveRDS("data/zm_data.RDS")

# region and fixation labels for graphs
region_labs <- c(
  `1` = "... to / the",
  `2` = "test / visit",
  `3` = "for the...",
  `ff` = "first fixation",
  `fp` = "first pass",
  `rp` = "regression path",
  `tt` = "total time"
)

# Visualization 1: All regions and fixation types
zm_data %>% 
  ggplot(aes(x = derivation, y = logRT, group = cond, fill = base))+
  geom_boxplot()+
  facet_grid(fixationtype~region, scales = "free", labeller = as_labeller(region_labs))+
  scale_x_discrete("Derivation", limits = c("underived", "derived"))


# Modeling ----------------------------------------------------------------

# set contrasts
contrasts(zm_data$base) = c(-0.5,0.5)
contrasts(zm_data$derivation) = c(0.5,-0.5)

### Is there an effect of derivation in region 1 for any dependent measure?
# Subsetting region 1
region_1 <- zm_data %>% 
  filter(region == 1)

reg_1_model <- function(df){
  lmer(logRT ~ derivation*base + (1|subj) + (1|item), data = df)
}

reg_1_mod <- region_1 %>% 
  group_by(fixationtype) %>% 
  nest() %>% 
  mutate(model = map(data, reg_1_model))

summary(reg_1_mod$model[[1]]) # ff
summary(reg_1_mod$model[[2]]) # fp
summary(reg_1_mod$model[[3]]) # rp
summary(reg_1_mod$model[[4]]) # tt

### Does mean rating effect reading times at region 2?
# Subsetting region 2
region_2 <- zm_data %>% 
  filter(region == 2)

rating_fn <- function(df){
  lmer(logRT ~ mean_rating + (1|subj) + (1|item), data = df)
}

rating_model <- region_2 %>%
  group_by(fixationtype) %>% 
  nest() %>% 
  mutate(model = map(data, rating_fn))

summary(rating_model$model[[1]]) # ff
summary(rating_model$model[[2]]) # fp
summary(rating_model$model[[3]]) # rp
summary(rating_model$model[[4]]) # tt

# resoundingly, yes, mean rating does affect the reading times.

region_2 %>% 
  ggplot(aes(x = mean_rating, y = logRT))+
  geom_point()+
  geom_smooth(method = lm)+
  facet_wrap(~fixationtype, labeller = as_labeller(region_labs))+
  ggtitle("Log Reading Time at Region 2\nPredicted by the Mean Rating\nfor Four Fixation Types")+
  xlab("Mean Rating")+
  ylab("Log Reading Time")


### Does derivation incur a reading time slowdown?
# model 1: logRT predicted by the interaction between derivation and
# base, with subject and item as random effects
der_model_fn1 <- function(df){
  lmer(logRT ~ derivation*base + (1|subj) + (1|item), data = df)
}

# model 2: logRT predicted by derivation nested within levels of base,
# with subject and item as random effects
der_model_fn2 <- function(df){
  lmer(logRT ~ base/derivation + (1|subj) + (1|item) + 0, data = df)
}

# model 3: logRT predicted by derivation nested within levels of base,
# with mean_rating as a fixed effect and subject as a random effect
der_model_fn3 <- function(df){
  lmer(logRT ~ base/derivation + mean_rating + (1|subj) + 0, data = df)
}

der_models <- region_2 %>%
  group_by(fixationtype) %>% 
  nest() %>% 
  mutate(
    model_1 = map(data, der_model_fn1),
    pred_1 = map2(data, model_1, add_predictions),
    model_2 = map(data, der_model_fn2),
    pred_2 = map2(data, model_2, add_predictions),
    model_3 = map(data, der_model_fn3),
    pred_3 = map2(data, model_3, add_predictions)
  )

# interaction between derivation and base
summary(der_models$model_1[[1]]) # ff
summary(der_models$model_1[[2]]) # fp
summary(der_models$model_1[[3]]) # rp
summary(der_models$model_1[[4]]) # tt

# derivation nested within levels of base
summary(der_models$model_2[[1]]) # ff
summary(der_models$model_2[[2]]) # fp
summary(der_models$model_2[[3]]) # rp
summary(der_models$model_2[[4]]) # tt

# models including mean_rating as a fixed effect
summary(der_models$model_3[[1]]) # ff
summary(der_models$model_3[[2]]) # fp
summary(der_models$model_3[[3]]) # rp
summary(der_models$model_3[[4]]) # tt


# Visualizations ----------------------------------------------------------

der_preds <- der_models  %>% 
  unnest(pred_1)

# model 1 derivation predictions
der_preds %>%
  ggplot(aes(x = derivation, y = pred, group = cond, fill = base))+
  geom_boxplot()+
  scale_x_discrete("Derivation", limits = c("underived", "derived"))+
  facet_wrap(~fixationtype, scales = "free", labeller = as_labeller(region_labs))+
  ylab("Log Reading Time")

# derivation predictions plus linear model
der_preds %>% 
  filter(fixationtype %in% c("fp", "rp")) %>% 
  ggplot(aes(x = derivation, y = pred, group = base, color = base))+
    geom_smooth(method = "lm")+
    facet_grid(~fixationtype, scales = "free", labeller = as_labeller(region_labs))+
    geom_jitter(alpha = 0.2)+
    scale_x_discrete("Derivation", limits = c("underived", "derived"))+
    labs(title = "First Pass and Regression Path \nLog Reading Times at the Critical Region \nControlling for Subject and Item")+
    ylab("Log Reading Time")+
    theme(legend.position = "right")

# model 2 derivation predictions
der_preds2 <- der_models  %>% 
  unnest(pred_2)

der_preds2 %>%
  ggplot(aes(x = derivation, y = pred, group = cond, fill = base))+
  geom_boxplot()+
  scale_x_discrete("Derivation", limits = c("underived", "derived"))+
  facet_wrap(~fixationtype, scales = "free", labeller = as_labeller(region_labs))+
  ylab("Log Reading Time")


