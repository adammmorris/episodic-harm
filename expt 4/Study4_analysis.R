# setup -------------------------------------------------------------------

# load packages with groundhog
require(groundhog)
pkgs = c('ggplot2', 'lme4', 'lmerTest',  'dplyr', 'afex', 'RColorBrewer', 'mediation', 'simr')
groundhog.library(pkgs, '2021-04-01')

# Only works in RStudio -- otherwise you have to set the path manually
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

se <- function(data) {return(sd(data)/sqrt(length(data)))}
se.prop <- function(data) {return(sqrt(mean(data)*(1-mean(data))/length(data)))}
dodge <- position_dodge(width=0.9)

theme_update(strip.background = element_blank(),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             panel.background = element_blank(),
             plot.background = element_blank(),
             axis.text=element_text(size=30, colour = "black"),
             axis.title=element_text(size=18, face = "bold"),
             axis.title.x = element_text(vjust = 0),
             legend.title = element_text(size = 24, face = "bold"),
             legend.text = element_text(size = 20),
             plot.title = element_text(size = 26, face = "bold", vjust = 1),
             panel.margin = unit(1.0, "lines"), 
             plot.margin = unit(c(0.5,  0.5, 0.5, 0.5), "lines"),
             axis.line = element_line(colour = "black", size = 2),
             axis.ticks = element_line(color = 'black', size = 3),
             axis.ticks.length = unit(.25, 'cm')
)


print.stats = function(m, m.std, param) {
  m.coef = summary(m)$coefficients[param,]
  m.std.coef = summary(m.std)$coefficients[param,]
  b = m.coef[1]
  t = m.coef[4]
  df = m.coef[3]
  p = m.coef[5]
  m.ci = suppressMessages(confint(m, param))
  m.std.ci = suppressMessages(confint(m.std, param))
  beta = m.std.coef[1]
  
  if (p >= .001) {
    p.string = paste0(', p = ', signif(p,2))
  } else {
    p.string = paste0(', p < .001')
  }
  
  paste0('b = ', signif(b,2), ' [', signif(m.ci[1], 2), ',', signif(m.ci[2], 2),
         '], β = ', signif(beta, 2), ' [', signif(m.std.ci[1], 2), ',',
         signif(m.std.ci[2], 2), '], t(', signif(df, 2), ') = ',
         signif(t, 2), p.string)
}

options(digits = 2)

getExcludeSubj = function(data, minRespLength, correctNumResp) {
  return((data %>% group_by(subject) %>%
            summarise(nResponses = length(response),
                      respLength = mean(nchar(encodeString(as.character(na.omit(response)))))) %>%
            filter(respLength < minRespLength | nResponses != correctNumResp))$subject)
}

# import data -------------------------------------------------------------
expt = '5'

data <- read.csv(paste0('Study', expt, '_data.csv')) %>% arrange(subject) %>%
  mutate(imagine = factor(imagine, labels = c("Control", "Imagine"), levels = c(0,1)),
         condition = factor(condition, labels = c("Positive","Negative"), levels = c(2,0)),
         rating_vivid = (rating_detail + rating_cohere) / 2,
         rating_check = rating_check == 1)

# exclusion
data = data %>%
  filter(!(subject %in% getExcludeSubj(data, 50, length(unique(story_id)) / 2)))

# make data frames grouped by subject & condition
df.bysubj <- data %>% group_by(condition, imagine, subject) %>%
  summarise(rating_will = mean(rating_will), rating_happy = mean(rating_happy),
            rating_check = mean(rating_check), rating_justif = mean(rating_justif),
            rating_vivid = mean(rating_vivid))
df.bycond <- df.bysubj %>% group_by(condition, imagine) %>%
  summarise(will = mean(rating_will), will.se = se(rating_will),
            happy = mean(rating_happy), happy.se = se(rating_happy),
            check = mean(rating_check), check.se = se(rating_check),
            justif = mean(rating_justif), justif.se = se(rating_justif),
            vivid = mean(rating_vivid), vivid.se = se(rating_vivid))
# manipulation check ------------------------------------------------------

# plotting variables
plt.labs = c("Justified", "Unjustified")
plt.vals = c("Darkgreen", "Red")

# Check out justification
df.bycond %>% dplyr::select(imagine, check, justif, justif.se) %>%
  mutate(justif.ci.low = justif - 1.96*justif.se,
         justif.ci.high = justif + 1.96*justif.se)

# Does justification rating differ across conditions?
model.justif = lmer_alt(rating_justif ~ condition + (1 | subject) +
                          (condition | story_id),
                        data = data %>% filter(imagine == 'Imagine'))
summary(rePCA(model.justif))
model.justif.r1 = lmer_alt(rating_justif ~ condition + (1 | subject) +
                          (condition || story_id),
                        data = data %>% filter(imagine == 'Imagine'))
summary(rePCA(model.justif.r1))
model.justif.r2 = lmer_alt(rating_justif ~ condition + (1 | subject) +
                             (1 | story_id),
                           data = data %>% filter(imagine == 'Imagine'))
summary(rePCA(model.justif.r2))

summary(model.justif.r2)
model.justif.r2.std = lmer_alt(scale(rating_justif) ~ condition + (1 | subject) +
                             (1 | story_id),
                           data = data %>% filter(imagine == 'Imagine'))
print.stats(model.justif.r2, model.justif.r2.std, 'conditionNegative')


# Does people's self-reported ability to follow the instructions differ across conditions?
model.check = lmer_alt(rating_check ~ condition + (1 | subject) +
                         (condition | story_id),
                       data = data %>% filter(imagine == 'Imagine'))
summary(rePCA(model.check))
model.check.r1 = lmer_alt(rating_check ~ condition + (1 | subject) +
                         (condition || story_id),
                       data = data %>% filter(imagine == 'Imagine'))
summary(rePCA(model.check.r1))
model.check.r2 = lmer_alt(rating_check ~ condition + (1 | subject) +
                            (1 | story_id),
                          data = data %>% filter(imagine == 'Imagine'))
summary(rePCA(model.check.r2))

summary(model.check.r2)
model.check.r2.std = lmer_alt(scale(rating_check) ~ condition + (condition | subject) +
                            (1 | story_id),
                          data = data %>% filter(imagine == 'Imagine'))
print.stats(model.check.r2, model.check.r2.std, 'conditionNegative')

# main test ---------------------------------------------------------------
# how does imagining (vs control) alter the likelihood of performing behavior in each condition?

## make plot
ggplot(df.bycond, aes(x = imagine, y = will, colour = condition, group = condition)) +
  geom_line(aes(), size = 1) +
  geom_point(aes(), size = 5) +
  geom_errorbar(aes(ymax = will + will.se, ymin = will - will.se), width = .1) +
  labs(x = "", y = "Reported likelihood of\nperforming behavior") +
  theme(axis.title.y = element_text(vjust = 1)) +
  guides(linetype = "none",
         colour = "none",
         shape = "none") +
  ylim(2, 6) +
  scale_colour_manual(labels = plt.labs, values = plt.vals) +
  scale_x_discrete(labels = c('Mismatched', 'Matched'))

## fit linear models

## test for interaction & main effects
model <- lmer_alt(rating_will ~ condition * imagine + (1 | subject) +
                    (condition * imagine | story_id), data = data)
summary(rePCA(model))
model.r1 <- lmer_alt(rating_will ~ condition * imagine + (1 | subject) +
                       (condition * imagine || story_id), data = data)
summary(rePCA(model.r1))
model.r2 <- lmer_alt(rating_will ~ condition * imagine + (1 | subject) +
                       (condition + imagine || story_id), data = data)
summary(rePCA(model.r2))
model.r3 <- lmer_alt(rating_will ~ condition * imagine + (1 | subject) +
                       (condition || story_id), data = data)
summary(rePCA(model.r3))
model.r4 <- lmer_alt(rating_will ~ condition * imagine + (1 | subject) +
                       (imagine || story_id), data = data)
summary(rePCA(model.r4))
model.r5 <- lmer_alt(rating_will ~ condition * imagine + (1 | subject) +
                       (condition | story_id), data = data)
summary(rePCA(model.r5))

summary(model.r5)
model.r5.std <- lmer_alt(scale(rating_will) ~ condition * imagine + (1 | subject) +
                           (condition | story_id), data = data)
print.stats(model.r5, model.r5.std, 'conditionNegative:imagineImagine')
print.stats(model.r5, model.r5.std, 'conditionNegative')
print.stats(model.r5, model.r5.std, 'imagineImagine')
powerSim(model.r5, test = fixed('conditionNegative:imagineImagine', 't'), nsim = 100)
# power 31% [22, 41]

## test for simple effects
# justified condition
model.pos <- lmer_alt(rating_will ~ imagine + (1 | subject) +
                        (imagine | story_id),
                      data = data %>% filter(condition == 'Positive'))
summary(rePCA(model.pos))

summary(model.pos)
model.pos.std <- lmer_alt(scale(rating_will) ~ imagine + (1 | subject) +
                            (imagine | story_id),
                          data = data %>% filter(condition == 'Positive'))
print.stats(model.pos, model.pos.std, 'imagineImagine')
powerSim(model.pos, test = fixed('imagineImagine', 't'), nsim = 100)
# power 57% [47, 67]

# unjustified
model.neg <- lmer_alt(rating_will ~ imagine + (1 | subject) +
                        (imagine | story_id),
                      data = data %>% filter(condition == 'Negative'))
summary(rePCA(model.neg))
model.neg.r1 <- lmer_alt(rating_will ~ imagine + (1 | subject) +
                        (imagine || story_id),
                      data = data %>% filter(condition == 'Negative'))
summary(rePCA(model.neg.r1))

summary(model.neg.r1)
model.neg.r1.std <- lmer_alt(scale(rating_will) ~ imagine + (1 | subject) +
                              (imagine || story_id), data = data %>% filter(condition == 'Negative'))
print.stats(model.neg.r1, model.neg.r1.std, 'imagineImagine')

# vividness ---------------------------------------------------------------
## is vividness similar across conditions?
ggplot(df.bycond %>% filter(imagine == 'Imagine'), aes(x = condition, y = vivid)) +
  geom_col(aes(), position = dodge) +
  geom_errorbar(aes(ymax = vivid + vivid.se, ymin = vivid - vivid.se), width = .1) +
  labs(x = "", y = "Vividness")

## how does vividness relate to willingness?
ggplot(data %>% filter(imagine == 'Imagine'), aes(x = rating_vivid, y = rating_will)) +
  geom_point(alpha = 0, aes(color = condition)) + geom_smooth(method = "lm", se = T, aes(color = condition)) +
  labs(x = "Vividness of simulation", y = "Likelihood of performing behavior") +
  guides(colour = guide_legend(title = 'Condition')) +
  scale_colour_manual(labels = plt.labs, values = plt.vals) +
  ylim(1, 7) +
  theme(panel.background = element_rect(fill = 'white'))

# need to center the vividness ratings
data.vivid = data %>% filter(imagine == 'Imagine') %>%
  mutate(rating_vivid = scale(rating_vivid))
model.vivid <- lmer_alt(rating_will ~ condition * rating_vivid +
                          (rating_vivid | subject) +
                          (condition * rating_vivid | story_id),
                        data = data.vivid)
summary(rePCA(model.vivid))
model.vivid.r1 <- lmer_alt(rating_will ~ condition * rating_vivid +
                             (rating_vivid || subject) +
                             (condition * rating_vivid || story_id),
                           data = data.vivid)
summary(rePCA(model.vivid.r1))
model.vivid.r2 <- lmer_alt(rating_will ~ condition * rating_vivid +
                             (1 | subject) +
                             (condition + rating_vivid || story_id),
                           data = data.vivid)
summary(rePCA(model.vivid.r2))
model.vivid.r3 <- lmer_alt(rating_will ~ condition * rating_vivid +
                             (1 | subject) +
                             (condition + rating_vivid | story_id),
                           data = data.vivid)
summary(rePCA(model.vivid.r3))

summary(model.vivid.r3)
model.vivid.r3.std <- lmer_alt(scale(rating_will) ~ condition * rating_vivid +
                                 (1 | subject) +
                                 (condition + rating_vivid | story_id),
                               data = data.vivid)
print.stats(model.vivid.r3, model.vivid.r3.std, 'conditionNegative:rating_vivid')
print.stats(model.vivid.r3, model.vivid.r3.std, 'rating_vivid')

# affect -----------------------------------------------
df.affect = df.bysubj %>% group_by(condition) %>%
  summarize(happy = mean(rating_happy), happy.se = se(rating_happy))

## plot
ggplot(df.bycond, aes(x = imagine, y = happy, colour = condition, group = condition)) +
  geom_line(aes(), size = 1) +
  geom_point(aes(), size = 5) +
  geom_errorbar(aes(ymax = happy + happy.se, ymin = happy - happy.se), width = .1) +
  labs(x = "", y = "How the story\nmade you feel") +
  theme(axis.title.y = element_text(vjust = 1)) +
  guides(linetype = guide_legend(title = "Condition"),
         colour = guide_legend(title = "Condition"),
         shape = guide_legend(title = "Condition")) +
  ylim(2, 5) +
  scale_colour_manual(labels = plt.labs, values = plt.vals) +
  scale_x_discrete(labels = c('Reading\nControl', 'Imagine')) +
  theme(panel.background = element_rect(fill = 'white'))

# test for effect of condition
model.affect <- lmer_alt(rating_happy ~ condition * imagine +
                           (1 | subject) +
                           (condition * imagine | story_id),
                         data = data)
summary(rePCA(model.affect))
model.affect.r1 <- lmer_alt(rating_happy ~ condition * imagine +
                              (1 | subject) +
                              (condition * imagine || story_id),
                            data = data)
summary(rePCA(model.affect.r1))
model.affect.r2 <- lmer_alt(rating_happy ~ condition * imagine +
                              (1 | subject) +
                              (condition + imagine || story_id),
                            data = data)
summary(rePCA(model.affect.r2))
model.affect.r3 <- lmer_alt(rating_happy ~ condition * imagine +
                              (1 | subject) +
                              (condition || story_id),
                            data = data)
summary(rePCA(model.affect.r3))
model.affect.r4 <- lmer_alt(rating_happy ~ condition * imagine +
                              (1 | subject) +
                              (imagine || story_id),
                            data = data)
summary(rePCA(model.affect.r4))
model.affect.r5 <- lmer_alt(rating_happy ~ condition * imagine +
                              (1 | subject) +
                              (condition | story_id),
                            data = data)
summary(rePCA(model.affect.r5))

summary(model.affect.r5)
model.affect.r5.std <- lmer_alt(scale(rating_happy) ~ condition * imagine +
                                  (1 | subject) +
                                  (condition | story_id),
                                data = data)
print.stats(model.affect.r5, model.affect.r5.std, 'conditionNegative')
print.stats(model.affect.r5, model.affect.r5.std, 'conditionNegative:imagineImagine')

## mediation analysis
detach(package:afex, unload=T)
detach(package:lmerTest, unload = T)
data.pos = data %>% filter(condition == 'Positive') %>%
  mutate(imagine.numeric = as.numeric(imagine))

model.med1 = lmer(rating_happy ~ imagine + (1 | subject), data = data.pos)
summary(rePCA(model.med1))

model.med2 = lmer(rating_will ~ imagine + rating_happy +
                    (rating_happy | subject), data = data.pos)
summary(rePCA(model.med2))
model.med2.r1 = lmer(rating_will ~ imagine + rating_happy +
                       (1 | subject), data = data.pos)
summary(rePCA(model.med2.r1))

model.med3 = mediate(model.med1, model.med2.r1, treat = 'imagine', mediator = 'rating_happy')
summary(model.med3)
plot(model.med3)

# this is a secondary mediation analysis specifically for expt. 5,
# testing whether, within the 'matched' condition, affect mediates the effect of justification
# condition on reported likelihood of harming
data.matched = data %>% filter(imagine == 'Imagine') %>%
  mutate(condition.numeric = as.numeric(condition))
model2.med1 = lmer(rating_happy ~ condition + (1 | subject), data = data.matched)
model2.med2 = lmer(rating_will ~ condition + rating_happy + (rating_happy | subject),
                  data = data.matched)
summary(rePCA(model2.med2))
model2.med2.r1 = lmer(rating_will ~ condition + rating_happy + (1 | subject),
                  data = data.matched)
summary(rePCA(model2.med2.r1))
model2.med3 = mediate(model2.med1, model2.med2.r1, treat = 'condition', mediator = 'rating_happy')
summary(model2.med3)
plot(model2.med3)

groundhog.library(c('lmerTest', 'afex'), '2021-04-01')


# save --------------------------------------------------------------------

save.image(paste0('Study', expt, '_output.rdata'))

