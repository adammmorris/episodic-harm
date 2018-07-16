# setup -------------------------------------------------------------------

# only works in Rstudio -- otherwise you have to set the path manually!
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

require(ggplot2)
require(lme4)
require(lmerTest)
require(dplyr)
require(RColorBrewer)
require(mediation)

theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black"),
             axis.text=element_text(size=20, colour = "black"), axis.title=element_text(size=18, face = "bold"), axis.title.x = element_text(vjust = 0),
             legend.title = element_text(size = 24, face = "bold"), legend.text = element_text(size = 20), plot.title = element_text(size = 26, face = "bold", vjust = 1),
             axis.title.y = element_text(vjust = 1))

se <- function(data) {return(sd(data)/sqrt(length(data)))}
dodge <- position_dodge(width=0.9)
getExcludeSubj = function(data, minRespLength, correctNumResp) {
  return((data %>% group_by(subject) %>%
            summarise(nResponses = length(response),
                      respLength = mean(nchar(encodeString(as.character(na.omit(response)))))) %>%
            filter(respLength < minRespLength | nResponses != correctNumResp))$subject)
}

# experiments 2a, 2b, 3, and 4, respectively
studies = c('rhyme', 'rhyme_justifboth', 'read', 'scrambled')

# import data -------------------------------------------------------------

# choose which study to analyze
# set to studies[c(1,2,3)] to get combined analysis
study = studies[c(2)]

# read in data
# 'imagine' column is trial type (imagine vs. control)
# 'condition' column is justification condition
data <- read.csv("data.csv") %>% arrange(subject) %>%
  filter(version %in% study) %>%
  mutate(imagine = factor(imagine, labels = c("Control", "Imagine"), levels = c(0,1)),
         condition = factor(condition, labels = c("Positive","Neutral","Negative"), levels = c(2,1,0)),
         rating_vivid = (rating_detail + rating_cohere) / 2)

# exclusion
data = data %>%
  filter(!(subject %in% getExcludeSubj(data, 50,
                                       ifelse(study == 'scrambled',
                                              length(unique(story_id)) / 2,
                                              length(unique(story_id))))))

# neutral condition is only used in expt. 2a
if (length(study) > 1 | !(study %in% c('rhyme', 'rhyme_justifboth'))) {
  data = data %>% filter(condition != 'Neutral')
}

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



# main test ---------------------------------------------------------------
# how does imagining (vs control) alter the likelihood of performing behavior in each condition?

## make plot
if (length(study) == 1 & study %in% c('rhyme', 'rhyme_justifboth')) {
  plt.labs = c("Justified", "Neutral", "Unjustified")
  plt.vals = c("Darkgreen", "Blue", "Red")
} else {
  plt.labs = c("Justified", "Unjustified")
  plt.vals = c("Darkgreen", "Red")
}
ggplot(df.bycond, aes(x = imagine, y = will, colour = condition, group = condition)) +
  geom_line(aes(), size = 1) +
  geom_point(aes(), size = 5) +
  geom_errorbar(aes(ymax = will + will.se, ymin = will - will.se), width = .1) +
  labs(x = "", y = "Likelihood of performing behavior") +
  theme(axis.title.y = element_text(vjust = 1)) +
  guides(linetype = guide_legend(title = "Condition"),
         colour = guide_legend(title = "Condition"),
         shape = guide_legend(title = "Condition")) +
  ylim(2, 7) +
  scale_colour_manual(labels = plt.labs, values = plt.vals) +
  theme(panel.background = element_rect(fill = 'white'))
  

## fit linear models

# test for interaction
model <- lmer(rating_will ~ condition * imagine + (1 + imagine | subject) + (1 + condition * imagine | story_id), data = data, REML = F)
model.noint <- lmer(rating_will ~ condition + imagine + (1 + imagine | subject) + (1 + condition * imagine | story_id), data = data, REML = F)
anova(model, model.noint) 

# test for main effects
model.onlyimag <- lmer(rating_will ~ imagine + (1 + imagine | subject) + (1 + condition * imagine | story_id), data = data, REML = F)
model.onlycond <- lmer(rating_will ~ condition + (1 + imagine | subject) + (1 + condition * imagine | story_id), data = data, REML = F)
anova(model.noint, model.onlycond)
anova(model.noint, model.onlyimag)

# test for simple effects
model.pos <- lmer(rating_will ~ imagine + (1 + imagine | subject) + (1 + imagine | story_id), data = data %>% filter(condition == 'Positive'), REML = F)
model.pos.null <- lmer(rating_will ~ 1 + (1 + imagine | subject) + (1 + imagine | story_id), data = data %>% filter(condition == 'Positive'), REML = F)
anova(model.pos, model.pos.null)

model.neut <- lmer(rating_will ~ imagine + (1 + imagine | subject) + (1 + imagine | story_id), data = data %>% filter(condition == 'Neutral'), REML = F)
model.neut.null <- lmer(rating_will ~ 1 + (1 + imagine | subject) + (1 + imagine | story_id), data = data %>% filter(condition == 'Neutral'), REML = F)
anova(model.neut, model.neut.null)

model.neg <- lmer(rating_will ~ imagine + (1 + imagine | subject) + (1 + imagine | story_id), data = data %>% filter(condition == 'Negative'), REML = F)
model.neg.null <- lmer(rating_will ~ 1 + (1 + imagine | subject) + (1 + imagine | story_id), data = data %>% filter(condition == 'Negative'), REML = F)
anova(model.neg, model.neg.null)

# get bayes factor for null
bf = anovaBF(rating_will ~ imagine + subject, data = data %>% filter(condition == 'Negative'), whichRandom="subject")
1 / exp(bf@bayesFactor$bf)

# affect -----------------------------------------------
df.affect = df.bysubj %>% group_by(condition) %>%
  summarize(happy = mean(rating_happy), happy.se = se(rating_happy))

## plot
ggplot(df.bycond, aes(x = imagine, y = happy, colour = condition, group = condition)) +
  geom_line(aes(), size = 1) +
  geom_point(aes(), size = 5) +
  geom_errorbar(aes(ymax = happy + happy.se, ymin = happy - happy.se), width = .1) +
  labs(x = "", y = "How the story made you feel") +
  theme(axis.title.y = element_text(vjust = 1)) +
  guides(linetype = guide_legend(title = "Condition"),
         colour = guide_legend(title = "Condition"),
         shape = guide_legend(title = "Condition")) +
  ylim(2, 7) +
  scale_colour_manual(labels = plt.labs, values = plt.vals)

# test for effect of condition
model.onaffect.cond <- lmer(rating_happy ~ condition + (1 | subject) + (1 + condition | story_id), data = data, REML = F)
model.onaffect.cond.null <- lmer(rating_happy ~ 1 + (1 | subject) + (1 + condition | story_id), data = data, REML = F)
anova(model.onaffect.cond, model.onaffect.cond.null)

# test for interaction on affect
model.onaffect <- lmer(rating_happy ~ condition * imagine + (1 | subject) + (1 + condition * imagine | story_id), data = data, REML = F)
model.onaffect.noint <- lmer(rating_happy ~ condition + imagine + (1 | subject) + (1 + condition * imagine | story_id), data = data, REML = F)
anova(model.onaffect, model.onaffect.noint)

# test for simple effects
model.pos <- lmer(rating_happy ~ imagine + (1 + imagine | subject) + (1 + imagine | story_id), data = data %>% filter(condition == 'Positive'), REML = F)
model.pos.null <- lmer(rating_happy ~ 1 + (1 + imagine | subject) + (1 + imagine | story_id), data = data %>% filter(condition == 'Positive'), REML = F)
anova(model.pos, model.pos.null)

model.neut <- lmer(rating_happy ~ imagine + (1 + imagine | subject) + (1 + imagine | story_id), data = data %>% filter(condition == 'Neutral'), REML = F)
model.neut.null <- lmer(rating_happy ~ 1 + (1 + imagine | subject) + (1 + imagine | story_id), data = data %>% filter(condition == 'Neutral'), REML = F)
anova(model.neut, model.neut.null)

model.neg <- lmer(rating_happy ~ imagine + (1 + imagine | subject) + (1 + imagine | story_id), data = data %>% filter(condition == 'Negative'), REML = F)
model.neg.null <- lmer(rating_happy ~ 1 + (1 + imagine | subject) + (1 + imagine | story_id), data = data %>% filter(condition == 'Negative'), REML = F)
anova(model.neg, model.neg.null)

## how does affect impact willingness?
# plot willingness vs affect
ggplot(data, aes(x = rating_happy, y = rating_will)) +
  geom_jitter() +
  geom_smooth(method = lm) +
  labs(x = "How the story made you feel", y = "Willingness") +
  theme(axis.title.y = element_text(vjust = 1))

# plot willingness vs affect, split by all
ggplot(data, aes(x = rating_happy, y = rating_will)) +
  geom_jitter() +
  geom_smooth(method = lm) +
  labs(x = "How the story made you feel", y = "Willingness") +
  theme(axis.title.y = element_text(vjust = 1)) +
  facet_wrap(~ imagine + condition)

## mediation analysis
# this is the main mediation analysis,
# testing whether affect mediates the effect of imagining justified harm on
# reported likelihood of harming.
detach(package:lmerTest, unload=T)
data.pos = data %>% filter(condition == 'Positive') %>% mutate(imagine.numeric = as.numeric(imagine))
model.med1 = lmer(rating_happy ~ imagine + (1 + imagine.numeric | subject), data = data.pos)
model.med2 = lmer(rating_will ~ imagine + rating_happy + (1 + imagine.numeric + rating_happy | subject), data = data.pos)
model.med3 = mediate(model.med1, model.med2, treat = 'imagine', mediator = 'rating_happy')
summary(model.med3)
plot(model.med3)
require(lmerTest)

# get coefficients
model.med1.coeff = lmer(rating_happy ~ imagine + (1 + imagine.numeric | subject), data = data.pos)
summary(model.med1.coeff)
model.med2.coeff = lmer(rating_will ~ rating_happy + (1 + rating_happy | subject), data = data.pos)
summary(model.med2.coeff)
model.med3.coeff = lmer(rating_will ~ imagine + (1 + imagine.numeric | subject), data = data.pos)
summary(model.med3.coeff)

# this is a secondary mediation analysis specifically for expt. 4,
# testing whether, within the 'matched' condition, affect mediates the effect of justification
# condition on reported likelihood of harming
detach(package:lmerTest, unload=T)
data.pos = data %>% filter(imagine == 'Imagine') %>% mutate(imagine.numeric = as.numeric(condition))
model.med1 = lmer(rating_happy ~ condition + (1 + imagine.numeric | subject), data = data.pos)
model.med2 = lmer(rating_will ~ condition + rating_happy + (1 + imagine.numeric + rating_happy | subject), data = data.pos)
model.med3 = mediate(model.med1, model.med2, treat = 'condition', mediator = 'rating_happy')
summary(model.med3)
plot(model.med3)
require(lmerTest)

# model.affect.simple = lmer(rating_will ~ rating_happy + (rating_happy | subject) + (rating_happy | story_id), data = data, REML = F)
# model.affect.simple.null = lmer(rating_will ~ 1 + (rating_happy | subject) + (rating_happy | story_id), data = data, REML = F)
# summary(model.affect.simple)
# anova(model.affect.simple,model.affect.simple.null)

# vividness ---------------------------------------------------------------
## is vividness similar across conditions?
ggplot(df.bycond, aes(x = condition, y = vivid)) +
  geom_col(aes(), position = dodge) +
  geom_errorbar(aes(ymax = vivid + vivid.se, ymin = vivid - vivid.se), width = .1) +
  labs(x = "", y = "Vividness")

model.onvivid <- lmer(rating_vivid ~ condition + (1 | subject) + (1 + condition | story_id),
                            data = data %>% filter(imagine == 'Imagine'), REML = F)
model.onvivid.null <- lmer(rating_vivid ~ 1 + (1 | subject) + (1 + condition | story_id),
                                 data = data %>% filter(imagine == 'Imagine'), REML = F)
anova(model.onvivid, model.onvivid.null)

## how does vividness relate to willingness?
ggplot(data %>% filter(imagine == 'Imagine'), aes(x = rating_vivid, y = rating_will)) +
  geom_point(alpha = 0, aes(color = condition)) + geom_smooth(method = "lm", se = T, aes(color = condition)) +
  labs(x = "Vividness of simulation", y = "Likelihood of performing behavior") +
  guides(colour = guide_legend(title = 'Condition')) +
  scale_colour_manual(labels = plt.labs, values = plt.vals) +
  ylim(1, 7) +
  theme(panel.background = element_rect(fill = 'white'))

# is there an interaction?
model.vivid <- lmer(rating_will ~ condition * rating_vivid + (1 + rating_vivid | subject) + (1 + condition * rating_vivid | story_id),
                     data = data %>% filter(imagine == 'Imagine'), REML = F)
model.vivid.noint <- lmer(rating_will ~ condition + rating_vivid + (1 + rating_vivid | subject) + (1 + condition * rating_vivid | story_id),
                           data = data %>% filter(imagine == 'Imagine'), REML = F)
anova(model.vivid, model.vivid.noint)

# main effect
model.vivid.null <- lmer(rating_will ~ condition + (1 + rating_vivid | subject) + (1 + condition * rating_vivid | story_id),
                          data = data %>% filter(imagine == 'Imagine'), REML = F)
anova(model.vivid.noint, model.vivid.null)

# simple effect
model.vivid.pos <- lmer(rating_will ~ rating_vivid + (1 + rating_vivid | subject) + (1 + rating_vivid | story_id),
                          data = data %>% filter(imagine == 'Imagine', condition == 'Positive'), REML = F)
model.vivid.pos.null <- lmer(rating_will ~  1 + (1 + rating_vivid | subject) + (1 + rating_vivid | story_id),
                         data = data %>% filter(imagine == 'Imagine', condition == 'Positive'), REML = F)
anova(model.vivid.pos, model.vivid.pos.null)

model.vivid.neg <- lmer(rating_will ~ rating_vivid + (1 + rating_vivid | subject) + (1 + rating_vivid | story_id),
                        data = data %>% filter(imagine == 'Imagine', condition == 'Negative'), REML = F)
model.vivid.neg.null <- lmer(rating_will ~  1 + (1 + rating_vivid | subject) + (1 + rating_vivid | story_id),
                             data = data %>% filter(imagine == 'Imagine', condition == 'Negative'), REML = F)
anova(model.vivid.neg, model.vivid.neg.null)
