# setup -------------------------------------------------------------------

# load packages with groundhog
require(groundhog)
pkgs = c('ggplot2', 'lme4', 'lmerTest',  'dplyr', 'afex', 'RColorBrewer', 'mediation', 'simr')
#groundhog.library(pkgs, '2021-04-01')
lapply(pkgs, require, character.only = T)

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

# read in data
# 'imagine' column is trial type (imagine vs. control)
# 'condition' column is justification condition
data <- read.csv(paste0('Study', expt, '_data.csv')) %>% arrange(subject) %>%
  mutate(imagine = factor(imagine, labels = c("Control", "Imagine"), levels = c(0,1)),
         condition = factor(condition, labels = c("Positive","Neutral","Negative"), levels = c(2,1,0)),
         rating_vivid = (rating_detail + rating_cohere) / 2,
         rating_justif = rating_justif2)

# exclusion
data = data %>%
  filter(!(subject %in% getExcludeSubj(data, 50, length(unique(story_id)))))

# make data frames grouped by subject & condition
df.bysubj <- data %>% group_by(condition, imagine, subject) %>%
  summarise(rating_will = mean(rating_will), rating_happy = mean(rating_happy),
            rating_justif = mean(rating_justif),
            rating_vivid = mean(rating_vivid))
df.bycond <- df.bysubj %>% group_by(condition, imagine) %>%
  summarise(will = mean(rating_will), will.se = se(rating_will),
            happy = mean(rating_happy), happy.se = se(rating_happy),
            #check = mean(rating_check), check.se = se(rating_check),
            justif = mean(rating_justif), justif.se = se(rating_justif),
            vivid = mean(rating_vivid), vivid.se = se(rating_vivid))


# manipulation check ------------------------------------------------------

# plotting variables
plt.labs = c("Justified", "Unjustified")
plt.vals = c("Darkgreen", "Red")

# Check out justification
df.bycond %>% dplyr::select(imagine, justif, justif.se) %>%
  mutate(justif.ci.low = justif - 1.96*justif.se,
         justif.ci.high = justif + 1.96*justif.se)

# Does justification rating differ across conditions?
model.justif = lmer_alt(rating_justif ~ condition + (1 | subject) +
                          (condition | story_id),
                        data = data %>% filter(imagine == 'Imagine'))
summary(rePCA(model.justif))
  

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
  scale_colour_manual(labels = plt.labs, values = plt.vals)

## test for interaction & main effects
model <- lmer_alt(rating_will ~ condition * imagine +
                    (imagine | subject) +
                    (condition * imagine | story_id), data = data)
summary(rePCA(model))
model.r1 <- lmer_alt(rating_will ~ condition * imagine +
                    (imagine | subject) +
                    (condition * imagine || story_id), data = data)
summary(rePCA(model.r1))
summary(model.r1)
model.r1.std <- lmer_alt(scale(rating_will) ~ condition * imagine +
                       (imagine | subject) +
                       (condition * imagine || story_id), data = data)
print.stats(model.r1, model.r1.std, 'conditionNegative:imagineImagine')
print.stats(model.r1, model.r1.std, 'conditionNegative')
print.stats(model.r1, model.r1.std, 'imagineImagine')

df.test = data %>% filter(condition %in% c('Positive', 'Negative'))
model.test <- lmer_alt(scale(rating_will) ~ condition * imagine + (imagine | subject) + (1 | story_id),
                  data = df.test)
summary(model.test)
fixef(model.test)['conditionNegative:imagineImagine'] = -.2

model.test.ext = extend(model.test, along = 'subject', )

model.test2 = lmer(scale(rating_will) ~ condition * imagine + (imagine | subject) + (1 | story_id),
                       data = df.test)

pwr.test = mixedpower(model.test2, df.test, c('condition', 'imagine'),
           simvar = 'subject', steps = c(100,200,500), critical_value = 2, n_sim = 100,
           SESOI = c(.18, -.5, .33, -.2), databased = F)
pwr.test

## test for simple effects
# in justified condition
model.pos <- lmer_alt(rating_will ~ imagine + (imagine | subject) +
                        (imagine | story_id),
                      data = data %>% filter(condition == 'Positive'))
model.pos.std <- lmer_alt(scale(rating_will) ~ imagine + (imagine | subject) +
                        (imagine | story_id),
                      data = data %>% filter(condition == 'Positive'))
summary(rePCA(model.pos))
summary(model.pos)
print.stats(model.pos, model.pos.std, 'imagineImagine')
pwr.pos = powerSim(model.pos.r2, test = fixed('imagineImagine', 't'), nsim = 100)
pwr.pos
# 86% power [77, 92]

model.pos.r2.std.pwr = model.pos.r2.std
fixef(model.pos.r2.std.pwr)['imagineImagine'] = .2
summary(model.pos.r2.std.pwr)
pwr.pos2 = powerSim(model.pos.r2.std.pwr, test = fixed('imagineImagine', 't'), nsim = 500)
pwr.pos2

model.pos.r2.std.pwr.ext = extend(model.pos.r2.std.pwr, along = 'subject', n = 28)
pwr.pos3 = powerSim(model.pos.r2.std.pwr.ext, test = fixed('imagineImagine', 't'), nsim = 500)
pwr.pos3

pwr.pos2.curve = powerCurve(model.pos.r2.std.pwr.ext, along = 'subject',
                            test = fixed('imagineImagine', 't'), breaks = c(100,500), nsim = 100)
pwr.pos2.curve2 = powerCurve(model.pos.r2.std.pwr, within = 'subject', test = fixed('imagineImagine', 't'),
                             breaks = c(100,500))
plot(pwr.pos2.curve2)

# unjustified condition
model.neg <- lmer_alt(rating_will ~ imagine + (imagine | subject) + (imagine | story_id), data = data %>% filter(condition == 'Negative'))
summary(rePCA(model.neg))
model.neg.std <- lmer_alt(scale(rating_will) ~ imagine + (imagine | subject) + (imagine | story_id), data = data %>% filter(condition == 'Negative'))
print.stats(model.neg, model.neg.std, 'imagineImagine')

bf = lmBF(rating_will ~ imagine + subject,
          data = data %>%
            filter(condition == 'Negative') %>%
            mutate(subject = factor(as.numeric(factor(subject)))),
          whichRandom="subject")
bf.null = lmBF(rating_will ~ subject,
               data = data %>%
                 filter(condition == 'Negative') %>%
                 mutate(subject = factor(as.numeric(factor(subject)))),
               whichRandom="subject")
bf.null / bf

# analyzing the role of vividness ---------------------------------------------------------------
## is vividness similar across conditions?
ggplot(df.bycond %>% filter(imagine == 'Imagine'), aes(x = condition, y = vivid)) +
  geom_col(aes(), position = dodge) +
  geom_errorbar(aes(ymax = vivid + vivid.se, ymin = vivid - vivid.se), width = .1) +
  labs(x = "", y = "Vividness")

## how does vividness relate to willingness?
ggplot(data %>% filter(imagine == 'Imagine'), aes(x = rating_vivid, y = rating_will)) +
  geom_point(alpha = 0, aes(color = condition)) + geom_smooth(method = "lm", se = T, aes(color = condition)) +
  labs(x = "Vividness of simulation", y = "Reported likelihood of\nperforming behavior") +
  guides(colour = "none") +
  scale_colour_manual(labels = plt.labs, values = plt.vals) +
  ylim(1, 7)

# need to center the vividness ratings
data.vivid = data %>% filter(imagine == 'Imagine') %>%
  mutate(rating_vivid = scale(rating_vivid))
model.vivid <- lmer_alt(rating_will ~ condition * rating_vivid +
                      (rating_vivid | subject) +
                      (condition * rating_vivid | story_id),
                    data = data.vivid)
summary(rePCA(model.vivid))

model.vivid.r1 <- lmer_alt(rating_will ~ condition * rating_vivid +
                          (rating_vivid | subject) +
                          (condition * rating_vivid || story_id),
                        data = data.vivid)
summary(rePCA(model.vivid.r1))

model.vivid.r2 <- lmer_alt(rating_will ~ condition * rating_vivid +
                             (rating_vivid | subject) +
                             (condition + rating_vivid || story_id),
                           data = data.vivid)
summary(rePCA(model.vivid.r2))
summary(model.vivid.r2)

model.vivid.r2.std <- lmer_alt(rating_will ~ condition * rating_vivid +
                             (rating_vivid | subject) +
                             (condition + rating_vivid || story_id),
                           data = data.vivid %>% mutate(rating_will = scale(rating_will),
                                                        rating_vivid = scale(rating_vivid)))
print.stats(model.vivid.r2, model.vivid.r2.std, 'conditionNegative:rating_vivid')
print.stats(model.vivid.r2, model.vivid.r2.std, 'rating_vivid')

# simple effects
model.vivid.pos <- lmer_alt(rating_will ~ rating_vivid +
                          (rating_vivid | subject) +
                          (rating_vivid | story_id),
                        data = data.vivid %>% filter(condition == 'Positive'))
summary(rePCA(model.vivid.pos))
model.vivid.pos.r1 <- lmer_alt(rating_will ~ rating_vivid +
                              (rating_vivid | subject) +
                              (rating_vivid || story_id),
                            data = data.vivid %>% filter(condition == 'Positive'))
summary(rePCA(model.vivid.pos.r1))
model.vivid.pos.r2 <- lmer_alt(rating_will ~ rating_vivid +
                                 (rating_vivid | subject) +
                                 (1 | story_id),
                               data = data.vivid %>% filter(condition == 'Positive'))
summary(rePCA(model.vivid.pos.r2))
summary(model.vivid.pos.r2)
model.vivid.pos.r2.std <- lmer_alt(rating_will ~ rating_vivid +
                                 (rating_vivid | subject) +
                                 (1 | story_id),
                               data = data.vivid %>% filter(condition == 'Positive') %>%
                                 mutate(rating_vivid = scale(rating_vivid)))
print.stats(model.vivid.pos.r2, model.vivid.pos.r2.std, 'rating_vivid')

model.vivid.neg <- lmer_alt(rating_will ~ rating_vivid +
                              (rating_vivid | subject) +
                              (rating_vivid | story_id),
                            data = data.vivid %>% filter(condition == 'Negative'))
summary(rePCA(model.vivid.neg))
model.vivid.neg.r1 <- lmer_alt(rating_will ~ rating_vivid +
                                 (rating_vivid || subject) +
                                 (rating_vivid || story_id),
                               data = data.vivid %>% filter(condition == 'Negative'))
summary(rePCA(model.vivid.neg.r1))
model.vivid.neg.r2 <- lmer_alt(rating_will ~ rating_vivid +
                                 (1 | subject) +
                                 (1 | story_id),
                               data = data.vivid %>% filter(condition == 'Negative'))
summary(rePCA(model.vivid.neg.r2))
model.vivid.neg.r2.std <- lmer_alt(rating_will ~ rating_vivid +
                                     (1 | subject) +
                                     (1 | story_id),
                                   data = data.vivid %>% filter(condition == 'Negative') %>%
                                     mutate(rating_vivid = scale(rating_vivid)))
print.stats(model.vivid.neg.r2, model.vivid.neg.r2.std, 'rating_vivid')

bf.vivid = lmBF(rating_will ~ rating_vivid + subject,
          data = data.vivid %>%
            filter(condition == 'Negative') %>%
            mutate(subject = factor(as.numeric(factor(subject)))),
          whichRandom="subject")
bf.vivid.null = lmBF(rating_will ~ subject,
               data = data.vivid %>%
                 filter(condition == 'Negative') %>%
                 mutate(subject = factor(as.numeric(factor(subject)))),
               whichRandom="subject")
bf.vivid.null / bf.vivid

# analyze the role of affect -----------------------------------------------
df.affect = df.bysubj %>% group_by(condition) %>%
  summarize(happy = mean(rating_happy), happy.se = se(rating_happy)) %>%
  mutate(happy.ci.low = happy - 1.96*happy.se,
         happy.ci.high = happy + 1.96*happy.se)
df.affect

## plot
ggplot(df.bycond, aes(x = imagine, y = happy, colour = condition, group = condition)) +
  geom_line(aes(), size = 1) +
  geom_point(aes(), size = 5) +
  geom_errorbar(aes(ymax = happy + happy.se, ymin = happy - happy.se), width = .1) +
  labs(x = "", y = "How the story\nmade you feel") +
  theme(axis.title.y = element_text(vjust = 1)) +
  guides(linetype = "none",
         colour = "none",
         shape = "none") +
  ylim(2, 5) +
  scale_colour_manual(labels = plt.labs, values = plt.vals)

data.test = data %>% filter(condition == 'Positive') %>%
  group_by(imagine, subject) %>%
  summarize(j = mean(rating_justif), a = mean(rating_happy)) %>%
  arrange(subject, imagine) %>%
  group_by(subject) %>%
  mutate(j.diff = j - lag(j), a.diff = a - lag(a)) %>%
  filter(imagine == 'Imagine')
ggplot(data.test, aes(x = j.diff, y = a.diff)) +
  geom_point()

# test for effect of condition
model.affect <- lmer_alt(rating_happy ~ condition * imagine +
                              (imagine | subject) +
                              (condition * imagine | story_id),
                            data = data)
summary(rePCA(model.affect))
model.affect.r1 <- lmer_alt(rating_happy ~ condition * imagine +
                           (imagine | subject) +
                           (condition * imagine || story_id),
                         data = data)
summary(rePCA(model.affect.r1))
model.affect.r2 <- lmer_alt(rating_happy ~ condition * imagine +
                              (imagine | subject) +
                              (condition + imagine || story_id),
                            data = data)
summary(rePCA(model.affect.r2))
model.affect.r2.std <- lmer_alt(scale(rating_happy) ~ condition * imagine +
                              (imagine | subject) +
                              (condition + imagine || story_id),
                            data = data)
print.stats(model.affect.r2, model.affect.r2.std, 'conditionNegative:imagineImagine')
summary(model.affect.r2)

## simple effects
model.affect.pos <- lmer_alt(rating_happy ~  imagine +
                           (imagine | subject) +
                           (imagine | story_id),
                         data = data %>% filter(condition == 'Positive'))
summary(rePCA(model.affect.pos))
model.affect.pos.std <- lmer_alt(scale(rating_happy) ~  imagine +
                               (imagine | subject) +
                               (imagine | story_id),
                             data = data %>% filter(condition == 'Positive'))
print.stats(model.affect.pos, model.affect.pos.std, 'imagineImagine')

model.affect.neg <- lmer_alt(rating_happy ~  imagine +
                               (imagine | subject) +
                               (imagine | story_id),
                             data = data %>% filter(condition == 'Negative'))
summary(rePCA(model.affect.neg))
model.affect.neg.r1 <- lmer_alt(rating_happy ~  imagine +
                               (imagine | subject) +
                               (imagine || story_id),
                             data = data %>% filter(condition == 'Negative'))
summary(rePCA(model.affect.neg.r1))
model.affect.neg.r2 <- lmer_alt(rating_happy ~  imagine +
                                  (imagine | subject) +
                                  (1 | story_id),
                                data = data %>% filter(condition == 'Negative'))
summary(rePCA(model.affect.neg.r2))
model.affect.neg.r2.std <- lmer_alt(scale(rating_happy) ~  imagine +
                                   (imagine | subject) +
                                   (1 | story_id),
                                 data = data %>% filter(condition == 'Negative'))
print.stats(model.affect.neg.r2, model.affect.neg.r2.std, 'imagineImagine')

bf = lmBF(rating_happy ~ imagine + subject,
          data = data %>%
            filter(condition == 'Negative') %>%
            mutate(subject = factor(subject)),
          whichRandom="subject")
bf.null = lmBF(rating_will ~ subject,
               data = data %>%
                 filter(condition == 'Negative') %>%
                 mutate(subject = factor(subject)),
               whichRandom="subject")
bf.null / bf

## mediation analysis
detach(package:afex, unload=T)
detach(package:lmerTest, unload = T)
data.pos = data %>% filter(condition == 'Positive') %>%
  mutate(imagine.numeric = as.numeric(imagine))

model.med1 = lmer(rating_happy ~ imagine + rating_justif + (imagine.numeric | subject), data = data.pos)
summary(rePCA(model.med1))

model.med2 = lmer(rating_will ~ imagine + rating_happy + rating_justif +
                    (imagine.numeric + rating_happy | subject), data = data.pos)
summary(rePCA(model.med2))

model.med3 = mediate(model.med1, model.med2, treat = 'imagine', mediator = 'rating_happy')
summary(model.med3)
plot(model.med3)
#groundhog.library(c('lmerTest', 'afex'), '2021-04-01')
require(lmerTest)
require(afex)

summary(model.med1)
summary(model.med2)

# demographics ------------------------------------------------------------
demo <- read.csv(paste0('Study', expt, '_demo.csv')) %>% arrange(subject)
mean(demo$gender == 'Female')
mean(demo$age)
for (i in 1:nrow(data)) {
  subj = data$subject[i]
}

# basic
test1 = data %>% group_by(subject) %>%
  summarize(rating_will = mean(rating_will), rating_justif = mean(rating_justif),
            rating_happy.m = mean(rating_happy), rating_happy.se = se(rating_happy))
for (i in 1:nrow(test1)) {
  subj = test1$subject[i]
  test1$gender[i] = demo$gender[demo$subject == subj]
  test1$age[i] = demo$age[demo$subject == subj]
}
cor.test(test1$age, test1$rating_will)
cor.test(test1$age, test1$rating_justif)
cor.test(test1$age, test1$rating_happy)
test1.gender = test1 %>% group_by(gender) %>%
  summarize(rating_will = mean(rating_will), rating_justif = mean(rating_justif),
            rating_happy.m = mean(rating_happy.m), rating_happy.se = mean(rating_happy.se))
test1.gender

# imagine effect
test2 = data %>%
  group_by(imagine, subject) %>%
  summarize(will = mean(rating_will)) %>%
  arrange(subject) %>%
  group_by(subject) %>%
  mutate(will.diff = will - lag(will)) %>%
  filter(imagine == 'Imagine')
for (i in 1:nrow(test2)) {
  subj = test2$subject[i]
  test2$gender[i] = demo$gender[demo$subject == subj]
  test2$age[i] = demo$age[demo$subject == subj]
}
cor.test(test2$age, test2$will.diff)
test2.gender = test2 %>% group_by(gender) %>%
  summarize(will.diff = mean(will.diff))
test2.gender

# imagine x justification
test3 = data %>% 
  group_by(imagine, condition, subject) %>%
  summarize(will = mean(rating_will)) %>%
  arrange(subject, condition) %>%
  group_by(subject) %>%
  mutate(will.diff = will - lag(will)) %>%
  filter(imagine == 'Imagine') %>%
  ungroup() %>%
  mutate(will.diff2 = will.diff - lag(will.diff)) %>%
  filter(condition == 'Negative')
for (i in 1:nrow(test3)) {
  subj = test3$subject[i]
  test3$gender[i] = demo$gender[demo$subject == subj]
  test3$age[i] = demo$age[demo$subject == subj]
}
cor.test(test3$age, test3$will.diff2)
test3.gender = test3 %>% group_by(gender) %>%
  summarize(will.diff2 = mean(will.diff2))
test3.gender

# save --------------------------------------------------------------------

save.image(paste0('Study', expt, '_output.rdata'))

