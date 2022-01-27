# setup -------------------------------------------------------------------

# load packages with groundhog
require(groundhog)
pkgs = c('ggplot2', 'lme4', 'lmerTest',  'dplyr', 'afex', 'RColorBrewer', 'mediation', 'simr', 'mixedpower')
#groundhog.library(pkgs, '2021-04-01')
lapply(pkgs, require, character.only = T)

# Only works in RStudio -- otherwise you have to set the path manually
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

se <- function(data) {return(sd(data, na.rm = T)/sqrt(sum(!is.na(data))))}
se.prop <- function(data) {return(sqrt(mean(data,na.rm = T)*(1-mean(data, na.rm = T))/sum(!is.na(data))))}
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


# choose which experiment -------------------------------------------------

expt = readline(prompt="Which study do you want to analyze? (Enter 1a, 1b, or 1comb [for combined]): ")

# import data --------------------------------------------------------------------

data.raw <- read.csv(paste0('Study', expt, '_data.csv'))
data <- arrange(data.raw, subject)
data$imagine <- factor(data$imagine, labels = c("Control", "Imagine"), levels = c(0,1))
data$story_ind <- factor(data$story, labels = 1:10)


# exclusion & data prep ---------------------------------------------------------------

exclude.subj <- (data %>% group_by(subject) %>%
  summarise(nResponses = length(response), respLength = mean(nchar(encodeString(as.character(na.omit(response)))))) %>%
  filter(nResponses != 10 | respLength < 50))$subject

df.filt = data %>% filter(!(subject %in% exclude.subj)) %>%
  mutate(justif.high = factor(rating_justif >= median(rating_justif), c(F,T), c('Low', 'High')),
         justif.high.numeric = as.numeric(justif.high) - 1,
         imagine.numeric = as.numeric(imagine) - 1,
         justif.highscale = factor(rating_justif >= 4, c(F,T), c('Low', 'High')),
         justif.bin = cut(rating_justif, 6),
         justif.round = round(rating_justif),
         justif.bin.num = as.numeric(justif.bin),
         justif.bin.centered = justif.bin.num - mean(justif.bin.num),
         rating_justif.centered = rating_justif - mean(rating_justif))

# effect of trial type (imagine vs. control) on reported likelihood of harming ---------------------------------------------------------

df.bytrialtypeandsubject = df.filt %>%
  group_by(imagine, subject) %>%
  summarize(justif = mean(rating_justif), will = mean(rating_will),
            justif.high = mean(justif.high.numeric))
df.bytrialtype = df.bytrialtypeandsubject %>%
  group_by(imagine) %>%
  summarize(justif.m = mean(justif), justif.se = se(justif),
            justif.high.m = mean(justif.high), justif.high.se = se.prop(justif.high),
            will.m = mean(will), will.se = se(will))

ggplot(df.bytrialtype, aes(x = imagine, y = will.m)) +
  geom_point(aes(), size = 5) +
  geom_errorbar(aes(ymax = will.m + will.se, ymin = will.m - will.se), width = .1) +
  labs(x = "Trial Type", y = "Reported likelihood of\nperforming behavior") +
  guides(linetype = guide_legend(title = "Condition"),
         colour = guide_legend(title = "Condition"),
         shape = guide_legend(title = "Condition")) +
  scale_y_continuous(limits = c(2.9,4), breaks = c(3,4)) +
  theme(axis.text.x = element_text(size = 18))

# identify model
m.will = lmer_alt(rating_will ~ imagine +
                    (imagine | subject) +
                    (imagine | story_ind),
                  data = df.filt)
summary(rePCA(m.will))

m.will.r1 = lmer_alt(rating_will ~ imagine +
                       (imagine || subject) +
                       (imagine || story_ind),
                     data = df.filt)
summary(rePCA(m.will.r1))

if (expt == '1a') {
  m.will.r2 = lmer_alt(rating_will ~ imagine +
                         (1 | subject) +
                         (imagine || story_ind),
                       data = df.filt)
  summary(rePCA(m.will.r2))
  
  # get stats
  summary(m.will.r2)
  m.will.r2.std = lmer_alt(scale(rating_will) ~ imagine +
                         (1 | subject) +
                         (imagine || story_ind),
                       data = df.filt)
  print.stats(m.will.r2, m.will.r2.std, 'imagineImagine')
  
  m.will.r2.pwr = lmer(scale(rating_will) ~ imagine +
                         (1 | subject) +
                         (imagine | story_ind),
                         data = df.filt)
  summary(m.will.r2.pwr)
  pwr.will = mixedpower(m.will.r2.pwr, df.filt, c('imagine'), simvar = 'subject', steps = 96,
                        critical_value = 2, n_sim = 500)
  pwr.will
} else if (expt == '1b') {
  # get stats
  summary(m.will.r1)
  m.will.r1.std = lmer_alt(scale(rating_will) ~ imagine +
                             (imagine || subject) +
                             (imagine || story_ind),
                           data = df.filt)
  print.stats(m.will.r1, m.will.r1.std, 'imagineImagine')
  pwr.will = powerSim(m.will.r1, test = fixed('imagineImagine', 't'), nsim = 100)
  pwr.will
  # power is 67% (57, 76)
  
  # redo w/o item rfx for prereg
  m.will.noitem = lmer_alt(rating_will ~ imagine +
                             (imagine | subject),
                           data = df.filt)
  m.will.noitem.pca = rePCA(m.will.noitem)
  summary(m.will.noitem.pca)
  
  summary(m.will.noitem)
  m.will.noitem.std = lmer_alt(scale(rating_will) ~ imagine +
                             (imagine | subject),
                           data = df.filt)
  print.stats(m.will.noitem, m.will.noitem.std, 'imagineImagine')
}

# by subject
m.will.subj = lm(will.diff ~ j2, data = df.bytrialtypeandsubject)
summary(m.will.subj)
m.will.subj.std = lm(scale(will.diff) ~ scale(j2), data = df.bytrialtypeandsubject)
summary(m.will.subj.std)

m.will.subj2 = lm(will.diff ~ j2 + j2.diff, data = df.bytrialtypeandsubject)
summary(m.will.subj2)

df.bytrialtypeandsubject = df.bytrialtypeandsubject %>%
  arrange(subject) %>%
  group_by(subject) %>%
  mutate(will.diff = will - lag(will), j2.diff = j2 - lag(j2)) %>%
  filter(imagine == 'Imagine')
ggplot(df.bytrialtypeandsubject, aes(x = j2, y = will.diff)) +
  geom_point() +
  geom_smooth(method = 'lm', color = 'black') +
  #geom_abline(intercept = m.will.subj.int, slope = m.will.subj.slope) +
  #geom_abline(intercept = m.will.subj.int + 1.96*m.will.subj.se, slope = m.will.subj.slope, linetype = 'dashed') +
  scale_x_continuous(limits = c(1,7), breaks = c(1,4,7)) +
  scale_y_continuous(limits = c(-3.75, 4), breaks = c(-3,0,3)) +
  labs(x = 'Average justification', y = 'Average effect of\nimagination on likelihood')

# by subject 2 (editor suggestion)

test = df.filt %>% 
  group_by(subject) %>%
  mutate(justif.high.subj = rating_justif > median(rating_justif))

subjlist.test = unique(test$subject)
subj.coef = numeric(length(subjlist.test))
for (i in 1:length(subjlist.test)) {
  subj = subjlist.test[i]
  test.cur = test %>% filter(subject == subj)
  test.m = lm(rating_will ~ imagine * justif.high.subj, data = test.cur)
  coefs = summary(test.m)$coefficients
  if (nrow(coefs) == 4) {
    subj.coef[i] = summary(test.m)$coefficients[4,1]
  } else {
    subj.coef[i] = NA
  }
}


df.bytrialtypeandsubject2 = test %>%
  group_by(imagine, justif.high.subj, subject) %>%
  summarize(will = mean(rating_will)) %>%
  arrange(subject, justif.high.subj) %>%
  group_by(subject) %>%
  mutate(will.diff = will - lag(will)) %>%
  filter(imagine == 'Imagine') %>%
  mutate(will.diff2 = will.diff - lag(will.diff)) %>%
  filter(justif.high.subj == T)
hist(df.bytrialtypeandsubject2$will.diff2)
test.m = mean(df.bytrialtypeandsubject2$will.diff2)
test.se = se(df.bytrialtypeandsubject2$will.diff2)
c(test.m - 1.96*test.se, test.m, test.m + 1.96*test.se)
mean(df.bytrialtypeandsubject2$will.diff2 > 0)
se.prop(df.bytrialtypeandsubject2$will.diff2 > 0)

## explore justification
# does it differ between conditions?
ggplot(df.bytrialtype, aes(x = imagine, y = j2.m)) +
  geom_point(aes(), size = 5) +
  geom_errorbar(aes(ymax = j2.m + j2.se, ymin = j2.m - j2.se), width = .1) +
  labs(x = "Trial Type", y = "Justification\nrating") +
  scale_y_continuous(limits = c(2.9,4), breaks = c(3,4)) +
  theme(axis.text.x = element_text(size = 18))

ggplot(df.bytrialtype, aes(x = imagine, y = j2.high.m)) +
  geom_point(aes(), size = 5) +
  geom_errorbar(aes(ymax = j2.high.m + j2.high.se, ymin = j2.high.m - j2.high.se), width = .1) +
  labs(x = "Trial Type", y = "% justification\nabove median") +
  theme(axis.text.x = element_text(size = 18))

# identify model
m.justif = lmer_alt(rating_justif ~ imagine +
                      (imagine | subject) +
                      (imagine | story_ind),
                    data = df.filt)
m.justif.pca = rePCA(m.justif)
summary(m.justif.pca)

if (expt == '2a') {
  m.justif.r1 = lmer_alt(rating_justif ~ imagine +
                        (imagine || subject) +
                        (imagine || story_ind),
                      data = df.filt)
  m.justif.r1.pca = rePCA(m.justif.r1)
  summary(m.justif.r1.pca)
  
  m.justif.r2 = lmer_alt(rating_justif ~ imagine +
                           (1 | subject) +
                           (1 | story_ind),
                         data = df.filt)
  m.justif.r2.pca = rePCA(m.justif.r2)
  summary(m.justif.r2.pca)
  
  # get stats
  summary(m.justif.r2)
  m.justif.r2.std = lmer_alt(scale(rating_justif) ~ imagine +
                           (1 | subject) +
                           (1 | story_ind),
                         data = df.filt)
  print.stats(m.justif.r2, m.justif.r2.std, 'imagineImagine')
} else if (expt == '2b' || expt == '2comb') {
  summary(m.justif)
  m.justif.std = lmer_alt(scale(rating_justif) ~ imagine +
                               (imagine | subject) +
                               (imagine | story_ind),
                             data = df.filt)
  print.stats(m.justif, m.justif.std, 'imagineImagine')
}

## moderation

df.bytrialtypeandjustif = df.filt %>%
  group_by(imagine, justif.high, subject) %>%
  summarize(will = mean(rating_will)) %>%
  group_by(imagine, justif.high) %>%
  summarize(will.m = mean(will), will.se = se(will))

ggplot(df.bytrialtypeandjustif, aes(x = justif.high, y = will.m, group = imagine, fill = imagine)) + 
  geom_col(position = dodge) +
  geom_errorbar(aes(ymin = will.m - will.se, ymax = will.m + will.se), width = .1, position = dodge) +
  labs(x = "Relative justification\n(split by median)", y = "Reported likelihood of\nperforming behavior") +
  scale_y_continuous(limits = c(0,5), breaks = c(0,5)) +
  guides(fill = guide_legend(title = "Trial Type")) +
  scale_fill_brewer() +
  theme(legend.position = 'bottom')

df.bytrialtypeandjustifbinned = df.filt %>%
  group_by(imagine, justif.bin, subject) %>%
  summarize(will = mean(rating_will)) %>%
  group_by(imagine, justif.bin) %>%
  summarize(will.m = mean(will), will.se = se(will)) %>%
  arrange(justif.bin) %>%
  group_by(justif.bin) %>%
  mutate(will.diff = will.m - lag(will.m), will.diff.se = sqrt(will.se^2 + lag(will.se)^2))

ggplot(df.bytrialtypeandjustifbinned %>% filter(imagine == 'Imagine'), aes(x = justif.bin, y = will.diff)) + 
  geom_col(position = dodge) +
  geom_errorbar(aes(ymin = will.diff - will.diff.se, ymax = will.diff + will.diff.se), width = .1, position = dodge) +
  labs(x = "Absolute justification", y = "Difference in likelihood\n(imagine - control)") +
  scale_y_continuous(breaks = c(0,0.5)) +
  guides(fill = guide_legend(title = "Trial Type")) +
  scale_fill_brewer() +
  geom_vline(xintercept = median(df.filt$rating_justif), color = 'red', linetype = 'dashed') +
  scale_x_discrete(labels = c('1-2', '2-3', '3-4', '4-5', '5-6', '6-7')) +
  theme(axis.text=element_text(size=18, colour = "black"),
axis.title=element_text(size=18, face = "bold"))

ggplot(df.bytrialtypeandjustifbinned, aes(x = justif.bin, y = will.m, fill = imagine)) + 
  geom_col(position = dodge) +
  geom_errorbar(aes(ymin = will.m - will.se, ymax = will.m + will.se), width = .1, position = dodge) +
  labs(x = "Absolute justification", y = "Difference in reported likelihoods\n(imagine - control)") +
  #scale_y_continuous(limits = c(-1,2), breaks = c(-1,2)) +
  guides(fill = guide_legend(title = "Trial Type")) +
  scale_fill_brewer() +
  geom_vline(xintercept = median(df.filt$rating_justif), color = 'red', linetype = 'dashed') +
  scale_x_discrete(labels = c('1-2', '2-3', '3-4', '4-5', '5-6', '6-7'))

ggplot(df.filt, aes(x = justif.bin)) +
  geom_histogram(stat='count') +
  scale_x_discrete(labels = c('1-2', '2-3', '3-4', '4-5', '5-6', '6-7')) +
  scale_y_discrete(breaks = NULL) +
  geom_vline(xintercept = median(df.filt$rating_justif), color = 'red', linetype = 'dashed') +
  labs(x = 'Justification', y = 'Number of\nTrials')
m = mean(df.filt$rating_justif > 4)
s = se.prop(df.filt$rating_justif > 4)
c(m-1.96*s,m+1.96*s)

# identify model
m.mod = lmer_alt(rating_will ~ justif.high * imagine +
               (justif.high * imagine | subject) +
               (justif.high * imagine | story_ind),
             data = df.filt)
m.mod.pca = rePCA(m.mod)
summary(m.mod.pca)

m.mod.r1 = lmer_alt(rating_will ~ justif.high * imagine +
                          (justif.high * imagine || subject) +
                          (justif.high * imagine || story_ind),
                        data = df.filt)
m.mod.r1.pca = rePCA(m.mod.r1)
summary(m.mod.r1.pca)

if (expt == '2a') {
  m.mod.r2 = lmer_alt(rating_will ~ justif.high * imagine +
                            (justif.high * imagine || subject) +
                            (justif.high + imagine || story_ind),
                          data = df.filt)
  m.mod.r2.pca = rePCA(m.mod.r2)
  summary(m.mod.r2.pca)
  
  m.mod.r3 = lmer_alt(rating_will ~ justif.high * imagine +
                               (justif.high * imagine || subject) +
                               (justif.high || story_ind),
                             data = df.filt)
  m.mod.uncorr.r3.pca = rePCA(m.mod.r3)
  summary(m.mod.uncorr.r3.pca)
  
  m.mod.r3.alt = lmer_alt(rating_will ~ justif.high * imagine +
                        (justif.high * imagine || subject) +
                        (imagine || story_ind),
                      data = df.filt)
  m.mod.r3.alt.pca = rePCA(m.mod.r3.alt)
  summary(m.mod.r3.alt.pca)
  
  m.mod.r3.recorr = lmer_alt(rating_will ~ justif.high * imagine +
                               (justif.high * imagine || subject) +
                               (justif.high | story_ind),
                             data = df.filt)
  m.mod.r3.recorr.pca = rePCA(m.mod.r3.recorr)
  summary(m.mod.r3.recorr.pca)
  
  # get stats
  summary(m.mod.r3.recorr)
  m.mod.r3.recorr.std = lmer_alt(scale(rating_will) ~ justif.high * imagine +
                               (justif.high * imagine || subject) +
                               (justif.high | story_ind),
                             data = df.filt)
  print.stats(m.mod.r3.recorr, m.mod.r3.recorr.std, 'justif.highHigh:imagineImagine')

  m.mod.r3.recorr.pwr = lmer(scale(rating_will) ~ justif.high * imagine +
                                   (justif.high * imagine | subject) +
                                   (justif.high | story_ind),
                                 data = df.filt)
  summary(m.mod.r3.recorr.pwr)
  pwr.mod = mixedpower(m.mod.r3.recorr.pwr, df.filt, c('justif.high', 'imagine'), simvar = 'subject', steps = 96,
                        critical_value = 2, n_sim = 500, databased = T)
  pwr.mod
  
  # excluding failed sims
  
  m.mod.r3.recorr.nofail = lmer_alt(rating_will ~ justif.high * imagine +
                               (justif.high * imagine || subject) +
                               (justif.high | story_ind),
                             data = df.filt %>% filter(obs.failed == F))
  summary(m.mod.r3.recorr.nofail)
  
} else if (expt == '2b') {
  m.mod.r2 = lmer_alt(rating_will ~ justif.high * imagine +
                        (justif.high + imagine || subject) +
                        (justif.high + imagine || story_ind),
                      data = df.filt)
  m.mod.r2.pca = rePCA(m.mod.r2)
  summary(m.mod.r2.pca)
  
  m.mod.r3 = lmer_alt(rating_will ~ justif.high * imagine +
                        (justif.high + imagine || subject) +
                        (justif.high || story_ind),
                      data = df.filt)
  m.mod.r3.pca = rePCA(m.mod.r3)
  summary(m.mod.r3.pca)
  
  m.mod.r4 = lmer_alt(rating_will ~ justif.high * imagine +
                        (justif.high + imagine || subject) +
                        (imagine || story_ind),
                      data = df.filt)
  m.mod.r4.pca = rePCA(m.mod.r4)
  summary(m.mod.r4.pca)
  
  m.mod.r5 = lmer_alt(rating_will ~ justif.high * imagine +
                        (justif.high + imagine | subject) +
                        (justif.high | story_ind),
                      data = df.filt)
  m.mod.r5.pca = rePCA(m.mod.r5)
  summary(m.mod.r5.pca)
  
  # stats
  summary(m.mod.r5)
  m.mod.r5.std = lmer_alt(scale(rating_will) ~ justif.high * imagine +
                        (justif.high + imagine | subject) +
                        (justif.high | story_ind),
                      data = df.filt)
  print.stats(m.mod.r5, m.mod.r5.std, 'justif.highHigh:imagineImagine')
  
  powerSim(m.mod.r5, test = fixed('justif.highHigh:imagineImagine', 't'), nsim = 100)
  # power is 78% [69%, 86%]
  
  # redo w/o item rfx for prereg
  m.mod.noitem = lmer_alt(rating_will ~ justif.high * imagine +
                     (justif.high * imagine | subject),
                   data = df.filt)
  m.mod.noitem.pca = rePCA(m.mod.noitem)
  summary(m.mod.noitem.pca)
  m.mod.noitem.r1 = lmer_alt(rating_will ~ justif.high * imagine +
                            (justif.high * imagine || subject),
                          data = df.filt)
  m.mod.noitem.r1.pca = rePCA(m.mod.noitem.r1)
  summary(m.mod.noitem.r1.pca)
  m.mod.noitem.r2 = lmer_alt(rating_will ~ justif.high * imagine +
                               (justif.high + imagine || subject),
                             data = df.filt)
  m.mod.noitem.r2.pca = rePCA(m.mod.noitem.r2)
  summary(m.mod.noitem.r2.pca)
  m.mod.noitem.r3 = lmer_alt(rating_will ~ justif.high * imagine +
                               (justif.high + imagine | subject),
                             data = df.filt)
  m.mod.noitem.r3.pca = rePCA(m.mod.noitem.r3)
  summary(m.mod.noitem.r3.pca)
  
  summary(m.mod.noitem.r3)
  m.mod.noitem.r3.std = lmer_alt(scale(rating_will) ~ justif.high * imagine +
                               (justif.high + imagine | subject),
                             data = df.filt)
  print.stats(m.mod.noitem.r3, m.mod.noitem.r3.std, 'justif.highHigh:imagineImagine')
}

# simple effects
m.mod.high = lmer_alt(rating_will ~ imagine +
                   (imagine | subject) +
                   (imagine | story_ind),
                 data = df.filt %>% filter(justif.high == 'High'))
summary(rePCA(m.mod.high))
m.mod.high.r1 = lmer_alt(rating_will ~ imagine +
                        (imagine || subject) +
                        (imagine || story_ind),
                      data = df.filt %>% filter(justif.high == 'High'))
summary(rePCA(m.mod.high.r1))
m.mod.high.r2 = lmer_alt(rating_will ~ imagine +
                           (1 | subject) +
                           (1 | story_ind),
                         data = df.filt %>% filter(justif.high == 'High'))
summary(rePCA(m.mod.high.r2))
summary(m.mod.high.r2)
m.mod.high.r2.std = lmer_alt(scale(rating_will) ~ imagine +
                           (1 | subject) +
                           (1 | story_ind),
                         data = df.filt %>% filter(justif.high == 'High'))
print.stats(m.mod.high.r2, m.mod.high.r2.std, 'imagineImagine')


m.mod.low = lmer_alt(rating_will ~ imagine +
                        (imagine | subject) +
                        (imagine | story_ind),
                      data = df.filt %>% filter(justif.high == 'Low'))
summary(rePCA(m.mod.low))
m.mod.low.r1 = lmer_alt(rating_will ~ imagine +
                           (imagine || subject) +
                           (imagine || story_ind),
                         data = df.filt %>% filter(justif.high == 'Low'))
summary(rePCA(m.mod.low.r1))
m.mod.low.r2 = lmer_alt(rating_will ~ imagine +
                           (imagine || subject) +
                           (1 | story_ind),
                         data = df.filt %>% filter(justif.high == 'Low'))
summary(rePCA(m.mod.low.r2))
summary(m.mod.low.r2)
m.mod.low.r2.std = lmer_alt(scale(rating_will) ~ imagine +
                          (imagine || subject) +
                          (1 | story_ind),
                        data = df.filt %>% filter(justif.high == 'Low'))
print.stats(m.mod.low.r2, m.mod.low.r2.std, 'imagineImagine')


# continuous: full
df.filt = df.filt %>% mutate(
  justif.mirrored = rating_justif - 4,
  justif.mirrored.abs = abs(justif.mirrored),
  justif.mirrored.abs.centered = justif.mirrored.abs - mean(justif.mirrored.abs),
  justif.mirrored.sq = justif.mirrored ^ 2,
  justif.mirrored.sq.centered = justif.mirrored.sq - mean(justif.mirrored.sq),
  justif.bin.mirrored = justif.bin.num - 3.5,
  justif.bin.mirrored.sq = justif.bin.mirrored ^ 2,
  justif.bin.mirrored.abs = abs(justif.mirrored),
  justif.bin.mirrored.centered = justif.bin.mirrored - mean(justif.bin.mirrored),
  justif.bin.mirrored.sq.centered = justif.bin.mirrored.sq - mean(justif.bin.mirrored.sq),
  justif.bin.mirrored.abs.centered = justif.bin.mirrored.abs - mean(justif.bin.mirrored.abs))

m.mod.cont = lmer_alt(rating_will ~ rating_justif.centered * justif.mirrored.sq.centered * imagine + 
                        (rating_justif.centered + imagine + justif.mirrored.sq.centered || subject) +
                        (rating_justif.centered * imagine + justif.mirrored.sq.centered || story_ind),
                      data = df.filt)
summary(rePCA(m.mod.cont))
summary(m.mod.cont)
m.mod.cont.std = lmer_alt(rating_will ~ rating_justif.centered * justif.mirrored.sq.centered * imagine + 
                        (rating_justif.centered + imagine + justif.mirrored.sq.centered || subject) +
                        (rating_justif.centered * imagine + justif.mirrored.sq.centered || story_ind),
                      data = df.filt %>% mutate(rating_will = scale(rating_will),
                                                rating_justif.centered = scale(rating_justif.centered),
                                                justif.mirrored.sq.centered = scale(justif.mirrored.sq.centered)))
print.stats(m.mod.cont, m.mod.cont.std, 'rating_justif.centered:justif.mirrored.sq.centered:imagineImagine')

ggplot(df.filt, aes(x = rating_justif, y = rating_will, color = imagine)) +
  geom_point() +
  stat_smooth(method='lm', formula = y ~ x)

# continuous: simple effects
m.mod.bin1 = lmer_alt(rating_will ~ imagine +
                   (imagine | subject) +
                   (imagine || story_ind),
                 data = df.filt %>% filter(rating_justif >= 1 & rating_justif < 2))
summary(rePCA(m.mod.bin1))
summary(m.mod.bin1)
m.mod.bin1.std = lmer_alt(scale(rating_will) ~ imagine +
                        (imagine | subject) +
                        (imagine || story_ind),
                      data = df.filt %>% filter(rating_justif >= 1 & rating_justif < 2))
print.stats(m.mod.bin1, m.mod.bin1.std, 'imagineImagine')

m.mod.bin2 = lmer_alt(rating_will ~ imagine +
                        (imagine | subject) +
                        (imagine | story_ind),
                      data = df.filt %>% filter(rating_justif >= 2 & rating_justif < 3))
summary(rePCA(m.mod.bin2))
summary(m.mod.bin2)
m.mod.bin2.std = lmer_alt(scale(rating_will) ~ imagine +
                        (imagine | subject) +
                        (imagine | story_ind),
                      data = df.filt %>% filter(rating_justif >= 2 & rating_justif < 3))
print.stats(m.mod.bin2, m.mod.bin2.std, 'imagineImagine')

m.mod.bin3 = lmer_alt(rating_will ~ imagine +
                        (imagine | subject) +
                        (imagine | story_ind),
                      data = df.filt %>% filter(rating_justif >= 3 & rating_justif < 4))
m.mod.bin3.r1 = lmer_alt(rating_will ~ imagine +
                        (imagine || subject) +
                        (imagine || story_ind),
                      data = df.filt %>% filter(rating_justif >= 3 & rating_justif < 4))
summary(m.mod.bin3.r1)
m.mod.bin3.r1.std = lmer_alt(scale(rating_will) ~ imagine +
                        (imagine || subject) +
                        (imagine || story_ind),
                      data = df.filt %>% filter(rating_justif >= 3 & rating_justif < 4))
print.stats(m.mod.bin3.r1, m.mod.bin3.r1.std, 'imagineImagine')


m.mod.bin4 = lmer_alt(rating_will ~ imagine +
                        (imagine | subject) +
                        (imagine | story_ind),
                      data = df.filt %>% filter(rating_justif >= 4 & rating_justif < 5))
m.mod.bin4.r1 = lmer_alt(rating_will ~ imagine +
                           (imagine || subject) +
                           (imagine || story_ind),
                         data = df.filt %>% filter(rating_justif >= 4 & rating_justif < 5))
summary(rePCA(m.mod.bin4.r1))
m.mod.bin4.r2 = lmer_alt(rating_will ~ imagine +
                           (1 | subject) +
                           (1 | story_ind),
                         data = df.filt %>% filter(rating_justif >= 4 & rating_justif < 5))
summary(rePCA(m.mod.bin4.r2))
summary(m.mod.bin4.r2)
m.mod.bin4.r2.std = lmer_alt(scale(rating_will) ~ imagine +
                           (1 | subject) +
                           (1 | story_ind),
                         data = df.filt %>% filter(rating_justif >= 4 & rating_justif < 5))
print.stats(m.mod.bin4.r2, m.mod.bin4.r2.std, 'imagineImagine')

m.mod.bin5.r1 = lmer_alt(rating_will ~ imagine +
                           (imagine || subject) +
                           (imagine || story_ind),
                         data = df.filt %>% filter(rating_justif >= 5 & rating_justif < 6))
summary(rePCA(m.mod.bin5.r1))
m.mod.bin5.r2 = lmer_alt(rating_will ~ imagine +
                           (imagine || subject) +
                           (1 | story_ind),
                         data = df.filt %>% filter(rating_justif >= 5 & rating_justif < 6))
summary(rePCA(m.mod.bin5.r2))
summary(m.mod.bin5.r2)
m.mod.bin5.r2.std = lmer_alt(scale(rating_will) ~ imagine +
                           (imagine || subject) +
                           (1 | story_ind),
                         data = df.filt %>% filter(rating_justif >= 5 & rating_justif < 6))
print.stats(m.mod.bin5.r2, m.mod.bin5.r2.std, 'imagineImagine')

m.mod.bin6.r1 = lmer_alt(rating_will ~ imagine +
                           (imagine || subject) +
                           (imagine || story_ind),
                         data = df.filt %>% filter(rating_justif >= 6 & rating_justif <= 7))
summary(rePCA(m.mod.bin6.r1))
m.mod.bin6.r2 = lmer_alt(rating_will ~ imagine +
                           (imagine || subject) +
                           (1 | story_ind),
                         data = df.filt %>% filter(rating_justif >= 6 & rating_justif <= 7))
summary(rePCA(m.mod.bin6.r2))
m.mod.bin6.r3 = lmer_alt(rating_will ~ imagine +
                           (imagine || subject),
                         data = df.filt %>% filter(rating_justif >= 6 & rating_justif <= 7))
summary(m.mod.bin6.r3)
m.mod.bin6.r3.std = lmer_alt(scale(rating_will) ~ imagine +
                           (imagine || subject),
                         data = df.filt %>% filter(rating_justif >= 6 & rating_justif < 7))
print.stats(m.mod.bin6.r3, m.mod.bin6.r3.std, 'imagineImagine')


bin.coefs = c(summary(m.mod.bin1.std)$coefficients[2,1],
              summary(m.mod.bin2.std)$coefficients[2,1],
              summary(m.mod.bin3.r1.std)$coefficients[2,1],
              summary(m.mod.bin4.r2.std)$coefficients[2,1],
              summary(m.mod.bin5.r2.std)$coefficients[2,1],
              summary(m.mod.bin6.r3.std)$coefficients[2,1])
bin.ses = c(summary(m.mod.bin1.std)$coefficients[2,2],
              summary(m.mod.bin2.std)$coefficients[2,2],
              summary(m.mod.bin3.r1.std)$coefficients[2,2],
              summary(m.mod.bin4.r2.std)$coefficients[2,2],
              summary(m.mod.bin5.r2.std)$coefficients[2,2],
              summary(m.mod.bin6.r3.std)$coefficients[2,2])
bin.df = data.frame(bin = 1:6, coef = bin.coefs, se = bin.ses)


ggplot(bin.df, aes(x = factor(bin), y = coef)) + 
  geom_col(position = dodge) +
  geom_errorbar(aes(ymin = coef - se, ymax = coef + se), width = .1, position = dodge) +
  labs(x = "Justification level", y = "Effect of imagination\non reported harm\n(std. beta)") +
  scale_y_continuous(breaks = c(-.25,0,0.25)) +
  guides(fill = guide_legend(title = "Trial Type")) +
  scale_fill_brewer() +
  geom_vline(xintercept = median(df.filt$rating_justif), color = 'red', linetype = 'dashed', size = 1.2) +
  scale_x_discrete(labels = c('1-2', '2-3', '3-4', '4-5', '5-6', '6-7')) +
  theme(axis.text=element_text(size=18, colour = "black"),
        axis.title=element_text(size=18, face = "bold"))

# extra analyses for SM on combined data
if (expt == '2comb') {
  ## moderation via justification, split by scale midpoint
  df.bytrialtypeandjustifscale = df.filt %>%
    group_by(imagine, justif.highscale, subject) %>%
    summarize(will = mean(rating_will)) %>%
    group_by(imagine, justif.highscale) %>%
    summarize(will.m = mean(will), will.se = se(will))
  
  ggplot(df.bytrialtypeandjustifscale, aes(x = justif.highscale, y = will.m, group = imagine, fill = imagine)) + 
    geom_col(position = dodge) +
    geom_errorbar(aes(ymin = will.m - will.se, ymax = will.m + will.se), width = .1, position = dodge) +
    labs(x = "Justification rating\n(split by scale midpoint)", y = "Reported likelihood of\nperforming behavior") +
    #scale_y_continuous(limits = c(0,5), breaks = c(0,5)) +
    guides(fill = guide_legend(title = "Trial Type")) +
    scale_fill_brewer()
  
  # identify interaction model
  m.mod.scale = lmer_alt(rating_will ~ justif.highscale * imagine +
                     (justif.highscale * imagine | subject) +
                     (justif.highscale * imagine | story_ind),
                   data = df.filt)
  summary(rePCA(m.mod.scale))
  m.mod.scale.r1 = lmer_alt(rating_will ~ justif.highscale * imagine +
                           (justif.highscale * imagine || subject) +
                           (justif.highscale * imagine || story_ind),
                         data = df.filt)
  summary(rePCA(m.mod.scale.r1))
  
  summary(m.mod.scale.r1)
  m.mod.scale.r1.std = lmer_alt(scale(rating_will) ~ justif.highscale * imagine +
                              (justif.highscale * imagine || subject) +
                              (justif.highscale * imagine || story_ind),
                            data = df.filt)
  print.stats(m.mod.scale.r1, m.mod.scale.r1.std, 'justif.highscaleHigh:imagineImagine')
  
  # high only
  m.mod.scale.high = lmer_alt(rating_will ~ imagine +
                           (imagine | subject) +
                           (imagine | story_ind),
                         data = df.filt %>% filter(justif.highscale == 'High'))
  summary(rePCA(m.mod.scale.high))
  m.mod.scale.high.r1 = lmer_alt(rating_will ~ imagine +
                                (imagine | subject) +
                                (imagine || story_ind),
                              data = df.filt %>% filter(justif.highscale == 'High'))
  summary(rePCA(m.mod.scale.high.r1))
  m.mod.scale.high.r2 = lmer_alt(rating_will ~ imagine +
                                   (imagine | subject) +
                                   (1 | story_ind),
                                 data = df.filt %>% filter(justif.highscale == 'High'))
  summary(rePCA(m.mod.scale.high.r2))
  
  summary(m.mod.scale.high.r2)
  m.mod.scale.high.r2.std = lmer_alt(scale(rating_will) ~ imagine +
                                   (imagine | subject) +
                                   (1 | story_ind),
                                 data = df.filt %>% filter(justif.highscale == 'High'))
  print.stats(m.mod.scale.high.r2, m.mod.scale.high.r2.std, 'imagineImagine')
  
  # low only
  m.mod.scale.low = lmer_alt(rating_will ~ imagine +
                                (imagine | subject) +
                                (imagine | story_ind),
                              data = df.filt %>% filter(justif.highscale == 'Low'))
  summary(rePCA(m.mod.scale.low))
  
  summary(m.mod.scale.low)
  m.mod.scale.low.std = lmer_alt(scale(rating_will) ~ imagine +
                               (imagine | subject) +
                               (imagine | story_ind),
                             data = df.filt %>% filter(justif.highscale == 'Low'))
  print.stats(m.mod.scale.low, m.mod.scale.low.std, 'imagineImagine')
  
  ## moderation via justification, continuous
  df.bytrialtypeandjustifbin = df.filt %>%
    group_by(imagine, justif.bin, subject) %>%
    summarize(will = mean(rating_will)) %>%
    group_by(imagine, justif.bin) %>%
    summarize(will.m = mean(will), will.se = se(will))
  
  ggplot(df.bytrialtypeandjustifbin, aes(x = justif.bin, y = will.m, group = imagine, fill = imagine)) + 
    geom_col(position = dodge) +
    geom_errorbar(aes(ymin = will.m - will.se, ymax = will.m + will.se), width = .1, position = dodge) +
    labs(x = "Justification rating\n(rounded to nearest integer)", y = "Reported likelihood of\nperforming behavior") +
    #scale_y_continuous(limits = c(0,5), breaks = c(0,5)) +
    guides(fill = guide_legend(title = "Trial Type")) +
    scale_fill_brewer()
  
  hist(as.numeric(df.filt$rating_justif))
  mean(df.filt$justif.bin == '(6,7.01]')
  mean(df.filt[df.filt$justif.bin == '(6,7.01]',]$rating_will > 6.5)
  
  # identify model
  # starting with a less-than-full rfx structure here, b/c very likely the full model won't converge
  m.mod.bin = lmer_alt(rating_will ~ 
                         rating_justif.centered * I(rating_justif.centered^2) * imagine +
                          (rating_justif.centered + I(rating_justif.centered^2) + imagine || subject) +
                          (rating_justif.centered + I(rating_justif.centered^2) + imagine || story_ind),
                       data = df.filt)
  summary(rePCA(m.mod.bin))
  m.mod.bin.r1 = lmer_alt(rating_will ~ 
                         rating_justif.centered * I(rating_justif.centered^2) * imagine +
                         (rating_justif.centered + I(rating_justif.centered^2) || subject) +
                         (rating_justif.centered + I(rating_justif.centered^2) + imagine || story_ind),
                       data = df.filt)
  summary(rePCA(m.mod.bin.r1))
  m.mod.bin.r2 = lmer_alt(rating_will ~ 
                            rating_justif.centered * I(rating_justif.centered^2) * imagine +
                            (rating_justif.centered + imagine || subject) +
                            (rating_justif.centered + I(rating_justif.centered^2) + imagine || story_ind),
                          data = df.filt)
  summary(rePCA(m.mod.bin.r2))
  m.mod.bin.r3 = lmer_alt(rating_will ~ 
                            rating_justif.centered * I(rating_justif.centered^2) * imagine +
                            (I(rating_justif.centered^2) + imagine || subject) +
                            (rating_justif.centered + I(rating_justif.centered^2) + imagine || story_ind),
                          data = df.filt)
  summary(rePCA(m.mod.bin.r3))
  m.mod.bin.r4 = lmer_alt(rating_will ~ 
                            rating_justif.centered * I(rating_justif.centered^2) * imagine +
                            (rating_justif.centered || subject) +
                            (rating_justif.centered + I(rating_justif.centered^2) + imagine || story_ind),
                          data = df.filt)
  summary(rePCA(m.mod.bin.r4))
  m.mod.bin.r5 = lmer_alt(rating_will ~ 
                            rating_justif.centered * I(rating_justif.centered^2) * imagine +
                            (I(rating_justif.centered^2) || subject) +
                            (rating_justif.centered + I(rating_justif.centered^2) + imagine || story_ind),
                          data = df.filt)
  summary(rePCA(m.mod.bin.r5))
  m.mod.bin.r6 = lmer_alt(rating_will ~ 
                            rating_justif.centered * I(rating_justif.centered^2) * imagine +
                            (imagine || subject) +
                            (rating_justif.centered + I(rating_justif.centered^2) + imagine || story_ind),
                          data = df.filt)
  summary(rePCA(m.mod.bin.r6))
  
  summary(m.mod.bin.r6)
  m.mod.bin.r6.std = lmer_alt(scale(rating_will) ~ 
                            rating_justif.centered * I(rating_justif.centered^2) * imagine +
                            (imagine || subject) +
                            (rating_justif.centered + I(rating_justif.centered^2) + imagine || story_ind),
                          data = df.filt %>% mutate(rating_justif.centered = scale(rating_justif.centered)))
  print.stats(m.mod.bin.r6, m.mod.bin.r6.std, 'rating_justif.centered:I(rating_justif.centered^2):imagineImagine')
  print.stats(m.mod.bin.r6, m.mod.bin.r6.std, 'rating_justif.centered:imagineImagine')
  print.stats(m.mod.bin.r6, m.mod.bin.r6.std, 'I(rating_justif.centered^2):imagineImagine')
}

## mediation

detach(package:afex, unload = T)

# model 1
model.med1 = glmer(justif.high.numeric ~ imagine +
                     (imagine.numeric | subject),
                   family = 'binomial',
                   data = df.filt)
summary(rePCA(model.med1))

# model 2
model.med2 = lmer(rating_will ~ imagine + justif.high.numeric +
                    (imagine.numeric + justif.high.numeric | subject),
                  data = df.filt)
summary(rePCA(model.med2))

if (expt == '2a') {
  model.med1.r1 = glmer(justif.high.numeric ~ imagine +
                          (1 | subject),
                        family = 'binomial',
                        data = df.filt)
  summary(rePCA(model.med1.r1))
  
  # model 2
  model.med2.r1 = lmer(rating_will ~ imagine + justif.high.numeric +
                         (justif.high.numeric | subject),
                       data = df.filt)
  summary(rePCA(model.med2.r1))
  
  model.med2.r2 = lmer(rating_will ~ imagine + justif.high.numeric +
                         (justif.high.numeric | subject),
                       data = df.filt)
  summary(rePCA(model.med2.r2))
  
  # mediate
  model.med3 = mediate(model.med1.r1, model.med2.r1, treat = 'imagine', mediator = 'justif.high.numeric')
  summary(model.med3)
  
  summary(model.med1.r1)
  summary(model.med2.r1)
} else if (expt == '2b') {
  model.med3 = mediate(model.med1, model.med2, treat = 'imagine', mediator = 'justif.high.numeric')
  summary(model.med3)
  summary(model.med1)
  summary(model.med2)
}
require(afex)


# ratings analysis --------------------------------------------------------

df.imagharm = df.filt %>% filter(imagine == 'Imagine')

# get ratings
df.justif.all = read.csv('Study1a_ratings.csv') %>% mutate(resp_num = resp_num + 1)

df.justif = df.justif.all %>% # b/c of zero-indexing 
  group_by(resp_num) %>%
  summarize(rating1.pctyes = mean(rating1 == 1), rating1.pctno = mean(rating1 == 2), rating1.pctfail = mean(rating1 == 3),
            rating2.mean = mean(rating2), rating2.se = se(rating2)) %>%
  arrange(resp_num) %>%
  mutate(rating.bin = cut(rating2.mean, breaks = c(1,2,3,4,5,6,7)),
         rating.high = rating2.mean >= 4, rating.high.med = rating2.mean > median(rating2.mean))

df.imagharm$obs.justif.pctyes = df.justif$rating1.pctyes
df.imagharm$obs.justif.pctno = df.justif$rating1.pctno
df.imagharm$obs.justif.pctfail = df.justif$rating1.pctfail
df.imagharm$obs.justif.mean = df.justif$rating2.mean
df.imagharm$obs.justif.se = df.justif$rating2.se
df.imagharm$obs.justif.high = df.justif$rating.high
df.imagharm$obs.justif.bin = df.justif$rating.bin

for (i in 1:nrow(df.justif.all)) {
  orig.row = df.justif.all$resp_num[i]
  df.justif.all$orig.justif[i] = df.imagharm$rating_justif[orig.row]
  df.justif.all$orig.will[i] = df.imagharm$rating_will[orig.row]
  df.justif.all$orig.subject[i] = df.imagharm$subject[orig.row]
  df.justif.all$orig.story[i] = df.imagharm$story_id[orig.row]
}

# for (i in 1:nrow(df.imagharm)) {
#   justif.row = which(df.justif$resp_num == i)
#   if (length(justif.row) > 0) {
#     df.imagharm$obs.justif.pctyes[i] = df.justif$rating1.pctyes[justif.row]
#     df.imagharm$obs.justif.pctno[i] = df.justif$rating1.pctno[justif.row]
#     df.imagharm$obs.justif.pctfail[i] = df.justif$rating1.pctfail[justif.row]
#     df.imagharm$obs.justif.mean[i] = df.justif$rating2.mean[justif.row]
#     df.imagharm$obs.justif.bin[i] = df.justif$rating2.bin[justif.row]
#     df.imagharm$obs.justif.se[i] = df.justif$rating2.se[justif.row]
#     df.imagharm$obs.justif.high[i] = df.justif$rating.high.med[justif.row]
#   }
# }

ggplot(df.imagharm, aes(x = rating_justif, y = obs.justif.mean)) +
  geom_point() +
  geom_smooth(method='lm')
m.ratings1.r2 = lmer_alt(rating_justif ~ obs.justif.mean + (1 | subject) +
           (1 | story_ind),
         data = df.imagharm)
summary(rePCA(m.ratings1.r2))
summary(m.ratings1.r2)
m.ratings1.r2.std = lmer_alt(rating_justif ~ obs.justif.mean + (1 | subject) +
                           (1 | story_ind),
                         data = df.imagharm %>% mutate(rating_justif = scale(rating_justif),
                                                       obs.justif.mean = scale(obs.justif.mean)))
print.stats(m.ratings1.r2, m.ratings1.r2.std, 'obs.justif.mean')

m.ratings2.r3 = lmer_alt(orig.justif ~ rating2 +
                           (rating2 | orig.story) + (rating2 | orig.subject),
                         data = df.justif.all)
summary(rePCA(m.ratings2.r3))
summary(m.ratings2.r3)

ggplot(df.imagharm, aes(x = obs.justif.bin)) +
  geom_histogram(stat='count') +
  scale_x_discrete(labels = c('1-2', '2-3', '3-4', '4-5', '5-6', '6-7')) +
  scale_y_discrete(breaks = NULL) +
  #geom_vline(xintercept = 3.5, color = 'black', size = 1.4) +
  labs(x = 'Observer-rated justification', y = 'Number of trials')
m = mean(df.imagharm$obs.justif.mean > 4)
s = se.prop(df.imagharm$obs.justif.mean > 4)
c(m-1.96*s,m,m+1.96*s)

m = mean(df.imagharm$obs.justif.mean)
s = se(df.imagharm$obs.justif.mean)
c(m-1.96*s,m,m+1.96*s)
t.test(df.imagharm$obs.justif.mean - 4)

m.ratings3 = lmer_alt(obs.justif.pctyes ~ 1 + (1 | subject) + (1 | story_ind),
                      data = df.imagharm %>% mutate(obs.justif.pctyes = obs.justif.pctyes - .5))
summary(m.ratings3)

m.ratings3 = lmer_alt(rating2 ~ 1 + (1 | subject) + (1 | orig.story) + (1 | orig.subject),
                      data = df.justif.all %>% mutate(rating2 = rating2 - 4))
summary(m.ratings3)

ggplot(df.imagharm, aes(x = obs.justif.mean, y = rating_will)) +
  geom_point() +
  geom_smooth(method='lm')

for (i in 1:nrow(df.filt)) {
  imagharm.row = which(df.imagharm$subject == df.filt$subject[i] & df.imagharm$story_id == df.filt$story_id[i])
  if (length(imagharm.row) > 0) {
    df.filt$obs.justif[i] = df.imagharm$obs.justif.mean[imagharm.row]
    df.filt$obs.pctfail[i] = df.imagharm$obs.justif.pctfail[imagharm.row]
  } else {
    df.filt$obs.justif[i] = df.filt$rating_justif[i]
  }
}

df.filt = df.filt %>% mutate(obs.failed = factor(obs.pctfail > .5, c(F,T), c(F,T)))

df.bytrialtypeandsubject.obs = df.filt %>%
  group_by(imagine, subject) %>%
  summarize(j2 = mean(rating_justif), j2.obs = mean(obs.justif), will = mean(rating_will)) %>%
  arrange(subject) %>%
  group_by(subject) %>%
  mutate(will.diff = will - lag(will), j2.diff = j2 - lag(j2)) %>%
  filter(imagine == 'Imagine')
ggplot(df.bytrialtypeandsubject.obs, aes(x = j2.obs, y = will.diff)) +
  geom_point() +
  geom_smooth(method = 'lm', color = 'black') +
  scale_x_continuous(limits = c(3,7.1), breaks = c(3,7)) +
  scale_y_continuous(limits = c(-3.75, 4), breaks = c(-3,0,3)) +
  labs(x = 'Average \nobserver-rated justification', y = 'Average effect of\nimagination on likelihood')

m.will.subj.obs = lm(will.diff ~ j2.obs, data = df.bytrialtypeandsubject.obs)
summary(m.will.subj.obs)

df.bytrialtypeandjustif = df.filt %>%
  group_by(imagine, obs.justif.high, subject) %>%
  summarize(will = mean(rating_will)) %>%
  group_by(imagine, obs.justif.high) %>%
  summarize(will.m = mean(will), will.se = se(will))

ggplot(df.bytrialtypeandjustif, aes(x = obs.justif.high, y = will.m, group = imagine, fill = imagine)) + 
  geom_col(position = dodge) +
  geom_errorbar(aes(ymin = will.m - will.se, ymax = will.m + will.se), width = .1, position = dodge) +
  labs(x = "Relative justification\n(split by median)", y = "Reported likelihood of\nperforming behavior") +
  scale_y_continuous(limits = c(0,5), breaks = c(0,5)) +
  guides(fill = guide_legend(title = "Trial Type")) +
  scale_fill_brewer() +
  theme(legend.position = 'bottom')

# save --------------------------------------------------------------------

save.image(paste0('Study', expt, '_output.rdata'))

# make stuff for follow-up

# stuff for making experiment - DELETE BEFORE POSTING
to.print <- '["'
responses <- df.imagharm$response

response <- as.character(responses[1])
to.print <- paste0(to.print, substr(response, gregexpr('"', response)[[1]][1], nchar(response)))
for (i in 2:length(responses)) {
  response <- as.character(responses[i])
  to.print <- paste(to.print, paste0('"', substr(response, gregexpr('"', response)[[1]][1], nchar(response))), sep = '",\n')
}
to.print <- paste0(to.print, '"];')
cat(to.print)

sample(responses, 25)

paste(as.character(df.imagharm$story_id), sep="' '", collapse=", ")
