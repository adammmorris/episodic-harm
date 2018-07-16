# setup -------------------------------------------------------------------

# only works in Rstudio -- otherwise you have to set the path manually!
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

require(ggplot2)
require(lme4)
require(lmerTest)
require(dplyr)
require(RColorBrewer)

theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black"),
             axis.text=element_text(size=20, colour = "black"), axis.title=element_text(size=18, face = "bold"), axis.title.x = element_text(vjust = 0),
             legend.title = element_text(size = 24, face = "bold"), legend.text = element_text(size = 20), plot.title = element_text(size = 26, face = "bold", vjust = 1),
             axis.title.y = element_text(vjust = 1))
ymin = 2
ymax = 6

se <- function(data) {return(sd(data)/sqrt(length(data)))}
dodge <- position_dodge(width=0.9)
getExcludeSubj = function(data, minRespLength, nResponses) {
  return((data %>% group_by(subject) %>%
            summarise(nResponses = length(response), respLength = mean(nchar(encodeString(as.character(na.omit(response)))))) %>%
            filter(respLength < minRespLength | nResponses != nResponses))$subject)
}


# import data -------------------------------------------------------------

data <- read.csv("data.csv") %>% arrange(subject) %>%
  mutate(imagine = factor(imagine, labels = c("Control", "Imagine"), levels = c(0,1)),
         condition = factor(condition, labels = c("Harm","Help"), levels = c(0,1)),
         rationalization = factor(rationalization, levels = c(0,1), labels = c("No", "Yes")),
         story_ind = factor(story, labels = 1:14))
data = data %>% filter(!(subject %in% getExcludeSubj(data, 50, length(unique(story_id)))))

df.bysubj <- data %>% group_by(condition, imagine, subject) %>%
  summarise(rating_will = mean(rating))

df.bycond <- df.bysubj %>% group_by(condition, imagine) %>%
  summarise(will = mean(rating_will), will.se = se(rating_will))


# main test ---------------------------------------------------------------
# how does imagining (vs control) alter the likelihood of performing behavior in each condition?

## plot result
ggplot(df.bycond, aes(x = imagine, y = will, colour = condition, group = condition)) +
  geom_line(aes(), size = 1) +
  geom_point(aes(), size = 5) +
  geom_errorbar(aes(ymax = will + will.se, ymin = will - will.se), width = .1) +
  labs(x = "", y = "Likelihood of performing behavior") +
  guides(linetype = guide_legend(title = "Condition"),
         colour = guide_legend(title = "Condition"),
         shape = guide_legend(title = "Condition")) +
  ylim(ymin, ymax)

## fit linear models
# first, test for interaction
model <- lmer(rating ~ condition * imagine + (1 + imagine | subject) + (1 + condition * imagine | story_ind), data = data, REML = F)
model.noint <- lmer(rating ~ condition + imagine + (1 + imagine | subject) + (1 + condition * imagine | story_ind), data = data, REML = F)
anova(model, model.noint)

# then, test for main effets
model.onlyimag <- lmer(rating ~ imagine + (1 + imagine | subject) + (1 + condition * imagine | story_ind), data = data, REML = F)
model.onlycond <- lmer(rating ~ condition + (1 + imagine | subject) + (1 + condition * imagine | story_ind), data = data, REML = F)
anova(model.noint, model.onlycond) # main effect of imagine (chisq(1) = 6.49, p = .01)
anova(model.noint, model.onlyimag) # main effect of condition (chisq(1) = 10.78, p = .001)

# finally, test for simple effects
model.help <- lmer(rating ~ imagine + (1 + imagine | subject) + (1 + imagine | story_ind), data = data %>% filter(condition == 'Help'), REML = F)
model.help.null <- lmer(rating ~ 1 + (1 + imagine | subject) + (1 + imagine | story_ind), data = data %>% filter(condition == 'Help'), REML = F)
anova(model.help, model.help.null)
model.harm <- lmer(rating ~ imagine + (1 + imagine | subject) + (1 + imagine | story_ind), data = data %>% filter(condition == 'Harm'), REML = F)
model.harm.null <- lmer(rating ~ 1 + (1 + imagine | subject) + (1 + imagine | story_ind), data = data %>% filter(condition == 'Harm'), REML = F)
anova(model.harm, model.harm.null)


# justifiability ----------------------------------------------------------
df.imagharm = data %>% filter(imagine == 'Imagine' & condition == 'Harm')

# stuff for making experiment
to.print <- '[Q"Q'
responses <- df.imagharm$response

response <- as.character(responses[1])
to.print <- paste0(to.print, substr(response, gregexpr('"', response)[[1]][1], nchar(response)))
for (i in 2:length(responses)) {
  response <- as.character(responses[i])
  to.print <- paste(to.print, paste0('Q"Q', substr(response, gregexpr('"', response)[[1]][1], nchar(response))), sep = 'Q"Q, ')
}
to.print <- paste0(to.print, 'Q"Q];')
cat(to.print)

paste(as.character(df.rating$story_ind), sep="' '", collapse=", ")

# get ratings
df.justif = read.csv('justif_ratings.csv') %>% group_by(resp_num) %>%
  summarize(rating1.pctyes = mean(rating1 == 0), rating1.pctno = mean(rating1 == 1), rating1.pctfail = mean(rating1 == 2),
            rating2.mean = mean(rating2), rating2.se = se(rating2)) %>%
  arrange(resp_num) %>%
  mutate(rating2.mean = rating2.mean + 1, rating.high = rating2.mean > 4)

# add to df.imagharm
df.imagharm$justif.pctyes = df.justif$rating1.pctyes
df.imagharm$justif.pctno = df.justif$rating1.pctno
df.imagharm$justif.pctfail = df.justif$rating1.pctfail
df.imagharm$justif.mean = df.justif$rating2.mean
df.imagharm$justif.se = df.justif$rating2.se
df.imagharm$justif = df.justif$rating.high

# check out data
ggplot(df.imagharm, aes(df.imagharm$justif.mean)) + geom_histogram() + xlim(1,7) +
  labs(x = "How morally justified (1-7)", y = "Count") +
  geom_vline(xintercept = 4, color = 'red', linetype = 2)
mean(df.imagharm$justif.mean > 4) # 78% over midpoint
t.test(df.imagharm$justif.mean - 4)
ggplot(df.imagharm, aes(df.imagharm$justif.pctyes)) + geom_histogram() + xlim(0,1)
mean(df.imagharm$justif.pctyes > .5) # 75% had >50% of ppl say yes

# test for relationship w/ likelihood ratings
ggplot(df.imagharm, aes(x = df.imagharm$justif.mean, y = df.imagharm$rating)) +
  geom_point(size = 2, stroke = 2) +
  geom_smooth(method='lm') +
  labs(x = "How morally justified (1-7)", y = "Likelihood of performing behavior")
model.justif1 = lmer(rating ~ justif.mean + (1 + justif.mean | subject) + (1 + justif.mean | story_ind), data = df.imagharm)
model.justif1.null = lmer(rating ~ 1 + (1 + justif.mean | subject) + (1 + justif.mean | story_ind), data = df.imagharm)
summary(model.justif1)
anova(model.justif1, model.justif1.null)

ggplot(df.imagharm, aes(x = df.imagharm$justif.pctyes, y = df.imagharm$rating)) +
  geom_point(size = 2, stroke = 2) +
  geom_smooth(method='lm') +
  xlim(0,1)
model.justif2 = lmer(rating ~ justif.pctyes + (1 | subject) + (0 + justif.pctyes | subject) + (1 + justif.pctyes | story_ind), data = df.imagharm)
summary(model.justif2)

# put together w/ control data
df.imagharm.collapsed = df.imagharm %>% group_by(justif, subject) %>% summarize(rating = mean(rating)) %>%
  group_by(justif) %>% summarize(rating.mean = mean(rating), rating.se = se(rating))
  
df.combined = data.frame(category = factor(), rating = numeric(), rating.se = numeric())
df.combined = rbind(df.combined, data.frame(category = 'Control',
                             rating = df.bycond$will[df.bycond$condition == 'Harm' & df.bycond$imagine == 'Control'],
                             rating.se = df.bycond$will.se[df.bycond$condition == 'Harm' & df.bycond$imagine == 'Control']))
df.combined = rbind(df.combined, data.frame(category = 'Imagine--unjustified',
                                      rating = df.imagharm.collapsed$rating.mean[df.imagharm.collapsed$justif == F],
                                      rating.se = df.imagharm.collapsed$rating.se[df.imagharm.collapsed$justif == F]))
df.combined = rbind(df.combined, data.frame(category = 'Imagine--justified',
                                      rating = df.imagharm.collapsed$rating.mean[df.imagharm.collapsed$justif == T],
                                      rating.se = df.imagharm.collapsed$rating.se[df.imagharm.collapsed$justif == T]))

ggplot(df.combined, aes(x = category, y = rating)) +
  geom_point(size = 2, stroke = 2) + 
  ylim(2,6) + #scale_y_continuous(breaks = c(2,4,6)) +
  geom_errorbar(width = .1, aes(ymin = rating - rating.se, ymax = rating + rating.se)) +
  labs(x = "Trial type", y = "Likelihood of performing behavior")

data$justif = numeric(nrow(data))
for (i in 1:nrow(data)) {
  if (data$condition[i] == 'Harm' & data$imagine[i] == 'Imagine') {
    data$justif[i] = as.numeric(df.imagharm$justif[df.imagharm$id == data$id[i]])
  } else {
    data$justif[i] = -1
  }
}

model.harm.justif <- lmer(rating ~ imagine + (1 + imagine | subject) + (1 + imagine | story_ind), data = data %>% filter(condition == 'Harm' & justif != 0), REML = F)
model.harm.justif.null <- lmer(rating ~ 1 + (1 + imagine | subject) + (1 + imagine | story_ind), data = data %>% filter(condition == 'Harm' & justif != 0), REML = F)
anova(model.harm.justif, model.harm.justif.null)

model.harm.unjustif <- lmer(rating ~ imagine + (1 + imagine | subject) + (1 + imagine | story_ind), data = data %>% filter(condition == 'Harm' & justif != 1), REML = F)
model.harm.unjustif.null <- lmer(rating ~ 1 + (1 + imagine | subject) + (1 + imagine | story_ind), data = data %>% filter(condition == 'Harm' & justif != 1), REML = F)
anova(model.harm.unjustif, model.harm.unjustif.null)