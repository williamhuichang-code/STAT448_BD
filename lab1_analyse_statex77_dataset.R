


##### Load the dataset
states <- as.data.frame(state.x77)



##### check cols and rows
# check for cols
length(states)
if (length(states) != 0) {
  print("Not empty")
}
# check for rows
nrow(states)



##### exploring the dataset
# dimensions
dim(states)
# structure
str(states)



##### histograms
hist(
  states$Income,
  col = "orange",
  main = "Histogram of Income",
  xlab = "Income"
)



##### boxplot
boxplot(
  states$Income,
  col = "orange",
  main = "Boxplot of Income",
  ylab = "Income (per capita)"
)


##### scatter plot
plot(
  states$Income,
  states$Illiteracy,
  pch = 16,
  col = "red",
  main = "Income vs. Illiteracy",
  xlab = "Illiteracy Rate (%)",
  ylab = "Income ($)"
)



##### compute correlation coefficient
cor(states$Income, states$Illiteracy)



##### build a linear model
model <- lm(Income ~ Illiteracy, data = states)



##### model summary
summary(model)


##### Assumptions of Linear Regression
# linearity
#   the relationship between the independent and dependent variables is linear
# independence
#   observations are independent
# homoscedasticity
#   the variance of residual is the same
#       for any value of the independent variable
# normality
#   the residuals of the model are normally distributed
# no multicollinearity
#   independent variables are not highly correlated with each other



##### plotting the regression line
# puts the columns of states into the search path
attach(states)
# train the model
model <- lm(Income ~ Illiteracy, data = states)
# plot the actual data
plot(Illiteracy, Income, pch = 16, col = "red")
# add the regression line
abline(lm(Income ~ Illiteracy), col = "blue")
# add predicted values as points
lines(Illiteracy, lm(Income ~ Illiteracy)$fitted.values, col="green", type='b')
# diagnostic tool: add a lowess smoother
lines(lowess(Illiteracy, Income), col="orange")



##### plotting the residuals
plot(model$fitted.values, model$residuals, pch = 16, col = 'red',
     main="Residuals vs Fitted Values",
     xlab="Fitted Values",
     ylab="Residuals")
abline(h=0, col="blue")



##### QQ plot of residuals
# linear regression assumes
#   ε∼N(0,σ2)
# slope β formula works purely from geometry (matrix manipulations)
#   no normality needed
# why we care about normal residuals
#   hypothesis testing for coefficients
#   this derived slope is statistically significant 0.95
# the theoretical quantile
#   derived from inverse normal CDF
qqnorm(model$residuals)
qqline(model$residuals, col = "blue")



##### histogram of the model residuals
hist(model$residuals)


##### random string here