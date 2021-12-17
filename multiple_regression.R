## Create data
mouse.data <- data.frame(
  size = c(1.4, 2.6, 1.0, 3.7, 5.5, 3.2, 3.0, 4.9, 6.3),
  weight = c(0.9, 1.8, 2.4, 3.5, 3.9, 4.4, 5.1, 5.6, 6.3),
  tail = c(0.7, 1.3, 0.7, 2.0, 3.6, 3.0, 2.9, 3.9, 4.0))

## Look at the data
mouse.data

## Focus on how well weight predicts size
## Draw a graph of the data to make sure the relationship makes sense
plot(mouse.data$weight, mouse.data$size, pch=16, cex=2)

## Do the regression (size = y-intercept+slope*weight)
simple.regression <- lm(size ~ weight, data=mouse.data)

## Look at the R^2, F-value and p-value
summary(simple.regression)
# According to R^2 and p-value, weight is an important factor

#Add least-squares fit line
abline(simple.regression, lwd=5, col="red")

## Let's do multiple regression by adding an extra term, tail length
## Draw a graph of the data to make sure the relationship make sense
plot(mouse.data)
## This graph is more complex because it shows the relationships between all
## of the columns

## Do the regression (size=y-intercept+ slope_1*weight + slope_2*tail )
multiple.regression <- lm(size ~ weight + tail, data=mouse.data)

## Look at the R^2, F-value and p-value
summary(multiple.regression)
## As seen, p-vale and R^2 is adjusted which indicated  a better prediction.
# And also coefficients tells that, using weight and tail is not better than using tail alone but,
# is better than using weight alone.
## So using tail values are more reasonable.