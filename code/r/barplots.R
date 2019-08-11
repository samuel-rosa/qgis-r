##Graphics=group
##Bar plot=name
##showplots
##Layer=vector
##Field=Field Layer

# Keep only the 10 first characters of the field name
field10chars <- substr(Field, start = 1, stop = 10)

# Generate bar plot
barplot(table(Layer[[field10chars]]), main = paste("Bar plot of", Field), xlab = paste(Field), ylab = "Number")
