##Graphics=group
##Bar plot=name
##showplots
##Layer=vector
##Field=Field Layer

# Keep only the 10 first characters of the field name
field10chars <- substr(Field, start = 1, stop = 10)

# Generate bar plot
lattice::barchart(table(Layer[[field10chars]]), 
                  main = paste("Bar plot of", Field),
                  horizontal = FALSE,
                  col = 'gray',
                  xlab = paste(Field),
                  ylab = "Frequency", 
                  panel = function (x, y, ...) {
                    lattice::panel.grid(v = 0, h = -1, lty = 'dotted')
                    lattice::panel.barchart(x, y, ...)
                  })
