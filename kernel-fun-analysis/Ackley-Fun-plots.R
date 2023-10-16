library(ggplot2)
library(plot3D)

## Ackley
obj.fn = makeAckleyFunction(dimensions = 2L)
print(obj.fn)
print(autoplot(obj.fn))
plot3D(obj.fn, length.out = 50L, contour = TRUE)

print(autoplot(makeAckleyFunction(1)))

## Branin
obj.fn = makeBraninFunction()
plot3D(obj.fn, length.out = 50L, contour = TRUE)

print(autoplot(makeBraninFunction()))
