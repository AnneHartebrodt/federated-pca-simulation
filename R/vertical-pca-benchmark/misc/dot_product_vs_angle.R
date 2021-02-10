require(REdaS)
require(ggplot2)
require(data.table)
options(digits = 12)

angle<- function(co){
  
theta <- acos(co)
a <- rad2deg(theta)
# angle can be Nan
# print(angle)
# return the canonical angle
if (a > 90)
 {
   return (np.abs(a - 180))
 } else{
  return (a)
}
}

# coarse grained plot
steps <- seq(0.0, 1.0, 0.05)
a<-sapply(steps, function(x) angle(x))
data<- as.data.table(list(a=a, b=steps))
ggplot(data, aes(steps, a))+geom_line()+ylab('Angle[degree]')+xlab('Dot product')


# fine grained plot
steps <- seq(0.9999, 1.0, 1e-8)

a<-sapply(steps, function(x) angle(x))
data<- as.data.table(list(a=a, b=steps))
ggplot(data, aes(steps, a))+geom_line()+ylab('Angle[degree]')+xlab('Dot product')+
  scale_x_continuous(n.breaks = 15)+theme(axis.text = element_text(angle=90))



# fine grained plot
steps <- seq(0.9999999, 1.0, 1e-10)

a<-sapply(steps, function(x) angle(x))
data<- as.data.table(list(a=a, b=steps))
ggplot(data, aes(steps, a))+geom_line()+ylab('Angle[degree]')+xlab('Dot product')+
  scale_x_continuous(n.breaks = 15)+theme(axis.text = element_text(angle=90))


