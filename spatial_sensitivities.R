### Plotting the spatial sensitivities of a gravimeter, a gradiometer and a synthetic gradiometer (2 gravimeters at vertical distance L)

require(grDevices) # for colours
require(viridis)

z <- seq(-5,5,0.01) # 
x <- seq(-5,5,0.01) # horizontal extend
y <- seq(-5,5,0.01) # horizontal extend

G=6.67e-11 #gravitational constant
m1= 1 #point mass
L1 = 1.2 # distance between the gravimeters

lvls <- seq(-40,-40,-1)
lvls2 <- seq(-20,-40,-1)


Tzz <- matrix(nrow= length(x), ncol= length(z),0)
TzzL1 <- matrix(nrow= length(x), ncol= length(z),0)
#TzzL2 <- matrix(nrow= length(x), ncol= length(z),0)
gzz <- matrix(nrow= length(x), ncol= length(z),0)

for (i in 1:length(z)){ 
  
  for (j in 1:length(x)){ 
    
    Tzz[j,i] = G*m1*(abs(x[j])^2-2*abs(z[i])^2)/(abs(x[j])^2+abs(z[i])^2)^(5/2)# gradiometer
    #2 gravimeters with distance L apart from each other 
    TzzL1[j,i]=G*m1*((abs(z[i])+L1/2)/((abs(z[i])+L1/2)^2+abs(x[j])^2)^(3/2)-(abs(z[i])-L1/2)/((abs(z[i])-L1/2)^2+abs(x[j])^2)^(3/2))/L1 
    #TzzL2[j,i]=G*m1*((abs(z[i])+L2/2)/((abs(z[i])+L2/2)^2+abs(x[j])^2)^(3/2)-(abs(z[i])-L2/2)/((abs(z[i])-L2/2)^2+abs(x[j])^2)^(3/2))/L2 
    gzz[j,i]=G*m1*abs(z[i])/(abs(z[i])^2+abs(x[j])^2)^(3/2) #gravimeter
  }
}

# par(mfrow=c(1,3))
lvls <- seq(-40,-20,1)
lvls2 <- lvls
par(mar=c(4,5,1,1), xpd=F)
p1 <- filled.contour(x,z,(log(abs(gzz*10))),levels = lvls,key.axes={axis(4, cex.axis=3)}, plot.axes = { contour(x,z,log(abs(gzz)), col=c("brown"),
                                                                           levels = lvls2, drawlabels = F, axes = FALSE, frame.plot = FALSE, add = TRUE); axis(1, cex.axis=3); axis(2, cex.axis=3) }, color = viridis,plot.title =title(main = " ",sub ="",xlab=c("x"), ylab=c(" z"),cex.lab=3))



p2 <- filled.contour(x,z,(log(abs(Tzz))),levels = lvls, color = viridis,plot.axes = { contour(x,z,log(abs(Tzz)), col=c("brown"), 
                                                                                                            levels = lvls2, drawlabels = F, axes = FALSE,
                                                                                                            frame.plot = FALSE, add = TRUE); axis(1,cex.axis=3); 
                                                                                        axis(2,cex.axis=3)}, plot.title = title(main="", sub="",xlab=c("x"), cex.lab=3,ylab=c(" z")))

p3 <- filled.contour(x,z,(log(abs(TzzL1))),levels = lvls,key.axes={axis(4, cex.axis=3)},plot.axes = { contour(x,z,log(abs(TzzL1)), col=c("brown"), levels = lvls2, 
                                                                           drawlabels = F, axes = FALSE, frame.plot = FALSE, add = TRUE); 
  axis(1, cex.axis=3, cex.lab=3); axis(2, cex.axis=3) },color = viridis,plot.title = title(main = "",sub="",xlab=c("x"),  cex.lab=3,ylab=c(" z")))



grid.arrange(p1, p2, p3, ncol=3)
