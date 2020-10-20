
# Influence of operator's body mass on vertical gravity gradient estimation (based on gravity measurements on two heights)
as a function of horizontal distance


## Gravity forward-modeling : prism code (after Okabe, 1979)
prism_grav_fun <- function(P,xc,yc,zc,delta_rho,G=6.67259e-11) {
  
  
  dx <- matrix(nrow =1, ncol = 2)
  dx[,1] <- P[1]-xc[1]
  dx[,2] <- P[1]-xc[2]
  dx
  
  dy <- matrix(nrow =1, ncol = 2)
  dy[,1] <- P[2]-yc[1]
  dy[,2] <- P[2]-yc[2]
  dy
  
  dz <- matrix(nrow =1, ncol = 2)
  dz[,1] <- P[3]-zc[1]
  dz[,2] <- P[3]-zc[2]
  dz
  
  delta_grav <- c(0)
  
  eps =  1e-6
  
  for (i in 1:2){
    for (j in 1:2){
      for (k in 1:2){
        
        r <- sqrt(dx[,i]^2+dy[,j]^2+dz[,k]^2)
        
        
        delta_grav <- delta_grav + (-G*delta_rho)*((-1)^i)*((-1)^j)*((-1)^k)*(dx[,i]*log(dy[,j]+r+eps)+dy[,j]*log(dx[,i]+r+eps)+2*dz[,k]*atan((dx[,i]+dy[,j]+r+eps)/dz[,k]))
      }
    }
  }
  delta_grav <-  delta_grav*10e8 # microGal
  return(delta_grav)
}
#  Test of a prism of my "size"

L1 = 0.261
L2 = 0.5925 + L1#m
L3 =  1.175 + L1#m

# distance between the heights
D1 <- 0.5925# m
D2 <- 0.5825# m
G=6.67e-11 #gravitational constant

# Influence of body on vertical gravity gradient estimation based on two gravity measurements

height <- 1.64 # m
distance <- seq(0,1,0.01)
len <- 0.20 # m
wid <- 0.20 # m
rho_human = 985 #kg/m3
#Point of observation pillar 3
P <- c(3.85,1.5,-L1)
GradP3H02 <- rep(0,length(distance))

for (i in 1:length(distance)){ 

xc<- c(P[1]+sqrt(distance[i]/2),P[1]+sqrt(distance[i]/2)+len)
yc<- c(P[2]+sqrt(distance[i]/2),P[2]+sqrt(distance[i]/2)+wid)
zc <- c(0,-height) #c(1,1.7)
  
P[3] <- -L1
P3H0_grav_me <- - prism_grav_fun(P,xc,yc,zc,delta_rho=rho_human,G=6.67259e-11)   
P[3] <- -L2
P3H1_grav_me <- - prism_grav_fun(P,xc,yc,zc,delta_rho=rho_human,G=6.67259e-11)   

P[3] <- -L3
P3H2_grav_me<- - prism_grav_fun(P,xc,yc,zc,delta_rho=rho_human,G=6.67259e-11)   


GradP3H02[i] <- -((P3H0_grav_me-P3H2_grav_me)/((D1+D2)*100))*1000

}

## Plotting

par(mar=c(6,6,1,1))
plot(distance,GradP3H02, cex.axis=2.5,type="l",lwd=2.5, ylab="", xlab="") #main="Influence on vertical gradient  \n (estimated from three measurement on pillar 3) \n of a box-shaped human (1.64 m 65 kg) \n as a function of distance", cex.main=0.8
grid(col="grey30", lty=3)
mtext(text = "VGG [E]", side=2,line=4, cex=2.5)
mtext(text = "distance [m]",side=1, line=4, cex=2.5)
