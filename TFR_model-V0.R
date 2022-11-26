####################### TFR MOdel 

a=b=1

nu=1   # N(t) ~ poisson(nu)

mu=2; theta = 0.5   # W~Gamma(mu,theta)


gamma = 0.2 ; eta = 1.5

lambda_0 = 2 ; gamma_0=0.5

# % Reliability function

##### HPP

R=function(t, a, b, nu, gamma, eta, theta, mu){
  
  F=function(s,t){
    (1/(1-a*theta*((s/eta)^gamma-(t/eta)^gamma)))^mu * 
      I(a*((s/eta)^gamma-(t/eta)^gamma) < 1/theta)
  }
  
  H=function(t){sapply(t, function(t) {
    integrate(F, lower = 0, upper = t, t=t)$value
  })}
  
  exp(-b*(t/eta)^gamma - nu*t + nu * H(t))
  
}

R(seq(0,1,0.01), a, b, nu, gamma, eta, theta, mu)
R(Inf, a, b, nu, gamma, eta, theta, mu)

Rt = function(t){R(t, a, b, nu, gamma, eta, theta, mu)}



plot(Rt,0,1)

plot(Rt,0,3)

######NHPP 


R=function(t, a, b, nu, gamma, eta, theta, mu){
  
  nu = function(t) lambda_0 + gamma_0 * t
  NU = function(t) lambda_0 * t + gamma_0 * t^2/2
  
  F=function(s,t){
    nu(s)*(1/(1-a*theta*((s/eta)^gamma-(t/eta)^gamma)))^mu *
      I(a*((s/eta)^gamma-(t/eta)^gamma) < 1/theta)
  }
  
  H=function(t){sapply(t, function(t) {
    integrate(F, lower = 0, upper = t, t=t)$value
  })}
  
  exp(-b*(t/eta)^gamma - NU(t) + H(t))
  
}

R(seq(0,1,0.01), a, b, nu, gamma, eta, theta, mu)
R(Inf, a, b, nu, gamma, eta, theta, mu)

Rt = function(t){R(t, a, b, nu, gamma, eta, theta, mu)}

plot(Rt,0,0.1)

plot(Rt,0,3)


############################### Maintenance Section 

###  %%%%%%   Condition Based Maintenance _ considering NHPP

k=1
tau=2


#### Preventive Probability 


Pp = function(n,m,k) {
  
  nu = function(t) lambda_0 + gamma_0 * t
  
  NU = function(t) lambda_0 * t + gamma_0 * t^2/2
  
  Gam = function(t) (t/eta)^gamma
  
  F0 = function(t) {
    nu(t)*(1/(1-a*theta*(Gam(t)-Gam((k+1)*tau))))^mu *
      I(a*(Gam(t)-Gam((k+1)*tau)) < 1/theta)
  }
  F0_int = integrate(F0, k*tau, (k+1)*tau)$value
  
  p=(F0_int^m/factorial(m)) * exp(-NU((k+1)*tau)) * ((NU(k*tau))^n)/factorial(n) *
    exp(b*(Gam(k*tau)-Gam((k+1)*tau))) * (1/(1-a*theta*(Gam(k*tau)-Gam((k+1)*tau))))^(n*mu) *
    I(a*(Gam(k*tau)-Gam((k+1)*tau)) < 1/theta) * 
    Rt(k*tau)
  return(p)
}
Pp(1,1,3)


pp=0
M=200
ppnstar = function(nstar,k){
  
  for(m in 1:M){
    for(n in 0:nstar){
      pp=pp+Pp(n,m,k)
    }
  }
  return(pp)
}


ppnstar(5,3)

nstar = 2

ppk <- function(k) {
  ppnstar(nstar = nstar , k)
}

ppk(1)
s1 = seq(0,10,0.1)
s2 = c()
for (i in 1:length(s1)) {
  s2[i] = ppk(s1[i])
}

plot(s1,s2,type = "l")

#### Corrective Probability 


# type_1

Pc= function(k){
  Rt(k*tau) - Rt((k+1)*tau)
}

Pc(0)

# type_2

Pc_0 <- function(k) {
  
  nu = function(t) lambda_0 + gamma_0 * t
  
  NU = function(t) lambda_0 * t + gamma_0 * t^2/2
  
  Gam = function(t) (t/eta)^gamma
  
  Mgf = function(t , theta, mu) (1/(1-theta*t))^mu
  
  F1 = function(t) {
    nu(t) * Mgf(Gam(t)-Gam((k+1)*tau), theta * a , mu) * 
      I(a*(Gam(t)-Gam((k+1)*tau))<1/theta)
  }
  
  pc0 = exp(integrate(F1 , k*tau , (k+1)*tau)$value -
              NU((k+1)*tau) + b *(Gam(k * tau) - Gam((k+1) * tau)) + 
              (Mgf(Gam(k * tau)-Gam((k+1)*tau), theta * a , mu) *
                 I(a*(Gam(k * tau)-Gam((k+1)*tau))<1/theta) *
                 NU(k*tau)))
  
  
  return((1-pc0) * Rt(k * tau))
  
}

Pc(0.2)
Pc_0(0.2)

plot(seq(0,1,0.1),Pc(seq(0,1,0.1)),type = "l")

s1 = seq(0,1,0.1)
s2 = c()
for (i in 1:length(s1)) {
  s2[i] = Pc_0(s1[i])
}
plot(s1,s2,type = "l")

############################# Lung rate Cost Function
K = 10
c_ins = 1 ; c_p = 2 ; c_c = 3

ETr = function(nstar){
  etr = 0
  for (k in 1:K) {
    etr = etr + k*tau * (ppnstar(nstar = nstar , k-1) + Pc(k-1))
  }
  return(etr)
} 

ETr(2)


ECr = function(nstar){
  ecr = 0
  for (k in 1:K) {
    ecr = ecr + ((k*c_ins + c_p)*ppnstar(nstar = nstar , k-1) + 
                   (k*c_ins + c_c)*Pc(k-1))
  }
  return(ecr)
} 
ECr(2)


L <- function(nstar) {
  l = ECr(nstar) / ETr(nstar)
  return(l)
}

tau = 2

nsta = seq(1,10,1)
LRCR = sapply(nsta , L)
LRCR

 which.min(LRCR)

plot(nsta, LRCR, type = "l")

plot(seq(tau,K*tau,tau), LRCR, type = "l")

library("plot3D")

scatter3D(x = nsta, y = seq(tau,K*tau,tau), z = LRCR,
          phi = 0, bty = "g",  type = "h", 
          ticktype = "detailed",
          pch = 19, cex = 0.5, xlab ="n*", ylab = "t", zlab = "L", zlim = NULL)



################################
###  %%%%%%   Imperfect Maintenance _ considering NHPP   
q=0.2


Rq=function(t, a, b, nu, gamma, eta, theta, mu){
  
  nu = function(t) lambda_0 + gamma_0 * t
  NU = function(t) lambda_0 * t + gamma_0 * t^2/2
  
  F=function(s,t){
    nu(s)*(1/(1-q*a*theta*((s/eta)^gamma-(t/eta)^gamma)))^mu * 
      I(q*a*((s/eta)^gamma-(t/eta)^gamma)<1/theta)
  }
  
  H=function(t){sapply(t, function(t) {
    integrate(F, lower = 0, upper = t, t=t)$value
  })}
  
  exp(-b*(t/eta)^gamma - NU(t) + H(t))
  
}

Rq(0:10, a, b, nu, gamma, eta, theta, mu)
Rq(Inf, a, b, nu, gamma, eta, theta, mu)

Rtq = function(t){R(t, a, b, nu, gamma, eta, theta, mu)}

plot(Rtq,0,0.1)

plot(Rtq,0,3)


########### Preventive Probability

Ppq = function(n,m,k) {
  
  nu = function(t) lambda_0 + gamma_0 * t
  
  NU = function(t) lambda_0 * t + gamma_0 * t^2/2
  
  Gam = function(t) (t/eta)^gamma
  
  F0 = function(t) {
    nu(t)*(1/(1-a*theta*(Gam(t)-Gam((k+1)*tau))))^mu *
      I(a*(Gam(t)-Gam((k+1)*tau))<1/theta)
  }
  F0_int = integrate(F0, k*tau, (k+1)*tau)$value
  
  p=(F0_int^m/factorial(m)) * exp(-NU((k+1)*tau)) * ((NU(k*tau))^n)/factorial(n) *
    exp(b*(Gam(k*tau)-Gam((k+1)*tau))) * (1/(1-q*a*theta*(Gam(k*tau)-Gam((k+1)*tau))))^(n*mu) *
    I(q*a*(Gam(k*tau)-Gam((k+1)*tau)) < 1/theta) * 
    Rtq(k*tau)
  return(p)
}
Ppq(1,1,3)


pp=0
M=200
ppnstarq = function(nstar,k){
  
  for(m in nstar+1:M){
    for(n in 0:nstar){
      pp=pp+Ppq(n,m,k)
    }
  }
  return(pp)
}


ppnstarq(10,3)

nstar = 5

ppkq <- function(k) {
  ppnstarq(nstar = nstar , k)
}

s1 = seq(0,10,0.1)
s2 = c()
for (i in 1:length(s1)) {
  s2[i] = ppkq(s1[i])
}

plot(s1,s2,type = "l")


#### Corrective Probability 



# type_1

Pcq= function(k){
  Rtq(k*tau) - Rtq((k+1)*tau)
}

Pcq(0)



# type_2

Pc_0q <- function(k) {
  
  nu = function(t) lambda_0 + gamma_0 * t
  
  NU = function(t) lambda_0 * t + gamma_0 * t^2/2
  
  Gam = function(t) (t/eta)^gamma
  
  Mgf = function(t , theta, mu) (1/(1-theta*t))^mu
  
  F1 = function(t) {
    nu(t) * Mgf(Gam(t)-Gam((k+1)*tau), theta * a , mu) *
      I(a*(Gam(t)-Gam((k+1)*tau)) < 1/theta)
  }
  
  pc0 = exp(integrate(F1 , k*tau , (k+1)*tau)$value - NU((k+1)*tau) +
              b *(Gam(k * tau) - Gam((k+1) * tau)) + 
              (Mgf(Gam(k * tau)-Gam((k+1)*tau), theta * a *q , mu) *
                 I(a*q*(Gam(k * tau)-Gam((k+1)*tau)))  *
                 NU(k*tau)))
  
  
  return((1-pc0) * Rtq(k * tau))
  
}

Pcq(10)
Pc_0q(10)

plot(seq(0,10,0.1),Pcq(seq(0,10,0.1)),type = "l")

s1 = seq(0,10,0.1)
s2 = c()
for (i in 1:length(s1)) {
  s2[i] = Pc_0q(s1[i])
}
plot(s1,s2,type = "l")


############################# Lung rate Cost Function
K = 10
c_ins = 1 ; c_p = 2 ; c_c = 3

ETrq = function(nstar){
  etr = 0
  for (k in 1:K) {
    etr = etr + k*tau * (ppnstarq(nstar = nstar , k-1) + Pcq(k-1))
  }
  return(etr)
} 

ETr(2)


ECrq = function(nstar){
  ecr = 0
  for (k in 1:K) {
    ecr = ecr + ((k*c_ins + c_p)*ppnstarq(nstar = nstar , k-1) + 
                   (k*c_ins + c_c)*Pcq(k-1))
  }
  return(ecr)
} 
ECr(2)


Lq <- function(nstar) {
  l = ECrq(nstar) / ETrq(nstar)
  return(l)
}


q= 0.001

nsta = seq(1,10,1)
LRCR = sapply(nsta , Lq)
LRCR

which.min(LRCR)

plot(nsta, LRCR, type = "l")

plot(seq(tau,K*tau,tau), LRCR, type = "l")

library("plot3D")

scatter3D(x = nsta, y = seq(tau,K*tau,tau), z = LRCR,
          phi = 0, bty = "g",  type = "h", 
          ticktype = "detailed",
          pch = 19, cex = 0.5, xlab ="n*", ylab = "t", zlab = "L", zlim = NULL)



