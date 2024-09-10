library(PLNmodels)
data(oaks)

myPLN <- PLN(Abundance ~ 0 + tree + offset(log(Offset)), data = oaks)
myZIPLN <- ZIPLN(Abundance ~ 0 + tree + offset(log(Offset)) | 0 + tree, data = oaks)

## bricolage car quelques points explosent, sans conséquences sur la convergence finale
obj_PLN <- -myPLN$optim_par$objective
plot(obj_PLN, ## optimisation nlopt de la vraisemblance profilée
     ylim = c(obj_PLN[1], # fonction objective sans les constantes
              tail(obj_PLN,1)), log = "xy")

plot(myZIPLN$optim_par$objective, type = "l", log = "xy") ## Les itérations du VEM: la vraisemblance estimée
