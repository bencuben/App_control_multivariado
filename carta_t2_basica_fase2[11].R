#Brahian Cano Urrego
#Estudiante de estadística universidad nacional
#Fecha: 17/05/2019
#APLICACION IDENTIFICAR FUENTES DE VARIACION
#PARA CARTA DE CONTROL T2 PARA DATOS NORMALES

#library(MSQC)
library(gtools)

# Definición de la función graficadora del T2 -----------------------------
#Argumentos:
#dat:Matriz,arreglo o data.frame con las observaciones provenientes de la fase 2
#X: El vector de medias teorico o hallado en la fase 1
#S: matriz de covarianzas de las p variables
#n: número de observaciones usadas en el HDS
#alpha: nivel de significancia para el cálculo del UCL

carta.T2.obs.faseI<-function(dat,HSD,alpha=0.05,Ksim=0.8){
  
  # Ploteo de la carta T2 ---------------------------------------------------
  
  #Número de observaciones usadas en la fase 1
  n <- nrow(HSD)
  #Vector de medias historico
  X<-apply(HSD,2,mean)
  #Matriz de varianzas historicas
  S <- var(HSD)
  #Número de observaciones usadas en la fase 2
  n1 <- nrow(dat)
  #Número de variables
  p <- ncol(dat)
  
  #cálculo del UCL teórico a partir de la media y varianza especificada
  UCL<-((p*(n+1)*(n-1)) / (n*(n-p)))*qf((1-alpha),p,n-p)
  
  #cálculo del T2 para cada vector de observaciones
  T2<-mahalanobis(dat,center=X,cov=S)
  
  #Indices 
  Observacion<-1:n1
  
  par(mar = c(4, 5, 3, 5))
  #Plot del T2 vs indice de la observación
  plot(Observacion,T2,las=1,type="l",xlim=c(0,n1*1.1+3),
       ylim=c(0,max(UCL,max(T2))+2),main=expression("Carta"*~T^2),
       ylab=expression(T^2),xlab="No. Observación",font=2)
  
  #Se añade el UCL al gráfico
  abline(h=UCL,lty=3)
  mtext(paste(" UCL=", round(UCL, 2)), side = 4, at = UCL, 
        las = 2,font=1,cex=0.8)
  #Se agregan los puntos al gráfico identificando cuales salen del UCL
  #además de guardar los valores fuera de los limites
  alarma <- NULL
  for(i in 1:n1){
    if(T2[i]>UCL){
      col <- "red"
      alarma <- c(alarma,i)
    }else{
      col <- "black"
    }
    points(Observacion[i],T2[i],pch=20,col=col)
    if(T2[i]>UCL)text(i,T2[i],labels=paste0("obs=",i,collapse = ""),pos=3,font=2,cex=0.7)
  }
  
  
  
  # Metodo MYT --------------------------------------------------------------
  
  print("--------------------------Método MYT-----------------------------")
  for(j in alarma){
    print(paste0("Para la observación ",j))
    
    #se crean dos vectores con el fin de identificar las posibles causas de alarma
    culpables<- NULL
    inocentes<-1:p
    
    #Donde ir guardando los resultados para luego mostrarlos
    T_i<-NULL
    
    #Para los incondicionales
    
    #ucl para los terminos individuales
    ucl <- ( (n +1) / (n) )* qf(1-alpha,1,n -1)
    contador<-1
    #T_i para terminos incondicionales
    for(i in 1:p){
      T_i[i]<-mahalanobis(dat[j,i],center=X[i],cov=S[i,i] )
      if(T_i[i]>ucl){
        culpables[contador] <- inocentes[i]
        contador <- contador+1
      }
    }
    
    #forma linda de mostrar los resultados
    texto<- paste(rep("T_",p),1:p,sep="")
    
    tabla<-cbind(T_i,ucl)
    rownames(tabla) <- texto
    
    #se quitan del procedimientos las variables que ya mostraron alarma
    if(length(culpables)>0){
      inocentes<-inocentes[-culpables]
    }
    if(length(inocentes>0)){
      
      
      #se chequea el subvector
      
      ucl <- ( ((length(inocentes))*(n +1 )*(n -1))
               /  (n*(n-length(inocentes)))  )* qf(1-alpha,length(inocentes),
                                                   n-length(inocentes))
      
      T_i<-mahalanobis(dat[j,inocentes],center=X[inocentes],cov=S[inocentes,inocentes] )
      
      # y se continua con las que siguen siendo inocentes para los casos por pares y etc
      
      if(T_i>=ucl){
        
        #La iteración se empieza en los subconjuntos de tamaño 2
        tamaño <- 2
        culpables_aux<-NULL
        
        while(length(inocentes)>0 ){
          
          T_i<-NULL
          texto<-NULL
          contador<-1
          #Definición de subconjuntos a tomar según el tamaño en el paso
          #auxiliar<-combn(inocentes,tamaño)
          auxiliar<-t(permutations(length(inocentes),tamaño,inocentes))
          #ucl calculado según la cantidad de terminos condicionantes
          ucl <- (((n +1)*(n -1))/(n*(n -(nrow(auxiliar)-1) -1)))*qf(1-alpha,1,
                                                                     n -(nrow(auxiliar)-1) -1)
          #atravez de todas las posibles combinaciones del tamaño especificado, se hallara su respectivo T2
          #a partir de la sustracción
          for(i in 1:ncol(auxiliar)){
            #se calcula del T_i condicional a partir de la igual de las sustraciones
            T_i[i]<-mahalanobis(dat[j,auxiliar[,i]],center=X[auxiliar[,i]],cov=S[auxiliar[,i],auxiliar[,i]] )-
              mahalanobis(dat[j,auxiliar[-1,i]],center=X[auxiliar[-1,i]],cov=S[auxiliar[-1,i],auxiliar[-1,i]] )
            
            #creación del T_i para identificarlo en la tabla
            texto[i] <- paste0(paste0("T_",auxiliar[1,i],"|"),
                               paste0(auxiliar[-1,i],collapse = ","))
            
            if(T_i[i]>ucl){
              culpables_aux<- unique(c(culpables_aux,auxiliar[,i]))
            }
          }
          
          #Se añaden lo culpables encontrados
          culpables <- c(culpables,culpables_aux)
          
          #Se agregan los resultados a la salida
          tabla_aux<-cbind(T_i,ucl)
          rownames(tabla_aux) <- texto
          
          tabla<- rbind(tabla,tabla_aux)
          
          #Si encontro nuevos culpables debemos quitarlos de la lista de inocentes
          if(length(culpables_aux)>0){
            #Se quita de los inocentes aquellos allados como culpables
            inocentes<-inocentes[-match(culpables_aux,inocentes)]
          }
          #En el caso de que nos quedaramos sin inocentes se debe para el ciclo
          if(length(inocentes)>0){#para cuando hayan inocentes
            #Calculo del subvector restante
            ucl <- ( ((length(inocentes))*(n +1 )*(n -1))
                     /  (n*(n-length(inocentes)))  )* qf(1-alpha,length(inocentes),
                                                         n-length(inocentes))
            
            T_i<-mahalanobis(dat[j,inocentes],center=X[inocentes],cov=S[inocentes,inocentes] )
            #Si el T2 con el subvector no es significativo, se termino el proceso
            if(T_i<=ucl){
              break
            }
            
          }else{#se para el ciclo si no hay inocentes
            break
          }
          #aumento del contador del tamaño de subconjuntos
          tamaño<- tamaño+1
        }
      }
      
    }
    
    # Se imprimen los resultados para cada alarma encontrada
    print(round(tabla,3))
    print(paste0("Las variables a las cuales se debe la alarma son: ",paste0(colnames(dat)[culpables] ,collapse=",") ) )
    cat("\n","")
  }
  
  
  # The Murphy out-of-control algorithm ----------------------------------------
  print("------------------Murphy out of control algorithm---------------")
  
  for (j in alarma){
    print(paste0("Para la observación ",j))
    # se crea un vecto nulo para alojar a los culpables en cada paso
    culpables <- NULL
    #Se empieza con todas las variables inocentes y se iran eliminando
    inocentes <- 1:p
    tabla<-NULL
    
    
    while(length(inocentes)>0 ){
      #culpable encontrado en cada paso
      culpables_aux<-NULL
      #vector donde se alojaran las diferencias
      T_i_diff <- NULL
      
      contador<-1
      #Para todos los que siguen siendo inocentes se prueba su distancia al T2
      #completo en presencia a las que ya han sido encontradas culpables
      for(i in inocentes){
        auxiliar<-c(i,culpables)
        auxiliar<-auxiliar[order(auxiliar)]
        T_i_diff[contador] <- T2[j]-mahalanobis(dat[j,auxiliar],center=X[auxiliar],
                                                cov=S[auxiliar,auxiliar] )
        contador <- contador+1
      }
      #controlador de los grados de libertad de la chi cuadrado
      tamaño <- length(auxiliar)
      
      #se haya la variable que esta generando el mínimo
      culpables_aux<-inocentes[which(T_i_diff==min(T_i_diff))[1]]
      
      #se agrega a la lista completa de culpables
      culpables<-c(culpables,culpables_aux)
      
      #Si encontro nuevos culpables debemos quitarlos de la lista de inocentes
      if(length(culpables_aux)>0){
        #Se quita de los inocentes aquellos hallados como culpables
        inocentes<-inocentes[-match(culpables_aux,inocentes)]
      }
      #Inserción en tabla de resultados
      if(length(inocentes>1)){
        tabla_aux<-cbind(T_diff=min(T_i_diff),valor_critico=qchisq(1-alpha,p-tamaño))
        rownames(tabla_aux)<-paste0("T_",paste0(culpables,collapse = ","))
        
        tabla<-rbind(tabla,tabla_aux)
        
      }
      
      #criterio de parada por chi cuadrado
      if(min(T_i_diff)<=qchisq(1-alpha,p-tamaño) ){
        break
      }
      
    }
    print(round(tabla,3))
    print(paste0("Las variables a las cuales se debe la alarma son: ",paste0(colnames(dat)[culpables] ,collapse=",") ) )
    cat("\n","")
  }
  
  #Médodo DFT-----------------------------------------------------------
  print("--------------------------Método DFT-----------------------------")
  #Para todos los casos de alarma
  for(j in alarma){
    #Se imprime la observación alarmante
    print(paste0("Para la observación ",j))
    
    #Cálcule del t_i para cada variable
    t= as.numeric((dat[j,]-X)/sqrt(diag(S)*(1+1/n)))
    
    #Cuantil para luego hallar el k ind(La confianza menor en cada caso)
    T.t.df=pt(t,n-1)
    Kind=abs(2*T.t.df-1)
    
    #Significancia de bonferroni
    Kbonf=(p+Ksim-1)/p
    #color=ifelse(as.vector(Kind>Kbonf),"red","white")
    Diagnóstico=ifelse(as.vector(Kind>Kbonf),"Variable altamente sospechosa","")
    
    #par(mar=c(5.1, 4.1, 4.1, 2.1))
    #barplot(Kind,space=0,names.arg=1:p,col=color,ylim=c(0,max(Kind)+0.3),xlab="i",ylab=expression(K[ind]))
    #legend("topleft","Variables altamente sospechosas",col=2,pch=15,bty="n",pt.cex=2)
    
    #Matriz con los resultados
    res=data.frame(t,Kind,Kbonf,Diagnóstico)
    
    cat("Nivel de confianza nominal",alpha, "\n")
    cat("Nivel de confianza simultáneo ",Ksim, "\n")
    #cat("UCL=",UCL, "\n")
    #cat("Estadístico T2=",T2[j], "\n")
    cat("\n")
    print(res)
  }
  
  return(list(T2=T2,alarma=alarma,p=p,n=n,alpha=alpha))
}



# Prueba de la función ----------------------------------------------------

#generación de un data frame de 100 obs y 2 variables
#de dos normales con medias y varianzas diferentes para fase 1
datos<- matrix(c(rnorm(50,15,4),rnorm(50,10,2),rnorm(50,20,4),rnorm(50,25,4)),nrow=50)

#Agregando nombres random
colnames(datos)<-c("che1","che2","che3","che4")

#Vector de medias simulado
#media<-apply(datos,2,mean)
#Matriz de varianzas simulada
#varianza<- var(datos)


#datos simulados para fase2 con la misma dist que la original
datos.prueba<- matrix(c(rnorm(50,15,4),rnorm(50,10,2),rnorm(50,20,4),rnorm(50,25,4)),nrow=50)
colnames(datos.prueba)<-c("che1","che2","che3","che4")

resultados <- carta.T2.obs.faseI(datos.prueba,datos,0.05)
#DFT.method(datos.prueba,media,varianza,50)

# Prueba 1 nelfy ------------------------------------------------------------


datos1<- read.csv2("DATOS1FASEI.csv")
#vector de medias
#media1<-apply(datos1,2,mean)
#Matriz de varianzas 
#varianza1<- var(datos1)

datos.prueba1<-data.frame(x1=26.73,x2= 18.00 ,x3=2.75)

resultados <- carta.T2.obs.faseI(datos.prueba1,datos1,0.05)

#DFT.method(datos.prueba,media,varianza,50)

# Prueba nelfy 2 ----------------------------------------------------------

datos2<- read.csv2("DATOS2FASEI.csv")
#vector de medias
#media2<-apply(datos2,2,mean)
#Matriz de varianzas 
#varianza2<- var(datos2)

datos.prueba2<-data.frame(x1=1.95,x2=4.56,x3=1.86)

resultados <- carta.T2.obs.faseI(datos.prueba2,datos2,0.05)

#DFT.method(datos.prueba2,media2,varianza2,46)

# Prueba 3 nelfy ----------------------------------------------------------

datos3<- read.csv2("DATOS3FASEI.csv")
#vector de medias
#media3<-apply(datos3,2,mean)
#Matriz de varianzas 
#varianza3<- var(datos3)

datos.prueba3<-data.frame(x1=16.33,x2=33.32,x3=3.86)

resultados <- carta.T2.obs.faseI(datos.prueba3,datos3,0.05)
#DFT.method(datos.prueba3,media3,varianza3,50)

# DFT ---------------------------------------------------------------------

# DFT.method=function(dat,X,S,n,Knom=0.95,Ksim=0.8){
#   print("--------------------------Método DFT-----------------------------")
#   #Para todos los casos de alarma
#   for(j in resultados$alarma){
#     #Se imprime la observación alarmante
#     print(paste0("Para la observación ",j))
#     
#     #Número de variables que se tienen
#     p=ncol(dat)
#     #UCL=(p*(n+1)*(n-1))/(n*(n-p))*qf(Knom,p,n-p)
#     #T2=mahalanobis(x,center=meanHDS,cov=COVHDS)
#     
#     #Cálcule del t_i para cada variable
#     t= as.numeric((dat[j,]-X)/sqrt(diag(S)*(1+1/n)))
#     
#     #Cuantil para luego hallar el k ind(confianza menor en cada caso)
#     T.t.df=pt(t,n-1)
#     Kind=abs(2*T.t.df-1)
#     Kbonf=(p+Ksim-1)/p
#     #color=ifelse(as.vector(Kind>Kbonf),"red","white")
#     Diagnóstico=ifelse(as.vector(Kind>Kbonf),"Variable altamente sospechosa","")
#     
#     #par(mar=c(5.1, 4.1, 4.1, 2.1))
#     #barplot(Kind,space=0,names.arg=1:p,col=color,ylim=c(0,max(Kind)+0.3),xlab="i",ylab=expression(K[ind]))
#     #legend("topleft","Variables altamente sospechosas",col=2,pch=15,bty="n",pt.cex=2)
#     
#     res=data.frame(t,Kind,Kbonf,Diagnóstico)
#     
#     #res=data.frame(i=1:p,nrow=1,as.vector(t),as.vector(Kind),as.vector(Kbonf),as.vector(Diagnóstico))
#     cat("Nivel de confianza nominal",Knom, "\n")
#     cat("Nivel de confianza simultáneo ",Ksim, "\n")
#     #cat("UCL=",UCL, "\n")
#     #cat("Estadístico T2=",T2[j], "\n")
#     cat("\n")
#     print(res)
#   }
#   
# }


#Prueba nelfy del pdf

xnew=data.frame(x1=15,x2=10,x3=20,x4=-5)
meanxref=c(0,0,0,0)
Sref=matrix(c(102.74,88.67,67.04,54.06,88.67,142.74,86.56,80.03,67.04,86.56,84.57,
69.42,54.06,80.03,69.42,99.06),ncol=4,byrow=T)

resultados<-carta.T2.obs.faseI(xnew,meanxref,Sref,n=40)
#DFT.method(dat=xnew,X=meanxref,COVHDS=Sref,n=40)

