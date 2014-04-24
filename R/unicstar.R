unicstar <-
function(datos_j){
#                                 
#	This function identifies the sample ties \eqn{(y^{*}_k)} in \eqn{\mathcal{S}_j} and computes their associated frequencies.
#   
#	Input:
#		datos_j		-	Data matrix with features and number of individuals in the sample 'S_j'
#						
#	Output. The function 'unicstar' produces an object list with two elements:
#		U		    - 	Number of sample ties in 'S_j'
#		unicstar	-	Data matrix with sample ties 'y_i' and associated frequencies 
#

y_i <- as.data.frame(datos_j)
y_i <- as.data.frame(datos_j[,"y_i"])
colnames(y_i) <- c("y_i")

# 	Identifying sample ties
ystar <- sort(as.matrix(unique(y_i)))
unicstar <- cbind(ystar,ystar)
colnames(unicstar) <- c("ystar","m")
dimUnicos <- dim(unicstar)

# 	Computing frequencies 
k <- 1
for(k in 1:dimUnicos[1]){
	ids <- which(datos_j[,"y_i"]==unicstar[k,"ystar"])
	unicstar[k,"m"] <- sum(datos_j[ids,"n_i"])
}
# 	Counting the number sample ties
U <- dimUnicos[1]

#	Output
salida <- list(U,unicstar)
return(salida)

#
#  --   END of unicstar.R  --
}
