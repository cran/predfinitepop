conteo <-
function(datos_j,N_j){
#
#	This function counts the number of individuals in \eqn{\mathcal{S}_j} and 
#	\eqn{widetilde{\mathcal{S}}_j}, with respect to the total number of individuals 
#	in a given planned domain \eqn{j}. 
#
#	Input:
#		"datos_j" 	-	This object should contain two columns labelled: "n_i" and "domplan".  
#						"n_i"     - Number of individuals in the ith group of individuals, and 
#						"domplan" - Column vector with categories for the planned domains
#		"N_j"  		-	This object should include two columns labelled: "N_j" and "domplan".
#						"N_j"     - Number of individuals in the jth planned domain, and 
#						"domplan" - Column vector with categories for the planned domains
#
#	Output:
#		M			-	Number of groups in "S_j"
#		CardS		-	Number of individuals in "S_j" (if "n_i"=1 for any "i", then Cards=M=nrow(datos_j))
#		CardNoS		-	Number of individuals out of sample "Stilde_j"
#

#	 Counting the number of individuals in the sample
CardS <- as.matrix(sum(datos_j[,"n_i"]))
colnames(CardS) <- c("CardS")

#	 Counting the number of individuals outside the sample
CardNoS <- as.matrix(N_j - CardS)
colnames(CardNoS) <- c("CardNoS")

#	Number of groups in the sample
M <- as.matrix(nrow(datos_j))
colnames(M) <- c("M")

#	Output
salida <- list(M,CardS,CardNoS)
return(salida)

#
#  --  END of conteo --
}
