pesosDP <-
function(unicos_sal,alphaDP){
#
#	This function computes the weights associated with the predictive distribution \eqn{\hat{G}_j} for a given planned domain \eqn{j}, 
#	using the sample ties \eqn{(y^{*}j)}, under the Dircihlet processes prior.
#   
#	Input:
#		unicos_sal	-	(U x 2) matrix with sample ties and associated frequencies (produced with 'unicos')
#		alphaDP		-	Positive scalar for the Dirichlet process
#                       
#	Details:
#		The matrix with sample ties 'unicos_sal' must include a column named 'm_k' for the frequencies of the sample ties.
#
#	Output:
#		rho		-	(U x 2) dimensional vector with weights associated with 'ystar' (the sample ties)
#		phi		-	Probability weight associated with 'G_{j0}' (the continuous part of 'Ghat_j')
#

# computing "rho" and "phi"
unicstar <- unicos_sal[[2]]
U <- unicos_sal[[1]]
rho <- NaN * unicstar[,"m"]
M <- sum(unicstar[,"m"])
for(l in 1:U){
	rho[l] <- unicstar[l,"m"]/(M+alphaDP)
	}
phi <- as.matrix(alphaDP/(M+alphaDP))
rho <- as.matrix(rho)
colnames(rho) <- c("rho")
colnames(phi) <- c("phi")

#	Output
salida <- list(rho,phi)
return(salida)

#
#  --  END of pesosDP.R  --
}
