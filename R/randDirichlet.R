randDirichlet <-
function(alpha,n){
#
#   This function simulates samples from the Dirichlet distribution with (px1) vector parameter \eqn{\alpha}.
#
#	Input:
#		alpha	-	p-dimensional vector with positive entries
#		n		-	Number of simulated replicates
#
#	Output:
#		randDir	-	(Nxp) matrix with 'n' simulated replicates
#
#	References:
#	-	"Non-Uniform Random Variate Generation", Berlin: Springer-Verlag. Devrog, Luz (1986)
#

#	Dimension "p"
p <- length(alpha)

#	Repository
randDir <- matrix(NaN,n,p)
pAux <- NaN * rep(1:p)

#	Simulation
for(i in 1:n)
{
	for(k in 1:p)
	{
		pAux[k] <- rgamma(1, shape = alpha[k], scale = 1)
	}
	randDir[i,] <- pAux / sum(pAux)
}

#	Output
salida <- randDir
return(salida)

#
#	--	END of randDirichlet.R --
}
