# source('C:/Ecology/Drive/Research/Andropogon/Andropogon - Pheno-geno Modeling/Code/TEMP.r')

	### plot
	########
	
	# colors
	cols <- cividis(101)
		
	dirCreate('./Figures & Tables/SNP Models - Genetic Diversity')

	### generazliation
	thold <- 0.9 # topmost quantile of psi above which is designated as range "core"
	
	
	highQuant <- quantile(lopodPredictVarP@data$psi_i, thold, na.rm=TRUE)
	
	# plot
	png(paste0('./Figures & Tables/SNP Models - Genetic Diversity/Genetic Diversity.png'), width=2 * 1200, height=1000, res=300)

		par(mfrow=c(1, 2), oma=0.1 * c(1, 4, 1, 26), mai=0 * c(1, 1, 1, 1), mar=rep(0, 4))

		### present
		###########

		# diversity
		plot(sqDiversity, col=cols, breaks=seq(0, 1, by=0.01), legend=FALSE, axes=FALSE, box=FALSE)
		plot(nam1, add=TRUE)

		# range core
		highQuant <- quantile(lopodPredictVarP$psi_i, thold, na.rm=TRUE)
		whichCore <- which(lopodPredictVarP@data$psi_i >= highQuant)
		highPred <- lopodPredictVarP[whichCore, ]
		highPred <- gUnaryUnion(highPred)
		plot(highPred, col=NA, border='white', lwd=1, add=TRUE)

		# fade by probability of occurrence
		plot(lopodPredictVarP, col=alpha('white', 1 - lopodPredictVarP$psi_i), border=NA, add=TRUE)
		
		# add populations
		popCol <- cols[round(100 * obsDiversity)]
		points(phenoGenoPops, pch=21, bg=popCol, cex=1.2)
		
		# title
		usr <- par('usr')
		x <- usr[1]
		y <- usr[4] - 0.1 * (usr[4] - usr[3])
		text(x, y, labels='a) Current climate', xpd=NA, cex=0.7, pos=4, font=2)
		
		# legend
	
		### future
		###########

		# diversity
		plot(futDiversity, col=cols, breaks=seq(0, 1, by=0.01), legend=FALSE, axes=FALSE, box=FALSE)
		plot(nam1, add=TRUE)

		# range core
		whichCore <- which(lopodPredictVarP@data$psi_ensembleGcm_rcp85_2070s >= highQuant)
		highPred <- lopodPredictVarP[whichCore, ]
		highPred <- gUnaryUnion(highPred)
		plot(highPred, col=NA, border='white', lwd=1.5, add=TRUE)

		# fade by probability of occurrence
		plot(lopodPredictVarP, col=alpha('white', 1 - lopodPredictVarP$psi_ensembleGcm_rcp85_2070s), border=NA, add=TRUE)
		
		# add populations
		popCol <- cols[round(100 * obsDiversity)]
		points(phenoGenoPops, pch=21, bg=popCol, cex=1.2)

		# title
		usr <- par('usr')
		x <- usr[1]
		y <- usr[4] - 0.1 * (usr[4] - usr[3])
		text(x, y, labels='b) 2070s RCP8.5', xpd=NA, cex=0.7, pos=4, font=2)
		
		# legend
		legendGrad(
			x='right',
			inset = -0.17,
			vert = TRUE,
			width = 0.12,
			height = 1,
			labels = c(0, 0.25, 0.5, 0.75, 1),
			labAdj = 0.35,
			cex = 0.65,
			col = cols,
			border = 'black',
			title = 'Diversity',
			titleAdj = c(0.5, 0.855),
			adjX = c(0, 0.35),
			# adjY = c(0.23, 0.775),
			adjY = c(0.13, 0.775),
			boxBg = par('bg'),
			boxBorder = NULL
			# swatches = list(list(swatchAdjY=c(0.13, 0.175), col=NA, border='orange', lwd=1.5, labels='Range\ncore'))
		)	

		title(sub=date(), cex.sub=0.35, line=-1, outer=TRUE)
		
	dev.off()
