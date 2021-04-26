library(reshape2)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(qtl2)
library(AlphaSimR)
library(magicdesign)
library(bindata)

### set working directory.
setwd("")

### functions needed.
{
  # function to create files in qtl2 format (input marker genotypes should be 0/1/2/NA).
  format.qtl2 <- function(map, geno.founder, geno.ril, filename){
    
    # get the number of founders.
    n <- ncol(geno.founder)
    
    # create the genetic map file.
    colnames(map) <- c("marker", "chr", "pos")
    
    # create the founder geno file.
    geno.founder[geno.founder==0] <- "A"
    geno.founder[geno.founder==1] <- "H"
    geno.founder[geno.founder==2] <- "B"
    geno.founder[is.na(geno.founder)] <- "-"
    geno.founder <- cbind(map[,1], geno.founder)
    colnames(geno.founder) <- c("marker", LETTERS[1:n])
    
    # create the RILs geno file.
    geno.ril[geno.ril==0] <- "A"
    geno.ril[geno.ril==1] <- "H"
    geno.ril[geno.ril==2] <- "B"
    geno.ril[is.na(geno.ril)] <- "-"
    geno.ril <- cbind(map[,1], geno.ril)
    colnames(geno.ril)[1] <- "marker"
    
    # export the files.
    write.csv(map, paste("gmap_", filename, ".csv", sep=""), quote=F, row.names=F)
    write.csv(geno.founder, paste("geno_founder_", filename, ".csv", sep=""), quote=F, row.names=F)
    write.csv(geno.ril, paste("geno_ril_", filename, ".csv", sep=""), quote=F, row.names=F)
    
  }
  
  # function to convert simulated marker data to founder calls.
  magic.hap <- function(hap, n){
    
    # get the marker indices that correspond to each founder allele.
    idx <- lapply(1:(n-1), FUN=function(x) seq(x, ncol(hap), n-1))
    
    # create a new haplotype matrix.
    out <- matrix(1, nrow=nrow(hap), ncol=ncol(hap)/(n-1))
    
    # convert the haplotypes into parental genotypes.
    for(i in 1:ncol(out)){
      for(j in 1:length(idx)){
        out[hap[,idx[[j]][i]]==1,i] <- j + 1
      }
    }
    
    return(out)
    
  }
  
  # function to simulate and count the recombinant haplotypes in each population.
  magic.posthoc <- function(crossfile, fcall, self, n.sim){
    
    # read the marker files using qtl2.
    crossfile <- read_cross2(crossfile)
    
    # PART 1.
    {
      
      # get gmap.
      gmap <- crossfile$gmap
      gmap <- lapply(1:length(gmap), FUN=function(x) unname(gmap[[x]]))
      
      # get cross info (xinfo).
      xinfo <- crossfile$cross_info
      
      # get the number of founder (n) and crossing generation (nx).
      n <- ncol(crossfile$cross_info)
      nx <- log(n, 2)
      
      # convert xinfo into xplan.
      xplan <- list(t(matrix(c(t(xinfo)), nrow=2)))
      for(i in 2:nx) xplan[[i]] <- matrix(1:nrow(xplan[[i-1]]), ncol=2, byrow=T)
      
      # add selfing generations.
      xplan2 <- list()
      for(i in 1:nx){
        if(self[i] > 0){
          xplan2 <- c(xplan2, xplan[i], replicate(self[i], list(cbind(1:nrow(xplan[[i]]), 1:nrow(xplan[[i]])))))
        } else {
          xplan2 <- c(xplan2, xplan[i])
        }
      }
      xplan <- xplan2
      
      # get the number of generations.
      n.gen <- length(xplan)
      
      # create founder haplotypes.
      n.chr <- length(gmap)
      n.marker <- sapply(1:n.chr, FUN=function(x) length(gmap[[x]]))
      gen.map <- lapply(1:n.chr, FUN=function(x) (gmap[[x]] - min(gmap[[x]]))/100)
      gen.map <- lapply(1:n.chr, FUN=function(x) sort(rep(unname(gen.map[[x]]), n-1)))
      fhap <- lapply(1:n.chr, FUN=function(x) matrix(rep(rbind(0, diag(n-1)), n.marker[x]), nrow=n, ncol=n.marker[x]*(n-1)))
      
      # objects to store the summary outputs from each round of simulation.
      out.rseg <- replicate(n.chr, vector())
      
      # simulate MAGIC population for n.sim-times.
      for(i in 1:n.sim){
        
        # create the founder population.
        founder <- newMapPop(genMap=gen.map,
                             haplotypes=fhap,
                             inbred=T,
                             ploidy=2L)
        SP <- SimParam$new(founder)
        sim <- newPop(founder, simParam=SP)
        
        for(j in 1:n.gen){
          sim <- makeCross(pop=sim,
                           crossPlan=xplan[[j]],
                           nProgeny=1,
                           simParam=SP)
        }
        
        # extract the haploid marker data for all chromosomes.
        hap <- lapply(1:n.chr, FUN=function(x) pullSegSiteHaplo(pop=sim, chr=x, simParam=SP))
        hap <- lapply(1:n.chr, FUN=function(x) hap[[x]][seq(1,nrow(hap[[x]]),2),,drop=F])
        hap <- lapply(1:n.chr, FUN=function(x) magic.hap(hap=hap[[x]], n=n))
        
        # keep the proportion of each founder in each simulation.
        for(j in 1:n.chr){
          
          # identify the marker segment (mseg) for each individual.
          mseg <- vector()
          for(k in 1:nrow(hap[[j]])){
            temp <- which(!(hap[[j]][k, -1] == hap[[j]][k, -ncol(hap[[j]])]))
            if(length(temp) > 0){
              temp <- cbind(k, c(1,temp+1), c(temp, ncol(hap[[j]])))
            } else {
              temp <- cbind(k, 1, ncol(hap[[j]]))
            }
            mseg <- rbind(mseg, temp)
          }
          
          # identify the length of each mseg.
          temp <- gmap[[j]]
          temp <- c(temp[1], (temp[-1] + temp[-length(temp)])/2, temp[length(temp)])
          mseg <- cbind(mseg, 
                        temp[mseg[,3] + 1] - temp[mseg[,2]])
          
          # identify the founder of each mseg.
          mseg <- cbind(mseg, 
                        sapply(1:nrow(mseg), FUN=function(x) hap[[j]][mseg[x,1], mseg[x,2]]))
          attr(mseg, "dimnames") <- NULL
          
          # calculate the number of recombinant haplotype.
          rhap <- cbind(sort(rep(1:n, n)), rep(1:n,n))
          rhap <- rhap[!(rhap[,1]==rhap[,2]), ]
          rseg <- matrix(0, nrow=nrow(hap[[j]]), ncol=nrow(rhap))
          for(k in 1:nrow(rhap)){
            temp <- mseg[mseg[-nrow(mseg),5]==rhap[k,1] & mseg[-1,5]==rhap[k,2] & mseg[-nrow(mseg),1]==mseg[-1,1], 1]
            temp1 <- unique(temp)
            temp2 <- unname(c(table(temp)))
            rseg[temp1, k] <- temp2
          }
          rseg <- colSums(rseg)/nrow(rseg)
          out.rseg[[j]] <- cbind(out.rseg[[j]], rseg)
          
        }
        
      }
      
      
      # remove unnecessary attributes.
      for(j in 1:n.chr) attr(out.rseg[[j]], "dimnames") <- NULL
      
    }
    
    # PART 2.
    {
      
      # get the real founder genotype.
      geno <- crossfile$founder_geno
      
      # create a new crossfile for simulated RIL.
      crossfile_sim <- crossfile
      
      # convert the true simulated RIL hap to marker based on founder genotype.
      # set heterozygous founder genotypes (2) to missing (NA).
      hap2geno <- replicate(n.chr, vector())
      for(i in 1:n.chr){
        for(j in 1:ncol(hap[[i]])){
          geno[[i]][geno[[i]] == 2] <- NA
          hap2geno[[i]] <- cbind(hap2geno[[i]], geno[[i]][hap[[i]][,j], j])
        }
        attr(hap2geno[[i]], "dimnames") <- attr(crossfile$geno[[i]], "dimnames")
      }
      
      # replace the real RIL genotype in crossfile_sim with simulated RIL genotype.
      crossfile_sim$geno <- hap2geno
      
      # calculate the genotype probability from the simulated RIL genotype.
      gp <- calc_genoprob(cross=crossfile_sim, error_prob=0.01, map_function="haldane", cores=0)
      
      # call the founder genotypes in the RILs with a minimum probability of 0.5.
      fcall_sim <- maxmarg(probs=gp, minprob=0.5001, cores=0)
      
      # identify the xo positions using qtl2 function.
      xo <- locate_xo(geno=fcall_sim, map=crossfile_sim$gmap, cores=0)
      
      # objects to store the outputs from each chromosome.
      out.rseg2 <- list()
      
      for(j in 1:n.chr){
        
        # get the fcall for j-th chr and set NA to 0.
        hap <- fcall_sim[[j]]
        hap[is.na(hap)] <- 0
        attr(hap, "dimnames") <- NULL
        
        # identify the marker segment (mseg) for each individual.
        mseg <- vector()
        for(k in 1:nrow(hap)){
          temp <- which(!(hap[k, -1] == hap[k, -ncol(hap)]))
          if(length(temp) > 0){
            temp <- cbind(k, c(1,temp+1), c(temp, ncol(hap)))
          } else {
            temp <- cbind(k, 1, ncol(hap))
          }
          mseg <- rbind(mseg, temp)
        }
        
        # identify the founder of each mseg.
        mseg <- cbind(mseg, 
                      sapply(1:nrow(mseg), FUN=function(x) hap[mseg[x,1], mseg[x,2]]))
        attr(mseg, "dimnames") <- NULL
        
        # remove any segment that is 0 (remove NA segment).
        mseg <- mseg[!(mseg[,4]==0),]
        
        # merge any segment that are the same (e.g. 5-0-5 becomes 5-5 previously, and now 5-5 becomes 5).
        temp1 <- which(mseg[-1, 4] == mseg[-nrow(mseg), 4] & mseg[-1, 1] == mseg[-nrow(mseg), 1])
        if(length(temp1) > 0){
          temp2 <- mseg[temp1 + 1, 3]
          mseg[temp1, 3] <- temp2
          mseg <- mseg[-c(temp1 + 1), ]
        }
        
        # identify the length of each mseg.
        mseg <- cbind(mseg, 0)
        for(k in unique(mseg[,1])){
          temp <- c(xo[[j]][[k]], max(gmap[[j]])) - c(0, xo[[j]][[k]])
          mseg[mseg[,1]==k, 5] <- temp
        }
        
        # rearrange mseg so the columns are the same as in the simulation (prevent confusion).
        mseg <- mseg[,c(1,2,3,5,4)]
        
        # calculate the number of recombinant haplotype.
        rhap <- cbind(sort(rep(1:n, n)), rep(1:n,n))
        rhap <- rhap[!(rhap[,1]==rhap[,2]), ]
        rseg <- matrix(0, nrow=nrow(hap), ncol=nrow(rhap))
        for(k in 1:nrow(rhap)){
          temp <- mseg[mseg[-nrow(mseg),5]==rhap[k,1] & mseg[-1,5]==rhap[k,2] & mseg[-nrow(mseg),1]==mseg[-1,1], 1]
          temp1 <- unique(temp)
          temp2 <- unname(c(table(temp)))
          rseg[temp1, k] <- temp2
        }
        rseg <- t(rseg)
        out.rseg2 <- c(out.rseg2, list(rseg))
        
      }
      
    }
    
    # PART 3.
    {
      # identify the xo positions using qtl2 function.
      xo <- locate_xo(geno=fcall, map=crossfile$gmap, cores=0)
      
      # objects to store the outputs from each chromosome.
      out.rseg3 <- list()
      
      for(j in 1:n.chr){
        
        # get the fcall for j-th chr and set NA to 0.
        hap <- fcall[[j]]
        hap[is.na(hap)] <- 0
        attr(hap, "dimnames") <- NULL
        
        # identify the marker segment (mseg) for each individual.
        mseg <- vector()
        for(k in 1:nrow(hap)){
          temp <- which(!(hap[k, -1] == hap[k, -ncol(hap)]))
          if(length(temp) > 0){
            temp <- cbind(k, c(1,temp+1), c(temp, ncol(hap)))
          } else {
            temp <- cbind(k, 1, ncol(hap))
          }
          mseg <- rbind(mseg, temp)
        }
        
        # identify the founder of each mseg.
        mseg <- cbind(mseg, 
                      sapply(1:nrow(mseg), FUN=function(x) hap[mseg[x,1], mseg[x,2]]))
        attr(mseg, "dimnames") <- NULL
        
        # remove any segment that is 0 (remove NA segment).
        mseg <- mseg[!(mseg[,4]==0),]
        
        # merge any segment that are the same (e.g. 5-0-5 becomes 5-5 previously, and now 5-5 becomes 5).
        temp1 <- which(mseg[-1, 4] == mseg[-nrow(mseg), 4] & mseg[-1, 1] == mseg[-nrow(mseg), 1])
        if(length(temp1) > 0){
          temp2 <- mseg[temp1 + 1, 3]
          mseg[temp1, 3] <- temp2
          mseg <- mseg[-c(temp1 + 1), ]
        }
        
        # identify the length of each mseg.
        mseg <- cbind(mseg, 0)
        for(k in unique(mseg[,1])){
          temp <- c(xo[[j]][[k]], max(gmap[[j]])) - c(0, xo[[j]][[k]])
          mseg[mseg[,1]==k, 5] <- temp
        }
        
        # rearrange mseg so the columns are the same as in the simulation (prevent confusion).
        mseg <- mseg[,c(1,2,3,5,4)]
        
        # calculate the number of recombinant haplotype.
        rhap <- cbind(sort(rep(1:n, n)), rep(1:n,n))
        rhap <- rhap[!(rhap[,1]==rhap[,2]), ]
        rseg <- matrix(0, nrow=nrow(hap), ncol=nrow(rhap))
        for(k in 1:nrow(rhap)){
          temp <- mseg[mseg[-nrow(mseg),5]==rhap[k,1] & mseg[-1,5]==rhap[k,2] & mseg[-nrow(mseg),1]==mseg[-1,1], 1]
          temp1 <- unique(temp)
          temp2 <- unname(c(table(temp)))
          rseg[temp1, k] <- temp2
        }
        rseg <- t(rseg)
        out.rseg3 <- c(out.rseg3, list(rseg))
        
      }
      
    }
    
    # return the simulation and real data output.
    return(list(out.rseg, out.rseg2, out.rseg3))
    
  }
  
  # function to plot recombinant haplotypes from magic.posthoc output (individual rec hap).
  magic.plot.rseg <- function(dat, n, n.sim=100, ylim=NULL){
    
    dat.plot1 <- dat[[1]][[1]]
    for(i in 2:length(dat[[1]])){
      dat.plot1 <- dat.plot1 + dat[[1]][[i]]
    }
    dat.plot1 <- data.frame(t(dat.plot1), sim=1:n.sim)
    dat.plot1 <- melt(dat.plot1, id.vars=c("sim"))
    colnames(dat.plot1)[2] <- "rhap"
    
    temp <- cbind(sort(rep(1:n, n)), rep(1:n,n))
    temp <- temp[!(temp[,1]==temp[,2]), ]
    temp <- paste(temp[,1], temp[,2], sep="_")
    
    dat.plot1$rhap <- factor(dat.plot1$rhap, labels=temp)
    
    dat.plot2 <- dat[[2]][[1]]
    for(i in 2:length(dat[[2]])){
      dat.plot2 <- dat.plot2 + dat[[2]][[i]]
    }
    dat.plot2 <- data.frame(value=rowSums(dat.plot2)/ncol(dat.plot2), rhap=1:length(temp), pop="sim")
    
    dat.plot3 <- dat[[3]][[1]]
    for(i in 2:length(dat[[3]])){
      dat.plot3 <- dat.plot3 + dat[[3]][[i]]
    }
    dat.plot3 <- data.frame(value=rowSums(dat.plot3)/ncol(dat.plot3), rhap=1:length(temp), pop="actual")
    
    dat.anno <- rbind(dat.plot2, dat.plot3)
    dat.anno$rhap <- as.factor(dat.anno$rhap)
    dat.anno$rhap <- factor(dat.anno$rhap, labels=temp)
    dat.anno$pop <- as.factor(dat.anno$pop)
    dat.anno$pop <- factor(dat.anno$pop, levels=c("actual", "sim"))
    
    if(is.null(ylim)) ylim <- max(c(dat.plot1$value, dat.anno$value))
    
    ggplot() +
      scale_x_discrete() +
      annotate("rect", xmin=seq(0.5, length(temp), 2), xmax=seq(1.5, length(temp), 2), ymin=-Inf, ymax=Inf, fill="#EEEEEE", color=NA) +
      annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#CCCCCC") +
      geom_boxplot(data=dat.plot1, aes(x=rhap, y=value), outlier.size=0.2, lwd=0.2) +
      geom_point(data=dat.anno, aes(x=rhap, y=value, color=pop), size=1, alpha=0.8) +
      theme(panel.background=element_blank(), panel.grid=element_blank()) +
      theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
      scale_y_continuous(limits=c(0, ylim)) +
      scale_color_manual(values=c("#FF5555", "#5555FF")) +
      guides(color=guide_legend(override.aes=list(size=2))) +
      xlab("recombinant haplotype") +
      ylab("mean count") +
      theme(legend.position="bottom")
    
  }
  
  # function to plot recombinant haplotypes for wheat-UK16 specifically (individual rec hap).
  magic.plot.rseg2 <- function(dat, n, n.sim=100){
    
    dat.plot1 <- dat[[1]][[1]]
    for(i in 2:length(dat[[1]])){
      dat.plot1 <- dat.plot1 + dat[[1]][[i]]
    }
    dat.plot1 <- data.frame(t(dat.plot1), sim=1:n.sim)
    dat.plot1 <- melt(dat.plot1, id.vars=c("sim"))
    colnames(dat.plot1)[2] <- "rhap"
    
    temp <- cbind(sort(rep(1:n, n)), rep(1:n,n))
    temp <- temp[!(temp[,1]==temp[,2]), ]
    temp <- paste(temp[,1], temp[,2], sep="_")
    
    dat.plot1$rhap <- factor(dat.plot1$rhap, labels=temp)
    dat.plot1$group <- 4
    dat.plot1$group[dat.plot1$rhap%in%temp[1:60]] <- 1
    dat.plot1$group[dat.plot1$rhap%in%temp[61:120]] <- 2
    dat.plot1$group[dat.plot1$rhap%in%temp[121:180]] <- 3
    
    dat.plot3 <- dat[[3]][[1]]
    for(i in 2:length(dat[[3]])){
      dat.plot3 <- dat.plot3 + dat[[3]][[i]]
    }
    dat.plot3 <- data.frame(value=rowSums(dat.plot3)/ncol(dat.plot3), rhap=1:length(temp))
    
    dat.plot3$rhap <- as.factor(dat.plot3$rhap)
    dat.plot3$rhap <- factor(dat.plot3$rhap, labels=temp)
    dat.plot3$group <- c(rep(1,60), rep(2,60), rep(3,60), rep(4,60))
    
    ylim <- max(c(dat.plot1$value, dat.plot3$value))
    
    ggplot() +
      scale_x_discrete() +
      annotate("rect", xmin=seq(0.5, length(temp)/4, 2), xmax=seq(1.5, length(temp)/4, 2), ymin=-Inf, ymax=Inf, fill="#EEEEEE", color=NA) +
      annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#CCCCCC") +
      geom_boxplot(data=dat.plot1, aes(x=rhap, y=value), outlier.size=0.2, lwd=0.2) +
      geom_point(data=dat.plot3, aes(x=rhap, y=value), color="#FF7777", size=1, alpha=0.8) +
      facet_wrap(vars(group), scales="free_x", nrow=4) +
      theme(strip.background=element_blank(), strip.text=element_blank()) +
      theme(panel.background=element_blank(), panel.grid=element_blank()) +
      theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
      scale_y_continuous(limits=c(0, ylim)) +
      xlab("recombinant haplotype") +
      ylab("mean count")
    
  }
  
  # function to plot recombinant haplotypes from magic.posthoc output (matrix plot).
  magic.plot.rseg3 <- function(dat, n, n.sim=100, lim=NULL){
    
    temp <- cbind(sort(rep(1:n, n)), rep(1:n,n))
    temp <- temp[!(temp[,1]==temp[,2]), ]
    colnames(temp) <- c("y", "x")
    
    dat.plot1 <- dat[[1]][[1]]
    for(i in 2:length(dat[[1]])){
      dat.plot1 <- dat.plot1 + dat[[1]][[i]]
    }
    dat.plot1 <- data.frame(temp,
                            value=rowSums(dat.plot1)/ncol(dat.plot1))

    dat.plot3 <- dat[[3]][[1]]
    for(i in 2:length(dat[[3]])){
      dat.plot3 <- dat.plot3 + dat[[3]][[i]]
    }
    dat.plot3 <- data.frame(temp, value=rowSums(dat.plot3)/ncol(dat.plot3))
    
    if(is.null(lim)) lim <- max(c(dat.plot1$value, dat.plot3$value))
    
    ggplot() +
      annotate("segment", x=0.5, xend=n+0.5, y=seq(0.5, n+0.5, 1), yend=seq(0.5, n+0.5, 1), color="#DDDDDD") +
      annotate("segment", x=seq(0.5, n+0.5, 1), xend=seq(0.5, n+0.5, 1), y=0.5, yend=n+0.5, color="#DDDDDD") +
      ggforce::geom_circle(data=dat.plot3, aes(x0=x, y0=y, r=sqrt(value/(lim*4*1.5)), color=value, fill=value)) +
      ggforce::geom_circle(data=dat.plot1, aes(x0=x, y0=y, r=sqrt(value/(lim*4*1.5))), alpha=0.5) +
      scale_y_reverse(breaks=1:n, expand=c(0,0)) +
      scale_x_continuous(breaks=1:n, expand=c(0,0), position="top") +
      scale_fill_gradient(low="#EEEEEE", high="#FF0000", limits=c(0, lim), name="mean count") +
      scale_color_gradient(low="#EEEEEE", high="#FF0000", limits=c(0, lim), name="mean count") +
      coord_fixed() +
      theme(panel.background=element_blank(), panel.grid=element_blank()) +
      theme(axis.ticks=element_blank(), axis.title=element_blank()) +
      theme(axis.text=element_text(size=8)) +
      theme(legend.title=element_text(size=8), legend.text=element_text(size=7)) +
      theme(legend.position="bottom")
      
  }
  
  # function to calculate mean and sd for total recombinations.
  magic.totalrec <- function(dat, g){
    
    a <- dat[[1]]
    a <- lapply(1:length(a), FUN=function(x) colSums(a[[x]]))
    b <- rowSums(sapply(1:length(a), FUN=function(x) a[[x]]))
    a.1 <- round(sapply(1:length(a), FUN=function(x) c(mean(a[[x]]), sd(a[[x]]))),2)
    a.2 <- round(sapply(1:length(a), FUN=function(x) c(mean(a[[x]]/g[x]), sd(a[[x]]/g[x]))),2)
    b.1 <- round(rbind(mean(b), sd(b)), 2)
    b.2 <- round(rbind(mean(b/sum(g)), sd(b/sum(g))), 2)
    out <- data.frame(Chr=c(names(g), "All"), t(cbind(a.1, b.1)), t(cbind(a.2, b.2)))
    colnames(out)[2:5] <- c("simpop_mean", "simpop_sd", "simpop_mean_per_M", "simpop_sd_per_M")
    
    a <- dat[[2]]
    a <- lapply(1:length(a), FUN=function(x) colSums(a[[x]]))
    b <- rowSums(sapply(1:length(a), FUN=function(x) a[[x]]))
    a.1 <- round(sapply(1:length(a), FUN=function(x) c(mean(a[[x]]), sd(a[[x]]))),2)
    a.2 <- round(sapply(1:length(a), FUN=function(x) c(mean(a[[x]]/g[x]), sd(a[[x]]/g[x]))),2)
    b.1 <- round(rbind(mean(b), sd(b)), 2)
    b.2 <- round(rbind(mean(b/sum(g)), sd(b/sum(g))), 2)
    out <- data.frame(out, t(cbind(a.1, b.1)), t(cbind(a.2, b.2)))
    colnames(out)[6:9] <- c("sim_mean", "sim_sd", "sim_mean_per_M", "sim_sd_per_M")
    
    a <- dat[[3]]
    a <- lapply(1:length(a), FUN=function(x) colSums(a[[x]]))
    b <- rowSums(sapply(1:length(a), FUN=function(x) a[[x]]))
    a.1 <- round(sapply(1:length(a), FUN=function(x) c(mean(a[[x]]), sd(a[[x]]))),2)
    a.2 <- round(sapply(1:length(a), FUN=function(x) c(mean(a[[x]]/g[x]), sd(a[[x]]/g[x]))),2)
    b.1 <- round(rbind(mean(b), sd(b)), 2)
    b.2 <- round(rbind(mean(b/sum(g)), sd(b/sum(g))), 2)
    out <- data.frame(out, t(cbind(a.1, b.1)), t(cbind(a.2, b.2)))
    colnames(out)[10:13] <- c("actual_mean", "actual_sd", "actual_mean_per_M", "actual_sd_per_M")
    
    return(out)
  }
  
  # function to calculate mean and sd for total recombinations, wheat-UK16 specifically.
  magic.totalrec2 <- function(dat, g){
    a <- dat[[1]]
    a <- lapply(1:length(a), FUN=function(x) colSums(a[[x]]))
    b <- rowSums(sapply(1:length(a), FUN=function(x) a[[x]]))
    a.1 <- round(sapply(1:length(a), FUN=function(x) c(mean(a[[x]]), sd(a[[x]]))),2)
    a.2 <- round(sapply(1:length(a), FUN=function(x) c(mean(a[[x]]/g[x]), sd(a[[x]]/g[x]))),2)
    b.1 <- round(rbind(mean(b), sd(b)), 2)
    b.2 <- round(rbind(mean(b/sum(g)), sd(b/sum(g))), 2)
    out <- data.frame(Chr=c(names(g), "All"), t(cbind(a.1, b.1)), t(cbind(a.2, b.2)))
    colnames(out)[2:5] <- c("simpop_mean", "simpop_sd", "simpop_mean_per_M", "simpop_sd_per_M")
    
    a <- dat[[3]]
    a <- lapply(1:length(a), FUN=function(x) colSums(a[[x]]))
    b <- rowSums(sapply(1:length(a), FUN=function(x) a[[x]]))
    a.1 <- round(sapply(1:length(a), FUN=function(x) c(mean(a[[x]]), sd(a[[x]]))),2)
    a.2 <- round(sapply(1:length(a), FUN=function(x) c(mean(a[[x]]/g[x]), sd(a[[x]]/g[x]))),2)
    b.1 <- round(rbind(mean(b), sd(b)), 2)
    b.2 <- round(rbind(mean(b/sum(g)), sd(b/sum(g))), 2)
    out <- data.frame(out, t(cbind(a.1, b.1)), t(cbind(a.2, b.2)))
    colnames(out)[6:9] <- c("actual_mean", "actual_sd", "actual_mean_per_M", "actual_sd_per_M")
    
    return(out)
  }
  
  # function to calculate the percent of rec-hap recovered.
  magic.rhap.recovery <- function(dat){
    a <- dat[[1]][[1]]
    for(i in 2:length(dat[[1]])) a <- a + dat[[1]][[i]]
    a <- rowSums(a)/ncol(a)
    
    b <- dat[[3]][[1]]
    for(i in 2:length(dat[[3]])) b <- b + dat[[3]][[i]]
    b <- rowSums(b)/ncol(b)
    
    print(data.frame(mean=round(mean(b/a),3), sd=round(sd(b/a),3)), row.names=F)
  }
  
  # function to simulate and calculate percent of rec-hap recovered for different marker density.
  # chr.len and marker.int are given in centiMorgan.
  magic.rhap.marker <- function(xinfo, self=c(0,0,4), f.cor, chr.len=200, marker.int=0.05){
    
    # limit to 8 founders.
    n <- 8
    
    # calculate the number of markers given chr.len and marker.int.
    n.marker <- floor(chr.len/marker.int)
    
    # create the marker positions.
    gmap <- seq(0, chr.len, marker.int)
    gmap <- gmap[1:n.marker]
    names(gmap) <- paste("SITE", formatC(1:n.marker, width=nchar(n.marker), format="d", flag="0"), sep="_")
    
    # create the founder marker data.
    geno.founder <- bindata::rmvbin(n=n.marker,
                                    margprob=rep(0.5, n),
                                    bincorr=f.cor)
    geno.founder <- t(geno.founder*2)
    rownames(geno.founder) <- LETTERS[1:n]
    colnames(geno.founder) <- paste("SITE", formatC(1:n.marker, width=nchar(n.marker), format="d", flag="0"), sep="_")
    
    # get the number of crossing generation (nx).
    nx <- log(n, 2)
    
    # convert xinfo into xplan.
    xplan <- list(t(matrix(c(t(xinfo)), nrow=2)))
    for(i in 2:nx) xplan[[i]] <- matrix(1:nrow(xplan[[i-1]]), ncol=2, byrow=T)
    
    # add selfing generations.
    xplan2 <- list()
    for(i in 1:nx){
      if(self[i] > 0){
        xplan2 <- c(xplan2, xplan[i], replicate(self[i], list(cbind(1:nrow(xplan[[i]]), 1:nrow(xplan[[i]])))))
      } else {
        xplan2 <- c(xplan2, xplan[i])
      }
    }
    xplan <- xplan2
    
    # get the number of generations.
    n.gen <- length(xplan)
    
    # create founder haplotypes for use in AlphaSimR.
    n.ind <- nrow(xplan[[length(xplan)]])
    gen.map <- sort(rep(unname(gmap)/100, n-1))
    fhap <- matrix(rep(rbind(0, diag(n-1)), n.marker), nrow=n, ncol=n.marker*(n-1))
    
    # simulate MAGIC population.
    founder <- newMapPop(genMap=list(gen.map),
                         haplotypes=list(fhap),
                         inbred=T,
                         ploidy=2L)
    SP <- SimParam$new(founder)
    sim <- newPop(founder, simParam=SP)
    
    for(j in 1:n.gen){
      sim <- makeCross(pop=sim,
                       crossPlan=xplan[[j]],
                       nProgeny=1,
                       simParam=SP)
    }
    
    # extract the haploid marker data for all chromosomes.
    hap <- pullSegSiteHaplo(pop=sim, chr=NULL, simParam=SP)
    hap <- hap[seq(1, nrow(hap), 2), , drop=F]
    hap <- magic.hap(hap=hap, n=n)
    
    # calculate the number of rec-hap from the true simulated data.
    out <- replicate(2, vector())
    {
      # identify the marker segment (mseg) for each individual.
      mseg <- vector()
      for(k in 1:nrow(hap)){
        temp <- which(!(hap[k, -1] == hap[k, -ncol(hap)]))
        if(length(temp) > 0){
          temp <- cbind(k, c(1,temp+1), c(temp, ncol(hap)))
        } else {
          temp <- cbind(k, 1, ncol(hap))
        }
        mseg <- rbind(mseg, temp)
      }
      
      # identify the founder of each mseg.
      mseg <- cbind(mseg, 
                    sapply(1:nrow(mseg), FUN=function(x) hap[mseg[x,1], mseg[x,2]]))
      attr(mseg, "dimnames") <- NULL
      
      # remove any segment that is 0 (remove NA segment).
      mseg <- mseg[!(mseg[,4]==0),]
      
      # merge any segment that are the same (e.g. 5-0-5 becomes 5-5 previously, and now 5-5 becomes 5).
      temp1 <- which(mseg[-1, 4] == mseg[-nrow(mseg), 4] & mseg[-1, 1] == mseg[-nrow(mseg), 1])
      if(length(temp1) > 0){
        temp2 <- mseg[temp1 + 1, 3]
        mseg[temp1, 3] <- temp2
        mseg <- mseg[-c(temp1 + 1), ]
      }
      
      # calculate the number of recombinant haplotype.
      rhap <- cbind(sort(rep(1:n, n)), rep(1:n,n))
      rhap <- rhap[!(rhap[,1]==rhap[,2]), ]
      rseg <- matrix(0, nrow=nrow(hap), ncol=nrow(rhap))
      for(k in 1:nrow(rhap)){
        temp <- mseg[mseg[-nrow(mseg),4]==rhap[k,1] & mseg[-1,4]==rhap[k,2] & mseg[-nrow(mseg),1]==mseg[-1,1], 1]
        temp1 <- unique(temp)
        temp2 <- unname(c(table(temp)))
        rseg[temp1, k] <- temp2
      }
      rseg <- t(rseg)
      rseg <- rowSums(rseg)/ncol(rseg)
      
      # add the results to out.
      out[[1]] <- rseg
    }
    
    
    # convert hap into geno for the RILs.
    geno <- matrix(0, nrow=nrow(hap), ncol=ncol(hap))
    for(j in 1:n){
      temp <- t(matrix(rep(geno.founder[j,], nrow(hap)), nrow=ncol(hap)))
      temp[!(hap == j)] <- 0
      geno <- geno + temp
    }
    rownames(geno) <- rownames(xinfo)
    colnames(geno) <- colnames(geno.founder)
    
    # create a crossfile.
    crossfile <- list(crosstype="riself8",
                      geno=list(geno + 1),
                      gmap=list(gmap),
                      founder_geno=list(geno.founder + 1),
                      is_x_chr=FALSE,
                      is_female=rep(FALSE, nrow(xinfo)),
                      cross_info=xinfo,
                      alleles=LETTERS[1:n])
    class(crossfile) <- "cross2"
    
    # create a list of marker subsets.
    temp <- marker.int
    while(max(temp) < 10) temp <- c(temp, 2*max(temp))
    idx <- lapply(temp, FUN=function(x) which(gmap%in%seq(0, chr.len, x)) )

    # loop to calculate the fcall and identify rhap for each data subset.
    for(i in 1:length(idx)){
      
      # subset the crossfile to smaller sets of markers.
      temp <- crossfile
      temp$geno[[1]] <- temp$geno[[1]][, idx[[i]] ]
      temp$gmap[[1]] <- temp$gmap[[1]][ idx[[i]] ]
      temp$founder_geno[[1]] <- temp$founder_geno[[1]][, idx[[i]] ]
      
      # calculate genotype probabilities.
      gp <- calc_genoprob(cross=temp, error_prob=0.01, map_function="haldane", cores=0)
      
      # call the founder at minprob of 0.5001.
      fcall <- maxmarg(probs=gp, minprob=0.5001, cores=0)[[1]]
      fcall[is.na(fcall)] <- 0
      attr(fcall, "dimnames") <- NULL
      
      # identify the marker segment (mseg) for each individual.
      mseg <- vector()
      for(k in 1:nrow(fcall)){
        temp <- which(!(fcall[k, -1] == fcall[k, -ncol(fcall)]))
        if(length(temp) > 0){
          temp <- cbind(k, c(1,temp+1), c(temp, ncol(fcall)))
        } else {
          temp <- cbind(k, 1, ncol(fcall))
        }
        mseg <- rbind(mseg, temp)
      }
      
      # identify the founder of each mseg.
      mseg <- cbind(mseg, 
                    sapply(1:nrow(mseg), FUN=function(x) fcall[mseg[x,1], mseg[x,2]]))
      attr(mseg, "dimnames") <- NULL
      
      # remove any segment that is 0 (remove NA segment).
      mseg <- mseg[!(mseg[,4]==0),]
      
      # merge any segment that are the same (e.g. 5-0-5 becomes 5-5 previously, and now 5-5 becomes 5).
      temp1 <- which(mseg[-1, 4] == mseg[-nrow(mseg), 4] & mseg[-1, 1] == mseg[-nrow(mseg), 1])
      if(length(temp1) > 0){
        temp2 <- mseg[temp1 + 1, 3]
        mseg[temp1, 3] <- temp2
        mseg <- mseg[-c(temp1 + 1), ]
      }
      
      # calculate the number of recombinant haplotype.
      rhap <- cbind(sort(rep(1:n, n)), rep(1:n,n))
      rhap <- rhap[!(rhap[,1]==rhap[,2]), ]
      rseg <- matrix(0, nrow=nrow(fcall), ncol=nrow(rhap))
      for(k in 1:nrow(rhap)){
        temp <- mseg[mseg[-nrow(mseg),4]==rhap[k,1] & mseg[-1,4]==rhap[k,2] & mseg[-nrow(mseg),1]==mseg[-1,1], 1]
        temp1 <- unique(temp)
        temp2 <- unname(c(table(temp)))
        rseg[temp1, k] <- temp2
      }
      rseg <- t(rseg)
      rseg <- rowSums(rseg)/ncol(rseg)
      
      # add the results to out.
      out[[2]] <- c(out[[2]], list(rseg))

    }
    
    return(out)
    
  }
  
  # function to simulate and compare the error rates at different minprob.
  magic.fminprob <- function(crossfile, self){
    
    # read the marker files using qtl2.
    crossfile <- read_cross2(crossfile)
    
    # get gen.map.
    gen.map <- crossfile$gmap
    
    # get cross info (xinfo).
    xinfo <- crossfile$cross_info
    
    # get the number of founder (n) and crossing generation (nx).
    n <- ncol(crossfile$cross_info)
    nx <- log(n, 2)
    
    # convert xinfo into xplan.
    xplan <- list(t(matrix(c(t(xinfo)), nrow=2)))
    for(i in 2:nx) xplan[[i]] <- matrix(1:nrow(xplan[[i-1]]), ncol=2, byrow=T)
    
    # add selfing generations.
    xplan2 <- list()
    for(i in 1:nx){
      if(self[i] > 0){
        xplan2 <- c(xplan2, xplan[i], replicate(self[i], list(cbind(1:nrow(xplan[[i]]), 1:nrow(xplan[[i]])))))
      } else {
        xplan2 <- c(xplan2, xplan[i])
      }
    }
    xplan <- xplan2
    
    # get the number of generations.
    n.gen <- length(xplan)
    
    # get the number of chromosomes.
    n.chr <- length(gen.map)
    
    # create a vector to indicate which chromosome a marker belongs to.
    chr <- vector()
    for(i in 1:n.chr) chr <- c(chr, rep(i, length(gen.map[[i]])))
    
    # create founder haplotypes for use in AlphaSimR.
    n.ind <- nrow(xplan[[length(xplan)]])
    n.marker <- sapply(1:n.chr, FUN=function(x) length(gen.map[[x]]))
    gen.map <- lapply(1:n.chr, FUN=function(x) (gen.map[[x]] - gen.map[[x]][1])/100)
    gen.map <- lapply(1:n.chr, FUN=function(x) sort(rep(gen.map[[x]], n-1)))
    fhap <- lapply(1:n.chr, FUN=function(x) matrix(rep(rbind(0, diag(n-1)), n.marker[x]), nrow=n, ncol=n.marker[x]*(n-1)))
    
    # simulate MAGIC population.
    founder <- newMapPop(genMap=gen.map,
                         haplotypes=fhap,
                         inbred=T,
                         ploidy=2L)
    SP <- SimParam$new(founder)
    sim <- newPop(founder, simParam=SP)
    
    for(j in 1:n.gen){
      sim <- makeCross(pop=sim,
                       crossPlan=xplan[[j]],
                       nProgeny=1,
                       simParam=SP)
    }
    
    # extract the haploid marker data for all chromosomes.
    hap <- pullSegSiteHaplo(pop=sim, chr=NULL, simParam=SP)
    hap <- hap[seq(1, nrow(hap), 2), , drop=F]
    hap <- magic.hap(hap=hap, n=n)
    
    # convert hap into geno for the RILs.
    fgeno <- do.call(cbind, crossfile$founder_geno)
    geno <- matrix(0, nrow=nrow(hap), ncol=ncol(hap))
    for(j in 1:n){
      temp <- t(matrix(rep(fgeno[j,], nrow(hap)), nrow=ncol(hap)))
      temp[!(hap == j)] <- 0
      geno <- geno + temp
    }
    
    # create "cross2" class object for qtl2.
    dat <- crossfile
    dat$geno <- lapply(1:n.chr, FUN=function(x) geno[, chr==x, drop=F])
    for(j in 1:n.chr){
      rownames(dat$geno[[j]]) <- rownames(crossfile$geno[[j]])
      colnames(dat$geno[[j]]) <- colnames(crossfile$geno[[j]])
    }
    
    # calculate genotype probabilities.
    gp <- calc_genoprob(cross=dat, error_prob=0.01, map_function="haldane", cores=0)
    
    # minprob for testing.
    mp <- seq(0.1, 1, 0.1)
    
    # call the founders at different minprob.
    out <- vector()
    for(j in 1:10){
      
      fcall <- maxmarg(probs=gp, minprob=mp[j], cores=0)
      fcall <- do.call(cbind, fcall)
      
      fout <- fcall == hap
      fout <- c(correct=sum(fout, na.rm=T), incorrect=sum(!fout, na.rm=T), missing=sum(is.na(fout)))
      out <- rbind(out, fout/sum(fout))
      
    }
    
    return(out)
    
  }
  
}

### Figure 1 - Distribution of currently available MAGIC population designs.
{
  dat <- read.csv("input/magic_pop_designs.csv", as.is=T)
  dat$n <- as.factor(dat$n)
  dat$design <- as.factor(dat$design)
  dat$design <- factor(dat$design, labels=c("Basic", "Partial", "Semi-structured", "Unstructured"))
  colnames(dat)[2] <- "Design"
  
  ggplot() +
    geom_bar(data=dat, aes(x=n, y=count, fill=Design), stat="identity") +
    theme(panel.background=element_blank(), panel.grid=element_blank()) +
    scale_y_continuous(breaks=0:26, labels=c(0,rep("",9),10,rep("",9),20,rep("",6)), expand=c(0,0)) +
    annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#AAAAAA") +
    theme(axis.ticks.y=element_line(color=c("#888888",rep("#CCCCCC",9),"#888888",rep("#CCCCCC",9),"#888888",rep("#CCCCCC",6)))) +
    theme(axis.ticks.y=element_line(size=0.25)) +
    theme(axis.ticks.x=element_line(color="#888888", size=0.25)) +
    theme(legend.position=c(0.8,0.8)) +
    theme(legend.title=element_text(size=10), legend.text=element_text(size=8)) +
    guides(fill=guide_legend(override.aes=list(size=0.25))) +
    xlab("Number of founders") +
    ylab("Number of MAGIC populations")
  
  ggsave(filename="output/Figure_1.svg",
         height=3.5,
         width=3,
         scale=4/3,
         units="in")
  
}

### Figure 2A, 2B, 4, S1, S2, S4 - Comparison of recombinant haplotypes from different designs.
{
  ## Part 1 - Data preparation.
  {
    
    # wheat UK elite data 1 (Mackay et al 2014, map from Gardner et al 2016).
    dat.UK <- read.csv("input/wheat_uk_elite.csv", as.is=T)
    dat.UK <- dat.UK[!(dat.UK$Marker=="Kukri_c60913_155"), ] # remove this duplicated marker.
    format.qtl2(map=dat.UK[, 1:3],
                geno.founder=as.matrix(dat.UK[, 4:11]),
                geno.ril=as.matrix(dat.UK[, 12:654]),
                filename="wheat_uk_elite_1")
    saveRDS(dat.UK, "input/wheat_uk_elite_1.RDS")
    
    # wheat DE elite data 1 (Sannemann et al 2018, map from Wang et al 2014).
    dat.DE <- read.csv("input/wheat_de_elite.csv", as.is=T)
    temp <- as.matrix(dat.DE[, 4:921]) # set imputed marker as NA.
    temp[!(temp%in%c(0,1,2))] <- NA
    dat.DE <- data.frame(dat.DE[,1:3], temp)
    format.qtl2(map=dat.DE[, 1:3],
                geno.founder=as.matrix(dat.DE[, 4:11]),
                geno.ril=as.matrix(dat.DE[, 12:921]),
                filename="wheat_de_elite_1")
    saveRDS(dat.DE, "input/wheat_de_elite_1.RDS")
    
    # wheat UK elite data 2 (map from Gardner et al 2016, common markers with DE).
    # keep only markers where the marker names and chromosomes match.
    temp <- merge(x=dat.DE[,1:3], y=dat.UK[,1:3], by=c("Marker", "Chr"), all.x=T, all.y=F, sort=F)
    temp <- temp$Marker[!is.na(temp$Pos.y)]
    temp <- dat.UK$Marker[dat.UK$Marker%in%temp]
    dat.UK2 <- dat.UK
    rownames(dat.UK2) <- dat.UK2$Marker
    dat.UK2 <- dat.UK2[temp, ]
    rownames(dat.UK2) <- NULL
    format.qtl2(map=dat.UK2[, 1:3],
                geno.founder=as.matrix(dat.UK2[, 4:11]),
                geno.ril=as.matrix(dat.UK2[, 12:654]),
                filename="wheat_uk_elite_2")
    saveRDS(dat.UK2, "input/wheat_uk_elite_2.RDS")
    
    # wheat DE elite data 2 (map from Gardner et al 2016, common markers with UK).
    # make sure the markers are in the same exact order as dat.UK2.
    dat.DE2 <- dat.DE
    rownames(dat.DE2) <- dat.DE2$Marker
    dat.DE2 <- dat.DE2[temp, ]
    rownames(dat.DE2) <- NULL
    dat.DE2$Pos <- dat.UK2$Pos
    format.qtl2(map=dat.DE2[, 1:3],
                geno.founder=as.matrix(dat.DE2[, 4:11]),
                geno.ril=as.matrix(dat.DE2[, 12:921]),
                filename="wheat_de_elite_2")
    saveRDS(dat.DE2, "input/wheat_de_elite_2.RDS")
    
    # cowpea data.
    dat.CP <- read.csv("input/cowpea.csv", as.is=T)
    temp1 <- cbind(substr(dat.CP$Allele, 1, 1), substr(dat.CP$Allele, 3, 3))
    temp1 <- cbind(paste(temp1[,1], temp1[,1], sep=""),
                   paste(temp1[,1], temp1[,2], sep=""),
                   paste(temp1[,2], temp1[,1], sep=""),
                   paste(temp1[,2], temp1[,2], sep=""))
    temp2 <- as.matrix(dat.CP[,5:317])
    temp2[temp2==temp1[,1]] <- 0
    temp2[temp2==temp1[,2]] <- 1
    temp2[temp2==temp1[,3]] <- 1
    temp2[temp2==temp1[,4]] <- 2
    temp2[temp2=="--"] <- NA
    temp2 <- matrix(as.numeric(temp2), nrow=nrow(temp2), ncol=ncol(temp2))
    colnames(temp2) <- colnames(dat.CP)[5:317]
    dat.CP <- data.frame(dat.CP[,1:3], temp2)
    dat.CP <- dat.CP[rowSums(is.na(as.matrix(dat.CP[,4:11])))==0, ]
    format.qtl2(map=dat.CP[,1:3],
                geno.founder=dat.CP[,4:11],
                geno.ril=dat.CP[,12:316],
                filename="cowpea")
    saveRDS(dat.CP, "input/cowpea.RDS")
    
    # tomato data.
    dat.TM <- read.csv("input/tomato.csv", as.is=T)
    temp1 <- cbind(dat.TM$AlleleX,
                   paste(dat.TM$AlleleX, dat.TM$AlleleY, sep="/"),
                   paste(dat.TM$AlleleY, dat.TM$AlleleX, sep="/"),
                   dat.TM$AlleleY)
    temp2 <- as.matrix(dat.TM[,6:251])
    temp2[temp2==temp1[,1]] <- 0
    temp2[temp2==temp1[,2]] <- 1
    temp2[temp2==temp1[,3]] <- 1
    temp2[temp2==temp1[,4]] <- 2
    temp2[temp2=="-"] <- NA
    temp2 <- matrix(as.numeric(temp2), nrow=nrow(temp2), ncol=ncol(temp2))
    colnames(temp2) <- colnames(dat.TM)[6:251]
    dat.TM <- data.frame(dat.TM[,1:3], temp2)
    dat.TM <- dat.TM[, c(1:3,5,4,7,6,8,9,11,10,12:249)]
    format.qtl2(map=dat.TM[,1:3],
                geno.founder=dat.TM[,4:11],
                geno.ril=dat.TM[,12:249],
                filename="tomato")
    saveRDS(dat.TM, "input/tomato.RDS")
    
  }
  
  ## Part 2 - Calling founders.
  {
    
    # some notes:
    # save the genotype probabilities in external drive (large file).
    # call the founder genotypes in the RILs with a minimum probability of 0.5001.
    
    # wheat UK elite 1.
    temp <- read_cross2("input/qtl2/wheat_uk_elite_1.zip")
    gp <- calc_genoprob(cross=temp, error_prob=0.01, map_function="haldane", cores=0)
    saveRDS(gp, "D:/gp_wheat_uk_elite_1.RDS")
    fcall <- maxmarg(probs=gp, minprob=0.5001, cores=0)
    saveRDS(fcall, "output/fcall_wheat_uk_elite_1.RDS")
    rm(temp, gp, fcall)
    
    # wheat UK elite 2.
    temp <- read_cross2("input/qtl2/wheat_uk_elite_2.zip")
    gp <- calc_genoprob(cross=temp, error_prob=0.01, map_function="haldane", cores=0)
    saveRDS(gp, "D:/gp_wheat_uk_elite_2.RDS")
    fcall <- maxmarg(probs=gp, minprob=0.5001, cores=0)
    saveRDS(fcall, "output/fcall_wheat_uk_elite_2.RDS")
    rm(temp, gp, fcall)
    
    # wheat DE elite 1.
    temp <- read_cross2("input/qtl2/wheat_de_elite_1.zip")
    gp <- calc_genoprob(cross=temp, error_prob=0.01, map_function="haldane", cores=0)
    saveRDS(gp, "D:/gp_wheat_de_elite_1.RDS")
    fcall <- maxmarg(probs=gp, minprob=0.5001, cores=0)
    saveRDS(fcall, "output/fcall_wheat_de_elite_1.RDS")
    rm(temp, gp, fcall)
    
    # wheat DE elite 2.
    temp <- read_cross2("input/qtl2/wheat_de_elite_2.zip")
    gp <- calc_genoprob(cross=temp, error_prob=0.01, map_function="haldane", cores=0)
    saveRDS(gp, "D:/gp_wheat_de_elite_2.RDS")
    fcall <- maxmarg(probs=gp, minprob=0.5001, cores=0)
    saveRDS(fcall, "output/fcall_wheat_de_elite_2.RDS")
    rm(temp, gp, fcall)
    
    # cowpea.
    temp <- read_cross2("input/qtl2/cowpea.zip")
    gp <- calc_genoprob(cross=temp, error_prob=0.01, map_function="haldane", cores=0)
    saveRDS(gp, "D:/gp_cowpea.RDS")
    fcall <- maxmarg(probs=gp, minprob=0.5001, cores=0)
    saveRDS(fcall, "output/fcall_cowpea.RDS")
    rm(temp, gp, fcall)
    
    # tomato.
    temp <- read_cross2("input/qtl2/tomato.zip")
    gp <- calc_genoprob(cross=temp, error_prob=0.01, map_function="haldane", cores=0)
    saveRDS(gp, "D:/gp_tomato.RDS")
    fcall <- maxmarg(probs=gp, minprob=0.5001, cores=0)
    saveRDS(fcall, "output/fcall_tomato.RDS")
    rm(temp, gp, fcall)
    
    # wheat UK diverse (16 founders).
    # we will call the founders with equivalent minprob of 0.5.
    # since STITCH outputs 0-2 instead of 0-1 prob, we will use cutoff of 1.
    # WARNING: THIS IS VERY SLOW.
    filename <- list.files(path="D:/MAGIC_HAPLOTYPES/", pattern="hap.txt")
    filename <- filename[-22] # ignore the unmapped markers.
    filename <- paste("D:/MAGIC_HAPLOTYPES/", filename, sep="")
    
    fcall <- replicate(21, vector())
    
    for(i in 1:21){
      mat <- t(as.matrix(read.table(file=filename[i], sep="\t", skip=1)[,-c(1:3)])) > 1
      attr(mat, "dimnames") <- NULL
      
      z1 <- seq(1, ncol(mat), 16)
      z2 <- seq(16, ncol(mat), 16)
      for(j in 1:length(z1)){
        fcall[[i]] <- cbind(fcall[[i]], mat[,z1[j]:z2[j]]%*%c(1:16))
      }
    }
    
    saveRDS(fcall, "output/fcall_wheat_uk_diverse.RDS")
    
  }
  
  ## Part 3 - Compare the simulated vs actual data for each population.
  {
    # wheat uk elite 1.
    fcall <- readRDS("output/fcall_wheat_uk_elite_1.RDS")
    wue.1 <- magic.posthoc(crossfile="input/qtl2/wheat_uk_elite_1.zip", fcall=fcall, self=c(0,0,4), n.sim=100)
    saveRDS(wue.1, "output/summary_wheat_uk_elite_1.RDS")
    ggplot_build(magic.plot.rseg(dat=wue.1, n=8))$layout$panel_scales_y[[1]] #3.69
    magic.plot.rseg(dat=wue.1, n=8, ylim=3.70)
    ggsave(filename="output/iplot_rseg_wheat_uk_elite_1.svg",
           width=7,
           height=3,
           units="in",
           scale=4/3)
    
    max(rowSums(do.call(cbind, wue.1[[1]]))/ncol(wue.1[[1]][[1]])) #3.50
    max(rowSums(do.call(cbind, wue.1[[3]]))/ncol(wue.1[[3]][[1]])) #2.15
    magic.plot.rseg3(dat=wue.1, n=8, lim=3.50)
    ggsave(filename="output/mplot_rseg_wheat_uk_elite_1.svg",
           width=3,
           height=3,
           units="in",
           scale=4/3)
    
    # wheat uk elite 2.
    fcall <- readRDS("output/fcall_wheat_uk_elite_2.RDS")
    wue.2 <- magic.posthoc(crossfile="input/qtl2/wheat_uk_elite_2.zip", fcall=fcall, self=c(0,0,4), n.sim=100)
    saveRDS(wue.2, "output/summary_wheat_uk_elite_2.RDS")
    ggplot_build(magic.plot.rseg(dat=wue.2, n=8))$layout$panel_scales_y[[1]] #3.31
    magic.plot.rseg(dat=wue.2, n=8, ylim=5.30)
    ggsave(filename="output/iplot_rseg_wheat_uk_elite_2.svg",
           width=7,
           height=3,
           units="in",
           scale=4/3)
    
    max(rowSums(do.call(cbind, wue.2[[1]]))/ncol(wue.2[[1]][[1]])) #3.09
    max(rowSums(do.call(cbind, wue.2[[3]]))/ncol(wue.2[[3]][[1]])) #1.52
    magic.plot.rseg3(dat=wue.2, n=8, lim=5.07)
    ggsave(filename="output/mplot_rseg_wheat_uk_elite_2.svg",
           width=3,
           height=3,
           units="in",
           scale=4/3)
    
    # wheat de elite 1.
    fcall <- readRDS("output/fcall_wheat_de_elite_1.RDS")
    wde.1 <- magic.posthoc(crossfile="input/qtl2/wheat_de_elite_1.zip", fcall=fcall, self=c(0,0,4), n.sim=100)
    saveRDS(wde.1, "output/summary_wheat_de_elite_1.RDS")
    ggplot_build(magic.plot.rseg(dat=wde.1, n=8))$layout$panel_scales_y[[1]] #3.54
    magic.plot.rseg(dat=wde.1, n=8, ylim=3.70)
    ggsave(filename="output/iplot_rseg_wheat_de_elite_1.svg",
           width=7,
           height=3,
           units="in",
           scale=4/3)
    
    max(rowSums(do.call(cbind, wde.1[[1]]))/ncol(wde.1[[1]][[1]])) #3.37
    max(rowSums(do.call(cbind, wde.1[[3]]))/ncol(wde.1[[3]][[1]])) #3.07
    magic.plot.rseg3(dat=wde.1, n=8, lim=3.50)
    ggsave(filename="output/mplot_rseg_wheat_de_elite_1.svg",
           width=3,
           height=3,
           units="in",
           scale=4/3)
    
    
    # wheat de elite 2.
    fcall <- readRDS("output/fcall_wheat_de_elite_2.RDS")
    wde.2 <- magic.posthoc(crossfile="input/qtl2/wheat_de_elite_2.zip", fcall=fcall, self=c(0,0,4), n.sim=100)
    saveRDS(wde.2, "output/summary_wheat_de_elite_2.RDS")
    ggplot_build(magic.plot.rseg(dat=wde.2, n=8))$layout$panel_scales_y[[1]] #5.29
    magic.plot.rseg(dat=wde.2, n=8, ylim=5.30)
    ggsave(filename="output/iplot_rseg_wheat_de_elite_2.svg",
           width=7,
           height=3,
           units="in",
           scale=4/3)
    
    max(rowSums(do.call(cbind, wde.2[[1]]))/ncol(wde.2[[1]][[1]])) #5.07
    max(rowSums(do.call(cbind, wde.2[[3]]))/ncol(wde.2[[3]][[1]])) #2.34
    magic.plot.rseg3(dat=wde.2, n=8, lim=5.07)
    ggsave(filename="output/mplot_rseg_wheat_de_elite_2.svg",
           width=3,
           height=3,
           units="in",
           scale=4/3)
    
    # cowpea.
    fcall <- readRDS("output/fcall_cowpea.RDS")
    cp <- magic.posthoc(crossfile="input/qtl2/cowpea.zip", fcall=fcall, self=c(0,0,7), n.sim=100)
    saveRDS(cp, "output/summary_cowpea.RDS")
    magic.plot.rseg(dat=cp, n=8)
    ggsave(filename="output/iplot_rseg_cowpea.svg",
           width=7,
           height=3,
           units="in",
           scale=4/3)
    magic.plot.rseg3(dat=cp, n=8)
    ggsave(filename="output/mplot_rseg_cowpea.svg",
           width=3,
           height=3,
           units="in",
           scale=4/3)
    
    # tomato.
    fcall <- readRDS("output/fcall_tomato.RDS")
    tm <- magic.posthoc(crossfile="input/qtl2/tomato.zip", fcall=fcall, self=c(0,0,3), n.sim=100)
    saveRDS(tm, "output/summary_tomato.RDS")
    magic.plot.rseg(dat=tm, n=8)
    ggsave(filename="output/iplot_rseg_tomato.svg",
           width=7,
           height=3,
           units="in",
           scale=4/3)
    magic.plot.rseg3(dat=tm, n=8)
    ggsave(filename="output/mplot_rseg_tomato.svg",
           width=3,
           height=3,
           units="in",
           scale=4/3)
    
    # wheat uk diverse.
    # we do not have genetic map, so we will simulate the best we can.
    ped <- read.csv("input/ped_wheat_uk_diverse.csv", as.is=T)
    chr.len <- read.csv("input/wheat_uk_elite.csv", as.is=T)
    chr.len <- round(sapply(unique(chr.len$Chr), FUN=function(x) max(chr.len$Pos[chr.len$Chr==x]) - min(chr.len$Pos[chr.len$Chr==x]) ), 1)/100
    
    out.rseg <- replicate(21, vector())
    for(i in 1:100){
      
      # get the founder haplotypes from a simulated population.
      hap <- magic.eval(ped=ped, marker.dist=0.005, chr.len=chr.len, n.sim=1, n.hap=1, keep=T)$hap
      temp <- paste("chr", 1:21, "_", sep="")
      hap <- lapply(temp, FUN=function(x) hap[,grepl(x, colnames(hap))])
      
      # keep the proportion of each founder in each simulation.
      for(j in 1:21){
        
        # identify the marker segment (mseg) for each individual.
        mseg <- vector()
        for(k in 1:nrow(hap[[j]])){
          temp <- which(!(hap[[j]][k, -1] == hap[[j]][k, -ncol(hap[[j]])]))
          if(length(temp) > 0){
            temp <- cbind(k, c(1,temp+1), c(temp, ncol(hap[[j]])))
          } else {
            temp <- cbind(k, 1, ncol(hap[[j]]))
          }
          mseg <- rbind(mseg, temp)
        }
        
        # identify the founder of each mseg.
        mseg <- cbind(mseg, 
                      sapply(1:nrow(mseg), FUN=function(x) hap[[j]][mseg[x,1], mseg[x,2]]))
        attr(mseg, "dimnames") <- NULL
        
        # calculate the number of recombinant haplotype.
        rhap <- cbind(sort(rep(1:16, 16)), rep(1:16, 16))
        rhap <- rhap[!(rhap[,1]==rhap[,2]), ]
        rseg <- matrix(0, nrow=nrow(hap[[j]]), ncol=nrow(rhap))
        for(k in 1:nrow(rhap)){
          temp <- mseg[mseg[-nrow(mseg),4]==rhap[k,1] & mseg[-1,4]==rhap[k,2] & mseg[-nrow(mseg),1]==mseg[-1,1], 1]
          temp1 <- unique(temp)
          temp2 <- unname(c(table(temp)))
          rseg[temp1, k] <- temp2
        }
        rseg <- colSums(rseg)/nrow(rseg)
        out.rseg[[j]] <- cbind(out.rseg[[j]], rseg)
        
      }
      
    }
    
    for(i in 1:length(out.rseg)) attr(out.rseg[[i]], "dimnames") <- NULL
    
    fcall <- readRDS("output/fcall_wheat_uk_diverse.RDS")
    
    out.rseg3 <- list()
    for(j in 1:21){
      
      # get the fcall for j-th chr.
      hap <- fcall[[j]]
      
      # identify the marker segment (mseg) for each individual.
      mseg <- vector()
      for(k in 1:nrow(hap)){
        temp <- which(!(hap[k, -1] == hap[k, -ncol(hap)]))
        if(length(temp) > 0){
          temp <- cbind(k, c(1,temp+1), c(temp, ncol(hap)))
        } else {
          temp <- cbind(k, 1, ncol(hap))
        }
        mseg <- rbind(mseg, temp)
      }
      
      # identify the founder of each mseg.
      mseg <- cbind(mseg, 
                    sapply(1:nrow(mseg), FUN=function(x) hap[mseg[x,1], mseg[x,2]]))
      attr(mseg, "dimnames") <- NULL
      
      # remove any segment that is 0 (remove NA segment).
      mseg <- mseg[!(mseg[,4]==0),]
      
      # merge any segment that are the same (e.g. 5-0-5 becomes 5-5 previously, and now 5-5 becomes 5).
      temp1 <- which(mseg[-1, 4] == mseg[-nrow(mseg), 4] & mseg[-1, 1] == mseg[-nrow(mseg), 1])
      if(length(temp1) > 0){
        temp2 <- mseg[temp1 + 1, 3]
        mseg[temp1, 3] <- temp2
        mseg <- mseg[-c(temp1 + 1), ]
      }
      
      # calculate the number of recombinant haplotype.
      rhap <- cbind(sort(rep(1:16, 16)), rep(1:16, 16))
      rhap <- rhap[!(rhap[,1]==rhap[,2]), ]
      rseg <- matrix(0, nrow=nrow(hap), ncol=nrow(rhap))
      for(k in 1:nrow(rhap)){
        temp <- mseg[mseg[-nrow(mseg),4]==rhap[k,1] & mseg[-1,4]==rhap[k,2] & mseg[-nrow(mseg),1]==mseg[-1,1], 1]
        temp1 <- unique(temp)
        temp2 <- unname(c(table(temp)))
        rseg[temp1, k] <- temp2
      }
      rseg <- t(rseg)
      out.rseg3 <- c(out.rseg3, list(rseg))
      
    }
    
    wud <- c(list(out.rseg), list(vector()), list(out.rseg3))
    saveRDS(wud, "output/summary_wheat_uk_diverse.RDS")
    
    magic.plot.rseg2(dat=wud, n=16)
    ggsave(filename="output/iplot_rseg_wheat_uk_diverse.svg",
           width=7,
           height=10,
           units="in",
           scale=4/3)
    magic.plot.rseg3(dat=wud, n=16)
    ggsave(filename="output/mplot_rseg_wheat_uk_diverse.svg",
           width=6,
           height=6,
           units="in",
           scale=4/3)
    
  }
  
  ## Misc
  {
    temp <- cbind(sort(rep(1:8,8)), rep(1:8,8))
    temp <- temp[!(temp[,1]==temp[,2]), ]
    temp <- paste(temp[,1], temp[,2], sep="")
    temp <- which(temp%in%c("12","21","34","43","56","65","78","87"))
    
    cp <- readRDS("output/summary_cowpea.RDS")
    a <- cp[[3]][[1]]
    for(i in 2:11) a <- a + cp[[3]][[i]]
    a <- rowSums(a)/ncol(a)
    mean(a) # 0.463
    sd(a) # 0.234
    mean(a[temp]) # 0.936
    sd(a[temp]) # 0.171
    mean(a[-temp]) # 0.384
    sd(a[-temp]) # 0.123
    
    tm <- readRDS("output/summary_tomato.RDS")
    a <- tm[[3]][[1]]
    for(i in 2:12) a <- a + tm[[3]][[i]]
    a <- rowSums(a)/ncol(a)
    mean(a) # 0.486
    sd(a) # 0.207
    mean(a[temp]) # 0.907
    sd(a[temp]) # 0.095
    mean(a[-temp]) # 0.416
    sd(a[-temp]) # 0.116
    
    wud <- readRDS("output/summary_wheat_uk_diverse.RDS")
    a <- wud[[3]][[1]]
    for(i in 2:21) a <- a + wud[[3]][[i]]
    a <- rowSums(a)/ncol(a)
    mean(a) # 0.878
    sd(a) # 0.102
  }
  
}

### Figure 2C, 2D - Histogram of counts of recombinant haplotypes for UK 2 and DE 2 datasets.
{
  # wheat UK elite 2.
  wue.2 <- readRDS("output/summary_wheat_uk_elite_2.RDS")
  a <- wue.2[[3]][[1]]
  for(i in 2:21) a <- a + wue.2[[3]][[i]]
  a <- rowSums(a)/ncol(a)
  mean(a) # 0.879
  sd(a) # 0.227
  
  wde.2 <- readRDS("output/summary_wheat_de_elite_2.RDS")
  b <- wde.2[[3]][[1]]
  for(i in 2:21) b <- b + wde.2[[3]][[i]]
  b <- rowSums(b)/ncol(b)
  mean(b) # 0.997
  sd(b) # 0.416
  
  temp <- cbind(sort(rep(1:8, 8)), rep(1:8, 8))
  temp <- temp[!(temp[,1]==temp[,2]), ]
  temp <- paste(temp[,1], temp[,2], sep="")
  temp <- which(temp%in%c("12","21","34","43","56","65","78","87"))
  mean(b[temp]) # 1.910
  sd(b[temp]) # 0.313
  mean(b[-temp]) # 0.845
  sd(b[-temp]) # 0.150
  
  temp <- c("0", rep("", 4), "0.5", rep("", 4), "1.0", rep("", 4), "1.5", rep("", 4), "2.0", rep("", 4), "2.5")
    
  ggplot() +
    annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
    geom_histogram(aes(x=a), binwidth=0.1, fill="#CCCCCC", color="#555555", lwd=0.1) +
    theme(panel.background=element_blank(), panel.grid=element_blank()) +
    xlab("mean count of recombinant haplotype per RIL") +
    scale_x_continuous(expand=c(0.02,0), limits=c(0, 2.55), breaks=seq(0, 2.5, 0.1), labels=temp) +
    scale_y_continuous(expand=c(0.02,0), limits=c(0, 18)) +
    theme(axis.title=element_text(size=9), axis.text=element_text(size=8)) +
    theme(axis.ticks=element_line(size=0.2))
  
  ggsave(filename="output/hist_rseg_uk_2.svg",
         width=3,
         height=2,
         units="in",
         scale=4/3)
  
  ggplot() +
    annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
    geom_histogram(aes(x=b), binwidth=0.1, fill="#CCCCCC", color="#555555", lwd=0.1) +
    theme(panel.background=element_blank(), panel.grid=element_blank()) +
    xlab("mean count of recombinant haplotype per RIL") +
    scale_x_continuous(expand=c(0.02,0), limits=c(0, 2.55), breaks=seq(0, 2.5, 0.1), labels=temp) +
    scale_y_continuous(expand=c(0.02,0), limits=c(0, 18)) +
    theme(axis.title=element_text(size=9), axis.text=element_text(size=8)) +
    theme(axis.ticks=element_line(size=0.2))
  
  ggsave(filename="output/hist_rseg_de_2.svg",
         width=3,
         height=2,
         units="in",
         scale=4/3)
}

### Table 2, S3, S4 - Total recombination breakpoints in each chromosome/all.
{
  # wheat uk elite 1.
  g <- readRDS("input/wheat_uk_elite_1.RDS")
  g <- sapply(unique(g$Chr), FUN=function(x) max(g$Pos[g$Chr==x]) - min(g$Pos[g$Chr==x]) )/100
  out.wue.1 <- magic.totalrec(dat=wue.1, g=g)
  write.csv(out.wue.1, "output/totalrec_wheat_uk_elite_1.csv", quote=F, row.names=F)
  
  # wheat uk elite 2.
  g <- readRDS("input/wheat_uk_elite_2.RDS")
  g <- sapply(unique(g$Chr), FUN=function(x) max(g$Pos[g$Chr==x]) - min(g$Pos[g$Chr==x]) )/100
  out.wue.2 <- magic.totalrec(dat=wue.2, g=g)
  write.csv(out.wue.2, "output/totalrec_wheat_uk_elite_2.csv", quote=F, row.names=F)
  
  # wheat de elite 1.
  g <- readRDS("input/wheat_de_elite_1.RDS")
  g <- sapply(unique(g$Chr), FUN=function(x) max(g$Pos[g$Chr==x]) - min(g$Pos[g$Chr==x]) )/100
  out.wde.1 <- magic.totalrec(dat=wde.1, g=g)
  write.csv(out.wde.1, "output/totalrec_wheat_de_elite_1.csv", quote=F, row.names=F)
  
  # wheat de elite 2.
  g <- readRDS("input/wheat_de_elite_2.RDS")
  g <- sapply(unique(g$Chr), FUN=function(x) max(g$Pos[g$Chr==x]) - min(g$Pos[g$Chr==x]) )/100
  out.wde.2 <- magic.totalrec(dat=wde.2, g=g)
  write.csv(out.wde.2, "output/totalrec_wheat_de_elite_2.csv", quote=F, row.names=F)
  
  # cowpea.
  g <- readRDS("input/cowpea.RDS")
  g <- sapply(unique(g$Chr), FUN=function(x) max(g$Pos[g$Chr==x]) - min(g$Pos[g$Chr==x]) )/100
  out.cp <- magic.totalrec(dat=cp, g=g)
  write.csv(out.cp, "output/totalrec_cowpea.csv", quote=F, row.names=F)
  
  # tomato.
  g <- readRDS("input/tomato.RDS")
  g <- sapply(unique(g$Chr), FUN=function(x) max(g$Pos[g$Chr==x]) - min(g$Pos[g$Chr==x]) )/100
  out.tm <- magic.totalrec(dat=tm, g=g)
  write.csv(out.tm, "output/totalrec_tomato.csv", quote=F, row.names=F)
  
  # wheat uk diverse.
  g <- readRDS("input/wheat_uk_elite_1.RDS")
  g <- sapply(unique(g$Chr), FUN=function(x) max(g$Pos[g$Chr==x]) - min(g$Pos[g$Chr==x]) )/100
  out.wud <- magic.totalrec2(dat=wud, g=g)
  write.csv(out.wud, "output/totalrec_wheat_uk_diverse.csv", quote=F, row.names=F)
}

### Table 1 - Percent rec-hap recovery.
{
  # wheat uk elite 1.
  wue.1 <- readRDS("output/summary_wheat_uk_elite_1.RDS")
  magic.rhap.recovery(dat=wue.1)
  #  mean    sd
  # 0.380 0.094
  
  # wheat uk elite 2.
  wue.2 <- readRDS("output/summary_wheat_uk_elite_2.RDS")
  magic.rhap.recovery(dat=wue.2)
  #  mean    sd
  # 0.296 0.077
  
  # wheat de elite 1.
  wde.1 <- readRDS("output/summary_wheat_de_elite_1.RDS")
  magic.rhap.recovery(dat=wde.1)
  #  mean    sd
  # 0.737 0.117
  
  # wheat de elite 2.
  wde.2 <- readRDS("output/summary_wheat_de_elite_2.RDS")
  magic.rhap.recovery(dat=wde.2)
  #  mean    sd
  # 0.331 0.058
  
  # cowpea.
  cp <- readRDS("output/summary_cowpea.RDS")
  magic.rhap.recovery(dat=cp)
  #  mean    sd
  # 0.674 0.207
  
  # tomato.
  tm <- readRDS("output/summary_tomato.RDS")
  magic.rhap.recovery(dat=tm)
  #  mean    sd
  # 0.410 0.106
  
  # wheat uk diverse.
  wud <- readRDS("output/summary_wheat_uk_diverse.RDS")
  magic.rhap.recovery(dat=wud)
  #  mean    sd
  # 0.799 0.092
  
}

### Figure 3, S3 - Proportion of unique vs identical rec-hap.
{
  ### 1 cM as identical.
  {
    # wheat uk elite 2.
    dat <- read_cross2("input/qtl2/wheat_uk_elite_2.zip")
    fcall <- readRDS("output/fcall_wheat_uk_elite_2.RDS")
    xo <- locate_xo(geno=fcall, map=dat$gmap, cores=0)
    
    out <- replicate(21, vector())
    for(j in 1:21){
      
      # get the fcall for j-th chr and set NA to 0.
      hap <- fcall[[j]]
      attr(hap, "dimnames") <- NULL
      
      # identify the recombinant haplotype for each individual.
      idx <- which(sapply(1:length(xo[[j]]), FUN=function(x) length(xo[[j]][[x]]) > 0))
      for(k in idx){
        temp <- hap[k,]
        temp <- temp[!is.na(temp)]
        temp <- rle(temp)$values
        out[[j]] <- rbind(out[[j]], data.frame(RIL=k, hap=paste(temp[-length(temp)], temp[-1], sep="_"), pos=xo[[j]][[k]]))
      }
    }
    
    out2 <- vector()
    for(j in 1:21){
      bin <- seq(min(dat$gmap[[j]]), max(dat$gmap[[j]]), length=round(diff(range(dat$gmap[[j]]))/1)+1)
      bin[length(bin)] <- bin[length(bin)] + 1
      bin <- cbind(bin[-length(bin)], bin[-1])
      for(k in unique(out[[j]]$hap)){
        temp <- out[[j]][out[[j]]$hap==k, ]
        temp <- sapply(1:nrow(temp), FUN=function(x) which(bin[,1] <= temp$pos[x] & bin[,2] > temp$pos[x]))
        temp <- sapply(1:nrow(bin), FUN=function(x) sum(temp==x))
        out2 <- rbind(out2, data.frame(hap=k, chr=j, bin, count=temp))
      }
    }
    
    saveRDS(out2, "output/hap_1cM_actual_UK.RDS")
    
    out3 <- out2$count[out2$count > 0]
    sum(out3 == 1)/length(out3) # 0.6523
    sum(out3 > 1)/length(out3) # 0.3477
    sum(out3[out3 == 1])/sum(out3) # 0.3667 
    sum(out3[out3 > 1])/sum(out3) # 0.6333
    
    temp <- table(out3)
    dat.plot <- data.frame(n.rhap=as.numeric(names(temp)), count=unname(c(temp)), pop="UK")
    
    
    # wheat de elite 2.
    dat <- read_cross2("input/qtl2/wheat_de_elite_2.zip")
    fcall <- readRDS("output/fcall_wheat_de_elite_2.RDS")
    xo <- locate_xo(geno=fcall, map=dat$gmap, cores=0)
    
    out <- replicate(21, vector())
    for(j in 1:21){
      
      # get the fcall for j-th chr and set NA to 0.
      hap <- fcall[[j]]
      attr(hap, "dimnames") <- NULL
      
      # identify the recombinant haplotype for each individual.
      idx <- which(sapply(1:length(xo[[j]]), FUN=function(x) length(xo[[j]][[x]]) > 0))
      for(k in idx){
        temp <- hap[k,]
        temp <- temp[!is.na(temp)]
        temp <- rle(temp)$values
        out[[j]] <- rbind(out[[j]], data.frame(RIL=k, hap=paste(temp[-length(temp)], temp[-1], sep="_"), pos=xo[[j]][[k]]))
      }
    }
    
    out2 <- vector()
    for(j in 1:21){
      bin <- seq(min(dat$gmap[[j]]), max(dat$gmap[[j]]), length=round(diff(range(dat$gmap[[j]]))/1)+1)
      bin[length(bin)] <- bin[length(bin)] + 1
      bin <- cbind(bin[-length(bin)], bin[-1])
      for(k in unique(out[[j]]$hap)){
        temp <- out[[j]][out[[j]]$hap==k, ]
        temp <- sapply(1:nrow(temp), FUN=function(x) which(bin[,1] <= temp$pos[x] & bin[,2] > temp$pos[x]))
        temp <- sapply(1:nrow(bin), FUN=function(x) sum(temp==x))
        out2 <- rbind(out2, data.frame(hap=k, chr=j, bin, count=temp))
      }
    }
    
    saveRDS(out2, "output/hap_1cM_actual_DE.RDS")
    
    out3 <- out2$count[out2$count > 0]
    sum(out3 == 1)/length(out3) # 0.5423
    sum(out3 > 1)/length(out3) # 0.4577
    sum(out3[out3 == 1])/sum(out3) # 0.1883 
    sum(out3[out3 > 1])/sum(out3) # 0.8117
    
    temp <- table(out3)
    dat.plot <- rbind(dat.plot,
                      data.frame(n.rhap=as.numeric(names(temp)), count=unname(c(temp)), pop="DE"))
    dat.plot$pop <- as.factor(dat.plot$pop)
    dat.plot$pop <- factor(dat.plot$pop, levels=c("UK", "DE"))
    
    
    
    ggplot() +
      annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
      geom_point(data=dat.plot, aes(x=n.rhap, y=count), size=0.5) +
      facet_wrap(vars(pop), nrow=2, scales="free_x") +
      theme(panel.background=element_blank(), panel.grid=element_blank()) +
      scale_x_continuous(limits=c(0, 140), breaks=seq(0, 140, 20), expand=c(0.01,0)) +
      theme(strip.background=element_blank(), strip.text=element_blank()) +
      theme(axis.title=element_text(size=9), axis.text=element_text(size=8)) +
      xlab("Number of identical recombinant haplotype")
    
    ggsave(filename="output/hap_1cM_actual_count.svg",
           height=3,
           width=5,
           units="in",
           scale=4/3)
    
    dat.plot2 <- data.frame(prop=c(0.3662, 0.6338, 0.1879, 0.8121),
                            type=c("unique", "non-unique", "unique", "non-unique"),
                            pop=c("UK", "UK", "DE", "DE"))
    dat.plot2$type <- as.factor(dat.plot2$type)
    dat.plot2$type <- factor(dat.plot2$type, levels=c("unique", "non-unique"))
    dat.plot2$pop <- as.factor(dat.plot2$pop)
    dat.plot2$pop <- factor(dat.plot2$pop, levels=c("UK", "DE"))
    
    ggplot() +
      annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
      geom_col(data=dat.plot2, aes(x=pop, y=prop, fill=type), position="stack") +
      theme(panel.background=element_blank(), panel.grid=element_blank()) +
      theme(legend.position="bottom", legend.title=element_blank(), legend.text=element_text(size=8)) +
      theme(axis.title=element_text(size=9), axis.text=element_text(size=8)) +
      scale_fill_discrete(guide=guide_legend(nrow=2))
    
    ggsave(filename="output/hap_1cM_actual_prop.svg",
           height=3,
           width=1.5,
           units="in",
           scale=4/3)
  }
  
  ### 10 cM as identical.
  {
    # wheat uk elite 2.
    dat <- read_cross2("input/qtl2/wheat_uk_elite_2.zip")
    fcall <- readRDS("output/fcall_wheat_uk_elite_2.RDS")
    xo <- locate_xo(geno=fcall, map=dat$gmap, cores=0)
    
    out <- replicate(21, vector())
    for(j in 1:21){
      
      # get the fcall for j-th chr and set NA to 0.
      hap <- fcall[[j]]
      attr(hap, "dimnames") <- NULL
      
      # identify the recombinant haplotype for each individual.
      idx <- which(sapply(1:length(xo[[j]]), FUN=function(x) length(xo[[j]][[x]]) > 0))
      for(k in idx){
        temp <- hap[k,]
        temp <- temp[!is.na(temp)]
        temp <- rle(temp)$values
        out[[j]] <- rbind(out[[j]], data.frame(RIL=k, hap=paste(temp[-length(temp)], temp[-1], sep="_"), pos=xo[[j]][[k]]))
      }
    }
    
    out2 <- vector()
    for(j in 1:21){
      bin <- seq(min(dat$gmap[[j]]), max(dat$gmap[[j]]), length=round(diff(range(dat$gmap[[j]]))/10)+1)
      bin[length(bin)] <- bin[length(bin)] + 1
      bin <- cbind(bin[-length(bin)], bin[-1])
      for(k in unique(out[[j]]$hap)){
        temp <- out[[j]][out[[j]]$hap==k, ]
        temp <- sapply(1:nrow(temp), FUN=function(x) which(bin[,1] <= temp$pos[x] & bin[,2] > temp$pos[x]))
        temp <- sapply(1:nrow(bin), FUN=function(x) sum(temp==x))
        out2 <- rbind(out2, data.frame(hap=k, chr=j, bin, count=temp))
      }
    }
    
    saveRDS(out2, "output/hap_10cM_actual_UK.RDS")
    
    out3 <- out2$count[out2$count > 0]
    sum(out3 == 1)/length(out3) # 0.3781
    sum(out3 > 1)/length(out3) # 0.6219
    sum(out3[out3 == 1])/sum(out3) # 0.1373 
    sum(out3[out3 > 1])/sum(out3) # 0.8627
    
    temp <- table(out3)
    dat.plot <- data.frame(n.rhap=as.numeric(names(temp)), count=unname(c(temp)), pop="UK")
    
    
    # wheat de elite 2.
    dat <- read_cross2("input/qtl2/wheat_de_elite_2.zip")
    fcall <- readRDS("output/fcall_wheat_de_elite_2.RDS")
    xo <- locate_xo(geno=fcall, map=dat$gmap, cores=0)
    
    out <- replicate(21, vector())
    for(j in 1:21){
      
      # get the fcall for j-th chr and set NA to 0.
      hap <- fcall[[j]]
      attr(hap, "dimnames") <- NULL
      
      # identify the recombinant haplotype for each individual.
      idx <- which(sapply(1:length(xo[[j]]), FUN=function(x) length(xo[[j]][[x]]) > 0))
      for(k in idx){
        temp <- hap[k,]
        temp <- temp[!is.na(temp)]
        temp <- rle(temp)$values
        out[[j]] <- rbind(out[[j]], data.frame(RIL=k, hap=paste(temp[-length(temp)], temp[-1], sep="_"), pos=xo[[j]][[k]]))
      }
    }
    
    out2 <- vector()
    for(j in 1:21){
      bin <- seq(min(dat$gmap[[j]]), max(dat$gmap[[j]]), length=round(diff(range(dat$gmap[[j]]))/10)+1)
      bin[length(bin)] <- bin[length(bin)] + 1
      bin <- cbind(bin[-length(bin)], bin[-1])
      for(k in unique(out[[j]]$hap)){
        temp <- out[[j]][out[[j]]$hap==k, ]
        temp <- sapply(1:nrow(temp), FUN=function(x) which(bin[,1] <= temp$pos[x] & bin[,2] > temp$pos[x]))
        temp <- sapply(1:nrow(bin), FUN=function(x) sum(temp==x))
        out2 <- rbind(out2, data.frame(hap=k, chr=j, bin, count=temp))
      }
    }
    
    saveRDS(out2, "output/hap_10cM_actual_DE.RDS")
    
    out3 <- out2$count[out2$count > 0]
    sum(out3 == 1)/length(out3) # 0.3163
    sum(out3 > 1)/length(out3) # 0.6837
    sum(out3[out3 == 1])/sum(out3) # 0.0688 
    sum(out3[out3 > 1])/sum(out3) # 0.9312
    
    temp <- table(out3)
    dat.plot <- rbind(dat.plot,
                      data.frame(n.rhap=as.numeric(names(temp)), count=unname(c(temp)), pop="DE"))
    dat.plot$pop <- as.factor(dat.plot$pop)
    dat.plot$pop <- factor(dat.plot$pop, levels=c("UK", "DE"))
    
    ggplot() +
      annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
      geom_point(data=dat.plot, aes(x=n.rhap, y=count), size=0.5) +
      facet_wrap(vars(pop), nrow=2, scales="free_x") +
      theme(panel.background=element_blank(), panel.grid=element_blank()) +
      scale_x_continuous(limits=c(0, 140), breaks=seq(0, 140, 20), expand=c(0.01,0)) +
      theme(strip.background=element_blank(), strip.text=element_blank()) +
      theme(axis.title=element_text(size=9), axis.text=element_text(size=8)) +
      xlab("Number of identical recombinant haplotype")
    
    ggsave(filename="output/hap_10cM_actual_count.svg",
           height=3,
           width=5,
           units="in",
           scale=4/3)
    
    dat.plot2 <- data.frame(prop=c(0.1370, 0.8630, 0.0686, 0.9314),
                            type=c("unique", "non-unique", "unique", "non-unique"),
                            pop=c("UK", "UK", "DE", "DE"))
    dat.plot2$type <- as.factor(dat.plot2$type)
    dat.plot2$type <- factor(dat.plot2$type, levels=c("unique", "non-unique"))
    dat.plot2$pop <- as.factor(dat.plot2$pop)
    dat.plot2$pop <- factor(dat.plot2$pop, levels=c("UK", "DE"))
    
    ggplot() +
      annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
      geom_col(data=dat.plot2, aes(x=pop, y=prop, fill=type), position="stack") +
      theme(panel.background=element_blank(), panel.grid=element_blank()) +
      theme(legend.position="bottom", legend.title=element_blank(), legend.text=element_text(size=8)) +
      theme(axis.title=element_text(size=9), axis.text=element_text(size=8)) +
      scale_fill_discrete(guide=guide_legend(nrow=2))
    
    ggsave(filename="output/hap_10cM_actual_prop.svg",
           height=3,
           width=1.5,
           units="in",
           scale=4/3)
  }
  
}

### Figure 5A, 5B - Error rates at different minprob.
{
  # compare different minprob for wheat UK data 2.
  dat1 <- magic.fminprob(crossfile="input/qtl2/wheat_uk_elite_2.zip", self=c(0,0,4))
  
  # compare different minprob for wheat DE data 2.
  dat2 <- magic.fminprob(crossfile="input/qtl2/wheat_de_elite_2.zip", self=c(0,0,4))
  
  # plot the difference.
  dat.plot <- rbind(data.frame(minprob=seq(0.1,1,0.1), dat1, Population="UK"),
                    data.frame(minprob=seq(0.1,1,0.1), dat2, Population="DE"))
  dat.plot <- melt(dat.plot, id.vars=c("minprob","Population"))
  dat.plot$Population <- factor(dat.plot$Population, levels=c("UK", "DE"))
  
  # save the data.
  saveRDS(dat.plot, "output/minprob.RDS")
  
  ggplot(data=dat.plot, aes(x=minprob, y=value, color=variable)) +
    geom_point() +
    geom_line() +
    facet_wrap(vars(Population), nrow=2) +
    theme(panel.background=element_blank(), panel.grid=element_blank()) +
    annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
    guides(color=guide_legend(title="Founder calls")) +
    ylab("proportion") +
    scale_x_continuous(breaks=seq(0,1,0.1)) +
    theme(axis.title=element_text(size=9), axis.text=element_text(size=8)) +
    theme(legend.title=element_text(size=8), legend.text=element_text(size=8))
  
  ggsave(filename="output/minprob.svg",
         height=2.5,
         width=3.5,
         units="in",
         scale=4/3)
}

### Figure 5C - Marker number vs recombinant haplotypes identified.
{
  # read the crossfile for wheat uk elite.
  xinfo <- read.csv("input/cross_info_sim.csv", row.names=1, as.is=T)
  xinfo <- as.matrix(xinfo)
  
  # get the founder marker data so we can estimate the correlations among founders.
  dat.UK <- readRDS("input/wheat_uk_elite_1.RDS")
  f.cor <- round(cor(as.matrix(dat.UK[, 4:11])) ,2)
  f.cor[f.cor < 0] <- 0
  
  # calculate the number of rec-hap recovered based on different subsets of markers.
  dat <- magic.rhap.marker(xinfo=xinfo, self=c(0,0,4), f.cor=f.cor, chr.len=200, marker.int=0.05)
  
  # save the output.
  saveRDS(dat, "output/rhap_marker.RDS")
  
  # plot the output.
  dat.plot <- data.frame(prop=colSums(sapply(1:9, FUN=function(x) dat[[2]][[x]]/dat[[1]]))/56,
                         x=c(0.05, 0.10, 0.20, 0.40, 0.80, 1.60, 3.20, 6.40, 12.80))
  
  ggplot() +
    geom_point(data=dat.plot, aes(x=x, y=prop), size=0.5) +
    theme(panel.background=element_blank(), panel.grid=element_blank()) +
    annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=NA, color="#DDDDDD") +
    ylab("proportion") +
    xlab("cM/marker") +
    theme(axis.title=element_text(size=9), axis.text=element_text(size=8)) +
    scale_x_continuous(expand=c(0.02,0), breaks=seq(0, 12, 2)) +
    scale_y_continuous(expand=c(0,0), limits=c(0,1))
  
  ggsave(filename="output/rhap_marker.svg",
         height=2,
         width=3.5,
         units="in",
         scale=4/3)
  
}

### Table 3, Figure 6, 7, Table S5 - Simulating 5 designs.
{
  # design 1 - full
  mpop.1 <- magic.eval(n=8, m=45, reps=c(1,1,3), self=c(0,0,4), balanced=T, chr.len=c(1.0, 1.5, 2, 2.5, 3.0), n.sim=100)
  
  # design 2 - partial, balanced
  mpop.2 <- magic.eval(n=8, m=1, reps=c(1,9,15), self=c(0,0,4), balanced=T, chr.len=c(1.0, 1.5, 2, 2.5, 3.0), n.sim=100)
  
  # design 3 - partial unbalanced
  mpop.3 <- magic.eval(n=8, m=7, reps=c(1,9,15), self=c(0,0,4), balanced=F, chr.len=c(1.0, 1.5, 2, 2.5, 3.0), n.sim=100)
  
  # design 4 - partial balanced + self
  mpop.4 <- magic.eval(n=8, m=1, reps=c(1,9,15), self=c(0,1,4), balanced=T, chr.len=c(1.0, 1.5, 2, 2.5, 3.0), n.sim=100)
  
  # design 5 - basic
  mpop.5 <- magic.eval(n=8, m=0, reps=c(1,4,4), self=c(0,0,0), addx=1, repx=15, selfx=4, chr.len=c(1.0, 1.5, 2, 2.5, 3.0), n.sim=100)

  # compare all designs with "interval".
  #temp <- cbind(sort(rep(1:8, 8)), rep(1:8, 8))
  #temp <- temp[!(temp[,1]==temp[,2]), ]
  temp <- cbind(c(1,1,1,3,3,3), c(2,3,5,2,4,8))
  p <- magic.plot(input=list(mpop.1, mpop.2, mpop.3, mpop.4, mpop.5), display="interval", fpair=temp)
  ggsave(plot=p,
         filename="output/sim_interval.svg",
         width=7,
         height=5,
         units="in",
         scale=4/3)
  
  # compare all designs with "whole".
  p <- magic.plot(input=list(mpop.1, mpop.2, mpop.3, mpop.4, mpop.5), display="whole", annotate=F)
  ggsave(plot=p,
         filename="output/sim_whole.svg",
         width=7,
         height=7,
         units="in",
         scale=4/3)
  
  # get summary information for all designs.
  magic.summary(input=list(mpop.1, mpop.2, mpop.3, mpop.4, mpop.5))
  #             design 1 design 2 design 3 design 4 design 5
  #founder             8        8        8        8        8
  #type             full  partial  partial  partial    basic
  #reps          1,1,3,0 1,9,15,0 1,9,15,0 1,9,15,0 1,4,4,15
  #self          0,0,4,0  0,0,4,0  0,0,4,0  0,1,4,0  0,0,0,4
  #cross      28,210,315 28,14,63 19,13,63 28,14,63 4,2,4,64
  #generation          7        7        7        8        8
  #RIL               945      945      945      945      960
  #funnel            315        7        7        7        1
  
  # save the results.
  saveRDS(mpop.1, "output/mpop1.RDS")
  saveRDS(mpop.2, "output/mpop2.RDS")
  saveRDS(mpop.3, "output/mpop3.RDS")
  saveRDS(mpop.4, "output/mpop4.RDS")
  saveRDS(mpop.5, "output/mpop5.RDS")
  
  # calculate mean and CV for total rec hap.
  temp <- colSums(mpop.1[[2]])
  mean(temp) #0.167
  sd(temp)/mean(temp) #0.075
  
  temp <- colSums(mpop.2[[2]])
  mean(temp) #0.167
  sd(temp)/mean(temp) #0.093
  
  temp <- colSums(mpop.3[[2]])
  mean(temp) #0.169
  sd(temp)/mean(temp) #0.101
  
  temp <- colSums(mpop.4[[2]])
  mean(temp) #0.186
  sd(temp)/mean(temp) #0.114
  
  temp <- colSums(mpop.5[[2]])
  mean(temp) #0.202
  sd(temp)/mean(temp) #0.271
  
  # calculate mean and CV for number of unique rec hap.
  temp <- colSums(mpop.1[[2]] > 0)
  mean(temp) #52.200
  sd(temp)/mean(temp) #0.035
  
  temp <- colSums(mpop.2[[2]] > 0)
  mean(temp) #49.930
  sd(temp)/mean(temp) #0.041
  
  temp <- colSums(mpop.3[[2]] > 0)
  mean(temp) #50.290
  sd(temp)/mean(temp) #0.041
  
  temp <- colSums(mpop.4[[2]] > 0)
  mean(temp) #47.810
  sd(temp)/mean(temp) #0.045
  
  temp <- colSums(mpop.5[[2]] > 0)
  mean(temp) #34.510
  sd(temp)/mean(temp) #0.171
}






