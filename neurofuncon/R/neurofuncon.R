neurofun2d = function(traces, position, all.stimuli.ids, nlayer, ggm.stimuli.ids,
                      num.neurons = 25, num.edges = 25, rho.seq = seq(5,0.001,-0.001),
                      Plot =TRUE, node.size = 4){
  # traces :: d x n matrix
  # position :: n x 3 data frame  >> n x 2 in future version
  # all.stimuli.ids :: vector of all stimuli ids
  # nlayer :: the layer we want to observe >> remove in future version
  # ggm.stimuli.ids :: vector of stimuli ids we want to observe
  # num.neurons :: number of nodes (neurons) to observe in the GGM
  # num.edges :: number of edges (connections) to observe in the GGM
  # rho.seq :: sequence of potential rho values
  # Plot :: boolean to create the plot
  # node.size :: integer value of the size of the neurons with the highest activation averages

  if(!require(glasso)) install.packages('glasso')
  library(glasso)
  if(!require(ggplot2)) install.packages('ggplot2')
  library(ggplot2)
  if(!require(ggmap)) install.packages('ggmap')
  library(ggmap)
  if(!require(ggpubr)) install.packages('ggpubr')
  library(ggpubr)

  # Setting up data
  stim.resp <- traces # stim.resp == stimulus response
  stim.istim <- data.frame(all.stimuli.ids)
  colnames(stim.istim) = "stim.id"

  stim.resp <- cbind(stim.resp, stim.istim) # Binding stimulus ids to stimuli response

  stim.id = ggm.stimuli.ids
  stim.idx <- which(all.stimuli.ids %in% ggm.stimuli.ids)    #time bins related to stimuli in "stim"
  layer.idx <- which(position[, 3] == nlayer)     #neurons on layer (z-dim) "nlayer"
  traces.select <- traces[stim.idx,layer.idx]     #subset of data about "stim" and "nlayer"
  layer <- position[layer.idx, ]     #positions of the selected neurons on layer (z-dim) "nlayer"

  #All activity in nlayer from stim.resp
  stim.resp.all <- stim.resp[,layer.idx]

  #All activity in nlayer from stimulus
  stim.resp.layer <- sqrt(stim.resp.all[stim.idx,])
  averages <- colMeans(stim.resp.layer)
  variances <- apply(stim.resp.layer, 2, var)

  layer <- data.frame(cbind(layer,averages,variances))


  #Using top num.neurons neurons with most activation (can change num.neurons based on how many we want)
  max.activation = head(sort(colSums(stim.resp.layer),decreasing=TRUE),n=num.neurons)
  max.activation.idx = double()
  for(i in 1:length(max.activation)){
    max.activation.idx = append(max.activation.idx,which(colSums(stim.resp.layer)==max.activation[i]))
  }

  max.activation.idx <- sort(unique(max.activation.idx))

  neurons.max.act = data.frame(matrix(ncol = 0, nrow = nrow(stim.resp.layer)))
  list.max.idx = list()

  for(idx in max.activation.idx){
    neurons.max.act = cbind(neurons.max.act, stim.resp.layer[,idx])
    title = paste0("V",idx)
    names(neurons.max.act)[ncol(neurons.max.act)] <- title
    list.max.idx <- append(list.max.idx, idx)
  }

  num.neurons <- length(max.activation.idx)

  # Making size and opacity of the points based on the top num.neurons values
  # and generating the location for the lines
  size <- double()
  alpha <- double()
  lines.loc <- data.frame(x=NULL, y=NULL)
  for (len in 1:nrow(layer)) {
    if(len %in% max.activation.idx){
      size <- append(size, layer[,"averages"][len]+0.5)
      alpha <- append(alpha, 1)
      lines.loc <- rbind(lines.loc, layer[len,1:2])
    }
    else{
      size <- append(size, 0.5)
      alpha <- append(alpha, 0.4)
    }
  }

  size <- node.size*(size/max(size))

  layer <- cbind(layer,size)
  layer <-cbind(layer,alpha)

  neurons.cov <- cov(neurons.max.act)
  neurons.cor <- cor(neurons.max.act)

  # Changing the rho based on a certain number of edges ####
  for(rho.value in rho.seq){
    theta <- glasso::glasso(neurons.cor, rho = rho.value)$wi
    edge <- theta!=0
    edge.num = (sum(edge)-num.neurons)/2
    if(edge.num >= num.edges){
      break
    }
  }

  # Determining if the edges have positive or negative correlation between neurons ####
  d = nrow(lines.loc)
  loc.p = data.frame(x1 = NA, x2 = NA, y1= NA, y2 = NA)
  loc.g = data.frame(x1 = NA, x2 = NA, y1= NA, y2 = NA)
  for(i in 1:(d-1)){
    for(j in (i+1):(d)){
      if(theta[i,j]<0){
        loc.p <- rbind(loc.p,c(lines.loc[i,1], lines.loc[j,1], lines.loc[i,2], lines.loc[j,2]))
      }
      else if(theta[i,j]>0){
        loc.g <- rbind(loc.g,c(lines.loc[i,1],lines.loc[j,1],lines.loc[i,2], lines.loc[j,2]))
      }
    }
  }
  loc.p <- loc.p[-1,]
  loc.g <- loc.g[-1,]

  num.edges <- nrow(loc.p)+nrow(loc.g)

  #If theta is negative, then partial correlation is positive & vice versa

  # Now plotting the GGM ####
  if(Plot == TRUE)
  {
    plt <- ggplot(layer, aes(layer[,1], layer[,2])) +
      geom_segment(data = loc.p, aes(x = x1, y=y1, xend=x2, yend=y2), color = "blue", linewidth = 0.75) +
      geom_segment(data = loc.g, aes(x = x1, y=y1, xend=x2, yend=y2), color = "red", linewidth = 0.75) +
      geom_point(aes(color = variances), alpha = alpha, size = size) +
      scale_color_viridis_c() +
      labs(title = paste0("Functional connectivity in layer ", nlayer, " for ", num.neurons, " neurons and ", num.edges," connections")) +
      theme(axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            text = element_text(family="Palatino"),
            plot.title = element_text(family="Palatino", face = "bold"),
            panel.background = element_rect(fill = "white",linewidth = 0.5, linetype = "solid", color="grey")) +
      annotate("text",x=max(layer[,1]),y=min(layer[,2]),label=paste0("rho = ", rho.value),family="Palatino") +
      coord_cartesian(clip="off")
    print(plt)
  }

  # Return function ####
  return(list(layer=nlayer,loction.positive=loc.p,location.negative=loc.g,layer=layer,rho=rho.value,true.num.edges=num.edges,true.num.neurons=num.neurons))
}

neurofun3d = function(traces, position, all.stimuli.ids, nlayers, ggm.stimuli.ids,
                      num.neurons = 25, num.edges = 25, rho.seq = seq(5,0.001,-0.001),
                      Plot = TRUE, node.size = 6){
  # traces :: d x n matrix
  # position :: n x 3 data frame  >> n x 2 in future version
  # all.stimuli.ids :: vector of all stimuli ids
  # nlayer :: the layer we want to observe >> remove in future version
  # ggm.stimuli.ids :: vector of stimuli ids we want to observe
  # num.neurons :: number of nodes (neurons) to observe in the GGM
  # num.edges :: number of edges (connections) to observe in the GGM
  # rho.seq :: sequence of potential rho values
  # Plot :: boolean to create the plot
  # node.size :: integer value of the size of the neurons with the highest activation averages

  if(!require(glasso)) install.packages('glasso')
  library(glasso)
  if(!require(ggplot2)) install.packages('ggplot2')
  library(ggplot2)
  if(!require(ggmap)) install.packages('ggmap')
  library(ggmap)
  if(!require(ggpubr)) install.packages('ggpubr')
  library(ggpubr)
  if(!require(colourvalues)) install.packages('colourvalues')
  library(colourvalues)
  if(!require(rgl)) install.packages('rgl')
  library(rgl)

  stim.resp <- traces
  stim.istim <- data.frame(all.stimuli.ids)
  colnames(stim.istim) = "stim.id"

  stim.resp <- cbind(stim.resp, stim.istim)

  stim.id = ggm.stimuli.ids
  stim.idx <- which(all.stimuli.ids %in% ggm.stimuli.ids)    #time bins related to stimuli in "stim"
  layer.idx <- which(position[, 3] %in% nlayers)     #neurons on layers (z-dim) "nlayers"
  traces.select <- traces[stim.idx,layer.idx]     #subset of data about "stim" and "nlayers"
  layer <- position[layer.idx, ]     #positions of the selected neurons on layer (z-dim) "nlayer"

  #All activity in nlayers from stim.resp
  stim.resp.all <- stim.resp[,layer.idx]

  #All activity in layer n from stimulus
  stim.resp.layer <- sqrt(stim.resp.all[stim.idx,])
  averages <- colMeans(stim.resp.layer)
  variances <- apply(stim.resp.layer, 2, var)
  layer <- data.frame(cbind(layer,averages,variances))


  #Using top n neurons with most activation (can change n based on how many we want)
  max.activation = head(sort(colSums(stim.resp.layer),decreasing=TRUE),n=num.neurons)
  max.activation.idx = double()
  for(i in 1:length(max.activation)){
    max.activation.idx = append(max.activation.idx,which(colSums(stim.resp.layer)==max.activation[i]))
  }

  max.activation.idx <- sort(unique(max.activation.idx))

  neurons.max.act = data.frame(matrix(ncol = 0, nrow = nrow(stim.resp.layer)))
  list.max.idx = list()

  for(idx in max.activation.idx){
    neurons.max.act = cbind(neurons.max.act, stim.resp.layer[,idx])
    title = paste0("V",idx)
    names(neurons.max.act)[ncol(neurons.max.act)] <- title
    list.max.idx <- append(list.max.idx, idx)
  }


  #Making size and opacity of the points based on the top n values
  #and generating the location for the lines
  size <- double()
  color <- double()
  lines.loc <- data.frame(x=NULL, y=NULL,z=NULL)
  for (len in 1:nrow(layer)) {
    if(len %in% max.activation.idx){
      size <- append(size, layer[,"averages"][len] + 0.5)
      color <- append(color, layer[,"variances"][len])
      lines.loc <- rbind(lines.loc, layer[len,1:3])
    }
  }

  size <- node.size*(size/max(size))

  lines.loc <- cbind(lines.loc,size)
  lines.loc<-cbind(lines.loc,color)

  neurons.cov <- cov(neurons.max.act)
  neurons.cor <- cor(neurons.max.act)

  #Changing the rho based on a certain number of edges

  for(rho.value in rho.seq){
    theta <- glasso::glasso(neurons.cor, rho = rho.value)$wi
    edge <- theta!=0
    edge.num = (sum(edge)-num.neurons)/2
    if(edge.num >= num.edges){
      break
    }
  }

  d = nrow(lines.loc)
  loc.p = data.frame(x1 = NA, y1= NA, z1=NA, x2 = NA,  y2 = NA, z2=NA)
  loc.g = data.frame(x1 = NA, y1= NA, z1=NA, x2 = NA,  y2 = NA, z2=NA)
  for(i in 1:(d-1)){
    for(j in (i+1):(d)){
      if(theta[i,j]<0){
        loc.p <- rbind(loc.p,c(lines.loc[i,1], lines.loc[i,2],lines.loc[i,3],lines.loc[j,1],lines.loc[j,2],lines.loc[j,3]))
      }
      else if(theta[i,j]>0){
        loc.g <- rbind(loc.g,c(lines.loc[i,1],lines.loc[i,2],lines.loc[i,3],lines.loc[j,1],lines.loc[j,2],lines.loc[j,3]))
      }
    }
  }
  loc.p <- loc.p[-1,]
  loc.g <- loc.g[-1,]

  if(Plot==TRUE){
    if(length(unique(layer[,3])) == 1){
      return(neurofun2d(traces,position,all.stimuli.ids,unique(layer[,3])[1],ggm.stimuli.ids,num.neurons,num.edges))
    }
    else
    {  # layers.avg <- data.frame(cbind(neuron_position, averages))
      with(layer, rgl::plot3d(layer[,1], layer[,2], layer[,3], col = colourvalues::color_values(variances), size=2,alpha=0.4,axes=FALSE,xlab="",ylab="",zlab=expression(paste("layer depth (",mu,"m)"))))
      for(lns in range(1,nrow(loc.p))){
        rgl::segments3d(x=as.vector(t(loc.p[,c(1,4)])),
                        y=as.vector(t(loc.p[,c(2,5)])),
                        z=as.vector(t(loc.p[,c(3,6)])),col="blue")
      }
      for(lns in range(1,nrow(loc.p))){
        rgl::segments3d(x=as.vector(t(loc.g[,c(1,4)])),
                        y=as.vector(t(loc.g[,c(2,5)])),
                        z=as.vector(t(loc.g[,c(3,6)])),col="red")
      }
      for(i in 1:nrow(lines.loc)) {
        rgl::points3d(x=lines.loc[i,1], y=lines.loc[i,2],z=lines.loc[i,3],size=lines.loc[i,4],color=colourvalues::color_values(lines.loc[i,5]))
      }
      rgl::axes3d(edge='z--',at = nlayers)
      rgl::box3d()
    }
  }
  return(list(nlayers=unique(layer[,3]),location.positives=loc.p,location.negative=loc.g,layers=layer,rho=rho.value,true.num.edges=edge.num))
}
