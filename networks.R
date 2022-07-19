.libPaths("C:/Program Files/R/R-4.1.2/library")
setwd('C:/Users/mihal/OneDrive')

# install.packages('igraph')
library(igraph)
m1 <- graph( ~'Lioudis'-'Ariadni'-'Michalis'-'Elewnora'-'Manos'-'Ariadni'-'Anastasia'-'Lioudis'-'Nikos'-'Michalis',
'Nikos'-'Manos'-'Anastasia'-'Nikos', 'Nikos'-'Ariadni',
'Michalis'-'Manos'-'Lioudis', 'Lioudis'-'Michalis'-'Anastasia', 'Elewnora'-'Nikos', 'Michalis'-'Katerina'-'Manos', 
'Ariadni'-'Katerina'-'Elewnora', 'Anastasia'-'Katerina'-'Lioudis', 'Nikos'-'Katerina'-'Manos',
'Michalis'-'Dimitra'-'Manos', 'Katerina'-'Dimitra'-'Nikos', 'Ariadni'-'Dimitra'-'Elewnora', 'Lioudis'-'Dimitra',
'Michalis'-'Maniou'-'Katerina', 'Lioudis'-'Maniou'-'Manos', 'Elewnora'-'Maniou',
'Michalis'-'Dorina'-'Maniou', 'Dimitra'-'Dorina'-'Manos', 'Michalis'-'Eve'-'Dorina', 'Maniou'-'Eve'-'Dimitra', 
'Michalis'-'Eve'-'Manos', 'Michalis'-'Yiannis'-'Dorina',
'Michalis'-'Marios'-'Ilias'-'Eva',  'Michalis'-'Eva'-'Marios', 'Michalis'-'Ilias',
'Lego'-'Dimitris'-'Michalis'-'Nektarios', 'Dimitris'-'Nektarios'-'Lego'-'Michalis'-'Tzou','Dimitris'-'Tzou'-'Lego',
'Marios'-'Mitsos'-'Michalis'-'Nick'-'Marios','Marios'-'Nick'-'Mitsos',
'Mariangela'-'Michalis'-'Irene'-'Mariangela', 'Panos'-'Irene'-'Andreas'-'Panos'-'Michalis'-'Andreas',
'Michalis'-'Chrysa'-'Mariangela', 'Michalis'-'Eleni'-'Lego','Nektarios'-'Eleni'-'Dimitris', 'Irene'-'Dimitris'

)
plot(m1,edge.arrow.size=1)

# Adjacency matrix
m1[]
E(m1)
V(m1)




# Add attributes to the network
V(m1)$name
V(m1)$gender <- c("male", "female", "male", "female", "male", "female", "male",'female','female','female','female','female',
'male','male','male','female','female','male','male','female','male','male','female','female','male','male','female','female')
vertex_attr(m1)

# Visualization
plot(m1, edge.arrow.size=.5, vertex.label.color="black", vertex.label.dist=1.5,
     vertex.color=c( "pink", "skyblue")[1+(V(m1)$gender=="male")] )


# Network descriptive statistics
vcount(m1)
# Number of nodes

ecount(m1)
# Number of edges

# Density: The proportion of present edges from all possible edges in thenetwork.
edge_density(m1)

# Note that it's the same as:
ecount(m1)/(vcount(m1)*(vcount(m1)-1)/2)


reciprocity(m1) 
# Always 1 for undirected graphs!


# Transitivity: ratio of triangles (direction disregarded) to connected triples
transitivity(m1, type="global")




# Diameter: A network diameter is the longest geodesic distance (length of the shortest path between two nodes) in the network.
diameter(m1)


# Node degree: The function degree() has a mode of in for in-degree, out for out-degree, and all or total for total degree.
deg <- degree(m1, mode="all")
plot(m1, vertex.size=deg*3) 



# Community detection
# Fast hierarchical agglomeration algorithm for finding community structure (Clauset, Newman, and Moore, 2004)
m <- make_graph(V(m1)$name)
fc <- cluster_fast_greedy(as.undirected(m))
plot(fc,m1)



# Visualizations
d = get.diameter(m1)
E(m1)$color = "grey"
E(m1)$width = 1
E(m1, path=d)$color = "red"
E(m1, path=d)$width = 2
V(m1)$color  = "SkyBlue2"
V(m1)[d]$color = "red"
coords = layout.fruchterman.reingold(m1)
plot(m1, layout=coords, vertex.label = NA, vertex.size=3)


hist(deg, breaks=1:vcount(m1)-1, main="Histogram of node degree")





library(network)

# install.packages('intergraph')
library(intergraph)
net.m1 <- asNetwork(m1)



# install.packages('ergm')
library(ergm)
# The ergm package fits the model using Monte Carlo maximum likelihood estimation (MCMLE).It implements maximum likelihood estimates of ERGMs to be calculated using Markov Chain Monte Carlo
# Fitting the Erdos-Renyi model
model1 <- ergm(net.m1 ~ edges)
model1
summary(model1)




# Metric of Gof for model1
control.ergm(seed = 2345)
m1gof <- gof(model1, GOF = ~ distance + espartners + triadcensus,
             verbose = TRUE, interval = 5e+4)
par(mfrow = c(2,2))
plot(m1gof, cex.lab=1.6, cex.axis=1.6, plotlogodds = TRUE)
# It's visible that the fit is really bad in this model!










# install.packages('latentnet')
library(latentnet)
# Fit latent cluster random effect models. The main idea is: Each node has a latent-unobserved position in a d-dimensional Euclidean space (e.g. d = 2). Each node is associated with one (unobserved) cluster. The model can also take into account covariates.
# Setting up the model with G=2 clusters
model.fit <- ergmm(net.m1 ~ euclidean(d = 2, G = 2), verbose = TRUE)
summary(model.fit)

# Posterior membership probabilities
attr(model.fit$sample, "Q")

par(mfrow = c(1,2))
plot(model.fit)
plot(model.fit, pie = TRUE, vertex.cex = 2.5)





# Setting up the model with G=3 clusters
model.fit.2 <- ergmm(net.m1 ~ euclidean(d = 2, G = 3), verbose = TRUE)
summary(model.fit.2)

# Posterior membership probabilities
attr(model.fit.2$sample, "Q")

par(mfrow = c(1,2))
plot(model.fit.2)
plot(model.fit.2, pie = TRUE, vertex.cex = 2.5)








# Setting up the model with G=4 clusters
model.fit.3 <- ergmm(net.m1 ~ euclidean(d = 2, G = 4), verbose = TRUE)
summary(model.fit.3)

# Posterior membership probabilities
attr(model.fit.3$sample, "Q")

par(mfrow = c(1,2))
plot(model.fit.3)
plot(model.fit.3, pie = TRUE, vertex.cex = 2.5)







