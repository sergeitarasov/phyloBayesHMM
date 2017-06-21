

set.seed(34344)
test.tree = tree.bisse(c(0.1, 0.1, 0.03, 0.03, 0.01, 0.07), x0=0, max.taxa = 50)
tree<-pbtree(n=4)
state4=c(0,1,1,0)
names(state4)=c("t2", "t3", "t4", "t1")

?tree.bisse

str(tree)
plot(tree)
plot(test.tree)
nodelabels()

two.state.part.like(test.tree, test.tree$tip.state, 0.02, 0.05)

## Arguments
M=matrix(c(-.02,.05,.02,-.05),2,2)

my.tree=tree
my.data=state4

my.tree=test.tree
my.data=test.tree$tip.state
forward.rate=0.02
backward.rate=0.05
elapsed.time=1

my.tree=tt
my.data=states
##########
str(my.tree)
plot(my.tree)
nodelabels()
str(test.tree)

## Matrix Exp
two.state.trans.prob = function(forward.rate, backward.rate, elapsed.time){
  total.rate = forward.rate + backward.rate

  return((matrix(c(rep(backward.rate,2),rep(forward.rate,2)),2,2) +
            matrix(c(forward.rate, -backward.rate, -forward.rate, backward.rate),2,2)*
            exp(-total.rate*elapsed.time))/total.rate)
}
###
str(tree)
str(my.tree )
plot(tree)
nodelabels()

tree $edge
tree $edge.length
my.tree $edge
my.tree $edge.length

test.tree$edge
##########
#two.state.part.like = function(my.tree, my.data, forward.rate, backward.rate){
  two.state.part.like = function(my.tree, my.data, M){

  ## reorder the edges in the "pruningwise" order
  my.tree = reorder(my.tree, order = "pr")

  if (!("phylo" %in% class(my.tree)))
    stop("Error: object \"my.tree\" is not of class \"phylo\"")

  if (is.null(my.tree$edge.length))
    stop("Error: tree \" my.tree\" must have branch lengths.")

  ## reorder data on tips to match the order of the my.tree phylo object
  if (!is.null(names(my.data))) {
    if(!any(is.na(match(names(my.data), my.tree$tip.label)))){
      my.data = my.data[my.tree$tip.label]
    }else{
      warning('the names of argument "my.data" and the names of the tip labels
did not match: the former were ignored in the analysis.')
    }
  }

  ## prepare transition probability matrices (this of course can and should be done in C++ as well)
 # prob.array = array(0, dim=c(2,2,length(my.tree$edge.length)))
  #for (i in 1:length(my.tree$edge.length)){
#    prob.array[,,i] = two.state.trans.prob(forward.rate, backward.rate, my.tree$edge.length[i])
 # }
  prob.array=exp_array_of_edges(my.tree$edge.length, M)

   ##return(partLike(my.tree$edge, my.tree$Nnode, my.data, 2, prob.array))
  return(Ln_Markov1(my.tree$edge, my.tree$Nnode, my.data, 2, prob.array))
  #return(partLike(my.tree$edge, my.tree$Nnode, my.data, 2, my.tree$edge.length, M))
}



############

ptm <- proc.time()
for (i in 1:10000){my_cube(2, my.tree$edge.length, M)}
proc.time() - ptm

ptm <- proc.time()
for (i in 1:10000){
  prob.array = array(0, dim=c(2,2,length(my.tree$edge.length)))
  for (i in 1:length(my.tree$edge.length)){
    prob.array[,,i] = two.state.trans.prob(forward.rate, backward.rate, my.tree$edge.length[i])
  }
}
proc.time() - ptm

ptm <- proc.time()
for (i in 1:10000){batch_exp(my.tree$edge.length, M))}
proc.time() - ptm



ptm <- proc.time()
for (i in 1:100000){two.state.part.like(test.tree, test.tree$tip.state, 0.02, 0.05)}
proc.time() - ptm

mk2.lik = make.mk2(test.tree, test.tree$tip.state)
tes=exp(mk2.lik(c(0.02,0.05),root=ROOT.BOTH))
sum(tes*.5)

ptm <- proc.time()
for (i in 1:100000){make.mk2(test.tree, test.tree$tip.state)}
proc.time() - ptm

ptm <- proc.time()
for (i in 1:100000){two.state.part.like(test.tree, test.tree$tip.state, M)}
proc.time() - ptm

two.state.part.like(test.tree, test.tree$tip.state, M)




ptm <- proc.time()
for (i in 1:100000){Ln_Markov1(my.tree$edge, my.tree$Nnode, my.data, my.tree$edge.length, M)}
proc.time() - ptm

ln1=Ln_Markov1(my.tree$edge, my.tree$Nnode, my.data, my.tree$edge.length, M)
ln1*2

Ln_Markov1(my.tree$edge, my.tree$Nnode, my.data, my.tree$edge.length, M, c(.5,.5))


[1,] 6.907984e-08 2.724959e-08
