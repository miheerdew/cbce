jaccard_sim <- function(A, B){
  length(intersect(A, B))/length(union(A,B))
}

# P(A|B)
OL_frac <- function(A, B) {
  length(unique(intersect(A,B)))/length(unique(B))
}

#' Represents P(c1|c2)
#'
#' @export
P <- function(c1, c2, wt.x, wt.y) {
  wt.x*OL_frac(c1$x, c2$x) + wt.y*OL_frac(c1$y, c2$y)
}

#Lift binary op to bipartite sets
setop <- function(c1, c2, op) {
 list(x=op(c1$x, c2$x), y=op(c1$y, c2$y) )
}

com_intersect <- function(c1, c2) setop(c1, c2, op=intersect)

com_union <- function(c1, c2) setop(c1, c2, op=union)

coms_union <- function(com) {
  Reduce(com_union, com, list(x=c(), y=c()))
}
com_unique <- function(c) list(x=unique(c$x), y=unique(c$x))

split_xy <- function(A, dx) {
  list(x=A[A <= dx], y=A[A > dx] - dx)
}

comms <- function(res, dx) {
  res$extract_res %>>%
    list.filter(exists("StableComm") && length(StableComm) > 0) %>>%
    list.map(split_xy(StableComm, dx))
}

comms_sim <- function(sim) {
  list.map(sim$bms, split_xy(., sim$dx))
}

OL_matrix <- function(com1, com2=com1, wt.x=0.5, wt.y=1 - wt.x) {
 com1 %>>%
    list.map(c1 ~
               list.mapv(com2, c2 ~ P(c2,c1, wt.x=wt.x, wt.y))
           ) %>>%
    list.rbind
}

#Returns the fraction of pairs of c found in com
pair_coverage <- function(c, com, com_weights=rep(1, length(com))) {
  c <- com_unique(c)
  Counts <- matrix(0, length(c$x), length(c$y))
  for (i in 1:length(com)) {
    c2 <- com_intersect(c, com[[i]])
    for (x in c2$x) {
      for (y in c2$y) {
       xindx <- which(c$x == x)
       yindx <- which(c$y == y)
       Counts[xindx, yindx] <- max(Counts[xindx, yindx], com_weights[i])
      }
    }
  }
  mean(Counts)
}

node_coverage <- function(c, com,
                          wt.com=rep(1, length(com)),
                          wt.x=0.5,
                          wt.y=1-wt.x) {
  ord <- order(wt.com, decreasing = TRUE)

  lx <- ly <- 0
  cx <- c$x
  cy <- c$y

  for(i in ord) {
    ix <- intersect(cx, com[[i]]$x)
    iy <- intersect(cy, com[[i]]$y)

    lx <- lx + wt.com[i]*length(ix)
    ly <- ly + wt.com[i]*length(iy)

    cx <- setdiff(cx, ix)
    cy <- setdiff(cy, iy)
  }
  wt.x*(lx/length(c$x)) + wt.y*(ly/length(c$y))
}

is_element <- function(c, com, wt.x=0.5, wt.y=1-wt.x) {
  max(list.mapv(com, similarity(c, ., wt.x = wt.x, wt.y = wt.y)))
}

similarity <- function(c1,c2, wt.x=0.5, wt.y=1-wt.x) {
   wt.x*jaccard_sim(c1$x, c2$x) + wt.y*jaccard_sim(c1$y, c2$y)
}
#For each comm in com1 find fraction of pairs found in com2

penalty_func1 <- function(nums){
  nums^2
}

#' Compare two decks of (bi-)communities
#'
#' Given two decks of bi-communities, compare these two decks (for similariy and differences)
#' using a collection of metrics.
#'
#' @param comms1,comms2 list. Deck of communites. For instance the object returned by \code{\link{comms}}
#' @param wt.x,wt.y numeric. The weightage to be given to X, Y sides while comparing approximate distances.
#'                   For instances if wt.x = 1, then two communites with same X set will be considered equal.
#' @param weights function. This is a weightage function which is used to assess if two communities are equal/contained
#'                to each other based on their (perhaps one way) jaccard distances. The default value is
#'                \code{function(nums) nums > .9 } which is a binary function, which asseses yes/no based on the theshold .9,
#'                but it is possible to have a smoother transition.
#' @param weak_coverage Boolean. For "coverage" metric, use the strong or weak version.
#'
#' @return A list containing following
#' \describe{
#' \item{report1, report2}{report1, report2 are lists with entries corresponding to each
#' community in comms1 and comms2 repectively. Each entry is a report which compares that community to the other pile.
#' This report has
#'  \describe{
#'  \item{coverage}{Coverage is a measure of "how much of this community is found in the other deck combined".
#'  Weak coverage, only looks at the X and Y sets, while strong coverage actually looks at number of pairs.}
#'  \item{outer_coverage}{(Actually, this metric is not very useful)
#'    This measures how well the community is covered by communities in the other deck that "contain it".
#'    Quotations were used in the previous sentence because the notion is approximately used. (see the weights parameter)
#'    This number will typically be close to 1 if the community is contained inside some community from the other deck.}
#'  \item{inner_coverage}{This measures how well the community is covered by communities in the other deck that are "contained in it".
#'      This number will typically be close to 0.}
#'  \item{closest_match}{This measures the closest match to the communities from the other deck }
#'  }
#' }
#' \item{OL12, OL21}{These are overlap matrices. \code{OL12[i,j]} gives the fraction of
#'  community i that is present in community j}
#' }
#' @export
compare <- function(comms1, comms2,
                    wt.x=0.5, wt.y=0.5,
                    weights=function(nums) nums > 0.9,
                    weak_coverage = TRUE) {

  OL <- function(cs1, cs2) {
    OL_matrix(cs1, cs2, wt.x=wt.x, wt.y=wt.y)
  }
  coverage <- function(c1, com2, wt.com=rep(1, length(com2)), weak=weak_coverage) {
    if(weak) {
      node_coverage(c1, com2, wt.com=wt.com, wt.x=wt.x, wt.y=wt.y)
    } else {
      pair_coverage(c1, com2, wt.com)
    }
  }

  OL12 <- OL(comms1, comms2)
  OL21 <- OL(comms2, comms1)
  # C12 <- cover(comms1, comms2)
  # C21 <- cover(comms2, comms1)

  report <- function(c1, ol12, ol21, comm2) {
    l1 <-
      list(
         #How well is it covered by comms?
         coverage=coverage(c1, comm2),

         #How well is it covered by comms it is contained in?
         outer_coverage=coverage(c1, comm2, wt.com = weights(ol12), weak = TRUE),

         #How well is it coverged by comms that are contained in it?
         inner_coverage=coverage(c1, comm2, wt.com = weights(ol21), weak = TRUE),

         #How closes is c1 to a community in comm1
         closest_match=is_element(c1, comm2, wt.x = wt.x, wt.y = wt.y)
        )
  }

  list(report1=list.stack(list.map(comms1, report(., OL12[.i,], OL21[,.i], comms2))),
       report2=list.stack(list.map(comms2, report(., OL21[.i,], OL12[,.i], comms1))),
       OL12=OL12,
       OL21=OL21)
}

#' Calculate the effective bimodules based on its OL_matrix
#'
#' @param OLmat The overlap matrix of the bimodule as returned by OL_matrix(comm)
#' @param weight_fun Function to determine when two communities are considered the same.
#' @param sym Boolean. If TRUE the equality is both way inclusion else just one way inclusion.
#' @return A count representing the effective number of communities
#' @export
effective_number <- function(OLmat, weight_fun=function(x) x > 0.95, sym=TRUE) {
  if(sym) {
    OLmat <- (OLmat + t(OLmat))/2
  }

  diag(OLmat) <- 0

  1:nrow(OLmat) %>>%
   list.mapv(1 - weight_fun(max(OLmat[., 1:.]))) %>>%
    sum
}

analyse <- function(com){
 OLx <- OL_matrix(com, com, wt.x = 1)
 OLy <- OL_matrix(com, com, wt.x = 0)
 OL <- OL_matrix(com, com, wt.x = 0.5)
 #Add effective number.
}