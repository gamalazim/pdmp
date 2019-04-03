#' @import methods
#' @import utils
#' @useDynLib pdmp

###################################################
## Start pedigrees as 'data.frames' with 'factors'
###################################################

## something that can be used the same way as is.na
is.unk <- function(x) {x == "0"}

#' Recode pedigree data using integer numbers
#'
#' \code{pedRecode} expects a pedigree and the column locations of the individual and the two parents.
#' The function then transforms the character IDs in the 3 columns into integer values. The function
#' recods IDs consistantly across the 3 columns so that a parent that also appears as an individual
#' will have a single integer ID throughout the pedigree.
#'
#' @param p data frame or matrix with at least 3 columns for the individual and two parents.
#' @param cfm vector of column names or positions of individual, paternal parent and maternal parent.
#'
#' @details The column of individuals determined by \code{p[,cfm[1]]} is coded from 1 to N. Any parent
#' who appears as an individual will take his/her integer value of the individuals' column. This
#' process of recoding is normally performed after structuring and sorting the pedigree. Usually
#' pedigree data would be free of any inconsistencies before performing this recoding. See 'Examples'.
#'
#' @return data frame for the recoded pedigree where columns \code{p[,cfm]} are replaced by integer
#' values consistent across the 3 columns.
#'
#' @seealso
#' \code{\link{pedCensor}}, \code{\link{pedSort}}
#'
#' @examples
#' ## Load pedigree data
#' data(pedigree)
#'
#' ## Restructure pedigree data in 'p01' and remove inconsistencies
#' p.clean <- pedCensor(p01)
#'
#' ## Sort in chronological order
#' p.sort <- pedSort(p.clean)
#'
#' ## Recode sorted pedigree
#' p.recode <- pedRecode(p.sort)
#'
#' @export
pedRecode <- function(p, cfm=1:3) {

	if(is.data.frame(p) || is.matrix(p)) {
		if(length(cfm) != 3 || max(cfm) > dim(p)[2])
			stop("incorrect dimensions for animal, sire and dam")
		if( pedExamine(p, signal=T) )
			stop("incorrect pedigree structure,\n... use function 'pedExamine' to reveal inconsistencies ...")

		if(is.matrix(p)) p <- as.data.frame(p)

		if(!is.factor(p[,cfm[2]])) p[,cfm[2]] <- factor(p[,cfm[2]])
		if(!is.factor(p[,cfm[3]])) p[,cfm[3]] <- factor(p[,cfm[3]])

		levels(p[,cfm[2]]) <- match(levels(p[,cfm[2]]), p[ ,cfm[1]], nomatch=0)
		levels(p[,cfm[3]]) <- match(levels(p[,cfm[3]]), p[ ,cfm[1]], nomatch=0)

		p[,cfm[1]] <- 1:nrow(p)
		for(i in c(cfm[2], cfm[3])) p[,i] <- as.integer(as.character(p[,i]))
		p
	}

	else {
		stop("argument must be a data.frame or matrix")
	}
}

#' Sort pedigree data in chronological order so parents appear before their progeny
#'
#' \code{pedSort} receives a consistent pedigree and the column locations of the individual and
#' the two parents. The function then performs sorting in chronological oredr without needing
#' a sorting factor such as birth dates. \code{pedSort} is super fast!
#'
#' @param p.org data frame or matrix with at least 3 columns for the individual and two parents.
#' @param cfm vector of column names or positions of individual, paternal parent and maternal parent.
#' @param recode logical for whether or not to output the sorted pedigree in integer-coded format.
#'
#' @details Given an input of consistent pedigree, chronological odering is performed on the pedigree.
#' The input pedigree is checked first by \code{pedExamine} for all inconsistencie and an error is
#' generated if inconsistencies are encountered. It is recommended to always run \code{pedCensor} before
#' \code{pedSort}. See \sQuote{Examples} below.
#'
#' @return data frame for the sorted pedigree
#'
#' @seealso
#' \code{\link{pedCensor}}, \code{\link{pedExamine}}
#'
#' @examples
#' ## Load pedigree data
#' data(pedigree)
#'
#' ## Restructure pedigree data in 'p01' and remove inconsistencies
#' p.clean <- pedCensor(p01)
#'
#' ## Sort in chronological order
#' p.sort <- pedSort(p.clean)
#'
#' ## Sort pedigree and generate integer-coded output
#' p.sort.recoded <- pedSort(p.org = p.clean, recode = TRUE)
#'
#' @export
pedSort <- function (p.org, cfm=1:3, recode=FALSE)
# input pedigree should be free of inconsistencies, pass pedExamine, etc
# recode=FALSE will not ouput a recoded pedigree
{
        p <- pedRecode(p.org, cfm)
        np <- setdiff(p[,cfm[1]], c(p[,cfm[2]], p[,cfm[3]]))
        dex <- 0:nrow(p)

        if( any(p[,cfm[2]]>p[,cfm[1]]) || any(p[,cfm[3]]>p[,cfm[1]]) ) {
                PEDSORTOUT <- .C("pedSortMain",
                as.integer(c(0,p[,cfm[1]])),
                as.integer(c(0,p[,cfm[2]])),
                as.integer(c(0,p[,cfm[3]])),
                as.integer(c(0,np)),
                rerank = as.integer(dex),
                as.integer(nrow(p)),
                as.integer(length(np)),
                Base = integer(1),
                PACKAGE = "pdmp")

                if(recode)
                        pedRecode(p.org[(PEDSORTOUT$rerank[-1]),])
                else
                        p.org[(PEDSORTOUT$rerank[-1]),]
        }

        else stop("Already Sorted Pedigree")

}


#' Examine pedigree data for specific inconsistencies
#' \code{pedExamine} receives a messy pedigree and the column locations of the individual and
#' the two parents. The function then performs specific integrity checks and outputs results as
#' a list of problematic individuals.
#'
#' @param p data frame or matrix with at least 3 columns for the individual and two parents.
#' @param cfm vector of column names or positions of individual, paternal parent and maternal parent.
#' @param recode logical for whether or not to output the sorted pedigree in integer-coded format.
#'
#' @details Given an input of consistent pedigree, chronological odering is performed on the pedigree.
#' The input pedigree is checked first by \code{pedExamine} for all inconsistencie and an error is
#' generated if inconsistencies are encountered. It is recommended to always run \code{pedCensor} before
#' \code{pedSort}. See \sQuote{Examples} below.
#'
#' @return data frame for the sorted pedigree
#'
#' @seealso
#' \code{\link{pedCensor}}, \code{\link{pedExamine}}
#'
#' @examples
#' ## Load pedigree data
#' data(pedigree)
#'
#' ## Perform all possible checks and print a 'LIST' of messed up records
#' LIST.messed <- pedExamine(p01)

#' ## Limit checking to 'misid'-entified records and print whole records. LIST$misid is output!! Still a LIST.
#' head(p01[ pedExamine(p01,examine="misid")$misid, ])
#'
#'           individual              Sire              Dam
#' #42   BR2CTY111545932    BR2CTY12338584  BR2CTY110892064
#' #46  BR2CTY1817309548    BR2CTY12347416  BR2CTY110788816
#' #67  BR2CTY1817821648    BR2CTY12388072  BR2CTY111532864
#' # ...
#'
#' #As shown, pedExamine identified BR2CTY111545932, BR2CTY1817309548, and BR2CTY1817821648 as misidentified IDs.
#' #We may investigate further the first individual to find out that the individual appeared as herself and her dam
#' #on record 101334. See below,
#'
#' p01[p01$Dam == 'BR2CTY111545932',]
#'
#' #              individual           Sire             Dam
#' #78     BR2CTY436082749120 BR2CTY12317968 BR2CTY111545932
#' #80     BR2CTY436090674292 BR2CTY12367660 BR2CTY111545932
#' #101334    BR2CTY111545932 BR2CTY12338584 BR2CTY111545932
#'
#' #Finally, note that misid will almost always have to be removed. To remove use,
#' x <- pedExamine(p01, examine="misid") # Limit examine to wrong ID's
#' p01.NOmisid <- p01[!is.element(p01$individual, x$misid),] # Remove wrong ID's

## pedExamine(ped, examine = 'misid')$misid
## character(0) ## indicates that misid's were correctly cleaned up!

#' @export
"pedExamine" <-
  function(p, cfm=1:3,
           examine = c("duplicate", "parent.record", "overlap", "misid"), signal=F) {

    LIST <- list()

    ## a list of duplicated individuals
    if(any(pmatch(examine, "duplicate", nomatch=0))) {
      LIST$duplicate <- p[duplicated(p[,cfm[1]]), cfm[1]][ ,drop=T]
      if(signal) return(0)
    }

    ## a list of parents that never appeared as individuals, ie in cfm[,1]
    if(any(pmatch(examine, "parent.record", nomatch=0))) {
      LIST$missing.parent.record <- setdiff(c(as.character(p[,cfm[2]]), as.character(p[,cfm[3]])), p[,cfm[1]])
      LIST$missing.parent.record <- LIST$missing.parent.record[!is.unk(LIST$missing.parent.record)]
      if(signal) return(0)
    }

    ## 'overlap' Base individuals that are not in the sire or dam columns
    ## Sometimes those are needed to be included any way
    if(any(pmatch(examine, "overlap", nomatch=0))) {
      base <- p[p[,cfm[1]] != 0 & p[,cfm[2]]==0 & p[,cfm[3]]==0, cfm[1]]
      LIST$overlap <- setdiff(base, c(as.character(p[,cfm[2]]), as.character(p[,cfm[3]])))
      if(signal) return(0)
    }

    ## identify row names of 0-id's for individuals; id == sire; id == dam; and known sire is also a known dam
    ## Sometimes needed as in the MGS pedigree
    if(any(pmatch(examine, "misid", nomatch=0))) {
      LIST$misid <- union( row.names(p[ p[,cfm[1]] == 0 |
        as.character(p[,cfm[1]])==as.character(p[,cfm[2]]) |
        as.character(p[,cfm[1]])==as.character(p[,cfm[3]]),]),
        row.names(p[(p[,cfm[3]] != 0 &
        is.element(as.character(p[,cfm[3]]), as.character(p[,cfm[2]]))),]))

      ## Controversial individuals, eg, <A B C> and <A B D>
      p.tmp <- unique(p[,cfm[1:2]])
      p.tmp <- p.tmp[!is.unk(p.tmp[,cfm[1]]) & !is.unk(p.tmp[,cfm[2]]),]
      dupFather <- p.tmp[duplicated(p.tmp[,cfm[1]]), cfm[1]][,drop=T]

      p.tmp <- unique(p[,cfm[c(1,3)]])
      p.tmp <- p.tmp[!is.unk(p.tmp[,cfm[1]]) & !is.unk(p.tmp[,cfm[2]]),]
      dupMother <- p.tmp[duplicated(p.tmp[,cfm[1]]), cfm[1]][,drop=T]

      dupFM <- union(as.character(dupFather), as.character(dupMother))
      LIST$misid <- union( row.names(p[is.element(p[,cfm[1]], dupFM),]), LIST$misid)

      if(signal) return(0)
    }

    ## Diagnose unsorted pedigrees
    message("\n Done ...\n Examine doesn't check for unsorted pedigree, this is not an error\n")

    LIST
  }



## Pedigree Stacking Function.
## Syntax:
##       pedStack(pedFile, c(ind, father, mother), list of vectors of 3 elements each)
##       pedStack(ped, c(1,2,3), list( c("s", "ss", "ds"),
##				       c("d", "sd", "dd"),
##				       c(4, 8, 9),
##				       c(5, 10, 11),
## 				       etc,) )

## This is not just stacking of anything. For example, the function wouldn't repeat
## animal records of sires and dams! So it is a smart/pedigree stacking not just
## stacking vectors over vectors!! Also notice that the function can utilize names
## and/or vector numbers.

#' Stack relatives in a consistent 3-column pedigree format

#' @export
pedStack <- function(p, v, ...) {
	for(i in 1:ncol(p))
		if(!is.factor(p[,i])) p[,i] <- as.factor(p[,i])

	dots <- (...)
	N <- length(dots)
	D <- p[,v]
	for(i in 1:N) {
		d <- p[,dots[[i]]]
		names(d) <- names(D)
		D <- rbind(D, d)
	}
	unique(D)
}

## Example:
##   For the most scarse pedigree data ever, I asked a producer to supply 'complete' pedigree data
##   for his cows! It is pretty typical that the herdman will give a list of cows born in his herd
##   along with their dams and service sires used to inseminate the dams. Then automatically
##   the grand sire, great grand sire, and great great grand sire are all known for a given service
##   bull.
##
##                           Dc  [Great Great Grand Sire]
##                            \  /
##                             \/
##                 Db  [Great Grand Sire]
##                  \  /
##                   \/
##           Da  [Grand Sire]
##            \  /
##             \/
##   Dam  [Service Bull]
##     \  /
##      \/
##   [Cow]


## read.table("~/GGA/stotz/block.data.1",sep=" ")->garbage
## cbind(garbage, Da=paste("Da", 1:nrow(garbage), sep="")) -> garbage  ## Dummy Variable Da
## cbind(garbage, Db=paste("Db", 1:nrow(garbage), sep="")) -> garbage  ## Dummy Variable Db
## cbind(garbage, Dc=paste("Dc", 1:nrow(garbage), sep="")) -> garbage  ## Dummy Variable Dc
## garbage$unk <- 0  ## column of unknowns ...

## x <- pedStack(garbage, c(1,3,4), list( c("V3", "V6", "Da"),
## c("V4", "unk", "unk"), c("V7", "Db", "V8"), c("Da","unk","unk"),
## c("V8", "V9", "Dc"), c("Db", "unk", "unk"), c("V9", "unk", "unk"),
## c("Dc", "unk", "unk")))


## pedEdit is specific and more technical, must understand its inner workings which is that l
## must be provided or that l <- pedExamine(p) has been run first. In this regard pedCensor is
## more general because it runs pedExamine first anyway! In addition pedEdit() edits selected or all
## items based on which of the edit arguments is set to FALSE, by default all edits are set to TRUE

## IMPORTANT: This function does particular edits and skips others (based on options supplied), therefore,
## inconsistencies can be created. Even if all edits are requested, pedEdit edits the supplied pedigree one
## time. If you are not sure about those inconsistencies in your particular pedigree data, only resort
## to pedCensor that iteratively calls pedEdit until all inconsistencies have been removed!!

## In certain cases you may need to call pedEdit iteratively based on special supplied list.
## Another problem with using pedEdit is that you 'can' confuse it if the 'l' provided is not the same as
## pedExamine(p)!! But that's of course a mistake of usage!!

#' Edit particular inconsistencies
#'
#' @export
"pedEdit" <-
  function(p, cfm=1:3, l,
           duplicate.rm=TRUE,
           parent.record.add=TRUE,
           overlap.rm=TRUE,
           misid.rm=TRUE) {

    if(length(l$duplicate) == 0) duplicate.rm <- FALSE
    if(length(l$missing.parent.record) == 0) parent.record.add <- FALSE
    if(length(l$overlap) == 0) overlap.rm <- FALSE
    if(length(l$misid) == 0) misid.rm <- FALSE

    if(duplicate.rm == TRUE) { ## order by k which is the amount of info (0, 1, or 2)
      k <- (p[,cfm[2]] != 0) + (p[,cfm[3]] != 0)
      p <- p[order(k,decreasing=TRUE),]
      y <- row.names(p[duplicated(p[,cfm[1]]),])
      p[y, cfm[1]] <- NA ## assign NA to individuals with less info, rm @ the end
    }

    if(misid.rm == TRUE) p[l$misid, cfm[1]] <- NA

    if(overlap.rm == TRUE) {
      if(any(is.na(p[,cfm[1]]))) {
        p <- p[!is.na(p[,cfm[1]]), ]
        ## re-run pedExamine on clean 'p' with examine='overlap'
        p[is.element(as.character(p[,cfm[1]]), as.character(pedExamine(p, cfm, examine='overlap')$overlap)), cfm[1]] <- NA
      }
      else
        p[is.element(as.character(p[,cfm[1]]), as.character(l$overlap)), cfm[1]] <- NA
    }

    if(parent.record.add == TRUE) {
      if(any(is.na(p[,cfm[1]]))) {
        p <- p[!is.na(p[,cfm[1]]), ]
        ## re-run pedExamine on clean 'p' with examine='parent.record'
        y <- data.frame(cbind(pedExamine(p, cfm, examine = 'parent.record')$missing.parent.record, 0, 0))
        if(dim(y)[2] >= 3) {
          names(y) <- names(p)
          p <- rbind(y, p)
        }
      }
      else {
        y <- data.frame(cbind(l$missing.parent.record, 0, 0))
        names(y) <- names(p)
        p <- rbind(y, p)
      }
    }
    p[!is.na(p[,cfm[1]]), ]
  }

## pedCensor(p, cfm=1:3, only.check=TRUE) is similar to pedExamine(p), HOWEVER, for specific examinations, you
## still need pedExamine, e.g., pedExamine(p, examine == 'overlap') !!
## Dec 2015: Added overlap.rm = TRUE/FALSE, sometimes those are needed
##

#' Check for all consistencies and fix them all
#'
#' @export
"pedCensor" <-
  function(p, cfm=1:3, only.check = FALSE, overlap.rm=TRUE) {
    LIST <- pedExamine(p,cfm)
    if(only.check == TRUE) LIST
    else {
      if(overlap.rm) {
        while(any(sapply(LIST, length) > 0)) {
          p <- pedEdit(p, cfm, LIST)
          LIST <- pedExamine(p,cfm)
        }
      }
      else{
        while(any(sapply(LIST[-3], length) > 0)) {
          p <- pedEdit(p, cfm, LIST, overlap.rm=F)
          LIST <- pedExamine(p,cfm)
        }
      }

      p
    }
  }

## pedBind binds two pedigree frames based on the first 3 col's (default) of x and y,
## by default names are taken from the first 3 col names of x unless otherwise supplied
## as specific int col numbers in arg *cfm or specific character names in arg Names.
## Dec 2015: Added overlap.rm = TRUE/FALSE, sometimes those are needed
##

#' Bind two pedigree frames to form a new frame without repetition or inconsistencies
#'
#' @export
"pedBind" <-
  function(x, y, xcfm=1:3, ycfm=1:3, Names=names(x)[xcfm], overlap.rm=TRUE) {
    names(y)[ycfm] <- names(x)[xcfm] <- Names
    if(overlap.rm) pedCensor(rbind(x[,xcfm], y[,ycfm]))
    else pedCensor(rbind(x[,xcfm], y[,ycfm]), overlap.rm=F)
  }



## To call for the 1st time don't provide any lpanc - just 'paternal.ancestors(p, i)'
#' Extract paternal ancestors
#'
#' @export
"paternal.ancestors" <-
  function(p, i, lpanc=c(0, i) ) {
    lpanc <- c(lpanc, x<-p[(p[,1]==i), 2])
    if( !any(duplicated(lpanc)) ) paternal.ancestors(p, x, lpanc) else lpanc
  }

## To call for the 1st time don't provide any lmanc - just 'maternal.ancestors(p, i)'
#' Extract maternal ancestors
#'
#' @export
"maternal.ancestors" <-
  function(p, i, lmanc=c(0, i) ) {
  lmanc <- c(lmanc, x<-p[(p[,1]==i), 3])
  if( !any(duplicated(lmanc)) )  maternal.ancestors(p, x, lmanc) else lmanc
  }

#' Identify loops in a pedigree frame
#' @export
"pedLoop" <-
  function(p, i) {
    x <- paternal.ancestors(p, i)
    y <- maternal.ancestors(p, i)
    LIST <- list()

    ifelse(x[length(x)]==0,
      LIST$lpanc <- data.frame(x[-1], c(x[-c(1,2)],0)),
      { LIST$lpanc <- data.frame(x[-1], c(x[-c(1,2)],NA)); warning("trouble"); } )

    ifelse(y[length(y)]==0,
      LIST$lmanc <- data.frame(y[-1], c(y[-c(1,2)],0)),
      { LIST$lmanc <- data.frame(y[-1], c(y[-c(1,2)],NA)); warning("trouble"); } )


    names(LIST$lpanc) <- names(p)[1:2]
    names(LIST$lmanc) <- names(p)[c(1,3)]
    LIST
}

## Return the two parents
## if i is a vector of IDs, the value returned is simply p[i,], i.e, pedigree
## section with i, s, and d ('i' being the input IDs).

#' Search for parents and return them for a vector of individuals
#'
#' @export
"parents" <-
  function(p, i) {
    if(length(i) == 1) {
      if(!is.unk(i)) {
        return ( sapply(p[which(i == p[,1]), 2:3], function(x) levels(factor(x)[,drop=T])) )
      }
      else return ( c('0', '0') )
    }
    if(length(i) > 1) {
      p[match(i, p[,1]),]
    }
}

#' Return the maternal grandsire of an individual
#' @export
"mgs" <-
  function(p, i) {
    parents(ped, parents(ped, i)[2])[1]
}

##
##
##   All non parents:
##   np <- setdiff(p[,cfm[1]], c(p[,cfm[2]], p[,cfm[3]]))
##   used aftre clearing all obvious inconsistencies ... after pedCensor at least!
##   pedigree frame doesn't have to be recoded by 'pedRecode()'

#' Return all ancestors of a subset of individuals
#'
#' @export
"ancestors" <-
  function(p, i, ng) { ## i is a character string e.g., 'HOUSA000015563369', '40', or 'Ahmad'
    # planc is the previous list of ancestors:
    RESULT <- planc <- i ## first on the list is the individual sent to 'ancestors()'
    for(g in 1:ng) { ## ng number of past geneations to search through
       ## Direct ancestors in a generation (e.g., in past gen. 1 = 2 parents then 4 grandparents, etc.):
       iig <- vector(mode = "character", length = (2^g))
       for(j in 1:(2^(g-1))) {
         iig[(2*j -1):(2*j)] <- parents(p, planc[j])
       }
       if(all(iig == "0")) break
       planc <- iig
       RESULT <- c(RESULT, iig)
    }
    RESULT
}

## Not efficient
## i=actual ID not a rank
#' Clear a loop in the pedigree
#' @export
"clearLoop" <-
  function(p, i, ng) {
     x <- ancestors(p, i, ng)
     y <- x[x != 0]
     TOUNK <- unique(y[duplicated(y)])

     planc <- i
     for(g in 1:ng) { ## ng number of past geneations to search through
       ## Direct ancestors in a generation (e.g., in past gen. 1 = 2 parents then 4 grandparents, etc.):
       iig <- vector(mode = "character", length = (2^g))
       for(j in 1:(2^(g-1))) {
         iig[(2*j -1):(2*j)] <- parents(p, planc[j])
         if( is.element(iig[(2*j -1)], TOUNK) ) p[(p[,1] == planc[j]), 2] <- "0"
         if( is.element(iig[(2*j)], TOUNK) ) p[(p[,1] == planc[j]), 3] <- "0"
       }
       if(all(iig == "0")) break
       planc <- iig
    }
    p
}


## Not efficiesnt
## loopi is a vector of individuals not their ranks in the pedigree ..

#' Clear loops in pedigree for a vector individuals
#'
#' @export
"clearLoops" <-
  function(p, loopi, ng) {
    for(i in as.character(loopi))
      p <- clearLoop(p, i, ng)
    p
}

## pedCensor and pedSort will be needed after the extraction...

#' Extract all relevant pedigree info for a subset of individuals
#'
#' @export
"pedExtract" <-
  function(p.org, sample) {
    p <- pedRecode(p.org)
    sample <- factor(sample)[,drop=T]
    levels(sample) <- match(levels(sample), p.org[ ,1], nomatch=0)
    xtract <- as.numeric(levels(sample))

    repeat {
       newxt <- setdiff(c(p[xtract, 2], p[xtract, 3]), xtract)
       newxt <- newxt[newxt != 0]
       if(length(newxt) > 0) xtract <- c(newxt, xtract)
       else break
    }
    p.org[xtract,]
}

##
## March 31, 2011
##
## Don't use pedCensor after this function many things will look like
## inconsistency to pedExamine such as individuals apearing as both sires and dams
## Can only use pedExamine(pedigree, examine=c("duplicate", "parent.record", "overlap"))
## which should not matter if i/s/d pedigree was censored before pedMgs.
##
## To extract bulls and get rid of females in the individual's column,
## use pedExtract(pedigree, sample = a list of bulls from data)
##
## Reasonably fast for 1.2 million pedigree, it took 3.7 seconds
#  system.time(p <- pedMgs(ped))
#    user  system elapsed
#   3.652   0.016   3.667

#' Create a male-only pedigree by replace all dams with their sires to get
#' a sire-maternal-grandsire pedigree
#'
#' @export
"pedMgs" <-
  function(ped) {
    p = ped
    levels(p[,3]) = match(levels(p[,3]), p[,1], nomatch=0)
    p[,3] = as.numeric(as.character(p[,3]))
    p[p[,3] != 0, 3] = as.character(p[p[p[,3] != 0 ,3] ,2])
    p[,3] = factor(p[,3])[,drop=T]
    p
}

#' Create the Numerator Relationship Matrix with \code{sparse} option
#'
#' @export
"kinship" <-
#
# Assumes consistent, sorted, and 1-n recoded pedigree
#
  function(p, cfm=1:3, sparse = FALSE) {

    ORD <- nrow(p)

    KINSHIPOUT <- .C("KinshipMain",
       as.integer(c(0,p[,cfm[2]])),
       as.integer(c(0,p[,cfm[3]])),
       as.integer(ORD),
       A = c(0, as.double(matrix(0, ORD, ORD))),
       PACKAGE = "pdmp" )

    matrix(KINSHIPOUT$A[-1], ORD, ORD, byrow=T)

}

#' Get the L factor of A, where A = L.D.t(L)
#'
#' Get the \code{L} factor of \code{A}, with \code{sparse} option, where \code{A = L.D.t(L)}
#' Assumes consistent, sorted, and 1-N recoded pedigree
#'
#' @export
"getL" <-
  function(p, cfm=1:3, sparse = FALSE) {

    ORD <- nrow(p)

    getLOUT <- .C("getLMain",
       as.integer(c(0,p[,cfm[2]])),
       as.integer(c(0,p[,cfm[3]])),
       as.integer(ORD),
       L = c(0, as.double(matrix(0, ORD, ORD))),
       PACKAGE = "pdmp" )

    matrix(getLOUT$L[-1], ORD, ORD, byrow=T)
}


##
#' Get Nnmber of paternal and maternal generations in a pedigree
#' @export
"numGen" <-
  function(p, anc = 1, inds) {
# p sorted and 1-n recoded pedigree
# anc = 1 (paternal) or 2 (maternal)
# inds: vector of individuals to calculate number of generations for

    anc = anc+1
    n = length(inds)
    RESULT <- numeric(n)
    for(k in 1:n) {
      tng = 0; j = inds[k];
      while(p[j, anc] != 0) {tng = tng + (p[j, anc] != 0); j = p[j, anc]}
      RESULT[k] = tng
      }
      RESULT
}

## To use, repeatedly get a random sample out of the nonparents and
## average their number of ancestors, see the example below:

# Time elapsed: about 60 seconds for 20 samples of 20 individuals per sample and a pedigree
# of over 600,000 individuals

# x <- numeric(20); for(i in 1:20) x[i] <- mean(numGen(p, anc=1, inds=sample(300000:600000, size=20)))
# x
# [1] 5.85 6.80 5.20 5.85 6.00 6.45 6.05 5.95 6.30 5.55 6.15 5.10 6.15 7.05 5.95
# [16] 5.85 6.40 6.10 4.95 6.35
# mean(x)
# [1] 6.0025
#


#' Return the inverse of the Numerator Relationship Matrix
#'
#' Given consistent pedigree, return the inverse of the Numerator Relationship Matrix either in
#' dense format (default) or in sparse format (\code{sparse = TRUE}).
#'
#' @export
"kinInv" <-
#
# Assumes consistent, sorted, and 1-n recoded pedigree
# Not complete yet --
  function(p, cfm=1:3, sparse = FALSE) {

    ORD <- nrow(p)
    iCoef <- diag( kinship(p = p, cfm = 1:3, sparse = FALSE) ) - 1

    KININVOUT <- .C("KininvMain",
       as.double(iCoef),
       as.integer(c(0,p[,cfm[2]])),
       as.integer(c(0,p[,cfm[3]])),
       as.integer(ORD),
       Ainv = c(0, as.double(matrix(0, ORD, ORD))),
       PACKAGE = "pdmp" )

    matrix(KININVOUT$Ainv[-1], ORD, ORD, byrow=T)

}


