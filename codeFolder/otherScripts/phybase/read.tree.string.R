`read.tree.string` <-
function (file = "", format = "nexus") 
{
    
    X <- scan(file = file, what = "character", sep = "\n", quiet = TRUE)
    X <- paste(unlist(strsplit(paste(X,collapse=""),split=";")),";",sep="")
    
    y1<-grep("[[]", X, ignore.case = TRUE)
    y2<-grep("[]]", X, ignore.case = TRUE)

    if(length(y1) > 0)
    {
		for(i in 1:length(y1))
		{
			if(y1[i] == y2[i])
			{
				X[y1[i]]<-gsub("[[].*[]]","",X[y1[i]])
			}
			else
			{
				X[y1[i]]<-gsub("[[].*","",X[y1[i]])
				X[y2[i]]<-gsub(".*[]]","",X[y2[i]])
				X[(y1[i]+1):(y2[i]-1)]<-""
			}
		}
    }
    X<-X[X!=""]

    if (length(grep("phylip",format,ignore.case=TRUE)))
    {
	tree<-X
	translation <- FALSE
    }
    else{
    i1 <- grep("BEGIN TREES;", X, ignore.case = TRUE)
    i2 <- grep("TRANSLATE", X, ignore.case = TRUE)
    translation <- FALSE
    if (length(i2) == 1 && length(i1) == 1) {
        if (i2 > i1) 
            translation <- TRUE
    }
    if (translation) {
        semico <- grep(";", X)
        end <- semico[semico > i2][1]
        x <- X[(i2 + 1):end]
        x <- unlist(strsplit(x, "[,; \t]"))
        x <- x[nzchar(x)]
        TRANS <- matrix(x, ncol = 2, byrow = TRUE)
        TRANS[, 2] <- gsub("['\"]", "", TRANS[, 2])
        nspecies <- dim(TRANS)[1]
        speciesname <- TRANS[, 2]
    }
    
    tree <- X[grep("=",X)]
    tree <- gsub(".*=", "", tree)
    }
    
    if (!translation) {
        speciesname <- species.name(tree[1])
    }
    string <- unlist(strsplit(tree[1], NULL))
    leftpar <- which(string == "(")
    rightpar <- which(string == ")")
    comma <- which(string == ",")
    if (length(leftpar) != length(rightpar)) 
        stop("The number of left parenthesis is NOT equal to the number of right  parenthesis")
    if (length(leftpar) == length(comma)) 
        rooted <- TRUE
    else if (length(leftpar) == (length(comma) - 1)) 
        rooted <- FALSE
    else 
	{
		print("The tree is not a binary tree!")
		rooted <- NA
	}

    z <- list(tree = "", names = "", root = TRUE)
    z$tree <- tree
    z$names <- speciesname
    z$root <- rooted
    z
}
