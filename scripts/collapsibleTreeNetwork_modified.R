# ------------------------------------------------------------------------------------ #
# ----------------------- Code for producing interactive trees ----------------------- #
# ------------------------------------------------------------------------------------ #

require(icd)

# NOTE: The main function in this file, collapsibleTreeNetwork2(), was adapted from the collapsibleTree package
# (version 0.1.6) created by Adeel Khan. The original version of this package is available via CRAN at the 
# following link: https://CRAN.R-project.org/package=collapsibleTree. 

# Expand icd9 codes into nice html format
explainer.zip <- function(df) collapse(c(rep("&emsp;&emsp;&emsp;&emsp;",df["indent"]),"&bull; ",df["icd9_decimal"],": ",df["explainer"],"<br/>"))

expand_html <- function(dat,digits=3,tr){
  # Dat is a data frame or list with three named entries:
  # outcomes (list of associated outcomes)
  # leaf (is that node a leaf) 
  # child (name of the node itself)
  # beta is the coefficient associated with all of the outcomes listed
  frmt <- paste0("%.",digits,"f")
  beta <- sprintf(frmt,exp(dat$beta))
  beta_ciu <- sprintf(frmt,exp(dat$beta_ciu))
  beta_cil <- sprintf(frmt,exp(dat$beta_cil))
  nsamp <- format(dat$nsamp,big.mark=",")
  chapslist <- sapply(icd9_chapters,function(chap) paste(chap[1],"to",chap[2]))
  chapsnames <- names(icd9_chapters)
  subchapslist <- sapply(icd9_sub_chapters,function(chap) paste(chap[1],"to",chap[2]))
  subchapsnames <- names(icd9_sub_chapters)
  if(dat$leaf){ # If this is a leaf node...
    codes <- dat$outcomes[[1]]
    which.chaps <- codes %in% chapslist
    which.subchaps <- codes %in% subchapslist
    which.codes <- !(which.chaps | which.subchaps)
    if(sum(which.chaps | which.subchaps | which.codes) != length(codes)){
      stop("stuffed up somewhere")
    }
    explainer <- character(length=length(codes))
    icd9_decimal <- character(length=length(codes))
    if(sum(which.chaps > 0)){
      explainer[which.chaps] <- sapply(codes[which.chaps], function(code) chapsnames[chapslist == code])
      icd9_decimal[which.chaps] <- codes[which.chaps]
    }
    if(sum(which.subchaps > 0)){
      explainer[which.subchaps] <- sapply(codes[which.subchaps], function(code) subchapsnames[subchapslist == code])
      icd9_decimal[which.subchaps] <- codes[which.subchaps]
    }
    explainer[which.codes] <- sapply(codes[which.codes],icd_explain,condense=F,brief=F)
    icd9_decimal[which.codes] <- as.character(icd_short_to_decimal(codes[which.codes]))
   
    # Now get the shape of the 'tree' for this set of nodes
    outcomes.sort <- order(as.numeric(dat$outcomes[[1]]))
    subtr <- induced_subgraph(tr,dat$outcomes[[1]][outcomes.sort])
    comp <- components(subtr)$membership
    html_out <- character()
    for(co in unique(comp)){
      # Get subtree formed by a connected component
      subtr.comp <- induced_subgraph(subtr,V(subtr)[comp == co])
      n.comp <- length(V(subtr.comp))
      if(n.comp > 1){
        # Find root node
        root.subtr <- V(subtr.comp)[ego_size(subtr.comp,mode="in",mindist=1,order=1)==0]
        # ordering of nodes using data.tree package
        subtr.comp.data <- data.tree:::FromDataFrameNetwork(as.data.frame(get.edgelist(subtr.comp)))
        codes.order <- subtr.comp.data$Get("name")
        codes.extract <- sapply(codes.order,function(code) which(codes==code))
        # Amount of indentation (based on distance from root)
        indent <- as.numeric(distances(subtr.comp,v=root.subtr,to=codes.order,mode="out"))
        # Put it all in a data frame
        paste.df <- data.frame(icd9_decimal=icd9_decimal[codes.extract],explainer=explainer[codes.extract],
                               indent=indent,row.names=NULL,stringsAsFactors = FALSE)
        # collapse
        html_out <- c(html_out,collapse(apply(paste.df,1,explainer.zip)))
      } else {
        html_out <- c(html_out,paste0("&bull; ",icd9_decimal[comp==co],": ",explainer[comp==co],"<br/>"))
      }
    }
    html_out <- paste0("Odds Ratio (95% Credible Interval) = ",beta," (",beta_cil,",",beta_ciu,")<br/>","n = ",nsamp,"<br/>","<p align=\"left\">ICD9 codes:<br/>",collapse(html_out),"</p>")
    
  } else {
    # If this is not a leaf node
    code <- dat$child
    if(code %in% chapslist){
      icd9_decimal <- code
      explainer <- chapsnames[chapslist==code]
    } else if(code %in% subchapslist){
      icd9_decimal <- code
      explainer <- subchapsnames[subchapslist==code]
    } else {
      icd9_decimal <- icd_short_to_decimal(code)
      explainer <- icd_explain(code,condense=F,brief=F)
    }
    html_out <- paste0(icd9_decimal,": ",explainer)
  }
  
  return(html_out)
}

expand_html_sims <- function(dat,digits=3,tr){
  # Dat is a data frame or list with three named entries:
  # outcomes (list of associated outcomes)
  # leaf (is that node a leaf) 
  # child (name of the node itself)
  # beta is the coefficient associated with all of the outcomes listed
  frmt <- paste0("%.",digits,"f")
  beta1 <- sprintf(frmt,exp(dat$beta1))
  beta2 <- sprintf(frmt,exp(dat$beta2))
  chapslist <- sapply(icd9_chapters,function(chap) paste(chap[1],"to",chap[2]))
  chapsnames <- names(icd9_chapters)
  subchapslist <- sapply(icd9_sub_chapters,function(chap) paste(chap[1],"to",chap[2]))
  subchapsnames <- names(icd9_sub_chapters)
  if(dat$leaf){ # If this is a leaf node...
    codes <- dat$outcomes[[1]]
    which.chaps <- codes %in% chapslist
    which.subchaps <- codes %in% subchapslist
    which.codes <- !(which.chaps | which.subchaps)
    if(sum(which.chaps | which.subchaps | which.codes) != length(codes)){
      stop("stuffed up somewhere")
    }
    explainer <- character(length=length(codes))
    icd9_decimal <- character(length=length(codes))
    if(sum(which.chaps > 0)){
      explainer[which.chaps] <- sapply(codes[which.chaps], function(code) chapsnames[chapslist == code])
      icd9_decimal[which.chaps] <- codes[which.chaps]
    }
    if(sum(which.subchaps > 0)){
      explainer[which.subchaps] <- sapply(codes[which.subchaps], function(code) subchapsnames[subchapslist == code])
      icd9_decimal[which.subchaps] <- codes[which.subchaps]
    }
    explainer[which.codes] <- sapply(codes[which.codes],icd_explain,condense=F,brief=F)
    icd9_decimal[which.codes] <- as.character(icd_short_to_decimal(codes[which.codes]))
    
    # Now get the shape of the 'tree' for this set of nodes
    outcomes.sort <- order(as.numeric(dat$outcomes[[1]]))
    subtr <- induced_subgraph(tr,dat$outcomes[[1]][outcomes.sort])
    comp <- components(subtr)$membership
    html_out <- character()
    for(co in unique(comp)){
      # Get subtree formed by a connected component
      subtr.comp <- induced_subgraph(subtr,V(subtr)[comp == co])
      n.comp <- length(V(subtr.comp))
      if(n.comp > 1){
        # Find root node
        root.subtr <- V(subtr.comp)[ego_size(subtr.comp,mode="in",mindist=1,order=1)==0]
        # ordering of nodes using data.tree package
        subtr.comp.data <- data.tree:::FromDataFrameNetwork(as.data.frame(get.edgelist(subtr.comp)))
        codes.order <- subtr.comp.data$Get("name")
        codes.extract <- sapply(codes.order,function(code) which(codes==code))
        # Amount of indentation (based on distance from root)
        indent <- as.numeric(distances(subtr.comp,v=root.subtr,to=codes.order,mode="out"))
        # Put it all in a data frame
        paste.df <- data.frame(icd9_decimal=icd9_decimal[codes.extract],explainer=explainer[codes.extract],
                               indent=indent,row.names=NULL,stringsAsFactors = FALSE)
        # collapse
        html_out <- c(html_out,collapse(apply(paste.df,1,explainer.zip)))
      } else {
        html_out <- c(html_out,paste0("&bull; ",icd9_decimal[comp==co],": ",explainer[comp==co],"<br/>"))
      }
    }
    html_out <- paste0("Scenario 1 Odds Ratio = ",beta1,"<br/>","Scenario 2 Odds Ratio = ",beta2,"<br/>","<p align=\"left\">ICD9 codes:<br/>",collapse(html_out),"</p>")
    
  } else {
    # If this is not a leaf node
    code <- dat$child
    if(code %in% chapslist){
      icd9_decimal <- code
      explainer <- chapsnames[chapslist==code]
    } else if(code %in% subchapslist){
      icd9_decimal <- code
      explainer <- subchapsnames[subchapslist==code]
    } else {
      icd9_decimal <- icd_short_to_decimal(code)
      explainer <- icd_explain(code,condense=F,brief=F)
    }
    html_out <- paste0(icd9_decimal,": ",explainer)
  }
  
  return(html_out)
}


collapsibleTreeNetwork2 <- function (df, inputId = NULL, attribute = "leafCount", aggFun = sum, 
                                     fill = "lightsteelblue", linkLength = NULL, fontSize = 10, 
                                     tooltip = TRUE, tooltipHtml = NULL, nodeSize = NULL, nodeSizeAggFun = identity,
                                     nodeSizeScaleFun = stats::median,
                                     collapsed = TRUE, zoomable = TRUE, width = NULL, height = NULL) 
{
  nodeAttr <- c("leafCount", "count")
  if (!is.data.frame(df)) 
    stop("df must be a data frame")
  if (sum(is.na(df[, 1])) != 1) 
    stop("there must be 1 NA for root in the first column")
  if (!is.character(fill)) 
    stop("fill must be a either a color or column name")
  if (!(attribute %in% c(colnames(df), nodeAttr))) 
    stop("attribute column name is incorrect")
  if (!is.null(tooltipHtml)) 
    if (!(tooltipHtml %in% colnames(df))) 
      stop("tooltipHtml column name is incorrect")
  if (!is.null(nodeSize)) 
    if (!(nodeSize %in% c(colnames(df), nodeAttr))) 
      stop("nodeSize column name is incorrect")
  root <- df[is.na(df[, 1]), ]
  tree <- df[!is.na(df[, 1]), ]
  node <- data.tree::FromDataFrameNetwork(tree)
  rootAttr <- root[-(1:2)]
  Map(function(value, name) node[[name]] <- value, rootAttr, 
      names(rootAttr))
  leftMargin <- nchar(node$name)
  rightLabelVector <- node$Get("name", filterFun = function(x) x$level == 
                                 node$height)
  rightMargin <- max(sapply(rightLabelVector, nchar))
  options <- list(hierarchy = 1:node$height, input = inputId, 
                  attribute = attribute, linkLength = linkLength, fontSize = fontSize, 
                  tooltip = tooltip, collapsed = collapsed, zoomable = zoomable, 
                  margin = list(top = 20, bottom = 20, left = (leftMargin * 
                                                                 fontSize/2) + 25, right = (rightMargin * fontSize/2) + 
                                  25))
  jsonFields <- NULL
  if (fill %in% colnames(df)) {
    node$Do(function(x) x$fill <- x[[fill]])
    jsonFields <- c(jsonFields, "fill")
  } else {
    options$fill <- fill
  }
  if (tooltip & is.null(tooltipHtml)) {
    if (is.numeric(df[[attribute]]) & substitute(aggFun) != 
        "identity") {
      t <- data.tree::Traverse(node, "pre-order")
      data.tree::Do(t, function(x) {
        x$WeightOfNode <- data.tree::Aggregate(x, attribute, 
                                               aggFun)
        x$WeightOfNode <- prettyNum(x$WeightOfNode, big.mark = ",", 
                                    digits = 3, scientific = FALSE)
      })
    }
    else {
      node$Do(function(x) x$WeightOfNode <- x[[attribute]])
    }
    jsonFields <- c(jsonFields, "WeightOfNode")
  }
  if (tooltip & !is.null(tooltipHtml)) {
    node$Do(function(x) x$tooltip <- x[[tooltipHtml]])
    jsonFields <- c(jsonFields, "tooltip")
  }
  if (!is.null(nodeSize)) {
    scaleFactor <- 10/data.tree::Aggregate(node, nodeSize, nodeSizeScaleFun)
    t <- data.tree::Traverse(node, "pre-order")
    if ( substitute(nodeSizeAggFun) != "identity" | nodeSize %in% nodeAttr) {
      data.tree::Do(t, function(x) {
        x$SizeOfNode <- data.tree::Aggregate(x, nodeSize, 
                                             sum)
        x$SizeOfNode <- round(sqrt(x$SizeOfNode * scaleFactor) * 
                                pi, 2)
      })
    } else {
      #node$Do(function(x) x$SizeOfNode <- x[[nodeSize]])
      data.tree::Do(t, function(x) {
        x$SizeOfNode <- x[[nodeSize]]
        x$SizeOfNode <- round(sqrt(x$SizeOfNode * scaleFactor) * 
                                pi, 2)
      })
    }
    options$margin$left <- options$margin$left + node$SizeOfNode - 
      10
    jsonFields <- c(jsonFields, "SizeOfNode")
  }
  if (is.null(jsonFields)) 
    jsonFields <- NA
  data <- data.tree::ToListExplicit(node, unname = TRUE, keepOnly = jsonFields)
  x <- list(data = data, options = options)
  htmlwidgets::createWidget("collapsibleTree", x, width = width, 
                            height = height, htmlwidgets::sizingPolicy(viewer.padding = 0))
}

