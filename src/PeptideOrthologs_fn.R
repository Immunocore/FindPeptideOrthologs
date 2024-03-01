# Functions -----------------------------------------------------------

#' Write to a log file
writeLog <- function(to_write,peptide){
  fname <- paste0("log/",peptide,".txt")
  write(x = paste0(to_write," - ",Sys.time()),file = fname,append = file.exists(fname)) 
}


#' Grep a genome for a peptide string
#' 
#' @param peps_in vector of peptides to search for
#' @param fname_subject file name of the genome to search
GrepAGenome <- function(peps_in,fname_subject){
  sysCall <- paste0("grep -nE '",paste0(peps_in,collapse="|"),"' ",fname_subject)
  
  a <- system(sysCall,intern=TRUE)
  if(is.null(attr(a,"status"))){
    return(a)
  }else{
    return(0)
  }
}
PrecisionGrep <- function(strng,trgt){
  sysCall <- paste0("echo '",trgt,"' | grep -aob '",strng,"'")
  out <- system(sysCall,intern=TRUE)
  if(length(out)>0){
    out_split <- t(sapply(out,function(x) unlist(strsplit(x,split=":"))))
    return(data.frame(pep = out_split[,2],pos = out_split[,1],stringsAsFactors = FALSE))
  }else{
    return(data.frame(pep = NA,pos = NA,stringsAsFactors = FALSE))
  }
}


#' Make wildcarded peptides
#'
#' @param peptide
#' @param n_mismatches
#' @return Vector of wildcarded peptide
make_wildcarded_peptide <- function(peptide,n_mismatches,min_mismatches = 1){
  peptide_split <- unlist(strsplit(peptide,split=""))
  peptide_length <- nchar(peptide)
  
  # Generate a list of wildcarded versions of the target peptide
  wildcardedPeps_list <- vector(mode="list",length=n_mismatches)
  for(i in min_mismatches:n_mismatches){
    modInd <- combn(x = 1:peptide_length,m = i,simplify = FALSE)
    wildcardedPeps_list[[i]] <- 
      lapply(modInd,function(x) {
        modified <- peptide_split
        modified[x] <- "."
        return(paste0(modified,collapse=""))
      })
  }
  wildcardedPeps <- unlist(wildcardedPeps_list)
  return(wildcardedPeps)
}


BloSumR <- function(pep_in,subject,root_score,outType="both"){
  
  split_subject <- unlist(strsplit(subject,split=""),use.names = FALSE)
  split_pep <- unlist(strsplit(pep_in,split=""),use.names = FALSE) 
  n_mismatches <- sum(split_pep != split_subject)
  
  pep_score <- 
    sum(
      BLOSUM62[
        cbind(fmatch(x = split_pep,table = blosum_names),
              fmatch(x = split_subject,table = blosum_names))
      ]
    )
  
  BLOSUM_prop <- pep_score / root_score
  if(outType=="both"){
    return(list(BLOSUM_prop = BLOSUM_prop,n_mismatches = n_mismatches))  
  }else{
    return(BLOSUM_prop=BLOSUM_prop)
  }
}


#' Core function to find orthologs
#' 
#' @param peptide Peptide to search for
#' @param fname_genome File name of reference genome to search
#' @param n_mismatches Number of mismatches that we can about
FindOrthologs <- function(peptide,fname_genome,n_mismatches = 4){
  # peptide <- "LLTQIGCTL";nMismatches <- n_mismatches;fname_genome <- fname_human
  
  tmp_fname <- paste0("tmp/tmp",round(runif(n = 1,min=10000000,max=89999999)),".fasta")
  
  splitPep <- unlist(strsplit(peptide,split=""))
  root_score <- 
    sum(
      BLOSUM62[cbind(fmatch(x = splitPep,table = blosum_names),
                     fmatch(x = splitPep,table = blosum_names))]
    )
  pepLen <- length(splitPep)
  
  writeLog("Making wildcarded peptides",peptide)
  wildcarded_peptides <- make_wildcarded_peptide(peptide,n_mismatches = n_mismatches)
  
  # Run the genome grep once, return the line numbers that contain any one of the peptides
  writeLog("Grepping wildcarded peptides against whole proteome",peptide) 
  matched_genome_lines <- GrepAGenome(peps_in = wildcarded_peptides,fname_subject = fname_genome)
  headerIndx <- grep(">",matched_genome_lines)
  if(length(headerIndx)>0){
    matched_genome_lines <- matched_genome_lines[-headerIndx]  
  }
  writeLines(matched_genome_lines,tmp_fname)
  
  # Run it again on just the lines that have matches, one peptide at a time
  writeLog("Grepping matched wildcarded peptides against subset of proteome",peptide) 
  matches2_tmp <- lapply(wildcarded_peptides,GrepAGenome,fname_subject = tmp_fname)
  matches2 <- 
    unlist(lapply(1:length(matches2_tmp),function(x) 
      if(matches2_tmp[[x]][1]!=0){
        sapply(matches2_tmp[[x]],function(y) 
          paste0(wildcarded_peptides[x],":",y)
        )
      }else{
        character(0)
      }
    ))
  writeLog("Collating outputs",peptide)
  
  match_tbl <- t(sapply(matches2,function(x) unlist(strsplit(x,split=":"))))
  match_outputs <- 
    do.call(rbind,
            # for(x in c(1:dim(match_tbl)[1])){
            lapply(1:dim(match_tbl)[1],function(x){ #
              grepout <- PrecisionGrep(strng = match_tbl[x,1],trgt = match_tbl[x,4])
              nFound <- dim(grepout)[1]
              ProtIndx <- rep(as.numeric(match_tbl[x,3]),nFound)
              return(cbind(grepout,ProtIndx=ProtIndx))
            }))
  rownames(match_outputs) <- NULL
  
  # Add gene/protein information back in
  writeLog("Getting gene info",peptide)
  
  gene_header <- genome_raw[as.numeric(match_outputs$ProtIndx) - 1]
  gene_id <- sapply(gene_header,function(z) unlist(strsplit(z,split="\\|"))[2])
  gene_name_raw <- sapply(gene_header,function(z) unlist(strsplit(z,split="\\|"))[3])
  gene_name <- sapply(gene_name_raw,function(x) gsub(pattern = "_HUMAN",replacement = "",x = x))
  ForTable <- cbind(gene_id,gene_name)
  
  # Add BLOSUM scores
  blosum_output_raw <- t(sapply(match_outputs[,1],function(x) unlist(BloSumR(pep_in = as.character(x),
                                                                             subject = peptide,
                                                                             root_score = root_score))))
  rownames(ForTable) <- NULL
  rownames(blosum_output_raw) <- NULL
  
  
  # Combine everything to a final output
  final_output <- cbind(ForTable,match_outputs[,-3],blosum_output_raw)
  
  final_output_pasted <- apply(final_output,1,paste0,collapse="")
  Duplicates <- which(duplicated(final_output_pasted))
  if(length(Duplicates)>0){
    final_output <- final_output[-Duplicates,]  
  }
  
  rownames(final_output) <- NULL
  writeLog("Returning final dataframe",peptide)
  return(as.data.frame(final_output))
}
