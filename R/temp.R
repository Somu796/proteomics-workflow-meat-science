
# from 01
#' Counts maximum number of gene ids
#'
#' Count maximum number of gene ids or accessions in a given cell having maximum number of gene id possible for a given column.
#'
#' @param gene_data A data.frame object of the lfq data.
#' @param colName A character string containing column name of data.frame object of the gene id/accession.
#' @param delimiter A character string containing a regular expression, separating one gene id from another in a cell of multiple gene id.
#' @param endsWith Logical, if TRUE it considers that delimiter is also in the end of the string in a cell.
#' @return count of maximum number of gene ids in the respective cell of a given column
#' @export
countAccession <-  function(gene_data, colName, delimiter, endsWith = FALSE){
  max_count <- 0
  for (i in 1:length(gene_data[[colName]])){
    # delimiter <- glue("\\{delimiter}")
    matches <- gregexpr(delimiter, gene_data[[colName]][i])[[1]]
    count <- ifelse(matches[1] == -1, 0, length(matches))

    if (endsWith){
      # print(count)
      if (max_count < count){
        max_count <- count
        index <- i
      }
    }else{
      count <- count + 1
      # print(count)
      if (max_count < count){
        max_count <- count
        index <- i
      }
    }

  }
  return(max_count)
}

#' Extract gene id from additional information(MODIFY).
#' If gene ids has additional information tag along with it, this function helps to clean or extract the gene id. eg: following pattern_from = "\\|.*" and patten_to = "", the later part of "A0A3Q1M2L8|A0A3Q1M2L8_BOVIN"   returns "A0A3Q1M2L8".
#'
#' @param gene_data A data.frame object of the lfq data.
#' @param colName A character string containing column name of data.frame object of the gene id/accession.
#' @param pattern_from A character string of regular expression that has to be replaced with
getAccession <- function(gene_data, ColName, pattern_from = "\\|.*", pattern_to = "", verbose = FALSE){
  accessionNo = c()
  for (i in 1:nrow(gene_data)){
    accessionNo = append(accessionNo, sub(pattern_from, pattern_to, gene_data[i,ColName]))
  }
  if (sum(is.na(accessionNo)) == 0 & verbose) print("No Na present in accession Numbers") else print("NA present in Accession Numbers")
  return(accessionNo)
}

#####################



#' Accesses information about the gene ids from uniprot (MODIFY)
#' 
#' Once a list of accession or gene id is provided, it try to access all related nformation from UniProt mentioned in the information. Dependent on httr
#' 
#' @param accession_list A list of uniprot accessions as a string.
#' @param information A string with all the information (separated with ",") has to be accessed from the UniProt
#' 
#' @export

getGeneUniProt <- function(accession_list,  information)
{
  
  message("Please wait we are processing your accessions ...")
  pb <- progress::progress_bar$new(total = length(accession_list))
  
  baseUrl <- "https://rest.uniprot.org/uniprotkb/search?query=accession:"
  
  ProteinInfoParsed_total = data.frame()
  for (ProteinAcc in accession_list)
  {
    #to see if Request == 200 or not
    Request <- tryCatch(
      {
        GET(paste0(baseUrl , ProteinAcc,"&format=tsv"))
      },error = function(cond)
      {
        message("Internet connection problem occurs and the function will return the original error")
        message(cond)
      }
    ) 
    #this link return information in tab formate (format = tab)
    ProteinName_url <- paste0(ProteinAcc,"&format=tsv&fields=",information)
    RequestUrl <- paste0(baseUrl , ProteinName_url)
    RequestUrl <- URLencode(RequestUrl)
    if (length(Request) == 0)
    {
      message("Internet connection problem occurs")
      return()
    }
    if (Request$status_code == 200){
      # parse the information in DataFrame
      ProteinDataTable <- tryCatch(read.csv(RequestUrl, header = TRUE, sep = '\t'), error=function(e) NULL)
      if (!is.null(ProteinDataTable))
      {
        ProteinDataTable <- ProteinDataTable[1,]
        ProteinInfoParsed <- as.data.frame(ProteinDataTable,row.names = ProteinAcc)
        # add Dataframes together if more than one accession
        ProteinInfoParsed_total <- rbind(ProteinInfoParsed_total, ProteinInfoParsed)
      }
      
    }else {
      HandleBadRequests(Request$status_code)
    }
    pb$tick()
    
  }
  return(ProteinInfoParsed_total)
}

#################

