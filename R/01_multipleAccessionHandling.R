
library(tidyverse)

# 1. Clean gene ids ----------------------------------------
# dependent libraraies stringr

#' Clean gene ids (Optional)
#'
#' Clean gene ids from additional infromations. Both in case singe and multiple gene ids in a cell, this is helpful.
#'
#' @param gene_data A data.frame object of the lfq data.
#' @param colName A character string containing column name of data.frame object of the gene id/accession.
#' @param pattern A character string of regular expression that indicates position of the gene id that has to be extracted.
#' 
#' @export
getCleanAccession <- function(gene_data, colName, pattern) {
  # Validate inputs
  if (!is.data.frame(gene_data)) {
    stop("The input data must be a data frame.")
  }
  if (!colName %in% names(gene_data)) {
    stop("The specified column does not exist in the data frame.")
  }
  if (!is.character(pattern) || length(pattern) != 1) {
    stop("Pattern must be a single character string.")
  }

  # Apply the cleaning function to the specified column
  gene_data <- gene_data %>%
    mutate(!!colName := sapply(.[[colName]], cleanExtractAccession, pattern = pattern))
  
  return(gene_data)
}

#' Clean gene ids (Optional) helper
#' 
#' takes a single input string and pattern, and clean the the cell.
#' 
#' @param input_string A single input string of raw accessions.
#' @param pattern A character string of regular expression that indicates position of the gene id that has to be extracted.
#' 
#' @export
cleanExtractAccession <- function(input_string, pattern) {
  matches <- str_extract_all(input_string, pattern)
  cleaned_result <- unique(unlist(matches))
  paste(cleaned_result, collapse = ",")
}


#' getCleanAccession(f_data, "Accession", "\\b[A-Za-z0-9]+(?=\\|)")


# 2. Get Gene Names --------------------

#' Counts maximum number of gene ids
#'
#' Count maximum number of gene ids or accessions in a given cell having maximum number of gene id possible for a given column.
#'
#' @param gene_data A data.frame object of the lfq data.
#' @param colName A character string containing column name of data.frame object of the gene id/accession.
#' @param delimiter A character string containing a regular expression, separating one gene id from another in a cell of multiple gene id.
#' @param endsWith Logical, if TRUE it considers that delimiter is also in the end of the string in a cell.
#' @return count of maximum number of gene ids in the respective cell of a given column
#' @param verbose_ A Logical input, to show extra information
#' 
#' @export
countAccession <-  function(gene_data, colName, delimiter, endsWith = FALSE, verbose_){
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


#' Accesses information about the gene ids from uniprot (helper)
#' 
#' Once a list of accession or gene id is provided, it try to access all related nformation from UniProt mentioned in the information.
#' 
#' @param accession_list A list of uniprot accessions as a string.
#' @param information A string with all the information (separated with ",") has to be accessed from the UniProt
#' @param verbose_ A Logical input, to show extra information
#' 
#' @export
getGeneUniProt <- function(accession_list, information = "accession,gene_primary,gene_synonym,organism_name,protein_name,sequence", verbose_) {
  message("Processing your accessions, please wait...")
  
  # Initialize progress bar and empty data frame
  progress_bar <- progress::progress_bar$new(total = length(accession_list))
  progress_bar$tick(0) # Show progress bar immediately
  
  combined_data <- data.frame()
  
  for (accession in accession_list) {
    # Fetch data for the current accession
    protein_data <- fetchUniprotData(accession, information)
    
    # If valid data is returned, add it to the combined data
    if (!is.null(protein_data)) {
      protein_data <- protein_data[1, ] # Use only the first row
      combined_data <- rbind(combined_data, as.data.frame(protein_data, row.names = accession))
    }
    progress_bar$tick() # Update progress bar
  }
  
  return(combined_data)
}

#' Accesses information about the gene ids from uniprot (helper)
#' 
#' This takes individual accession, check if the URL is valid and returns the corresponding result
#' 
#' @param accession Individual accession as a string
#' @param information A string
#' @param base_url (default) base link of uniprot
#' @param verbose_ A Logical input, to show extra information
#' 
#' @export
# Function to fetch data for a single accession
fetchUniprotData <- function(accession, information, base_url = "https://rest.uniprot.org/uniprotkb/search?query=accession:", verbose_) {
  # Construct the query URL
  query_url <- paste0(base_url, accession, "&format=tsv&fields=", information)
  query_url <- URLencode(query_url) # Ensure the URL is safe
  
  # Make the API request and handle errors
  response <- tryCatch(
    {
      httr::GET(query_url)
    },
    error = function(e) {
      message("Error during API request for accession: ", accession)
      return(NULL)
    }
  )
  
  # Check if the request was successful
  if (is.null(response) || response$status_code != 200 && verbose_) {
    message("Failed to fetch data for accession: ", accession, 
            ". HTTP status: ", ifelse(is.null(response), "NULL", response$status_code))
    return(NULL)
  }
  
  # Parse the response into a data frame
  protein_data <- tryCatch(
    {
      read.csv(query_url, header = TRUE, sep = "\t")
    },
    error = function(e) {
      if(verbose_){
        message("Error parsing data for accession: ", accession)
      }
      
      return(NULL)
    }
  )
  
  return(protein_data)
}


#' Accessing the merged gene ids
#'
#' There are a lot of the gene id that has been merged with a new gene id, that can be retrieved using this package.
#' 
#' @param Accession A string or a vector of string containing gene ids/ accesssions
#' 
#' @export
getMergedGeneNames <- function(Accession) {
  new_accession_list <- c()
  base_url <- "https://www.uniprot.org/uniprot/"
  query_url <- paste0(base_url, Accession, ".json")
  for (i in query_url) {
    response <- GET(i)
    if (status_code(response) != 200) {
      data <- fromJSON(content(response, as = "text", encoding = "UTF-8"))
      if (!is.null(data$path)) {
        new_accession = sub(".*/", "", data$path)
        
      } else {
        new_accession = NA
      }
    } else {
      new_accession = NA
    }
    new_accession_list = append(new_accession_list, new_accession)
  }
  merged_accession = data.frame(Accession = Accession, new_Accession = new_accession_list)
  return(merged_accession)
}

#' Accesses information about the gene ids from uniprot of a cell having multiple entries
#' 
#' Accesses information about the gene ids from uniprot of a cell having multiple entries
#'
#' @param gene_data A data.frame object, should be the lfq data
#' @param colName A string of column name that contains the gene ids
#' @param delimiter A string of regular expressions through which the gene ids are pasted in a cell
#' @param information A string mentioning all the information that has to be retrieved from UniProt
#' @param countAccession A function to count accessions. Refer ?countAccession() for details.
#' @param getGeneUniProt A function to get the gene informations from UniProt. Refer to ?getGeneUniProt() for details.
#' @param verbose_ A Logical input, to show extra information
#' 
#' @export
getGeneInformationforMultipleAccession <- function(gene_data, colName, delimiter, information, countAccession, getGeneUniProt, verbose_) {
  
  UniprotNames = data.frame(matrix(nrow = 0, ncol = length(unlist(strsplit(information, ","))))) # Dummy data.frame to collect information from uniprot
  colnames(UniprotNames) = unlist(strsplit(information, ",")) # renaming the dummy data.frame colnames with the information header
  
  NA_GeneData = gene_data
  
  # Handling multiple accession number:
  no_accession = sum(grepl("Accession", colnames(gene_data)))
  
  # First layer of retrieving the gene details
  for (n in 1:(no_accession-1)){
    print(n)
    UniprotNames_dummy = getGeneUniProt(as.vector(NA_GeneData[[paste0(colName,".", n)]]), information) #Enriching Gene Information from uniprot using Accession Number
    colnames(UniprotNames_dummy) = str_to_title(unlist(strsplit(information, ","))) #Changing column to more readable names and capitalize first letter
    rownames(UniprotNames_dummy) = NULL #Rownames converted to NULL
    
    UniprotNames = rbind(UniprotNames, UniprotNames_dummy)
    
    ## If Accession.1 UniProt is NA (UniProt --> if NA)
    NA_geneprimary = c()
    
    NA_geneprimary = which(is.na(UniprotNames$Gene_primary)) # This will return entries where Accession.1 had NA, Accession.2, Accession.n
    
    ## Filtering gene_data according to Accession.1 UniProt information is NA (if NA --> Accession.n+1)
    NA_GeneData = c()
    NA_GeneData = gene_data[gene_data[[paste0(colName,".", n)]] %in% UniprotNames$Accession[NA_geneprimary],] 
    
    ## Filtering Accession numbers which are not NA in Acc.2 (Accession.n+1 --> if NOT NA)
    NA_GeneData  = NA_GeneData[!is.na(NA_GeneData[[paste0(colName,".", n+1)]]), ]
    
    if(nrow(NA_GeneData)==0){
      print("Accession retrieved (along with the modified entries).")
      break
    }
    
    ## Dropping UniProt Entries (Accession.n+1 --> else drop entries in UniProt) 
    UniprotNames = UniprotNames %>%
      filter(Accession %in% setdiff(UniprotNames$Accession, NA_GeneData[[paste0("Accession.", n)]]))
    print(n)
  }
  
  # merged genes doesn't get replaced with new genes so for that we have to re retrieve the gene information
  merged_Accession_UniprotNames <-  UniprotNames %>% 
    filter(Protein_name == "merged")
  
  if (nrow(merged_Accession_UniprotNames)>0){
    
    print("Please wait we are replacing merged entries with the new information...")
    
    # enriching "merged" accessions with UniProt
    merged_Accession_UniprotNames_information = getMergedGeneNames(merged_Accession_UniprotNames$Accession) #collecting merged_accession -> new_accession
    
    merged_Accession_UniprotNames_information_dummmy = getGeneUniProt(as.vector(merged_Accession_UniprotNames_information[["new_Accession"]]), information)  #collecting details corresponding the merged accession number
    colnames(merged_Accession_UniprotNames_information_dummmy) = str_to_title(unlist(strsplit(information, ","))) #Changing column to more readable names and capitalize first letter
    rownames(merged_Accession_UniprotNames_information_dummmy) = NULL #Rownames converted to NULL
    
    merged_Accession_UniprotNames_information <- merged_Accession_UniprotNames_information %>% 
      full_join(merged_Accession_UniprotNames_information_dummmy, by = join_by(new_Accession == Accession))
    
    # Replacing merged merged_accession_uniprot with UniProtNam
    UniprotNames <- UniprotNames %>% 
      filter(Protein_name != "merged") # removing the (Protein_names == merged) entries from UniprotNames
    
    UniprotNames <- rbind(UniprotNames, merged_Accession_UniprotNames_information %>% 
                            select(-new_Accession))
  }
  
  return(UniprotNames)
}



#' Populate the cell having multiple entries with single valid gene id
#' 
#' Populate the cell having multiple entries with single valid gene id
#'
#' @param gene_data A data.frame object, should be the lfq data
#' @param colName A string of column name that contains the gene ids
#' @param delimiter A string of regular expressions through which the gene ids are pasted in a cell
#' @param information A string mentioning all the information that has to be retrieved from UniProt
#' @param countAccession A function to count accessions. Refer ?countAccession() for details.
#' @param getGeneUniProt A function to get the gene informations from UniProt. Refer to ?getGeneUniProt() for details.
#' @param getGeneInformationforMultipleAccession A function that works on multiple Accession returns the uniprot Information
#' @param verbose_ A Logical input, to show extra information
#' 
#' @export
populateGeneNamesfromMultipleAccession <- function(gene_data, 
                                                   colName, 
                                                   delimiter, 
                                                   information = "accession,gene_primary,gene_synonym,organism_name,protein_name,sequence", 
                                                   countAccession = countAccession, 
                                                   getGeneUniProt = getGeneUniProt,
                                                   getGeneInformationforMultipleAccession = getGeneInformationforMultipleAccession,
                                                   verbose_ = FALSE){
  
  # Checkpoint for valid inputs
  if(nchar(colName)<1 && !is.data.frame(gene_data) && (colName %in% colnames(gene_data))){
    if(nchar(colName)<1){
      message("Error: Column name can't be empty")
    }
    
    if(!is.data.frame(gene_data)){
      message("Error: gene_data has to be a data.frame object")
    }
    
    if(colName %in% colnames(gene_data)){
      message("Error: Give column name is not in gene_data data.frame")
    }
  
    break
  }
  
  # converting the number of columns
  gene_data = gene_data %>%
    separate_wider_delim(colName,
                         delimiter,
                         names = paste0(
                           colName, ".",
                           1:countAccession(gene_data, colName, delimiter, endsWith = FALSE)
                         ),
                         too_few = "align_start")
  
  # accessing information from UniProt
  UniprotNames <- getGeneInformationforMultipleAccession(gene_data, colName, delimiter, information, countAccession, getGeneUniProt, verbose_)
  
  if(nrow(UniprotNames) == 0){
    message("No uniprot information is retrieved")
    break
  }
  
  # merging uniprot information with the gene_data object
  gene_data = gene_data  %>%
    # converting all the accession as individual entries
    pivot_longer(
      cols = starts_with(paste0(colName, ".")), 
      names_to = paste0(colName, "_no"), #creating extra column "Accession_no" to saving number Accession."1, 2, 3, 4"
      names_prefix = paste0(colName, "."),#to remove "Accession." from Accession.1, Accession.2,... 
      values_to = colName,
      values_drop_na = TRUE
    ) %>%
    relocate(c(colName, paste0(colName, "_no")), .before = 1)
  
  if(verbose_){
    accession = c()
    entries_count = c()
    multiple_entries = c()
    j = 0 # to ignore sum(Accession.9) - sum(accession.10)
    
    for (i in max(gene_data$Accession_no):min(gene_data$Accession_no)){
      accession = append(accession, i) 
      entries_count = ifelse(j>0, 
                             glue("Only {sum(gene_data$Accession_no == i) - sum(gene_data$Accession_no == (i+1))} entries has {i} accessions."),
                             glue("Only {sum(gene_data$Accession_no == i)} entries has {i} accessions."))
      
      multiple_entries = append(multiple_entries, entries_count)
      j = j+1
    }
    
    print(multiple_entries)
  }
  
  gene_data = gene_data %>%
    select(-"Accession_no")
  
  enriched_gene_data =  UniprotNames %>%
    left_join(gene_data, by= join_by(Accession)) #%>% 
    # relocate(c("GeneName.Desc"), .before = "Organism_name") 
  
  enriched_gene_data_01 = enriched_gene_data[!(is.na(enriched_gene_data$Gene_primary) & enriched_gene_data$Protein_name == "deleted"), ]
  
  if(nrow(enriched_gene_data_01) != nrow(enriched_gene_data)){
    removed_gene_list <- enriched_gene_data[(is.na(enriched_gene_data$Gene_primary) & enriched_gene_data$Protein_name == "deleted"), ] %>% 
      pull(colName)
    print(paste0("List of accessions that has been deleted from UniProt(", (nrow(enriched_gene_data) - nrow(enriched_gene_data_01)),"):",paste(unlist(removed_gene_list), collapse = ", ")))
  }
  
  return(enriched_gene_data_01)
}              
  

# column_names start with X fix that
  
