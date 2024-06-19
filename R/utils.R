xy2ThetaR<-function(x){
  r<-apply(x,1,function(a) sqrt(a[1]^2+a[2]^2))
  theta <- -apply(x,1,function(a) atan(a[2]/a[1]))
  return(data.frame( Theta=theta,R=r))

}

draw.polygon<-function(){
  p<-locator(1)
  pret<-p
  while (p$x>0){
    p2<-locator(1)
    if (p2$x>0){
      segments(p$x,p$y,p2$x,p2$y, col="red")
      pret$x<-c(pret$x,p2$x)
      pret$y<-c(pret$y,p2$y)
    }
    p<-p2
  }
  segments(pret$x[1],pret$y[1],tail(pret$x,1),tail(pret$y,1), col="red")
  return(pret)
}

parse_api_url <- function(url){
  ## get protocol (default = "https://")
  protocolless_url <- gsub("(^http://|^https://)(.*)$", "\\2",url)
  brapi_protocol <- gsub("(^http://|^https://)(.*)$", "\\1",url)
  brapi_protocol <- ifelse(brapi_protocol == url, "https://", brapi_protocol)

  ## get base url and port (default = 443)
  db_split <- strsplit(gsub("([^/]*).*", "\\1",protocolless_url), ":")
  brapi_db <- db_split[[1]][1]
  brapi_port <- ifelse(is.na(db_split[[1]][2]),80,as.numeric(db_split[[1]][2]))

  ## brapi api path (default = "/")
  brapi_apipath <- ifelse(grepl("/.*",protocolless_url),gsub("[^/]*/(.*)", "\\1", protocolless_url),"/")

  return(
    list(
      brapi_protocol = brapi_protocol,
      brapi_db = brapi_db,
      brapi_port = brapi_port,
      brapi_apipath = brapi_apipath
    )
  )
}

# brapi_get_variants <- function (con = NULL, referenceDbId = "", start="", variantDbId="", variantSetDbId="", page = 0, pageSize = 1000) {
#   usedArgs <- brapirv2:::brapi_usedArgs(origValues = FALSE)
#   brapirv2::brapi_checkCon(con = usedArgs[["con"]], verbose = FALSE)
#   brapirv2:::brapi_checkArgs(usedArgs, reqArgs = "")
#   callurl <- brapirv2:::brapi_GET_callURL(usedArgs = usedArgs,
#                                           callPath = "/variants", reqArgs = "", packageName = "BrAPI-Genotyping",
#                                           callVersion = 2)
#   try({
#     resp <- brapirv2:::brapi_GET(url = callurl, usedArgs = usedArgs)
#     cont <- httr::content(x = resp, as = "text", encoding = "UTF-8")
#     out <- brapirv2:::brapi_result2df(cont, usedArgs)
#   })
#   class(out) <- c(class(out), "brapi_get_variants")
#   brapirv2:::brapi_serverinfo_metadata(cont)
#   return(out)
# }
# brapi_get_variantsets <- function (con = NULL, studyDbId = "", studyName = "", variantDbId = "", variantSetDbId = "", page = 0, pageSize = 1000) {
#   usedArgs <- brapirv2:::brapi_usedArgs(origValues = FALSE)
#   brapirv2::brapi_checkCon(con = usedArgs[["con"]], verbose = FALSE)
#   brapirv2:::brapi_checkArgs(usedArgs, reqArgs = "")
#   callurl <- brapirv2:::brapi_GET_callURL(usedArgs = usedArgs,
#                                           callPath = "/variantsets", reqArgs = "", packageName = "BrAPI-Genotyping",
#                                           callVersion = 2)
#   try({
#     resp <- brapirv2:::brapi_GET(url = callurl, usedArgs = usedArgs)
#     cont <- httr::content(x = resp, as = "text", encoding = "UTF-8")
#     out <- brapirv2:::brapi_result2df(cont, usedArgs)
#   })
#   class(out) <- c(class(out), "brapi_get_variantsets")
#   brapirv2:::brapi_serverinfo_metadata(cont)
#   return(out)
# }
#
# brapi_get_calls <- function (con = NULL, variantDbId = "", variantSetDbId = "",variantName = "", page = 0, pageSize = 1000) {
#   usedArgs <- brapirv2:::brapi_usedArgs(origValues = FALSE)
#   brapirv2::brapi_checkCon(con = usedArgs[["con"]], verbose = FALSE)
#   brapirv2:::brapi_checkArgs(usedArgs, reqArgs = "")
#   callurl <- brapirv2:::brapi_GET_callURL(usedArgs = usedArgs,
#                                           callPath = "/calls", reqArgs = "", packageName = "BrAPI-Genotyping",
#                                           callVersion = 2)
#   try({
#     resp <- brapirv2:::brapi_GET(url = callurl, usedArgs = usedArgs)
#     cont <- httr::content(x = resp, as = "text", encoding = "UTF-8")
#     out <- brapirv2:::brapi_result2df(cont, usedArgs)
#   })
#   class(out) <- c(class(out), "brapi_get_calls")
#   brapirv2:::brapi_serverinfo_metadata(cont)
#   return(out)
# }
#
# brapi_get_calls <- function (con = NULL, variantDbId = "", variantSetDbId = "",variantName = "", page = 0, pageSize = 1000) {
#   usedArgs <- brapirv2:::brapi_usedArgs(origValues = FALSE)
#   brapirv2::brapi_checkCon(con = usedArgs[["con"]], verbose = FALSE)
#   brapirv2:::brapi_checkArgs(usedArgs, reqArgs = "")
#   callurl <- brapirv2:::brapi_GET_callURL(usedArgs = usedArgs,
#                                           callPath = "/calls", reqArgs = "", packageName = "BrAPI-Genotyping",
#                                           callVersion = 2)
#   try({
#     resp <- brapirv2:::brapi_GET(url = callurl, usedArgs = usedArgs)
#     cont <- httr::content(x = resp, as = "text", encoding = "UTF-8")
#     out <- brapirv2:::brapi_result2df(cont, usedArgs)
#   })
#   class(out) <- c(class(out), "brapi_get_calls")
#   brapirv2:::brapi_serverinfo_metadata(cont)
#   return(out)
# }
# brapi_get_variants_variantDbId_calls <- function (con = NULL, variantDbId = "", page = 0, pageSize = 1000) {
#   usedArgs <- brapirv2:::brapi_usedArgs(origValues = FALSE)
#   brapirv2::brapi_checkCon(con = usedArgs[["con"]], verbose = FALSE)
#   brapirv2:::brapi_checkArgs(usedArgs, reqArgs = "")
#   callurl <- brapirv2:::brapi_GET_callURL(usedArgs = usedArgs,
#                                           callPath = "/variants/{variantDbId}/call", reqArgs = "variantDbId", packageName = "BrAPI-Genotyping",
#                                           callVersion = 2)
#   try({
#     resp <- brapirv2:::brapi_GET(url = callurl, usedArgs = usedArgs)
#     cont <- httr::content(x = resp, as = "text", encoding = "UTF-8")
#     out <- brapirv2:::brapi_result2df(cont, usedArgs)
#   })
#   class(out) <- c(class(out), "brapi_get_variants_variantDbId_calls")
#   brapirv2:::brapi_serverinfo_metadata(cont)
#   return(out)
# }
#
#
# bmsapi_Get_sample_list_search<-function(con, crop, programUUID, exactMatch=F, searchString="%25"){
#   url<-paste0(bmscon_geturl(con), "/crops/",crop,"/sample-lists/search?programUUID=",programUUID,"&exactMatch=",tolower(as.character(exactMatch)),"&searchString=",searchString)
#   #httr::with_config(
#   #  config = config(ssl_verifypeer = FALSE, ssl_verifyhost = FALSE),
#   #  req<-curl_fetch_memory(url,handle = h2)
#   #)
#   req <- httr::GET(url = url,
#                    httr::timeout(25),
#                    httr::add_headers(Authorization = paste("Bearer", con$token),`X-Auth-Token`=con$token))
#   res<-jsonlite::fromJSON(rawToChar(req$content))
#   return(data.table(res))
# }

brapi_get_calls_fast <- function(con = NULL,
                                    callSetDbId = '',
                                    variantDbId = '',
                                    variantSetDbId = '',
                                    expandHomozygotes = NA,
                                    unknownString = '',
                                    sepPhased = '',
                                    sepUnphased = '',
                                    pageToken = '',
                                    page=0,
                                    pageSize = 1000) {
  ## Create a list of used arguments
  usedArgs <- brapirv2:::brapi_usedArgs(origValues = FALSE)
  if (exists(usedArgs[["sepPhased"]]) && usedArgs[["sepPhased"]] %in% c("|", "/")) {
    usedArgs[["sepPhased"]] <- paste0("%",
                                      toupper(charToRaw(usedArgs[["sepPhased"]])))
  }
  if (exists(usedArgs[["sepUnphased"]]) && usedArgs[["sepUnphased"]] %in% c("|", "/")) {
    usedArgs[["sepUnphased"]] <- paste0("%",
                                        toupper(charToRaw(usedArgs[["sepUnphased"]])))
  }
  if (exists(usedArgs[["unknownString"]]) && usedArgs[["unknownString"]] %in% c("|", "/")) {
    usedArgs[["unknownString"]] <- paste0("%",
                                          toupper(
                                            charToRaw(
                                              usedArgs[["unknownString"]])))
  }
  ## Check if BrAPI server can be reached given the connection details
  brapi_checkCon(con = usedArgs[["con"]], verbose = FALSE)
  ## Check validity of used and required arguments
  brapirv2:::brapi_checkArgs(usedArgs, reqArgs = "")
  ## Obtain the call url
  callurl <- brapirv2:::brapi_GET_callURL(usedArgs = usedArgs,
                                          callPath = "/calls",
                                          reqArgs = "",
                                          packageName = "BrAPI-Genotyping",
                                          callVersion = 2.0)

  try({
    ## Make the call and receive the response
    resp <- brapirv2:::brapi_GET(url = callurl, usedArgs = usedArgs)
    ## Extract the content from the response object in human readable form
    cont <- httr::content(x = resp, as = "text", encoding = "UTF-8")
    ## Convert the content object into a data.frame
    out <- data.table(tidyr::unnest(jsonlite::fromJSON(cont)$result$data,cols = "genotypeMetadata", names_sep = "."))
  })
  ## Set class of output
  class(out) <- c(class(out), "brapi_get_calls")
  ## Show pagination information from metadata
  brapirv2:::brapi_serverinfo_metadata(cont)
  return(out)
}


brapi_post_search_callsets_fast <- function(con = NULL,
                                       callSetDbIds = '',
                                       callSetNames = '',
                                       germplasmDbIds = '',
                                       germplasmNames = '',
                                       page = 0,
                                       pageSize = 1000,
                                       sampleDbIds = '',
                                       sampleNames = '',
                                       variantSetDbIds = '') {
  ## Create a list of used arguments
  usedArgs <- brapirv2:::brapi_usedArgs(origValues = FALSE)
  ## Check if BrAPI server can be reached given the connection details
  brapi_checkCon(con = usedArgs[["con"]], verbose = FALSE)
  ## Check validity of used and required arguments
  brapirv2:::brapi_checkArgs(usedArgs, reqArgs = "")
  ## Obtain the call url
  callurl <- brapirv2:::brapi_POST_callURL(usedArgs = usedArgs,
                                           callPath = "/search/callsets",
                                           reqArgs = "",
                                           packageName = "BrAPI-Genotyping",
                                           callVersion = 2.0)
  ## Build the Body
  callbody <- brapirv2:::brapi_POST_callBody(usedArgs = usedArgs,
                                             reqArgs = "")

  try({
    ## Make the call and receive the response
    resp <- brapirv2:::brapi_POST(url = callurl, body = callbody, usedArgs = usedArgs)
    ## Message about call status
    if (httr::status_code(resp) == 200) {
      message(paste0("Immediate Response.", "\n"))
    } else if (httr::status_code(resp) == 202) {
      message(paste0("Saved or Asynchronous Response has provided a searchResultsDbId.", "\n"))
      message(paste0("Use the GET /search/callsets/{searchResultsDbId} call to retrieve the paginated output.", "\n"))
    } else {
      stop(paste0("The POST /search/callsets call resulted in Server status, ", httr::http_status(resp)[["message"]]))
    }
    ## Extract the content from the response object in human readable form
    cont <- httr::content(x = resp, as = "text", encoding = "UTF-8")
    ## Convert the content object into a data.frame
    out <- data.table(tidyr::unnest(jsonlite::fromJSON(cont)$result$data, cols = "externalReferences", names_sep = "."))
  })
  ## Set class of output
  class(out) <- c(class(out), "brapi_post_search_callsets")
  ## Show pagination information from metadata
  brapirv2:::brapi_serverinfo_metadata(cont)
  return(out)
}

brapi_post_search_samples_fast <- function(con = NULL,
                                      commonCropNames='',
                                      externalReferenceIds='',
                                      externalReferenceSources='',
                                      germplasmDbIds='',
                                      germplasmNames='',
                                      observationUnitDbIds='',
                                      page=0,
                                      pageSize=1000,
                                      plateDbIds='',
                                      plateNames='',
                                      programDbIds='',
                                      programNames='',
                                      sampleDbIds='',
                                      sampleGroupDbIds='',
                                      sampleNames='',
                                      studyDbIds='',
                                      studyNames='',
                                      trialDbIds='',
                                      trialNames='') {
  ## Create a list of used arguments
  usedArgs <- brapirv2:::brapi_usedArgs(origValues = FALSE)
  ## Check if BrAPI server can be reached given the connection details
  brapi_checkCon(con = usedArgs[["con"]], verbose = FALSE)
  ## Check validity of used and required arguments
  brapirv2:::brapi_checkArgs(usedArgs, reqArgs = "")
  ## Obtain the call url
  callurl <- brapirv2:::brapi_POST_callURL(usedArgs = usedArgs,
                                           callPath = "/search/samples",
                                           reqArgs = "",
                                           packageName = "BrAPI-Genotyping",
                                           callVersion = 2.0)
  ## Build the Body
  callbody <- brapirv2:::brapi_POST_callBody(usedArgs = usedArgs,
                                             reqArgs = "")

  try({
    ## Make the call and receive the response
    resp <- brapirv2:::brapi_POST(url = callurl, body = callbody, usedArgs = usedArgs)
    ## Message about call status
    if (httr::status_code(resp) == 200) {
      message(paste0("Immediate Response.", "\n"))
    } else if (httr::status_code(resp) == 202) {
      message(paste0("Saved or Asynchronous Response has provided a searchResultsDbId.", "\n"))
      message(paste0("Use the GET /search/samples/{searchResultsDbId} call to retrieve the paginated output.", "\n"))
    } else {
      stop(paste0("The POST /search/samples call resulted in Server status, ", httr::http_status(resp)[["message"]]))
    }
    ## Extract the content from the response object in human readable form
    cont <- httr::content(x = resp, as = "text", encoding = "UTF-8")
    ## Convert the content object into a data.frame
    out <- data.table(tidyr::unnest(jsonlite::fromJSON(cont)$result$data,cols = "externalReferences", names_sep = "."))
  })
  ## Set class of output
  class(out) <- c(class(out), "brapi_post_search_samples")
  ## Show pagination information from metadata
  brapirv2:::brapi_serverinfo_metadata(cont)
  return(out)
}
