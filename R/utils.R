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
