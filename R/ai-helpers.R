## =========================================================================
## AI / LLM helpers (module description, image generation)
## =========================================================================

#' List available local Ollama models
#' @export
ai.get_ollama_models <- function(models = NULL,
                                 size = NULL) {

  available.models <- system("ollama list | tail -n +2 | cut -d' ' -f 1", intern=TRUE)

  models.sizes <- system("ollama list | tail -n +2 | tr -s ' ' | cut -d' ' -f 3", intern=TRUE)
  models.sizes <- as.numeric(models.sizes)
  models.sizes <- ifelse(models.sizes < 100, models.sizes, models.sizes/1000)
  names(models.sizes) <- available.models

  if (!is.null(models) && !any(models == "OLLAMA_MODELS")) {
    available.models <- intersect(models,available.models)
  }

  msize <- models.sizes[available.models] 
  if (!is.null(size) && size == "S") {
    available.models <- available.models[which(msize <= 3)]
  }

  if (!is.null(size) && size == "M") {
    sel <- which( msize > 3 & msize <= 6)
    available.models <- available.models[sel]
  }

  if (!is.null(size) && size == "L") {
    msize <- models.sizes[available.models] 
    available.models <- available.models[which(msize > 6)]
  }
  
  return(available.models)

}

#' Ask an LLM a question, dispatching to the configured engine
#' @export
ai.ask <- function(question,
                   model,
                   engine = c("ellmer", "tidyprompt")[2]) {

  if (model == "ellmer" && grepl("grok", model)) model <- "tidyprompt"

  if (! engine %in% c("ellmer", "tidyprompt"))
    stop("[WGCNAplus::ai.ask] Error. 'engine' must be 'ellmer' or 'tidyprompt'")
  
  if (engine == "ellmer")
    resp <- ai.ask_ellmer(question = question, model = model, prompt = NULL) 

  if (engine == "tidyprompt")
    resp <- ai.ask_tidyprompt(question = question, model = model) 

  return(resp)

}

#' Ask an LLM a question using the ellmer package
#' @export
ai.ask_ellmer <- function(question,
                          model = DEFAULT_LLM,
                          prompt = NULL) {

  chat <- NULL

  if (inherits(model, "Chat")) {
    chat <- model
  } else if (is.character(model)) {
    if (model %in% ai.get_ollama_models() || grepl("^ollama:", model) ) {
      model1 <- sub("^ollama:", "", model)
      chat <- ellmer::chat_ollama(model = model1, system_prompt = prompt)
    } else if (grepl("^gpt|^openai:",model) && Sys.getenv("OPENAI_API_KEY") != "") {
      message("warning: using remote GPT model:", model)
      model1 <- sub("^openai:", "", model)
      key <- Sys.getenv("OPENAI_API_KEY")
      chat <- ellmer::chat_openai(model = model1, system_prompt = prompt, api_key = key)
    } else if (grepl("^grok|^xai:",model) && Sys.getenv("XAI_API_KEY") != "") {
      model1 <- sub("^xai:","",model)
      key <- Sys.getenv("XAI_API_KEY")
      chat <- ellmer::chat_openai(model = model1, system_prompt = prompt,
        api_key = key, base_url = "https://api.x.ai/v1/")
    } else if (grepl("^groq:",model) && Sys.getenv("GROQ_API_KEY") != "") {
      model1 <- sub("groq:", "", model)
      key <- Sys.getenv("GROQ_API_KEY")
      chat <- ellmer::chat_groq(model = model1, system_prompt = prompt, api_key = key)
    } else if (grepl("^gemini|^google:",model) && Sys.getenv("GEMINI_API_KEY") != "") {
      model1 <- sub("^google:","",model)
      key <- Sys.getenv("GEMINI_API_KEY")
      chat <- ellmer::chat_google_gemini(model = model1, system_prompt = prompt, api_key = key)
    }
  }

  if (is.null(chat)) {
    message("ERROR. could not create model ", model)
    return(NULL)
  }

  . <- chat$chat(question, echo = FALSE)

  chat$last_turn()@text

}

#' Ask an LLM a question using the tidyprompt package
#' @export
ai.ask_tidyprompt <- function(question,
                              model,
                              verbose = 0) {

  llm <- NULL
  if (model %in% ai.get_ollama_models() || grepl("^ollama:", model) ) {
    model1 <- sub("^ollama:", "", model)
    prms <- list(model = model1)
    llm <- tidyprompt::llm_provider_ollama(parameters = prms)
  } else if (grepl("^remote:", model) ) {
    remotesrv <- Sys.getenv("OLLAMA_REMOTE")
    if (remotesrv == "") message("error: please set OLLAMA_REMOTE")
    if (remotesrv != "") {
      model1 <- sub("^remote:", "", model)    
      if (verbose > 0) {
        message("connecting to remote ollama server = ", remotesrv)
        message("remote model = ", model1)        
      }
      prms <- list(model = model1)
      url <- paste0("http://", remotesrv, "/api/chat")
      llm <- tidyprompt::llm_provider_ollama(parameters = prms, url = url)
    }
  } else if (grepl("^groq:", model)) {
    model2 <- sub("groq:", "", model)
    prms <- list(model = model2)
    llm <- tidyprompt::llm_provider_groq(parameters = prms)
  } else if (grepl("^grok|^xai:", model)) {
    model2 <- sub("^xai:", "", model)
    prms <- list(model = model2)
    llm <- tidyprompt::llm_provider_xai(parameters = prms)
  } else if (grepl("^gpt-|^openai:", model)) {
    model2 <- sub("^openai:", "", model)
    prms <- list(model = model2)
    llm <- tidyprompt::llm_provider_openai(parameters = prms)
  } else if (grepl("^gemini-|^google:", model)) {
    model2 <- sub("^google:", "", model)
    prms <- list(model = model2)
    key <- Sys.getenv("GEMINI_API_KEY")
    llm <- tidyprompt::llm_provider_google_gemini(parameters = prms, api_key = key)
  }

  if (is.null(llm)) {
    message("warning. unsupported model: ", model)
    return(NULL)
  }
  
  if (verbose > 0) {
    message("model = ", model)
    message("question = ", question)
  }

  resp <- NULL

  resp <- question |>
    tidyprompt::send_prompt(
      llm_provider = llm,
      clean_chat_history = TRUE,
      verbose = FALSE,
      return_mode = "only_response"
    )

  resp <- sub("<think>.*</think>", "", resp)

  return(resp)

}

#' Generate image with Gemini (aka Nano Banana).
#' Note this model handles very large prompts correctly.
#' @export
ai.create_image_gemini <- function(prompt,
                                   model = "gemini-2.5-flash-image",
                                   api_key = Sys.getenv("GEMINI_API_KEY"),
                                   format = c("file","base64","raw")[1], 
                                   filename = "image.png",
                                   aspectRatio = "16:9",
                                   imageSize = "1K",
                                   base_url = "https://generativelanguage.googleapis.com/v1beta") {

  assertthat::assert_that(assertthat::is.string(prompt), assertthat::noNA(prompt))
  assertthat::assert_that(assertthat::is.string(model), assertthat::noNA(model))
  assertthat::assert_that(assertthat::is.string(api_key), assertthat::noNA(api_key))
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required for ai.create_image_gemini()")
  }
  `%>%` <- dplyr::`%>%`

  if (nchar(api_key) == 0) {
    stop("GEMINI_API_KEY environment variable is not set", call. = FALSE)
  }

  message("calling gemini image generation...")
  model <- sub("^google:", "", model)
  url <- glue::glue("{base_url}/models/{model}:generateContent")

  headers <- c(`x-goog-api-key` = api_key, `Content-Type` = "application/json")

  body <- list(
    contents = list(list(parts = list(list(text = prompt)))),
    generationConfig = list(
      responseModalities = list("TEXT", "IMAGE"),
      imageConfig =  list(aspectRatio = aspectRatio, imageSize = imageSize)      
    )
  )

  if (grepl("gemini-2.5",model)) {
    body$generationConfig$imageConfig <- list(aspectRatio = aspectRatio)
  }
  
  response <- httr::POST(
    url = url,
    httr::add_headers(.headers = headers),
    body = jsonlite::toJSON(body, auto_unbox = TRUE),
    encode = "raw"
  )

  httr::http_type(response)
  if (httr::http_type(response) != "application/json") {
    stop("Gemini API returned unexpected content type", call. = FALSE)
  }

  parsed <- response %>%
    httr::content(as = "text", encoding = "UTF-8") %>%
    jsonlite::fromJSON(flatten = TRUE)

  httr::http_error(response)
  if (httr::http_error(response)) {
    error_msg <- if (!is.null(parsed$error$message)) parsed$error$message else "Unknown error"
    stop(paste0("Gemini API request failed [", httr::status_code(response), "]: ", error_msg), call. = FALSE)
  }
  
  parts <- parsed$candidates$content.parts  
  b64 <- NULL
  mimetype <- NULL
  for (part in parts) {
    if (!is.null(part$inlineData.data)) {
      b64 <- part$inlineData.data
      b64 <- head(b64[!is.na(b64)],1)
      mimetype <- part$inlineData.mimeType
      mimetype <- head(mimetype[!is.na(mimetype)],1)
      break()
    }
  }
  
  if (is.null(b64) || length(b64)==0) stop("No image data found in response")

  if (format == "file") {
    raw_image <- base64enc::base64decode(b64)    
    filetype <- sub("jpeg", "jpg", sub("image/","",mimetype))
    filename2 <- paste0(sub("[.].*$", "", filename), ".", filetype)
    writeBin(raw_image, filename2)
    message("Saved image to: ", filename2)
    return(invisible(filename2))
  }

  if (format == "raw") {
    raw_image <- base64enc::base64decode(b64)    
    return(invisible(raw_image))
  }

  if (format == "base64") return(invisible(b64))

  stop("return error")

}

#' Generate image with OpenAI's dallE.
#' Note: limitation of the prompt is about 1000 characters. 
#' @export
ai.create_image_openai <- function (prompt,
                                    model = NULL, 
                                    size = c("512x512","1024x1024","256x256")[1], 
                                    format = c("file","base64","raw"),
                                    filename = "image.png",
                                    api_key = Sys.getenv("OPENAI_API_KEY"),
                                    base_url = "https://api.openai.com/v1",
                                    user = NULL,
                                    organization = NULL)  {
  
  format <- match.arg(format)
  assertthat::assert_that(assertthat::is.string(prompt), assertthat::noNA(prompt))
  assertthat::assert_that(assertthat::is.string(size), assertthat::noNA(size))
  assertthat::assert_that(assertthat::is.string(format), assertthat::noNA(format))

  if (!is.null(user)) {
    assertthat::assert_that(assertthat::is.string(user), assertthat::noNA(user))
  }

  assertthat::assert_that(assertthat::is.string(api_key), assertthat::noNA(api_key))

  if (!is.null(organization)) {
    assertthat::assert_that(assertthat::is.string(organization), assertthat::noNA(organization))
  }

  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required for ai.create_image_openai()")
  }
  `%>%` <- dplyr::`%>%`

  model <- sub("^openai:", "", model)
  
  if (grepl("api.x.ai", base_url, fixed = TRUE)) {
    message("calling grok ($0.07 per image)")    
    if (is.null(model)) model <- "grok-2-image-1212"
    size <- NULL
  } else if (grepl("api.openai.com", base_url, fixed = TRUE)) {
    message("calling openai ($0.05 per image)")    
  } else {
    stop("invalid base_url =", base_url)
  }
  
  url <- glue::glue("{base_url}/images/generations")
  headers <- c(Authorization = paste("Bearer", api_key), `Content-Type` = "application/json")

  body <- list()
  body[["model"]] <- model
  body[["prompt"]] <- prompt
  body[["n"]] <- 1
  body[["response_format"]] <- "b64_json"
  body[["size"]] <- size
  body[["user"]] <- user
  response <- httr::POST(url = url, httr::add_headers(.headers = headers), body = body, encode = "json")

  httr::http_type(response)
  if (httr::http_type(response) != "application/json") {
    paste("OpenAI API probably has been changed. Please check online documentation.") %>% stop()
  }

  parsed <- response %>% httr::content(as = "text", encoding = "UTF-8") %>% jsonlite::fromJSON(flatten = TRUE)
  if (httr::http_error(response)) {
    error_msg <- parsed$error
    if(is.list(error_msg)) error_msg <- parsed$error$message
    paste0("OpenAI API request failed [", httr::status_code(response), 
      "]:\n\n", error_msg) %>% stop(call. = FALSE)
  }

  b64 <- parsed$data[['b64_json']]
  if (is.null(b64) || length(b64) == 0) stop("No image data found in parsed response")

  if (format == "file") {
    raw_image <- base64enc::base64decode(b64)    
    writeBin(raw_image, filename)
    message("Saved image to: ", filename)
    return(invisible(filename))
  }

  if (format == "raw") {
    raw_image <- base64enc::base64decode(b64)    
    return(invisible(raw_image))
  }

  if (format == "base64") return(invisible(b64))

  stop("return error")
  
}

#' Generate image with Grok (which uses Flux).
#' Note: limitation of the prompt is about 1000 characters.
#' @export
ai.create_image_grok <- function(prompt,
                                 model = "grok-2-image-1212", 
                                 format = c("file","base64","raw")[1], 
                                 api_key = Sys.getenv("XAI_API_KEY"),
                                 base_url = "https://api.x.ai/v1",
                                 filename = "image.png",
                                 user = NULL,
                                 organization = NULL) {

  model <- sub("^grok:", "", model)
  ai.create_image_openai(
    prompt,
    size = "default",
    format = format,
    filename = filename,
    model = model,
    base_url = base_url,
    api_key = api_key,
    user = user,
    organization = organization
  )

}
