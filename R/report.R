## ======================================================================
## Report and diagram generation functions for WGCNAplus
## Extracted from playbase/R/pgx-R
## ======================================================================

#' Create report
#' @param wgcna A WGCNA result object.
#' @param ai_model LLM model name string.
#' @param graph Optional graph object.
#' @param annot Annotation table or NULL.
#' @param multi Use multi-omics mode.
#' @param ntop Number of top entries.
#' @param topratio Ratio threshold for top modules.
#' @param psig P-value significance threshold.
#' @param format Output format: "markdown" or "html".
#' @param verbose Verbosity level.
#' @param progress Optional Shiny progress object.
#' @return List with report, diagram, and descriptions.
#' @export
create_report <- function(wgcna,
                          ai_model,
                          graph = NULL,
                          annot = NULL,
                          multi = FALSE,
                          ntop = 100,
                          topratio = 0.85,
                          psig = 0.05,
                          format = "markdown",
                          verbose = 1,
                          progress = NULL) {

  if (is.null(ai_model)) ai_model <- ""
  if (is.null(graph) && !is.null(wgcna$graph)) graph <- wgcna$graph
  
  if (!multi) {
    layers <- list(gx = wgcna)
  } else if (!is.null(wgcna$layers)) {
    layers <- wgcna$layers
  } else {
    layers <- wgcna
  }

  ## get top modules (most correlated with some phenotype)
  top.modules <- getTopModules(layers, topratio = topratio, kx = 4, multi = TRUE)

  if (is.null(annot) && !is.null(layers[[1]]$annot)) annot <- layers[[1]]$annot
  
  if (is.null(annot)) message("[create_report] No annot supplied. This is recommended.")

  ## Step 1. Describe modules with LLM. We can use one LLM model or more.
  if (verbose) message("Extracting top modules...")
  out <- describeModules(
    layers,
    modules = top.modules,
    ntop = ntop,
    annot = annot,
    psig = psig,
    experiment = wgcna$experiment,
    verbose = verbose,
    model = ""
  )
  descriptions_prompts <- out$questions
  descriptions <- out$answers

  ## Step 2: Make consensus summary from the descriptions.
  summaries <- list()
  summaries_prompts <- list()
  results <- NULL
  if (ai_model != "") {
    if (verbose) message("Simmering modules...")
    for (k in names(descriptions)) {
      ss <- descriptions[[k]]
      q2 <-  paste("Following are descriptions of a certain WGCNA module by one or more LLMs. Create a consensus conclusion out of the independent descriptions. Describe the underlying biology, relate correlated phenotypes and mention key genes, proteins or metabolites. Just answer, no confirmation, use 1-2 paragraphs. Use prose as much as possible, do not use tables or bullet points.\n\n", ss)
      cc <- ai.ask(q2, model=ai_model)
      summaries[[k]] <- cc
      summaries_prompts[[k]] <- q2
    }
    results <- summaries
  } else {
    if (verbose) message("Skipping module summaries...")
    results <- descriptions
  }

  ## addd compute setttings
  if (!is.null(wgcna$settings)) {
    results[['compute_settings']] <- paste0(names(wgcna$settings),'=',wgcna$settings,collapse='; ')
  }

  ## collate all results
  all.results <- lapply(names(results), function(me)
    paste0("================= ",me," =================\n\n", results[[me]],"\n"))
  all.results <- paste(all.results, collapse="\n")

  ## Step 3: Make detailed report. We concatenate all summaries and
  ## ask a (better) LLM model to create a report.
  if (verbose) message("Baking full report...")
  qq = diagram = report = NULL
  if (ai_model == "") {
    report <- all.results
  } else {
    qq <- "These are the results of a WGCNA analysis. There are descriptions of the most relevant modules. Create a detailed report for this experiment. Give a detailed interpretation of the underlying biology by connecting WGCNA modules into biological functional programs, referring to key genes, proteins or metabolites. Build an cross-module integrative biological narrative. Suggest similarity to known diseases and possible therapies. Add a discussion and conclusion. Omit abstract, future directions, limitations, or references. Add a short paragraph describing methods and compute settings at the end. Important: format like a scientific article, use prose as much as possible, minimize the use of tables and bullet points. For long tables show at least the top 5, and at most top 10, up and down entries. Do not inject any inline code. Only write if there was evidence in the source text."

    if (multi) qq <- gsub("WGCNA","multiomics WGCNA",qq)

    xx <- wgcna$experiment
    pp <- paste("You are a biologist interpreting results from a WGCNA analysis for this experiment:",  xx, ".\n\n")
    qq <- paste(pp, qq)

    if (format == "markdown") qq <- paste(qq, "Format text and sections as markdown.")
    if (tolower(format) == "html") qq <- paste(qq, "Format text and sections as HTML.")
    qq <- paste(qq, "\n\nnHere are the results: <results>",all.results,"\n</results>")

    ## Ask LLM
    report <- ai.ask(qq, model = ai_model)
    report <- gsub("^```html|```$","",report)

    ## Step 4: Create diagram from report
    if (verbose) message("Mashing up diagram...")
    diagram <- create_diagram(report, graph = graph, ai_model = ai_model)

  }

  return(
    list(
      descriptions_prompts = descriptions_prompts,
      descriptions = descriptions,
      summaries_prompts = summaries_prompts,
      summaries = summaries,
      report_prompt = qq,
      report = report,
      diagram = diagram
    )
  )

}


#' Describe WGCNA modules using LLM
#' @param wgcna A WGCNA result object.
#' @param ntop Number of top entries.
#' @param psig P-value significance threshold.
#' @param annot Annotation table or NULL.
#' @param modules Module names to describe.
#' @param experiment Experiment description string.
#' @param verbose Verbosity level.
#' @param model LLM model name or NULL.
#' @param docstyle Style for LLM description.
#' @param numpar Max paragraphs for LLM output.
#' @param level Feature level: "gene" or "geneset".
#' @return List with prompt, questions, and answers.
#' @export
describeModules <- function(wgcna,
                            ntop = 50,
                            psig = 0.05,
                            annot = NULL,
                            modules = NULL,
                            experiment = "",
                            verbose = 1, 
                            model = getOption("WGCNAplus.default_llm"),
                            docstyle = "detailed summary",
                            numpar = 2,
                            level = "gene")  {
  
  if (is.null(annot)) {
    message("[describeModules] WARNING. user annot table is recommended.")
  }

  top <- getTopTables(wgcna, annot=annot, ntop=ntop,
    psig=psig, level=level, rename="gene_title")

  if (is.null(modules)) {
    modules <- union(names(top$genes), names(top$sets))
  }

  if (is.null(experiment)) experiment <- ""

  if (length(modules)==0) {
    info("[describeModules] warning: empty module list!")
    return(NULL)
  }

  ## If no LLM is available we do just a manual summary
  model <- setdiff(model, c("", NA))
  if (is.null(model) || length(model) == 0) {
    desc <- list()
    for (m in modules) {

      ss=gg=pp="<none>"

      if (!is.null(top$genes[[m]])) {
        gg <- paste( top$genes[[m]], collapse=', ')
      }

      if (!is.null(top$sets[[m]])) {
        ss <- paste( sub(".*:","",top$sets[[m]]), collapse='; ')
      }

      if (m %in% names(top$pheno)) {
        pp <- paste( top$pheno[[m]], collapse='; ')
      }

      d <- ""
      if (!is.null(pp)) d <- paste(d, "**Correlated phenotypes**:", pp, "\n\n")

      if (!is.null(gg) && gg != "") {
        d <- paste(d, "**Key genes**:", gg, "\n\n")
      }

      if (!is.null(ss) && ss != "") {
        d <- paste(d, "**Top enriched gene sets**:", ss, "\n\n")
      }

      desc[[m]] <- d

    }

    res <- list(prompt = NULL, questions = NULL, answers = desc)

    return(res)

  }

  prompt <- paste("Give a", docstyle, "of the main overall biological function of the following top enriched genesets belonging to module <MODULE>. Discuss the possible relationship with phenotypes <PHENOTYPES> of this experiment about \"<EXPERIMENT>\". Use maximum", numpar, "paragraphs. Do not use any bullet points. \n\nHere is list of enriched gene sets: <GENESETS>\n")

  if (verbose > 1) cat(prompt)

  desc <- list()
  questions <- list()
  for (k in modules) {
    if (verbose > 0) message("Describing module ", k)

    ss=gg=pp=""
    if (length(top$sets[[k]])>0) {
      ss <- sub( ".*:","", top$sets[[k]] ) ## strip prefix
      ss <- paste(ss, collapse=';')
    } else {
      ss <- "[no significant genesets]"
    }

    if (k %in% names(top$pheno)) {
      pp <- paste0("'",top$pheno[[k]],"'")
      pp <- paste( pp, collapse=';')
    }

    q <- prompt

    if (length(top$genes[[k]])>0) {
      gg <- paste( top$genes[[k]], collapse=';')
      q <- paste(q, "\nHere is the list of key genes/proteins/metabolites: <KEYGENES>\n")
    }

    q <- sub("<MODULE>", k, q)
    q <- sub("<PHENOTYPES>", pp, q)
    q <- sub("<EXPERIMENT>", experiment, q)
    q <- sub("<GENESETS>", ss, q)
    q <- sub("<KEYGENES>", gg, q)

    answer <- ""
    for (m in model) {
      if (verbose > 0) message("  ...asking LLM model ", m)
      a <- ai.ask(q, model = m)
      a <- paste0(a, "\n\n[AI generated using ", m, "]\n")
      if (length(model) > 1) a <- paste0("\n-------------------------------\n\n", a)
      answer <- paste0(answer, a)
    }

    desc[[k]] <- answer
    questions[[k]] <- q
  }

  res <- list(prompt = prompt, questions = questions, answers = desc)
  
  return(res)

}

#' Correct common DOT diagram issues
#' @param diagram DOT diagram string.
#' @return Corrected DOT diagram string.
#' @keywords internal
correct_dot_diagram <- function(diagram) {

  ## force as digraph
  diagram <- sub("^graph","digraph",diagram)

  ## crazy arrows
  diagram <- gsub("-x->","->",diagram)

  ## remove anything after DOT last curly bracket
  diagram <- sub("\\}\n.*","}\n",diagram)
  diagram <- gsub("\\[solid\\]","[style=solid]",diagram)
  diagram <- gsub("\\[dashed\\]","[style=dashed]",diagram)

  # avoid these problematic colors
  diagram <- gsub("lightgreen","palegreen", diagram)  ## avoid
  diagram <- gsub("lightorange","lightsalmon", diagram)  ## avoid
  diagram <- gsub("fillcolor=black","fillcolor=lightgrey", diagram)  ## avoid
  diagram <- gsub("#000000","#AAAAAA",diagram)  ## no black

  # match 3- or 6-digit hex color with replace with quoted version
  diagram <- gsub("(?<!['\"])\\b(#(?:[0-9A-Fa-f]{3}){1,2})\\b(?!['\"])",
    "\"\\1\"", diagram, perl = TRUE )

  return(diagram)

}


#' Create DOT string object from graph object for sending to
#' LLM. Clean unneeded attributes.
#' @param graph An igraph object.
#' @return DOT format string.
#' @keywords internal
#' @export
graph2dot <- function(graph) {

  aa <- names(igraph::vertex_attr(graph))
  aa <- intersect(aa, c("layer","fc","value"))
  for (a in aa)  graph <- igraph::delete_vertex_attr(graph, a)

  bb <- names(igraph::edge_attr(graph))
  bb <- intersect(bb, c("rho","connection_type"))
  for (b in bb)  graph <- igraph::delete_edge_attr(graph, b)

  file <- tempfile(fileext = ".dot")
  igraph::write_graph(graph, file=file, format="dot")
  dot <- readChar(file, file.info(file)$size)
  unlink(file)

  return(dot)

}

#' Change layout of DOT string object
#' @param dot DOT format string.
#' @param dir Direction: "TB" or "LR".
#' @return Modified DOT string.
#' @keywords internal
#' @export
dot.rankdir <- function(dot, dir) {

  if (dir == "TB") dot <- sub("rankdir=LR","rankdir=TB", dot)
  if (dir == "LR") dot <- sub("rankdir=TB","rankdir=LR", dot)

  return(dot)

}

#' Create module diagram in DOT format
#' @param wgcna_report Report text string.
#' @param ai_model LLM model name string.
#' @param graph Optional graph template.
#' @param rankdir DOT layout direction.
#' @param correct Apply DOT corrections.
#' @param double.check Validate DOT with DiagrammeR.
#' @return DOT diagram string.
#' @export
create_diagram <- function(wgcna_report,
                           ai_model,
                           graph = NULL,
                           rankdir = "TB",
                           correct = TRUE,
                           double.check = TRUE) {

  if (!is.null(graph)) {
    ## If we pass a graph (from e.g. Lasagna) we tell the LLM to use
    ## the graph as starting point or template. This constrains the
    ## connections and minimizes 'hallucinations'
    message("[create_diagram] using graph template...")
    dot <- graph2dot(graph)
    qq <- paste0("Create a directed diagram connecting modules according to the following WGCNA report. Use the given undirected graph as starting point and use known scientific information to infer connectivity and directionality. All modules must be connected with at least one other module. Annotate modules with main biological function and key features (gene, proteins or metabolites). Add extra nodes for inferred intermediate phenotypes. Determine which phenotypes are causal and which phenotypes are observed effects. Suggest cause and effect relations that explain phenotypes and modules. Group modules with same biological functions. Give just the code in clean DOT format.

Layout in TB direction. Do not use any special characters, without headers or footer text. Do not use subgraphs. Do not use hexadecimal color coding. Use solid lines for positive regulation, use dashed lines for negative regulation. Annotate modules with module name, biological function and key gene/protein or metabolite. Use rectangular shapes for module nodes, use oval shapes for phenotype nodes. Color fill nodes matching the module names with light palette so we can still read well the text. Never use black for fill. Again, do not fill any nodes with black, use grey instead. Color phenotype nodes lightyellow. ")
    qq <- paste(qq,
      "\n\n<report>", wgcna_report, "</report>",
      "\n\n<dot>", dot, "</dot>" )

  } else {
    ## If we do not pass a graph (from e.g. Lasagna) we let the LLM
    ## connect the modules itself based on its external knowledge.
    message("[create_diagram] no graph template...")
    qq <- paste0("Create a diagram connecting modules in the following WGCNA report. Annotate modules with main biological function and key features (gene, proteins or metabolites). Add phenotype nodes. Suggest cause and effect relations that explain the phenotypes. Group modules with same biological functions. Give just the code in clean DOT format. Layout in TB direction. Do not use any special characters, without headers or footer text. Do not use subgraphs. Do not use hexadecimal color coding. Use solid lines for positive regulation, use dashed lines for negative regulation. Color fill nodes matching the module names with light palette so we can still read well the text. Never use black for fill. Again, do not fill any nodes with black, use grey instead. Color phenotype nodes lightyellow.")
    qq <- paste(qq, "\n\n<report>", wgcna_report, "</report>")
  }

  aa <- ai.ask(qq, model = ai_model)

  ## cleanup
  aa <- gsub("```","",aa)
  diagram <- gsub("mermaid\n|dot\n","",aa)
  diagram <- gsub("&","and",diagram)
  diagram <- dot.rankdir(diagram, dir=rankdir)
  if (correct) diagram <- correct_dot_diagram(diagram)

  if (double.check) {
    if (!requireNamespace("DiagrammeR", quietly = TRUE)) stop("Package 'DiagrammeR' is required")
    if (!requireNamespace("DiagrammeRsvg", quietly = TRUE)) stop("Package 'DiagrammeRsvg' is required")
    code.error <- TRUE
    ntry <- 1
    while (code.error && ntry <= 5) {
      ## check valid code
      dg <- DiagrammeR::grViz(diagram)
      out <- try(DiagrammeRsvg::export_svg(dg))
      code.error <- inherits(out, "try-error")
      if (code.error) {
        ## try to correct
        diagram <- ai.ask(paste("Please double check the following DOT diagram code and correct. If the code is correct, do not change anything. Just return the corrected clean code:",diagram), model = ai_model)
        ntry <- ntry + 1
      }
    }
  }
  return(diagram)
}


#' Given the output from create_report() this function creates
#' an infographic by calling Gemini3. The genAI model is given the
#' report and asked to adhere to the included or external given
#' diagram (in DOT format).
#' @param report Report text string.
#' @param diagram Optional DOT diagram string.
#' @param prompt Optional custom prompt.
#' @param model Gemini model name.
#' @param filename Output PNG filename.
#' @return Path to output image file.
#' @export
create_infographic <- function(report,
                               diagram = NULL,
                               prompt = NULL,
                               model = "gemini-3-pro-image-preview",
                               filename = "infographic.png",
                               api_key = Sys.getenv("GEMINI_API_KEY")) {

  prompt <- paste(prompt, "\nCreate a graphical abstract according to the given diagram and information in the WGCNA report. Use scientific visual style like Nature journals. Illustrate biological concepts with small graphics. \n\n", report, "\n---------------\n\n", diagram)

  outfile <- try(
    ai.create_image_gemini(
      prompt = prompt,
      model = model,
      format = "file",
      filename = filename,
      api_key = api_key
    )
  )

  if (inherits(outfile,"try-error")) return(NULL)

  return(invisible(outfile))

}

#' Create infographic for a single module
#' @param rpt Report object from create_report.
#' @param module Module name string.
#' @param prompt Optional custom prompt.
#' @param model Gemini model name.
#' @param filename Output PNG filename.
#' @return Path to output image file.
#' @export
create_module_infographic <- function(rpt,
                                      module,
                                      prompt = NULL,
                                      model = "gemini-3-pro-image-preview",
                                      filename = "module-infographic.png",
                                      api_key = Sys.getenv("GEMINI_API_KEY")) {

  if (!module %in% names(rpt$summaries)) {
    stop(paste("module", module, "not in report summaries"))
  }

  mm <- paste0("**", module, "**: ", rpt$summaries[[module]])

  prompt <- paste(prompt, "Create an infographic summarizing the biological narrative of the following WGCNA module. Use scientific visual style like Nature journals. Illustrate biological concepts with small graphics. Match the background with the name of the module with a very light shade. Include the module name in the title or image. \n\n", mm)

  outfile <- ai.create_image_gemini(
    prompt,
    model,
    filename = filename,
    api_key = api_key
  )
  message("saving to ", outfile)

  return(invisible(outfile))

}
