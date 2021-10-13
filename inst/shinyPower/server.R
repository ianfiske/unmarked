if (exists(".SHINY_MODEL")) {
  mod <- .SHINY_MODEL
} else {
  object <- get(".SHINY_MODEL", envir = unmarked:::.shiny_env)
}

coefs <- unmarked:::check_coefs(NULL, mod, TRUE)

inline_wrap <- function(f, ...){
  out <- f(...)
  div(style='display:inline-block; width: 100px; vertical-align:top', out)
}
options(unmarked_shiny=TRUE)
#on.exit(options(unmarked_shiny=FALSE))

get_coef_ui <- function(coefs, nulls=FALSE){

  out <- c(list(h3("Coefficient values")))
  parbase <- "coef_"
  if(nulls){
    out <- c(list(h3("Null hypotheses")))
    parbase <- "null_"
  }

  for (i in 1:length(coefs)){
    pars <- coefs[[i]]
    submod_name <- names(coefs)[i]
    inps <- lapply(1:length(pars), function(x){
      par_name <- names(pars)[x]
      inp_name <- paste0(parbase,submod_name,"_",par_name)
      inline_wrap(numericInput, inputId=inp_name, label=par_name,
                  value=0, step=0.01)
    })
    out <- c(out, list(h4(submod_name)), inps)
  }
  out
}

get_coefs <- function(input, nulls=FALSE){
  parbase <- "coef_"
  if(nulls) parbase <- "null_"
  pass <- reactiveValuesToList(input)
  pass$shinymanager_where <- NULL
  inp_sub <- pass[grepl(parbase,names(pass), fixed=TRUE)]
  inp_sub <- pass[!is.na(names(inp_sub))]
  names(inp_sub) <- gsub(parbase, "", names(inp_sub))
  submods <- gsub("_(.*)$","",names(inp_sub))
  pars <- gsub("^(.*)_","",names(inp_sub))
  out <- lapply(unique(submods), function(x){
    vals <- unlist(inp_sub[which(submods==x)])
    names(vals) <- pars[which(submods==x)]
    vals
  })
  names(out) <- unique(submods)
  out
}

get_design_ui <- function(input, default, name){
  #out <- c(list(h4(paste0("Number of ",name))))
  #inp <- reactiveValuesToList(input)
  nval <- input[[paste0("ndesign_",name)]]
  inps <- lapply(1:nval, function(x){
      inp_name <- paste0("design_",name,"_",x)
      inline_wrap(numericInput, inputId=inp_name, label=NULL,
                  value=100, min=1, step=1)
    })
  inps
  #c(out, inps)
}

get_design <- function(input){
  pass <- reactiveValuesToList(input)
  pass$shinymanager_where <- NULL
  inp_sub <- unlist(pass[grepl("design_sites_",names(pass),fixed=TRUE)])
  if(length(inp_sub)==1 & inp_sub[1] == numSites(mod@data)) return(NULL)
  #list(M=inp_sub, J=obsNum(mod@data), T=1)
  expand.grid(M=inp_sub, J=obsNum(mod@data), T=1)
}

run_analysis <- function(mod, coefs, alpha, nsim, nulls, design){
  if(!is.null(design) && any(sapply(design, length) > 1)){
    out <- unmarkedPowerList(mod, coefs, design, alpha, nsim)
  } else {
    out <- powerAnalysis(mod, coefs=coefs, alpha=alpha, nsim=nsim, nulls=nulls,
                         design=design)
  }
  out
}

function(input, output, session){

  #res_auth <- secure_server(
  #  check_credentials = check_credentials(credentials)
  #)

  #output$auth_output <- renderPrint({
  #  reactiveValuesToList(res_auth)
  #})

  options(unmarked_shiny_session=session)
  output$plot <- renderPlot(plot(mod))
  output$coef_ui <- renderUI(get_coef_ui(coefs))
  output$coefs <- renderPrint(get_coefs(input))
  output$null_ui <- renderUI(get_coef_ui(coefs, nulls=TRUE))
  output$nulls <- renderPrint(get_coefs(input, nulls=TRUE))
  output$mod <- renderUI(HTML(paste0("<b>Model:</b> ","mod")))
  output$class <- renderUI(HTML(paste0("<b>Type:</b>&nbsp&nbsp&nbsp",
                                       class(mod)[1])))
  output$sites <- renderUI(HTML(paste0("<b>Sites:</b>&nbsp&nbsp",
                                       numSites(mod@data))))
  output$design_sites <- renderUI(get_design_ui(input,numSites(mod@data),"sites"))

  observeEvent(input$run, {
    coefs <- isolate(get_coefs(input))
    nulls <- isolate(get_coefs(input, nulls=TRUE))
    design <- isolate(get_design(input))
    alpha <- isolate(input$alpha)
    nsims <- isolate(input$nsims)
    output$summary <- renderTable(
                      summary(run_analysis(mod, coefs, alpha, nsims, nulls, design))
                      #summary(powerAnalysis(mod, coefs=coefs, alpha=alpha,
                      #                      nsim=nsims, nulls=nulls, design=design))
                      )
  })
}



