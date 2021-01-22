

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(
  withMathJax(),
  # Application title
  titlePanel(HTML("<b> <font size='6'> Evaluate the trustworthiness of your study!</font> </b>")),
  
  fluidRow(
    column(8, h4("This application is based on the research program strategy (RPS) as explained in Krefeld-Schwalb, Witte & Zenker (2018). Please follow instructions to your left (under Input 1 and 2), then read off results to your right.")
           ,HTML("<b> <font size='5'> Input 1: Simulate likelihood ratios prior to data collection:</font> </b>"))),
  
  fluidRow(
    column(3,
           HTML("<font size='3'>Instructions: Indicate the desired power (1-&beta;) and level of significance (&alpha;) along with the hypothesized effect size. Read off results from Figs. 1 and 2, and Table 1 to your right.</font>"),
           numericInput("alpha",
                        "$$ \\alpha $$",
                        min = .001,
                        max = .999,
                        step = .01,
                        value =.05),
           
           numericInput("d",
                        "effect size (d)",
                        min = 0,
                        max = 1,
                        step = .01,
                        value =.2),
           
           numericInput("pow",
                        "Power $$ (1-\\beta) $$",
                        min = .001,
                        max = .999,
                        step = .01,
                        value =.95),
           
           h4("Research program strategy (RPS):"),
           h5("1. Preliminary Discovery"),
           h6(" $$p \\leq \\alpha, \\alpha \\leq .05,  unknown \\beta$$"),
           h5("2. Substantial Discovery"), 
           h6(" $$p \\leq \\alpha, \\alpha \\leq .05,  known \\beta$$"),
           h5("3. Preliminary Falsification"),
           h6("$$\\frac{L(d>0|x)}{L(d=0|x)} > \\frac{1-\\beta}\\alpha$$"),
           h5("4. Substantial Falsification"),
           h6("$$\\frac{L(d>\\delta|x)}{L(d=0|x)} > \\frac{1-\\beta}\\alpha$$"),
           h5("5. Preliminary Verification"),
           h6("$$\\frac{L(d=\\delta|x)}{L(d=0|x)} > \\frac{1-\\beta}\\alpha$$"),
           h5("6. Substantial Verification"),
           h6("$$\\frac{L(d=\\delta|x)}{L(d=0|x)} > \\frac{1-\\beta}\\alpha \\cap  \\frac{L(d|x)}{L(d=\\delta|x)} > \\frac{pdf(P50|d)}{pdf(P95|d)} > 4$$"),
           HTML("<font size='3'>Note. &alpha; = probability of type I error, &beta; = probability of type II error, N = estimated group size to achieve the desired power; ncp = non-centrality-parameter of the t-distribution representing H1.</font>")
    )
    , 
    
    column(9,
           
           fluidRow(
             column(3,plotOutput("Tplot", height = "400px", width= "100%"),
                    h5("Figure 1. Distribution of t-values for given H0 (white) and H1 (grey).")),
             column(4,plotOutput("plotPos", height = "400px", width= "100%"),
                    h5("Figure 2. Proportion of LR as a function of criteria.")),
             column(3,
                    tableOutput("tablePosL")
                    ,h5("Table 1. Proportion of correct positive results."))),
           
           fluidRow(
             column(3,plotOutput("Tplot0",  width= "100%"),
                    h5("Figure 3. Distribution of t-values for given H0 (white) and H1, adapted to d
                        observed (grey)")),
             column(4,plotOutput("plotneg", width= "100%")
                    ,h5("Figure 4. Proportion of LR as a function of selected criteria.")),
             column(3,
                    tableOutput("tableNegL")
                    ,h5("Table 2. Proportion of false positive results."))))
  ),
  tags$hr(style="border-color: grey;"),
  fluidRow(
    column(9, HTML("<b> <font size='5'> Input 2: Simulate likelihood ratios post-hoc </font> </b>"))),
  column(3,
         
         "Indicate the sample size N per group along with the observed t-value to post-hoc estimate the power and simulate likelihood ratios. Read off results from Figs. 6 and Table 4 to your right",
         numericInput("alphaemp",
                      "$$ \\alpha $$",
                      min = .001,
                      max = .999,
                      step = .01,
                      value =.05),
         numericInput("Temp",
                      "t-value (empirical)",
                      min = 0,
                      max = 100,
                      step = .1,
                      value =3),
         numericInput("N",
                      "N per group (empirical)",
                      min = 1,
                      max = 10000,
                      step = 1,
                      value =30)
  ),
  fluidRow(
    column(8,
           column(4, plotOutput("TplotEmp", width= "100%"),
                  h5("Figure 6. Distributions of t-values for given H1(=t empirical) und H0")
           ),
           column(4, 
                  tableOutput("tablePosLEmp")
                  ,h5("Table 3. Proportion of correct positive results.")
                  ,tableOutput("tableNegLEmp")
                  ,h5("Table 4. Proportion of false positive results.")))),
  h6("2017 Antonia Krefeld-Schwalb, antonia.krefeld-schwalb@unige.ch"))


# Define server logic required to draw a histogram
server <- function(input, output) {
  
  Nsampl=100
  L10 =Lplus00 = rep(NA, Nsampl)
  p=pgive1=L0x=Lplus0x=L0x=L1x=LEmpx= ddeviate= Tval = rep(NA, Nsampl)
  p = rep(NA, Nsampl)
  
  ### Samples from H1, given alpha, d and power
  
  samples <- reactive({
    
    #Estimate Sample Size
    Nest = ((qnorm(1- input$alpha)+qnorm(input$pow))/(input$d/(sqrt(2))))^2
    df = 2*Nest -2
    #Non-centrality parameter
    ncp = input$d*sqrt(Nest)/sqrt(2)
    #Draw Nsampl studies
    for (ss in 1:Nsampl){
      #Distribution of Treatment condition  
      x1 =rnorm(Nest,input$d, 1)
      meanx1 = mean(x1)
      
      #Distribution of Control condition  
      x0 =rnorm(Nest,0, 1)
      meanx0 = mean(x0)
      
      #Calculate t-values
      Tval[ss] = (meanx1-meanx0)/sqrt(2/Nest)
      
      #Calculate PCDF values D|L(d=0|x)
      p[ss] = 1-pt(Tval[ss], df)
      
      L0x[ss] = dt(Tval[ss], df, 0)
      Lplus0x[ss] = 1-pt(Tval[ss], df, ncp)
      
      L1x[ss] = dt(Tval[ss], df, ncp)
      LEmpx[ss] = dt(Tval[ss], df, Tval[ss])
      
      ddeviate[ss] = LEmpx[ss]/L1x[ss]
    }
    
    samp = list(Tval = Tval, p=p, L1x=L1x, L0x=L0x, Lplus0x=Lplus0x, LEmpx=LEmpx, ddeviate=ddeviate, Nest = Nest, ncp =ncp)
    samp
  })
  
  
  samples30 <- reactive({
    samp = samples()
    #Estimate Sample Size
    Nest = samp$Nest *.5
    df = 2*Nest -2
    #Non-centrality parameter
    ncp = input$d*sqrt(Nest)/sqrt(2)
    #Draw Nsampl studies
    for (ss in 1:Nsampl){
      #Distribution of Treatment condition  
      x1 =rnorm(Nest,input$d, 1)
      meanx1 = mean(x1)
      
      #Distribution of Control condition  
      x0 =rnorm(Nest,0, 1)
      meanx0 = mean(x0)
      
      #Calculate t-values
      Tval[ss] = (meanx1-meanx0)/sqrt(2/Nest)
      
      #Calculate PCDF values D|L(d=0|x)
      p[ss] = 1-pt(Tval[ss], df)
      
      L0x[ss] = dt(Tval[ss], df, 0)
      Lplus0x[ss] = 1-pt(Tval[ss], df, ncp)
      
      L1x[ss] = dt(Tval[ss], df, ncp)
      LEmpx[ss] = dt(Tval[ss], df, Tval[ss])
      
      ddeviate[ss] = LEmpx[ss]/L1x[ss]
    }
    
    samp30 = list(L1x=L1x, L0x=L0x, Lplus0x=Lplus0x)
    samp30
  })
  
  
  
  ### Samples from differetn H relative to H1, given alpha, and Nest d
  sumLR <- reactive({
    samp <- samples()
    
    del = c(0,.25*input$d, .5*input$d, .75*input$d,
            min(1,1.25*input$d),min(1,1.5*input$d), min(1,1.75*input$d), min(1,2*input$d),1)
    probLR <- matrix(NA, nrow= length(del), ncol = 2)
    #Estimate Sample Size
    Nest = ((qnorm(1- input$alpha)+qnorm(input$pow))/(input$d/(sqrt(2))))^2
    df = 2*Nest -2
    #Non-centrality parameter
    for (dd in 1:length(del))
    {
      ncp = del[dd]*sqrt(Nest)/sqrt(2)
      #Draw Nsampl studies
      for (ss in 1:Nsampl){
        #Distribution of Treatment condition  
        x1 =rnorm(Nest,del[dd], 1)
        meanx1 = mean(x1)
        
        #Distribution of Control condition  
        x0 =rnorm(Nest,0, 1)
        meanx0 = mean(x0)
        
        #Calculate t-values
        Tval[ss] = (meanx1-meanx0)/sqrt(2/Nest)
        
        #Calculate PCDF values D|L(d=0|x)
        p[ss] = 1-pt(Tval[ss], df)
        
        L0x[ss] = dt(Tval[ss], df, 0)
        Lplus0x[ss] = 1-pt(Tval[ss], df, ncp)
        
        L1x[ss] = dt(Tval[ss], df, ncp)
        LEmpx[ss] = dt(Tval[ss], df, Tval[ss])
        
        ddeviate[ss] = LEmpx[ss]/L1x[ss]
      }
      probLR[dd,1:2] = c(mean(Lplus0x/L0x >=input$pow/input$alpha),
                         mean(L1x/L0x>=input$pow/input$alpha))
    }
    probLR <- rbind(probLR[1:4,], 
                    c(mean(samp$Lplus0x/samp$L0x >=input$pow/input$alpha), mean(samp$L1x/samp$L0x>=input$pow/input$alpha)),
                    probLR[5:9,])
    
    pLR = list(probLR = probLR, d = c(del[1:4], input$d, del[5:9]))
    pLR
  })
  
  ### Samples from H0, given alpha, and power and d assumed
  samples0 <- reactive({
    Nest = ((qnorm(1- input$alpha)+qnorm(input$pow))/(input$d/(sqrt(2))))^2
    
    df =2*round(Nest) -2
    
    for (ss in 1:Nsampl){
      x0 =rnorm(round(Nest),0, 1)
      x0b =rnorm(round(Nest),0, 1)
      
      meanx0 = mean(x0)
      meanx0b = mean(x0b)
      
      demp = (mean(x0)-mean(x0b))/(sqrt((sd(x0)^2+sd(x0b)^2)/2))
      #ncp = demp*sqrt(round(Nest))/sqrt(2)
      ncp = input$d*sqrt(round(Nest))/sqrt(2)
      
      Tval[ss] = (meanx0-meanx0b)/sqrt(2/round(Nest))
      Tval0 = 0
      
      p[ss] = 1-pt(Tval[ss], df)
      L0x[ss] = dt(Tval[ss], df)
      Lplus0x[ss] = 1-pt(Tval[ss], df, ncp)
      L1x[ss] =dt(Tval[ss], df, ncp)
    }
    samp0 =list(p=p,Tval = Tval,  L0x=L0x, Lplus0x=Lplus0x, L1x=L1x,Nest=Nest, ncp = ncp)
    samp0
  })
  
  
  samplesEmp <- reactive({
    
    demp = sqrt(2/input$N)*input$Temp
    pow = pnorm(sqrt(input$N)*demp/sqrt(2)-qnorm(1-input$alphaemp))
    df = 2*input$N -2
    
    for (ss in 1:Nsampl){
      #Distribution of Treatment condition  
      x1 =rnorm(input$N,demp, 1)
      meanx1 = mean(x1)
      
      #Distribution of Control condition  
      x0 =rnorm(input$N,0, 1)
      meanx0 = mean(x0)
      
      #Calculate t-values
      Tval[ss] = (meanx1-meanx0)/sqrt(2/input$N)
      
      #Calculate PCDF values D|L(d=0|x)
      p[ss] = 1-pt(Tval[ss], df)
      
      L0x[ss] = dt(Tval[ss], df, 0)
      Lplus0x[ss] = 1-pt(Tval[ss], df, input$Temp)
      
      L1x[ss] = dt(Tval[ss], df, input$Temp)
      
    }
    
    
    samp = list(d=demp, pow = pow,Tval = Tval, p=p, L1x=L1x, L0x=L0x, Lplus0x=Lplus0x)
    samp
  })
  
  
  samplesEmp0 <- reactive({
    sampEmp = samplesEmp()
    df =2*round(input$N) -2
    for (ss in 1:Nsampl){
      x0 =rnorm(round(input$N),0, 1)
      x0b =rnorm(round(input$N),0, 1)
      
      meanx0 = mean(x0)
      meanx0b = mean(x0b)
      
      demp = (mean(x0)-mean(x0b))/(sqrt((sd(x0)^2+sd(x0b)^2)/2))
      #ncp = demp*sqrt(round(input$N))/sqrt(2)
      ncp = demp*sqrt(round(input$N))/sqrt(2)
      
      Tval[ss] = (meanx0-meanx0b)/sqrt(2/round(input$N))
      Tval0 = 0
      
      p[ss] = 1-pt(Tval[ss], df)
      L0x[ss] = dt(Tval[ss], df)
      Lplus0x[ss] = 1-pt(Tval[ss], df, ncp)
      L1x[ss] =dt(Tval[ss], df, ncp)
    }
    samp0 =list(p=p,Tval = Tval,  L0x=L0x, Lplus0x=Lplus0x, L1x=L1x, ncp = ncp)
    samp0
  })
  
  
  
  output$Tplot <- renderPlot(res = 100, {
    samp <- samples()
    hist(rt(100, samp$Nest-1), breaks = 10, xlim=c(-max(samp$Tval), max(samp$Tval)), ylim = c(0,Nsampl/2), xlab = c(""), 
         main = "")
    hist(samp$Tval, breaks = 10, xlim=c(-max(samp$Tval), max(samp$Tval)), ylim = c(0,Nsampl/2), col = "grey", add = TRUE, xlab = c(""))
    legend(-max(samp$Tval)-.5, 50, c(paste("d =",input$d), paste("N =", round(samp$Nest)),  paste("ncp =", round(samp$ncp,2))), 
           bty = "n",
           cex = 1)
    
  })
  
  output$verification <- renderText({
    
    if (input$alpha <= .05 & input$pow >=.95 ){msg = sprintf('<font color="%s">%s</font>','green',"You are in the context of justification!")
    msg
    } else {
      msg <- sprintf('<font color="%s">%s</font>','red',"You are in the context of discovery!") 
      msg
    }}
  )
  
  output$verificationemp <- renderText({
    samp <- samplesEmp()
    if (input$alphaemp <= .05 & samp$pow >=.95 ){msg = sprintf('<font color="%s">%s</font>','green',"You are in the context of justification!")
    msg
    } else {
      msg <- sprintf('<font color="%s">%s</font>','red',"You are in the context of discovery!") 
      msg
    }}
  )
  
  output$TplotEmp <- renderPlot(res = 100, {
    samp <- samplesEmp()
    
    hist(rt(100, input$N-1), breaks = 10, xlim=c(-max(samp$Tval), max(samp$Tval)), ylim = c(0,Nsampl/2), xlab = c(""), 
         main = "")
    hist(samp$Tval, breaks = 10, xlim=c(-max(samp$Tval), max(samp$Tval)), ylim = c(0,Nsampl/2), col = "grey", add = TRUE, xlab = c(""))
    legend(-max(samp$Tval)-.5, 50, c(paste("1-beta =",round(samp$pow,2)), paste("d =",round(samp$d,2)), paste("ncp =", round(input$Temp,2))), 
           bty = "n",
           cex = 1)
    
  })
  
  output$Tplot0 <- renderPlot(res = 100, {
    samp <- samples0()
    sampT <- samples()
    hist(rt(100, samp$Nest-1), breaks = 10, xlim=c(-max(sampT$Tval), max(sampT$Tval)), ylim = c(0,Nsampl/2), xlab = c(""), 
         main = "")
    hist(samp$Tval, breaks = 10, xlim=c(-max(sampT$Tval), max(sampT$Tval)), ylim = c(0,Nsampl/2), col = "grey", add = TRUE, xlab = c(""))
    legend(-max(sampT$Tval)-.5, 50, c(paste("True d =", 0), paste("N =", round(samp$Nest)),  paste("ncp(emp) =", round(samp$ncp,2))), 
           bty = "n",
           cex = 1)
    
  })
  
  
  output$pLRplot <- renderPlot(res = 100, {
    pLR <- sumLR()
    barplot(pLR$probLR[,1], xlab = "True d", ylab= paste('P(LR) >=', round(input$pow/input$alpha,2)),
            names.arg = pLR$d, las = 1)
  })
  
  
  output$plotPos <- renderPlot(res = 100, {  
    samp = samples()
    par(xpd = T)
    barplot(c(
      sum(samp$L1x/samp$L0x >= input$pow/input$alpha),
      sum(samp$L1x/samp$L0x < input$pow/input$alpha & samp$L1x/samp$L0x >= 3),
      #sum(samp$L1x/samp$L0x < 3 & samp$L1x/samp$L0x >= 1),
      sum(samp$L1x/samp$L0x < 3)),
      names.arg = c(paste("LR >=", round(input$pow/input$alpha,2)), paste(input$pow/input$alpha,">LR>=3"), "3>LR"),
      ylim=c(0,Nsampl))
    
    legend(0,-19,bty = "n",legend= c( "L(d=delta|x)/L(d=0|x)","L(d>delta|x)/L(d=0|x)"), pch =c(16,16), col = c("grey", rgb(1,0,0,alpha=0.3)))
    barplot(c(sum(samp$Lplus0x/samp$L0x >= input$pow/input$alpha),
              sum(samp$Lplus0x/samp$L0x < input$pow/input$alpha & samp$Lplus0x/samp$L0x >= 3), 
              #sum(samp$Lplus0x/samp$L0x < 3 &samp$Lplus0x/samp$L0x >= 1),
              sum(samp$Lplus0x/samp$L0x < 3)),
            ylim=c(0,Nsampl),
            col = rgb(1,0,0,alpha=0.3),
            add = TRUE)
  })
  
  output$plotneg <- renderPlot(res = 100, {     
    
    samp = samples0()
    par(xpd = TRUE)
    barplot(c(
      sum(samp$L1x/samp$L0x >= input$pow/input$alpha),
      sum(samp$L1x/samp$L0x < input$pow/input$alpha & samp$L1x/samp$L0x >= 3),
      # sum(samp$L1x/samp$L0x < 3 & samp$L1x/samp$L0x >= 1),
      sum(samp$L1x/samp$L0x < 3)),
      names.arg = c(paste("LR >=", round(input$pow/input$alpha,2)), 
                    paste(round(input$pow/input$alpha,2),">LR>=3"), 
                    "1>LR"),
      ylim=c(0,Nsampl))
    legend(0,-19,bty = "n", legend= c( "L(d=delta|x)/L(d=0|x)","L(d>delta|x)/L(d=0|x)"), pch =c(16,16), col = c("grey", rgb(1,0,0,alpha=0.3)))
    
    barplot(c(sum(samp$Lplus0x/samp$L0x >= input$pow/input$alpha),
              sum(samp$Lplus0x/samp$L0x < input$pow/input$alpha & samp$Lplus0x/samp$L0x >= 3), 
              #sum(samp$Lplus0x/samp$L0x < 3 &samp$Lplus0x/samp$L0x >= 1),
              sum(samp$Lplus0x/samp$L0x < 3)),
            ylim=c(0,Nsampl),
            col = rgb(1,0,0,alpha=0.3),
            add = TRUE)
    
  })
  
  
  
  
  output$tablePosL <- renderTable({
    samp = samples()
    samp30 = samples30()
    df = 2*samp$Nest-2
    pow = pnorm(sqrt(samp$Nest*1.5)*input$d/sqrt(2)-qnorm(1-input$alpha))
    
    diff = mean(samp$Lplus0x/samp$L0x >=input$pow/input$alpha) - round(sum(samp$ddeviate[samp$L1x/samp$L0x >=input$pow/input$alpha] < dt( qt(.5, df,samp$ncp),df,samp$ncp)/dt( qt(.95, df,samp$ncp),df,samp$ncp))/100,2)
    
    LR = c(mean(samp$Lplus0x/samp$L0x >=input$pow/input$alpha),
           mean(samp$L1x/samp$L0x>=input$pow/input$alpha),
           round(sum(samp$ddeviate[samp$L1x/samp$L0x >=input$pow/input$alpha] < dt( qt(.5, df,samp$ncp),df,samp$ncp)/dt( qt(.95, df,samp$ncp),df,samp$ncp))/100,2),
           diff,
           mean((log(samp$L1x/samp$L0x)+log(samp30$L1x/samp30$L0x))>= log(pow/(1-pow)))
    )
    
    
    
    LR = data.frame(cbind(c("4. Substantial Falsification", 
                            "5. Preliminary Verification", 
                            "6. Substantial Verification",
                            "False Negatives: Substantial Falsification - Substantial Verification", 
                            paste("Substantial Verification if N'= N +N/2",  round(samp$Nest*.5),", 1-beta=",round(pow,2)))
                          ,LR))
    
    colnames(LR) = c("Step in RPS","Proportion")
    LR
  })
  
  
  output$tablePosLEmp <- renderTable({
    samp = samplesEmp()
    
    LR = c(#mean(samp$Lplus0x/samp$L0x >=99),
      mean(samp$Lplus0x/samp$L0x >=samp$pow/input$alphaemp),
      #mean(samp$L1x/samp$L0x>99),
      mean(samp$L1x/samp$L0x>=samp$pow/input$alphaemp))
    
    LR = data.frame(cbind(c("4. Substantial Falsification",  "5. Preliminary Verification")
                          ,LR))
    colnames(LR) = c("Step in RPS", "Proportion")
    LR
  })
  
  output$tableNegLEmp <- renderTable({
    
    samp = samplesEmp0()
    sampP = samplesEmp()
    
    
    
    LR = c(#mean(samp$Lplus0x/samp$L0x >=99),
      mean(samp$Lplus0x/samp$L0x >=sampP$pow/input$alphaemp),
      #mean(samp$L1x/samp$L0x>=99),
      mean(samp$L1x/samp$L0x>=sampP$pow/input$alphaemp))
    
    LR = data.frame(cbind(c("4. Substantial Falsification", "5. Preliminary Verification")
                          ,LR))
    colnames(LR) = c("Step in RPS", "Proportion")
    LR 
    
  })
  
  
  output$LRextrem <- renderTable({
    samp <-  samples()
    samp0 <-  samples0()
    LR = data.frame(cbind(min(samp$L1x/samp$L0x),quantile(samp0$L1x/samp0$L0x, .99)))
    colnames(LR) = c("min(L(d|x)/L(d=0|x)) if H1 is true", "max(L(d|x)/L(d=0|x)) if H0 is true")
    LR
  })
  
  output$tableNegL <- renderTable({
    samp = samples0()
    LR = c(mean(samp$Lplus0x/samp$L0x>=input$pow/input$alpha),
           mean(samp$L1x/samp$L0x>=input$pow/input$alpha))
    
    LR = data.frame(cbind(c("4. Substantial Falsification", "5. Preliminary Verification")
                          ,LR))
    colnames(LR) = c("Step in RPS", "Proportion")
    LR
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)

