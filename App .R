#Libraries needed
library(shiny)
library(ggplot2)
library(gridExtra)
library(matrixStats)

# Define UI for application that draws a histogram
ui <- 
  
  fluidPage(headerPanel("Poisson Process"),
            
            #Define user input parameters   
            
            sidebarPanel(
              
              
              #helpText('Events occurrence plot'), 
              numericInput("lambda", "Enter P.P rate (lambda)",0.5),
              numericInput("timee", "Enter time period",20),
              numericInput("Trials","Enter Number of trials for simulations", 500),
              
              
              HTML("<p style='color:blue;'>S(i) distribution parameters:</p>"),
              numericInput("Sdistribution", "Enter a number as i",1),
              tags$em("To prove this property simulated CDF of S(i) has been campared with CDF of Gamma distribution [i,lambda]"),
              
              HTML("<p style='color:blue;'>Poisson Thinning parameters:</p>"),
              numericInput("EventType1", "Enter Probability of Event Type1",0.6),
              tags$em("Probability of event Type2=1-(Probability of event Type1)"),
              
              HTML("<p style='color:blue;'>Compound P.P parameters:</p>"),
              numericInput("lambdacompound", "Enter a number as rate of exponential distribution",5),
              
              
              HTML("<p style='color:blue;'>Supperposition parameters:</p>"),
              numericInput("superposition1", "Enter rate for P.P1",3),
              numericInput("superposition2", "Enter rate for P.P2",3)
             
      ),
      
      
    
    mainPanel(plotOutput("hist"),plotOutput("hist2"),plotOutput("hist3"),plotOutput("hist4"),plotOutput("hist5")),tableOutput("table1"),
    tags$footer("R Shiny application created by:", tags$b("Vahid Azizi and Samira KarimZadeh"),tags$br() ,
            "Department of Industrial and Manufacturing Systems Engineering",tags$br(),
            "Iowa State University",
            align = "left", style = "
            bottom:0;
            width:100%;
            height:0px;   /* Height of the footer */
            color: black;
            padding: 10px;
            z-index: 1000;"))

  
server<-function(input, output){
  
  output$hist<-renderPlot({
    
    x=c(0)
    y=c(0)
    while(sum(x)<input$timee){
      a=rexp(1,rate=input$lambda)
      x=c(x,a)
      y=c(y,1)
      
    }
    
    x<- x[1:(length(x)-1)]
    y <- c(0,y[1:(length(y)-1)])
    
    x=as.numeric(cumsum(x))
    y=as.numeric(cumsum(y))
    limitt=as.numeric( c(0,input$timee))

    plot(stepfun(x,y), main = "Events occurrence based on P.P",xlab = "Time",ylab="N(t)",xlim = limitt,col = "red", do.points = F)
    
    })
  
  
  
  
 
  
  output$hist2<-renderPlot({
    
    SimlistS=rep(0,input$Trials)
    
    for (i in 1:input$Trials){
      
      lambda2=input$lambda
      Stimes=rexp(input$Sdistribution,lambda2)
      Sfinal=sum(Stimes)
      SimlistS[i]=Sfinal
    }
    
    steps=100
    x<-rep(0,steps)
    y<-rep(0,steps)
    Gamma<-rep(0,steps)
    
    for (i in 2:steps){
      
      x[i]=(i*input$timee)/steps
      yy=(length(which(SimlistS<x[i])))/input$Trials
      y[i]=yy
      Gamma[i]=pgamma(x[i], input$Sdistribution, lambda2 ,lower.tail = TRUE)
      
    }

    timeunits=c(x,x)
    y=as.numeric(y)
    SimulatedValues =y
    
    Probability <-c(SimulatedValues,Gamma )
    Legend <-c(rep("Simulated CDF of S(i)",steps),rep("CDF of Gamma", steps))
    mydata2 <-data.frame(x, Probability)
    
     plot =  ggplot(data=mydata2, aes(x=timeunits, y=Probability, group=Legend,colour =Legend)) +
        geom_line()+
        geom_point()+
        labs(title = "S(i) Distribution")+
        xlab("Time") + ylab("Probability")+
        theme(plot.title = element_text( face="bold", size=14))

      print(plot)

  })
  
  output$hist3<-renderPlot({
    
    simlist1=c(0)
    simlist2=c(0)
    lambda2=input$lambda*input$timee
    
    for (i in 1:input$Trials){
      
      p=runif(1,0,1)
      numEvent=rpois(1,lambda2)
      
      if (p<=input$EventType1)
      {
        
        simlist1=c(simlist1,numEvent)}
      
      else
      {
        
        simlist2=c(simlist2,numEvent)}
      
    }
    
    total<-c(simlist1,simlist2)
    maxx<-max(total)
    simlist111<-rep(0,maxx+1)
    simlist222<-rep(0,maxx+1)
    exact<-rep(0,maxx+1)
    
    for (i in 1:(maxx+1)){
      
      w=length(which(simlist1==i-1))
      q=length(which(simlist2==i-1))
      simlist111[i]=w
      simlist222[i]=q
      exact[i]=dpois(i-1,lambda2)
    }
    
    probsimlist111=simlist111/length(total)
    probsimlist222=simlist222/length(total)
    EventNum<-c(0,c(1:maxx))
    Type1 <-probsimlist111
    Type2 <-probsimlist222
    Type1PlusType2<-probsimlist111+probsimlist222
    Probability <-c(Type1PlusType2, exact )
    Legend <-c(rep("Type1PlusType2",maxx+1),rep("Exact",maxx+1))
    mydata <-data.frame(EventNum, Probability)
    
    
    
    plot =  ggplot(data=mydata, aes(x=EventNum, y=Probability , group=Legend,colour =Legend ,fill = Legend)) +
      
      geom_bar(stat = "identity", position = "dodge")+
      xlab("Number of Events") + ylab("Probability") +
      ggtitle("Poisson Thinning") +
      theme(plot.title = element_text( face="bold", size=14))
    
    
    
    print(plot)
    
  })
  
  output$hist4<-renderPlot({
    
    repy<-rep(0,10* input$lambdacompound * input$timee)
    xsum<-rep(0,10* input$lambdacompound * input$timee)

    simlist_compound = numeric(input$Trials)
    for(i in 1:input$Trials){
      
      N=rpois(1,input$lambda*input$timee)
      x=rexp(N,rate=input$lambdacompound)
      repy[N]=repy[N]+1
      xsum[N]=xsum[N]+sum(x)
      simlist_compound[i] = sum(x)
      
    }
    #confidence interval
   # mean_sim = sum(xsum)/input$Trials
   # v = var(xsum)
    mean_sim = mean(simlist_compound)
    v = var(simlist_compound)
    interval = 1.96*sqrt(v)/sqrt(length(xsum))
    
    index=max(which(repy>0))
    
    repy<- repy[1:index]
    xsum<- xsum[1:index]
    expectedx<-rep(0,length(xsum))
    
    #calculating expected valu of X conditionening on value of N=n
    
    for (j in 1:length(xsum)){
      if (repy[j]!=0){
        expectedx[j]=(xsum[j]/(repy[j]))
        }
    }
    
    h=input$lambda*input$timee*(1/input$lambdacompound)
    exactmeanYi<-rep(h,length(xsum))
    pp=sum(xsum)/input$Trials
    meanYi<-rep(pp,index)
    xx=length(repy)
    steps2=c(1:length(repy))
    m=length(steps2)
    timeunits2=c(steps2,steps2,steps2)
    exactmeanYi=as.numeric(exactmeanYi)
    meanYi =as.numeric(meanYi)
    expectedx=as.numeric(expectedx)
    Expexted <-c(expectedx,meanYi,exactmeanYi )
    
    Legend <-c(rep("Simulated E(X|N=n)",m),rep("Simulated E(X)", m),rep("Exact E(X)", m))
    mydata3 <-data.frame(steps2, Expexted)

 
        
    plot = ggplot(data=mydata3, aes(x=timeunits2, y=Expexted, group=Legend, colour =Legend)) +
      geom_line()+
      geom_point()+
      labs(title = "Compound Poisson Process")+
      xlab("Number of Events") + ylab("Value")+
      theme(plot.title = element_text( face="bold", size=14))

    pasteinterval = paste("95% C.L. E(x) is ", round(mean_sim,2),"+-",round(interval,3) , "and Exact value for E(x) is ", exactmeanYi)   
    plot = gridExtra::grid.arrange(plot,bottom = pasteinterval)
    plot = gridExtra::grid.arrange(plot,bottom = "")    
    print(plot)

  })
  
  
  output$hist5<-renderPlot({  
    
    simlist_sum = matrix(nrow = input$timee, ncol = input$Trials)
    for(i in 1: input$Trials){
      for(j in 1:input$timee){
        
        #simlist_sum[j,i] = rpois(1, (input$superposition1 + input$superposition2)*j)
        simlist_sum[j,i] = rpois(1, input$superposition1 *j)+rpois(1,  input$superposition2*j)
      }
    }
    
    Nt=rowMeans(simlist_sum)
    Ns = numeric(input$timee)
    for(j in 1:input$timee){
      Ns [j] =  (input$superposition1+input$superposition2)*j
    }
    Event_time<-c(1:input$timee)
    Event_Num <-c(Nt,Ns)
    Legend <-c(rep("E[N(t)]", length(Nt)),rep("E[N1(t)+N2(t)]", length(Ns)))
    mydata4 <-data.frame(Event_time, Event_Num)
    t_unit = c(Event_time,Event_time)
    
    ######Confidence interval 
    m_sum = rowMeans(simlist_sum)
    sd_sum = apply(simlist_sum,1,sd)
    interval_sum = 1.96*sd_sum/sqrt(input$Trials)
    pasteinterval = paste("95% C.L. for E[N(t)] is as table below")
    
    plot = ggplot(data = mydata4, aes(x=t_unit ,y  =Event_Num, group=Legend,colour =Legend)) +
      geom_line()+
      geom_point()+
      labs(title = "Super position")+
      xlab("Time") + ylab("E[N(t)]")+
      theme(plot.title = element_text( face="bold", size=14))
    
      
    
    
    plot = gridExtra::grid.arrange(plot,bottom = pasteinterval)
    plot = gridExtra::grid.arrange(plot,bottom = "")
    print(plot)
    in_interval = logical(input$timee)
    for(i in 1:input$timee){
      in_interval[i] =abs( Nt[i] - Nt[i] ) <= interval_sum[i]
      
    }

    superpos_table = matrix(c("MeanExact" , round(Ns,3),"MeanSim" , round(m_sum,3) , "Interval (+-)" , round(interval_sum,3), "Fits Interval" , in_interval),nrow = 4, byrow = T)
    colnames(superpos_table )=c("Time" , c(1:input$timee))
    
    output$table1<-renderTable(superpos_table)
    
  })
  

}


shinyApp(ui=ui, server=server)