#' @title Outlier detection for single or grouped time series
#'
#' @description This function provides anomaly signals (even a graphical visualization) when there is a 'jump' in a single time series, or the 'jump' is too much different respect those ones of grouped similar time series.
#'
#' @param dati data frame (grouped time series: phenomenon+date+group+values) or vector (single time series)
#' @param graph logical value for graphical analysis (default=FALSE)
#' @param linear_analysis logical value for linear analysis (default=FALSE)
#' @param val_test_limit value for outlier detection sensitiveness (default=5 ; max=10)
#' @param save_out logical value for saving detected outliers (default=FALSE)
#' @param out_out a character file name for saving outliers in csv form, delimited with ";" and using ',' as decimal separator  (default out.csv)
#' @param pdf_out a character file name for saving graphic analysis in pdf file (default=out.pdf)
#' @param r_out rows number of graphs (default=3)
#' @param c_out cols number of graphs (default=2)
#' @param first_line value for first dotted line in graphic analysis (default=1)
#' @param pace_line value for pace in dotted line in graphic analysis (default=6)
#' @return Data frame of possible outliers in a triad. Output record: rows/time.2/series/y1/y2/y3/test(AV)/AV/ n/median(AV)/mad(AV)/madindex(AV). Where time.2 is the center of the triad y1, y2, y3; test(AV) is the number to compare with 5 to detect outlier; n is the number of observations of the group ....
#' @export
#' @examples
#' ## we can start with data without outliers but structured with co-movement between groups
#'data("dati")
#'## first column for phenomenon
#'## 2° col for time written in ordered numbers or strings
#'## 3° col for group classification variable
#'## 4° col for values
#'str(dati)
#'#######################################
#'## a data frame without any outlier
#'#######################################
#'out=wash.out(dati)
#'out   ## empity data frame
#'length(out[,1])  ## no row
#'## we can add two outliers
#'####  time=3 temperature value=0
#'dati[99,4]=  0
#'## ... and then for "rain" phenomenon!
#'####  time=3 rain value=37
#'dati[118,4]=  37
#'#######################################
#'##   data.frame with 2 fresh outliers
#'#######################################
#'out=wash.out(dati)
#'##  all "three terms" time series
#'## let's take a look at anomalous time series
#'out
#'## ... the same but we save results in a file....
#'## If we don't specify a name, out.csv  is the default
#'out=wash.out(dati,save_out=TRUE,out_out="tabel_out.csv")
#'out
#'## we put the parameter from 5 to 10, using this upper one  to capture
#'##       only  particularly anomalous outliers
#'out=wash.out(dati, val_test_limit = 10)
#'out
#'## save plots and outliers in a pdf file "out.pdf" as a default
#'out=wash.out(dati, val_test_limit = 10, graph=TRUE)
#'out
#'## we can make the usual analysis for groups but we can also use that one
#'## reserved for every single time series
#'## (linear_analysis): two files for saved outliers (out.csv and linout.csv)
#'##  and for graph display in two pdf files (out.pdf and linout.pdf)
#'out=wash.out(dati,val_test_limit=5,save_out=TRUE,linear_analysis=TRUE,graph=TRUE)
#'out
#'## out return only the linear analysis...
#'## ... in this case we lose the co-movement information an we run the risk
#'##     of finding too much variance in a single time series
#'##     and detecting not too much likely outliers
#'##########################################################
#'##  single time series analysis
#'##########################################################
#'data(ts)
#'str(ts)
#'sts= ts$dati
#'plot(sts,type="b",pch=20,col="red")
#'## a time series with a variability and an increasing trend
#'## sts is a vector and linear analysis is the default one
#'out=wash.out(sts)
#'out
#'## we find no outlier
#'out=wash.out(sts,val_test_limit=5,linear_analysis=TRUE,graph=TRUE)
#'out
#'## no outlier
#'## We can add an outlier with limited amount
#'sts[5]=sts[5]*2
#'plot(sts,type="b",pch=20,col="red")
#'out=wash.out(sts,val_test_limit=5)
#'out
#'## test is over 5 for a bit
#'out=wash.out(sts,val_test_limit=5,save_out=TRUE,graph=TRUE)
#'out
#'data(ts)
#'sts= ts$dati
#'sts[5]=sts[5]*3
#'## we can try a greater value to put an outlier of a certain importance
#'plot(sts,type="b",pch=20,col="blue")
#'out=wash.out(sts,val_test_limit=5,save_out=TRUE,graph=TRUE)
#'out
#'## washer procedure identify three triads of outliers values
#'system("rm *.csv *.pdf")
#'

wash.out = function( dati                      ,
                   #   p      t     i     y
                   # dati structure:  phenom./date/series/values/... other
                   graph=FALSE                 ,
                   linear_analysis=FALSE       ,
                   val_test_limit = 5          ,
                   save_out=FALSE              ,
                   out_out="out.csv"           ,
                   pdf_out="out.pdf"           ,
                   r_out =3                    ,
                   c_out=2                     ,
                   first_line =1               ,
                   pace_line = 6
)
{
  ## start function code
  ## sub function recall


  washer2.AV = function( dati ) #   p      t     i     y
  {          # dati structure:  phenom./date/series/values/... other
    # example:    Phenomenon     Time    Zone    Value    ...
    #             -----------  --------   --     -----  --------
    #             Temperature  20091231   A1      20.1    ...
    #             Temperature  20091231   A2      21.0    ...
    #                             ...
    #             Rain         20081231   B1     123.0    ...
    #                             ...
    ###############################################################################################
    AV      =  function(y) {   # y matrix 3 columns (y1 y2 y3) and n rows
      AV=array(0,length(y[,1]))
      100*(2*y[,2]-y[,1]-y[,3])/(stats::median(y[,1]+y[,2]+y[,3])+ y[,1]+y[,2]+y[,3]) }
    # output array AV
    ###############################################################################################
    test.AV =  function(AV) {  # AV array n rows
      t(rbind(test.AV=abs(AV-stats::median(AV))/stats::mad(AV),AV=AV,n=length(AV),median.AV=stats::median(AV),mad.AV=stats::mad(AV) ,
              madindex.AV=stats::mad(AV)*1000/150  ))     }
    # col      1      2   3        5          6         7
    # output: test / AV / n /  median(AV) / mad(AV) / madindex
    ################################################################################################
    if (min(dati[,4])> 0) {
      dati=dati[which(!is.na(dati[,4])),]
      dati=dati[order(dati[,1],dati[,3],dati[,2]),]
      fen=rownames( table(dati[,1]) )
      nfen=length(fen)
      out= NA
      for ( fi in 1:nfen)
      { print(c("phenomenon:",fi) ,quote=FALSE)
        time=rownames( table(dati[which(fen[fi]==dati[,1]),2]) )
        n=length(time)
        for ( i in 2:(n-1) )
        { datiy= dati[dati[,1] == fen[fi] & dati[,2] %in% c(time[i-1],time[i],time[i+1]),c(2,3,4)]
        y1=stats::reshape(datiy,timevar=colnames(datiy)[1],idvar=colnames(datiy)[2],direction="wide" )
        y1=y1[!is.na(apply(y1[,c(2:4)],1,sum)),]
        y=y1[,c(2,3,4)]
        colnames(y)=c("t.1","t.2","t.3")
        out=rbind(out,data.frame(fen=fen[fi],t.2=time[i],
                                 series=y1[,1],y=y,test.AV(AV(y))))
        }
      }
      rownames(out)=(1:length(out[,1])-1)
      washer2.AV=out[2:length(out[,1]),]
      # col      1      2      3     4  5  6      7    8  9     10         11       12
      # output: rows /time.2/series/y1/y2/y3/test(AV)/AV/ n /median(AV)/mad(AV)/madindex(AV)
      # end function washer2.AV
    } else print(" . . . zero or negative y:  t r a n s l a t i o n   r e q u i r e d !!!")
  }


  ####################################################################################################





  td_w= function(dati, phen="T")
  { ## serie storica a gruppi di tre
    ## ----------------------------------
    time=rep(0,(length(dati)-2)*3)
    zone=rep(0,(length(dati)-2)*3)
    value=rep(0,(length(dati)-2)*3)
    out = data.frame(phen,time,zone,value)

    k=1
    for (i in 1:(length(dati)-2))
    {for( j in 1:3)
    {
      out$zone[k+j-1]  = i
      out$time[k+j-1]  = j
      out$value[k+j-1] = dati[i+j-1]
    }
      k=k+3
    }
    return(out)
  }


  translate = function( dati )
  {
    min_var=min(dati[,4],na.rm = TRUE)
    if (min_var>0) min_var = 0.0001
    dati[,4]= dati[,4] - min_var+0.0001
    return(dati)
  }

  inv_translate = function (out,dati)
  {
    min_var=min(dati[,4],na.rm = TRUE)
    if (min_var>0) min_var = 0.0001

    out$y.t.1 = out$y.t.1 +min_var -0.0001
    out$y.t.2 = out$y.t.2 +min_var -0.0001
    out$y.t.3 = out$y.t.3 +min_var -0.0001
    return(out)
  }







  graph_f = function(dati,out,pdf_out,r_out,c_out,first_line,pace_line,val_test_limit)
  {

    #requireNamespace(gplots)
    grDevices::pdf(pdf_out)

    graphics::par(mfrow=c(r_out,c_out),cex.main=.6)

    a= out[out[,"test.AV"]> val_test_limit,]
    a= a[order(paste(a$fen,a$series)),]

    series_old = 0
    fen_old =  0


    for (i in 1:length(a[,1]))


    {
      #######################
      series=a[i,"series"]
      fen=a[i,"fen"]
      #######################

      if( !(series == series_old & fen == fen_old))

      {

        elenco=a[a$series==a[i,"series"] & a$fen==a[i,"fen"],c(1,3,2,4,5,6,7)]
        elenco[4:7]=round(elenco[4:7],2)
        #    names(elenco)[1]= "settore"
        #    names(elenco)[2] = "area"
        if (length(elenco[,1]) > 25) elenco = elenco[1:25,]
        gplots::textplot(elenco,cex = .50 , show.rownames = FALSE)
        graphics::title("Outlier(s) list \n (max 25)")

        ##
        cond  = (  dati[,3]== series & dati[,1]== fen )
        ss <- dati[ cond , c(2,4)]
        ss=ss[order(ss[,1]),]
        #ss[,1]=as.integer(ss[,1])
        #ss


        x=1:length(ss[,1])
        time=ss[,1]
        y=ss[,2]
        graphics::plot(x,y,xaxt="n",type="l",pch=20, cex=0.5,col="blue",
             main=paste('Series:',a[i,"series"],' phen. ',fen),
             lwd=0.4)
        graphics::points(x,y,col="darkblue",pch=20, cex=0.5)
        graphics::axis(1, at=seq(first_line, length(ss[,1]), by = pace_line),
             labels=time[seq(first_line, length(ss[,1]), by = pace_line)],
             las=2,cex.axis=0.8)

        graphics::abline(h=0, col = "blue",lty = 3)
        graphics::abline(v=seq(first_line-0.1, length(ss[,1]), by = pace_line), col = "blue",lty = 3)



        series_old = series
        fen_old =  fen

      }



      out_t2=which(ss[,1]==a[i,"t.2"])
      graphics::abline(v=out_t2-1.1,col="red", lwd=0.1)
      graphics::abline(v=out_t2+1.1,col="red", lwd =0.1)


    }





    grDevices::dev.off() #

  }











  if(is.vector(dati))
  { linear_analysis=TRUE
  if(length(dati)<10)
  {
    print("No linear analisys: too few data!")
    return(0)
  }
  if (min(dati)>=0 ) min_var = 0.0001
  else min_var = min(dati)
  dati= dati - min_var+0.0001
  out= washer2.AV(td_w(dati))
  out=out[out[,7]>val_test_limit,]
  if( length(out[,1])==0) print("NO outlier!!!")
  dum=out[,2]
  out[,2]=out[,3]+1
  out[,3]=dum
  out$y.t.1 = out$y.t.1 +min_var -0.0001
  out$y.t.2 = out$y.t.2 +min_var -0.0001
  out$y.t.3 = out$y.t.3 +min_var -0.0001
  dati= dati +min_var-0.0001
  if(save_out & length(out[,1])>0){
    utils::write.csv2(out[out[,7]>val_test_limit,],file=out_out)
    print(paste("File '",out_out,"' saved!",sep=""))
    print(paste("Dir -> ",getwd()))
  }

  if(graph & length(out[,1])>0)
  {
    dati=data.frame(phen=rep("T",length(dati)),data=1:length(dati),series=rep(2,length(dati)),dati=dati)
    #print(dati)
    #print(out)
    graph_f(dati,out,paste("lin",pdf_out,sep=""),r_out=2,c_out=1,first_line,pace_line,val_test_limit)
    print(paste("File'",paste("lin",pdf_out,sep=""),"' saved!",sep=""))
    print(paste("Dir -> ",getwd()))
  }
  return(out[out[,7]>val_test_limit,])
  }
  else  if(is.data.frame(dati))
  {
    out= washer2.AV(translate(dati))
    out=out[out[,7]>val_test_limit,]
    if( length(out[,1])==0) print("NO outlier!!!")
    else  out=inv_translate(out,dati)
    if(save_out & length(out[,1])>0){
      utils::write.csv2(out[out[,7]>val_test_limit,],file=out_out)
      print(paste("File '",out_out,"' saved!",sep=""))
      print(paste("Dir -> ",getwd()))
    }
    if(graph & length(out[,1])>0)
    {
      graph_f(dati,out,pdf_out,r_out,c_out,first_line,pace_line,val_test_limit)
      print(paste("File '",pdf_out,"' saved!",sep=""))
      print(paste("Dir -> ",getwd()))
    }
    if(linear_analysis)
    {
      fen_series=dati[,c(1,3)]
      datiw=dati[,c(1,2,3,4)]
      datiw[,1] = paste(datiw[,1],datiw[,3], sep="")
      fen_series[,3]=paste("p_",datiw[,1],sep="")
      datiw[,3] <- NULL
      names(datiw)[1]="col"
      names(datiw)[2]="DATA"
      names(datiw)[3]="p"

      datiw=stats::reshape(datiw,v.names="p",idvar="DATA",timevar="col",direction="wide",sep="_")
      n_col = length(datiw[1,])
      print(paste ("cols:",n_col))
      n_row = length(datiw[,1])
      if(n_row<10)
      {
        print("No linear analisys: too few data!")
        return(out[out[,7]>val_test_limit,])
      }
      print(paste ("rows:",n_row))
      mesi = datiw$DATA[c(-1,-n_row)]
      ###############################################################
      ##  first elaboration out of cycle
      ###############################################################
      min_var= min(datiw[,2],na.rm = TRUE)
      datiw[,2]= datiw[,2] - min_var +0.0001
      out=washer2.AV(td_w(datiw[,2]))
      out$fen= colnames(datiw)[2]
      out$t.2=mesi
      out$y.t.1 = out$y.t.1 +min_var -0.0001
      out$y.t.2 = out$y.t.2 +min_var -0.0001
      out$y.t.3 = out$y.t.3 +min_var -0.0001
      datiw[,2] = datiw[,2] +min_var -0.0001

      ###############################################################
      for (i in 3:n_col)
      ###############################################################
      {
        min_var= min(datiw[,i],na.rm = TRUE)
        datiw[,i]= datiw[,i] - min_var +0.0001
        out0=washer2.AV(td_w(datiw[,i]))
        out0$fen= colnames(datiw)[i]
        out0$t.2=mesi
        out0$y.t.1 = out0$y.t.1 +min_var -0.0001
        out0$y.t.2 = out0$y.t.2 +min_var -0.0001
        out0$y.t.3 = out0$y.t.3 +min_var -0.0001
        datiw[,i] = datiw[,i]   +min_var -0.0001
        out=rbind(out,out0)
        print(paste("Time series: ",i))

      }
      ###############################################################

      for (i in 1:length(out[,1]))
      {pos=match(out[i,1],fen_series[,3])
      as.character(fen_series[pos,1]) -> out[i,1]
      as.character(fen_series[pos,2]) -> out[i,3]
      }
    }
    out=out[out[,7]>val_test_limit,]
    if(linear_analysis & save_out & length(out[,1])>0){
      utils::write.csv2(out[out[,7]>val_test_limit,],file=paste("lin",out_out,sep=""))
      print(paste("File 'lin",out_out,"' saved!",sep=""))
      print(paste("Dir -> ",getwd()))
    }
    if(linear_analysis & graph & length(out[,1])>0)
    {
      graph_f(dati,out,paste("lin",pdf_out,sep=""),r_out,c_out,first_line,pace_line,val_test_limit)
      print(paste("File '",paste("lin",pdf_out,sep=""),"' saved!",sep=""))
      print(paste("Dir -> ",getwd()))
    }


  }
  else print ("not compliant data ")

  #dev.off() #
  return(out)
  ## end function
}
