###### Getting started:
# Copy - Paste the whole collection function in the command windows of R (CRAN)
# If necessary, download the two required libraries: tcltk and combinat
#
# Launch SIMIL (if not already launched) with the command: SIMIL.GUI()
#
###### Usage:
#
# Structure provides outputs without populations nomenclature, that's why
#
# YOU MUST PROVIDE A '.pop' file that contains populations names, see example files!
#
# The '.pop' file is a tab-delimited file and can be edited with a simple text editor.
# It must have a filename that contains the first three letters
# of the related Run of structure.
#
# e.g. Run = CspAlps_AC_Csp!_f.14
# '.pop' file = CspAlps_pops.pop or Csp_pops.pop or Csp.pop
#

require(tcltk)
require(combinat)

############## SIMIL calculation
SIMIL=function(...){
  COMBIN=function(data1.p,data2.p,permu=t(array(unlist(permn(nc)),dim=c(nc,gamma(nc+1)))),...){
    funfun=function(X=1,...){
    ranks=permu[X,]
    sum(abs(data1.p-data2.p[,order(ranks)]))
    }
    nc=ncol(data2.p)
    delta=sapply(c(1:gamma(nc+1)),funfun)
    data2.p=data2.p[,order(permu[which.min(delta),])]
    data2.p
    }

  #### Import Structure Runs
  speciesA=IMPSTRUCT.P(choose.files(default = "", caption = "Select Structure Run Species A",multi = TRUE, index = nrow(Filters)))
  speciesB=IMPSTRUCT.P(choose.files(default = "", caption = "Select Structure Run Species B",multi = TRUE, index = nrow(Filters)))

  ### Determine Overlaps and cut datasets
  popsA=speciesA$CodePop
  popsB=speciesB$CodePop
  inters.A=speciesA[match(intersect(popsA,popsB),popsA),colnames(speciesA)=='P']
  inters.B=speciesB[match(intersect(popsA,popsB),popsB),colnames(speciesB)=='P']

  if(is.null(dim(inters.A))){
    stop(call='These runs have K=1 groups, function stopped!')
    }

  n.overlap=nrow(inters.A)
  n.speciesA=nrow(speciesA)
  n.speciesB=nrow(speciesB)
  n.overall=n.speciesA+n.speciesB-n.overlap

  if((100*n.overlap/n.overall)<10){
    usr.choice=winDialog('yesno','Less than 10% of overlap between SpeciesA and SpeciesB. Continue anyway?')
    if(usr.choice=='NO') stop(call='Function stopped.')
    }

  ### find best homology between matrices
  inters.B=try(COMBIN(inters.A,inters.B),silent=T)

    Delta=abs(inters.A-inters.B)
    Delta=1-mean(apply(Delta,1,mean))

    ### Final results
    GSS.i=Delta
    GSS.s=Delta*n.overlap/min(n.speciesA,n.speciesB)
    GSS.o=Delta*n.overlap/n.overall
    GSS.m=mean(c(GSS.s,GSS.o))

    output=list(GSS=data.frame('GSS.unweighted'=GSS.i,'GSS.smallest'=GSS.s
              ,'GSS.overall'=GSS.o,'GSS.mean'=GSS.m
              ,'PercentOverlap'=100*n.overlap/n.overall),
         SpB.shuffled=inters.B)
    assign('GSS.out',output,env = .GlobalEnv)
    output$GSS
  }
  
############## SIMIL calculation
SIMIL.R=function(pathA,pathB){
  COMBIN=function(data1.p,data2.p,permu=t(array(unlist(permn(nc)),dim=c(nc,gamma(nc+1)))),...){
    funfun=function(X=1,...){
    ranks=permu[X,]
    sum(abs(data1.p-data2.p[,order(ranks)]))
    }
    nc=ncol(data2.p)
    delta=sapply(c(1:gamma(nc+1)),funfun)
    data2.p=data2.p[,order(permu[which.min(delta),])]
    data2.p
    }

  #### Import Structure Runs
  speciesA=IMPSTRUCT.P(pathA)
  speciesB=IMPSTRUCT.P(pathB)

  ### Determine Overlaps and cut datasets
  popsA=speciesA$CodePop
  popsB=speciesB$CodePop
  inters.A=speciesA[match(intersect(popsA,popsB),popsA),colnames(speciesA)=='P']
  inters.B=speciesB[match(intersect(popsA,popsB),popsB),colnames(speciesB)=='P']

  if(is.null(dim(inters.A))){
    stop(call='These runs have K=1 groups, function stopped!')
    }

  n.overlap=nrow(inters.A)
  n.speciesA=nrow(speciesA)
  n.speciesB=nrow(speciesB)
  n.overall=n.speciesA+n.speciesB-n.overlap

  if((100*n.overlap/n.overall)<10){
    usr.choice=winDialog('yesno','Less than 10% of overlap between SpeciesA and SpeciesB. Continue anyway?')
    if(usr.choice=='NO') stop(call='Function stopped.')
    }

  ### find best homology between matrices
  inters.B=try(COMBIN(inters.A,inters.B),silent=T)

    Delta=abs(inters.A-inters.B)
    Delta=1-mean(apply(Delta,1,mean))

    ### Final results
    GSS.i=Delta
    GSS.s=Delta*n.overlap/min(n.speciesA,n.speciesB)
    GSS.o=Delta*n.overlap/n.overall
    GSS.m=mean(c(GSS.s,GSS.o))

    output=list(GSS=data.frame('GSS.unweighted'=GSS.i,'GSS.smallest'=GSS.s
              ,'GSS.overall'=GSS.o,'GSS.mean'=GSS.m
              ,'PercentOverlap'=100*n.overlap/n.overall),
         SpB.shuffled=inters.B)
    assign('GSS.out',output,env = .GlobalEnv)
    output$GSS
  }
  
############## Importation Individuals
IMPSTRUCT.I=function(output){
  dd=readLines(output)
  ### Main
  Indiv=grep('individuals',dd)
  Clust=grep('clusters',dd)

  ### Find Number of individuals
  N.Indiv=dd[min(Indiv)]
  N.Indiv=as.numeric(unlist(strsplit(N.Indiv, split=' ')))
  N.Indiv=N.Indiv[is.na(N.Indiv)!=T]

  ### Find Number of Clusters
  N.Clust=dd[min(Clust)]
  N.Clust=as.numeric(unlist(strsplit(N.Clust, split=' ')))
  N.Clust=N.Clust[is.na(N.Clust)!=T]

  ### Find Estimated Ln
  EstLn=grep('Estimated Ln Prob of Data',dd,value=T)
  EstLn.val=as.numeric(unlist(strsplit(EstLn, split=' ')))
  EstLn.val=EstLn.val[is.na(EstLn.val)!=T]

  ### Find individuals clustering
  ## Define Limits of table
  test.pos=sort(c(Indiv,Clust))
  test.chk=c(0,diff(test.pos))

  if(any(test.chk==1)){
    START=1+test.pos[test.chk==1]
    }

  STOP=START+N.Indiv-1

  ## Build Table
  table.vect=dd[START:STOP]
  table.mat=data.frame('list'=table.vect)

  ## Produce first line
    line.split=unlist(strsplit(as.vector(table.mat[2,1]),split=' '))
    line.probs=as.numeric(line.split[is.na(as.numeric(line.split))!=T])
    line.labels=line.split[is.na(as.numeric(line.split))==T]
    line.labels=line.labels[line.labels!='']
    line.labels=line.labels[line.labels!=':']

    line.tot=c(line.labels,line.probs)

  ## colnames
    l.labs=rep('lab',length(line.labels))
    l.probs=rep('P',length(line.probs))


  table.final=as.data.frame(matrix(0,STOP-START,length(line.tot)))
  for(i in 1:nrow(table.mat)){
    line.split=unlist(strsplit(as.vector(table.mat[i,1]),split=' '))

    line.probs=as.numeric(line.split[is.na(as.numeric(line.split))!=T])
    line.labels=line.split[is.na(as.numeric(line.split))==T]
    line.labels=line.labels[line.labels!='']
    line.labels=line.labels[line.labels!=':']

    table.final[i,1:length(line.labels)]=line.labels
    table.final[i,(length(line.labels)+1):ncol(table.final)]=line.probs
    }
table.final=cbind(rep(EstLn.val,nrow(table.final)),table.final)
colnames(table.final)=c('EstLn',l.labs,l.probs)
table.final
}

IMPORTMULTI.I=function(listfiles,path=getwd(),best=T){
  listfiles=sort(listfiles)

  if(length(listfiles)==1){
    f.best=1
    } else {
      ##### Explore Dataset
      expl=NULL
      for(i in 1:length(listfiles)){
        ### Open Structure Outputs
        dd=readLines(listfiles[i])
        ### Main
        Clust=grep('clusters',dd)
        ### Find Number of Clusters
        N.Clust=dd[min(Clust)]
        N.Clust=as.numeric(unlist(strsplit(N.Clust, split=' ')))
        N.Clust=N.Clust[is.na(N.Clust)!=T]
        ### Find Estimated Ln
        EstLn=grep('Estimated Ln Prob of Data',dd,value=T)
        EstLn.val=as.numeric(unlist(strsplit(EstLn, split=' ')))
        EstLn.val=EstLn.val[is.na(EstLn.val)!=T]
        st=cbind(i,N.Clust,EstLn.val)

        expl=rbind(expl,st)
        }
      stats=data.frame(K=1:max(expl[,2]),
		       ML=aggregate(expl[,3], by=list(expl[,2]), max)[,2], 
		       sd=aggregate(expl[,3], by=list(expl[,2]), sd)[,2])
      write.table(stats, file=paste(path,"summary_ML.txt",sep='/'), sep='\t')
      k.best=aggregate(expl[,3],by=list(expl[,2]),max)
      k.best=k.best[k.best[,1]!=1,]
      GetIt=function(x,...) expl[expl[,3]==x,1]
      f.best=unlist(sapply(k.best$x,GetIt))
      }

  ## Import and save f.bests
  dir.create(paste(path,'EXPORTS.IND',sep='/'), showWarnings = F, recursive = FALSE)
  if(best==T){
    for(i in 1:length(f.best)){
      dd=IMPSTRUCT.I(listfiles[f.best[i]])
      write.table(dd,paste(path,'EXPORTS.IND',paste(basename(listfiles[f.best[i]]),'.txt',sep=''),sep='/'),sep='\t')
      }
    } else {
      for(i in 1:length(listfiles)){
      dd=IMPSTRUCT.I(listfiles[i])
      write.table(dd,paste(path,'EXPORTS.IND',paste(basename(listfiles[i]),'.txt',sep=''),sep='/'),sep='\t')
      }
    }
  }



############## Importation of populations
#### Single Run importation:
# listfiles=choose.files(default = "", caption = "Select files",multi = TRUE, index = nrow(Filters))
# IMPSTRUCT.P(listfiles)

#### Multiple Runs importation:
# IMPORTMULTI.P()

IMPSTRUCT.P=function(output){
  ### Open Population Names file
  target=substr(basename(output),start=1,stop=3)
  allfiles=dir(dirname(output))
  txtfiles=allfiles[which(regexpr('.pop',allfiles)>0)]
  pops=0
  if(length(txtfiles)!=0){
    pops.f=which(regexpr(target,txtfiles)>0)
    popsinfos=paste(dirname(output),txtfiles[pops.f],sep='/')
    pops=read.table(popsinfos)
    }

  ### Open Structure Outputs
  dd=readLines(output)

  ### Main
  Clust=grep('clusters',dd)

  ### Find Number of Clusters
  N.Clust=dd[min(Clust)]
  N.Clust=as.numeric(unlist(strsplit(N.Clust, split=' ')))
  N.Clust=N.Clust[is.na(N.Clust)!=T]

  ### Find Estimated Ln
  EstLn=grep('Estimated Ln Prob of Data',dd,value=T)
  EstLn.val=as.numeric(unlist(strsplit(EstLn, split=' ')))
  EstLn.val=EstLn.val[is.na(EstLn.val)!=T]

  ### Find Populations clustering
  ## Define Limits of table
  START=grep('Given',dd)+3
  STOP=grep('Allele-freq.',dd)-3

  ## Build Table
  table.vect=dd[START:STOP]
  table.mat=data.frame('list'=table.vect)

  ## Produce first line
    line.split=unlist(strsplit(as.vector(table.mat[2,1]),split=' '))
    line.probs=as.numeric(line.split[is.na(as.numeric(line.split))!=T])
    line.labels=line.split[is.na(as.numeric(line.split))==T]
    line.labels=line.labels[line.labels!='']
    line.labels=line.labels[line.labels!=':']

    line.tot=c(line.labels,line.probs)

  ## colnames
    l.labs=rep('lab',length(line.labels))
    l.probs=rep('P',length(line.probs)-1)
    l.inds=c('NbIndivs')


  table.final=as.data.frame(matrix(0,STOP-START,length(line.tot)))
  for(i in 1:nrow(table.mat)){
    line.split=unlist(strsplit(as.vector(table.mat[i,1]),split=' '))

    line.probs=as.numeric(line.split[is.na(as.numeric(line.split))!=T])
    line.labels=line.split[is.na(as.numeric(line.split))==T]
    line.labels=line.labels[line.labels!='']
    line.labels=line.labels[line.labels!=':']

    table.final[i,1:length(line.labels)]=line.labels
    table.final[i,(length(line.labels)+1):ncol(table.final)]=line.probs
    }
if(pops==0){
  pops=sub(":","",table.final$V1)
  table.final=cbind(pops,rep(EstLn.val,nrow(table.final)),table.final)
  colnames(table.final)=c('CodePop','EstLn',l.labs,l.probs,l.inds)
  } else {
  table.final=cbind(pops,rep(EstLn.val,nrow(table.final)),table.final)
  colnames(table.final)=c('CodePop','EstLn',l.labs,l.probs,l.inds)
  }
  table.final
}

IMPORTMULTI.P=function(listfiles,path=getwd(),best=T){
  listfiles=sort(listfiles)
  if(length(listfiles)==1){
    f.best=1
    } else {
    ##### Explore Dataset
    expl=NULL
    for(i in 1:length(listfiles)){
      ### Open Structure Outputs
      dd=readLines(listfiles[i])
      ### Main
      Clust=grep('clusters',dd)
      ### Find Number of Clusters
      N.Clust=dd[min(Clust)]
      N.Clust=as.numeric(unlist(strsplit(N.Clust, split=' ')))
      N.Clust=N.Clust[is.na(N.Clust)!=T]
      ### Find Estimated Ln
      EstLn=grep('Estimated Ln Prob of Data',dd,value=T)
      EstLn.val=as.numeric(unlist(strsplit(EstLn, split=' ')))
      EstLn.val=EstLn.val[is.na(EstLn.val)!=T]
      st=cbind(i,N.Clust,EstLn.val)

      expl=rbind(expl,st)
      }
    stats=data.frame(K=1:max(expl[,2]),
		      ML=aggregate(expl[,3], by=list(expl[,2]), max)[,2], 
		      sd=aggregate(expl[,3], by=list(expl[,2]), sd)[,2])
    write.table(stats, file=paste(path,"summary_ML.txt",sep='/'), sep='\t')
    k.best=aggregate(expl[,3],by=list(expl[,2]),max)
    k.best=k.best[k.best[,1]!=1,]
    GetIt=function(x,...) expl[expl[,3]==x,1]
    f.best=unlist(sapply(k.best$x,GetIt))
    }

  ## Import and save f.bests
  dir.create(paste(path,'EXPORTS.POP',sep='/'), showWarnings = F, recursive = FALSE)
  if(best==T){
    for(i in 1:length(f.best)){
      dd=IMPSTRUCT.P(listfiles[f.best[i]])
      write.table(dd,paste(path,'EXPORTS.POP',paste(basename(listfiles[f.best[i]]),'.txt',sep=''),sep='/'),sep='\t')
      }
    } else {
      for(i in 1:length(listfiles)){
      dd=IMPSTRUCT.P(listfiles[i])
      write.table(dd,paste(path,'EXPORTS.POP',paste(basename(listfiles[i]),'.txt',sep=''),sep='/'),sep='\t')
      }
    }
}


############## Graphical User Interface
numbfiles=0
exp.dir=getwd()

"HANDLE.GUI" <-
function(...){

    options.apply=function(...){
      path=tclvalue(exp.dir)

      if(tclvalue(each.dum)==1) best.id=F
      if(tclvalue(best.dum)==1) best.id=T

      if(tclvalue(ind.dum)==1) IMPORTMULTI.I(listfiles,path,best.id)
      if(tclvalue(pop.dum)==1) IMPORTMULTI.P(listfiles,path,best.id)

      tkdestroy(OPTIONS)
      cat(paste('Files Successfully Produced in Exporting Directory'))
      tkdestroy(OPTIONS)
      }

    choose.apply=function(...){
      assign("exp.dir",tkchooseDirectory(),env=.GlobalEnv)
      tkdestroy(OPTIONS)
      HANDLE.GUI()
      }

    choose.f=function(...){
      listfiles=choose.files()
      assign("numbfiles",length(listfiles),env=.GlobalEnv)
      assign("listfiles",listfiles,env=.GlobalEnv)
      tkdestroy(OPTIONS)
      HANDLE.GUI()
      }

    quit.apply=function(...){
      tkdestroy(OPTIONS)
      }

    ## GUI
    OPTIONS <- tktoplevel()
    tktitle(OPTIONS) <- "Handle STRUCTURE outputs"

    ind.dum=tclVar(init="2")
    pop.dum=tclVar(init="2")

    each.dum=tclVar(init="2")
    best.dum=tclVar(init="2")

    frame2 <- tkframe(OPTIONS, relief = "groove", borderwidth = 2)
    browse.file <- tkbutton(frame2, text = "Select Files", command = choose.f,padx = 20, font = "arial 10")
    button.ind <- tkcheckbutton(frame2, text="Individual Probabilities", 'onvalue'=1, 'offvalue'=0, 'variable'=ind.dum)
    button.pop <- tkcheckbutton(frame2, text="Population Probabilities", 'onvalue'=1, 'offvalue'=0, 'variable'=pop.dum)
    nfiles <- tklabel(frame2,text=paste(numbfiles,'Files Imported',sep=' '))
    tkgrid(browse.file)
    tkgrid(button.ind)
    tkgrid(button.pop)
    tkgrid(nfiles, sticky = "w")
    tkpack(frame2)

    button.each <- tkcheckbutton(frame2, text="Export Each Treated Run", 'onvalue'=1, 'offvalue'=0, 'variable'=each.dum)
    button.best <- tkcheckbutton(frame2, text="Export Most Likely Run (one export per K value)", 'onvalue'=1, 'offvalue'=0, 'variable'=best.dum)
    #exp.gui <- tklabel(frame2,text=paste('Exporting Directory:',tclvalue(exp.dir),sep=' '))
    browse.but <- tkbutton(frame2, text = "Select Export Directory", command = choose.apply,padx = 10, font = "arial 10")
    tkgrid(browse.but)
    tkgrid(button.each,sticky = "w")
    tkgrid(button.best,sticky = "w")
    tkpack(frame2)

    frame3 <- tkframe(OPTIONS, relief = "groove", borderwidth = 1)
    save.but <- tkbutton(frame3, text = "Apply", command = options.apply,padx = 10, font = "arial 10")
    quit.but <- tkbutton(frame3, text = "Quit", command = quit.apply,padx = 10, font = "arial 10")
    tkgrid(save.but,quit.but,sticky = "w")
    #tkgrid(exp.gui)
    tkpack(frame3,side='bottom')
    }


############## Running the functions collection
SIMIL.GUI=function(...){
  winMenuAdd("SIMIL")
  winMenuAddItem('SIMIL','Handle STRUCTURE Outputs', 'HANDLE.GUI()')
  winMenuAddItem('SIMIL','Compute GSS', 'SIMIL()')
  }
  
SIMIL.GUI()
