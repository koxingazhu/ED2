#==========================================================================================#
#==========================================================================================#
#     Leave these commands at the beginning.  They will refresh the session.               #
#------------------------------------------------------------------------------------------#
rm(list=ls())
graphics.off()
#==========================================================================================#
#==========================================================================================#



#==========================================================================================#
#==========================================================================================#
#      Here is the user defined variable section.                                          #
#------------------------------------------------------------------------------------------#

#----- Paths. -----------------------------------------------------------------------------#
# here           = "thispath"    # Current directory.
# there          = "thatpath"    # Directory where analyses/history are 
# srcdir         = "thisrscpath" # Source  directory.
# outroot        = "thisoutroot" # Directory for figures
# pov.incs       = "thispovincs" # Path with POV-Ray include files

# 
here           = "F:\\@@@ED2\\@3D Tree community"    # Current directory.
there          = "F:\\@@@ED2\\@3D Tree community\\h5output"    # Directory where analyses/history are 
srcdir         = "F:\\@@@ED2\\@3D Tree community\\R-utils" # Source directory.
outroot        = "F:\\@@@ED2\\@3D Tree community\\figures" # Directory for figures
# pov.incs       = "F:\\3D_Tree_community\\povray-3.7-stable\\povray-3.7-stable\\distribution\\include\\" # Path with POV-Ray include files
pov.incs       = "F:\\@@@ED2\\@3D Tree community\\povray-3.7-stable\\povray-3.7-stable\\distribution\\include" # Path with POV-Ray include files


# here           = "/data/users/zhuyu/3D_tree/@3D_Tree_community"    # Current directory.
# there          = "/data/users/zhuyu/3D_tree/@3D_Tree_community/h5output"    # Directory where analyses/history are 
# srcdir         = "/data/users/zhuyu/3D_tree/@3D_Tree_community/R-utils" # Source directory.
# outroot        = "/data/users/zhuyu/3D_tree/@3D_Tree_community/figures" # Directory for figures
# # pov.incs       = "F:\\3D_Tree_community\\povray-3.7-stable\\povray-3.7-stable\\distribution\\include\\" # Path with POV-Ray include files
# pov.incs       = "/usr/local/share/povray-3.6/include" # Path with POV-Ray include files



#------------------------------------------------------------------------------------------#


#----- Time options. ----------------------------------------------------------------------#
monthbeg    = 1   # First month to use
yearbeg     = 1999    # First year to consider
yearend     = 1999    # Maximum year to consider
reload.data = FALSE         # Should I reload partially loaded data?
pov.month   = 5            # Months for POV-Ray plots
pop.scale   = 1.0          # Scaling factor to REDUCE displayed population.
sasmonth    = sequence(12)
#------------------------------------------------------------------------------------------#



#----- Name of the simulations. -----------------------------------------------------------#
myplaces       = c("harvard")
#------------------------------------------------------------------------------------------#



#----- Plot options. ----------------------------------------------------------------------#
depth          = 1200                   # PNG resolution, in pixels per inch
depth          = 6000
paper          = "letter"               # Paper size, to define the plot shape
ptsz           = 16                     # Font size.
ibackground    = 0           # Background settings (check load_everything.r)
ibackground    = 1
#------------------------------------------------------------------------------------------#


#------ Miscellaneous settings. -----------------------------------------------------------#
slz.min        = -5.0         # The deepest depth that trees access water.
idbh.type      = 1   # Type of DBH class
                              # 1 -- Every 10 cm until 100cm; > 100cm
                              # 2 -- 0-10; 10-20; 20-35; 35-50; 50-70; > 70 (cm)
                              # 3 -- 0-10; 10-35; 35-70; > 70 (cm)
myklight=1
myallom=0

klight         = myklight     # Weighting factor for maximum carbon balance
iallom         = myallom      # Allometry

# isoil.hydro    = myslhydro    # Soil hydrology method
#------------------------------------------------------------------------------------------#

pov.dbh.min=10
pov.gap.line=20
pov.n.line=20
pov.gap.area = pov.gap.line^2
pov.nxy.patch = pov.n.line^2
pov.total.area = pov.gap.area * pov.nxy.patch



#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#      NO NEED TO CHANGE ANYTHING BEYOND THIS POINT UNLESS YOU ARE DEVELOPING THE CODE...  #
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#
#==========================================================================================#



#----- Loading some packages and scripts. -------------------------------------------------#
# source(file.path(srcdir,"load.everything.r"))
source(file.path(srcdir,"locations.r"))
source(file.path(srcdir,"plotsize.r"))
source(file.path(srcdir,"timeutils.r"))
source(file.path(srcdir,"monthly.template.r"))

# Are you sure you want to proceed [y|N]? N
# Error in eval(ei, envir) : Missing packages!!!

#------------------------------------------------------------------------------------------#


#----- Avoid unecessary and extremely annoying beeps. -------------------------------------#
options(locatorBell=FALSE)
#------------------------------------------------------------------------------------------#


# 
# #----- Load observations. -----------------------------------------------------------------#
# obsrfile = paste(srcdir,"LBA_MIP.v9.RData",sep="/")
# load(file=obsrfile)
# #------------------------------------------------------------------------------------------#



#----- Define plot window size ------------------------------------------------------------#
size = plotsize(proje=FALSE,paper=paper)
#------------------------------------------------------------------------------------------#



#---- Create the main output directory in case there is none. -----------------------------#
if (! file.exists(outroot)) dir.create(outroot)
#------------------------------------------------------------------------------------------#


place="harvard"
where = place
ici = tolower(where)

   #----- Retrieve default information about this place and set up some variables. --------#
   thispoi = locations(where,here=there,yearbeg=1999,yearend=1999,monthbeg=5,daybeg=1
                       ,filetype="E",fullonly=FALSE) # 
 
   inpref  = thispoi$pathin
   bnpref  = basename(inpref)
   outmain = paste(outroot,place,sep="/")
   outpref = paste(outmain,"povray",sep="/")
   lieu    = thispoi$lieu
   iata    = thispoi$iata
   suffix  = thispoi$iata
   yeara   = thispoi$yeara
   yearz   = thispoi$yearz
   meszz   = thispoi$monz
   #---------------------------------------------------------------------------------------#



   #----- Create the directories in case they don't exist. --------------------------------#
   if (! file.exists(outmain)) dir.create(outmain)
   if (! file.exists(outpref)) dir.create(outpref)
   #---------------------------------------------------------------------------------------#




   #----- Print a banner to entretain the user. -------------------------------------------#
   cat(" + Post-processing output from ",lieu,"...","\n")
   #---------------------------------------------------------------------------------------#


   
   #---------------------------------------------------------------------------------------#
   #  PATCH -- patch level variables, we save as lists because the dimensions vary.    #
   #---------------------------------------------------------------------------------------#
   patch                = list()
   patch$ipa            = list()
   patch$area           = list()
   
   NPATCHES_GLOBAL <- h5read("F:\\@@@ED2\\@3D Tree community\\h5output\\harvard\\analy\\harvard\\analy-E-1999-05-00-000000-g01.h5", "NPATCHES_GLOBAL")
   patch$ipa = sequence(NPATCHES_GLOBAL)
   
   AREA     = h5read("F:\\@@@ED2\\@3D Tree community\\h5output\\harvard\\analy\\harvard\\analy-E-1999-05-00-000000-g01.h5", "AREA")
   areasi     = h5read("F:\\@@@ED2\\@3D Tree community\\h5output\\harvard\\analy\\harvard\\analy-E-1999-05-00-000000-g01.h5", "AREA_SI")
   npatches   = h5read("F:\\@@@ED2\\@3D Tree community\\h5output\\harvard\\analy\\harvard\\analy-E-1999-05-00-000000-g01.h5", "SIPA_N")
   areapa      = AREA * rep(areasi,times=npatches)
   patch$area = areapa
   #---------------------------------------------------------------------------------------#
   
   #----- Cohort level, we save as lists because the dimensions vary. ---------------------#
   cohort                = list()
   cohort$ipa            = list()
   cohort$ico            = list()
   cohort$pft            = list()
   cohort$dbh            = list()

   PACO_N <- h5read("F:\\@@@ED2\\@3D Tree community\\h5output\\harvard\\analy\\harvard\\analy-E-1999-05-00-000000-g01.h5", "PACO_N")
   ncohorts    = PACO_N
   ipaconow    = rep(sequence(NPATCHES_GLOBAL),times=PACO_N)
   cohort$ipa  = ipaconow
   icoconow    = unlist(sapply(X = PACO_N, FUN = sequence))
   cohort$ico  = icoconow
   cohort$pft <- h5read("F:\\@@@ED2\\@3D Tree community\\h5output\\harvard\\analy\\harvard\\analy-E-1999-05-00-000000-g01.h5", "PFT")
   cohort$dbh <- h5read("F:\\@@@ED2\\@3D Tree community\\h5output\\harvard\\analy\\harvard\\analy-E-1999-05-00-000000-g01.h5", "DBH")

   NPLANT <- h5read("F:\\@@@ED2\\@3D Tree community\\h5output\\harvard\\analy\\harvard\\analy-E-1999-05-00-000000-g01.h5", "NPLANT")
   HITE <- h5read("F:\\@@@ED2\\@3D Tree community\\h5output\\harvard\\analy\\harvard\\analy-E-1999-05-00-000000-g01.h5", "HITE")
   
   pftconow = cohort$pft
   nplantconow = NPLANT
   areaconow  = rep(areapa,times=ncohorts)
   cohort$nplant = nplantconow * areaconow
   
   #---------------------------------------------------------------------------------------#

# 
#    years   = sort(unique(datum$year))
    years=1
#    pclabs  = paste("y",sprintf("%4.4i",years),"m",sprintf("%2.2i",pov.month),sep="")
#    pcwhens = paste(mon2mmm(pov.month,cap1=TRUE),years,sep="-")
    pcwhen = paste(mon2mmm(pov.month,cap1=TRUE),years,sep="-")
    pcout  = paste(sprintf("%4.4i",years),sprintf("%2.2i",pov.month),sep="-")
#    pcloop  = sequence(length(pclabs))


      w = 1
      #----- Grab this time. --------------------------------------------------------------#

      #------------------------------------------------------------------------------------#
      cat (" + Creating POV-Ray image for ",pcwhen,"...","\n")


      #----- Copy some patch variables to local variables. --------------------------------#
      ipa      = patch$ipa
      areapa   = patch$area
      #------------------------------------------------------------------------------------#



      #----- Copy some cohort variables to local variables. -------------------------------#
      ipaco    = cohort$ipa
      icoco    = cohort$ico
      nplantco = round(cohort$nplant * pov.total.area * pop.scale)
      dbhco    = cohort$dbh   
      pftco    = cohort$pft 
      #------------------------------------------------------------------------------------#


      #----- Remove small plants to reduce the clutter. -----------------------------------#
      keep     = is.finite(dbhco) & dbhco >= pov.dbh.min
      ipaco    = ipaco   [keep]
      icoco    = icoco   [keep]
      nplantco = nplantco[keep]
      dbhco    = dbhco   [keep]
      pftco    = pftco   [keep]
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #     Find the tree coordinates.                                                     #
      #------------------------------------------------------------------------------------#
         #----- Determine the quadricules where trees of each patch can go. ---------------#
         npatches    = length(areapa)
         nquadpa.1st = round( pov.nxy.patch * areapa )
         totquad     = sum(nquadpa.1st)
         off         = pov.nxy.patch - totquad
         nfix        = sample( c( rep( x = 0        , times = npatches-abs(off))
                                , rep( x = sign(off), times = abs(off)         )
                                )#end c
                             )#end sample
         nquadpa     = nquadpa.1st + nfix


         #---------------------------------------------------------------------------------#
         #      Create a list with the quadricules for each patch.  We double each count   #
         # to make sure each patch gets at least 2 numbers (so we can safely use 'sample'. #
         #---------------------------------------------------------------------------------#
         ipa.quad    = unlist(mapply(FUN=rep,x=ipa,each=2*nquadpa))
         quad        = split(x=rep(sample(pov.nxy.patch),each=2),f=ipa.quad)
         names(quad) = paste("patch",sprintf("%3.3i",sort(unique(ipa))),sep="_")
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Expand population variables.                                                #
         #---------------------------------------------------------------------------------#
         ipaco  = unlist(mapply(FUN=rep,x=ipaco,each=nplantco))
         dbhco  = unlist(mapply(FUN=rep,x=dbhco,each=nplantco))
         pftco  = unlist(mapply(FUN=rep,x=pftco,each=nplantco))
         nco    = sum(nplantco)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Randomly assign the cohorts to the patches where they belong.               #
         #---------------------------------------------------------------------------------#
         nco.pa.pop  = table(ipaco)
         lab         = paste("patch",sprintf("%3.3i",as.integer(names(nco.pa.pop))),sep="_")
         idx         = match(lab,names(quad))
         nco.pa      = rep(0,times=length(quad))
         nco.pa[idx] = nco.pa.pop
         quadco      = unlist( mapply( FUN      = sample
                                     , x        = quad
                                     , size     = nco.pa
                                     , MoreArgs = list(replace=TRUE)
                                     )#end mapply
                             )#end unlist
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #     Assign coordinates for all cohorts.                                         #
         #---------------------------------------------------------------------------------#
         pov.patch.xmax <- 380
         pov.patch.ymax <- 380
         
         # xco    = ( sample(x=0.1*(sequence(10*pov.patch.xmax)-0.5),size=nco,replace=TRUE)
         #          + pov.x0[quadco] )
         # yco    = ( sample(x=0.1*(sequence(10*pov.patch.ymax)-0.5),size=nco,replace=TRUE)
         #          + pov.y0[quadco] )
         
         
         
         xco    = ( sample(x=0.1*(sequence(10*pov.patch.xmax)-0.5),size=nco,replace=TRUE)-206
                    )
         yco    = ( sample(x=0.1*(sequence(10*pov.patch.ymax)-0.5),size=nco,replace=TRUE)-206
                    )
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #     Make the labels.                                                            #
         #---------------------------------------------------------------------------------#
         povplant = rbind(       "//----- The plants. ---------------------------------//"
                         , rbind( paste("plant(", unlist( mapply( FUN = paste
                                                                , sprintf("%2i"  ,iallom)
                                                                , sprintf("%7.2f",dbhco )
                                                                , sprintf("%2i"  ,pftco )
                                                                , sprintf("%7.2f",xco   )
                                                                , sprintf("%7.2f",yco   )
                                                                , MoreArgs = list(sep=",")
                                                                )#end mapply
                                                        )#end unlist
                                       ,")"
                                       ,sep=""
                                       )#end paste
                                )#rbind
                         ,       "//---------------------------------------------------//"
                         )#end rbind
         #---------------------------------------------------------------------------------#
      #------------------------------------------------------------------------------------#


      #------------------------------------------------------------------------------------#
      #      Copy the POV-ray template file to the working directory and append the set-   #
      # tings for this time.                                                               #
      #------------------------------------------------------------------------------------#
      # povscript = file.path(here,place,"polygon.pov")
      povscript = file.path(here,place,"povray_template.pov")
      dummy     = file.copy(from=file.path(srcdir,"povray_template.pov"),to=povscript)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Decide the colour of the text depending on the background.                    #
      #------------------------------------------------------------------------------------#
      if (ibackground == 0){
         pigment = "       pigment  { color rgb <0. ,0. ,0. >}   "
      }else if (ibackground == 1){
         pigment = "       pigment  { color rgb <1. ,1. ,1. >}   "
      }else if (ibackground == 2){
         pigment = "       pigment  { color rgb <1. ,1. ,1. >}   "
      }#end if
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #      Append the title and time stamp.                                              #
      #------------------------------------------------------------------------------------#
      lesim    = unlist(strsplit(lieu,split="\n"))
      povtitle = NULL
      for (n in sequence(length(lesim))){
         yt       = sprintf("%7.1f",200 - 20 * (n-1))
         povtitle = rbind( povtitle
                         ,       "//----- The header. --------------------------------//"
                         , paste("text { ttf \"cyrvetic.ttf\" \"",lesim[n],"\" 5,0",sep="")
                         ,               pigment
                         ,       "       finish{ ambient 1.0 diffuse 0.0}"
                         ,       "       scale     <   12.0,   12.0,    0.1>"
                         ,       "       rotate    <   28.8,    0.0,    0.0>"
                         ,       "       rotate    <    0.0,   45.0,    0.0>"
                         , paste("       translate < -150.0,", yt,",  200.0>",sep="")
                         ,       "     }// end text"
                         ,       "//--------------------------------------------------//"
                         ,       " "
                         ,       " "
                         )#end rbind
      }#end for (n in 1:lesim)
      povstamp = rbind(       "//----- The time stamp. ----------------------------//"
                      , paste("text { ttf \"cyrvetic.ttf\" \"",pcwhen,"\" 5,0",sep="")
                      ,        pigment
                      ,       "       finish{ ambient 1.0 diffuse 0.0}"
                      ,       "       scale     <   12.0,   12.0,    0.1>"
                      ,       "       rotate    <   28.8,    0.0,    0.0>"
                      ,       "       rotate    <    0.0,   45.0,    0.0>"
                      ,       "       translate <  100.0,  180.0, -150.0>"
                      ,       "     }// end text"
                      ,       "//--------------------------------------------------//"
                      ,       " "
                      ,       " "
                      )#end rbind
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Append the title, stamp, and plants...                                        #
      #------------------------------------------------------------------------------------#
      write(x = povtitle, file = povscript, ncolumns = 1, append = TRUE)
      write(x = povstamp, file = povscript, ncolumns = 1, append = TRUE)
      write(x = povplant, file = povscript, ncolumns = 1, append = TRUE)
      #------------------------------------------------------------------------------------#

      # For Windows, using POV-Ray 3.7
      
      # 
      # #------------------------------------------------------------------------------------#
      # #     Call POV-Ray (must use system, though...)                                      #
      # #------------------------------------------------------------------------------------#
      # povray  = Sys.which("povray")
      # outtemp = file.path(tempdir(),paste(place,".png",sep=""))
      # outfile = file.path(outpref,paste(bnpref,"-",pcout,".png",sep=""))
      # 
      # povopts = paste("-D"
      #                ,"-V"
      #                ,"+UA"
      #                ,paste0("+L",pov.incs)
      #                ,paste0("+W",round(size$width*depth ))
      #                ,paste0("+H",round(size$height*depth))
      #                ,paste0("+O",outtemp)
      #                ,sep = " "
      #                )#end paste
      # dummy   = system( command       = paste(povray,povopts,povscript,sep=" ")
      #                 , intern        = TRUE
      #                 , ignore.stdout = TRUE
      #                 , ignore.stderr = TRUE
      #                 )#end system
      # dummy   = file.copy(from=outtemp,to=outfile,overwrite=TRUE)
      # dummy   = file.remove(outtemp,povscript)
      # #------------------------------------------------------------------------------------#

   #---------------------------------------------------------------------------------------#

#==========================================================================================#
#==========================================================================================#
