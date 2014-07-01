Raghav Kunnawalkam Elayavalli
June 4th 2014
CERN 

RAA: 
This directory is created with the wisdom obtained attending quark matter 2014. this comes after seeing people suffer during the build up by making tons of plots with corrections after corrections and with small changes etc... all this would have killed me if i had my previous set up. hence the stark reality and the need for change. The main idea with the analysis setup which i must follow is the following steps: 
    1) Reading macro
       this macro reads in the data files and creates useful histograms like the pt and eta and what ever etc... so that we can avoid rerunning over all the files again and again. saves us a lot of time. also create root files which contain those histograms. 
    2) Analysis macro 
       As the name suggests we will do the analysis steps here. like unfolding and calculating the RAA etc... the final result will be stored as raw text files of histograms which will help us to read them in later and also while sharing the output as txt files and also save them as root files. do both in the output directory.       
    3) Plotting macro
       this is the major change that ive noticed is that I should separate from the previous macros. people are extremely picky about the plots and so it would help a lot to have a plotting macro that plots things quick and fast. 

Another important point I realize as very important is that we need to have a standard setup for data files (txt files). just like the way HEP data is saved and transmitted. people know how to read that and it is also easier to move and to debug later. 

So here is the analysis setup.(in terms of the directory) 

                      Home - rkunnawa
                            WORK
                            RAA
                      CMSSW_5_3_8_HI/src/
               Macro        Headers          Output

like i mentioned before we will have the three macros followed by the root logon file. 
Please include the comments for each file (as you have been doing quite well so far) and good luck!!!  

June 12th
  Finished adding the read macros for data and mc. almost done with the unfolding/analysis macros. 

June 18th
  Ive added a macro fakeCheck.C which looks at the fake jets. it runs on cordor and you have to run it like so 
  $ source condor_submit.sh Nfiles nfilesperjob
  finished the analysis macro. Only have to check it with the corresponding plots which are freaking numerous to make.  
  I still need to add the plotting macro. 

July 1st 
  Ive added the plotting macro and the NLO comparison macro from the previous setup. need to add the macros from Gabie as soon as i get it from her. 
  finish all the rest of the macros needed for the analysis.
