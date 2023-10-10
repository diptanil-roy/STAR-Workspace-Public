#!/bin/bash

mkdir -p /gpfs01/star/pwg_tasks/jetcorr03/zweighted_sample

for (( bin = 0; bin < 10; ++bin ));
do
	
zlow=`echo "scale=1; $bin/10" | bc`
zhigh=`echo "scale=1; $(( bin + 1 ))/10" | bc`

outputdir="/gpfs01/star/pwg_tasks/jetcorr03/zweighted_sample/ZBin_$bin"
mkdir -p $outputdir
rm -rf $outputdir/*
mkdir -p $outputdir/out
mkdir -p $outputdir/log
mkdir -p $outputdir/gen
mkdir -p $outputdir/trash


cat << EOT > SubmissionScript_$bin.xml
<?xml version="1.0" encoding="utf-8" ?> 
<!DOCTYPE note [
<!ENTITY PRODNAME   "pythia8.200GeV.pp.HFjets.ZBin.$bin">
<!ENTITY RECOCHAIN  "P16idAuAu200hft">

<!ENTITY NEVENTS    "200">
<!ENTITY STRIDE     "100">

<!ENTITY NFILES     "1000">
<!ENTITY WORKINGDIR "/gpfs01/star/pwg/droy1/STAR-Workspace/D0Analysis/SimulationAuAu/">
<!ENTITY DATASET    "filelist:&WORKINGDIR;/simulation.y2014x.pythia8.HFjets.input.3_inf.ZBin_0$bin.list" >

<!ENTITY MUTOPICODIR "/gpfs01/star/pwg/droy1/STAR-Workspace/D0Analysis/MuDstToPicoDst/">


<!ENTITY OUTDIR     "$outputdir/out/" >
<!ENTITY LOGDIR     "$outputdir/log/" >
<!ENTITY JOBFILES   "$outputdir/gen/" > 

<!ENTITY GEOMETRY      "y2014x"> 
<!ENTITY TIMESTAMP     "sdt20140530"> 
<!ENTITY DBVTIMESTAMP  "DbV20150316"> 

<!ENTITY RECLIB     "SL16d_embed">
<!ENTITY SIMLIB     "dev">
<!ENTITY SIMULATE      "false">

]>

<job  name               = "&PRODNAME;" 
      maxFilesPerProcess = "1"
      fileListSyntax     = "paths"  
      filesPerHour       = "0.1"
      simulateSubmission = "&SIMULATE;">

  <!--logControl="UCM"-->
  <stdout   URL="file:&JOBFILES;\$JOBID_\$FILEBASENAME.log" />
  <stderr   URL="file:&JOBFILES;\$JOBID_\$FILEBASENAME.err" />
  <input    URL="&DATASET;" nFiles="&NFILES;" />

  <Generator><Location>&JOBFILES;</Location></Generator>
  <ResourceUsage><Priority>100</Priority></ResourceUsage>

  <command>       

    echo "JobStartTime="\`/bin/date\`
#    setenv datevalue \`cat \${SUBMITTINGDIRECTORY}/date.txt\`    

 
    #########################################################
    #  Check that the SRM script can be found      
    #########################################################
    
       source /star/u/starreco/.cshrc
       unsetenv NODEBUG
    #setenv DB_SERVER_LOCAL_CONFIG /afs/rhic.bnl.gov/star/packages/conf/dbLoadBalancerLocalConfig_nightly.xml

    setup 32b
    starver &SIMLIB;    

    set production = &PRODNAME;
    set site       = rcf

    echo "SUMS_prodTag=\${production}"
    echo "SUMS_site=\${site}"
    echo "SUMS_prodType=reco"

    cd \$SCRATCH

    #########################################################
    # Transfer files needed to run the simulation
    #########################################################
    mkdir StRoot/StarGenerator/macros/ -p
    cd StRoot/StarGenerator/macros/ 
    cp -R &WORKINGDIR;/StRoot/StarGenerator/macros/* .
    cd -
    cp -R &WORKINGDIR;/starsim_zbin.C .
    cp -R &WORKINGDIR;/runBfc.C .
    cp -r /star/u/droy1/Y2019/STAR/FastJet .

	ln -s /star/u/droy1/Y2019/STAR/FastJet/fastjet-install/include/fastjet
	ln -s /star/u/droy1/Y2019/STAR/FastJet/fastjet-install/include/siscone

	setenv FASTJET FastJet/fastjet-install

    cp -R &WORKINGDIR;/.\$STAR_HOST_SYS .

    # Check to see that everything is setup and making sense
    ls -la

    rehash

cat &lt;&lt;+++
    ######################################################### 
    # Run starsim \`date\` 
    #########################################################
+++
    which starsim   

    set runfile   = \$INPUTFILE0
#   set runnumber = \`basename \$runfile:t .txt\`
    set runnumber = \`echo \$runfile:t | awk -F "_" '// { print \$1 }'\`
    set chunk     = \`echo \$runfile:t | awk -F "_" '// { print \$2 }'\`

    echo RUNNUMBER = \${runnumber}

    set   EMBEDDING_START_TIME = \`date\`
    # Run starsim and create 10 output fzd files

    cp \${runfile} .

    set finalrunfile = \`echo \$runfile | awk -F "/" '// { print \$NF }'\`

    echo \${finalrunfile} 
    echo \${finalrunfile} &gt; mc_\${runnumber}_\${chunk}.log
    echo SIMULATION_START_TIME = \`date\` &gt;&gt; mc_\${runnumber}_\${chunk}.log
    root4star -q -b starsim_zbin.C\(&NEVENTS;,\${runnumber},\"\${finalrunfile}\",1,&STRIDE;,3,0,$zlow,$zhigh\) &gt;&gt; mc_\${runnumber}_\${chunk}.log
    echo SIMULATION_END_TIME = \`date\` &gt;&gt; mc_\${runnumber}_\${chunk}.log

    set fzdin=\`ls *.fzd\`

    echo \${fzdin}

    cp *.fzd &OUTDIR;

    # Run reconstruction
    starver &RECLIB;
      mv .\$STAR_HOST_SYS meh

    which starsim 

    echo RECONSTRUCTION_START_TIME = \`date\` &gt; rc_\${runnumber}_\${chunk}.log
    root4star -q -b runBfc.C\(&NEVENTS;,\"\${fzdin}\"\)             &gt;&gt; rc_\${runnumber}_\${chunk}.log
    echo RECONSTRUCTION_END_TIME = \`date\` &gt;&gt; rc_\${runnumber}_\${chunk}.log

    cp *.MuDst.root &OUTDIR;

    set   EMBEDDING_END_TIME = \`date\`

    echo EMBEDDING_START_TIME = \${EMBEDDING_START_TIME} &gt; stats_\${runnumber}_\${chunk}.log
    starver &SIMLIB;
    cp -R &MUTOPICODIR;/.\$STAR_HOST_SYS .
    cp -R &MUTOPICODIR;/genDst.C .

    set mudstfile=\`ls *.MuDst.root\`

    root4star -l -b -q genDst.C\(-1,\"picoDst,PicoVtxMode:PicoVtxVpdOrDefault,PicoCovMtxMode:PicoCovMtxWrite\",\"\${mudstfile}\"\) &gt;&gt; stats_\${runnumber}_\${chunk}.log

    echo EMBEDDING_END_TIME = \${EMBEDDING_END_TIME} &gt;&gt; stats_\${runnumber}_\${chunk}.log
    grep SIMU*_TIME mc*_\${chunk}.log &gt;&gt; stats_\${runnumber}_\${chunk}.log 
    grep RECO*_TIME rc*_\${chunk}.log &gt;&gt; stats_\${runnumber}_\${chunk}.log
    ls -l --human-readable *.root *.fzd &gt;&gt; stats_\${runnumber}_\${chunk}.log

    cp *.picoDst*.root &OUTDIR;

    if ( &NEVENTS; &lt; 11 ) then
       cp *.fzd &OUTDIR;
       cp *.log*        &LOGDIR;   
     else
       gzip --best *.log
       cp *.log*        &LOGDIR;
     endif

   </command>


</job>

EOT

star-submit SubmissionScript_$bin.xml

mv SubmissionScript_$bin.xml $outputdir
mv *.dataset $outputdir/trash/.
mv *.session.xml $outputdir/trash/.

done
