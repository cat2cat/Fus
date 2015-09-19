#########################################################################
## Use of this script to run PRADA(RNAseq)                             ##
## Example: sh FlowControl6.sh -SARC.id 93                             ##
## arg1:SARC is cancertype(Main folder name)                           ##                                                         ##
## arg2:Number of jobs for PRADA-preprocess to be submited             ##
## arg3:Disk limit 85% for safe                                        ##
## arg4: download speed                                                ##
## arg5:config for PRADA version (walltime)                            ##
## arg6:every 30 mins to loop over all samples                         ##
#########################################################################

#!/usr/bin/bash
TODAY=$(date)
HOST=$(hostname)
echo "-----------------------------------------------------"
echo "Date: $TODAY  Host:$HOST"
echo "-----------------------------------------------------"

while true; do

## Argument followed command line
Nosubmit=${2:-15}    ## each loop
SpaceLimGT=95        ## If space exceed 95%,stop downloading

filename=`echo "$1"  | sed 's/-//'`        ## analysisID
Cancertype=`echo $filename | sed 's/.id//'| sed 's/.All//'`

DiskLimit=${3:-85}                         ## Scratch Space Percent Limit 
#speed=${4:-5}                             ## download speed
PBSadj=${5:-a}                             ## PRADA conf.txt files with different parameters
TimeInterval=${6:-15}                      ## Interval of sleep 15mins 
duration=${7:-7200}                        ## Repeat from the beginning after 120mins running
 

ScDir="/scratch/genomic_med/xhu3"
ReDir="/rsrch1/genomic_med/xhu3"
 
ScSPACE=$(df -h $ScDir | grep -oP "[0-9]{1,3}(?=%)")
ReSPACE=$(df -h $ReDir | grep -oP "[0-9]{1,3}(?=%)")
 

## My files for each cancer type 
Maindir="/scratch/genomic_med/xhu3/fusion/$Cancertype"

rsrDirfastq="/rsrch1/genomic_med/xhu3/fusion/$Cancertype.FastqEnd"
rsrDirtar="/rsrch1/genomic_med/xhu3/fusion/$Cancertype.tar"
rsrDirRNAseq="/rsrch1/genomic_med/xhu3/fusion/$Cancertype.RNAseq"

PRADAdir="/scratch/genomic_med/xhu3/PRADA/pyPRADA_1.2"

#################################################################################################
##  Download Samples,if scrach space exceed 95% limit, save the files in rsrch1 instead        ##
#################################################################################################

date1=$(date +"%s")
#echo $date1 

declare -i sample=0
COUNTER=1
cd $Maindir

 while read -r line; do
   ## Check each sample from the list of the cancer type
   folder=$line
   ((sample++))
   
   #Script to check if EOB has passed
   hour=`date +%k`
   if [ $hour -gt 7 ] && [ $hour -lt 18 ];then
      time="DAY"
   else
      time="NIGHT"
   fi
   echo -e "\nIt is about $hour o'clock,It is $time."
   
   if [ "$time" == "NIGHT" ];then
      child=8
   else
      child=5
   fi
   speed=${4:-$child}       ## download speed
   
    ## Moninor realtime PBS running status
    Rjobs=`qstat -u xhu3 -f | grep 'job_state = R'| wc -l`
    Qjobs=`qstat -u xhu3 -f | grep 'job_state = Q'| wc -l`
 
    BAMJobsByTime=`qstat -u xhu3 -f | grep 'resources_used.walltime'| awk '{gsub("resources_used.walltime = ", "");print}'|cut -d':' -f1 | awk '$1>=3 && $1<=25 {s+=1} END {print s+0}'`
 
    # echo -e "\n$Rjobs PBS jobs are running."
    # echo -e "$Qjobs PBS jobs are in queue."
    # echo -e "$BAMJobsByTime PBS jobs are using temp space.\n"

   echo -e "\nStart No#$sample Sample from $Cancertype list analysisID : $folder @`date`"

   ## First check whether files are available from rsrch1
   FilersrDLstatus="$rsrDirRNAseq/$folder.downloadStatus.txt"
   FilersrDL="$rsrDirRNAseq/$folder"
   PanCancerDL="/rsrch1/rists/cghub/TCGA/RNA/$Cancertype/$folder"

   cd $Maindir/$folder  

################################################################################ 
## Unzip and rename the files                                                ###
################################################################################

 tarfile=`find -name '*tar'`
 targzFile=`find -name '*tar.gz'`
 EndfastqNumber=`find -name '*end*.fastq' | wc -l`
 FastqNumber=`find *.fastq | wc -l`

 EndfastqNumberEnd=`find -name '*end*.fastq'`

 echo "Number of End Fastq = $EndfastqNumber"
     
 if [ $EndfastqNumber == 2 ]; then
     echo "Both Endfastq file Exist, their file size =";
     ls $Maindir/$folder/*.fastq | du -h ; 
     if [ -f "$tarfile" ] || [ -f "$targzFile" ]; then
         echo "Moving Tar file to rsrch1..."
         mv $Maindir/$folder/*tar*  $rsrDirtar
     else  
         echo "One tar file for $folder is already moved to rsrch1."
     fi  
 elif [ $FastqNumber == 1 ]; then
     ls $Maindir/$folder/*fastq >> $Maindir/BADtar.summary.txt
     mkdir $rsrDirfastq/$folder
     echo "Moving out single end fastq file, Don't stop...."
     mv $Maindir/$folder/*gz     $rsrDirfastq/$folder
     mv $Maindir/$folder/*fastq  $rsrDirfastq/$folder
     # awk '!/$folder/' $Maindir/$Cancertype.id > temp && mv temp $Maindir/$Cancertype.id
     continue
 else 
    echo "Need to Unzip and format original files."
    # if [ "$gtStatus" == "Downloaded" ] && [ -f "$targzFile" ] ;then
    if [ -f "$targzFile" ] ;then

         echo "Unzip tar.gz file for $folder,Don't interupt........"   
         tar -zxvf $Maindir/$folder/*tar*
         
         ## Rename the fastq file
         fastqNumber=`find -name '*fastq' | wc -l`
             echo $fastqNumber format Check for zipped fastq files
         if [ $fastqNumber == 2 ]; then
             echo "fastq file rename with end!!";
             NosepString=`ls *tar.gz | sed 's/[^_]//g'  | awk '{ print length }'`
             NosepStringFAS=`ls *1.fastq | sed 's/[^.]//g'  | awk '{ print length }'`
             if [ $NosepStringFAS == 2 ]; then
                fa=`ls *1.fastq | cut -d/ -f8 |cut -d. -f1 | sort -u`;
                fa1=`ls *1.fastq | cut -d/ -f8 |cut -d. -f1-2 | sort -u`;
                fa2=`ls *2.fastq | cut -d/ -f8 |cut -d. -f1-2 | sort -u`; 
                echo $fa1.fastq  $fa2.fastq $fa.1.fastq
                mv $Maindir/$folder/$fa1.fastq  $Maindir/$folder/$fa.end1.fastq
                mv $Maindir/$folder/$fa2.fastq  $Maindir/$folder/$fa.end2.fastq
                
             elif [ $NosepStringFAS == 1 ]; then

                fa1=`ls $Maindir/$folder/*1.fastq | cut -d/ -f8 |cut -d. -f1 | sort -u`;
                fa=`echo $fa1 | cut -d_ -f1-$NosepString`       ## Edit depends on sample center
                echo $fa1 $fa
                mv $fa1.fastq  $fa.end1.fastq
 
                fa2=`ls $Maindir/$folder/*2.fastq | cut -d/ -f8 |cut -d. -f1 | sort -u`;
                fa=`echo $fa2 | cut -d_ -f1-$NosepString`       ## Edit depends on sample center
                echo  $fa2  $fa
                mv $Maindir/$folder/$fa2.fastq  $Maindir/$folder/$fa.end2.fastq
             fi
              
             ## Change format
             f1=`ls *_*end1*`
             mv -v "$f1" $(echo "$f1" | tr '_' '.');
             f2=`ls *_*end2*`
             mv -v "$f2" $(echo "$f2" | tr '_' '.');
          fi
   
    elif [ "$gtStatus" == "Downloaded" ] && [ -f "$tarfile" ] ;then
         echo "Unzip tar file for $folder ........"
         tar -xvf $Maindir/$folder/*tar 
         echo "Unzip gz file for $folder, Don't interupt........ " 
         gunzip *fastq.gz
 
         ## Rename the fastq file
         fastqNumber=`find -name '*fastq' | wc -l`
         echo $fastqNumber fastq files in analysis folder 
         if [ $fastqNumber == 2 ]; then
             echo "fastq file rename with end!!";
             NosepString=`ls *fastq | sed 's/[^_]//g'  | awk '{ print length }'| head -1`
             # NosepString=$[$NosepString -1]

             fa1=`ls $Maindir/$folder/*1.fastq | cut -d/ -f8 |cut -d. -f1 | sort -u`;
             fa=`echo $fa1 | cut -d_ -f1-$NosepString`       ## Edit depends on sample center
             echo $fa1 $fa
             mv $fa1.fastq  $fa.end1.fastq
 
             fa2=`ls $Maindir/$folder/*2.fastq | cut -d/ -f8 |cut -d. -f1 | sort -u`;
             fa=`echo $fa2 | cut -d_ -f1-$NosepString`       ## Edit depends on sample center
             echo  $fa2  $fa
             mv $Maindir/$folder/$fa2.fastq  $Maindir/$folder/$fa.end2.fastq
 
             ## Change format
             f1=`ls *_*end1*`
             mv -v "$f1" $(echo "$f1" | tr '_' '.');
             f2=`ls *_*end2*`
             mv -v "$f2" $(echo "$f2" | tr '_' '.');
         
          elif [ $fastqNumber == 4 ]; then
             echo "combine fastq files."
             cat *1.fastq > $folder.end1.fastq
             cat *2.fastq > $folder.end2.fastq
             
          fi
    else 
         echo "Fail to unzip original file, tar file is already moved to rsrch1!"
    fi               
 fi
  
 ## QC check whether end fastq files are Good to go
 EndfastqNumber=`find -name '*end*.fastq' | wc -l` 

 #######################################################
 ## Prepare PBS files for PRARA-preprocess         #####
 #######################################################
     Prepbsfile=`find -name '*pbs' | wc -l`
     echo "Number of PBS files in $folder = $Prepbsfile"

     if [ $Prepbsfile == 0 ] && [ $EndfastqNumber == 2 ] ; then
         echo "Generate PBS file for PRADA preprocess:"
         fa1=`ls $Maindir/$folder/*end1.fastq | cut -d/ -f8 | sort -u`;
         NosepString=`ls $fa1 | sed 's/[^.]//g'  | awk '{ print length }'`
         NosepString=$[$NosepString -1]
         echo Separate Strings in end fastq file is $NosepString
         fa=`echo $fa1 | cut -d. -f1-$NosepString`
         echo This is input files for PRADA preprocess: $fa1 $fa

         $PRADAdir/prada-preprocess-bi -conf $PRADAdir/conf.$PBSadj.txt -inputdir $Maindir/$folder -sample $fa  -tag $fa -step 2_e1_1 -pbs $folder -outdir $Maindir/$folder -submit no
      
      elif [ $EndfastqNumber == 1 ] || [ $EndfastqNumber == 0 ] ; then
         echo "This sample only has single end sequence, NOT applicable for PRADA!"
      else 
         echo "PBS file for PRADA-preprocess exsit"
      fi
    
 #################################################################################
 ## Submit PBS job of PRADA preprocess depending on space                     ####
 #################################################################################
 
  if grep -q $folder "$Maindir/Prepbs.summary.txt"; then
      echo -e "Already Submited PBS file for PRADA preprocess in this folder $folder"
   
  else  
      echo "PBS file for PRADA preprocess in folder $folder is NOT run yet"
      # echo  $Rjobs $Qjobs
      if [ $ScSPACE -lt $DiskLimit ] || [ $ReSPACE -lt 99 ]; then  
         echo -e "Enough HPC space, Scratch is "$ScSPACE"% and Research1 is "$ReSPACE"%. $Rjobs jobs are running and $Qjobs jobs are in Queue. \n"
       
         ## Determine number of jobs to be submited at this round
         if [[ $Rjobs -lt 30 || $Qjobs -lt 5 || $BAMJobsByTime -lt 8 ]] && [ $ScSPACE -lt 55 ];then
             SubmitJobNoLimit=25
         elif [[ $Rjobs -lt 130 || $Qjobs -lt 20 || $BAMJobsByTime -lt 8 ]] && [ $ScSPACE -lt $DiskLimit ];then
             SubmitJobNoLimit=$Nosubmit
         else 
             SubmitJobNoLimit=3
         fi  
        
         if [ $COUNTER -lt $SubmitJobNoLimit ] && [ $EndfastqNumber == 2 ] ;then  ## count (-1) of the job to be submitted each time
             echo -e "Just submitting No $COUNTER PBS job for PRADA preprocess, this is my jobID.\n"
             msub $Maindir/$folder/$folder.pbs
             ls $Maindir/$folder/$folder.pbs >> $Maindir/Prepbs.summary.txt
             COUNTER=$[$COUNTER +1]
             sleep $[5] 
          elif [ $EndfastqNumber == 1 ] || [ $EndfastqNumber == 0 ]; then
             echo "No PBS file for $folder,This sample was ignored!"
          else 
             echo -e "This moment just submited ($COUNTER-1) PBS jobs for PRADA preprocess, the rest of PBS preprocess jobs are on hold until next one hour. \n" 
          fi 
         
        else 
          echo -e "\n#### HPC disk space warning: ####\n$(date)\n"
          echo -e "\n#### $Cancertype at $(date) ####\n Total Scratch is "$ScSPACE"% full and Total Research is "$ReSPACE"% full.\n";
          df -sh /scratch/genomic_med/xhu3/;
          df -sh /rsrch1/genomic_med/xhu3/;
        fi  
     fi 
   
#############################################################################################
##  After PRADA preprocess finished, Move fastq files to release the space               ####
##  Submit PBS job PRADA-fusion call                                                     ####
#############################################################################################
 PRADAsteps="PIPELINE" 
 logfile="$Maindir/$folder/$folder.log" 
 fastqEndfile=`find -name '*end1.fastq'`
 #echo $fastqEndfile  YYYYYYYY
 PBSfus="$Maindir/$folder/$folder.fusion.pbs"
  
 if [ -f "$logfile" ];then
 Result=`tail -1 $Maindir/$folder/$folder.log | cut -c1-8`
  
   if [ "$Result" == "$PRADAsteps" ]; then
      echo "PRADA preprocess has finished."
        
      if [ ! -f "$PBSfus" ];then
         echo "Checking ReadLength of BAM file and Generating PBS files for fusion call"
         
         ## Check read length of BAM files
         ReadLen=`samtools view $Maindir/$folder/*.bam | head -n 100000 | cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' | sort | uniq -c | cut -c9-10`
         echo -e "ReadLength of recalibrated BAM file is ==$ReadLen=="

         if [[ $ReadLen -gt 45 && $ReadLen -lt 55 ]]; then
            fuslen=40
         elif [[ $ReadLen -gt 70 && $ReadLen -lt 80 ]]; then
            fuslen=60
         else
            echo -e "Read Length Special!!!"
         fi

         sh $Maindir/fusionPBS.generate.sh $Cancertype $folder $fuslen 
         ls $Maindir/$folder/*/*rpkm* >> $Maindir/RPKM.summary.txt  
      else 
         echo "PRADA-fusion PBS file is already generated for $folder."
      fi 
 
      if [ -f "$fastqEndfile" ];then 
         echo "$fastqEndfile exsit, Moving fastq end files $folder to rsrch1,Don't interupt......";
         mkdir $rsrDirfastq/$folder
         ##  mv  $Maindir/$folder/*fastq  $rsrDirfastq/$folder
         rm -rf  $Maindir/$folder/*fastq
      else  
         echo "end.fastq files in $folder NOT exsit any more."
      fi  
      
     ## Fusion PBS job submit if PRADA preprocess finished
     if grep -q $folder "$Maindir/Fuspbs.summary.txt"; then
        echo -e "Already run FUSION PBS file $folder"
     else  
        if [[ $Qjobs -lt 80 && $ScSPACE -lt 97 ]] ; then
           echo "Submit PBS file $folder for fusion call"
           msub $Maindir/$folder/$folder.fusion.pbs
           ls $Maindir/$folder/$folder.fusion.pbs >> $Maindir/Fuspbs.summary.txt
           sleep $[15]
        else 
           echo "PBS file $folder for PRADA-fusion call later..."
        fi  
      fi  
    else  
      echo "PRADA preprocess steps have NOT finished yet"
    fi  
    
  else  
    echo "PRADA preprocess logfile NOT generated yet!!!"
  fi   

 ###################################################################################
 ## move PRADA result to rsrch1                                                 ####
 ###################################################################################

 fusionlogfile="$Maindir/$folder/$folder.fusion.log"
 fusioncallfile="$Maindir/$folder/fusion/$folder.fus.summary.txt"
 
 if [ -f "$fusionlogfile" ];then
    ResultFUS=`tail -1 $Maindir/$folder/$folder.fusion.log | cut -c1-9`
    if [ "$ResultFUS" == "step done" ]; then
        if [ -f "$fusioncallfile" ]; then
           echo $fusioncallfile exist.
           ls $Maindir/$folder/fusion/*fus.summary.txt >> $Maindir/FUSION.summary.txt
           echo -e "FUSION step completed, Moving the resulted files(BAM,RPKM,FUSION,etc) to rsrch1, Don't interupt.....\n"
           mv $Maindir/$folder  $rsrDirRNAseq
           mkdir $Maindir/$folder                              ## Debug GTdownload
           cp $rsrDirRNAseq/$folder/*pbs  $Maindir/$folder     ## Debug generate PBS files
           cp $rsrDirRNAseq/$folder/*log  $Maindir/$folder     ## Debug summary
        else
           echo "FUSION summary for $folder finished and is already moved to rsrch1."
        fi 
    else
        echo "FUSION step has NOT finished yet."
    fi
  fi
  
  ## Count current time
  date2=$(date +"%s")
  diff=$(($date2-$date1))
  echo "$(($diff / 60)) minutes and $(($diff % 60)) seconds elapsed."
  if [ $diff -gt $duration ]; then
     break
  fi
  done < "$filename"
  
  declare i sampleSize
  declare i NOrpkm
  sampleSize=`wc -l $Maindir/$Cancertype.id | cut -c1-3`
  # NOrpkm=`ls $Maindir/*/*/*rpkm* | wc -l | cut -c1-3`
  NOrpkm=`wc -l $Maindir/RPKM.summary.txt | cut -c1-3`
  NOfusion=`wc -l $Maindir/FUSION.summary.txt | cut -c1-3`
  echo -e "\nTotal sample size is $sampleSize, only $NOrpkm RPKM and $NOfusion FUSION results generated so far."
  
  echo -e "\nWaiting for another round of loop over all samples in analysisID list next hours.......\n"

  sleep $[ 10*$TimeInterval ]   ## minimum 30mins loop once
done 
 
## 07-16-2015 end ##

