#!/usr/local/bin/perl516


#my $RCS_ID = '$Id$ '

#############################################################################
#
# Name        : get_sample_stat.pl
# Purpose     : 
#
# Author      : Lily Bazak
# Created     : Sep 2013 
#
#############################################################################

######## runnning command ################

# perl516 ~/TCGA/scripts/get_sample_stat.pl LIHC_path_barcode_list ~/TCGA/LIHC/01/edit_out/ LIHC_01.stat get_sample_stat.log

########################


####  Packages to use  ######################################################

use strict;
use FileHandle;
use Bio::Perl;
use Bio::SeqIO::fasta;
use Bio::AlignIO;
use Bio::Tools::Run::StandAloneBlast;
use Bio::AlignIO::bl2seq;
use Data::Dumper;

####  Global variables  #####################################################

my $Prog = "get_sample_stat.pl";

####  The main section  #####################################################

# Analyze arguments

my $argc = @ARGV;
my $i;
for ($i=0;$i<scalar(@ARGV);$i++) {
    print "parameter $i: ",$ARGV[$i],"\n";
}

my ($inFiles) = $ARGV[0]; #path_barcode_list
# my ($inDir) = $ARGV[1];   #directory of ALU outputs
my ($outStat) = $ARGV[1]; #output file
my ($err) = $ARGV[2];     #error file

my $inFilesFh = new FileHandle ("$inFiles","r") || die "$Prog: FATAL: Cannot open input file '$inFiles': $!";
my $outFh = new FileHandle ("$outStat","w") || die "$Prog: FATAL: Cannot open output file '$outStat': $!";
my $errFh = new FileHandle ("$err","w") || die "$Prog: FATAL: Cannot open output file '$err': $!";


#read input file names and get stat

$outFh->print("SampleCode\tnum_of_fragments\tnum_of_uniq_align_fragments\tPercentage_of_uniq_align_fragments\tnum_of_reads_that_were_aligned_to_alu\tPercentage_of_alu_reads_from_all_the_reads\tPercentage_of_alu_reads_from_the_uniq_align_reads\tnum_of_Covered_alu\t");
$outFh->print("Percentage_of_covered_alu_(total_alu_is_1175329)\tnum_of_Covered_Alu_in_refseqs\tPercentage_of_covered_alu_in_refseq\tnum_of_refseq_with_covered_Alu\tnum_of_all_mm_in_Alu\tnum_of_mm_without_SNP\tnum_of_all_mm_in_Alu_-_after_rate_sites2\t");
$outFh->print("num_of_Alu_with_dominant_mm_-_after_clustering\tnum_of_edited_Alu_with_dominant_AG_or_TC\tPercentage_of_edited_alu_of_all_clustered_alu's\tPercentage_of_edited_alu_of_all_covered_alu's\tnum_of_noise_Alu_with_dominant_GA_or_CT\t");
$outFh->print("num_of_edited_Alu_in_refseqs\tPercentage_of_edited_alu_in_refseq\tnum_of_refseq_with_edited_Alu\tnum_of_all_sites_for_Alu_with_dominant_type\tnum_of_edited_Alu_Sites\tPercentage_of_edited_sites\tnum_of_noise_Alu_Sites_(GA_or_CT)\tnum_of_A\tnum_of_I\tIndex(I/A)\t");
$outFh->print("num_of__edited_Alu_in_lincRNA\tnum_of_edited_sites_in_Alu_in_lincRNA\tnum_of_edited_Alu_in_lincRNA_exons\tnum_of_edited_sites_in_Alu_in_lincRNA_exons\t");
$outFh->print("num_of_cvrApos\tnum_of_cvrCpos\tnum_of_cvrGpos\tnum_of_cvrTpos\tnum_of_cvrAs\tnum_of_cvrCs\tnum_of_cvrGs\tnum_of_cvrTs\tmedA\tmedC\tmedG\tmedT\tnum_of_Ads\tnum_of_ins\teditIndex\tnum_of_insWoSnp\teditIndexWoSnp\tnum_of_refCont\tnum_of_mmCont\tcontIndex\tnum_of_mmContWoSnp\tcontIndexWoSnp\n");

while (my $line = $inFilesFh->getline()) { 
    chomp $line;
    #/private2/common/Data/DataSets/InHouse/TCGA/LIHC/RNA-Seq/01/55574a1b-9f06-4758-aa9f-a1d4aa3c8f89/       TCGA-BC-A10Q-01A-11R-A131-07
    
    #in m barcode list: /home/alu/dvirlea/ALU/HpC/in_fastq/C10_1.txt.fastq	C10_s.bam
    
    
    #my ($path,$barcode) = split(/\s/,$line);
    #orig - my ($path,$barcode) = split(/\t/,$line);
    my ($opath,$barcode) = split(/\t/,$line);
    #$errFh->print("path:$path barcode:$barcode\n");
    #my (@dirs) = split(/\//,$path);
    #my $UUID = $dirs[scalar @dirs - 1];
    
    #$errFh->print("$barcode\n"); 
  
    #$outFh->print("$barcode\t$UUID\t");
    $outFh->print("$barcode\t");
    
    my $sampleCntDir = $opath . $barcode . "/counts/";
    #orig - my $sampleCntDir = $inDir . "/" . $barcode . "/counts/";

    #count fastq reads
    my $fastqF = $sampleCntDir . $barcode . "_fastq.cnt_x4";
    open( my $fastqFh, $fastqF);
    my $firstLine = <$fastqFh>;
    close ($fastqFh);
    my ($fastqL,$path) = split(/\s/,$firstLine);
    #orig - my $fastqR = $fastqL/2; #count reads (not pairs) because in the fastq file the number of reads is X4 and by dividing we get X2- number of not paired reads
    my $fastqR = $fastqL; #Reads are already counted through star
    #print "$fastqR\n";
    $errFh->print("$fastqR\t");
    $outFh->print("$fastqR\t");


    #count all uniq aligned reads in bam file 
    my $bamUnq = $sampleCntDir . $barcode . "_unq.cnt";
    open( my $bamUnqFh, $bamUnq);
    my $firstBLine = <$bamUnqFh>;
    # chomp $cntBamUnq;
    close ($bamUnqFh);
    my ($cntBamUnq,$path) = split(/\s/,$firstBLine);
    $errFh->print("$cntBamUnq\t");
    $outFh->print("$cntBamUnq\t");
    eval
    {
    #uniqe aligned reads percentage
    my $bamUnqPer = $cntBamUnq/$fastqR*100;
    $errFh->printf("%d\t%d\t%.2f\t",$cntBamUnq,$fastqR,$bamUnqPer);
    $outFh->printf("%.2f\t",$bamUnqPer); #%2 rount the result 2 digits after the decimal point

    #count aligned reads in bam file, after remove PCR duplicates
    #my $rmDupBam = $sampleCntDir . $barcode . "_rmdup.bam_cnt";
    #open( my $rmDupBamFh, $rmDupBam);
    #my $cntRmDupBam = <$rmDupBamFh>;
    #chomp $cntRmDupBam;
    #close ($rmDupBamFh);
    #$outFh->print("$cntRmDupBam\t");
 
    #count bam aligned reads to Alu
    my $aluBam = $sampleCntDir . $barcode . "_alu.bam_cnt";
    open( my $aluBamFh, $aluBam);
    my $cntAluBam = <$aluBamFh>;
    chomp $cntAluBam;
    close ($aluBamFh);
    $outFh->print("$cntAluBam\t");

    #alu aligned reads percentage
    my $aluInitPer = $cntAluBam/$fastqR*100;
    my $aluUnqPer = $cntAluBam/$cntBamUnq*100;
    $outFh->printf("%.2f\t%.2f\t",$aluInitPer,$aluUnqPer);

    #alu aligned reads - after remove PCR duplicates
    #my $aluWoDup = $sampleCntDir . $barcode . "_alu_rmdup.bam_cnt";
    #open( my $aluWoDupFh, $aluWoDup);
    #my $cntAluWoDup = <$aluWoDupFh>;
    #chomp $cntAluWoDup;
    #close ($aluWoDupFh);
    #$outFh->print("$cntAluWoDup\t");


    #count Covered Alu
    my $cvrAlu = $sampleCntDir . $barcode . "_cvr_alu_list.cnt";
    open( my $cvrAluFh, $cvrAlu);
    my $firstLine = <$cvrAluFh>;
    close ($cvrAluFh);
    my ($cvrAluCnt,$path) = split(/\s/,$firstLine);
    $outFh->print("$cvrAluCnt\t");

    #covered Alu percenatge
    my $allAluCnt = 1175329;
    my $cvrAluPer = $cvrAluCnt/$allAluCnt*100;
    $outFh->printf("%.2f\t",$cvrAluPer);

    #count Covered Alu in refseqs
    my $refDir = $opath . $barcode . "/refseqs/";
    $cvrAlu = $refDir . "/intersect_exp_alu_in_refseqs.cnt_alu";
    open($cvrAluFh, $cvrAlu);
    my $cntAlu = <$cvrAluFh>;
    chomp $cntAlu;
    close ($cvrAluFh);
    $outFh->print("$cntAlu\t");

    #covered Alu in refseq percentage
    my $aluRefseqPer = $cntAlu/$cvrAluCnt*100;
    $outFh->printf("%.2f\t",$aluRefseqPer);
    
    #count refseq with covered Alu
    my $refseq = $refDir . "/intersect_exp_alu_in_refseqs.cnt_refseq";
    open(my $refseqFh, $refseq);
    my $cntRef = <$refseqFh>;
    chomp $cntRef;
    close ($refseqFh);
    $outFh->print("$cntRef\t");
    
    #count all mm in Alu
    my $aluMm = $sampleCntDir . $barcode . "_alu_mm.out_cnt";
    open( my $aluMmFh, $aluMm);
    my $firstLine = <$aluMmFh>;
    close ($aluMmFh);
    my ($aluMmCnt,$path) = split(/\s/,$firstLine);
    $outFh->print("$aluMmCnt\t");

    #count all covered A bases - AA and AG
    my $aluA = $sampleCntDir . $barcode . "_alu_mpileup.cnt_alu_aa_ag";
    open( my $aluAFh, $aluA);
    my $firstLine = <$aluAFh>;
    close ($aluAFh);
    my ($aluAACnt,$aluAGCnt,$aluACnt) = split(/\s/,$firstLine);
    
    #count all covered T bases - TT and TC
    my $aluT = $sampleCntDir . $barcode . "_alu_mpileup.cnt_alu_tt_tc";
    open( my $aluTFh, $aluT);
    my $firstLine = <$aluTFh>;
    close ($aluTFh);
    my ($aluTTCnt,$aluTCCnt,$aluTCnt) = split(/\s/,$firstLine);

    #calc sample editing level - avg of AG and TC - before sites filtering
    my $editBase = ($aluAGCnt+$aluTCCnt)/2;
    my $allBase = ($aluAACnt+$aluTTCnt)/2;
    my $editLevel = $editBase/($allBase+$editBase)*100;
    #$errFh->printf("%s %d %d %.2f\n",$barcode,$allBase,$editBase,$editLevel);
    #$outFh->printf("%.2f\t",$editLevel);

    #count all mm in Alu - excluding SNPs
    $aluMm = $sampleCntDir . $barcode . "_alu_mm.wo_snp_cnt";
    open( $aluMmFh, $aluMm);
    $firstLine = <$aluMmFh>;
    close ($aluMmFh);
    ($aluMmCnt,$path) = split(/\s/,$firstLine);
    $outFh->print("$aluMmCnt\t");

    #count all mm in Alu - after rate_sites2
    $aluMm = $sampleCntDir . $barcode . "_alu_mm.rate_site2_cnt";
    open( $aluMmFh, $aluMm);
    $firstLine = <$aluMmFh>;
    close ($aluMmFh);
    ($aluMmCnt,$path) = split(/\s/,$firstLine);
    $outFh->print("$aluMmCnt\t");

    #count all Alu with mm - after rate_sites2
    #my $alu = $sampleCntDir . $barcode . "_alu_mm.rate_site2_cnt_alu";
    #open( my $aluFh, $alu);
    #$firstLine = <$aluFh>;
    #close ($aluFh);
    #my ($aluCnt,$path) = split(/\s/,$firstLine);
    #$outFh->print("$aluCnt\t");

    #count all Alu with dominant mm - after clustering
    my $alu = $sampleCntDir . $barcode . "_alu_mm.out_stat_cnt";
    open( my $aluFh, $alu);
    $firstLine = <$aluFh>;
    close ($aluFh);
    my ($clusAluCnt,$path) = split(/\s/,$firstLine);
    $outFh->print("$clusAluCnt\t");

    #count edited Alu with dominant AG or TC - after clustering
    $alu = $sampleCntDir . $barcode . "_alu_mm.edit_alu_list_cnt";
    open( $aluFh, $alu);
    $firstLine = <$aluFh>;
    close ($aluFh);
    my ($editAluCnt,$path) = split(/\s/,$firstLine);
    $outFh->print("$editAluCnt\t");

    #edited alu percentage of all clustered alu's (all Alu with dominant mm)
    my $editAluPer = $editAluCnt/$clusAluCnt*100;
    $outFh->printf("%.2f\t",$editAluPer);
    
    #edited alu percentage of all covered alu's
    my $editAluCvrPer = $editAluCnt/$cvrAluCnt*100;
    $outFh->printf("%.2f\t",$editAluCvrPer);

    #count noise Alu with dominant GA or CT - after clustering
    $alu = $sampleCntDir . $barcode . "_alu_mm.noise_alu_list_cnt";
    open( $aluFh, $alu);
    $firstLine = <$aluFh>;
    close ($aluFh);
    my ($noiseAluCnt,$path) = split(/\s/,$firstLine);
    $outFh->print("$noiseAluCnt\t");

    #count edited Alu in refseqs
    $cvrAlu = $refDir . "/intersect_edit_alu_in_refseqs.cnt_alu";
    open($cvrAluFh, $cvrAlu);
    my $cntRefAlu = <$cvrAluFh>;
    chomp $cntRefAlu;
    close ($cvrAluFh);
    $outFh->print("$cntRefAlu\t");

    #edited alu in refseq percentage
    my $editAluRefPer = $cntRefAlu/$editAluCnt*100;
    $outFh->printf("%.2f\t",$editAluRefPer);
    
    #count refseq with edited Alu
    $refseq = $refDir . "/intersect_edit_alu_in_refseqs.cnt_refseq";
    open(my $refseqFh, $refseq);
    $cntRef = <$refseqFh>;
    chomp $cntRef;
    close ($refseqFh);
    $outFh->print("$cntRef\t");


    #count all sites for Alu with dominant type
    my $sites = $sampleCntDir . $barcode . "_alu_sites_all.cnt";
    open( my $sitesFh, $sites);
    $firstLine = <$sitesFh>;
    close ($sitesFh);
    my ($siteCnt,$path) = split(/\s/,$firstLine);
    $outFh->print("$siteCnt\t");

    #count edited Alu Sites (AG or TC)
    $sites = $sampleCntDir . $barcode . "_alu_edit_sites.cnt";
    open( $sitesFh, $sites);
    $firstLine = <$sitesFh>;
    close ($sitesFh);
    my ($editSiteCnt,$path) = split(/\s/,$firstLine);
    $outFh->print("$editSiteCnt\t");

    #edited sites percentage
    my $editSitePer = $editSiteCnt/$siteCnt*100;
    $outFh->printf("%.2f\t",$editSitePer);

    #count noise Alu Sites (GA or CT)
    $sites = $sampleCntDir . $barcode . "_alu_noise_sites.cnt";
    open( $sitesFh, $sites);
    $firstLine = <$sitesFh>;
    close ($sitesFh);
    my ($noiseSiteCnt,$path) = split(/\s/,$firstLine);
    $outFh->print("$noiseSiteCnt\t");

    #sample editing level (Inosines/Adenosines)
    my $sampleDir = $opath . $barcode . "/";
    my $edit = $sampleDir . $barcode . "_all_sample_edit";
    open(my $editFh,$edit);
    $firstLine = <$editFh>;
    chomp $firstLine;
    close($editFh);
    my ($adn,$ins,$edit) = split(/\t/,$firstLine);
    $outFh->print("$adn\t$ins\t$edit\t");

    #lincRNAs
    my $samplelncRNADir = $opath . $barcode . "/lincRNA/";

    #count edited Alu in lincRNA
    my $aluLincRNA = $samplelncRNADir . "intersect_merge_lncRNA_edit_alu_strand_cnt";
    open( my $aluLincRNAFh, $aluLincRNA);
    $firstLine = <$aluLincRNAFh>;
    my ($cntAlu,$path) = split(/\s/,$firstLine);
    close ($aluLincRNAFh);
    $outFh->print("$cntAlu\t");

    #count edited sites in Alu in lincRNA
    my $siteLincRNA = $samplelncRNADir . "intersect_merge_lncRNA_edit_sites_strand_cnt";
    open( my $siteLincRNAFh, $siteLincRNA);
    $firstLine = <$siteLincRNAFh>;
    my ($cntSites,$path) = split(/\s/,$firstLine);;
    close ($siteLincRNAFh);
    $outFh->print("$cntSites\t");

    #count edited Alu in lincRNA exons
    $aluLincRNA = $samplelncRNADir . "intersect_merge_lncRNA_exons_edit_alu_strand_cnt";
    open( $aluLincRNAFh, $aluLincRNA);
    $firstLine = <$aluLincRNAFh>;
    ($cntAlu, $path) = split(/\s/,$firstLine);
    close ($aluLincRNAFh);
    $outFh->print("$cntAlu\t");

    #count edited sites in Alu in lincRNA exons
    $siteLincRNA = $samplelncRNADir . "intersect_merge_lncRNA_exons_edit_sites_strand_cnt";
    open( $siteLincRNAFh, $siteLincRNA);
    $firstLine = <$siteLincRNAFh>;
    ($cntSites, $path) = split(/\s/,$firstLine);
    close ($siteLincRNAFh);
    $outFh->print("$cntSites\t");


    #new index details (Auguest 2015)
    my $sampleNewIndexDir = $opath . $barcode . "/edit_index/";
    
    my $sampleNewDetails = $sampleNewIndexDir . $barcode . "_all_alu_edit_sample_level_plus.out";
  #  print "$sampleNewDetails\n";
    open(my $sampleNewDetailsFh, $sampleNewDetails);
    $firstLine = <$sampleNewDetailsFh>;
   # print "$firstLine";
    chomp($firstLine);
    close($sampleNewDetailsFh);
    $outFh->print("$firstLine");
    
    
    $outFh->print("\n");
    } || do
    {
    $outFh->print("\n");
    print "Something is wrong with sample $barcode!\n";
    
    };

}


0;

# End of main of "get_sample_stat.pl"

####  Functions section  ####################################################


#############################################################################

# End of script "get_sample_stat.pl"


