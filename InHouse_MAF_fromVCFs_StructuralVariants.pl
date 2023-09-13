#!/usr/bin/perl
use strict;
use warnings;

#define variables
my $Line;
my $file;
my %vars=();
my %contiglengths=();
my $count=0;
my $binlength=10000;

my($path1, $vcfsfile, $ID, $genesfile)=@ARGV;

#open each vcf file and add variants into %vars
open INPUT2, $vcfsfile or die "Cannot open $vcfsfile\n";
        lineloop2: foreach $file (<INPUT2>){
                chomp $file;
		$count++;
                open INPUT, $file or die "Cannot open $file\n";
	lineloop: foreach $Line (<INPUT>){
                        chomp $Line;
				##Store genome contig lengths
				if($Line=~/\#\#contig\=\<ID\=(\S+?)\,length\=(d+?)\>/){
					my $cont=$1;
					my $totlength=$2;
					if(!exists $contiglengths{$cont}){$contiglengths{$cont}=$totlength; print "Head\: $cont\t$totlength\n";}
					next lineloop;
					}#if head line matches contig name/ID
				if($Line=~/^\#\#/){
                                next lineloop;
                                }#skip other header lines
		my @linesplit2=split(/\t/,$Line); ##split information lines
			if($Line=~/^\#CHROM/){
			my $additionalIDs=(scalar(@linesplit2))-10;
			$count+=$additionalIDs;
			next lineloop;
			}####CHROM header
		my $chr=$linesplit2[0];
		my $SVstartPOS=$linesplit2[1];
		my $SVendPOS=$SVstartPOS;
		my $SVtype="NA";
		if($linesplit2[7]=~/END=(\d+?)\;/){$SVendPOS=$1;}
		if($linesplit2[7]=~/SVTYPE=(\S+?)\;/){$SVtype=$1;}
		if(!exists $contiglengths{$chr}){$contiglengths{$chr}=249000000; print "plus\: $chr\t249000000\n";}
		#Populate %vars with each per 'bin' position represented by SV
		if(!exists $vars{$chr}{$SVtype}){
			for(my $bs=0; $bs<=$contiglengths{$chr}; $bs+=$binlength){
			my $be=$bs+$binlength-1;
				$vars{$chr}{$SVtype}{$bs}{$be}{'total'}=0;
				$vars{$chr}{$SVtype}{$bs}{$be}{'hets'}=0;
				$vars{$chr}{$SVtype}{$bs}{$be}{'homs'}=0;
				$vars{$chr}{$SVtype}{$bs}{$be}{'gain'}=0;
				$vars{$chr}{$SVtype}{$bs}{$be}{'loss'}=0;
				}#for each bin
			}#initialise each per chr/contig structural variant type in %vars

	##Count up number of genotypes per bin
        my $bsp=0;
	for(my $p=$SVstartPOS;$p<=($SVendPOS+$binlength);$p+=$binlength){		#for each per 'bin' POS of SV
	   binloop: for(my $bs=$bsp; $bs<=$contiglengths{$chr}; $bs+=$binlength){
			my $be=$bs+$binlength-1;
			if($bs>$p){last binloop;}
			if($be<($p+1)){$bsp=$bs; next binloop;}
			if(($bs+1)<=$p and $p<=($be+1)){
			$bsp=$bs;
			##loop through each genotype
			loop: for(my $i=9;$i<scalar(@linesplit2);$i++){
				$vars{$chr}{$SVtype}{$bs}{$be}{'total'}+=2;

				##Add gain or loss counts for e.g CNV or LOH calls
				if($linesplit2[8]=~/^RC\:BC\:CN/){
					my @fields=split(/\:/,$linesplit2[$i]);
					if($fields[2]==2){next loop;}
					if($fields[2]==0){$vars{$chr}{$SVtype}{$bs}{$be}{'loss'}+=2;}
					if($fields[2]==1){$vars{$chr}{$SVtype}{$bs}{$be}{'loss'}+=1;}
					if($fields[2]==3){$vars{$chr}{$SVtype}{$bs}{$be}{'gain'}+=1;}
					if($fields[2]>=4){$vars{$chr}{$SVtype}{$bs}{$be}{'gain'}+=2;}
					}
				##Add gain or loss counts for e.g. CNV or LOH calls
				if($linesplit2[8]=~/^GT\:RC\:BC\:CN/){
					my @fields=split(/\:/,$linesplit2[$i]);
					if($fields[3]==2){next loop;}
					if($fields[3]==0){$vars{$chr}{$SVtype}{$bs}{$be}{'loss'}+=2;}
					if($fields[3]==1){$vars{$chr}{$SVtype}{$bs}{$be}{'loss'}+=1;}
					if($fields[3]==3){$vars{$chr}{$SVtype}{$bs}{$be}{'gain'}+=1;}
					if($fields[3]>=4){$vars{$chr}{$SVtype}{$bs}{$be}{'gain'}+=2;}
					}

				if($linesplit2[$i]=~/0\/0/ or $linesplit2[$i]=~/\.\/\./){next loop;} #increment total called genotypes
				if($linesplit2[$i]=~/([0123456789])\/([123456789])/){
					my $first=$1;
					my $second=$2;
					if($first eq $second){
						$vars{$chr}{$SVtype}{$bs}{$be}{'homs'}++;
						next loop;
						} #homozygote for variant ... ref genos 0/0 'looped' out line 39
					if($first ne $second){
						$vars{$chr}{$SVtype}{$bs}{$be}{'hets'}++;
						next loop;
						} #heterozygote
						}# if genotype non-REF
					}#for loop - genotypes
				}##if $p within or equal to $bs and $be

			}#for each bin
		}#for each position + binsize increments

	}#lineloop
        print "$file\n";
        close INPUT;
}#foreach vcf file
close INPUT2;


#open output file and print header
my %rarevars=();
my $output_file = $path1."/InHouse_MAFs_StructuralVariants_".$ID."_".$count."x_participants.txt";
open(OUT, ">>$output_file") || die "Cannot open file \"$output_file\" to write to!\n";
print OUT "\#CHROM\tPOS_BIN_START\tPOS_BIN_END\tSV_TYPE\tTOT_CALLED_GENOS\tHET\tHOM\tCN_LOSS\tCN_GAIN\tMAF_ALL\n";
	foreach my $c (sort keys %vars){
		foreach my $r (sort keys %{$vars{$c}}){
			foreach my $sp (sort {$a<=>$b} keys %{$vars{$c}{$r}}){
				foreach my $ep (sort {$a<=>$b} keys %{$vars{$c}{$r}{$sp}}){
				my $maf="\.";
				my $totalgenotypes=0;
				my $count2=$count*2;
				if($count2<$vars{$c}{$r}{$sp}{$ep}{'total'}){$count2=$vars{$c}{$r}{$sp}{$ep}{'total'};}
				unless($vars{$c}{$r}{$sp}{$ep}{'total'}==0){
					if(($r ne "CNV") and ($r ne "LOH")){$maf=($vars{$c}{$r}{$sp}{$ep}{'hets'}+(2*$vars{$c}{$r}{$sp}{$ep}{'homs'}))/$count2;}
					if(($r eq "CNV") or ($r eq "LOH")){$maf=($vars{$c}{$r}{$sp}{$ep}{'loss'}+$vars{$c}{$r}{$sp}{$ep}{'gain'})/$count2;}
					$totalgenotypes=$vars{$c}{$r}{$sp}{$ep}{'total'}/2;
					#Add rare SV bin coordinates to %rarevars
					if($maf<0.001){$rarevars{$c}{$r}{$sp}{$ep}=$maf;}
					#Print cumulative InHouse SVs to file
					print OUT "$c\t$sp\t$ep\t$r\t$totalgenotypes\t$vars{$c}{$r}{$sp}{$ep}{'hets'}\t$vars{$c}{$r}{$sp}{$ep}{'homs'}\t$vars{$c}{$r}{$sp}{$ep}{'loss'}\t$vars{$c}{$r}{$sp}{$ep}{'gain'}\t$maf\n";
					}#unless total=0
				}#per end pos
			}#per start pos
		}#per SV type
	}#per chr
close OUT;

##Open candidate/selected genes file
my %candG=();
open INPUT3, $genesfile or die "Cannot open $genesfile\n";
	lineloop3: foreach $Line (<INPUT3>){
                        chomp $Line;
			$Line=~s/\"//g;
			my @linesplit3=split(/\t/,$Line);
			if(!exists $candG{$linesplit3[0]}{$linesplit3[1]}{$linesplit3[2]}){
				$candG{$linesplit3[0]}{$linesplit3[1]}{$linesplit3[2]}=$linesplit3[3];
				} ##Genesymbol=ls3;ENSGID=ls4
		}
close INPUT3;

##Re-read through filtered SV files and Print rare SVs to total file
my $output_file2 = $path1."/Rare_StructuralVariants_".$ID."_".$count."x_participants.txt";
open(OUT2, ">>$output_file2") || die "Cannot open file \"$output_file2\" to write to!\n";
print OUT2 "Part_ID\tVariant\tGENOTYPE_CN\tGeneSymbol\tIHMAF\tPOPN\tCONSEQUENCE\tCHROM\tSTARTPOS\tENDPOS\tLENGTH\(Mb\)\n";
#open each vcf file and add variants into %vars
open INPUT2, $vcfsfile or die "Cannot open $vcfsfile\n";
        lineloop2: foreach $file (<INPUT2>){
                chomp $file;
                open INPUT, $file or die "Cannot open $file\n";
		my $ID="NA";
	lineloop3: foreach $Line (<INPUT>){
                        chomp $Line;
				if($Line=~/^\#\#/){
                                next lineloop3;
                                }#skip other header lines
		my @linesplit2=split(/\t/,$Line); ##split information lines
			if($Line=~/^\#CHROM/){
			$ID=$linesplit2[9];
			next lineloop3;
			}####CHROM header
		my $chr=$linesplit2[0];
		my $SVstartPOS=$linesplit2[1];
		my $SVendPOS=$SVstartPOS;
		my $SVtype="NA";
		if($linesplit2[7]=~/END=(\d+?)\;/){$SVendPOS=$1;}
		if($linesplit2[7]=~/SVTYPE=(\S+?)\;/){$SVtype=$1;}
		#Find if present within positional bins of %rarevars
		if(!exists $rarevars{$chr}{$SVtype}){next lineloop3;}
		#loop through each position of rare SVs per chr+type
		my $genoCN=$linesplit2[9];
		if($linesplit2[8]=~/^RC\:BC\:CN$/){$genoCN=~s/\d+\:\d+\:/CN/;}
		if($linesplit2[8]=~/^RC\:BC\:CN\:MCC$/){$genoCN=~s/\d+\:\d+\:/CN/;$genoCN=~s/\:\d+//;}
		if($linesplit2[8]=~/^GT/){$genoCN=~s/\:\S+//;}
		foreach my $sp (sort {$a<=>$b} keys %{$rarevars{$chr}{$SVtype}}){
			foreach my $ep (sort {$a<=>$b} keys %{$rarevars{$chr}{$SVtype}{$sp}}){
				if(($sp<=$SVstartPOS and $ep>=$SVstartPOS) or ($sp<=$SVendPOS and $ep>=$SVendPOS) or ($sp>=$SVstartPOS and $ep<=$SVendPOS)){
					##Add cand genes
					my $CGenes=".";
					if(exists $candG{$chr}){
					foreach my $cst (keys %{$candG{$chr}}){
						foreach my $ced (keys %{$candG{$chr}{$cst}}){
							if(($cst<=$SVstartPOS and $ced>=$SVstartPOS) or ($cst<=$SVendPOS and $ced>=$SVendPOS) or ($cst>$SVstartPOS and $ced<$SVendPOS)){
								$CGenes=$CGenes.", ".$candG{$chr}{$cst}{$ced};
								}#if coords match to the structural variant - add genename to $CGenes
							}#per gene end
						}#per gene start
					}#if chrom name in %CGenes
					$CGenes=~s/\.\, //;
					my $lengMB=($SVendPOS+1-$SVstartPOS)/1000000;
					my $SVtype2=$SVtype;
					if($SVtype=~/BND/){$SVtype2=$SVtype."_".$linesplit2[4];}
				print OUT2 "$ID\t$chr\_$SVstartPOS\_$SVendPOS\_$lengMB\_Mb\t$genoCN\t$CGenes\t$rarevars{$chr}{$SVtype}{$sp}{$ep}\tInHouse$count\t$SVtype2\t$chr\t$SVstartPOS\t$SVendPOS\t$lengMB\n";
				next lineloop3;
				}#if matches within position bin
			}#per end pos
		}#per start pos

	}#lineloop
        close INPUT;
}#foreach vcf file
close INPUT2;
close OUT2;

exit;
