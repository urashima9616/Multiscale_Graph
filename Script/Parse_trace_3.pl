#!/usr/bin/perl
use POSIX;
use strict;

my $input_trace_file = $ARGV[0];
my $output_trace_file = $input_trace_file."_ref";
my $src=-1;
my $row_num=0;
my $burst_cnt=0;
my $burst_num=0;
my $avg=0;
my $temp_cnt=0;
my $pkt_num=0;
my $ratio=0;
my @argv;
my $max_num=0;
my $len=0;
my $first=1;
  open(OUTPUT,">$output_trace_file")||die "Cannot open the output file $output_trace_file:$!\n";
    open(TRACE_FILE,"$input_trace_file")||die "Cannot open the trace file $input_trace_file:$!\n";
   while (<TRACE_FILE>)
   {

     my $line = $_;
     chomp ($line);
     @argv=split(/\s+/, $line);
     $len=$argv[1]+2 if $first==1 ;
     #print $len,"\n";

     #print $argv[1],"\n";
     #print  "$_\t" foreach @argv, "\n";
     if ($first==0 and $len>0 ){ 
      my $id=$argv[1];
      my $weight=$argv[2];
      $max_num+=$weight;
    
    }
    
    $len=$len-1 if $first==0 ;
    $first=0;
   } 
  
      close (TRACE_FILE);
$first=1;
open(TRACE_FILE,"$input_trace_file")||die "Cannot open the trace file $input_trace_file:$!\n";
   while (<TRACE_FILE>)
   {

     my $line = $_;
     chomp ($line);
     @argv=split(/\s+/, $line);
     $len=$argv[1]+2 if $first==1 ;
     $max_num=$max_num/$len if $first==1;
     print OUTPUT $len," ",$max_num,"\n" if $first==1;

     #print $len,"\n";

     #print $argv[1],"\n";
     #print  "$_\t" foreach @argv, "\n";
     if ($first==0 and $len>0 ){ 
      my $id=$argv[1];
      my $weight=$argv[2];
      my $ndep=$argv[3];
      
      # print OUTPUT $argv[3+$_]," ",$id," ","\{'weight',$weight\}\n";
    print OUTPUT $id," ",$weight,"\n";
        
      #print "\n";
    }
    
    $len=$len-1 if $first==0 ;
    $first=0;
   } 

   close (TRACE_FILE);


   close (OUTPUT);
 

 
