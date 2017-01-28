#!/usr/bin/perl
use POSIX;
use strict;

my $input_trace_file = $ARGV[0];
my $input_reference_file= $ARGV[1];
my $output_trace_file = $input_trace_file."final.dot";
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
my $max=1;
my $ref_line=0;
	  open(OUTPUT,">$output_trace_file")||die "Cannot open the output file $output_trace_file:$!\n";
    open(TRACE_FILE,"$input_trace_file")||die "Cannot open the trace file $input_trace_file:$!\n";
    open(REF_FILE,"$input_reference_file") ||die "Cannot open the ref file $input_trace_file:$!\n";
    $ref_line=<REF_FILE>;
    chomp($ref_line);
    @argv=split(/\s+/,$ref_line);
    $len=$argv[0];
    $max=$argv[1];  
    print $len," ",$max,"\n";
    close(REF_FILE);
    while (<TRACE_FILE>)
   {
     my $line = $_;
     my $line_ori=$line;
     print OUTPUT $line_ori if ($len==0 or $first==1);
     chomp ($line);
     @argv=split(/;/, $line);

     if ($first==0 and $len>0){ 
        my $id=$argv[0];
        open(REF_FILE,"$input_reference_file");
        while(<REF_FILE>){
          $ref_line = $_;
          chomp ($ref_line);
            my @ref_argv=split(/\s+/, $ref_line);
          if($ref_argv[0]==$id){
            #print $ref_argv[0],"\n";
            my $ratio=$ref_argv[1]/$max;
            #$ratio=4 if $ratio>4 ;
            $ratio=log(10*$ratio)/3 if $ratio>0;
            if($ratio>0){
              print OUTPUT "$id node [shape=point,weight=$ratio, height=$ratio];\n";
            }
            else{
              print OUTPUT "$id node [shape=doublecircle];\n";
            }
            last;
          }
      }
      $len=$len-1;
      close(REF_FILE);


     }

    
    
    $first =0 if $first==1 ;
  
   } 
  
   
   
    close (TRACE_FILE);
    close (OUTPUT);
 

 
