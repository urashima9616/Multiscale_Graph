#! c:/Perl64/bin/perl -w
  use strict;
  open(HANDLE, "test.txt") or die "cannot open\n";
  my $a;
  my $max=0;
  while(defined($a=<HANDLE>)){
 if($a=~/^Delay from.*[A-Za-z]+([1-9]+)$/){
 print $1, "\n";
  my $temp=$1;
  $max=$temp if $temp>$max ;
   }
 
  }
 
  print $max;
  <>;
