#!/usr/bin/perl -w
#===============================================================================
#
#         FILE: getpmf.pl
#
#        USAGE: ./getpmf.pl gwham_output  
#
#  DESCRIPTION: Extract the PMFs from the output of gwham
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Dejun Lin (), dejun.lin@gmail.com
# ORGANIZATION: University of Rochester Medical Center
#      VERSION: 1.0
#      CREATED: 04/13/2015 02:40:41 PM
#     REVISION: ---
#===============================================================================

use FileHandle;

my $wham = shift;
my $fh = FileHandle->new($wham,"r");

my $flag = 0;
while(<$fh>) {
  if(/^#\s*Bin\s*Vals\s*PMF\s*RhoNormalized\s*$/) {
    $flag = 1;
    print;
    next;
  }
  next unless $flag;
  if($flag && /^#/) { last; }
  print;
}
