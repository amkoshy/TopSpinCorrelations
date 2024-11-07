#!/usr/bin/perl

use strict;
use warnings;

my $index = 0;
my @in = read_input($ARGV[0] || die "No input");

my %configs = (
'modelvariation_toprapidity_02.py' => qq(weightVariable = 'toprapidity', slope = 0.2/2.5, weight1x = 1, maxWeight = 100),
'modelvariation_toprapidity_04.py' => qq(weightVariable = 'toprapidity', slope = 0.4/2.5, weight1x = 1, maxWeight = 100),
'modelvariation_toprapidity_06.py' => qq(weightVariable = 'toprapidity', slope = 0.6/2.5, weight1x = 1, maxWeight = 100),
'modelvariation_toprapidity_08.py' => qq(weightVariable = 'toprapidity', slope = 0.8/2.5, weight1x = 1, maxWeight = 100),
'modelvariation_toprapidity_m02.py' => qq(weightVariable = 'toprapidity', slope = -0.2/2.5, weight1x = 1, maxWeight = 100),
'modelvariation_toprapidity_m04.py' => qq(weightVariable = 'toprapidity', slope = -0.4/2.5, weight1x = 1, maxWeight = 100),
'modelvariation_toprapidity_m06.py' => qq(weightVariable = 'toprapidity', slope = -0.6/2.5, weight1x = 1, maxWeight = 100),
'modelvariation_toprapidity_m08.py' => qq(weightVariable = 'toprapidity', slope = -0.8/2.5, weight1x = 1, maxWeight = 100),
'modelvariation_toppt_100_02.py' => qq(weightVariable = 'toppt', slope = 0.2/100, weight1x = 100, maxWeight = 100),
'modelvariation_toppt_100_04.py' => qq(weightVariable = 'toppt', slope = 0.4/100, weight1x = 100, maxWeight = 100),
'modelvariation_toppt_100_06.py' => qq(weightVariable = 'toppt', slope = 0.6/100, weight1x = 100, maxWeight = 100),
'modelvariation_toppt_100_08.py' => qq(weightVariable = 'toppt', slope = 0.8/100, weight1x = 100, maxWeight = 100),
'modelvariation_toppt_100_1.py' => qq(weightVariable = 'toppt', slope = 1.0/100, weight1x = 100, maxWeight = 100),
'modelvariation_toppt_100_m02.py' => qq(weightVariable = 'toppt', slope = -0.2/100, weight1x = 100, maxWeight = 100),
'modelvariation_toppt_100_m04.py' => qq(weightVariable = 'toppt', slope = -0.4/100, weight1x = 100, maxWeight = 100),
'modelvariation_toppt_100_m06.py' => qq(weightVariable = 'toppt', slope = -0.6/100, weight1x = 100, maxWeight = 100),
'modelvariation_toppt_100_m08.py' => qq(weightVariable = 'toppt', slope = -0.8/100, weight1x = 100, maxWeight = 100),
'modelvariation_toppt_100_m1.py' => qq(weightVariable = 'toppt', slope = -1.0/100, weight1x = 100, maxWeight = 100),
);

for my $fn (keys %configs) {
    mv($fn, $configs{$fn});
}

sub read_input {
    my $fileName = shift;
    open my $FH, '<', $fileName or die $!;
    my @res = <$FH>;
    1 while ($res[++$index] !~ /PLACE_MODEL_VARIATION_HERE/ && $index < @res);
    #print "index=$index\n";
    return @res;
}

sub mv {
    my ($fn, $py) = @_;
    #print "$py, $fn\n";
    open my $OUT, '>', $fn or die $!;
    print $OUT @in[0..$index-1];
print $OUT <<END;
        from TopAnalysis.TopUtils.EventWeightDileptonModelVariation_cfi import *
        process.eventWeightDileptonModelVariation = eventWeightDileptonModelVariation.clone(
            $py
        )
        process.topsequence += process.eventWeightDileptonModelVariation
        process.maxEvents.input = int(157.5 * 1143.221)
END
    print $OUT @in[$index+1 .. $#in];
}


