#!/usr/bin/perl -w
#$ -S /usr/bin/perl

# The job is located in the current
# working directory.
#$ -cwd

# we don't indicate the queue, let the queuing system take care of that
#$ -q vangogh

# error messages in (there are none if it's ok, so I prefer to reuse the same file):
#$ -e job.perl.error

# If i'm right there won't be standard output.
#$ -o job.perl.output

# I make inputfiles, though may be not necessary, so I can afterwards and during check what it has done

@temperatures = (3.045, 3.035); 

@ratios = (1.8);

$seed = 822;

foreach $temp (@temperatures)
{
    for($i = 0; $i <@ratios; $i++)
    {
	$x = 1; $z = $ratios[$i];
        $max = $x > $z ? $x : $z;        	
        
	$filename = sprintf("ANI%dx%f.%f",$x,$z,$temp);       
	$outfilename = sprintf("ANIOUT%dx%f.%f",$x,$z,$temp);
	open(FILE,">>".$filename);
	for($j = 1; $j*$max<=30; $j+=4)
	{
	    printf(FILE "  %d   %d   %f  1 100  2000000   4 %d   1046\n%f  0.00000\n", $j*$x,$j*$x, $j*$z,$seed++,$temp);
	}
	# this inputline will halt to program w1n.
	printf(FILE " 0.00 \n 0 0  0   000000 0 00 00000\n");
	close(FILE);
	# file made, we now pipe it to w1n
	@output = `ani < $filename`; 
	# this will have taken some time
	open(OFILE, ">>".$outfilename);
	print OFILE @output;
	close(OFILE);
    }    
} 
