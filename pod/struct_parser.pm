package struct_parser;
#!/usr/bin/env perl
# Author: Yulong Li <liyulong12@mails.ucas.ac.cn>
# process STRUCTURE and ADMIXTURE results
use strict;
use warnings;
use POSIX;
use List::MoreUtils qw{uniq};
use Archive::Zip;
use File::Basename 'basename';
use Sort::Key::Natural qw( natsort );
use constant V    => '231219';
use constant BIG  => 99999;
sub parse_struct {
	#
	# reps => K => {'run'=>$run, 'estLnPk'=>$LnK, 'mea_ln'=>$mea_ln, 'var_ln'=>$var_ln,
	# 'Q'=>$clsts, 'indv'=>$$indvs, 'nloci'=>$loci}
	#
	# parse_struct(\@files, $in_path, $grp, $out_path)
	#
	# return ($clsts, \@K, \@files, $repss, $ma)
	#
	
	my $files   = shift;
	my $in_path = shift;
	my $grp     = shift;
	my $out_path= shift; # Path used to output errors.
	my @files   = @$files;
	my ($indivs, $loci, @K, $ma, $clsts, $repss); 
	
	my ($pops, $tmp, $N);
	if ($grp) {
		if ($grp =~ /^(\d+,)*\d+,?$/g) {
			# newpop size.
			my @parts = split(',', $grp);
			my $line = 1;
			for (my $i=0; $i<@parts; $i++) {
				my $size = $parts[$i];
				$N += $size;
				for (my $j=$line; $j<=$N; $j++) {
					$pops->{$j} = 'pop'. ($i+1);
					$line++; # line number.
				}
			}
			
		} else {
			($pops, $tmp) = parse_popmap($grp);
			$N = scalar(keys %$pops);
		}
	}
	
	#
	#
	#
	foreach my $file (@files) {
		my ($tclsts, $k, $indv, $LnK, $mea_ln, $var_ln);
		my $h;
		open(my $in_fh, "$in_path/$file") or die "$!";
		
		while(<$in_fh>) {
			chomp;
			$_ =~ s/\r//g; # This is important.
		
			if (/(\d+)\s+individuals/) {
				$indivs = $1;
			}
			if (/(\d+)\s+loci/) {
				$loci =$1;
			} 
			if (/(\d+)(\s+)populations(\s+)assumed/) {
				$k = $1;
				push @K, $k;
			}
		
			if (/(.+?)\s+=\s(.+)/) {
			#
			# Estimated Ln Prob of Data   = -9092.6
			# Mean value of ln likelihood = -7201.2
			# Variance of ln likelihood   = 3782.8
			# Mean value of alpha         = 0.0772
			#
				$LnK    = $2 if $1 eq 'Estimated Ln Prob of Data';
				$mea_ln = $2 if $1 eq 'Mean value of ln likelihood';
				$var_ln = $2 if $1 eq 'Variance of ln likelihood';
			}
		
			if (/(\d+)\s+.+\s+(\(\d+\))\s+(\d+)?\s+:\s+(.+)/) {
			#
			# Inferred ancestry of individuals:
			# 		Label (%Miss) Pop:  Inferred clusters
			# 1     JM01    (0)    1 :  0.122 0.113 0.110 0.023 0.115 0.122 0.110 0.129 0.132 0.024 
			# 2     JM02    (0)    1 :  0.124 0.167 0.118 0.006 0.110 0.107 0.118 0.101 0.127 0.021 
			# 3     JM03    (0)    1 :  0.105 0.143 0.134 0.008 0.130 0.107 0.137 0.103 0.116 0.015 
			# 4     JM04    (0)    1 :  0.130 0.161 0.116 0.008 0.107 0.101 0.118 0.119 0.114 0.024 
			#	
				my $pop;
				
				if ($grp) {
					#
					# group option
					#
					
					if ($N ne $indivs) {
						my $info = "PopMap (group size) error! Number of individuals are not consistent with that in file: $file\n";
                        $info   .= ", in popmap: $N, in file: $indivs\n";
						$info   .= ", maybe you should convert your file to UNIX format...\n";
						die $info;
					}
					$pop = $pops->{$1};
				} else {
					if (!$3) {
						# error, no pop indicator
						my $info = "No population indicator in file: $file, you should turn on the PopData option (check the box 'Putative population origin for each individual') in STRUCTURE!";
						$info   .= "\nOr upload a popmap file, or input a vector of new pop sizes!";
						die $info;
					}
					$pop = 'pop'.$3; # predefined pop.
				}
				
				$indv++;
				my $dat = $4;
				$dat    =~ s/^\s+|\s+$//g;
				
				my @parts = split(/\s+/, $dat);
				my $len   = @parts;
				
				if ($k < $len) {
					#
					# remove extra data after Q-matrix if exists
					#
					splice (@parts, $k-$len);
					
				} elsif ($k > $len) {
					#
					#
					#
					die "The number of columns (in Q-matrix) is smaller than K.";
				}
				
				die "K is not equal to clusters in STRUCTURE Q-matrix\n" if $k ne @parts;
			
				for (my $i=0; $i<@parts; $i++) {
					#each membership (=ancestry) coefficient (Q).
					#
					# return array of pops containing each clusters.
					#
					push @{$clsts->{$file}->{$i+1}->{$pop}}, $parts[$i]; # array of each cluster of each pop.
					push @{$ma->{$file}->[$i]}, $parts[$i];
				}
			}	 
			next;
		}
	
        if (!$k) {
			my $info = "Can not get k in $file, ";
            $info   .= "please check your file.\n";
			die $info;
        }
        if (!$indv || $indv != $indivs) {
			my $info = "Individuals are not identical in $file";
			$info   .= ", expected: $indivs, I get: $indv, \n ";
            $info   .= "please check your file.\n";
			die $info;
        }
		close $in_fh;
		
		my @f_parts = split(/_/, $file);
		my $run = $f_parts[$#f_parts-1];
		push @{$repss->{$k}}, {'run'=>$run, 'estLnPk'=>$LnK, 'mea_ln'=>$mea_ln, 'var_ln'=>$var_ln,
		'Q'=>$clsts, 'indv'=>$indivs, 'nloci'=>$loci};
	}
	
	@K = uniq @K;
	@K = sort{$a<=>$b} @K;
	
	return ($clsts, \@K, \@files, $repss, $ma);	
}

sub parse_admix {
	#
	# parse_admix(files, $in_path, popmap, $out_path, $grp, convert?)
	#
	# return ($clsts, \@K, \@files, $pops, $ma, $likelihood)
	#
	my $files   = shift;
	my $in_path = shift;
	my $pop     = shift; # this also hold new pop size.
	my $out_path= shift; # Path used to output errors.
	my $log     = shift;
	my $grp     = shift; # now useless.
	my $convert = shift;
	
	my ($cnt, @files, $clsts, $reps, @K, @order, $pops, $N, $tmp, $ma);
	my ($popsize, $indivs_p, $qtable);
	@files = @$files;
	
	#parse popmap
	
	if ($pop =~ /^(\d+,)*\d+,?$/g) {
		# newpop size.
		my @parts = split(',', $pop);
		my $line = 1;
        $popsize = scalar(@parts);
		for (my $i=0; $i<$popsize; $i++) {
			my $size = $parts[$i];
			$N += $size;
			for (my $j=$line; $j<=$N; $j++) {
				$pops->{$j} = 'pop' . ($i+1);
				push @order, 'pop' . ($i+1);
				$line++; # line number.
			}
		}
	} else {
		($pops, $tmp, $popsize) = parse_popmap($pop);
		@order = @$tmp;
		$N = scalar(keys %$pops);
	}
	
	if (defined $convert) {
		# Convert to labels for plot.
		open (my $out_fh, ">$in_path/labels.txt") or die "$!";
		for(my $i=0; $i<@order; ++$i) {
			print $out_fh $i+1, "\t", $order[$i], "\n";
		}
		close $out_fh;
	}
	
	foreach my $file (natsort @files) {
    
		my ($tclsts, $k, $popname, $popord, $indivs, $out_fh, @ind_table_q);
		open(my $in_fh, "$in_path/$file") or die "$!";
		
		if (defined $convert) {
			open ($out_fh, ">$in_path/$file.strq") or die "$!";
		}
		while(<$in_fh>) {
			#
			#
			#
			$_ =~ s/[\r\n]|^\s+|\s+$//g; # This is important, remove space and tab lines.
			next if /^#|^$/;
			
			$indivs++;
			
			if (not defined $pops->{$.}) {
				my $info = "PopMap (group size) error! Number of individuals are not consistent with that in file: $file,
                please check your popmap file.\n";
				die $info;
			}
			
			my @parts = split(/\s+/);
			$k = @parts;
			push @K, $k;
			
			if (!$popname || $popname ne $pops->{$.}) {
				$popord += 1; # Order of Pops.
				$popname = $pops->{$.};
				
			}
			
			for (my $i=0; $i<$k; $i++) {
					#each membership (=ancestry) coefficient (Q).
					#
					# return array of pops containing each clusters.
					#
					push @{$clsts->{$file}->{$i+1}->{$popname}}, $parts[$i]; # array of each cluster of each pop.
					push @{$ma->{$file}->[$i]}, $parts[$i];
					
			}
			
            my $tmp = join(' ', $., $., '(0)', $popord, ':').' '.$_."\n";
            push @ind_table_q, $tmp;
			if (defined $convert) {
				# Convert to structure Q format.
					
				print $out_fh $tmp;
					
			}
		}
		close $in_fh;
		close $out_fh if defined $convert;
		
		if ($N != $indivs) {
			my $info = "PopMap (group size) error! Number of individuals are not consistent with that in file: $file\n";
            $info   .= ", in popmap: $N, in file: $indivs\n";
            $info   .= ", please check your files...\n";
			die $info;
		}
        if ($indivs_p && $indivs_p != $indivs) {
            die "Number of individuals is not consistent in $file.\n";
        } 
        $indivs_p = $indivs;
        
        # store q tables.
        (my $a, my $pop_table_q) = calc_popq($clsts->{$file}, \@order, $popsize);
        push @{$qtable->{'indq'}->{$k}}, \@ind_table_q;
        push @{$qtable->{'popq'}->{$k}}, \@$pop_table_q;
        push @{$qtable->{'file'}->{$k}}, $file;
        
	}
	@K = uniq @K;
	@K = sort{$a<=>$b} @K;
	
	#
	# parse logs likelihood
	#
	my $likelihood;
	if ($log) {
		
		($likelihood, my $flag) = parse_logs($in_path, scalar(@files));
		
		if ($flag == 0) {
			my $info = "Number of logs are not consistent with Q files.\n";
			die $info;	
		}
		if ($flag == 2) {
			my $info = "Log files are wrong!\n";
			die $info;	
		}
	}
    
	# 1:popmap. 2: number of files. 
    # 3:number of K. 4: number of individuals. 
    push my @summary, ($pops, \@files, \@K, $indivs_p);
	
    return ($clsts, \@K, \@files, $pops, $ma, $likelihood, \@summary, $qtable);
}

sub calc_deltak {
	#
	# calc_deltak (K, repss, flag).
	# flag: pass mean and var likelihood directly.
	# rerurn and a success flag 
	my $K     = shift;
	my $repss = shift;
	my $flag  = shift;
	my @K     = @{$K};
	my $delta;

	
	if (!$flag) {
		# calculate Mean LnPK and std. 
		for (my $i=0; $i<@K; ++$i) {
			my $k = $K[$i];
			my (@iter, @LnPks);
			map {push @iter, $_->{'run'}; push @LnPks, $_->{'estLnPk'}} @{$repss->{$k}};
			my ($avg, $sum, $var)  = mean_and_var(\@LnPks);
			$delta->{$k}->{'estLnPk'} = $avg;
			$delta->{$k}->{'std'}  = $var;
			$delta->{$k}->{'iter'} = @iter;
		}
	} else {
		# just pass mean and var
		my $ll   = $repss;
		my $reps = $flag;
		
		for (my $i=0; $i<@K; ++$i) {
			my $k = $K[$i];
			$delta->{$k}->{'estLnPk'} = $ll->{$k}->{'ll'};
			$delta->{$k}->{'std'}     = $ll->{$k}->{'llstd'};
			$delta->{$k}->{'iter'}    = $reps->{$k};
		}
	}
	#
	# calculate L''(K)
	# lnpk   => L'(K)  = L(K) -L(K-1)
	# lnppk  => L''(K) = abs(L(K+1) - 2L(K) + L(K-1))
	# deltak => DeltaK = L''(K)/stdv
	#
    
    my $re_flg = 1; # return flag.
	
    for (my $i=1; $i<scalar(@K); ++$i) {

		my $k = $K[$i];
		my $this = $delta->{$k}->{'estLnPk'};
		my $var  = $delta->{$k}->{'std'};
		my $befr = $delta->{$k-1}->{'estLnPk'};
        if (defined $this && defined $befr) { 
            $delta->{$k}->{'lnpk'} = $this - $befr; # can from second to the last one.
        }
        if ($i != scalar(@K) - 1) {
			my $aftr = $delta->{$k+1}->{'estLnPk'};
            if (!(defined $this && defined $befr && defined $aftr)||($var < 0.000001)) {
                #
                # var is too small.
                # update 2018.01.26
                # not continuious
                # update 2019.01.17
                $re_flg = 0;
                next;
                
            }
			my $lnkk = abs($aftr - 2*$this + $befr);
            
			$delta->{$k}->{'lnppk'} = $lnkk;
			$delta->{$k}->{'deltak'} = $lnkk/$var;
		}
        
	}
	
	return ($delta, $re_flg);

}

sub calc_medk {
	#
	# calc_medk(clsts, thred, K)
	#
	
	#
	# step1. Calculate mean and median of each pop.
	#
	my $clsts = shift;
	my $thred = shift;
	my $K     = shift;
	my $reps;
	
	for my $run (keys %$clsts) {
		
		my $tclsts;
		my $k_run = scalar(keys %{$clsts->{$run}});
		for my $k (keys %{$clsts->{$run}}) {
            # assert each cluster.
			my $t->{'med'} = 0;
			$t->{'avg'} = 0;
			foreach my $pop (keys %{$clsts->{$run}->{$k}}) {
				my @tmp = sort{$b<=>$a} @{$clsts->{$run}->{$k}->{$pop}};	
				my $med         = median(\@tmp);
				my ($avg, $sum) = mean(\@tmp);
				#
				# may be faster here?
				#
				if ($med > $thred) {
					$t->{'med'}++;
				} 
				if ($avg > $thred) {
					$t->{'avg'}++;
				}
                if ($t->{'med'} >= 1 && $t->{'avg'} >= 1) {
                    last;
                }
			}
			$tclsts->{'med'} ++ if $t->{'med'} >= 1; # True clusters of median.
			$tclsts->{'avg'} ++ if $t->{'avg'} >= 1; # True clusters of mean.
			$tclsts->{'med'} = 0 if not defined $tclsts->{'med'};
			$tclsts->{'avg'} = 0 if not defined $tclsts->{'avg'};
		}
		#
		# for each reps.
		#
		push @{$reps->{$k_run}->{'avg'}}, $tclsts->{'avg'};
		push @{$reps->{$k_run}->{'med'}}, $tclsts->{'med'};
	
	}
	
	#
	#
	# Step3. Calculate Median and max for each K.
	#
	my @K = @{$K};
	my $maxk = $K[$#K];
	my $mink = $K[0];
	my (@medmed, @medmean, @maxmed, @maxmean, $median, $max, $iter);
	for (my $i=0; $i<@K; ++$i) {
		
		my $k = $K[$i];
		$iter->{$k} = @{$reps->{$k}->{'avg'}}; # Num of reps of each K.
		my @tmp1 = sort{$b<=>$a} @{$reps->{$k}->{'avg'}}; #means of reps.
		my @tmp2 = sort{$b<=>$a} @{$reps->{$k}->{'med'}}; #medians of reps.
		$median->{$k}->{'avg'} = upmedian(\@tmp1);
		$median->{$k}->{'med'} = upmedian(\@tmp2);
		$max->{$k}->{'avg'}    = $tmp1[0];
		$max->{$k}->{'med'}    = $tmp2[0];
		push @medmed, $median->{$k}->{'med'}; 
		push @medmean, $median->{$k}->{'avg'};
		push @maxmed, $max->{$k}->{'med'};
		push @maxmean, $max->{$k}->{'avg'};
	}	
	
	#
	# Step4. Calculate max across all K.
	#
	
	my @max1 = sort{$b<=>$a} @medmed;
	my @max2 = sort{$b<=>$a} @medmean;
	my @max3 = sort{$b<=>$a} @maxmed;
	my @max4 = sort{$b<=>$a} @maxmean;
	my $med  = {'medmed'=>\@medmed, 'medmean'=>\@medmean, 'maxmed'=>\@maxmed, 'maxmean'=>\@maxmean};
	my $medk = {'medmed'=>$max1[0], 'medmean'=>$max2[0], 'maxmed'=>$max3[0], 'maxmean'=>$max4[0]};
	return ($med, $medk, $iter);
	#return (\@medmed, \@medmean, \@maxmed, \@maxmean);
}

sub choos_k {
	#
	# algorithm of fastStructure
	# nil Raj, Matthew Stephens, and Jonathan K. Pritchard. fastSTRUCTURE: Variational Inference of 
	# Population Structure in Large SNP Data Sets , (Genetics) June 2014 197:573-589
	# we compute the ancestry contribution of each model component as the mean admixture proportion over all samples.
	# The number of relevant model components is then the minimum number of populations that have a cumulative ancestry contribution of at least 99.99%,
	#
	# usage: choos_k(in_path, ma)
	#
	my $in_path = shift;
	my $flag    = shift; # Flag for pass Q.
	my ($cnt, @files, @bestK);
	if ($flag) {
		#
		# $flag->{'run'}; run of Q file.
		# K: clusters
		# N: num of individuals.
		# ma: KxN matrix of Q.
		#
		foreach my $run (keys %{$flag}) {
			my (@Q, @C, $N, $best);
			my @ma = @{$flag->{$run}};
			for my $i (@ma) {
				my $sum = 0;
				map {$sum += $_} @$i;
				push @Q, $sum;
			}
			$N = @{$ma[0]};
			@Q = sort{$b<=>$a} @Q;
			push @C, $Q[0];
			for (my $i=1; $i<@Q; $i++) {
				push @C, $Q[$i] + $C[$i-1];
			}
			$best = 1;
			for (my $i=0; $i<@C; $i++) {
				$best++ if ($C[$i] < $N - 1);
			}
			push @bestK, $best;
			$cnt->{'f'}++;
		}	
	} else {
		#
		# read from files.
		#
		opendir (D, $in_path) or die "$!";
		while ((my $file = readdir(D))) {
		
			next if $file !~ /.+\.Q$|.+\.meanQ$/; # Requier .Q for admixture or .meanQ for fastStructure.
			push (@files, $file);
			$cnt->{'f'}++; # Num of files..
		}
	
	
		foreach my $file (@files) {
		
			open(my $in_fh, "$in_path/$file") or die "$!";
			my (@Q, @C, $N, $best);
			while(<$in_fh>) {
				next if /^#|^$/;
				chomp;
				my @parts = split;
				for (my $k=0; $k<@parts; $k++) {
					$Q[$k] += $parts[$k];
				}
				$N++; # Num of individuals.
			}
			@Q = sort{$b<=>$a} @Q;
			push @C, $Q[0];
			for (my $i=1; $i<@Q; $i++) {
				push @C, $Q[$i] + $C[$i-1];
			}
			$best = 1;
			for (my $i=0; $i<@C; $i++) {
				$best++ if ($C[$i] < $N - 1);
			}
			push @bestK, $best;	
		}
	}
	
	my (@K, %bin, $best);
	map {$bin{$_}++} @bestK;
	@K = sort{$a<=>$b} keys %bin;

	$best = $K[0]; # start at the minimum K.
	for (my $i=1; $i<@K; $i++) {
		$best = $K[$i] if $bin{$K[$i]} > $bin{$best}; # find the minimum K that comes at the most time.
	}
	return $best;
	
}

sub parse_f {
	#
	# return files that match the indicated suffix.
	#
	my $in_path = shift;
	my $fmt     = shift;
	my ($cnt, @files);
	$fmt = $fmt eq 'admixture' ? '.+\.Q$|.+\.meanQ' : $fmt;
	$fmt = $fmt eq 'structure' ? '_f' : $fmt;
	opendir (D, $in_path) or die "$!";
	while ((my $file = readdir(D))) {
		next if $file =~ /^\./;
        next if $file !~ /$fmt/;
		push (@files, $file);
		$cnt->{'f'}++; # Num of files..
	}
	return \@files;
}

sub parse_popmap {
	#
	# parse PopMap, return pops->{indiv} and the order of pops.
	# update 2018.02.20
    #
	my $pop = shift;
	my (@order, $pops, $popsize, $popid);
	open(my $in_fh, "$pop") or die "No PopMap file!";
	my $i=0;
    my $j=0;
	while (<$in_fh>) {
		$_ =~ s/[\r\n]|^\s+|\s+$//g; # This is important, and remove space or tab lines.
		next if /^#|^$/;
		my @part = split;
        my $name = $part[1];
		$i++; # individual order. 
        if (!$popsize->{$name}) {
            $j++; # pop id.
            push @order, $name;
            $popid->{$name} = $j; # convert pop name to pop id.
        }
        $popsize->{$name}++;
        $pops->{$i} = $popid->{$name};
	}
	close $in_fh;
	return($pops, \@order, $popsize);
}

sub check_format {
	#
	# Check file format
	# check_format(fmt, in_path, popmap?)
	#
	my $fmt     = shift;
	my $in_path = shift;
	my $pop     = shift;
	my ($in_file, $indiv, $test, $msg);
	if ($fmt eq 'str') {
		my @files;
		opendir (D, $in_path) or die "$!";
		while ((my $file = readdir(D))) {
			next if $file =~ /^\./;
			next if $file !~ /_f$/;
			push (@files, $file);
		}
		foreach my $file (@files) {
			$test = 0;
			open(my $in_fh, "$in_path/$file") or die "$!";
			while(<$in_fh>) {
				chomp;
				$_ =~ s/\r//g; 
				if (/(\d+)(\s+)populations(\s+)assumed/) {
					$test = 1
				}
				if (/(\d+)\s+individuals/) {
					$indiv = $1;
				}
			}
			close $in_fh;
			if (!$test) {
				#raise error...
				$msg = "error! $file is not a correct STRUCTURE file.\n";
				last;
			}
		}
		
	} elsif ($fmt eq 'admix') {
		my @files;
		opendir (D, $in_path) or die "$!";
		while ((my $file = readdir(D))) {
			next if $file =~ /^\./;
			next if $file !~ /.+\.Q$|.+\.meanQ/;
			push (@files, $file);
		}
		foreach my $file (@files) {
			$test = 0;
			open(my $in_fh, "$in_path/$file") or die "$!";
			$indiv = 0;
			while(<$in_fh>) {
				chomp;
				$_ =~ s/\r//g;
				my @parts = split;
				map {$test = 1 if /\d+/} @parts;
				$indiv++;
			}
			close $in_fh;
			if (!$test) {
				#raise error...
				$msg = "error! $file is not a correct Q-matrix (ADMIXTURE fastStructure) file.\n";
				last;
			}
		}
	} elsif ($fmt eq 'popmap') {
		my $num;
		($msg, $indiv) = check_format('admix', $in_path);
		open(my $in_fh, "$pop") or die "No PopMap file!";
		while (<$in_fh>) {
			chomp;
            $_ =~ s/\r//g;
            next if /^#|^$/;
			$num++;
		}
		close $in_fh;
		if ($num != $indiv) {
			#print $num, "\n", $indiv, "\n";
			$msg .= "Number of individuals in PopMap file are not identical to that in Q-matrix (ADMIXTURE fastStructure) file.\n";
		}
	}
	
	return($msg, $indiv);
	
}

sub parse_logs {
	#
	# parse_logs(in_path, cnt)
	# cnt: number of Q files.
	# return mean likelihood and a flag.
	#
	my $in_path = shift;
	my $cnt     = shift;
	my @logs    = @{parse_f($in_path, '.log')};
	if (@logs == $cnt) {
		my $likelihood;
		foreach my $file (@logs) {
			open(my $in_fh, "$in_path/$file") or die "$!";
			my $k;
			if ($file =~ /.+\.(\d+)\.log/) {
				$k = $1;
			}
			
			while (<$in_fh>) {
				chomp;
				if (/Marginal\s+Likelihood\s+=\s+(.+)/) {
					#
					# fastSTRUCTURE
					#
					push @{$likelihood->{$k}->{'ll'}}, $1;
				} 
				if (/^Loglikelihood:\s+(.+)/) {
					#
					# admixture
					#
					push @{$likelihood->{$k}->{'ll'}}, $1;
				}
				if (/^CV\s+error\s+\(K=(\d+)\):\s+(.+)/) {
					#
					# admixture CV error
					#
                    $k = $1; # add k if available.
					push @{$likelihood->{$k}->{'cv'}}, $2;	
				}
				if (/^CV\s+error\s+=\s+(.+),(.+)/) {
					#
					# fastSTRUCTURE CV error
					#
					push @{$likelihood->{$k}->{'cv'}}, $1;	
				}
			}
		
		}
		
		my $ll;
		foreach my $k (keys %{$likelihood}) {
			
			if ($likelihood->{$k}->{'ll'}) {
				my @lll = @{$likelihood->{$k}->{'ll'}};
				($ll->{$k}->{'ll'}, my $tmp, $ll->{$k}->{'llstd'}) = mean_and_var(\@lll);
			} else {
				# log is wrong.
				return(0,2);
			}
			if ($likelihood->{$k}->{'cv'}) {
				# if cv error is available.
				my @cv = @{$likelihood->{$k}->{'cv'}};
				($ll->{$k}->{'cv'}, my $tmp, $ll->{$k}->{'cvstd'}) = mean_and_var(\@cv);
			}
		}
		return ($ll, 1);
	} else {
		#
		# files are not consistent.
		#
		
		return(0,0);
		
	}
}

sub upmedian {
	#
	#calculate the ceil median 
	#the array must be sorted!
	my $num = shift;
	my @tmp = @{$num};
	my $cnt = @tmp;
	my $med = 0;
	
	if ($cnt%2 == 0) {
		#
		# return ceil if the median is not a int.
		#
		# the index of array start at 0 !!!
		#
		$med += $tmp[$cnt/2-1];
		$med += $tmp[$cnt/2];
		$med  = ceil($med/2);
	} else {
		$med  = $tmp[($cnt-1)/2];
	}
	return $med;
}

sub median {

	#the array must be sorted!
	my $num = shift;
	my @tmp = @{$num};
	my $cnt = @tmp;
	my $med = 0;
	
	if ($cnt%2 == 0) {
		#
		# the index of array start at 0 !!!
		#
		
		$med += $tmp[$cnt/2] + $tmp[$cnt/2-1];
		$med /= 2;
	} else {
		$med  = $tmp[($cnt-1)/2];
	}
	return $med;
}

sub mean {
	
	my $num = shift;
	my @tmp = @{$num};
	my $cnt = @tmp;
	my ($average, $sum);
	$average = 0;
	$sum     = 0;
	map {$sum += $_} @tmp;
	$average = $sum/$cnt;
	return($average, $sum);
}

sub mean_and_var {
	
	my $num = shift;
	my @tmp = @{$num};
	my $cnt = @tmp;
	my ($average, $sum, $var, $ss);
	$average = 0;
	$sum     = 0;
	$var     = 0;
	$ss      = 0;
	map {$sum += $_} @tmp;
	$average = $sum/$cnt;
	map {$ss += ($_ - $average)**2} @tmp;
	if ($cnt == 1) {
		$var = 0;
	} else {
		$var = sqrt($ss/($cnt-1));
	}
	return($average, $sum, $var);
}

sub parse_struct_new {
	#
	# reps => K => {'run'=>$run, 'estLnPk'=>$LnK, 'mea_ln'=>$mea_ln, 'var_ln'=>$var_ln,
	# 'Q'=>$clsts, 'indv'=>$$indvs, 'nloci'=>$loci}
	#
	# parse_struct(\@files, $in_path, $grp, $out_path)
	#
	# return ($clsts, \@K, \@files, $repss, $ma)
	#
    # new versions, only retain needed information.
	
	my $files   = shift;
	my $in_path = shift;
	my $grp     = shift;
	my $out_path= shift; # Path used to output errors.
	my @files   = @$files;
	my ($indivs, $indivs_p, $loci, @K, $ma, $clsts, $repss, $qtable);
    my ($pops, $order, $N, $popsize);
    my $pop_flag = 1;
    mkdir("$in_path/converted");
	
	if ($grp) {
		if ($grp =~ /^(\d+,)*\d+,?$/g) {
			# newpop size.
			my @parts = split(',', $grp);
            my @order;
			my $line = 1;
			for (my $i=0; $i<@parts; $i++) {
				my $size = $parts[$i];
                push @order, $i+1;
                $popsize->{$i+1} = $size;
				$N += $size;
				for (my $j=$line; $j<=$N; $j++) {
					$pops->{$j} = $i+1;
					$line++; # line number.
				}
			}
			$order = \@order;
		} else {
			($pops, $order, $popsize) = parse_popmap($grp);
			$N = scalar(keys %$pops);
		}
	}
	
	foreach my $file (natsort @files) {
		my ($tclsts, $k, $indv, $LnK, $mea_ln, $var_ln);
		my ($h, $ln_info, $before, @output, $out_fh, $out_fh_grp, @pop_table_q, @ind_table_q);
		my $first = 1;my $pop_table = '';
        my $has_run_par = 0;
        open(my $in_fh, "$in_path/$file") or die "$!";
		if (!$grp&&$first) {open($out_fh_grp, ">$out_path/PopMap_for_admixture") or die "$!";}
		while(<$in_fh>) {
            next if /^$/;
            my $pop_table_cnt = 0;
            if (/Run parameters/) {
                $has_run_par = 1;
            } 
            if (/(\d+)\s+individuals/) {
                $indivs = $1;
            }
            if (/(\d+)\s+loci/) {
                $loci =$1;
            } 
            if (/(\d+)(\s+)populations(\s+)assumed/) {
                $k = $1;
                push @K, $k;
            }        
            if (/Proportion of membership of each pre-defined/) {
                # population tables.
                $pop_table .= "--------------------------------------------\r\n";
                $pop_table .= $_;                
                my $tmp = <$in_fh>;
                while ($tmp !~/----------/) {
                    $pop_table .= $tmp;
                    $pop_table_cnt++;
                    push @pop_table_q, $tmp if $pop_table_cnt > 5; # store in clsts.
                    $tmp = <$in_fh>;
                    die "Too large file or the Format is not correct !!!" if $pop_table_cnt > BIG;
                }
                $pop_table .= "--------------------------------------------\r\n";
            } else {
                # other part before Q-matrix.
                if ($_ !~/Inferred ancestry of individuals/) {
                    $before .= $_;
                }
            }
            
			if (/(.+?)\s+=\s(.+)/) {
			#
			# Estimated Ln Prob of Data   = -9092.6
			# Mean value of ln likelihood = -7201.2
			# Variance of ln likelihood   = 3782.8
			# Mean value of alpha         = 0.0772
			#
				#print join("\t", $1, $2), "\n";
				$LnK    = $2 if $1 eq 'Estimated Ln Prob of Data';
				$mea_ln = $2 if $1 eq 'Mean value of ln likelihood';
				$var_ln = $2 if $1 eq 'Variance of ln likelihood';
                $ln_info .= $_;
			} 
            
            if (/Inferred ancestry of individuals/) {
                #
                # get Q-matrix.
                #

                #
                if (!$has_run_par) {
                    my $info = "Can not get Run parameters for STRUCTURE file: $file, ";
                       $info.= "please check your file.\n";
                    die $info;
                }
                if (!$k) {
                    my $info = "Can not get k in $file, ";
                    $info   .= "please check your file.\n";
                    die $info;
                }
                if (not defined $LnK) {
                    my $info = "Can not get likelihood values for STRUCTURE file: $file, ";
                    $info   .= "please check your file.\n";
                    die $info;
                }
                if ($indivs_p && $indivs_p != $indivs) {
                     die "Number of individuals is not consistent with previous files, in $file.\n";
                }
                if ($grp) {
				
                    if ($N != $indivs) {
                        my $info = "PopMap (group size) error! Number of individuals are not consistent with that in file: $file\n";
                        $info   .= ", in popmap: $N, in file: $indivs\n";
                        $info   .= ", please check your popmap file...\n";
                        die $info;
                    }
                }
                #
                $indivs_p  = $indivs;
                my $header = <$in_fh>;
                unless ($header =~ /^\s+/) {$header = <$in_fh>;} # any line before the true header was skipped.
                my $has_lab= 0;
                my $has_pop= 0;
                $header    = (split(/:/,$header))[0];
                my @parts  = split(/\s+/, $header);
                foreach my $idx (@parts) {
                    $has_pop = 1 if $idx eq 'Pop';
                    $has_lab = 1 if $idx eq 'Label';
                }
                $header = " Label (%Miss) Pop :  Inferred clusters\r\n";
                push @output, ($before, $pop_table, "\r\n".$_.$header); # output for converted files.
                for (my $i=0; $i<$indivs; $i++) {
                    $indv++; # number of individuals passed.
                    my $pop = '';
                    my $f   = <$in_fh>;
                    $f =~ s/[\n\r]//g; # This is important.
                    my $idx = rindex($f, ':');
                    my $id  = substr($f, 0, $idx);
                    my $q   = substr($f, $idx+1);
                    die "Number of individuals is not identical in $file, expected $indivs, but get $indv.\n" if $f =~ /^$/;
                    #
                    # get label and pop info, if not add it.
                    #
                    $id =~ s/^\s+|\s+$//g;
                    my @indId = split(/\s+/, $id);
                    my $seq   = $indId[0];
                    if (!$has_lab) {
                        # no individual labels, add sequential numbers.
                        unshift (@indId, $seq); 
                    } else {
                        $indId[1] = $seq; # replace label by seq number.
                    }
                    
                    if ($grp) {
					
                        $pop = $pops->{$seq};
                        if (!$has_pop) {
                            push (@indId, $pop); # no pop indicator, add popmap.
                        } else {
                            $indId[-1] = $pop;
                        }   
                    } else {
                        if (!$has_pop) {
                            # warning, no pop indicator
                            if ($pop_flag) {
                                my $info = "No population indicator in file: $file, you should turn on the PopData option (check the box 'Putative population origin for each individual') in STRUCTURE!";
                                $info   .= "\nOr upload a popmap file, or input a vector of new pop sizes!";
                                warn $info;
                                $pop_flag = 0;
                            }
                            push (@indId, 1); # no pop indicator, add 1 for all.
                        }
                        $pop = $indId[-1]; # predefined pop.
                        $pops->{$seq} = $pop;
                        if ($first) {print $out_fh_grp "$pop\n";}
                    }
                    
                    $q =~ s/^\s+|\s+$//g;
                    
                    my @qparts = split(/\s+/, $q);
                    my $len   = @qparts;
				
                    if ($k < $len) {
                        # remove extra data after Q-matrix if exists
                        splice (@qparts, $k-$len);
                    } elsif ($k > $len) {
                        die "The number of columns (in Q-matrix) is smaller than K.";
                    }
                    
                    for (my $i=0; $i<@qparts; $i++) {
                        push @{$clsts->{$file}->{$i+1}->{$pop}}, $qparts[$i]; # array of each cluster of each pop.
                        push @{$ma->{$file}->[$i]}, $qparts[$i];
                    }
                    my $tmp = join(" ", @indId) . ' :  ' . join(" ", @qparts) . "\r\n";
                    push @output, $tmp;
                    push @ind_table_q, $tmp;
                }
                
                last;
            }
    	}
        
        if (scalar(@output) == 0) {
            # nothing left in _f.
            my $info = "Can not get any Q-matrix for STRUCTURE file: $file, ";
               $info.= "please check your file.\n";
            die $info;
        }
        
        $first = 0;
		close $in_fh;
        if ($grp) {
            # recalculate population q.
            ($output[1], my $pop_table_q) = calc_popq($clsts->{$file}, $order, $popsize);
            @pop_table_q = @$pop_table_q;
        }
        push @{$qtable->{'popq'}->{$k}}, \@pop_table_q; 
        push @{$qtable->{'indq'}->{$k}}, \@ind_table_q;
        push @{$qtable->{'file'}->{$k}}, $file;
              
		my @f_parts = split(/_/, $file);
		my $run = $f_parts[-2];
		#print join("\t", $run, $var_ln, $LnK, $mea_ln), "\n";
		push @{$repss->{$k}}, {'run'=>$run, 'estLnPk'=>$LnK, 'mea_ln'=>$mea_ln, 'var_ln'=>$var_ln,
		'Q'=>$clsts, 'indv'=>$indivs, 'nloci'=>$loci};
        # output converted file.
        open($out_fh, ">$in_path/converted/$file") or die "$!";
		print $out_fh @output;
        close $out_fh;
	}
	
	@K = uniq @K;
	@K = sort{$a<=>$b} @K;
    
    # 1:popmap. 2: number of files. 
    # 3:number of K. 4: number of individuals. 
    push my @summary, ($pops, \@files, \@K, $indivs);
	
    return ($clsts, \@K, \@files, $repss, $ma, $pop_flag, \@summary, $qtable);	
}

sub calc_pi {
	my ($clsts, $prob) = @_;
	# The parsimony method
	# 1. the replicate runs with Pr[X|K] values smaller than the mean were discarded.
	my ($mea, $tmp) = mean(\values(%$prob));
	my ($nr, $ACSw, $ACSb, $SPS, $SPSp, $SPSs);
	# 2. an assignment quality score is calculated for each of the nr retained replicate runs.
	#    ACSw = 1/nw * sum_i(sum_j(sum_k(Qi*Qj))) within one cluster.
	#    ACSb = 1/nb * sum_j(sum_i(sum_k(Qi*Qj))) between clusters.
	#    SPS = ACSw - ACSb.
    foreach my $run (keys %$clsts) {
		next if $prob->{$run} < $mea;
		$nr++;
		#$clsts->{$run}
	}
	# 3. the harmonic mean of the sizes of well‐defined clusters is calculated.
	#    n5p = ()
	#    HMCS = count(n5p)/sum(1/n5p)

	# 4. the overall strength of the inferred population structure for a given run is calculated.
	#    SPSp = SPS/HMCS
	#    SPSs = max(sps_p) among nr runs.

	# 5. the total number of clusters whose sizes are not larger than 5 (as determined in step 3) across the nr
	# retained runs, n5−, is calculated.
	#    n5m = ()

	# 6. the consistency of ancestry assignments cross runs is calculated.
	#    CAA = 1/(n*(n-1)/2) * sum(score1+score2)
	#    score1 = 1 (m=0)  or 0 (m>0)
	#    score2 = 1 (m=nr) or 0 (m<nr)
	# 7. the overall assignment quality for an assumed number of k populations is measured by the parsimony index
    #    PI = SPSs + CAA - 2*n5m/(k*nr) - 1/(2*k)
	#
}

sub calc_popq {

    my $clsts   = shift;
    my $tmp     = shift;
    my $popsize = shift; # size.
    my @order   = @$tmp;
    my ($pop_table, @a, $popid, $c, @pop_table_q);
    foreach my $k (sort{$a<=>$b}  keys %$clsts) {
        # k from 1 to max.
        push @a, $k;
        foreach my $pop (keys %{$clsts->{$k}}) {
            # pop from 1 to max.
            my @clst = @{$clsts->{$k}->{$pop}}; # array of one cluster for a population.
            my ($mean,$avg) = mean(\@clst);
            push @{$pop_table->{$pop}}, (sprintf "%.3f", $mean); # population is sequential.
        }
    }
    my $ksize   = scalar(@a);
    map{$c .= sprintf "%2s%-5d", '  ', $_;} @a;
    my $b       = $ksize*5 + ($ksize-1)*2;
    my $output .= "--------------------------------------------\r\n";
       $output .= "Proportion of membership of each pre-defined\r\n";
       $output .= " population in each of the ".$ksize." clusters\r\n\r\n";
       $output .= (sprintf "%-9s", 'Given').(sprintf "%-*s", $b, 'Inferred Clusters')."     Number of\r\n";
       $output .= (sprintf "%-9s", ' Pop').$c."  Individuals\r\n\r\n";
    foreach (@order) {
        $popid++;
        my $tmp = (sprintf "%-9s", " $popid:").join("  ", @{$pop_table->{$popid}})."       $popsize->{$_}\r\n";
        push @pop_table_q, $tmp;
        $output.= $tmp;
    }
    $output .= "--------------------------------------------\r\n";
    return ($output, \@pop_table_q);
}

sub rm_ghost {
    #
    # remove ghost clusters
    # rm_ghost($clsts, $method, $thred, $out, $ma)
    # 180911 @Yulong
    # 
    my $clsts  = shift;
    my $method = shift;
    my $thred  = shift;
    my $out    = shift;
    my $ma     = shift; # $ma->{$run}[K-1] = Q;
    my $K      = {};
    my @Ks_new;
    if (!-d $out) {mkdir($out)}
    for my $run (keys %$clsts) {
		
		my (@Q, @ks);
		my $k_old = scalar(keys %{$clsts->{$run}});
		for my $k (keys %{$clsts->{$run}}) {
            # assert each cluster.
			my $t->{'med'} = 0;
			$t->{'avg'} = 0;
			foreach my $pop (keys %{$clsts->{$run}->{$k}}) {
				my @tmp = sort{$b<=>$a} @{$clsts->{$run}->{$k}->{$pop}};
                if ($method eq 'med') {
                    my $med = median(\@tmp);
                    if ($med > $thred) {
                        $t->{'med'}++;
                    }
                    last if $t->{'med'} >= 1;
                }
                
                if ($method eq 'avg') {
                    my ($avg, $sum) = mean(\@tmp);
                    if ($avg > $thred) {
                        $t->{'avg'}++;
                    }
                    last if $t->{'avg'} >= 1;
                }
                
			}
            
            if ($t->{'med'} >= 1 || $t->{'avg'} >= 1) {
                # ghost
                push @ks, $k;
            } 

		}
        
        for my $k (sort {$a<=>$b} @ks){
            --$k;
            my $i = 0;
            foreach (@{$ma->{$run}->[$k]}) {
                $Q[$i] .= $_ . "\t";
                $i++;
            }
        }
        
		#
		# for each reps.
		#
        my $k_new = scalar(@ks); push @Ks_new, $k_new;
        next if $k_new == 0;
        open(my $out_fh, ">$out/$run.$k_old.$k_new.Q") or die "$!";
        print $out_fh join("\n", @Q), "\n";
        close $out_fh;
        $K->{$k_old}->{$k_new}++;
	}
    @Ks_new = sort{$b<=>$a} @Ks_new;
    $K->{'max'} = $Ks_new[0];
    $K->{'min'} = $Ks_new[-1];
    return $K;
}

sub compress_files {
    #
    # compress files according to the original K.
    # compress_files($in_path, $out_path, $thred, $K)
    #
    my $in_path  = shift;
    my $out_path = shift;
    my $thred    = shift;
    my $K        = shift; # what K you want?
    my $ks;
    if (!-d $out_path) {mkdir($out_path)}
    my @files = <$in_path/*.Q>;
    foreach my $file (@files) {
        if ($file =~ /.*\.(\d+)\.(\d+)\.Q$/) {
            my $k_old = $1;
            my $k_new = $2;
            if (defined $K && $K eq 'new') {
                push @{$ks->{$k_new}}, $file;
            } else {
                push @{$ks->{$k_old}}, $file;
            }
        }
    }
    foreach my $k (keys %$ks) {
    
        my $zip = Archive::Zip->new();
        map {my $f=basename($_); $zip->addFile($_, "K=$k/$f");} @{$ks->{$k}};
        $zip->writeToFileNamed("$out_path/$k.$thred.zip");
    }
    my $zip = Archive::Zip->new();
    $zip->addTreeMatching( $in_path, "K=ALL", "\.Q");
    $zip->writeToFileNamed("$out_path/ALL.$thred.zip");
}

1
