#!C:/Strawberry/perl/bin/perl.exe
#Rna-Rna Interaction Using Genetic Algorithm
use warnings;
use strict;
use Text::Balanced;
main();

sub main{
	## input : 2 Rna sequences
	## output : Predicted the intermolecular and interamolecular bonding with minimum free energy
	open (my $in,"RnaRnaSeq.txt") or die "Can't read $!";  
	my $RnaRna = <$in>;
	$RnaRna = uc($RnaRna);
	chomp($RnaRna);
	my $Rna1 = $RnaRna;
	my $Rna2 = $RnaRna;
	$Rna1 =~ s/&.*//;
	$Rna2 =~ s/.*&//;
	
	#  Predict secondary structures of RNA1 and RNA2 
	my $strRna1 = Predict_secondary_structure($Rna1);
	my $strRna2 = Predict_secondary_structure($Rna2);
	
	#  optimal and suboptimal solution of predict secondary structures of RNA1 and RNA2 
	my @subOptRna1 = Optimal_secondary_structure($Rna1);
	my @subOptRna2 = Optimal_secondary_structure($Rna2);
	
	# Bit matrix to generate all possible intermolecular stems between RNA1 and RNA2 
	my @RnaRnaStems = RnaRnaStems($Rna1,$Rna2);
	#print_Stems_inFile(\@RnaRnaStems,$Rna1,$Rna2);
	
	my $populationSize = 1;
	my @solutions = ();
	for(my $i=0;$i<$populationSize;++$i){
		my $rndRna1 = $subOptRna1[rand @subOptRna1];
		my $rndRna2 = $subOptRna2[rand @subOptRna2];
		my %solution = (rna1=>$rndRna1, rna2=>$rndRna2, stems=>[] );
		Mutation(\%solution,\@RnaRnaStems,$Rna1,$Rna2);
		push @solutions, {%solution};
		
	}
}

# Predict Secondary Structure
sub Predict_secondary_structure{
	my($RNA)=@_;
	my($STR);	
	open(filehandelO3,'>',"secondary_structure.txt");
	print(filehandelO3	$RNA,"\n","@\n");
	close(filehandelO3);
	$STR= qx/RNAfold.exe <secondary_structure.txt/;
	$STR=~ s/.*\n//;
	$STR=~ s/ \(.*\)\n//;
	return($STR);
}
# End Predict Secondary Structure

# Calculate Energy Solution
sub CompleteEnergy{
	my($SP,$RNA)=@_;
	my($prog);
	my($ES);	
	open(filehandelO3,'>',"solutionstem.txt");
	print(filehandelO3	$RNA,"\n",$SP,"\n","@\n");
	close(filehandelO3);
	$ES= qx/RNAeval.exe <solutionstem.txt/;
	#print $ES;
	#<STDIN>;
	$ES=~ s/.*\n//;
	$ES=~ s/.* //;
	$ES=~ s/\(//;
	$ES=~ s/\)//;
	return($ES);
	
}
# End Calculate Energy Solution

# Predict Optimal and Suboptimal Secondary Structure
sub Optimal_secondary_structure{
	my($RNA)=@_;
	print $RNA,"\n\n";
	my($STR);	
	open(filehandelO3,'>',"optimal_secondary_structure.txt");
	print(filehandelO3	$RNA,"\n","@\n");
	close(filehandelO3);
	$STR= qx/RNAsubopt.exe <optimal_secondary_structure.txt/;
	$STR=~ s/.*\n//;
	my @Sub_optimal = split /\n/, $STR;
	for(my $i=0;$i<scalar(@Sub_optimal);$i++){
		$Sub_optimal[$i] =~ s/ .*//;
		#print $Sub_optimal[$i],"\n";
	}
	return(@Sub_optimal);
}
# End Predict Optimal and Suboptimal Secondary Structure

# Create Dot Plot
sub DotPlot {
	my($RNA1,$RNA2)=@_;
	my @seq1 = split('',$RNA1);
	my @seq2 = split('',$RNA2);
	my @matrix = ();
	for(my $i=0;$i<scalar(@seq1);$i++)
	{
		for(my $j=0;$j<scalar(@seq2);$j++)
		{
			if(($seq1[$i] eq "G" and $seq2[$j] eq "C") or ($seq1[$i]eq"A" and $seq2[$j]eq"U") or ($seq1[$i]eq"C" and $seq2[$j]eq"G") or ($seq1[$i]eq"U" and $seq2[$j]eq"A"))
			{
				$matrix[$i][$j]="//";
			}
			else
			{
				$matrix[$i][$j]=" ";
			}
		}
	}
	return @matrix;
}
#End Create Dot Plot

# Declare Rna1-Rna2 Stems
sub RnaRnaStems {
	my($RNA1,$RNA2)=@_;
	my @matrix = DotPlot($RNA1,$RNA2);
	my @seq1 = split('',$RNA1);
	my @seq2 = split('',$RNA2);
	my @stems = ();
	# first declare stems with length 3
	# u1 : initial nucletide position
	# u2 : final nucletide position
	# u3 : length of the stem
	for(my $i=0;$i<scalar(@seq1);++$i){
		for(my $j=scalar(@seq2)-1;$j>=0;--$j){
			if($matrix[$i][$j] eq "//"){
				my $start=$i;
				my $end=$j;
				my $k;
				my $min = scalar(@seq1);
				if(scalar(@seq2)<scalar(@seq1)){
					$min = scalar(@seq2);
				}
				for($k=1;$k<$min;++$k){
					if(($i+$k)<scalar(@seq1) and ($j-$k)>=0){
						if($matrix[$i+$k][$j-$k] eq " "){
						
							if(abs($i+$k-1-($j-$k+1)+1) < 3){
								$k=0;
							}
							last;
						}
					}
					else{
						last;
					}
				}
				my $length=$k;
				if($length>=3){
					push @stems, {u1=>$start, u2=>$end, u3=>$length };
				}
			}
		}
	}
	my @sorted =  sort { $a->{u3} <=> $b->{u3} } @stems;
	@stems = @sorted;
	return @stems;
}
#End Declare Stems

# Print rna1-rna2 stems in file
sub print_Stems_inFile {
	my @stems = @{$_[0]};
	my $RNA1 = $_[1];
	my $RNA2 = $_[2];
	my @seq1 = split('',$RNA1);
	my @seq2 = split('',$RNA2);
	open(my $out,'>',"D:/Stems.txt") or die "Can't open file for writing: $!"; 
	for(my $i=0;$i<scalar(@stems);$i++)
	{
		print $out $i," => ";
		print $out $stems[$i]{u1}," ",$stems[$i]{u2}," ",$stems[$i]{u3}," ";
		for(my $j=0;$j<$stems[$i]{u3};$j++){
			print $out $seq1[$stems[$i]{u1}+$j];
		}
		print $out " ";
		for(my $j=$stems[$i]{u3}-1;$j>=0;$j--){
			print $out $seq2[$stems[$i]{u2}-$j];
		}
		print $out "\n";
	}
	close $out or die "Failed to close file: $!";
}
#End Print rna1-rna2 stems in file

# Mutation
sub Mutation {
	my %new_solution = %{$_[0]};
	my @stems = @{$_[1]};
	# Rna1
	my $seqRNA1 = $_[2];
	my $strRNA1 = $new_solution{rna1};
	print $strRNA1,"\n";
	my @strRna1 = split('',$strRNA1);
	
	# Rna 2
	my $seqRNA2 = $_[3];
	my $strRNA2 = $new_solution{rna2};
	print $strRNA2,"\n";
	my @strRna2 = split('',$strRNA2);
	
	#compute original energy
	my $original_energy = CompleteEnergy($strRNA1,$seqRNA1) + CompleteEnergy($strRNA2,$seqRNA2);
	
	# Random pick a stem with length≥ 4 and energy < 0 
	my %rnd_stem =  %{$stems[rand @stems]};
	my $StrRnaRna = structure_Rna_Rna(\%rnd_stem,$seqRNA1,$seqRNA2);
	my $rnd_stem_energy = RnaRnaEnergy($StrRnaRna,$seqRNA1,$seqRNA2);
	#my $rnd_stem_energy = calculate_Stem_Energy(\%rnd_stem,$seqRNA1,$seqRNA2);
	while($rnd_stem{u3}<4 or $rnd_stem_energy>=0){
		%rnd_stem = %{$stems[rand @stems]};
		$StrRnaRna = structure_Rna_Rna(\%rnd_stem,$seqRNA1,$seqRNA2);
		$rnd_stem_energy = RnaRnaEnergy($StrRnaRna,$seqRNA1,$seqRNA2);
		#$rnd_stem_energy = calculate_Stem_Energy(\%rnd_stem,$seqRNA1,$seqRNA2);
	}
	print $rnd_stem{u1},",",$rnd_stem{u2},",",$rnd_stem{u3},"\n";
	print $rnd_stem_energy;
	
	# Insert the picked stem to the stem list of Snew 
	my @stem_list = @{$new_solution{stems}};
	push @stem_list,\%rnd_stem;
	
	# change RNA1 and RNA2 structure because of the collision between internal bonding and intermolecular bonding
	my @new_strRna1 = Remove_parentheses_rna1(\@strRna1,\%rnd_stem);
	my @new_strRna2 = Remove_parentheses_rna2(\@strRna2,\%rnd_stem);

	my $new_strRNA1 = join('',@new_strRna1);
	my $new_strRNA2 = join('',@new_strRna2);
	print $new_strRNA1,"\n";
	print $new_strRNA2,"\n";
	
	#  Interaction Length = number of interacting nucleotides + loop sizes
	my  $Interaction_Length = Interaction_Length(\%rnd_stem,\@new_strRna1,\@new_strRna2);
	print $Interaction_Length,"\n";
	
	#while($Interaction_Length < 65){
		my $new_energy = CompleteEnergy($new_strRNA1,$seqRNA1) + CompleteEnergy($new_strRNA2,$seqRNA2) + Energy_Model($rnd_stem_energy,$Interaction_Length);
		print "origin : ",$original_energy,"\n";
		print "new :",$new_energy,"\n";
		if( $new_energy > $original_energy ){
			$new_strRNA1 = $strRNA1;
			$new_strRNA2 = $strRNA2;
			@new_strRna1 = split('',$strRNA1);
			@new_strRna2 = split('',$strRNA2);
			print $new_strRNA1,"\n";
			print $new_strRNA2,"\n";
		}
		# randomly expand to the left and right by adding one nearby stem until the maximum interaction length is reached.
		
	#}
}
# End Mutation

# Stack Energy

sub stack_energy {
	my %stack_energy;
	$stack_energy{"AA"}{"UU"} = -0.9;
	$stack_energy{"AC"}{"UG"} = -2.2;
	$stack_energy{"AG"}{"UC"} = -2.1;
	$stack_energy{"AG"}{"UU"} = -0.6;
	$stack_energy{"AU"}{"UA"} = -1.1;
	$stack_energy{"AU"}{"UG"} = -1.4;
	$stack_energy{"CA"}{"GU"} = -2.1;
	$stack_energy{"CC"}{"GG"} = -3.3;
	$stack_energy{"CG"}{"GC"} = -2.4;
	$stack_energy{"CG"}{"GU"} = -1.4;
	$stack_energy{"CU"}{"GA"} = -2.1;
	$stack_energy{"CU"}{"GG"} = -2.1;
	$stack_energy{"GA"}{"CU"} = -2.4;
	$stack_energy{"GA"}{"UU"} = -1.3;
	$stack_energy{"GC"}{"CG"} = -3.4;
	$stack_energy{"GC"}{"UG"} = -2.5;
	$stack_energy{"GG"}{"CC"} = -3.3;
	$stack_energy{"GG"}{"CU"} = -1.5;
	$stack_energy{"GG"}{"UC"} = -2.1;
	$stack_energy{"GG"}{"UU"} = -0.5;
	$stack_energy{"GU"}{"CA"} = -2.2;
	$stack_energy{"GU"}{"CG"} = -2.5;
	$stack_energy{"GU"}{"UA"} = -1.4;
	$stack_energy{"GU"}{"UG"} = 1.3;
	$stack_energy{"UA"}{"AU"} = -1.3;
	$stack_energy{"UA"}{"GU"} = -1;
	$stack_energy{"UC"}{"AG"} = -2.4;
	$stack_energy{"UC"}{"GG"} = -1.5;
	$stack_energy{"UG"}{"AC"} = -2.1;
	$stack_energy{"UG"}{"AU"} = -1;
	$stack_energy{"UG"}{"GC"} = -1.4;
	$stack_energy{"UG"}{"GU"} = 0.3;
	$stack_energy{"UU"}{"AA"} = -0.9;
	$stack_energy{"UU"}{"AG"} = -1.3;
	$stack_energy{"UU"}{"GA"} = -0.6;
	$stack_energy{"UU"}{"GG"} = -0.5;
	return %stack_energy;
}

# End Stack Energy

# Calculate stem energy

sub calculate_Stem_Energy {
	my %stem = %{$_[0]};
	my $RNA1 = $_[1];
	my $RNA2 = $_[2];
	my @seq1 = split('',$RNA1);
	my @seq2 = split('',$RNA2);
	my %stack_energy = stack_energy();
	my $sum_energy = 0;
	for(my $j=0;$j<$stem{u3}-1;++$j){
		foreach my $top (sort { $a cmp $b} keys %stack_energy){
			foreach my $bottom ( keys %{$stack_energy{$top}}){
				my $pair1 = join('',$seq1[$stem{u1}+$j],$seq1[$stem{u1}+$j+1]);
				my $pair2 = join('',$seq2[$stem{u2}-$j],$seq2[$stem{u2}-$j-1]);
				if($pair1 eq $top and $pair2 eq $bottom)
				{
					$sum_energy = $sum_energy + $stack_energy{$top}{$bottom};
				}
			}
		}	
	}
	return $sum_energy;
}

# End Calculate stem energy

# Remove Parantheses rna1 according to rna1-rna2 stem

sub Remove_parentheses_rna1 {
	my @str = @{$_[0]};
	my %stem = %{$_[1]};
	my @remove = ();
	for(my $i=0;$i<$stem{u3};$i++){
		if($str[$stem{u1}+$i] eq "\("){
			$str[$stem{u1}+$i] = "\(@";
			push @remove,$stem{u1}+$i;
		}
		if($str[$stem{u1}+$i] eq "\)"){
			$str[$stem{u1}+$i] = "\)@";
			push @remove,$stem{u1}+$i;
		}
		
	}
	my @temp = ();
	for(my $i=0;$i<scalar(@str);$i++){
		if($str[$i] eq "\("){
			push @temp,$i;
		}
		if($str[$i] eq "\(@"){
			push @temp,($i."@");
		}
		if($str[$i] eq "\)"){
			if($temp[-1] =~ m/(\d@)/ ){
				$temp[-1] =~ s/@//g;
				push @remove,$i;
			}
			pop @temp;
		}
		if($str[$i] eq "\)@"){
			push @remove,$temp[-1];
			pop @temp;
		}
	}
	for(my $i=0;$i<scalar(@remove);$i++){
		$str[$remove[$i]] = ".";
	}
	return @str;
}
# End Remove Parantheses rna1 according to rna1-rna2 stem

# Remove Parantheses rna2 according to rna1-rna2 stem
sub Remove_parentheses_rna2 {
	my @str = @{$_[0]};
	my %stem = %{$_[1]};
	my @remove = ();
	for(my $i=0;$i<$stem{u3};$i++){
		if($str[$stem{u2}-$i] eq "\("){
			$str[$stem{u2}-$i] = "\(@";
			push @remove,$stem{u2}-$i;
		}
		if($str[$stem{u2}-$i] eq "\)"){
			$str[$stem{u2}-$i] = "\)@";
			push @remove,$stem{u2}-$i;
		}
		
	}
	my @temp = ();
	for(my $i=0;$i<scalar(@str);$i++){
		if($str[$i] eq "\("){
			push @temp,$i;
		}
		if($str[$i] eq "\(@"){
			push @temp,($i."@");
		}
		if($str[$i] eq "\)"){
			if($temp[-1] =~ m/(\d@)/ ){
				$temp[-1] =~ s/@//g;
				push @remove,$i;
			}
			pop @temp;
		}
		if($str[$i] eq "\)@"){
			push @remove,$temp[-1];
			pop @temp;
		}
	}
	for(my $i=0;$i<scalar(@remove);$i++){
		$str[$remove[$i]] = ".";
	}
	return @str;
}
# End Remove Parantheses rna2 according to rna1-rna2 stem

# Calculate Energy Solution
sub RnaRnaEnergy{
	my($SP,$RNA1,$RNA2)=@_;
	my($ES);	
	open(filehandelO3,'>',"solutionstem.txt");
	print(filehandelO3	$RNA1,"&",$RNA2,"\n",$SP,"\n","@\n");
	close(filehandelO3);
	$ES= qx/RNAeval.exe <solutionstem.txt/;
	#print $ES;
	#<STDIN>;
	$ES=~ s/.*\n//;
	$ES=~ s/.* //;
	$ES=~ s/\(//;
	$ES=~ s/\)//;
	return($ES);
	
}
# End Calculate Energy Solution

sub structure_Rna_Rna {
	my %rnd_stem = %{$_[0]};
	my $seqRNA1 = $_[1];
	my $seqRNA2 = $_[2];
	my @strRnaRna = ();
	for(my $i=0;$i<length($seqRNA1);++$i){
		push @strRnaRna,".";
	}
	for(my $i=0;$i<$rnd_stem{u3};$i++){
		$strRnaRna[$rnd_stem{u1}+$i] = "\(";
	}
	push @strRnaRna,"&";
	for(my $i=length($seqRNA1)+1;$i<length($seqRNA1)+1+length($seqRNA2);++$i){
		push @strRnaRna,".";
	}
	for(my $i=0;$i<$rnd_stem{u3};$i++){
		$strRnaRna[length($seqRNA1)+1+$rnd_stem{u2}-$i] = "\)";
	}
	return join('',@strRnaRna);
}

sub Interaction_Length {
	my %rnd_stem = %{$_[0]};
	my @strRna1 = @{$_[1]};
	my @strRna2 = @{$_[2]};
	my $k = 0;
	my $start_loop1 = $rnd_stem{u1}-$k;
	my $end_loop1 = $rnd_stem{u1}+$k;
	while($strRna1[$rnd_stem{u1}-$k] eq "."){
		$start_loop1 = $rnd_stem{u1}-$k;
		if($rnd_stem{u1}-$k+1>0){
			$k++;
		}
		else{
			last;
		}
	}
	$k = 0;
	while($strRna1[$rnd_stem{u1}+$k] eq "."){
		$end_loop1 = $rnd_stem{u1}+$k;
		if($rnd_stem{u1}+$k+1<scalar(@strRna1)){
			$k++;
		}
		else{
			last;
		}
	}
	$k = 0;
	my $start_loop2 = $rnd_stem{u2}-$k;
	my $end_loop2 = $rnd_stem{u2}+$k;
	while($strRna2[$rnd_stem{u2}-$k] eq "."){
		$start_loop2 = $rnd_stem{u2}-$k;
		if($rnd_stem{u2}-$k+1>0){
			$k++;
		}
		else{
			last;
		}
	}
	$k = 0;
	while($strRna2[$rnd_stem{u2}+$k] eq "."){
		$end_loop2 = $rnd_stem{u2}+$k;
		if($rnd_stem{u2}+$k+1<scalar(@strRna2)){
			$k++;
		}
		else{
			last;
		}
	}

	return  $rnd_stem{u3} + $end_loop1-$start_loop1+1 + $end_loop2-$start_loop2+1 ;
}

sub Energy_Model {
	my($rnd_stem_energy,$Interaction_Length) = @_;
	my $r = $rnd_stem_energy/$Interaction_Length;
	my $p0=0;
	my $p1=0;
	my $p2=0;
	my $p3=0;
	my $p4=0;
	if($Interaction_Length>=20){
		$p0 = -12;
	}
	if($Interaction_Length>=12){
		$p0 = -12;
	}
	if($Interaction_Length>=-3){
		$p0 = 8;
	}
	if($Interaction_Length<=4){
		$p0 = 12;
	}
	if($Interaction_Length<=8){
		$p0 = 12;
	}
	
	return $r*$rnd_stem_energy + $p0 + $p1 + $p2 + $p3 + $p4;
}