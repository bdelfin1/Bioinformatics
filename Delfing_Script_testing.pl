#!/usr/bin/perl -w
use strict;
my %pdbFiles; #stores raw file names
my %moleculeNames; #stores PDB IDs of proteins, taken from the file names
my %close; #atoms in active site
my @filecontents; #array of pdb file contents
my $line;
my @atom;
my $atom;
my @het;
my $het; 
my $dist; #distance from ligand
my $filename = '/Users/Hubble_Space_Telescope/Desktop/somePymolStuff/testing3.pml'; # file to write to
#currently it uses only one file and rewrites it every time. 
#I wasn't sure if this would be better or worse than generating a random filename unique to each run,
#but it makes it a lot easier for testing purposes, instead of flooded my computer with these files.
for (my $i=0; $i <= 5; $i++) {
	$pdbFiles{$i} = '';
}

my $count = 0;
my $temp = '';
while ($count <= 5 && $temp ne 'stop') {
	print "Enter the location of a .pdb file (up to 6, 1 at a time) or type \'stop\':  ";
	chomp($temp = <STDIN>);
	if ($temp ne 'stop'){
		$pdbFiles{$count} = $temp; # grab pdb file from user
		#$pdbFiles{$count} =~ s/\s+$//; # trim right side of name, basically chomp, but I think it has been more consistent
		if ($pdbFiles{$count} =~ /\/([\d\w]{4})\./){ #grabs 4 char between '/' and '.', aka the /(2hck).pdb bit
			$moleculeNames{$count} = $1;
		}
		$count++;
	}
}

open(my $fh, '>', $filename) or die "Could not open file '$filename' $!"; 
my $fetchLine = "fetch ";

foreach my $key ( keys %moleculeNames ){ #list the molecule names to be fetched
	$fetchLine .= $moleculeNames{$key};
	$fetchLine .= " ";
}

print $fh "$fetchLine, async=0\n hi ev\n sh cartoon\n bg white \n"; #begin pymol script
#A few lines of the beginning, hides everything but cartoon; white bg is nice for screenshots
print $fh "set cartoon_side_chain_helper, on\n";
print $fh "set stick_transparency, .7\n"; #makes side chains 70% transparent
#ligand transparency adjusted separately later

my $ligColor; #determines color for ligand, might end up going with default color instead.

foreach my $key ( keys %moleculeNames ){
	#print "key: $key, value: $moleculeNames{$key}\n";
	open (INPUT, "< $pdbFiles{$key}") || die "Could not open file\n";
	print "looking in  $pdbFiles{$key}\n";
	chomp(@filecontents = <INPUT>);
	foreach $line (@filecontents) {
        if (($line =~ /^ATOM/) && ($line =~ /\bCA/)) {
       # if ($line =~ /^ATOM/)  {
                push(@atom, $line);
        } elsif (($line =~ /^HETATM/) && !($line =~ /HOH/)) { #excluding water seems to have improved results, less clutter
                push(@het, $line);
        }
		}
	my $count = 0; #used for labeling selections
	foreach $atom (@atom) {
        my $x = substr($atom, 31, 8);
        my $y = substr($atom, 39, 8);
        my $z = substr($atom, 47, 8);
        my $resname = substr($atom, 17, 3);
        my $chainid = substr($atom, 21, 1);
        my $resnum  = substr($atom, 22, 4);

        foreach $het (@het) {
                my $a = substr($het, 31, 8);
                my $b = substr($het, 39, 8);
                my $c = substr($het, 47, 8);
        			my $hetresname = substr($het, 17, 3);
       			my $hetchainid = substr($het, 21, 1);
     			   	my $hetresnum  = substr($het, 22, 4);
        	
                my $dist = sqrt( ($x - $a)**2 + ($y - $b)**2 + ($z - $c)**2 );
                if ($dist <= 5) {
        				  $close{"$hetresname:$hetchainid:$resname:$chainid:$resnum"}  = $dist;
        	   	
        	   	      	  $resnum =~ s/^\s+//;
        	            print $fh "sele sele_$count, $chainid/$resname`$resnum/ and $moleculeNames{$key}\n";
        	         	  print $fh "sele sele_$count, sele_$count expand 5\n\n"; # 5 was $dist
        	           	   # select the predicted active site (and the ligand for now, recolored later)
        	         	   $count++;
        	   
              	}
                			
        }
	}

	my $i = 0;
	
	print $fh "sele binding_$moleculeNames{$key}, (";
	while ($i <  $count) { #making binding site selection out of the temp ones, to reduce clutter
			print $fh "sele_$i, ";
			$i++;
	}
	print $fh ")\n";

	$i = 0;
	print $fh "delete (";  #getting rid of the temporary selections
	while ($i <  $count) {
			print $fh "sele_$i, ";
			$i++;
	}
	print $fh ")\n";
	 #No idea why, but sele_0 keeps sticking around.
	#the previous loop SHOULD delete it, and DOES list it to be deleted, but it doesn't actually happen.
	
	print $fh "sele nonbind_$moleculeNames{$key}, !binding_$moleculeNames{$key} and $moleculeNames{$key}\n"; # make a new selection with everything not in predicted active site

	print $fh "color gray40, nonbind_$moleculeNames{$key}\n"; # recolor things not in the selection
	print $fh "desele\n";
	print $fh "sele ligand_$moleculeNames{$key}, binding_$moleculeNames{$key} and hetatm and $moleculeNames{$key}\n"; #selects the hetatms from the selection colored in foreach loop
	print $fh "sele ligand_$moleculeNames{$key}, ligand_$moleculeNames{$key} expand 2\n";
	print $fh "sh sticks, ligand_$moleculeNames{$key} and $moleculeNames{$key}\n"; #show ligand as sticks
	print $fh "set_bond stick_transparency, 0, ligand_$moleculeNames{$key}\n"; #This turns off transparency for the ligands
	print $fh "sh sticks, binding_$moleculeNames{$key} and $moleculeNames{$key}\n"; #show sticks for side chains
	
	$ligColor = 8 * ($key + 1);	#gives color value 8-48 
	print $fh "color $ligColor$ligColor, ligand_$moleculeNames{$key}\n";
	#assigns values 88, 1616, 2424, etc to ligands to color them uniquely
	 #Messing with ligand colors, might use default color instead. 
	print $fh "desele\n";	
	print $fh "hi ev, nonbind_$moleculeNames{$key}\n"; #hiding nonbinding by default, since transparency is weird.
	#possible way to make new objects from selection: 
	#create <object name>, <selection>

}

foreach my $key ( keys %moleculeNames ){
	if ($key != 0) {
		print $fh "align binding_$moleculeNames{$key}, binding_$moleculeNames{0}, object=align_$moleculeNames{$key}\n";
		#aligns each molecule with the first molecule entered		
		#the alignment object created here can be used to view the sequence alignment of that molecule in the sequence viewer.
		print $fh "hi cgo, align_$moleculeNames{$key}\n"; 
	} #fixed problems with align, cealign had worse output.
}
print $fh "set seq_view, 1\n"; #turns on the sequence viewer
print $fh "orient binding_$moleculeNames{0}\n";
print $fh "zoom binding_$moleculeNames{0}, 15 \n"; #trying to come up with a good starting viewpoint; this seems okay. maybe room to improve
system('open /Users/Hubble_Space_Telescope/Desktop/somePymolStuff/testing3.pml');	 #runs pymol script.
close $fh; 		
					
					

