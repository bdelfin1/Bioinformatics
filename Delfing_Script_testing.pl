#!/usr/bin/perl -w
use strict;

my @filecontents;
my $line;
my @atom;
my $atom;
my @het;
my $het;
my $dist;
my $name;
my $distance;
print "Enter the location of a .pdb file:  ";
my $pdbFile = <STDIN>; # grab pdb file from user
$pdbFile =~ s/\s+$//; # trim right side of name, basically chomp, but I think it has been more consistent
my $filename = '/Users/Hubble_Space_Telescope/Desktop/somePymolStuff/testing3.pml'; # file to write to
my %close;
print "Enter another location of a .pdb file:  ";
my @filecontents2;
my $line2;
my @atom2;
my $atom2;
my @het2;
my $het2;
my $dist2;
my $name2;
my $distance2;
my %close2;
my $pdbFile2 = <STDIN>; # grab pdb file from user
$pdbFile2 =~ s/\s+$//; # trim right side of name, basically chomp, but I think it has been more consistent


#
# Example PDB file format for atom type (PERL numbering start w/ 0)
#
#          1         2         3         4         5         6         7         8
#012345678901234567890123456789012345678901234567890123456789012345678901234567890
#
#ATOM    145  N   VAL A  25      32.433  16.336  57.540  1.00 11.92      A1   N
#ATOM    146  CA  VAL A  25      31.132  16.439  58.160  1.00 11.85      A1   C
#HETATM15786  U   IUM A1102      18.227  57.101 113.666  0.30 28.22           U
#
#Record Format
#
#COLUMNS      DATA TYPE        FIELD      DEFINITION
#------------------------------------------------------
# 1 -  6      Record name      "ATOM    "
# 7 - 11      Integer          serial     Atom serial number.
#13 - 16      Atom             name       Atom name.
#17           Character        altLoc     Alternate location indicator.
#18 - 20      Residue name     resName    Residue name.
#22           Character        chainID    Chain identifier.
#23 - 26      Integer          resSeq     Residue sequence number.
#27           AChar            iCode      Code for insertion of residues.
#31 - 38      Real(8.3)        x          Orthogonal coordinates for X in Angstroms
#39 - 46      Real(8.3)        y          Orthogonal coordinates for Y in Angstroms
#47 - 54      Real(8.3)        z          Orthogonal coordinates for Z in Angstroms
#55 - 60      Real(6.2)        occupancy  Occupancy.
#61 - 66      Real(6.2)        tempFactor Temperature factor.
#77 - 78      LString(2)       element    Element symbol, right-justified.
#79 - 80      LString(2)       charge     Charge on the atom.


open (INPUT, "< $pdbFile") || die "Could not open file\n";
chomp(@filecontents = <INPUT>);
#my $moleculeName = substr $pdbFile, -8, 4; #method 1, grabs 8 char from end, until 4 from end (.pdb)
my $moleculeName;
if ($pdbFile =~ /\/([\d\w]{4})\./){ #method 2, grabs 4 char between '/' and '.', aka the /(2hck).pdb bit
	$moleculeName = $1;
} #I think the binding_$moleculeName way works fine if you don't chomp filename, but I'm not sure if there's a reason not to do this way
my $moleculeName2;
if ($pdbFile2 =~ /\/([\d\w]{4})\./){ #method 2, grabs 4 char between '/' and '.', aka the /(2hck).pdb bit
	$moleculeName2 = $1;
} #I think the binding_$moleculeName way works fine if you don't chomp filename, but I'm not sure if there's a reason not to do this way
printf "$moleculeName $moleculeName2\n"; #just prints the molecule name in command line; not essential

open(my $fh, '>', $filename) or die "Could not open file '$filename' $!";
print $fh "fetch $moleculeName $moleculeName2, async=0\n hi ev\n sh cartoon\n bg white \n"; #begin pymol script
#A few lines of the beginning, hides everything but cartoon; white bg is nice for screenshots

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
        	            print $fh "sele sele_$count, $chainid/$resname`$resnum/\n";
        	         	   print $fh "sele sele_$count, sele_$count expand 5\n\n"; # 5 was $dist
        	           	   # select the predicted active site (and the ligand for now, recolored later)
        	         	   $count++;
        	   
              	}
                			
        }
}

#close INPUT;
#######second time; probably make this a function later
open (INPUT2, "< $pdbFile2") || die "Could not open file 2\n";
chomp(@filecontents2 = <INPUT2>);
#my $moleculeName = substr $pdbFile, -8, 4; #method 1, grabs 8 char from end, until 4 from end (.pdb)
#my $moleculeName;
#if ($pdbFile2 =~ /\/([\d\w]{4})\./){ #method 2, grabs 4 char between '/' and '.', aka the /(2hck).pdb bit
#	$moleculeName2 = $2;
#} #I think the binding_$moleculeName way works fine if you don't chomp filename, but I'm not sure if there's a reason not to do this way
#printf "$moleculeName\n"; #just prints the molecule name in command line; not essential

foreach $line (@filecontents2) {
        if (($line =~ /^ATOM/) && ($line =~ /\bCA/)) {
       # if ($line =~ /^ATOM/)  {
                push(@atom2, $line);
        } elsif (($line =~ /^HETATM/) && !($line =~ /HOH/)) { #excluding water seems to have improved results, less clutter
                push(@het2, $line);
        }
}

my $count2 = 0; #used for labeling selections
foreach $atom2 (@atom2) {
        my $x2 = substr($atom2, 31, 8);
        my $y2 = substr($atom2, 39, 8);
        my $z2 = substr($atom2, 47, 8);
        my $resname2 = substr($atom2, 17, 3);
        my $chainid2 = substr($atom2, 21, 1);
        my $resnum2  = substr($atom2, 22, 4);

        foreach $het2 (@het2) {
                my $a2 = substr($het2, 31, 8);
                my $b2 = substr($het2, 39, 8);
                my $c2 = substr($het2, 47, 8);
        			my $hetresname2 = substr($het2, 17, 3);
       			my $hetchainid2 = substr($het2, 21, 1);
     			   	my $hetresnum2  = substr($het2, 22, 4);
        	
                my $dist2 = sqrt( ($x2 - $a2)**2 + ($y2 - $b2)**2 + ($z2 - $c2)**2 );
                if ($dist2 <= 5) {
        				  $close2{"$hetresname2:$hetchainid2:$resname2:$chainid2:$resnum2"}  = $dist2;
        	   	
        	   	      	  $resnum2 =~ s/^\s+//;
        	            print $fh "sele sele_B_$count2, $chainid2/$resname2`$resnum2/\n";
        	         	   print $fh "sele sele_B_$count2, sele_B_$count2 expand 5\n\n"; # 5 was $dist
        	           	   # select the predicted active site (and the ligand for now, recolored later)
        	         	   $count2++;
        	   
              	}
                			
        }
}




my $i = 0;
print $fh "set cartoon_side_chain_helper, on\n";
print $fh "sele binding_$moleculeName, ("; 
while ($i <  $count) {
			print $fh "sele_$i, ";
			$i++;
}
print $fh ")\n";

$i = 0;
print $fh "delete ("; 
while ($i <  $count) {
			print $fh "sele_$i, ";
			$i++;
}
print $fh ")\n";


print $fh "sele nonbind_$moleculeName, !binding_$moleculeName and $moleculeName\n"; # make a new selection with everything not in predicted active site



print $fh "color gray40, nonbind_$moleculeName\n"; # recolor things not in the selection
print $fh "desele\n";
print $fh "sele ligand_$moleculeName, binding_$moleculeName and hetatm and $moleculeName\n"; #selects the hetatms from the selection colored in foreach loop
print $fh "sele ligand_$moleculeName, ligand_$moleculeName expand 2\n";
print $fh "sh sticks, binding_$moleculeName and $moleculeName\n";
print $fh "color red, ligand_$moleculeName\n"; #red seems to be a good color to make it stick out well
print $fh "set transparency, 0.7, nonbind_$moleculeName\n";
print $fh "desele\n";
#########second molecule
print $fh "sele binding_$moleculeName2, ("; 
$i = 0;
while ($i <  $count2) {
			print $fh "sele_B_$i, ";
			$i++;
}
print $fh ")\n";

$i = 0;
print $fh "delete ("; 
while ($i <  $count2) {
			print $fh "sele_B_$i, ";
			$i++;
}
print $fh ")\n";
print $fh "sele nonbind_$moleculeName2, !binding_$moleculeName2 and $moleculeName2 \n";
print $fh "color black, nonbind_$moleculeName2\n"; # recolor things not in the selection

print $fh "desele\n";
print $fh "sele ligand_$moleculeName2, binding_$moleculeName2 and hetatm and $moleculeName2\n"; #selects the hetatms from the selection colored in foreach loop
print $fh "sele ligand_$moleculeName2, ligand_$moleculeName2 expand 2\n";
print $fh "sh sticks, binding_$moleculeName2 and $moleculeName2\n";
print $fh "color red, ligand_$moleculeName2\n"; #red seems to be a good color to make it stick out well
print $fh "set transparency, 0.7, nonbind_$moleculeName2\n"; #idk why this isn't working
print $fh "desele\n";
print $fh "cealign binding_$moleculeName, binding_$moleculeName2\n";

print $fh "orient binding_$moleculeName\n";
print $fh "zoom binding_$moleculeName, 15 \n"; #trying to come up with a good starting viewpoint; this seems okay. maybe room to improve
system('open /Users/Hubble_Space_Telescope/Desktop/somePymolStuff/testing3.pml');	 #runs pymol script.
close $fh;

	
#while (($name, $distance) = each %close) {
#printf "%s %5.3f\n",$name, $distance;

#				}

########test
#my $moleculeName3 = '2hck';
#my $testnum = '3';
#my $testname = $moleculeName.$testnum;
#print "This is a test: "."$testname"."\n";
####end test					
					
					

