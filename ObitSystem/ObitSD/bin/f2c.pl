# This script parses a fortran file and converts (sort of ) to c.
# The intended translation from from AIPSish fortran to OBIT c
# with considerable hand work later.  This is merely an asistance tool.
# First input argument is the name of the input file
# Second input argument is the name of the output file.
#

require 5;

# mention subroutines
sub SwallowFile;  # Read input file and convert to an array of strings
sub BarfC;        # Write output
sub Lower;        # All to lower case
sub Semicolon;    # Adds semicolons
sub Comments;     # convert comments
sub Subscripts;   # convert subscripts
sub Substitute;   # common substitutions
sub Loops;        # do loops

# begin main block
{

# set default input file
    $SrcFile = "AIPS.FOR";

# use first argument if it exists as input file name
#printf ("number of args $#ARGV $ARGV[0]\n");
    if ($#ARGV>=0) {$SrcFile=$ARGV[0]};
    
# Set output file name
    $OutFile =  "Obit.c";
    if ($#ARGV>=1) {$OutFile=$ARGV[1]};
    
# initialize file
    @FileLines = 0;
    $NumLines = 0;
    
# Swallow file
    SwallowFile;

# Various operations
    Lower;
    Comments;
    Semicolon;
    Subscripts;
    Substitute;
    Loops;

# Write output
    BarfC;
   
# any cleanup operations
  EXIT:
} # end main block

#------------- Subroutines used ------------------------
# Get next argument from $line reading from SRC if necessary
sub SwallowFile {
  local $line;

# initialize file
  @FileLines = 0;
  $NumLines = 0;

printf ("reading file %s\n", $SrcFile); # debug

# open input file
  open (SRC, "<$SrcFile") || 
      die "Cannot open input $SrcFile\n";

# loop over file parsing
  while ($line = <SRC>)   {
      $FileLines[$NumLines]= $line;
      $NumLines++;
#     printf ("$line");

} # end loop over file

  close (SRC); # close file

} # end SwallowFile


# routine to write output file
sub BarfC {
    local $i;
    
printf ("writing file %s\n", $OutFile); # debug

    open (OUT, ">$OutFile") || 
	die "Cannot open $OutFile\n";
    
# loop over file parsing
    for ($i=0; $i<$NumLines; $i++) {
	printf (OUT $FileLines[$i]);
    } # end loop over file
    
    close (OUT); # close file
} # end BarfC

# routine to convert all to lower case
sub Lower {
    local $i;
    
# loop over file 
    for ($i=0; $i<$NumLines; $i++) {
	$FileLines[$i] = lc($FileLines[$i]);
} # end loop over file
    
} # end lower

# routine to convert comments
sub Comments {
    local $i;
    
# loop over file 
    for ($i=0; $i<$NumLines; $i++) {
	if (substr ($FileLines[$i],0,1) eq "c") {
	    chomp($FileLines[$i]); # remove newline
	    $FileLines[$i] = "/*".substr($FileLines[$i],1)." */\n";
	    # trim initial blanks
	    $FileLines[$i] =~s/\/\* */\/\* /;
	}
    } # end loop over file
    
} # end Comments

# routine to add semicolons
sub Semicolon {
    local $i;
    
    $/ = "\n";
# loop over file 
    for ($i=0; $i<$NumLines; $i++) {
	if (substr ($FileLines[$i],0,2) ne "/*") { # not in comment
	    if ((substr ($FileLines[$i+1],5,1) eq " ") ||
		(substr ($FileLines[$i+1],0,2) eq "/*")) {
		chomp($FileLines[$i]); # remove newline
		$FileLines[$i] = $FileLines[$i].";\n";
	    }
# remove continuation character 
	    if ((substr ($FileLines[$i],5,1) ne " ") &&
		(substr ($FileLines[$i],0,2) ne "/*")) {
		$FileLines[$i] = substr($FileLines[$i],6);
	    }
	} # end not a comment
    } # end loop over file
    
} # end Semicolon

# routine to convert subscripts
sub Subscripts {
    local $i, $j, $k, $n, $m, @subsrpt, $temp, @args, $decl;
    
# loop over file 
    for ($i=0; $i<$NumLines; $i++) {
	if (substr ($FileLines[$i],0,2) ne "/*") { # not in comment

	    # Is this a declaration?
	    $decl = ($FileLines[$i] =~ /integer/) ||
		($FileLines[$i] =~ /real/) ||
		($FileLines[$i] =~ /logical/) ||
		($FileLines[$i] =~ /double precision/);

	    # singly subscripted
	    # Require nonblank before the ( to separate subscripts 
	    # and call arguments
	    @subsrpt = ($FileLines[$i] =~ /\w\({1}\w{1,5}[\)]/g);
	    $n = @subsrpt;
	    if ($n > 0) {
		for ($j=0; $j<$n; $j++) {
		    # drop initial character
		    $subsrpt[$j] = substr ($subsrpt[$j] ,1);
		    # remove ( and )
		    $temp = $subsrpt[$j];
		    $temp =~ s/\(//;
		    $temp =~ s/\)//;
		    # reformat 
		    @args = split(/,/, $temp);
		    # subtract one if not declaration
		    if (!$decl) {
			$m = @args;
			for ($k=0; $k<$m; $k++) {
			    if (substr (@args[$k],0,1) =~ /\d/) {@args[$k] = @args[$k]-1;}
			    else {@args[$k] = @args[$k]."-1";}
			}
		    }
		    # escape ( and )
		    $subsrpt[$j] =~ s/\(/\\(/;
		    $subsrpt[$j] =~ s/\)/\\)/;
		    # replace
		    $FileLines[$i] =~ s/$subsrpt[$j]/[@args[0]]/;
		}
	    } # end loop over subscripts

	    # Doubly subscripted
	    @subsrpt = ($FileLines[$i] =~ /\w\({1}\w{1,5},\w{1,5}[\)]/g);
	    $n = @subsrpt;
	    if ($n > 0) {
		for ($j=0; $j<$n; $j++) {
		    # drop initial character
		    $subsrpt[$j] = substr ($subsrpt[$j] ,1);
		    # remove ( and )
		    $temp = $subsrpt[$j];
		    $temp =~ s/\(//;
		    $temp =~ s/\)//;
		    # reformat 
		    @args = split(/,/, $temp);
		    # subtract one if not declaration
		    if ($decl == null) {
			$m = @args;
			for ($k=0; $k<$m; $k++) {
			    if (substr (@args[$k],0,1) =~ /\d/) {@args[$k] = @args[$k]-1;}
			    else {@args[$k] = @args[$k]."-1";}
				}
		    }
		    # escape ( and )
		    $subsrpt[$j] =~ s/\(/\\(/;
		    $subsrpt[$j] =~ s/\)/\\)/;
		    # replace
		    $FileLines[$i] =~ s/$subsrpt[$j]/[@args[1]][@args[0]]/;
		}
	    } # end loop over subscripts

	    # triply subscripted
	    @subsrpt = ($FileLines[$i] =~ /\w\({1}\w{1,5},\w{1,5},\w{1,5}[\)]/g);
	    $n = @subsrpt;
	    if ($n > 0) {
		for ($j=0; $j<$n; $j++) {
		    # drop initial character
		    $subsrpt[$j] = substr ($subsrpt[$j] ,1);
		    # remove ( and )
		    $temp = $subsrpt[$j];
		    $temp =~ s/\(//;
		    $temp =~ s/\)//;
		    # reformat 
		    @args = split(/,/, $temp);
		    # subtract one if not declaration
		    if (!$decl) {
			$m = @args;
			for ($k=0; $k<$m; $k++) {
			    if (substr (@args[$k],0,1) =~ /\d/) {@args[$k] = @args[$k]-1;}
			    else {@args[$k] = @args[$k]."-1";}
			}
		    }
		    # escape ( and )
		    $subsrpt[$j] =~ s/\(/\\(/;
		    $subsrpt[$j] =~ s/\)/\\)/;
		    # replace
		    $FileLines[$i] =~ s/$subsrpt[$j]/[@args[1]][@args[0]]/;
		}
	    } # end loop over subscripts

	} # end not a comment
    } # end loop over file
    
} # end subscripts

# Do common substitions
sub Substitute {
    local $i, $ind, $routine, @parts, $temp, $label, $start;
    
    $start = 0;  # haven't started a routine
    $/ = "\n";
# loop over file 
    for ($i=0; $i<$NumLines; $i++) {
	if (substr ($FileLines[$i],0,2) ne "/*") { # not in comment
	    # change subroutine to void and get name 
	    if ($FileLines[$i] =~ / * subroutine */) {
		$FileLines[$i] =~s/ *subroutine/  void/;
		$temp = substr($FileLines[$i], 4+index ($FileLines[$i], "void"));
		$temp =~ s/^\s+//; # remove leading blanks
		@parts = split(/ /,$temp);
		$routine = $parts[0];
		chomp ($routine);     # remove newline
		$routine =~ s/;//g;   # remove any semicolons
		$routine =~ s/\s+$//; # remove trailing blanks
		$start = 1;
	    }

	    # replace semicolon at end of routine call arguments with {
	    if ($start && $FileLines[$i] =~ /;/) {
		$start = 0;
		$FileLines[$i] =~ s/;/ \n\{/; 
	    }

	    # macros
	    $FileLines[$i] =~s/ max *\(/ MAX \(/;
	    $FileLines[$i] =~s/ min *\(/ MIN \(/;

	    # declarations
	    $FileLines[$i] =~s/ *integer/  gint/;
	    $FileLines[$i] =~s/ *logical/  gboolean/;
	    $FileLines[$i] =~s/ *real/  gfloat/;
	    $FileLines[$i] =~s/ *double precision/  gdouble/;
	    
	    # if/then/else/end if
	    $FileLines[$i] =~s/ then *;*/ {/;
	    $FileLines[$i] =~s/ else *;*/  } else {/;
	    $FileLines[$i] =~s/ end if *;*/  } /;
	    $FileLines[$i] =~s/ end *;*/  } \/\* end of routine $routine \*\/ /;

	    # comparison operators ,logicals
	    $FileLines[$i] =~s/\.eq\./ == /g;
	    $FileLines[$i] =~s/\.ne\./ != /g;
	    $FileLines[$i] =~s/\.gt\./ > /g;
	    $FileLines[$i] =~s/\.ge\./ >= /g;
	    $FileLines[$i] =~s/\.lt\./ < /g;
	    $FileLines[$i] =~s/\.le\./ <= /g;
	    $FileLines[$i] =~s/\.not\./!/g;
	    $FileLines[$i] =~s/\.and\./ && /g;
	    $FileLines[$i] =~s/\.or\./ || /g;
	    $FileLines[$i] =~s/\.true\./TRUE/g;
	    $FileLines[$i] =~s/\.false\./FALSE/g;

	    # deal with go to
	    if ($FileLines[$i] =~ / * go to */) {
		$FileLines[$i] =~s/go to/goto/; # respell goto
		$label = substr($FileLines[$i], 4+index ($FileLines[$i], "goto"));
		chomp ($label);     # remove newline
		$label =~ s/;//g;   # remove any semicolons
		$label =~ s/^\s+//; # remove leading blanks
		# Can't do if too complex
		if (!($FileLines[$i] =~ /\(/)) {
		    $FileLines[$i] =~ s/$label/L$label:/;
		}
	    }

	   # change labels
	   if (substr($FileLines[$i],0,5) =~ /\d/) {
	       $label = substr($FileLines[$i],0,5);
	       $label =~ s/^\s+//; # remove leading blanks
	       $label =~ s/\s+$//; # remove trailing blanks
	       $FileLines[$i] =~ s/$label/L$label:/;
	   }

	} # end not a comment
    } # end loop over file
    
} # end Substitute


# Do loops => for loops
sub Loops {
    local $i, $label, @parts, @parts2, $var, $lower, $upper, $temp;
    
    $/ = "\n";
# loop over file 
    for ($i=0; $i<$NumLines-1; $i++) {
	if (substr ($FileLines[$i],0,2) ne "/*") { # not in comment
	    # do
	    if ($FileLines[$i] =~ / +do +/) {
		$temp = $FileLines[$i];
		chomp($temp);      # remove newline
		$temp =~ s/;//g;   # remove any semicolons
		$temp =~ s/^\s+//; # remove leading blanks
		@parts = split(/=/, $temp);  # parse line
		@parts2 = split(/ /, $parts[0]);
		$label = $parts2[1]; 
		$var   = $parts2[2];
		@parts2 = split(/,/, $parts[1]);
		$lower = $parts2[0];
		$upper = $parts2[1];
		#printf ("for ($var=$lower; $var<=$upper; $var++) { /* loop $label */\n");
		$FileLines[$i] = 
		    "      for ($var=$lower; $var<=$upper; $var++) { /* loop $label */\n";
	    }

	    # continue -> end } 
	    if ($FileLines[$i] =~ / *continue */) {
		$label = substr ($FileLines[$i],0,6);
		$FileLines[$i] = "      ".substr ($FileLines[$i],6,80);
		$FileLines[$i] =~s/continue/} \/\* end loop $label \*\//;
	    }
	} # end not a comment
    } # end loop over file
    
} # end Substitute

