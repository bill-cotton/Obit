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
sub Roots;        # Process routines
sub doDoc($);     # Beginning of routine documentation
sub GetArgType($);# Determine routine argument type
sub GetArgDoc($); # Determine routine Documentation
sub Obitize($);   # Determine routine Documentation

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
    
# prototype
    @ProtoLines = 0;
    $NumProto   = 0;

# Swallow file
    SwallowFile;

# Various operations
    Lower;
    Comments;
    Semicolon;
    Substitute;
    Subscripts;
    Loops;
    Roots;

# Write output
    BarfC;
   
# any cleanup operations
  EXIT:
} # end main block

#------------- Subroutines used ------------------------
# Get next argument from $line reading from SRC 
# Puts markers at the beginning of routines
sub SwallowFile {
  local $line, $newRoutine, $routineName, @temp, $number;

# initialize file
  @FileLines = 0;
  $NumLines = 0;
  $number = 0;     # How many routines

printf ("reading file %s\n", $SrcFile); # debug

# open input file
  open (SRC, "<$SrcFile") || 
      die "Cannot open input $SrcFile\n";

# loop over file parsing
  $newRoutine = 0;
  $routineName = "None";
  while ($line = <SRC>)   {
      # is the the beginning of a subroutine?
      if ($line =~ /^\s{6}SUBROUTINE\s{1}/) {
	  $newRoutine = 1;
	  # Label end of old if any
	  if ($number > 0 ) {
	      $FileLines[$NumLines++]= "C *****f2c: end of input routine $routineName\n";
	  }
	  $number++;
	  # extract routine name
	  @temp = split(" ",$line);
	  $routineName = @temp[1];
	  # mark beginning of new routine
	  $FileLines[$NumLines++]= "C *****f2c: beginning of input routine $routineName\n";
      }
      # Insert routine info before first comment in new routine
      if (($newRoutine==1) && ($line =~ /^C/)) {
	  $newRoutine = 0;
	  $FileLines[$NumLines++]= "C Routine translated from the AIPSish $SrcFile/$routineName\n";
	  #print $FileLines[$NumLines++];
      }

      # Make includes into comments
      $line =~ s/\s{6}INCLUDE/C      INCLUDE/;
      
      # If a continuation line, append to previous
      if ($line =~ /^\s{5}\S/) {
	  $line = substr($line,6);  # drop beginning and continuation marker
	  $line =~s/^\s*//;  # trim initial blanks
	  chomp($FileLines[$NumLines-1]); # remove old newline
	  $FileLines[$NumLines-1] = "$FileLines[$NumLines-1] $line";
      } else { # Start new line 
	  $FileLines[$NumLines++]= $line;
      }

} # end loop over file

  # mark end of last routine
  if ($number > 0 ) {
      $FileLines[$NumLines]= "C *****f2c: end of input routine $routineName\n";
      $NumLines++;
  }
  close (SRC); # close file

} # end SwallowFile


# routine to write output file
sub BarfC {
    local $i;
    
printf ("writing file %s\n", $OutFile); # debug

    open (OUT, ">$OutFile") || 
	die "Cannot open $OutFile\n";

# Add prototypes at top
    printf (OUT "/* Prototype for functions in $OutFile */\n");
    for ($i=0; $i<$NumProto; $i++) {
	printf (OUT $ProtoLines[$i]);
    }
    printf (OUT "/* End prototype for functions in $OutFile */\n\n");
    
# loop over file
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
	if (!($FileLines[$i] =~ /^[Cc]/)) { # don't do all-comment lines
	    $FileLines[$i] = lc($FileLines[$i]);
	}
} # end loop over file
    
} # end lower

# routine to convert comments
sub Comments {
    local $i;
    
# loop over file 
    for ($i=0; $i<$NumLines; $i++) {
	if (substr ($FileLines[$i],0,1) eq "C") {
	    chomp($FileLines[$i]); # remove newline
	    $FileLines[$i] = "/*".substr($FileLines[$i],1)." */\n";
	    # trim initial blanks
	    $FileLines[$i] =~s/^\/\*\s{35,80}/\/\* /;
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
    my @parts;
    
# loop over file 
    for ($i=0; $i<$NumLines; $i++) {
	#print $FileLines[$i];
	if (substr ($FileLines[$i],0,2) ne "/*") { # not in comment

	    # Is this a declaration?
	    $decl = ($FileLines[$i] =~ / *olong /) ||
		($FileLines[$i] =~ / *ofloat /) ||
		($FileLines[$i] =~ / *gboolean /) ||
		($FileLines[$i] =~ / *gchar /) ||
		($FileLines[$i] =~ / *odouble /);

	    # special handling for character strings 
	    if ($FileLines[$i] =~ / *gchar/) {
		@subsrpt = ($FileLines[$i] =~ /\w*\*\d*/g);
		$n = @subsrpt;
		if ($n > 0) {
		    for ($j=0; $j<$n; $j++) {
			# get string dimension
			$temp = $subsrpt[$j];
			@parts = split(/\*/,$temp);
			# add one for null
			$temp = @parts[1];
			if (substr ($temp,0,1) =~ /\d/) {$temp = $temp+1;}
			else {$temp = $temp."\+1";}
			# replace
			$FileLines[$i] =~ s/$parts[0]\*$parts[1]/$parts[0]\[$temp\]/;
			#print $FileLines[$i];
		    }
		}
	    } # end special handling for character declarations

	    # singly subscripted
	    # Require nonblank before the ( to separate subscripts 
	    # and call arguments
	    @subsrpt = ($FileLines[$i] =~ /\w\({1}\w{1,20}[\)]/g);

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

	    # Do it again for nested subscripts
	    @subsrpt = ($FileLines[$i] =~ /\w\({1}\w{1,20}[\)]/g);
	    $n = @subsrpt;
	    if ($n > 0) {
		for ($j=0; $j<$n; $j++) {
		    # drop initial character
		    $subsrpt[$j] = substr ($subsrpt[$j] ,1);
		    # replace * with 1
		    $subsrpt[$j] =~ s/\*/\\*/g;
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
		    $subsrpt[$j] =~ s/\(/\\(/g;
		    $subsrpt[$j] =~ s/\)/\\)/g;
		    # escape [ and ]
		    $subsrpt[$j] =~ s/\[/\\[/g;
		    $subsrpt[$j] =~ s/\]/\\]/g;
		    # replace
		    $FileLines[$i] =~ s/$subsrpt[$j]/[@args[0]]/;
		}
	    } # end loop over subscripts


	    # Doubly subscripted
	    @subsrpt = ($FileLines[$i] =~ /\w\({1}\w{1,20},\w{1,20}[\)]/g);
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
	    @subsrpt = ($FileLines[$i] =~ /\w\({1}\w{1,20},\w{1,20},\w{1,20}[\)]/g);
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
	#print $FileLines[$i];
	if (substr ($FileLines[$i],0,2) ne "/*") { # not in comment
	    # change subroutine to void and get name 
	    if ($FileLines[$i] =~ / * subroutine */) {
		$FileLines[$i] =~s/ *subroutine/void/;
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

	    # Comment out "call"
	    $FileLines[$i] =~s/\scall\s/\/\*CALL\*\/ /;

	    # macros
	    $FileLines[$i] =~s/ max *\(/ MAX \(/g;
	    $FileLines[$i] =~s/ min *\(/ MIN \(/g;
	    $FileLines[$i] =~s/ pi / G_PI /g;
	    $FileLines[$i] =~s/dg2rad/DG2RAD/g;
	    $FileLines[$i] =~s/rad2dg/RAD2DG/g;
	    $FileLines[$i] =~s/as2rad/AS2RAD/g;
	    $FileLines[$i] =~s/rad2as/RAD2AS/g;
	    $FileLines[$i] =~s/aipsish/AIPSish/g;
            # quotes
	    $FileLines[$i] =~s/\'/\"/g;
            # numbers - use e rather than d for exponent
	    $FileLines[$i] =~s/(\d)d(\d)/$1e$2/g;
	    $FileLines[$i] =~s/(\d)d-/$1e-/g;

	    # Common functions without a blank before (
	    $FileLines[$i] =~s/([ \+\-\*\/])sin\(/$1sin \(/g;
	    $FileLines[$i] =~s/([ \+\-\*\/])cos\(/$1cos \(/g;
	    $FileLines[$i] =~s/([ \+\-\*\/])tan\(/$1tan \(/g;
	    $FileLines[$i] =~s/([ \+\-\*\/])atan2\(/$1atan2 \(/g;
	    $FileLines[$i] =~s/([ \+\-\*\/])abs\(/$1abs \(/g;
	    $FileLines[$i] =~s/([ \+\-\*\/])sqrt\(/$1sqrt \(/g;
	    $FileLines[$i] =~s/([ \+\-\*\/])MAX\(/$1MAX \(/g;
	    $FileLines[$i] =~s/([ \+\-\*\/])MIN\(/$1MIN \(/g;

	    # declarations
	    $FileLines[$i] =~s/ *integer/  olong/;
	    $FileLines[$i] =~s/ *logical/  gboolean/;
	    $FileLines[$i] =~s/ *real/  ofloat/;
	    $FileLines[$i] =~s/ *double precision/  odouble/;
	    $FileLines[$i] =~s/ *character/  gchar/;
	    
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

	    # deal with GO TO
	    if ($FileLines[$i] =~ / * go to */) {
		$FileLines[$i] =~s/go to/goto/; # respell goto
		$label = substr($FileLines[$i], 4+index ($FileLines[$i], "goto"));
		chomp ($label);     # remove newline
		$label =~ s/;//g;   # remove any semicolons
		$label =~ s/^\s+//; # remove leading blanks
		# Can't do if too complex
		if (!($label =~ /\(/)) {
		    $FileLines[$i] =~ s/$label/L$label/;
		}
		#print $FileLines[$i],$label,"\n";
	    }

	   # change labels
	   if ((!($FileLines[$i]=~/^\/\*/)) && (substr($FileLines[$i],0,5)) =~ /\d/) {
	       $label = substr($FileLines[$i],0,5);
	       $label =~ s/^\s+//; # remove leading blanks
	       $label =~ s/\s+$//; # remove trailing blanks
	       #print $i, $FileLines[$i];
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
    
} # end Loops

# Process routines
sub Roots {
    local @rootLines, $NumRLines, $routine;
    local @newRoot, @newFile, $numNewFile;
    my $i, $j;
    
    # init
    $routine = "None";
    $numNewFile = 0;
# loop over file 
    for ($i=0; $i<$NumLines; $i++) {
	# look for start of routine 
	if ($FileLines[$i] =~ /^\/\* \*\*\*\*\*f2c: beginning of input routine/) { 
	    # initialize routine
	    @rootLines = 0;
	    $NumRLines = 0;
	    $routine = substr($FileLines[$i],40);
	    $routine =~ s/^\s+//;      # remove leading blanks
	    $routine =~ s/\s*\*\/.*//; # remove stuff at end
	    chomp($routine);
	    # end begining of routine
	} elsif ($FileLines[$i] =~ /^\/\* \*\*\*\*\*f2c: end of input routine/) {
	    # process routine
	    doDoc($routine);   # Documentation and interface, writes @newRoot

	    Obitize($routine); # Obit specific transformations, modifies @newRoot

	    # Copy to output
	    #print "$#newRoot line 30 ",$newRoot[30],"\n";
	    for ($j=0; $j<=$#newRoot; $j++) {
		$newFile[$numNewFile++] = $newRoot[$j];
	    }

	    #print $newFile[$numNewFile-1];
	} elsif ($routine ne "None") {    # In routine
	    $rootLines[$NumRLines++] = $FileLines[$i];
	} else {                          # Line before first routine copy to output
	    $newFile[$numNewFile++] = $FileLines[$i];
	}
    } # end loop over file

    # Copy $newFile back to $FileLines
    print " Output files now has $numNewFile lines\n";
    for ($i=0; $i<$numNewFile; $i++) {
	$FileLines[$i] = $newFile[$i];
    }
    $NumLines = $numNewFile;
    
} # end Roots

# Header documentation
# argument =  routine name, 
# Uses externals $NumRLines, @rootLines
# out output puts full routine in @newRoot
# function prototype modifies externals $NumProto, @ProtoLines
sub  doDoc($) {
    my $i, $name, @routine, $nLines, $temp, @rootArg, $out;
    my @argType, $numType, @argDoc, $numDoc;
    
    # Local copies of arguments
    $name    = $_[0];
    $nLines  = $NumRLines;

    # Get list of arguments
    @rootArg = 0;
    $temp = $rootLines[0];
    $temp =~ s/.*\(//;  # trim beginning
    $temp =~ s/\).*//s;  # trim end
    @rootArg = split(", ", $temp);
    #print $name, " ",$#rootArg," ",$nLines,"\n";

    # Get types of arguments
    @argType = 0;
    $numType = 0;
    for ($i=0; $i<=$#rootArg; $i++) {
	@argType[$numType++] = GetArgType($rootArg[$i]);
	#print "now $i $name $rootArg[$i] @argType[$numType-1]\n";
    }

    # add type to first line
    for ($i=0; $i<=$#rootArg; $i++) {
	$rootLines[0] =~ s/([\( ,])$rootArg[$i]/$1@argType[$i] $rootArg[$i]/;
    }

    # Add to prototypes
    @ProtoLines[$NumProto] = $rootLines[0];
    @ProtoLines[$NumProto] =~ s/\s*\n\{/;/;
    $NumProto++;

    # Get Descriptions of arguments
    @argDoc = 0;
    $numDoc = 0;
    for ($i=0; $i<=$#rootArg; $i++) {
	@argDoc[$numDoc++] = GetArgDoc($rootArg[$i]);
	#print "how $i $name $rootArg[$i] @argDoc[$numDoc-1]";
    }

    # Doxygen header
    $out = 0;
    @newRoot = 0;
    $newRoot[$out++] = "/**\n";

    # Copy non-copyright comments until first non comment
    $i = 0;
    while (($rootLines[++$i] =~ /^\/\*/) && !($rootLines[$i] =~/\s*inputs\:/i)) {
	if (!($rootLines[$i] =~ /^\/\*[ -][;-]/)) { # Copyright?
	    $temp = substr ($rootLines[$i], 2); # drop beginning comment delimit
	    chomp($temp);                       # remove newline
	    $temp =~ s/^\s*//;                  # trim initial whitespace
	    $temp =~ s/\*\///g;                 # drop end comment delimit
	    $newRoot[$out] = " * $temp \n";
	    $out++;
	    #print $newRoot[$out];
	}
    }

    for ($i=0; $i<=$#rootArg; $i++) {
	$newRoot[$out++] = " * \\param $rootArg[$i]   @argDoc[$i]\n";
	$out++;
    }
    $newRoot[$out++] = " */\n";

    #copy rest of routine 
    for ($i=0; $i<$nLines; $i++) {
	$newRoot[$out++] = $rootLines[$i];
    }
} # end doDoc

# Determine routine argument type
# Uses externals $NumRLines, @rootLines
# Also modify routine for output scalars (add '*')
# Returns type as e.g. "ofloat*" or "olong" or "gchar[21]"
sub GetArgType($) {
    my $out, $i, $arg, isOut;

    # Local copies of arguments
    $arg    = $_[0];
 
    $out = "gpointer";

    #print "arg = $arg \n";
    # loop over file looking for $arg in declaration
    for ($i=0; $i<$NumRLines; $i++) {
	if ($rootLines[$i] =~ /[ ,]$arg[\[ ,\(\;]/) {
	    if ($rootLines[$i] =~ /\s*Obit/) { # Obit variable
		$out = $rootLines[$i];
		chomp ($out);
		#$out =~ /(\s*)(Obit)(\w*)(.*)/Obit$3\*/;
	    } elsif ($rootLines[$i] =~ /\s*olong /) { # integer
		if ($rootLines[$i] =~ /[ ,]$arg[\[\(]/) {
		    $out = "olong*";
		} else {
		    $out = "olong";
		}
		break;
	    } elsif ($rootLines[$i] =~ /\s*gchar /) { # string
		if ($rootLines[$i] =~ /[ ,]$arg[\[\(]/) {
		    $out = "gchar*";
		} else {
		    $out = "gchar";
		}
		break;
	    } elsif ($rootLines[$i] =~ /\s*ofloat /) { # float
		if ($rootLines[$i] =~ /[ ,]$arg[\[\(]/) {
		    $out = "ofloat*";
		} else {
		    $out = "ofloat";
		}
		break;
	    } elsif ($rootLines[$i] =~ /\s*gboolean /) { # boolean
		if ($rootLines[$i] =~ /[ ,]$arg[\[\(]/) {
		    $out = "gboolean*";
		} else {
		    $out = "gboolean";
		}
		break;
	    } elsif ($rootLines[$i] =~ /\s*odouble /) { # double
		if ($rootLines[$i] =~ /[ ,]$arg[\[\(]/) {
		    $out = "odouble*";
		} else {
		    $out = "odouble";
		}
		break;
	    }
	    #print "found $rootLines[$i]\n";
	}
    } # end loop looking for declarations

    # If it's already a pointer type we're done
    if ($out =~ /\*/) {return($out);}

    $isOut = 0;
    # loop over file looking for $arg in output
    for ($i=0; $i<$NumRLines; $i++) {
	if (!($rootLines[$i] =~/^\/\*/ ) && ($rootLines[$i] =~ /\s* $arg.*=/)) {
	    $isOut = 1;
	    $out =~ s/(\w*)/$1\*/;
	    break;
	}
    }
    if ($isOut==0) {return($out);} # done ?

    # Add '*' to variable name
    # loop over file
    for ($i=0; $i<$NumRLines; $i++) {
	if (!($rootLines[$i] =~/^\/\*/ ) && ($rootLines[$i] =~ /\s*[ \+\*\-\/]$arg[; \+\*\-\/]/)) {
	    $rootLines[$i] =~ s/(\s*[ \+\*\-\/])($arg)([; \+\*\-\/])/$1\(\*$arg\)$3/g;
	}
    }
 
    return ($out);
} # end GetArgType

# Determine routine argument documentation string
# Uses externals $NumRLines, @rootLines
# Also modify routine for output scalars (add '*')
# Returns type as e.g. "ofloat*" or "olong" or "gchar[21]"
sub GetArgDoc($) {
    my $out, $i, $arg, $ip, $temp;

    # Local copies of arguments
    $arg    = $_[0];
 
    $out = "Unavailable\n";

    #print "arg = $arg \n";
    # loop over file looking for $arg in header comments
    $i = 0;
    while ($rootLines[++$i] =~ /^\/\*/) {  # Until first non comment
	if ($rootLines[$i] =~ /^\/\*\s*$arg[\(\[ ]/i) { # $arg first non blank
	    $out = $rootLines[$i];
	    $out =~ s/(^\/\*\s*$arg[\(\[ ])(.*)(\*\/)/$2/i;
	    $out =~ s/^\s*[IRDC].{3}//;   # remove type/size 
	    chomp($out);                  # remove newline
	    $out =~ s/^\s*\)//;
	    # Copy any other
	    $ip = index(lc($rootLines[$i]), $arg);
	    #print "$ip $rootLines[$i+1]\n";
	    while ($rootLines[++$i] =~ /^\/\*\s{$ip}/)  {
		$temp = substr($rootLines[$i], $ip);
		$temp =~ s/\*\///;
		chomp($temp);
		$out = $out."\n * $temp";
		#print "$out \n";
 	    }
	    break;
	}
   } # end loop looking for documentation

    return ($out);
} # end GetArgDoc

# Obit specific transformations
# Uses externals @newRoot, $NumProto, @ProtoLines
sub Obitize($) {
    my $routine, $nRoot, $i, $doErr, $inExec, $temp;

    # Local copies of arguments
    $routine  = $_[0];
    $nRoot = $#newRoot + 1;   # no. entries in @newRoot

    # Replace ierr with ObitErr?
    $doErr = (@ProtoLines[$NumProto-1] =~ /olong[\* ]*ierr/);
    if ($doErr) {
	@ProtoLines[$NumProto-1] =~ s/olong[\* ]*ierr/ObitErr\* err/;
    } else {     # nothing else for now
	return;
    }

    # loop over routine
    $inExec = 0;  # Gotten into executable code?
    for ($i=0; $i<$nRoot; $i++) {

	# modify call, ierr => err
	$newRoot[$i] =~ s/olong[\* ]*ierr/ObitErr\* err/;
	$newRoot[$i] =~ s/(\W)(ierr)(\W)/$1err$3/;

	# Startup at boundry between declarations and executable
	if ((!$inExec) && (@newRoot[$i] =~ /^\s{3,6}/)) {
	    $inExec = 1;
	    @newRoot[$i-1] = @newRoot[$i-1]."  gchar *routine = \"$routine\";\n\n";
	    if ($doErr) {
		@newRoot[$i-1] = @newRoot[$i-1]."  /* Error checks */\n".
		    "  g_assert(ObitErrIsA(err));\n".
		    "  if (err->error) return ;  /* previous error? */\n\n";
	    }
	}

	# test on error condition
	$newRoot[$i] =~ s/(\Werr)(\s*)(!=)/$1->error$2$3/;

	# update changing error condition
	$newRoot[$i] =~ s/\(\*err\)\s*=/err->error =/;

	# AIPS DMSG.INC -> local declaration
	$newRoot[$i] =~ s/\/\*      INCLUDE 'INCS:DMSG.INC' \*\//  gchar msgtxt[81];/;

	# string assignment of msgtxt into c
	$newRoot[$i] =~ s/(^.*)(\smsgtxt\s*)(=)(\s*)(".*")/$1strncpy \(msgtxt,$5,80\)/;

	# partial modification of write into msgtxt
	$newRoot[$i] =~ s/(^.*)(write\s*\(msgtxt,)(\s*\d*)(\))(.*)(;)/$1g_snprintf \(msgtxt,80,"format $3",$5\);/;

	# replace MSGWRT with message
	if ($newRoot[$i] =~ /\/\*CALL\*\/ msgwrt \([01234]\)/) {
	    $newRoot[$i] = "    Obit_log_error(err, OBIT_InfoErr, msgtxt);\n";
	}
	# replace MSGWRT with warning message
	if ($newRoot[$i] =~ /\/\*CALL\*\/ msgwrt \([5]\)/) {
	    $newRoot[$i] = "    Obit_log_error(err, OBIT_InfoWarn, msgtxt);\n";
	}
	# replace MSGWRT with error message
	if ($newRoot[$i] =~ /\/\*CALL\*\/ msgwrt \([6789]\)/) {
	    $newRoot[$i] = "    Obit_log_error(err, OBIT_Error, msgtxt);\n";
	}

    } # end loop over routine
} # end Obitize
