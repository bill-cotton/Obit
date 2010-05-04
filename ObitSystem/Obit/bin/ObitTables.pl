# This script parses the LaTeX documentation file and produces
# header and c source files for Obit Tables.
# Output files for Table type "XX" are:
# ObitTableXX.c           - source code for class
# ObitTableXX.h           - Class definition header
# ObitTableXXRowDef.h     - Class Row structure definition header
# ObitTableXXRowClassDef.h- Class Row Class structure definition header
# ObitTableXXDef.h        - Class structure header
# ObitTableXXClassDef.h   - ClassInfo structure header
# TableXX.inc             - Python class specific interface (in python dir)
#
# Keys on the following macroes:
# tablename: Argument is table type, if this is 'XX', class ObitTableXX
#            will be generated.  If XX is two characters this is an AIPS XX table, 
#            otherwise FITS only.
# tablecolindex: This defines zero or more symbols to be used in keyword names
#                (first arg which should be lower case)  and the label for
#                the relevant data column (second argument).  The 1-rel
#                column number of the specified column is determined and this value
#                is substituted for the first argument in any keyword names.
# tablekey:  One per application keyword, arguments are
#             1: name
#             2: type 'I'=integer, 'E'=real, 'D'=double, 'A'=character
#             3: Variable name
#             4: Default value if any
#             5: Range if structural keyword
#             6: Description string
# 
# tablecol:  One per table column, arguments are
#             1: name
#             2: units
#             3: datatype (character)
#                I=16 bit integer, 
#                J=32 bit integer, 
#                E=real, 
#                D=double, 
#                C=Complex (32 bit each), 
#                M=double complex (64 bit each), 
#                A=character, 
#                L=logical, 
#                X=bit array, 
#                B=byte array, 
#             4: dimension string (1-5 dimensions)
#             5: variable name
#             6: description
#
# NOTES:
#  For some characters the standard escape does not appear to work, in these cases,
#  resort to using the ascii values:
#  '\' => \x5c
#  '%' => "%%" , \x25 works no better than '%'

require 5;

# mention subroutines
sub ParseTableDoc;     # Parse the next table description
sub WriteClass;        # write class files
sub GetNextArg;        # parse next argument
sub GetTypeCode;       # Get Obit type code
sub DclTypeCode;       # Get variable declaration type code
sub RowType;           # Get row data name of correct type
sub GetColumnCount;    # Get OBIT column element count from dimension string
sub WriteClassDef;     # Write ClassInfo definition file
sub WriteRowDef;       # Write Class Row structure definition
sub WriteRowClassDef;  # Write Class Row structure definition
sub WriteDef;          # Write Class Structure definition file
sub WriteHeader;       # Write Class definition file
sub WriteSource;       # Write Class function definition file
sub WritePy;           # Write Python Class interface file
sub Copyright;         # Write Copyright
sub NoQuote;           # Return a string without quotes or blanks

# begin main block
{

# initialize descriptor information
$Extname = "Uninitialized"; # empty string or "AIPS "
$ExtType = "";

# set default documentation file
$DocFile = "doc/OBITdoc.tex";

# use first argument if it exists as input file name
#printf ("number of args $#ARGV $ARGV[0]\n");
if ($#ARGV>=0) {$DocFile=$ARGV[0]};

#printf ("reading file %s\n", $DocFile); # debug

# open input file
if (!open (SRC, "<$DocFile")) {
    printf("Cannot open input tex file %s\n", $DocFile);
    goto EXIT;
}

# Parse document file 
ParseTableDoc;

# Write any remaining class definition files
WriteClass;

# any cleanup operations
EXIT:
} # end main block

#------------- Subroutines used ------------------------
# Parse the next table description
sub ParseTableDoc  {
    local @args;

# loop over file parsing
    while ($line = <SRC>)   {
	# Start new cycle on "\tablename{" in line
	# '\' is \x5c which is the only way it can be said
	if ($line=~"\\x5ctablename{") {
	    WriteClass; # write old and clear
	    #print ($line);
	    # get table type
	    $Extname = GetNextArgument();
	    # GODDAMN VLBA assholes
	    $ExtnameClean = $Extname;
	    $ExtnameClean =~ s/-/_/;
            # Is this AIPS (2 char) or only FITS, (History special)?
            if ((length($Extname)==2) || ($Extname eq "History")) {
		$ExtType = "AIPS ";
	    } else {
		$ExtType = "";
	    }
	} # end start new definition
	
	# Swallow table introduction
	elsif ($line=~"\\x5ctableintro\\x5b") {
	    $done = 0;
	    while (($line = <SRC>) && !$done)  {
		$TableIntro[$NumIntro] =  $line;
		# Remove } and {
		$TableIntro[$NumIntro] =~  s/\{//;
		$TableIntro[$NumIntro] =~  s/\}//;
		$NumIntro++;
		$done = ($line =~ /\}/);
	    }
	} # end new table introduction
	
	# Any keyword column indices
	elsif ($line=~"\\x5ctablecolindex\\x5b") {
	    #print ($line);
	    $IndexSymb[$NumIndex] =  GetNextArgument();
	    #print ($IndexSymb[$NumIndex],"\n");
	    $IndexCol[$NumIndex] =  GetNextArgument();
	    #print ($IndexCol[$NumIndex],"\n");
	    $NumIndex++;
	} # end new keyword entry
	
	# Table keyword entry
	elsif ($line=~"\\x5ctablekey\\x5b") {
	    $KeyWord[$NumKey] =  GetNextArgument();
	    $KeyWord[$NumKey] =~ s/ +\"$/\"/; # remove trailing blanks
	    $KeyType[$NumKey] =  GetNextArgument();
	    # use int for integers
	    if ($KeyType[$NumKey] eq "I") {$KeyType[$NumKey] = "K"};
	    $KeyName[$NumKey] =  GetNextArgument();
	    $KeyDefV[$NumKey] =  GetNextArgument();
	    $KeyRang[$NumKey] =  GetNextArgument();
	    $KeyComm[$NumKey] =  GetNextArgument();
	    $KeyNumCol[$NumKey] =  0;
	    $NumKey++;
	} # end new keyword entry
	
	# Table column entry
	elsif ($line=~"\\x5ctablecol\\x5b") {
	    $ColName[$NumCol] =  GetNextArgument();
	    # NO $ColName[$NumCol] =~ s/ +\"$/\"/; # remove trailing blanks
	    $ColUnit[$NumCol] =  GetNextArgument();
	    $ColUnit[$NumCol] =~ s/ +\"$/\"/; # remove trailing blanks
	    $ColType[$NumCol] =  GetNextArgument();
	    $ColDim[$NumCol]  =  GetNextArgument();
	    $ColVName[$NumCol]=  GetNextArgument();
	    $ColComm[$NumCol] =  GetNextArgument();
	    # Split Col Name into root and suffix parts
	    @args = split(/\#/, $ColName[$NumCol]);
	    $ColName[$NumCol] = @args[0];
	    $ColSuff[$NumCol] = NoQuote(@args[1]);
	    $NumCol++;
	} # end new keyword entry
    } # end loop over file
	
    
} # end ParseTableDoc 

# Get next argument from $line reading from SRC if necessary
sub GetNextArgument {
  local $lvalue, $start, $stop;
  $lvalue = $line;
# trim line to next '{'
  $start = index ($line, "{");
  # See if need to try next line
  while ($start<0) {
    $line = <SRC>;
    $start = index ($line, "{");
  }
  $start = $start+1;
  $stop = index ($line, "}");
  $lvalue = substr ($line, $start, $stop-$start);
  # remove any '\'
  $lvalue =~ s/\x5c//g;
  #remove leading and trailing blanks
  $lvalue =~ s/^ +//;
  $lvalue =~ s/ +$//;
  # remove from $line
  $line = substr ($line, $stop+1);
  return $lvalue;
} # end GetNextArgument

# Return a string without quotes or blanks
sub NoQuote {
  local $temp = $_[0];
  $temp =~ s/\"//g;
  $temp =~ s/ //g;
  return $temp;
} # end NoQuote

#Get row data name of correct type
sub RowType {
  local $ltype;
  local $larg = $_[0];
  if ($larg eq 'I') {$ltype = 'siRow';}
  elsif ($larg eq 'J') {$ltype = 'iRow';}
  elsif ($larg eq 'K') {$ltype = 'iRow_oint';}
  elsif ($larg eq 'E') {$ltype = 'fRow';}
  elsif ($larg eq 'D') {$ltype = 'dRow';}
  elsif ($larg eq 'A') {$ltype = 'cRow';}
  elsif ($larg eq 'L') {$ltype = 'lRow';}
  elsif ($larg eq 'C') {$ltype = 'cxRow';}
  elsif ($larg eq 'M') {$ltype = 'dcxRow';}
  elsif ($larg eq 'X') {$ltype = 'iRow';}
  elsif ($larg eq 'B') {$ltype = 'bRow';}
  else {$ltype = 'Unknown type';} # probably should be fixed in documentation
  return $ltype;
} # end RowType

#Get Obit type code
sub GetTypeCode {
  local $ltype;
  local $larg = $_[0];
  if ($larg eq 'I') {$ltype = 'OBIT_short';}
  elsif ($larg eq 'J') {$ltype = 'OBIT_oint';}
  elsif ($larg eq 'K') {$ltype = 'OBIT_oint';}
  elsif ($larg eq 'E') {$ltype = 'OBIT_float';}
  elsif ($larg eq 'D') {$ltype = 'OBIT_double';}
  elsif ($larg eq 'A') {$ltype = 'OBIT_string';}
  elsif ($larg eq 'L') {$ltype = 'OBIT_bool';}
  elsif ($larg eq 'C') {$ltype = 'OBIT_complex';}
  elsif ($larg eq 'M') {$ltype = 'OBIT_dcomplex';}
  elsif ($larg eq 'X') {$ltype = 'OBIT_bits';}
  elsif ($larg eq 'B') {$ltype = 'OBIT_byte';}
  else {$ltype = 'Unknown type';} # probably should be fixed in documentation
  return $ltype;
} # end GetTypeCode

# Get declaration type code, 
# first argument is document type code
sub DclTypeCode {
  local $ltype;
  local $larg = $_[0];
  if ($larg eq 'I') {$ltype = 'gshort';}
  elsif ($larg eq 'J') {$ltype = 'oint';}
  elsif ($larg eq 'K') {$ltype = 'oint';}
  elsif ($larg eq 'E') {$ltype = 'ofloat';}
  elsif ($larg eq 'D') {$ltype = 'odouble';}
  elsif ($larg eq 'A') {$ltype = 'gchar';}
  elsif ($larg eq 'L') {$ltype = 'gboolean';}
  elsif ($larg eq 'C') {$ltype = 'OBIT_complex';}
  elsif ($larg eq 'M') {$ltype = 'OBIT_dcomplex';}
  elsif ($larg eq 'X') {$ltype = 'oint';}
  elsif ($larg eq 'B') {$ltype = 'guint8';}
  else {$ltype = 'Unknown type';} # probably should be fixed in documentation
  return $ltype;
} # end DclTypeCode

#Get OBIT column element count from dimension string, -1=indexed
sub GetColumnCount {
  local $lcount, $f1, $f2, $f3, $f4, $f5;
  local $larg = $_[0];

  # parse dimension string
  ($f1,$f2,$f3,$f4,$f5) = split(',',  substr($larg, 1, length($larg)-2));
 
 # return -1 if any dimension is not a constant digit
  $lcount = -1;
  #printf (" dim %s %s %s %s %s %s \n",$larg, $f1, $f2, $f3, $f4, $f5);#debug
  if ((substr ($f1,0,1) =~ /\D/)) {return $lcount;}
  if ((substr ($f2,0,1) =~ /\D/)) {return $lcount;}
  if ((substr ($f3,0,1) =~ /\D/)) {return $lcount;}
  if ((substr ($f4,0,1) =~ /\D/)) {return $lcount;}
  if ((substr ($f5,0,1) =~ /\D/)) {return $lcount;}
#        printf (" dim %s %d %d %d %d %d \n",$larg, $f1, $f2, $f3, $f4, $f5);
  if ($f1 <= 1) {$f1 = 1};
  if ($f2 <= 1) {$f2 = 1};
  if ($f3 <= 1) {$f3 = 1};
  if ($f4 <= 1) {$f4 = 1};
  if ($f5 <= 1) {$f5 = 1};
#        printf (" dim %d %d %d %d %d \n",$f1, $f2, $f3, $f4, $f5);
  $lcount = $f1 * $f2 * $f3 * $f4 * $f5;
  return $lcount;
} # end GetColumnCount

# Any adjustment of keywords for indexed symbols
sub AdjustKeywords {
  local $i, $j, $icol;
  #print ("AdjustKeywords NumIndex =",$NumIndex,".\n");
  if ($NumIndex<=0) {return;}  # only if something to do
  $icol = -1;
  for ($i=0; $i<$NumIndex; $i++) {
    for ($j=0; $j<$NumCol; $j++) {
      if ($ColName[$j]=~$IndexCol[$i]) {
	$icol = $j+1;  # 1-rel column number for FITS
	break;
      }
    } # end loop over columns
    #print ("AdjustKeywords Column ",$IndexCol[$i]," is ",$icol,".\n");
    # adjust any keywords
    if ($icol>0) {
      for ($j=0; $j<$NumKey; $j++) {
	#print ($KeyWord[$j], " ",$IndexSymb[$i], " ",$icol," to ");
	$KeyWord[$j] =~ s/$IndexSymb[$i]/$icol/;
	#print ($KeyWord[$j]," \n");
      } # end loop adjusting keywords
    } # End make adjustment
  } # end loop over indexed symbols
} # end adjust keywords

# routine to write class files
sub WriteClass {
    local $k;
  if ($Extname ne "Uninitialized") {

    # Any adjustment of keywords for index symbols
    AdjustKeywords;

# DEBUG
#    for ($k=0; $k<$NumKey; $k++) {
#	print ("Keyword $k  $KeyWord[$k] $KeyType[$k] $KeyName[$k] def: $KeyDefV[$k] range: $KeyRang[$k]\n    $KeyComm[$k]\n");
#    }
#    for ($k=0; $k<$NumCol; $k++) {
#	print ("Column $k  $ColName[$k] $ColSuff[$k] $ColUnit[$k] $ColType[$k] $ColDim[$k] $ColVName[$k] \n    $ColComm[$k]\n");
#    }

    # Ignore the XX table (template in LaTeX document)
    if ($Extname eq "XX") {return;}

    # Say what's happening
    print ("Processing class ObitTable$ExtnameClean \n");

    # ClassInfo definition file
    $filename = "include/ObitTable".$ExtnameClean."ClassDef.h";
    $backupname = "$filename.backup";
    if (-f $filename) {rename ($filename, $backupname);}
    open (OUT, ">$filename") || 
      die "Cannot open $filename\n";
    WriteClassDef;
    close OUT;

    # Class Structure definition file
    $filename = "include/ObitTable".$ExtnameClean."Def.h";
    $backupname = "$filename.backup";
    if (-f $filename) {rename ($filename, $backupname);}
    open (OUT, ">$filename") || 
      die "Cannot open $filename\n";
    WriteDef;
    close OUT;
 
    # Class Row Structure definition file
    $filename = "include/ObitTable".$ExtnameClean."RowDef.h";
    $backupname = "$filename.backup";
    if (-f $filename) {rename ($filename, $backupname);}
    open (OUT, ">$filename") || 
	die "Cannot open $filename\n";
    WriteRowDef;
    close OUT;

    # Class Row Class Structure definition file
    $filename = "include/ObitTable".$ExtnameClean."RowClassDef.h";
    $backupname = "$filename.backup";
    if (-f $filename) {rename ($filename, $backupname);}
    open (OUT, ">$filename") || 
	die "Cannot open $filename\n";
    WriteRowClassDef;
    close OUT;

   # Class definition file
    $filename = "include/ObitTable".$ExtnameClean.".h";
    $backupname = "$filename.backup";
    if (-f $filename) {rename ($filename, $backupname);}
    open (OUT, ">$filename") || 
      die "Cannot open $filename\n";
    WriteHeader;
    close OUT;
 
  # Class function definition file
    $filename = "src/ObitTable".$ExtnameClean.".c";
    $backupname = "$filename.backup";
    if (-f $filename) {rename ($filename, $backupname);}
    open (OUT, ">$filename") || 
      die "Cannot open $filename\n";
    WriteSource;
    close OUT;

  # Write Python interface
    $filename = "python/Table".$ExtnameClean.".inc";
    $backupname = "$filename.backup";
    if (-f $filename) {rename ($filename, $backupname);}
    open (OUT, ">$filename") || 
      die "Cannot open $filename\n";
    WritePy;
    close OUT;

  } # end of write file section

  # initialize descriptor information
  $Extname = "Uninitialized"; # No longer any information
  $NumKey = 0;     # Number of keywords
  $NumCol = 0;     # Number of columns
  $NumIndex= 0;    # Number of keyword column indices.
  $NumIntro= 0;    # Number lines of introduction.
  @KeyWord = 0;    # Keyword names
  @KeyType = 0;    # Keyword types
  @KeyDefV = 0;    # Keyword default value
  @KeyRang = 0;    # Keyword range
  @KeyName = 0;    # Keyword variable names
  @KeyNumCol = 0;  # Number of columns indexed by each keyword.
  @KeyComm = 0;    # Keyword descriptions
  @ColName = 0;    # Column titles
  @ColUnit = 0;    # Column physical units
  @ColType = 0;    # Column data type code
  @ColDim  = 0;    # Column dimensionality
  @ColVName= 0;    # Column variable name roots
  @ColComm = 0;    # Column description
  @ColSuff = 0;    # Column name suffix (keyword variable name )
  @IndexSymb=0;    # Column index symbol array
  @IndexCol =0;    # Column label value array
  @TableIntro=0;   # Table introduction section
} # end WriteClass

# Write ClassInfo definition file
sub WriteClassDef {
    Copyright;
    printf (OUT "/* Define the basic components of the ObitTable$ExtnameClean ClassInfo structure */\n");
    printf (OUT "/* This is intended to be included in a classInfo structure definition*/\n");
    printf (OUT "\#include \"ObitTableClassDef.h\"  /* Parent class ClassInfo definition file */\n");
    printf (OUT "/** Function pointer to convert. */\n");
    printf (OUT "ObitTableConvertFP ObitTable$ExtnameClean\Convert;\n");
    printf (OUT "/** Function pointer to read a row. */\n");
    printf (OUT "ObitTableReadRowFP ObitTable$ExtnameClean\ReadRow;\n");
    printf (OUT "/** Function pointer to write a row. */\n");
    printf (OUT "ObitTableWriteRowFP ObitTable$ExtnameClean\WriteRow;\n");
} # end WriteClassDef

# Write ClassInfo definition file
sub WriteRowClassDef {
    Copyright;
    printf (OUT "/* Define the basic components of the ObitTable$ExtnameClean\Row ClassInfo structure */\n");
    printf (OUT "/* This is intended to be included in a classInfo structure definition*/\n");
    printf (OUT "\#include \"ObitTableRowClassDef.h\"  /* Parent class ClassInfo definition file */\n");
} # end WriteRowClassDef

# Write class Row Definition file , Count total number of Columns
sub WriteRowDef {
    local $i, $j, $k, $lo, $hi, $TC, $count, $want, $colname, $comm;

    Copyright;
   printf (OUT "/*  Define the basic components of the ObitTableRow structure       */\n");
   printf (OUT "/*  This is intended to be included in a class structure definition   */\n");
   printf (OUT "/**\n");
   printf (OUT " * \x5cfile ObitTable$ExtnameClean\RowDef.h\n");
   printf (OUT " * ObitTable$ExtnameClean\Row structure members for derived classes.\n");
   printf (OUT " */\n");
   printf (OUT "\#include \"ObitTableRowDef.h\"  /* Parent class definitions */\n");

    $TotalColumns = 0;
    # Unindexed, scalar columns
    for ($k=0; $k<$NumCol; $k++) {
	$count = GetColumnCount($ColDim[$k]);
	#print ("column $k count $count \n");
	if (($count==1) && ($ColSuff[$k] eq "")) {
	    printf (OUT "/** $ColComm[$k] */\n");
	    $TC = DclTypeCode($ColType[$k]);
	    printf (OUT "$TC  $ColVName[$k];\n");
	    $TotalColumns++;
	}
    }

    # Unindexed, array columns - done as simple pointers
    for ($k=0; $k<$NumCol; $k++) {
	$count = GetColumnCount($ColDim[$k]);
	#print ("column $k count $count \n");
	if (($count!=1) && ($ColSuff[$k] eq "")) {
	    printf (OUT "/** $ColComm[$k] */\n");
	    $TC = DclTypeCode($ColType[$k])."*";
	    printf (OUT "$TC  $ColVName[$k];\n");
	    $TotalColumns++;
	}
    }

    # indexed - loop over keywords
    for ($k=0; $k<$NumKey; $k++) {
	# it must have a range to be useful
	($lo, $hi) = split(',', substr($KeyRang[$k], 1));
	if (($lo>=1) && ($hi>=$lo)) {
	    # Loop over index values
	    for ($j=$lo; $j<=$hi; $j++) {

		# Loop over columns
		for ($i=0; $i<$NumCol; $i++) {
		    $colname =  $ColVName[$i]."$j";
		    #print (" col $i $colname $ColSuff[$i] $KeyWord[$k] \n");
		    $count = GetColumnCount($ColDim[$i]);

		    # do we want this one? indexed by keyword $k
		    $want = $ColSuff[$i] eq NoQuote($KeyWord[$k]);
		    #$comm =  NoQuote($KeyWord[$k]);
		    #print (" :$ColSuff[$i]: and :$comm:\n");

		    # prepare Comment - replace index name with value
		    $comm =  $ColComm[$i];
		    $comm =~ s/$ColSuff[$i]/$j/g;

		    if (($count==1) && ($want)) {
			# Scalar, indexed columns
			printf (OUT "/** $comm */\n");
			$TC = DclTypeCode($ColType[$i]);
			printf (OUT "$TC  $colname;\n");
			$TotalColumns++;

		    } elsif (($count!=1) && ($want)){
			# array indexed columns - done as simple pointers
			printf (OUT "/** $comm */\n");
			$TC = DclTypeCode($ColType[$i])."*";
			printf (OUT "$TC  $colname;\n");
			$TotalColumns++;
		    }
		} # end loop over columns
	
	    } # end loop over indices
	}
    } # end loop over keywords

    # Add _status column
    printf (OUT "/** status 0=normal, 1=modified, -1=flagged */\n");
    printf (OUT "olong  status;\n");

   #print ("Total number of colums for $Extname table is $TotalColumns\n");

} # end WriteRowDef

# Write Class Structure definition file, count types of columns
sub WriteDef {
    local $i, $j, $k, $l, $TC, $TX, $comm, $want;

    Copyright;

    printf (OUT "/*  Define the basic components of the ObitTable$ExtnameClean  structure          */\n");
    printf (OUT "/*  This is intended to be included in a class structure definition   */\n");
    printf (OUT "/**\n");
    printf (OUT " * \x5cfile ObitTable$ExtnameClean\Def.h\n");
    printf (OUT " * ObitTable$ExtnameClean structure members for derived classes.\n");
    printf (OUT " */\n");
    printf (OUT "\#include \"ObitTableDef.h\"  /* Parent class definitions */\n");

    # Write keywords
    for ($k=0; $k<$NumKey; $k++) {
	printf (OUT "/** $KeyComm[$k] */\n");
	$TC = DclTypeCode($KeyType[$k]);
	$TX = ""; # special case for strings
	if ($KeyType[$k] eq "A") {$TX = "[MAXKEYCHARTABLE".$ExtnameClean."]";}
	printf (OUT "$TC  $KeyName[$k]$TX;\n");
    }


    # Number of non indexed columns
    $NumColNoindex = 0;

    # Unindexed, scalar columns
    for ($k=0; $k<$NumCol; $k++) {
	$count = GetColumnCount($ColDim[$k]);
	#print ("column $k count $count \n");
	if (($count==1) && ($ColSuff[$k] eq "")) {
	    printf (OUT "/** Column offset for $ColComm[$k] in table record */\n");
	    printf (OUT "olong  $ColVName[$k]Off;\n");
	    printf (OUT "/** Physical column number for $ColComm[$k] in table record */\n");
	    printf (OUT "olong  $ColVName[$k]Col;\n");
	    $NumColNoindex++;
	}
    }

    # Unindexed, array columns
    for ($k=0; $k<$NumCol; $k++) {
	$count = GetColumnCount($ColDim[$k]);
	#print ("column $k $ColVName[$k] count $count $ColDim[$k] $ColSuff[$k] \n");
	if (($count!=1) && ($ColSuff[$k] eq "")) {
	    printf (OUT "/** Column offset for $ColComm[$k] in table record */\n");
	    printf (OUT "olong  $ColVName[$k]Off;\n");
	    printf (OUT "/** Physical column number for $ColComm[$k] in table record */\n");
	    printf (OUT "olong  $ColVName[$k]Col;\n");
	    $NumColNoindex++;
	}
    }

    # indexed - loop over keywords
    for ($k=0; $k<$NumKey; $k++) {
	# it must have a range to be useful
	($lo, $hi) = split(',', substr($KeyRang[$k], 1));
	if (($lo>=1) && ($hi>=$lo)) {
	    # Loop over index values
	    for ($j=$lo; $j<=$hi; $j++) {

		# Loop over columns
		for ($i=0; $i<$NumCol; $i++) {
		    $colname =  $ColVName[$i]."$j";
		    #print (" col $i $colname $ColSuff[$i] $KeyWord[$k] \n");

		    # do we want this one? indexed by keyword $k
		    $want = $ColSuff[$i] eq NoQuote($KeyWord[$k]);

		    # prepare Comment - replace index name with value
		    $comm =  $ColComm[$i];
		    $comm =~ s/$ColSuff[$i]/$j/g;

		    if ($want) {
			printf (OUT "/** Column offset for $comm in table record */\n");
			printf (OUT "olong  $colname\Off;\n");
			printf (OUT "/** Physical column number for $comm in table record */\n");
			printf (OUT "olong  $colname\Col;\n");
			# How many columns indexed by this keyword
			if ($j==$k) {$KeyNumCol[$k]++;}
		    }
		} # end loop over columns
		
	    } # end loop over indices
	}
    } # end loop over keywords
} # end WriteDef

# Write Class definition file
sub WriteHeader {
    local $k, $temp, TC;

    Copyright;

    printf (OUT "\#ifndef OBITTABLE".$ExtnameClean."_H \n");
    printf (OUT "\#define OBITTABLE".$ExtnameClean."_H \n");
    printf (OUT "\n");
    printf (OUT "\#include \"Obit.h\"\n");
    printf (OUT "\#include \"ObitErr.h\"\n");
    printf (OUT "\#include \"ObitTable.h\"\n");
    printf (OUT "\#include \"ObitData.h\"\n");
    printf (OUT "\n");
    printf (OUT "/*-------- Obit: Merx mollis mortibus nuper ------------------*/\n");
    printf (OUT "/**\n");
    printf (OUT " * \x5cfile ObitTable".$ExtnameClean.".h\n");
    printf (OUT " * ObitTable".$ExtnameClean." class definition.\n");
    printf (OUT " *\n");
    printf (OUT " * This class is derived from the \#ObitTable class.\n");
    printf (OUT " *\n");
    # Write introduction section from document
    for ($k = 0; $k<$NumIntro; $k++) {
	printf (OUT " * $TableIntro[$k]");
   }
    printf (OUT " *\n");
    printf (OUT " * This class contains tabular data and allows access.\n");
    printf (OUT " * \"".$ExtType.$ExtnameClean."\" table\n");
    printf (OUT " * An ObitTable".$ExtnameClean." is the front end to a persistent disk resident structure.\n");
    if ($ExtType eq "AIPS ") {
      printf (OUT " * Both FITS (as Tables) and AIPS cataloged data are supported.\n");
    } else {
      printf (OUT " * Only FITS (as Tables) are supported.\n");
    }
    printf (OUT " *\n");
    printf (OUT " * \x5csection TableDataStorage Table data storage\n");
    printf (OUT " * In memory tables are stored in a fashion similar to how they are \n");
    printf (OUT " * stored on disk - in large blocks in memory rather than structures.\n");
    printf (OUT " * Due to the word alignment requirements of some machines, they are \n");
    printf (OUT " * stored by order of the decreasing element size: \n");
    printf (OUT " * double, float long, int, short, char rather than the logical order.\n");
    printf (OUT " * The details of the storage in the buffer are kept in the \n");
    printf (OUT " * \#ObitTable".$ExtnameClean."Desc.\n");
	printf (OUT " *\n");
    printf (OUT " * In addition to the normal tabular data, a table will have a \"\_status\"\n");
    printf (OUT " * column to indicate the status of each row.\n");
    if ($ExtType eq "AIPS ") {
      printf (OUT " * The status value is read from and written to (some modification) AIPS \n");
      printf (OUT " * tables but are not written to externally generated FITS tables which\n");
      printf (OUT " * don't have these colummns.  It will be written to Obit generated tables\n");
      printf (OUT " * which will have these columns.\n");
      printf (OUT " * Status values:\n");
      printf (OUT " * \x5cli status = 0 => normal\n");
      printf (OUT " * \x5cli status = 1 => row has been modified (or created) and needs to be\n");
      printf (OUT " *                   written.\n");
      printf (OUT " * \x5cli status = -1 => row has been marked invalid.\n");
    }
    printf (OUT " *\n");
    printf (OUT " * \x5csection ObitTable".$ExtnameClean."Specification Specifying desired data transfer parameters\n");
    printf (OUT " * The desired data transfers are specified in the member ObitInfoList.\n");
    if ($ExtType eq "AIPS ") {
      printf (OUT " * There are separate sets of parameters used to specify the FITS or AIPS \n");
      printf (OUT " * data files.\n");
    }
    printf (OUT " * In the following an ObitInfoList entry is defined by \n");
    printf (OUT " * the name in double quotes, the data type code as an \#ObitInfoType enum \n");
    printf (OUT " * and the dimensions of the array (? => depends on application).\n");
    if ($ExtType eq "AIPS ") {
      printf (OUT " * To specify whether the underlying data files are FITS or AIPS\n");
      printf (OUT " * \x5cli \"FileType\" OBIT_int (1,1,1) OBIT_IO_FITS or OBIT_IO_AIPS \n");
      printf (OUT " * which are values of an \#ObitIOType enum defined in ObitIO.h.\n");
    }
    printf (OUT " *\n");
    printf (OUT " * The following apply to both types of files:\n");
    printf (OUT " * \x5cli \"nRowPIO\", OBIT_int, Max. Number of visibilities per \n");
    printf (OUT " *     \"Read\" or \"Write\" operation.  Default = 1.\n");
    printf (OUT " *\n");
    printf (OUT " * \x5csubsection TableFITS FITS files\n");
    printf (OUT " * This implementation uses cfitsio which allows using, in addition to \n");
    printf (OUT " * regular FITS images, gzip compressed files, pipes, shared memory \n");
    printf (OUT " * and a number of other input forms.\n");
    printf (OUT " * The convenience Macro \#ObitTable".$ExtnameClean."SetFITS simplifies specifying the \n");
	printf (OUT " * desired data.\n");
    printf (OUT " * Binary tables are used for storing visibility data in FITS.\n");
    printf (OUT " * For accessing FITS files the following entries in the ObitInfoList \n");
    printf (OUT " * are used:\n");
    printf (OUT " * \x5cli \"FileName\" OBIT_string (?,1,1) FITS file name.\n");
    printf (OUT " * \x5cli \"TabName\"  OBIT_string (?,1,1) Table name (e.g. \"AIPS CC\").\n");
    printf (OUT " * \x5cli \"Ver\"      OBIT_int    (1,1,1) Table version number\n");
    printf (OUT " *\n");
    if ($ExtType eq "AIPS ") {
      printf (OUT " * \subsection ObitTable".$ExtnameClean."AIPS AIPS files\n");
      printf (OUT " * The ObitAIPS class must be initialized before accessing AIPS files; \n");
      printf (OUT " * this uses \#ObitAIPSClassInit.\n");
      printf (OUT " * For accessing AIPS files the following entries in the ObitInfoList \n");
      printf (OUT " * are used:\n");
      printf (OUT " * \x5cli \"Disk\" OBIT_int (1,1,1) AIPS \"disk\" number.\n");
      printf (OUT " * \x5cli \"User\" OBIT_int (1,1,1) user number.\n");
      printf (OUT " * \x5cli \"CNO\"  OBIT_int (1,1,1) AIPS catalog slot number.\n");
      printf (OUT " * \x5cli \"TableType\" OBIT_string (2,1,1) AIPS Table type\n");
      printf (OUT " * \x5cli \"Ver\"  OBIT_int    (1,1,1) AIPS table version number.\n");
    }
    printf (OUT " *\n");
    printf (OUT " * \x5csection ObitTable".$ExtnameClean."access Creators and Destructors\n");
    printf (OUT " * An ObitTable".$ExtnameClean." can be created using newObitTable".$ExtnameClean."Value which attaches the \n");
    printf (OUT " * table to an ObitData for the object.  \n");
    printf (OUT " * If the output ObitTable".$ExtnameClean." has previously been specified, including file information,\n");
    printf (OUT " * then ObitTable".$ExtnameClean."Copy will copy the disk resident as well as the memory \n");
    printf (OUT " * resident information.\n");
    printf (OUT " *\n");
    printf (OUT " * A copy of a pointer to an ObitTable".$ExtnameClean." should always be made using the\n");
    printf (OUT " * ObitTable".$ExtnameClean."Ref function which updates the reference count in the object.\n");
    printf (OUT " * Then whenever freeing an ObitTable".$ExtnameClean." or changing a pointer, the function\n");
    printf (OUT " * ObitTable".$ExtnameClean."Unref will decrement the reference count and destroy the object\n");
    printf (OUT " * when the reference count hits 0.\n");
    printf (OUT " *\n");
    printf (OUT " * \x5csection ObitTable".$ExtnameClean."Usage I/O\n");
    printf (OUT " * Visibility data is available after an input object is \"Opened\"\n");
    printf (OUT " * and \"Read\".\n");
    printf (OUT " * I/O optionally uses a buffer attached to the ObitTable".$ExtnameClean." or some external\n");
    printf (OUT " * location.\n");
    printf (OUT " * To Write an ObitTable".$ExtnameClean.", create it, open it, and write.\n");
    printf (OUT " * The object should be closed to ensure all data is flushed to disk.\n");
    printf (OUT " * Deletion of an ObitTable".$ExtnameClean." after its final unreferencing will automatically\n");
    printf (OUT " * close it.\n");
    printf (OUT " */\n");
    printf (OUT "\n");
    printf (OUT "/*--------------Class definitions-------------------------------------*/\n");
    printf (OUT "\n");
    printf (OUT "/** Number of characters for Table keyword */\n");
    printf (OUT " #define MAXKEYCHARTABLE".$ExtnameClean." 24\n");
    printf (OUT "\n");
    printf (OUT "/** ObitTable".$ExtnameClean." Class structure. */\n");
    printf (OUT "typedef struct {\n");
    printf (OUT "\#include \"ObitTable".$ExtnameClean."Def.h\"   /* this class definition */\n");
    printf (OUT "} ObitTable".$ExtnameClean.";\n");
    printf (OUT "\n");
    printf (OUT "/** ObitTable".$ExtnameClean."Row Class structure. */\n");
    printf (OUT "typedef struct {\n");
    printf (OUT "\#include \"ObitTable".$ExtnameClean."RowDef.h\"   /* this class row definition */\n");
    printf (OUT "} ObitTable".$ExtnameClean."Row;\n");
    printf (OUT "\n");
    printf (OUT "/*----------------- Macroes ---------------------------*/\n");
    printf (OUT "/** \n");
    printf (OUT " * Macro to unreference (and possibly destroy) an ObitTable".$ExtnameClean."\n");
    printf (OUT " * returns an ObitTable".$ExtnameClean."*.\n");
    printf (OUT " * in = object to unreference\n");
    printf (OUT " */\n");
    printf (OUT "\#define ObitTable".$ExtnameClean."Unref(in) ObitUnref (in)\n");
    printf (OUT "\n");
    printf (OUT "/** \n");
    printf (OUT " * Macro to reference (update reference count) an ObitTable".$ExtnameClean.".\n");
    printf (OUT " * returns an ObitTable".$ExtnameClean."*.\n");
    printf (OUT " * in = object to reference\n");
    printf (OUT " */\n");
    printf (OUT "\#define ObitTable".$ExtnameClean."Ref(in) ObitRef (in)\n");
    printf (OUT "\n");
    printf (OUT "/** \n");
    printf (OUT " * Macro to determine if an object is the member of this or a \n");
    printf (OUT " * derived class.\n");
    printf (OUT " * Returns TRUE if a member, else FALSE\n");
    printf (OUT " * in = object to reference\n");
    printf (OUT " */\n");
    printf (OUT "\#define ObitTable".$ExtnameClean."IsA(in) ObitIsA (in, ObitTable".$ExtnameClean."GetClass())\n");
    printf (OUT "\n");
    printf (OUT "/** \n");
    printf (OUT " * Macro to unreference (and possibly destroy) an ObitTable".$ExtnameClean."Row\n");
    printf (OUT " * returns an ObitTable".$ExtnameClean."Row*.\n");
    printf (OUT " * in = object to unreference\n");
    printf (OUT " */\n");
    printf (OUT "\#define ObitTable".$ExtnameClean."RowUnref(in) ObitUnref (in)\n");
    printf (OUT "\n");
    printf (OUT "/** \n");
    printf (OUT " * Macro to reference (update reference count) an ObitTable".$ExtnameClean."Row.\n");
    printf (OUT " * returns an ObitTable".$ExtnameClean."Row*.\n");
    printf (OUT " * in = object to reference\n");
    printf (OUT " */\n");
    printf (OUT "\#define ObitTable".$ExtnameClean."RowRef(in) ObitRef (in)\n");
    printf (OUT "\n");
    printf (OUT "/** \n");
    printf (OUT " * Macro to determine if an object is the member of this or a \n");
    printf (OUT " * derived class.\n");
    printf (OUT " * Returns TRUE if a member, else FALSE\n");
    printf (OUT " * in = object to reference\n");
    printf (OUT " */\n");
    printf (OUT "\#define ObitTable".$ExtnameClean."RowIsA(in) ObitIsA (in, ObitTable".$ExtnameClean."RowGetClass())\n");
    printf (OUT "\n");
    printf (OUT "/*---------------Public functions---------------------------*/\n");
    printf (OUT "/*----------------Table Row Functions ----------------------*/\n");
    printf (OUT "/** Public: Row Class initializer. */\n");
    printf (OUT "void ObitTable".$ExtnameClean."RowClassInit (void);\n");
    printf (OUT "\n");
    printf (OUT "/** Public: Constructor. */\n");
    printf (OUT "ObitTable".$ExtnameClean."Row* newObitTable".$ExtnameClean."Row (ObitTable".$ExtnameClean." *table);\n");
    printf (OUT "\n");
    printf (OUT "/** Public: ClassInfo pointer */\n");
    printf (OUT "gconstpointer ObitTable".$ExtnameClean."RowGetClass (void);\n");
    printf (OUT "\n");
    printf (OUT "/*------------------Table Functions ------------------------*/\n");
    printf (OUT "/** Public: Class initializer. */\n");
    printf (OUT "void ObitTable".$ExtnameClean."ClassInit (void);\n");
    printf (OUT "\n");
    printf (OUT "/** Public: Constructor. */\n");
    printf (OUT "ObitTable".$ExtnameClean."* newObitTable".$ExtnameClean." (gchar* name);\n");
    printf (OUT "\n");
    printf (OUT "/** Public: Constructor from values. */\n");
    printf (OUT "ObitTable".$ExtnameClean."* \n");
    printf (OUT "newObitTable".$ExtnameClean."Value (gchar* name, ObitData *file, olong *ver,\n");
    printf (OUT "  		     ObitIOAccess access,\n");
    # Add structural keywords to call sequence
    $temp =     "                  ";
    for ($k=0; $k<$NumKey; $k++) {
	if ($KeyRang[$k] ne "") {
	    $TC = DclTypeCode($KeyType[$k]);
	    $temp = $temp." $TC $KeyName[$k],";
	}
    }
    printf (OUT "  $temp\n");
    printf (OUT "		     ObitErr *err);\n");
    printf (OUT "\n");
    printf (OUT "/** Public: Class initializer. */\n");
    printf (OUT "void ObitTable".$ExtnameClean."ClassInit (void);\n");
    printf (OUT "\n");
    printf (OUT "/** Public: ClassInfo pointer */\n");
    printf (OUT "gconstpointer ObitTable".$ExtnameClean."GetClass (void);\n");
    printf (OUT "\n");
    printf (OUT "/** Public: Copy (deep) constructor. */\n");
    printf (OUT "ObitTable".$ExtnameClean."* ObitTable".$ExtnameClean."Copy  (ObitTable".$ExtnameClean." *in, ObitTable".$ExtnameClean." *out, \n");
    printf (OUT "			   ObitErr *err);\n");
    printf (OUT "\n");
    printf (OUT "/** Public: Copy (shallow) constructor. */\n");
    printf (OUT "ObitTable".$ExtnameClean."* ObitTable".$ExtnameClean."Clone (ObitTable".$ExtnameClean." *in, ObitTable".$ExtnameClean." *out);\n");
    printf (OUT "\n");
    printf (OUT "/** Public: Convert an ObitTable to an ObitTable".$ExtnameClean." */\n");
    printf (OUT "ObitTable".$ExtnameClean."* ObitTable".$ExtnameClean."Convert  (ObitTable *in);\n");
    printf (OUT "\n");
    printf (OUT "/** Public: Create ObitIO structures and open file */\n");
    printf (OUT "ObitIOCode ObitTable".$ExtnameClean."Open (ObitTable".$ExtnameClean." *in, ObitIOAccess access, \n");
    printf (OUT "			  ObitErr *err);\n");
    printf (OUT "\n");
    printf (OUT "/** Public: Read a table row */\n");
    printf (OUT "ObitIOCode \n");
    printf (OUT "ObitTable".$ExtnameClean."ReadRow  (ObitTable".$ExtnameClean." *in, olong i".$ExtnameClean."Row, ObitTable".$ExtnameClean."Row *row,\n");
    printf (OUT "		     ObitErr *err);\n");
    printf (OUT "\n");
    printf (OUT "/** Public: Init a table row for write */\n");
    printf (OUT "void \n");
    printf (OUT "ObitTable".$ExtnameClean."SetRow  (ObitTable".$ExtnameClean." *in, ObitTable".$ExtnameClean."Row *row,\n");
    printf (OUT "		     ObitErr *err);\n");
    printf (OUT "\n");
    printf (OUT "/** Public: Write a table row */\n");
    printf (OUT "ObitIOCode \n");
    printf (OUT "ObitTable".$ExtnameClean."WriteRow  (ObitTable".$ExtnameClean." *in, olong i".$ExtnameClean."Row, ObitTable".$ExtnameClean."Row *row,\n");
    printf (OUT "		     ObitErr *err);\n");
    printf (OUT "\n");
    printf (OUT "/** Public: Close file and become inactive */\n");
    printf (OUT "ObitIOCode ObitTable".$ExtnameClean."Close (ObitTable".$ExtnameClean." *in, ObitErr *err);\n");
    printf (OUT "\n");
    printf (OUT "/*----------- ClassInfo Structure -----------------------------------*/\n");
    printf (OUT "/**\n");
    printf (OUT " * ClassInfo Structure.\n");
    printf (OUT " * Contains class name, a pointer to any parent class\n");
    printf (OUT " * (NULL if none) and function pointers.\n");
    printf (OUT " */\n");
    printf (OUT "typedef struct  {\n");
    printf (OUT "\#include \"ObitTable".$ExtnameClean."ClassDef.h\"\n");
    printf (OUT "} ObitTable".$ExtnameClean."ClassInfo; \n");
    printf (OUT "\n");
    printf (OUT "/**\n");
    printf (OUT " * ClassInfo Structure For Table".$ExtnameClean."Row.\n");
    printf (OUT " * Contains class name, a pointer to any parent class\n");
    printf (OUT " * (NULL if none) and function pointers.\n");
    printf (OUT " */\n");
    printf (OUT "typedef struct  {\n");
    printf (OUT "\#include \"ObitTable".$ExtnameClean."RowClassDef.h\"\n");
    printf (OUT "} ObitTable".$ExtnameClean."RowClassInfo; \n");
    printf (OUT "\#endif /* OBITTABLE".$ExtnameClean."_H */ \n");
} # end WriteHeader

# Write Copyright
sub Copyright {
    printf (OUT "/* \x24Id:  \x24   */\n");
    printf (OUT "/* DO NOT EDIT - file generated by ObitTables.pl                      */\n");
    printf (OUT "/*--------------------------------------------------------------------*/\n");
    printf (OUT "/*;  Copyright (C)  2010                                              */\n");
    printf (OUT "/*;  Associated Universities, Inc. Washington DC, USA.                */\n");
    printf (OUT "/*;                                                                   */\n");
    printf (OUT "/*;  This program is free software; you can redistribute it and/or    */\n");
    printf (OUT "/*;  modify it under the terms of the GNU General Public License as   */\n");
    printf (OUT "/*;  published by the Free Software Foundation; either version 2 of   */\n");
    printf (OUT "/*;  the License, or (at your option) any later version.              */\n");
    printf (OUT "/*;                                                                   */\n");
    printf (OUT "/*;  This program is distributed in the hope that it will be useful,  */\n");
    printf (OUT "/*;  but WITHOUT ANY WARRANTY; without even the implied warranty of   */\n");
    printf (OUT "/*;  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    */\n");
    printf (OUT "/*;  GNU General Public License for more details.                     */\n");
    printf (OUT "/*;                                                                   */\n");
    printf (OUT "/*;  You should have received a copy of the GNU General Public        */\n");
    printf (OUT "/*;  License along with this program; if not, write to the Free       */\n");
    printf (OUT "/*;  Software Foundation, Inc., 675 Massachusetts Ave, Cambridge,     */\n");
    printf (OUT "/*;  MA 02139, USA.                                                   */\n");
    printf (OUT "/*;                                                                   */\n");
    printf (OUT "/*;Correspondence about this software should be addressed as follows: */\n");
    printf (OUT "/*;         Internet email: bcotton\@nrao.edu.                         */\n");
    printf (OUT "/*;         Postal address: William Cotton                            */\n");
    printf (OUT "/*;                         National Radio Astronomy Observatory      */\n");
    printf (OUT "/*;                         520 Edgemont Road                         */\n");
    printf (OUT "/*;                         Charlottesville, VA 22903-2475 USA        */\n");
    printf (OUT "/*--------------------------------------------------------------------*/\n");
} # end Copyright 

# Write preamble to source file
sub WritePreamble {
    printf (OUT "\#include \"ObitTable".$ExtnameClean.".h\"\n");
    printf (OUT "\#include \"ObitTableList.h\"\n");
    printf (OUT "\#include \"ObitData.h\"\n");
    printf (OUT "\n");
    printf (OUT "/*----------------Obit:  Merx mollis mortibus nuper ------------------*/\n");
    printf (OUT "/**\n");
    printf (OUT " * \x5cfile ObitTable".$ExtnameClean.".c\n");
    printf (OUT " * ObitTable".$ExtnameClean." class function definitions.\n");
    printf (OUT " *\n");
    printf (OUT " * This class is derived from the \#ObitTable class.\n");
    printf (OUT " */\n");
    printf (OUT "\n");
    printf (OUT "/** name of the class defined in this file */\n");
    printf (OUT "static gchar *myClassName = \"ObitTable".$ExtnameClean."\";\n");
    printf (OUT "\n");
    printf (OUT "/**  Function to obtain parent Table ClassInfo - ObitTable */\n");
    printf (OUT "static ObitGetClassFP ObitParentGetClass = ObitTableGetClass;\n");
    printf (OUT "\n");
    printf (OUT "/** name of the Row class defined in this file */\n");
    printf (OUT "static gchar *myRowClassName = \"ObitTable".$ExtnameClean."Row\";\n");
    printf (OUT "\n");
    printf (OUT "/**  Function to obtain parent TableRow ClassInfo */\n");
    printf (OUT "static ObitGetClassFP ObitParentGetRowClass = ObitTableRowGetClass;\n");
    printf (OUT "\n");
    printf (OUT "/*--------------- File Global Variables  ----------------*/\n");
    printf (OUT "/*----------------  Table Row  ----------------------*/\n");
    printf (OUT "/**\n");
    printf (OUT " * ClassInfo structure ObitTableClassInfo.\n");
    printf (OUT " * This structure is used by class objects to access class functions.\n");
    printf (OUT " */\n");
    printf (OUT "static ObitTable".$ExtnameClean."RowClassInfo myRowClassInfo = {FALSE};\n");
    printf (OUT "\n");
    printf (OUT "/*------------------  Table  ------------------------*/\n");
    printf (OUT "/**\n");
    printf (OUT " * ClassInfo structure ObitTable".$ExtnameClean."ClassInfo.\n");
    printf (OUT " * This structure is used by class objects to access class functions.\n");
    printf (OUT " */\n");
    printf (OUT "static ObitTable".$ExtnameClean."ClassInfo myClassInfo = {FALSE};\n");
    printf (OUT "\n");
    printf (OUT "/*---------------Private function prototypes----------------*/\n");
    printf (OUT "/** Private: Initialize newly instantiated Row object. */\n");
    printf (OUT "void  ObitTable".$ExtnameClean."RowInit  (gpointer in);\n");
    printf (OUT "\n");
    printf (OUT "/** Private: Deallocate Row members. */\n");
    printf (OUT "void  ObitTable".$ExtnameClean."RowClear (gpointer in);\n");
    printf (OUT "\n");
    printf (OUT "/** Private: Initialize newly instantiated object. */\n");
    printf (OUT "void  ObitTable".$ExtnameClean."Init  (gpointer in);\n");
    printf (OUT "\n");
    printf (OUT "/** Private: Deallocate members. */\n");
    printf (OUT "void  ObitTable".$ExtnameClean."Clear (gpointer in);\n");
    printf (OUT "\n");
    printf (OUT "/** Private: update table specific info */\n");
    printf (OUT "static void ObitTable".$ExtnameClean."Update (ObitTable".$ExtnameClean." *in, ObitErr *err);\n");
    printf (OUT "\n");
    printf (OUT "/** Private: copy table keywords to descriptor info list */\n");
    printf (OUT "static void ObitTable".$ExtnameClean."DumpKey (ObitTable".$ExtnameClean." *in, ObitErr *err);\n");
    printf (OUT "\n");
    printf (OUT "/** Private: Set Class function pointers */\n");
    printf (OUT "static void ObitTable".$ExtnameClean."ClassInfoDefFn (gpointer inClass);\n");
    printf (OUT "\n");
    printf (OUT "/** Private: Set Row Class function pointers */\n");
    printf (OUT "static void ObitTable".$ExtnameClean."RowClassInfoDefFn (gpointer inClass);\n");
    printf (OUT "/*----------------------Public functions---------------------------*/\n");
    printf (OUT "\n");
    printf (OUT "/*------------------  Table Row ------------------------*/\n");
} # end WritePreamble

sub WriteTableRow {
    printf (OUT "/**\n");
    printf (OUT " * Constructor.\n");
    printf (OUT " * If table is open and for write, the row is attached to the buffer\n");
    printf (OUT " * Initializes Row class if needed on first call.\n");
    printf (OUT " * \x5cparam name An optional name for the object.\n");
    printf (OUT " * \x5creturn the new object.\n");
    printf (OUT " */\n");
    printf (OUT "ObitTable".$ExtnameClean."Row* newObitTable".$ExtnameClean."Row (ObitTable".$ExtnameClean." *table)\n");
    printf (OUT "{\n");
    printf (OUT "  ObitTable".$ExtnameClean."Row* out;\n");
    printf (OUT "  odouble   *dRow;\n");
    printf (OUT "  oint      *iRow;\n");
    printf (OUT "  gshort    *siRow;\n");
    printf (OUT "  ofloat    *fRow;\n");
    printf (OUT "  gchar     *cRow;\n");
    printf (OUT "  gboolean  *lRow;\n");
    printf (OUT "  guint8    *bRow;\n");
    printf (OUT "\n");
    printf (OUT "  /* Class initialization if needed */\n");
    printf (OUT "  if (!myRowClassInfo.initialized) ObitTable".$ExtnameClean."RowClassInit();\n");
    printf (OUT "\n");
    printf (OUT "  /* allocate/init structure */\n");
    printf (OUT "  out = g_malloc0(sizeof(ObitTable".$ExtnameClean."Row));\n");
    printf (OUT "\n");
    printf (OUT "  /* initialize values */\n");
    printf (OUT "  out->name = g_strdup(\"Table".$ExtnameClean." Row\");\n");
    printf (OUT "\n");
    printf (OUT "  /* set ClassInfo */\n");
    printf (OUT "  out->ClassInfo = (gpointer)&myRowClassInfo;\n");
    printf (OUT "\n");
    printf (OUT "  /* initialize other stuff */\n");
    printf (OUT "  ObitTable".$ExtnameClean."RowInit((gpointer)out);\n");
    printf (OUT "  out->myTable   = (ObitTable*)ObitTableRef((ObitTable*)table);\n");
    printf (OUT "\n");
    printf (OUT "  /* If writing attach to buffer */\n");
    printf (OUT "  if ((table->buffer) && (table->myDesc->access != OBIT_IO_ReadOnly) &&\n");
    printf (OUT "      (table->myStatus != OBIT_Inactive)) {\n");
    printf (OUT "    /* Typed pointers to row of data */  \n");
    printf (OUT "    dRow  = (odouble*)table->buffer;\n");
    printf (OUT "    iRow  = (oint*)table->buffer;\n");
    printf (OUT "    siRow = (gshort*)table->buffer;\n");
    printf (OUT "    fRow  = (ofloat*)table->buffer;\n");
    printf (OUT "    cRow  = (gchar*)table->buffer;\n");
    printf (OUT "    lRow  = (gboolean*)table->buffer;\n");
    printf (OUT "    bRow  = (guint8*)table->buffer;\n");
    printf (OUT "  \n");
    printf (OUT "    /* Set row pointers to buffer */\n");
    # Unindexed, array columns
    for ($k=0; $k<$NumCol; $k++) {
	$count = GetColumnCount($ColDim[$k]);
	if (($count!=1) && ($ColSuff[$k] eq "")) {
	    $TC = RowType($ColType[$k]);
	    printf (OUT "    out->$ColVName[$k] = $TC + table->$ColVName[$k]Off;\n");
	}
    }

    # indexed - loop over keywords
    for ($k=0; $k<$NumKey; $k++) {
	# it must have a range to be useful
	($lo, $hi) = split(',', substr($KeyRang[$k], 1));
	if (($lo>=1) && ($hi>=$lo)) {
	    # Loop over index values
	    for ($j=$lo; $j<=$hi; $j++) {

		# Loop over columns
		for ($i=0; $i<$NumCol; $i++) {
		    $colname =  $ColVName[$i]."$j";
		    $count = GetColumnCount($ColDim[$i]);
		    #print (" col $i $colname $ColSuff[$i] $KeyWord[$k] \n");

		    # do we want this one? indexed by keyword $k
		    $want = $ColSuff[$i] eq NoQuote($KeyWord[$k]);
		    
		    if (($count==1) && ($want)) {
			# scalar indexed columns
			$colname = $ColVName[$i]."$j";
			$TC = RowType($ColType[$i])."[table->$colname\Off]";
			printf (OUT "    if (table->$KeyName[$k]>=$j) out->$colname = $TC;\n");
		    } elsif (($count!=1) && ($want)){
			# array indexed columns
			$colname = $ColVName[$i]."$j";
			$TC = RowType($ColType[$i]);
			printf (OUT "    if (table->$KeyName[$k]>=$j) out->$colname = $TC + table->$colname\Off;\n");
		    }
		} # end loop over columns
		
	    } # end loop over indices
	}
    } # end loop over keywords
    printf (OUT "  } /* end attaching row to table buffer */\n");
    printf (OUT "\n");
    printf (OUT " return out;\n");
    printf (OUT "} /* end newObitTable".$ExtnameClean."Row */\n");
    printf (OUT "\n");
    printf (OUT "/**\n");
    printf (OUT " * Returns ClassInfo pointer for the Row class.\n");
    printf (OUT " * \x5creturn pointer to the Row class structure.\n");
    printf (OUT " */\n");
    printf (OUT "gconstpointer ObitTable".$ExtnameClean."RowGetClass (void)\n");
    printf (OUT "{\n");
    printf (OUT "  /* Class initialization if needed */\n");
    printf (OUT "  if (!myRowClassInfo.initialized) ObitTable".$ExtnameClean."RowClassInit();\n");
    printf (OUT "  return (gconstpointer)&myRowClassInfo;\n");
    printf (OUT "} /* end ObitTable".$ExtnameClean."RowGetClass */\n");
    printf (OUT "\n");
} # End WriteTableRow

sub WriteConstructor {
    local $k, $temp, $TC, $want, $count, $f1, $f2, $f3, $f4, $f5, $colname;
    printf (OUT "/*------------------  Table  ------------------------*/\n");
    printf (OUT "/**\n");
    printf (OUT " * Constructor.\n");
    printf (OUT " * Initializes class if needed on first call.\n");
    printf (OUT " * \x5cparam name An optional name for the object.\n");
    printf (OUT " * \x5creturn the new object.\n");
    printf (OUT " */\n");
    printf (OUT "ObitTable".$ExtnameClean."* newObitTable".$ExtnameClean." (gchar* name)\n");
    printf (OUT "{\n");
    printf (OUT "  ObitTable".$ExtnameClean."* out;\n");
    printf (OUT "\n");
    printf (OUT "  /* Class initialization if needed */\n");
    printf (OUT "  if (!myClassInfo.initialized) ObitTable".$ExtnameClean."ClassInit();\n");
    printf (OUT "\n");
    printf (OUT "  /* allocate/init structure */\n");
    printf (OUT "  out = g_malloc0(sizeof(ObitTable".$ExtnameClean."));\n");
    printf (OUT "\n");
    printf (OUT "  /* initialize values */\n");
    printf (OUT "  if (name!=NULL) out->name = g_strdup(name);\n");
    printf (OUT "  else out->name = g_strdup(\"Noname\");\n");
    printf (OUT "\n");
    printf (OUT "  /* set ClassInfo */\n");
    printf (OUT "  out->ClassInfo = (gpointer)&myClassInfo;\n");
    printf (OUT "\n");
    printf (OUT "  /* initialize other stuff */\n");
    printf (OUT "  ObitTable".$ExtnameClean."Init((gpointer)out);\n");
    printf (OUT "\n");
    printf (OUT " return out;\n");
    printf (OUT "} /* end newObitTable".$ExtnameClean." */\n");
    printf (OUT "\n");
    printf (OUT "/**\n");
    printf (OUT " * Returns ClassInfo pointer for the class.\n");
    printf (OUT " * \x5creturn pointer to the class structure.\n");
    printf (OUT " */\n");
    printf (OUT "gconstpointer ObitTable".$ExtnameClean."GetClass (void)\n");
    printf (OUT "{\n");
    printf (OUT "  /* Class initialization if needed */\n");
    printf (OUT "  if (!myClassInfo.initialized) ObitTable".$ExtnameClean."ClassInit();\n");
    printf (OUT "\n");
    printf (OUT "  return (gconstpointer)&myClassInfo;\n");
    printf (OUT "} /* end Obit".$ExtnameClean."GetClass */\n");
    printf (OUT "\n");
    printf (OUT "/**\n");
    printf (OUT " * Constructor from values.\n");
    printf (OUT " * Creates a new table structure and attaches to the TableList of file.\n");
    printf (OUT " * If the specified table already exists then it is returned.\n");
    printf (OUT " * Initializes class if needed on first call.\n");
    printf (OUT " * Forces an update of any disk resident structures (e.g. AIPS header).\n");
    printf (OUT " * \x5cparam name   An optional name for the object.\n");
    printf (OUT " * \x5cparam file   ObitData which which the table is to be associated.\n");
    printf (OUT " * \x5cparam ver    Table version number. 0=> add higher, value used returned\n");
    printf (OUT " * \x5cparam access access (OBIT_IO_ReadOnly, means do not create if it doesn't exist.\n");
   # Add structural keywords to call description
    $temp = "		     ";
    for ($k=0; $k<$NumKey; $k++) {
	if ($KeyRang[$k] ne "") {
	    printf (OUT " * \x5cparam $KeyName[$k] $KeyComm[$k]\n");
	}
    }
    printf (OUT " * \x5cparam err Error stack, returns if not empty.\n");
    printf (OUT " * \x5creturn the new object, NULL on failure.\n");
    printf (OUT " */\n");
    printf (OUT "ObitTable".$ExtnameClean."* newObitTable".$ExtnameClean."Value (gchar* name, ObitData *file, olong *ver,\n");
printf (OUT " 	                    ObitIOAccess access,\n");

    # Add structural keywords to call sequence
    $temp = "		    ";
    for ($k=0; $k<$NumKey; $k++) {
	if ($KeyRang[$k] ne "") {
	    $TC = DclTypeCode($KeyType[$k]);
	    $temp = $temp." $TC $KeyName[$k],";
	}
    }
    printf (OUT "  $temp\n");
    printf (OUT "		     ObitErr *err)\n");
    printf (OUT "{\n");
    printf (OUT "  ObitTable".$ExtnameClean."* out=NULL;\n");
    printf (OUT "  ObitTable *testTab=NULL;\n");
    printf (OUT "  ObitTableDesc *desc=NULL;\n");
    printf (OUT "  ObitTableList *list=NULL;\n");
    printf (OUT "  ObitInfoList  *info=NULL;\n");
    printf (OUT "  gboolean exist, optional;\n");
    printf (OUT "  olong colNo, i, ncol, highVer;\n");
    printf (OUT "  ObitIOCode retCode;\n");
    # Drop any "IDI_" prefix
    if ($ExtnameClean=~"IDI_") {
	printf (OUT "  gchar *tabType = \"".$ExtType.substr($Extname,4)."\";\n");
    } else {
	printf (OUT "  gchar *tabType = \"".$ExtType.$Extname."\";\n");
    }
    printf (OUT "  gchar *routine = \"newObitTable".$ExtnameClean."Value\";\n");
    printf (OUT "\n");
    printf (OUT " /* error checks */\n");
    printf (OUT "  g_assert(ObitErrIsA(err));\n");
    printf (OUT "  if (err->error) return NULL;\n");
    printf (OUT "  g_assert (ObitDataIsA(file));\n");
    printf (OUT "\n");
    printf (OUT "  /* Class initialization if needed */\n");
    printf (OUT "  if (!myClassInfo.initialized) ObitTable".$ExtnameClean."ClassInit();\n");
    printf (OUT "\n");
    printf (OUT "  /* Check if the table already exists */\n");
    printf (OUT "  /* Get TableList */\n");
    printf (OUT "  list = ((ObitData*)file)->tableList;\n");
    printf (OUT "  info = ((ObitData*)file)->info;\n");
    printf (OUT "\n");
    printf (OUT "  /* Get highest version number if not specified */\n");
    printf (OUT "  if (*ver==0) { \n");
    # Drop any "IDI_" prefix
    if ($ExtnameClean=~"IDI_") {
	 printf (OUT "    highVer = ObitTableListGetHigh (list, \"".$ExtType.substr($Extname,4)."\");\n");
     } else {
	 printf (OUT "    highVer = ObitTableListGetHigh (list, \"".$ExtType.$Extname."\");\n");
     }
    printf (OUT "    if (access==OBIT_IO_ReadOnly) *ver = highVer;\n");
    printf (OUT "    else if (access==OBIT_IO_ReadWrite) *ver = highVer;\n");
    printf (OUT "    else if (access==OBIT_IO_WriteOnly) *ver = highVer+1;\n");
    printf (OUT "  }\n");
    printf (OUT "  /* See if it already exists */\n");
    printf (OUT "  exist = FALSE;\n");
    printf (OUT "  if (*ver>0) { /* has to be specified */\n");
    printf (OUT "    exist = ObitTableListGet(list, tabType, ver, &testTab, err);\n");
    printf (OUT "    if (err->error) /* add traceback,return */\n");
    printf (OUT "      Obit_traceback_val (err, routine,\"\", out);\n");
    printf (OUT "  \n");
    printf (OUT "    /* if readonly, it must exist to proceed */\n");
    printf (OUT "    if ((access==OBIT_IO_ReadOnly) && !exist) return out;\n");
    printf (OUT "    if (testTab!=NULL) { /* it exists, use it if is an ObitTable".$ExtnameClean." */\n");
    printf (OUT "      if (ObitTable".$ExtnameClean."IsA(testTab)) { /* it is an ObitTable".$ExtnameClean." */\n");
    printf (OUT "	out = ObitTableRef(testTab);\n");
    printf (OUT "      } else { /* needs conversion */\n");
    printf (OUT " 	out = ObitTable".$ExtnameClean."Convert(testTab);\n");
    printf (OUT "	/* Update the TableList */\n");
    printf (OUT "	ObitTableListPut(list, tabType, ver, (ObitTable*)out, err);\n");
    printf (OUT "	if (err->error) /* add traceback,return */\n");
    printf (OUT "	  Obit_traceback_val (err, routine,\"\", out);\n");
    printf (OUT "      }\n");
    printf (OUT "      testTab = ObitTableUnref(testTab); /* remove reference */\n");
    printf (OUT "      return out; /* done */\n");
    printf (OUT "    }\n");
    printf (OUT "  } /* end of test for previously existing table */\n");
    printf (OUT "  \n");
    printf (OUT "  /* If access is ReadOnly make sure one exists */\n");
    printf (OUT "  if (access==OBIT_IO_ReadOnly) { \n");
    # Drop any "IDI_" prefix
    if ($Extname=~"IDI_") {
	printf (OUT "    highVer = ObitTableListGetHigh (list, \"".$ExtType.substr($Extname,4)."\");\n");
     } else {
	printf (OUT "    highVer = ObitTableListGetHigh (list, \"".$ExtType.$Extname."\");\n");
    }
    printf (OUT "    if (highVer<=0) return out;\n");
    printf (OUT "  }\n");
    printf (OUT "  \n");
    printf (OUT "  /* create basal table */\n");
    printf (OUT "  testTab = newObitDataTable ((ObitData*)file, access, tabType,\n");
    printf (OUT "			       ver, err);\n");
    printf (OUT "  if (err->error) Obit_traceback_val (err, routine,\"\", out);\n");
    printf (OUT "  \n");
    printf (OUT "  /* likely need to convert */\n");
    printf (OUT "  if (ObitTable".$ExtnameClean."IsA(testTab)) { \n");
    printf (OUT "    out = ObitTableRef(testTab);\n");
    printf (OUT "  } else { /* needs conversion */\n");
    printf (OUT "    out = ObitTable".$ExtnameClean."Convert(testTab);\n");
    printf (OUT "  }\n");
    printf (OUT "  testTab = ObitTableUnref(testTab); /* remove reference */\n");
    printf (OUT "\n");
    printf (OUT "  /* Update the TableList */\n");
    printf (OUT "  ObitTableListPut(list, tabType, ver, (ObitTable*)out, err);\n");
    printf (OUT "  if (err->error) Obit_traceback_val (err, routine,\"\", out);\n");
    printf (OUT "\n");
    printf (OUT "  /* if it previously existed merely return it */\n");
    printf (OUT "  if (exist) return out; \n");
    printf (OUT "\n");
    printf (OUT "  /* set ClassInfo */\n");
    printf (OUT "  out->ClassInfo = (gpointer)&myClassInfo;\n");
    printf (OUT "\n");
    printf (OUT "  /* Set values */\n");

    # Write structural keywords
    for ($k=0; $k<$NumKey; $k++) {
	if ($KeyRang[$k] ne "") {
	    printf (OUT "  out->$KeyName[$k] = MAX (0, $KeyName[$k]);\n");
	}
    }

    # Write nonstructural keyword defaults
    for ($k=0; $k<$NumKey; $k++) {
	if (($KeyRang[$k] eq "") && ($KeyDefV[$k] ne "")){
	    # Strings aways different
	    if ($KeyType[$k] eq "A") {
		# String
		printf (OUT "  strncpy (out->$KeyName[$k], $KeyDefV[$k], MAXKEYCHARTABLE$ExtnameClean );\n");
	    } else {
		# Numeric/boolean
		printf (OUT "  out->$KeyName[$k] = $KeyDefV[$k];\n");
	    }
	}
    }

    printf (OUT "\n");
    printf (OUT "  /* initialize descriptor */\n");
    printf (OUT "  desc = out->myDesc;\n");
    printf (OUT "  /* How many columns actually in table? */\n");

    # number of columns actually in this table
    $temp =     "ncol = $NumColNoindex";
    for ($k=0; $k<$NumKey; $k++) {
	if ($KeyRang[$k] ne "") {
	    $temp = $temp." + out->$KeyName[$k]*$KeyNumCol[$k]";
	}
    }
    printf (OUT "  $temp ;\n");
    printf (OUT "  desc->FieldName = g_malloc0((ncol+1)*sizeof(gchar*));\n");
    printf (OUT "  desc->FieldUnit = g_malloc0((ncol+1)*sizeof(gchar*));\n");
    printf (OUT "  desc->type      = g_malloc0((ncol+1)*sizeof(ObitInfoType));\n");
    printf (OUT "  desc->dim       = g_malloc0((ncol+1)*sizeof(gint32*));\n");
    printf (OUT "  for (i=0; i<ncol+1; i++) \n");
    printf (OUT "    desc->dim[i] = g_malloc0(MAXINFOELEMDIM*sizeof(gint32));\n");
    printf (OUT "\n");
    printf (OUT "  desc->TableName = g_strdup(tabType);\n");
    printf (OUT "  desc->sort[0] = 0;\n");
    printf (OUT "  desc->sort[1] = 0;\n");
    printf (OUT "  colNo = 0;\n");
    printf (OUT "\n");
    printf (OUT "  /* Define Columns */\n");

    # Define columns - MUST be in same order as AIPS
    $done1 = 0;
    $keylim = $NumKey;
    if ($keylim<=0) {$keylim = 1;}
    # loop over keywords
    for ($k=0; $k<$keylim; $k++) {
	# it must have a range to be useful for indexing
	if ($k<$NumKey) {($lo, $hi) = split(',', substr($KeyRang[$k], 1));}
	else {$lo=0; $hi=0;} # no keywords - still need one pass
	if ((($lo>=1) && ($hi>=$lo)) || ($done1==0)) {
	    # Loop over index values
	    for ($jndex=$lo; $jndex<=$hi; $jndex++) {

		# Loop over columns
		for ($icol=0; $icol<$NumCol; $icol++) {

		    # Unindexed, scalar columns
		    $count = GetColumnCount($ColDim[$icol]);
		    # Need this one here? Only on first keyword loop
		    #print ("DEBUG icol $icol $ColName[$icol] count $count done1 $done1 lo $lo ColSuff .$ColSuff[$icol]. \n");
		    if (($done1==0) && ($count==1) && ($ColSuff[$icol] eq "")) {
			printf (OUT "  desc->FieldName[colNo] = g_strdup($ColName[$icol]);\n");
			printf (OUT "  desc->FieldUnit[colNo] = g_strdup($ColUnit[$icol]);\n");
			$TC = GetTypeCode($ColType[$icol]);
			printf (OUT "  desc->type[colNo] = $TC;\n");
			printf (OUT "  for (i=0; i<MAXINFOELEMDIM; i++) desc->dim[colNo][i] = 1;\n");
			printf (OUT "  colNo++;\n");
		    }

		    # Unindexed, array columns
		    $count = GetColumnCount($ColDim[$icol]);
		    if (($done1==0) && ($count!=1) && ($ColSuff[$icol] eq "")) {
			# parse dimension string
			$temp = substr($ColDim[$icol], 1);
			$temp =~ s/\)//g; # no end parenthese
			($f1,$f2,$f3,$f4,$f5) = split(',',  $temp);
			# Drop optional, zero length columns
			if (index($ColComm[$icol],"[OPTIONAL]")>=0) {
			    printf (OUT "  optional = TRUE;\n");
			} else {
			    printf (OUT "  optional = FALSE;\n");
			}
			printf (OUT "  if (($f1 > 0) || (!optional)) {\n");
			printf (OUT "    desc->FieldName[colNo] = g_strdup($ColName[$icol]);\n");
			printf (OUT "    desc->FieldUnit[colNo] = g_strdup($ColUnit[$icol]);\n");
			$TC = GetTypeCode($ColType[$icol]);
			printf (OUT "    desc->type[colNo] = $TC;\n");
			printf (OUT "    for (i=0; i<MAXINFOELEMDIM; i++) desc->dim[colNo][i] = 1;\n");
			if (($f1 > 1) || ((substr ($f1,0,1) =~ /\D/))) {
			    # Bits arrays are special - packed into oints
			    if ($ColType[$icol] eq "?") {  # but doesn't matter here
				printf (OUT "    desc->dim[colNo][0] = 1 + ($f1 - 1) / sizeof(oint);\n");
			    } else {
				printf (OUT "    desc->dim[colNo][0] = $f1;\n");
			    }
			}
			if (($f2 > 1) || ((substr ($f2,0,1) =~ /\D/))) {
			    printf (OUT "    desc->dim[colNo][1] = $f2;\n");}
			if (($f3 > 1) || ((substr ($f3,0,1) =~ /\D/))) {
			    printf (OUT "    desc->dim[colNo][2] = $f3;\n");}
			if (($f4 > 1) || ((substr ($f4,0,1) =~ /\D/))) {
			    printf (OUT "    desc->dim[colNo][3] = $f4;\n");}
			if (($f5 > 1) || ((substr ($f5,0,1) =~ /\D/))) {
			    printf (OUT "    desc->dim[colNo][4] = $f4;\n");}
			printf (OUT "    colNo++;\n");
			printf (OUT "  }\n");
		    }
		
		    # indexed 
		    $colname =  $ColVName[$icol]."$jndex";
		    $count = GetColumnCount($ColDim[$icol]);
		    #print (" col $icol $colname $ColSuff[$icol] $KeyWord[$k] \n");
		    
		    # do we want this one? indexed by keyword $k
		    $want = $ColSuff[$icol] eq NoQuote($KeyWord[$k]);
		    $want = $want && (($lo>=1) && ($hi>=$lo));
		    
		    if (($count==1) && ($want)) {
			# scalar indexed columns
			$colname = $ColName[$icol]."$jndex\"";
			printf (OUT "  if ($KeyName[$k]>=$jndex) {\n");
			printf (OUT "    desc->FieldName[colNo] = g_strdup($colname);\n");
			printf (OUT "    desc->FieldUnit[colNo] = g_strdup($ColUnit[$icol]);\n");
			$TC = GetTypeCode($ColType[$icol]);
			printf (OUT "    desc->type[colNo] = $TC;\n");
			printf (OUT "    for (i=0; i<MAXINFOELEMDIM; i++) desc->dim[colNo][i] = 1;\n");
			printf (OUT "    colNo++;\n");
			printf (OUT "  }\n");

		    } elsif (($count!=1) && ($want)){
			# array indexed columns
			# parse dimension string
			$temp = substr($ColDim[$icol], 1);
			$temp =~ s/\)//g; # no end parenthese
			($f1,$f2,$f3,$f4,$f5) = split(',',  $temp);
			printf (OUT "  if ($KeyName[$k]>=$jndex) {\n");
			$colname = $ColName[$icol]."$jndex\"";
			printf (OUT "    desc->FieldName[colNo] = g_strdup($colname);\n");
			printf (OUT "    desc->FieldUnit[colNo] = g_strdup($ColUnit[$icol]);\n");
			$TC = GetTypeCode($ColType[$icol]);
			printf (OUT "    desc->type[colNo] = $TC;\n");
			printf (OUT "    for (i=0; i<MAXINFOELEMDIM; i++) desc->dim[colNo][i] = 1;\n");
			if (($f1 > 1) || ((substr ($f1,0,1) =~ /\D/))) {
			    # Bits arrays are special - packed into oints
			    if ($ColType[$k] eq "X") {
				printf (OUT "    desc->dim[colNo][0] = 1 + ($f1 - 1) / sizeof(oint);\n");
			    } else {
				printf (OUT "    desc->dim[colNo][0] = $f1;\n");
			    }
			}
			if (($f2 > 1) || ((substr ($f2,0,1) =~ /\D/)))  {
			    printf (OUT "    desc->dim[colNo][1] = $f2;\n");}
			if (($f3 > 1) || ((substr ($f3,0,1) =~ /\D/)))  {
			    printf (OUT "    desc->dim[colNo][2] = $f3;\n");}
			if (($f4 > 1) || ((substr ($f4,0,1) =~ /\D/)))  {
			    printf (OUT "    desc->dim[colNo][3] = $f4;\n");}
			if (($f5 > 1) || ((substr ($f5,0,1) =~ /\D/)))  {
			    printf (OUT "    desc->dim[colNo][4] = $f4;\n");}
			printf (OUT "    colNo++;\n");
			printf (OUT "  }\n");
		    }

		    
		} # end loop over columns
		$done1 = 1;  # Should have written first set
	    } # end loop over indices
	} # end if keyword has index range
    } # end loop over keywords
    printf (OUT "  /* Add _status column at end */\n");
    printf (OUT "  desc->FieldName[colNo] = g_strdup(\"_status\");\n");
    printf (OUT "  desc->FieldUnit[colNo] = g_strdup(\"        \");\n");
    printf (OUT "  desc->type[colNo] = OBIT_long;\n");
    printf (OUT "  for (i=0; i<MAXINFOELEMDIM; i++) desc->dim[colNo][i] = 1;\n");
    printf (OUT "  \n");
    printf (OUT "  /* number of fields */\n");
    printf (OUT "  desc->nfield = colNo + 1;\n");
    printf (OUT "\n");
    printf (OUT "  /* initialize descriptor keywords */\n");
    printf (OUT "  ObitTable".$ExtnameClean."DumpKey (out, err);\n");
    printf (OUT " \n");
    printf (OUT "  /* index table descriptor */\n");
    printf (OUT "  ObitTableDescIndex (desc);\n");
    printf (OUT "\n");
    printf (OUT "  /* Open and Close to fully instantiate */\n");
    printf (OUT "  retCode = ObitTable".$ExtnameClean."Open(out, OBIT_IO_WriteOnly, err);\n");
    printf (OUT "  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */\n");
    printf (OUT "    Obit_traceback_val (err, routine, out->name, out);    \n");
    printf (OUT "  \n");
    printf (OUT "  retCode = ObitTable".$ExtnameClean."Close(out, err);\n");
    printf (OUT "  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */\n");
    printf (OUT "    Obit_traceback_val (err, routine, out->name, out); \n");
    printf (OUT "\n");
    printf (OUT "  /* Force update of disk resident info */\n");
    printf (OUT "  retCode = ObitIOUpdateTables (((ObitData*)file)->myIO, info, err);\n");
    printf (OUT "  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */\n");
    printf (OUT "    Obit_traceback_val (err, routine, out->name, out); \n");
    printf (OUT "  \n");
    printf (OUT " return out;\n");
    printf (OUT "} /* end newObitTable".$ExtnameClean."Value */\n");
    printf (OUT "\n");
} # end WriteConstructor

# Write copy routines
sub WriteCopy {
    local $k;
    printf (OUT "/**\n");
    printf (OUT " * Convert an ObitTable to an ObitTable".$ExtnameClean.".\n");
    printf (OUT " * New object will have references to members of in.\n");
    printf (OUT " * \x5cparam in  The object to copy, will still exist afterwards \n");
    printf (OUT " *            and should be Unrefed if not needed.\n");
    printf (OUT " * \x5creturn pointer to the new object.\n");
    printf (OUT " */\n");
    printf (OUT "ObitTable".$ExtnameClean."* ObitTable".$ExtnameClean."Convert (ObitTable* in)\n");
    printf (OUT "{\n");
    printf (OUT "  ObitTable".$ExtnameClean." *out;\n");
    printf (OUT "\n");
    printf (OUT "  /* error check */\n");
    printf (OUT "  g_assert(ObitTableIsA(in));\n");
    printf (OUT "\n");
    printf (OUT "  /* create basic object */\n");
    printf (OUT "  out = newObitTable".$ExtnameClean."(in->name);\n");
    printf (OUT "\n");
    printf (OUT "  /* Delete structures on new */\n");
    printf (OUT "  out->info   = ObitInfoListUnref(out->info);\n");
    printf (OUT "  out->thread = ObitThreadUnref(out->thread);\n");
    printf (OUT "  out->myDesc = ObitTableDescUnref(out->myDesc);\n");
    printf (OUT "  out->mySel  = ObitTableSelUnref(out->mySel);\n");
    printf (OUT "  \n");
    printf (OUT "  /* Reference members of in */\n");
    printf (OUT "  \n");
    printf (OUT "  out->info   = ObitInfoListRef(in->info);\n");
    printf (OUT "  out->thread = ObitThreadRef(in->thread);\n");
    printf (OUT "  out->myDesc = ObitTableDescRef(in->myDesc);\n");
    printf (OUT "  out->mySel  = ObitTableSelRef(in->mySel);\n");
    printf (OUT "\n");
    printf (OUT "  /* Remember who I am */\n");
    printf (OUT " out->tabType = g_strdup(in->tabType);\n");
    printf (OUT " out->tabVer  = in->tabVer;\n");
    printf (OUT "  /* Secret reference to host */ \n");
    printf (OUT " out->myHost  = in->myHost;\n");
    printf (OUT "\n");
    printf (OUT "  return out;\n");
    printf (OUT "} /* end ObitTable".$ExtnameClean."Convert */\n");
    printf (OUT "\n");
    printf (OUT "\n");
    printf (OUT "/**\n");
    printf (OUT " * Make a deep copy of input object.\n");
    printf (OUT " * Copies are made of complex members including disk files; these \n");
    printf (OUT " * will be copied applying whatever selection is associated with the input.\n");
    printf (OUT " * Objects should be closed on input and will be closed on output.\n");
    printf (OUT " * In order for the disk file structures to be copied, the output file\n");
    printf (OUT " * must be sufficiently defined that it can be written.\n");
    printf (OUT " * The copy will be attempted but no errors will be logged until\n");
    printf (OUT " * both input and output have been successfully opened.\n");
    printf (OUT " * ObitInfoList and ObitThread members are only copied if the output object\n");
    printf (OUT " * didn't previously exist.\n");
    printf (OUT " * Parent class members are included but any derived class info is ignored.\n");
    printf (OUT " * \x5cparam in  The object to copy\n");
    printf (OUT " * \x5cparam out An existing object pointer for output or NULL if none exists.\n");
    printf (OUT " * \x5cparam err Error stack, returns if not empty.\n");
    printf (OUT " * \x5creturn pointer to the new object.\n");
    printf (OUT " */\n");
    printf (OUT "ObitTable".$ExtnameClean."* ObitTable".$ExtnameClean."Copy (ObitTable".$ExtnameClean." *in, ObitTable".$ExtnameClean." *out, ObitErr *err)\n");
    printf (OUT "{\n");
    printf (OUT "  gchar *routine = \"ObitTable".$ExtnameClean."Copy\";\n");
    printf (OUT "\n");
    printf (OUT "  /* Class initialization if needed */\n");
    printf (OUT "  if (!myClassInfo.initialized) ObitTable".$ExtnameClean."ClassInit();\n");
    printf (OUT "\n");
    printf (OUT " /* error checks */\n");
    printf (OUT "  g_assert(ObitErrIsA(err));\n");
    printf (OUT "  if (err->error) return NULL;\n");
    printf (OUT "  g_assert (ObitIsA(in, &myClassInfo));\n");
    printf (OUT "  if (out) g_assert (ObitIsA(out, &myClassInfo));\n");
    printf (OUT "\n");
    printf (OUT "  /* Use parent class to copy */\n");
    printf (OUT "  out = (ObitTable".$ExtnameClean."*)ObitTableCopy ((ObitTable*)in, (ObitTable*)out, err);\n");
    printf (OUT "  if (err->error) /* add traceback,return */\n");
    printf (OUT "    Obit_traceback_val (err, routine,in->name, out);\n");
    printf (OUT "\n");
    printf (OUT "  /* Copy this class  info */\n");

    # Copy keywords
    for ($k=0; $k<$NumKey; $k++) {
	if ($KeyType[$k] eq "A") { # Strings more difficult
	    printf (OUT "  strncpy (out->$KeyName[$k], in->$KeyName[$k], MAXKEYCHARTABLE$ExtnameClean );\n");
	} else {
	    printf (OUT "  out->$KeyName[$k] = in->$KeyName[$k];\n");
	}
    }

    printf (OUT "  /* Update class specific info */\n");
    printf (OUT "  ObitTable".$ExtnameClean."Update (out, err);\n");
    printf (OUT "    \n");
    printf (OUT "  return out;\n");
    printf (OUT "} /* end ObitTable".$ExtnameClean."Copy */\n");
    printf (OUT "\n");
} # end WriteCopy

# Write Open
sub WriteOpen {
    printf (OUT "/**\n");
    printf (OUT " * Initialize structures and open file.\n");
    printf (OUT " * The image descriptor is read if OBIT_IO_ReadOnly or \n");
    printf (OUT " * OBIT_IO_ReadWrite and written to disk if opened OBIT_IO_WriteOnly.\n");
    printf (OUT " * After the file has been opened the member, buffer is initialized\n");
    printf (OUT " * for reading/storing the table unless member bufferSize is <0.\n");
    printf (OUT " * If the requested version (\"Ver\" in InfoList) is 0 then the highest\n");
    printf (OUT " * numbered table of the same type is opened on Read or Read/Write, \n");
    printf (OUT " * or a new table is created on on Write.\n");
    printf (OUT " * The file etc. info should have been stored in the ObitInfoList:\n");
    if ($ExtType eq "AIPS ") {
      printf (OUT " * \x5cli \"FileType\" OBIT_long scalar = OBIT_IO_FITS or OBIT_IO_AIPS \n");
      printf (OUT " *               for file type (see class documentation for details).\n");
    } else {
      printf (OUT " * \x5cli \"FileType\" OBIT_long scalar = OBIT_IO_FITS (AIPS not allowed) \n");
      printf (OUT " *               for file type (see class documentation for details).\n");
    }
    printf (OUT " * \x5cli \"nRowPIO\" OBIT_long scalar = Maximum number of table rows\n");
    printf (OUT " *               per transfer, this is the target size for Reads (may be \n");
    printf (OUT " *               fewer) and is used to create buffers.\n");
    printf (OUT " * \x5cparam in Pointer to object to be opened.\n");
    printf (OUT " * \x5cparam access access (OBIT_IO_ReadOnly,OBIT_IO_ReadWrite,\n");
    printf (OUT " *               or OBIT_IO_WriteOnly).\n");
    printf (OUT " *               If OBIT_IO_WriteOnly any existing data in the output file\n");
    printf (OUT " *               will be lost.\n");
    printf (OUT " * \x5cparam err ObitErr for reporting errors.\n");
    printf (OUT " * \x5creturn return code, OBIT_IO_OK=> OK\n");
    printf (OUT " */\n");
    printf (OUT "ObitIOCode ObitTable".$ExtnameClean."Open (ObitTable".$ExtnameClean." *in, ObitIOAccess access, \n");
    printf (OUT "			  ObitErr *err)\n");
    printf (OUT "{\n");
    printf (OUT "  ObitIOCode retCode = OBIT_IO_SpecErr;\n");
    printf (OUT "  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};\n");
    printf (OUT "  olong nRowPIO;\n");
    printf (OUT "  gchar *routine = \"ObitTable".$ExtnameClean."Open\";\n");
    printf (OUT "\n");
    printf (OUT "  /* Class initialization if needed */\n");
    printf (OUT "  if (!myClassInfo.initialized) ObitTable".$ExtnameClean."ClassInit();\n");
    printf (OUT "\n");
    printf (OUT "  /* error checks */\n");
    printf (OUT "  g_assert (ObitErrIsA(err));\n");
    printf (OUT "  if (err->error) return retCode;\n");
    printf (OUT "  g_assert (ObitIsA(in, &myClassInfo));\n");
    printf (OUT "\n");
    printf (OUT "   /* Do one row at a time */\n");
    printf (OUT "   nRowPIO = 1;\n");
    printf (OUT "   ObitInfoListPut(in->info, \"nRowPIO\", OBIT_long, dim, (gconstpointer)&nRowPIO, err);\n");
    printf (OUT "   if (err->error) /* add traceback,return */\n");
    printf (OUT "     Obit_traceback_val (err, routine, in->name, retCode);\n");
    printf (OUT "   \n");
    printf (OUT "   /* use parent class open */\n");
    printf (OUT "   retCode = ObitTableOpen ((ObitTable*)in, access, err);\n");
    printf (OUT "   if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */\n");
    printf (OUT "     Obit_traceback_val (err, routine, in->name, retCode);\n");
    printf (OUT "   \n");
    printf (OUT "   /* Update class specific info */\n");
    printf (OUT "   ObitTable".$ExtnameClean."Update (in, err);\n");
    printf (OUT "   \n");
    printf (OUT "   return retCode;\n");
    printf (OUT "} /* end ObitTable".$ExtnameClean."Open */\n");
    printf (OUT "\n");
} # end WriteOpen

sub WriteReadRow {
    local $k, $temp, $TC, $want, $count, $f1, $f2, $f3, $f4, $f5, $colname;
    printf (OUT "/**\n");
    printf (OUT " * Read a table row and return an easily digested version.\n");
    printf (OUT " * Scalar values are copied but for array values, pointers \n");
    printf (OUT " * into the data array are returned.\n");
    printf (OUT " * \x5cparam in       Table to read\n");
    printf (OUT " * \x5cparam i".$ExtnameClean."Row   Row number, -1 -> next\n");
    printf (OUT " * \x5cparam row      Table Row structure to receive data\n");
    printf (OUT " * \x5cparam err ObitErr for reporting errors.\n");
    printf (OUT " * \x5creturn return code, OBIT_IO_OK=> OK\n");
    printf (OUT " */\n");
    printf (OUT "ObitIOCode \n");
    printf (OUT "ObitTable".$ExtnameClean."ReadRow  (ObitTable".$ExtnameClean." *in, olong i".$ExtnameClean."Row, ObitTable".$ExtnameClean."Row *row,\n");
    printf (OUT "		     ObitErr *err)\n");
    printf (OUT "{\n");
    printf (OUT "  ObitIOCode retCode = OBIT_IO_SpecErr;\n");
    printf (OUT "  odouble   *dRow;\n");
    printf (OUT "  oint      *iRow;\n");
    printf (OUT "  gshort    *siRow;\n");
    printf (OUT "  ofloat    *fRow;\n");
    printf (OUT "  gchar     *cRow;\n");
    printf (OUT "  gboolean  *lRow;\n");
    printf (OUT "  guint8    *bRow;\n");
    printf (OUT "  gchar *routine = \"ObitTable".$ExtnameClean."ReadRow\";\n");
    printf (OUT "  \n");
    printf (OUT "  /* error checks */\n");
    printf (OUT "  g_assert (ObitErrIsA(err));\n");
    printf (OUT "  if (err->error) return retCode;\n");
    printf (OUT "  g_assert (ObitIsA(in, &myClassInfo));\n");
    printf (OUT "\n");
    printf (OUT "  if (in->myStatus == OBIT_Inactive) {\n");
    printf (OUT "    Obit_log_error(err, OBIT_Error,\n");
    printf (OUT "		   \"".$ExtType.$ExtnameClean." Table is inactive for  %%s \", in->name);\n");
    printf (OUT "    return retCode;\n");
    printf (OUT " }\n");
    printf (OUT "\n");
    printf (OUT "  /* read row i".$ExtnameClean."Row */\n");
    printf (OUT "  retCode = ObitTableRead ((ObitTable*)in, i".$ExtnameClean."Row, NULL,  err);\n");
    printf (OUT "  if (err->error) \n");
    printf (OUT "    Obit_traceback_val (err, routine, in->name, retCode);\n");
    printf (OUT "\n");
    printf (OUT "  /* Typed pointers to row of data */  \n");
    printf (OUT "  dRow  = (odouble*)in->buffer;\n");
    printf (OUT "  iRow  = (oint*)in->buffer;\n");
    printf (OUT "  siRow = (gshort*)in->buffer;\n");
    printf (OUT "  fRow  = (ofloat*)in->buffer;\n");
    printf (OUT "  cRow  = (gchar*)in->buffer;\n");
    printf (OUT "  lRow  = (gboolean*)in->buffer;\n");
    printf (OUT "  bRow  = (guint8*)in->buffer;\n");
    printf (OUT "  \n");
    printf (OUT "  /* Copy scalar fields, for arrays only set pointer*/\n");
    # Unindexed, scalar columns
    for ($k=0; $k<$NumCol; $k++) {
	$count = GetColumnCount($ColDim[$k]);
	if (($count==1) && ($ColSuff[$k] eq "")) {
	    $TC = RowType($ColType[$k])."[in->$ColVName[$k]Off]";
	    printf (OUT "  row->$ColVName[$k] = $TC;\n");
	}
    }

    # Unindexed, array columns
    for ($k=0; $k<$NumCol; $k++) {
	$count = GetColumnCount($ColDim[$k]);
	if (($count!=1) && ($ColSuff[$k] eq "")) {
	    $TC = RowType($ColType[$k]);
	    printf (OUT "  row->$ColVName[$k] = $TC + in->$ColVName[$k]Off;\n");
	}
    }

    # indexed - loop over keywords
    for ($k=0; $k<$NumKey; $k++) {
	# it must have a range to be useful
	($lo, $hi) = split(',', substr($KeyRang[$k], 1));
	if (($lo>=1) && ($hi>=$lo)) {
	    # Loop over index values
	    for ($j=$lo; $j<=$hi; $j++) {

		# Loop over columns
		for ($i=0; $i<$NumCol; $i++) {
		    $colname =  $ColVName[$i]."$j";
		    $count = GetColumnCount($ColDim[$i]);
		    #print (" col $i $colname $ColSuff[$i] $KeyWord[$k] \n");

		    # do we want this one? indexed by keyword $k
		    $want = $ColSuff[$i] eq NoQuote($KeyWord[$k]);
		    
		    if (($count==1) && ($want)) {
			# scalar indexed columns
			$colname = $ColVName[$i]."$j";
			$TC = RowType($ColType[$i])."[in->$colname\Off]";
			printf (OUT "  if (in->$KeyName[$k]>=$j) row->$colname = $TC;\n");
		    } elsif (($count!=1) && ($want)){
			# array indexed columns
			$colname = $ColVName[$i]."$j";
			$TC = RowType($ColType[$i]);
			printf (OUT "  if (in->$KeyName[$k]>=$j) row->$colname = $TC + in->$colname\Off;\n");
		    }
		} # end loop over columns
		
	    } # end loop over indices
	}
    } # end loop over keywords
    # _status column
    printf (OUT "  row->status = iRow[in->myDesc->statusOff];\n");
    printf (OUT "\n");
    printf (OUT "  return retCode;\n");
    printf (OUT "} /*  end ObitTable".$ExtnameClean."ReadRow */\n");
    printf (OUT "\n");
} # end WriteReadRow

sub WriteSetRow {
    local $k, $temp, $TC, $want, $count, $f1, $f2, $f3, $f4, $f5, $colname;
    printf (OUT "/**\n");
    printf (OUT " * Attach an ObitTableRow to the buffer of an ObitTable.\n");
    printf (OUT " * This is only useful prior to filling a row structure in preparation .\n");
    printf (OUT " * for a WriteRow operation.  Array members of the Row structure are .\n");
    printf (OUT " * pointers to independently allocated memory, this routine allows using .\n");
    printf (OUT " * the table IO buffer instead of allocating yet more memory..\n");
    printf (OUT " * This routine need only be called once to initialize a Row structure for write..\n");
    printf (OUT " * \x5cparam in  Table with buffer to be written \n");
    printf (OUT " * \x5cparam row Table Row structure to attach \n");
    printf (OUT " * \x5cparam err ObitErr for reporting errors.\n");
    printf (OUT " */\n");
    printf (OUT "void \n");
    printf (OUT "ObitTable".$ExtnameClean."SetRow  (ObitTable".$ExtnameClean." *in, ObitTable".$ExtnameClean."Row *row,\n");
    printf (OUT "		     ObitErr *err)\n");
    printf (OUT "{\n");
    printf (OUT "  odouble   *dRow;\n");
    printf (OUT "  oint      *iRow;\n");
    printf (OUT "  gshort    *siRow;\n");
    printf (OUT "  ofloat    *fRow;\n");
    printf (OUT "  gchar     *cRow;\n");
    printf (OUT "  gboolean  *lRow;\n");
    printf (OUT "  guint8    *bRow;\n");
    printf (OUT "  \n");
    printf (OUT "  /* error checks */\n");
    printf (OUT "  g_assert (ObitErrIsA(err));\n");
    printf (OUT "  if (err->error) return;\n");
    printf (OUT "  g_assert (ObitIsA(in, &myClassInfo));\n");
    printf (OUT "  g_assert (ObitIsA(row, &myRowClassInfo));\n");
    printf (OUT "\n");
    printf (OUT "  if (in->myStatus == OBIT_Inactive) {\n");
    printf (OUT "    Obit_log_error(err, OBIT_Error,\n");
    printf (OUT "		   \"".$ExtnameClean." Table is inactive for  %%s \", in->name);\n");
    printf (OUT "    return;\n");
    printf (OUT " }\n");
    printf (OUT "\n");
    printf (OUT "  /* Typed pointers to row of data */  \n");
    printf (OUT "  dRow  = (odouble*)in->buffer;\n");
    printf (OUT "  iRow  = (oint*)in->buffer;\n");
    printf (OUT "  siRow = (gshort*)in->buffer;\n");
    printf (OUT "  fRow  = (ofloat*)in->buffer;\n");
    printf (OUT "  cRow  = (gchar*)in->buffer;\n");
    printf (OUT "  lRow  = (gboolean*)in->buffer;\n");
    printf (OUT "  bRow  = (guint8*)in->buffer;\n");
    printf (OUT "  \n");
    printf (OUT "  /* Set row pointers to buffer */\n");
    # Unindexed, array columns
    for ($k=0; $k<$NumCol; $k++) {
	$count = GetColumnCount($ColDim[$k]);
	if (($count!=1) && ($ColSuff[$k] eq "")) {
	    $TC = RowType($ColType[$k]);
	    printf (OUT "  row->$ColVName[$k] = $TC + in->$ColVName[$k]Off;\n");
	}
    }

    # indexed - loop over keywords
    for ($k=0; $k<$NumKey; $k++) {
	# it must have a range to be useful
	($lo, $hi) = split(',', substr($KeyRang[$k], 1));
	if (($lo>=1) && ($hi>=$lo)) {
	    # Loop over index values
	    for ($j=$lo; $j<=$hi; $j++) {

		# Loop over columns
		for ($i=0; $i<$NumCol; $i++) {
		    $colname =  $ColVName[$i]."$j";
		    $count = GetColumnCount($ColDim[$i]);
		    #print (" col $i $colname $ColSuff[$i] $KeyWord[$k] \n");

		    # do we want this one? indexed by keyword $k
		    $want = $ColSuff[$i] eq NoQuote($KeyWord[$k]);
		    
		    if (($count==1) && ($want)) {
			# scalar indexed columns
			$colname = $ColVName[$i]."$j";
			$TC = RowType($ColType[$i])."[in->$colname\Off]";
			printf (OUT "  if (in->$KeyName[$k]>=$j) row->$colname = $TC;\n");
		    } elsif (($count!=1) && ($want)){
			# array indexed columns
			$colname = $ColVName[$i]."$j";
			$TC = RowType($ColType[$i]);
			printf (OUT "  if (in->$KeyName[$k]>=$j) row->$colname = $TC + in->$colname\Off;\n");
		    }
		} # end loop over columns
		
	    } # end loop over indices
	}
    } # end loop over keywords
    printf (OUT "\n");
    printf (OUT "} /*  end ObitTable".$ExtnameClean."SetRow */\n");
    printf (OUT "\n");
} # end WriteSetRow

# Write WriteRow
sub WriteWriteRow {
    local $k, $temp, $TC, $want, $count, $f1, $f2, $f3, $f4, $f5, $colname;
    printf (OUT "/**\n");
    printf (OUT " * Write a table row.\n");
    printf (OUT " * Before calling this routine, the row structure needs to be initialized\n");
    printf (OUT " * and filled with data. The array members of the row structure are  \n");
    printf (OUT " * pointers to independently allocated memory.  These pointers can be set to the \n");
    printf (OUT " * correct table buffer locations using ObitTable".$ExtnameClean."SetRow  \n");
    printf (OUT " * \x5cparam in       Table to read\n");
    printf (OUT " * \x5cparam i".$ExtnameClean."Row   Row number, -1 -> next\n");
    printf (OUT " * \x5cparam row Table Row structure containing data\n");
    printf (OUT " * \x5cparam err ObitErr for reporting errors.\n");
    printf (OUT " * \x5creturn return code, OBIT_IO_OK=> OK\n");
    printf (OUT " */\n");
    printf (OUT "ObitIOCode \n");
    printf (OUT "ObitTable".$ExtnameClean."WriteRow  (ObitTable".$ExtnameClean." *in, olong i".$ExtnameClean."Row, ObitTable".$ExtnameClean."Row *row,\n");
    printf (OUT "		      ObitErr *err)\n");
    printf (OUT "{\n");
    printf (OUT "  ObitIOCode retCode = OBIT_IO_SpecErr;\n");
    printf (OUT "  gshort    *siRow;\n");
    printf (OUT "  odouble   *dRow;\n");
    printf (OUT "  oint      *iRow, i;\n");
    printf (OUT "  ofloat    *fRow;\n");
    printf (OUT "  gchar     *cRow;\n");
    printf (OUT "  gboolean  *lRow;\n");
    printf (OUT "  guint8    *bRow;\n");
    printf (OUT "  gchar *routine = \"ObitTable".$ExtnameClean."WriteRow\";\n");
    printf (OUT "  \n");
    printf (OUT "\n");
    printf (OUT "  /* error checks */\n");
    printf (OUT "  g_assert (ObitErrIsA(err));\n");
    printf (OUT "  if (err->error) return retCode;\n");
    printf (OUT "  g_assert (ObitIsA(in, &myClassInfo));\n");
    printf (OUT "\n");
    printf (OUT "  if (in->myStatus == OBIT_Inactive) {\n");
    printf (OUT "    Obit_log_error(err, OBIT_Error,\n");
    printf (OUT "		   \"".$ExtType.$ExtnameClean." Table is inactive for %%s \", in->name);\n");
    printf (OUT "    return retCode;\n");
    printf (OUT " }\n");
    printf (OUT "\n");
    printf (OUT "  /* Typed pointers to row of data */  \n");
    printf (OUT "  dRow  = (odouble*)in->buffer;\n");
    printf (OUT "  siRow = (gshort*)in->buffer;\n");
    printf (OUT "  iRow  = (oint*)in->buffer;\n");
    printf (OUT "  fRow  = (ofloat*)in->buffer;\n");
    printf (OUT "  cRow  = (gchar*)in->buffer;\n");
    printf (OUT "  lRow  = (gboolean*)in->buffer;\n");
    printf (OUT "  bRow  = (guint8*)in->buffer;\n");
    printf (OUT "  \n");
    printf (OUT "  /* Make full copy of all data */\n");
    # Unindexed, scalar columns
    for ($k=0; $k<$NumCol; $k++) {
	$count = GetColumnCount($ColDim[$k]);
	if (($count==1) && ($ColSuff[$k] eq "")) {
	    $TC = RowType($ColType[$k])."[in->$ColVName[$k]Off]";
	    printf (OUT "  $TC = row->$ColVName[$k];\n");
	}
    }

    # Unindexed, array columns
    for ($k=0; $k<$NumCol; $k++) {
	$count = GetColumnCount($ColDim[$k]);
	if (($count!=1) && ($ColSuff[$k] eq "")) {
	    $TC = RowType($ColType[$k])."[in->$ColVName[$k]Off+i]";
	    printf (OUT "  if (in->$ColVName[$k]Col >= 0) { \n");
	    printf (OUT "    for (i=0; i<in->myDesc->repeat[in->$ColVName[$k]Col]; i++) \n");
	    printf (OUT "      $TC = row->$ColVName[$k]\[i];\n");
	    printf (OUT "  } \n");
	}
    }

    # indexed - loop over keywords
    for ($k=0; $k<$NumKey; $k++) {
	# it must have a range to be useful
	($lo, $hi) = split(',', substr($KeyRang[$k], 1));
	if (($lo>=1) && ($hi>=$lo)) {
	    # Loop over index values
	    for ($j=$lo; $j<=$hi; $j++) {

		# Loop over columns
		for ($i=0; $i<$NumCol; $i++) {
		    $colname =  $ColVName[$i]."$j";
		    $count = GetColumnCount($ColDim[$i]);
		    #print (" col $i $colname $ColSuff[$i] $KeyWord[$k] \n");

		    # do we want this one? indexed by keyword $k
		    $want = $ColSuff[$i] eq NoQuote($KeyWord[$k]);
		    
		    if (($count==1) && ($want)) {
			# scalar indexed columns
			$colname = $ColVName[$i]."$j";
			$TC = RowType($ColType[$i])."[in->$colname\Off]";
			printf (OUT "  if (in->$KeyName[$k]>=$j) $TC = row->$colname;\n");
		    } elsif (($count!=1) && ($want)){
			# array indexed columns
			$colname = $ColVName[$i]."$j";
			printf (OUT "  if (in->$KeyName[$k]>=$j)\n");
			$TC = RowType($ColType[$i])."[in->$colname\Off+i]";
			printf (OUT "    if (in->$colname\Col >= 0) { \n");
			printf (OUT "      for (i=0; i<in->myDesc->repeat[in->$colname\Col]; i++) \n");
			printf (OUT "        $TC = row->$colname\[i];\n");
			printf (OUT "    } \n");
		    }
		} # end loop over columns
		
	    } # end loop over indices
	}
    } # end loop over keywords
    printf (OUT "\n");
    printf (OUT "  /* copy status */\n");
    printf (OUT "  iRow[in->myDesc->statusOff] = row->status;\n");
    printf (OUT "   \n");
    printf (OUT "  /* Write one row */\n");
    printf (OUT "  in->myDesc->numRowBuff = 1;\n");
    printf (OUT " \n");
    printf (OUT "  /* Write row i".$ExtnameClean."Row */\n");
    printf (OUT "  retCode = ObitTableWrite ((ObitTable*)in, i".$ExtnameClean."Row, NULL,  err);\n");
    printf (OUT "  if (err->error) \n");
    printf (OUT "    Obit_traceback_val (err, routine,in->name, retCode);\n");
    printf (OUT "\n");
    printf (OUT "  return retCode;\n");
    printf (OUT "} /*  end ObitTable".$ExtnameClean."WriteRow */\n");
    printf (OUT "\n");
    printf (OUT "/**\n");
    printf (OUT " * Shutdown I/O.\n");
    printf (OUT " * \x5cparam in Pointer to object to be closed.\n");
    printf (OUT " * \x5cparam err ObitErr for reporting errors.\n");
    printf (OUT " * \x5creturn error code, OBIT_IO_OK=> OK\n");
    printf (OUT " */\n");
    printf (OUT "ObitIOCode ObitTable".$ExtnameClean."Close (ObitTable".$ExtnameClean." *in, ObitErr *err)\n");
    printf (OUT "{\n");
    printf (OUT "  ObitIOCode retCode = OBIT_IO_SpecErr;\n");
    printf (OUT "  gchar *routine = \"ObitTable".$ExtnameClean."Close\";\n");
    printf (OUT "\n");
    printf (OUT "  /* error checks */\n");
    printf (OUT "  g_assert (ObitErrIsA(err));\n");
    printf (OUT "  if (err->error) return retCode;\n");
    printf (OUT "  g_assert (ObitIsA((Obit*)in, &myClassInfo));\n");
    printf (OUT "  /* Something going on? */\n");
    printf (OUT "  if (in->myStatus == OBIT_Inactive) return OBIT_IO_OK;\n");
    printf (OUT "\n");
    printf (OUT "  /* Update keywords on descriptor if not ReadOnly*/\n");
    printf (OUT "  if (in->myDesc->access != OBIT_IO_ReadOnly) \n");
    printf (OUT "    ObitTable".$ExtnameClean."DumpKey (in, err);\n");
    printf (OUT "  if (err->error) \n");
    printf (OUT "    Obit_traceback_val (err, routine, in->name, retCode);\n");
    printf (OUT "\n");
    printf (OUT "  /* Close */\n");
    printf (OUT "  retCode = ObitTableClose ((ObitTable*)in, err);\n");
    printf (OUT "  if (err->error) \n");
    printf (OUT "    Obit_traceback_val (err, routine, in->name, retCode);\n");
    printf (OUT "\n");
    printf (OUT "  return retCode;\n");
    printf (OUT "} /* end ObitTable".$ExtnameClean."Close */\n");
    printf (OUT "\n");
} # end WriteWriteRow

# Write TableRow private routines
sub WriteTableRowPrivate {
    local $k, $temp, $TC, $want, $count, $f1, $f2, $f3, $f4, $f5, $colname;
    printf (OUT "/*---------------Private functions--------------------------*/\n");
    printf (OUT "/*----------------  Table".$ExtnameClean." Row  ----------------------*/\n");
    printf (OUT "/**\n");
    printf (OUT " * Creates empty member objects, initialize reference count.\n");
    printf (OUT " * Parent classes portions are (recursively) initialized first\n");
    printf (OUT " * \x5cparam inn Pointer to the object to initialize.\n");
    printf (OUT " */\n");
    printf (OUT "void ObitTable".$ExtnameClean."RowInit  (gpointer inn)\n");
    printf (OUT "{\n");
    printf (OUT "  ObitClassInfo *ParentClass;\n");
    printf (OUT "  ObitTable".$ExtnameClean."Row *in = inn;\n");
    printf (OUT "\n");
    printf (OUT "  /* error checks */\n");
    printf (OUT "  g_assert (in != NULL);\n");
    printf (OUT "\n");
    printf (OUT "  /* recursively initialize parent class members */\n");
    printf (OUT "  ParentClass = (ObitClassInfo*)(myRowClassInfo.ParentClass);\n");
    printf (OUT "  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) \n");
    printf (OUT "    ParentClass->ObitInit (inn);\n");
    printf (OUT "\n");
    printf (OUT "  /* set members in this class */\n");
    printf (OUT "  /* Set array members to NULL */\n");
    # Unindexed, array columns
    for ($k=0; $k<$NumCol; $k++) {
	$count = GetColumnCount($ColDim[$k]);
	if (($count!=1) && ($ColSuff[$k] eq "")) {
	    printf (OUT "  in->$ColVName[$k] = NULL;\n");
	}
    }

    # indexed - loop over keywords
    for ($k=0; $k<$NumKey; $k++) {
	# it must have a range to be useful
	($lo, $hi) = split(',', substr($KeyRang[$k], 1));
	if (($lo>=1) && ($hi>=$lo)) {
	    # Loop over index values
	    for ($j=$lo; $j<=$hi; $j++) {

		# Loop over columns
		for ($i=0; $i<$NumCol; $i++) {
		    $colname =  $ColVName[$i]."$j";
		    $count = GetColumnCount($ColDim[$i]);
		    #print (" col $i $colname $ColSuff[$i] $KeyWord[$k] \n");

		    # do we want this one? indexed by keyword $k
		    $want = $ColSuff[$i] eq NoQuote($KeyWord[$k]);
		    
		    if (($count!=1) && ($want)){
			# array indexed columns
			$colname = $ColVName[$i]."$j";
			printf (OUT "  in->$colname = NULL;\n");
		    }
		} # end loop over columns
		
	    } # end loop over indices
	}
    } # end loop over keywords
    printf (OUT "\n");
    printf (OUT "} /* end ObitTable".$ExtnameClean."RowInit */\n");
    printf (OUT "\n");
    printf (OUT "/**\n");
    printf (OUT " * Deallocates member objects.\n");
    printf (OUT " * Does (recursive) deallocation of parent class members.\n");
    printf (OUT " * For some reason this wasn't build into the GType class.\n");
    printf (OUT " * \x5cparam  inn Pointer to the object to deallocate.\n");
    printf (OUT " *           Actually it should be an ObitTable".$ExtnameClean."Row* cast to an Obit*.\n");
    printf (OUT " */\n");
    printf (OUT "void ObitTable".$ExtnameClean."RowClear (gpointer inn)\n");
    printf (OUT "{\n");
    printf (OUT "  ObitClassInfo *ParentClass;\n");
    printf (OUT "  ObitTable".$ExtnameClean."Row *in = inn;\n");
    printf (OUT "\n");
    printf (OUT "  /* error checks */\n");
    printf (OUT "  g_assert (ObitIsA(in, &myRowClassInfo));\n");
    printf (OUT "\n");
    printf (OUT "  /* delete this class members */\n");
    printf (OUT "  /* Do not free data array pointers as they were not malloced */\n");
    printf (OUT "  \n");
    printf (OUT "  /* unlink parent class members */\n");
    printf (OUT "  ParentClass = (ObitClassInfo*)(myRowClassInfo.ParentClass);\n");
    printf (OUT "  /* delete parent class members */\n");
    printf (OUT "  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) \n");
    printf (OUT "    ParentClass->ObitClear (inn);\n");
    printf (OUT "  \n");
    printf (OUT "} /* end ObitTable".$ExtnameClean."RowClear */\n");
    printf (OUT "\n");

    printf (OUT "/**\n");
    printf (OUT " * Initialize global ClassInfo Structure.\n");
    printf (OUT " */\n");
    printf (OUT "void ObitTable".$ExtnameClean."RowClassInit (void)\n");
    printf (OUT "{\n");
    printf (OUT "  if (myRowClassInfo.initialized) return;  /* only once */\n");
    printf (OUT "  \n");
    printf (OUT "  /* Set name and parent for this class */\n");
    printf (OUT "  myRowClassInfo.ClassName   = g_strdup(myRowClassName);\n");
    printf (OUT "  myRowClassInfo.ParentClass = ObitParentGetRowClass();\n");
    printf (OUT "\n");
    printf (OUT "  /* Set function pointers */\n");
    printf (OUT "  ObitTable".$ExtnameClean."RowClassInfoDefFn ((gpointer)&myRowClassInfo);\n");
    printf (OUT " \n");
    printf (OUT "  myRowClassInfo.initialized = TRUE; /* Now initialized */\n");
    printf (OUT " \n");
    printf (OUT "} /* end ObitTable".$ExtnameClean."RowClassInit */\n");
    printf (OUT "\n");
    printf (OUT "/**\n");
    printf (OUT " * Initialize global ClassInfo Function pointers.\n");
    printf (OUT " */\n");
    printf (OUT "static void ObitTable".$ExtnameClean."RowClassInfoDefFn (gpointer inClass)\n");
    printf (OUT "{\n");
    printf (OUT "  ObitTable".$ExtnameClean."RowClassInfo *theClass = (ObitTable".$ExtnameClean."RowClassInfo*)inClass;\n");
    printf (OUT "  ObitClassInfo *ParentClass = (ObitClassInfo*)myRowClassInfo.ParentClass;\n");
    printf (OUT "\n");
    printf (OUT "  if (theClass->initialized) return;  /* only once */\n");
    printf (OUT "\n");
    printf (OUT "  /* Check type of inClass */\n");
    printf (OUT "  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myRowClassInfo));\n");
    printf (OUT "\n");
    printf (OUT "  /* Initialize (recursively) parent class first */\n");
    printf (OUT "  if ((ParentClass!=NULL) && \n");
    printf (OUT "      (ParentClass->ObitClassInfoDefFn!=NULL))\n");
    printf (OUT "    ParentClass->ObitClassInfoDefFn(theClass);\n");
    printf (OUT "\n");
    printf (OUT "  /* function pointers defined or overloaded this class */\n");
    printf (OUT "  theClass->ObitClassInit = (ObitClassInitFP)ObitTable".$ExtnameClean."RowClassInit;\n");
    printf (OUT "  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitTable".$ExtnameClean."RowClassInfoDefFn;\n");
    printf (OUT "  theClass->ObitGetClass  = (ObitGetClassFP)ObitTable".$ExtnameClean."RowGetClass;\n");
    printf (OUT "  theClass->newObit         = NULL;\n");
    printf (OUT "  theClass->newObitTableRow = (newObitTableRowFP)newObitTable".$ExtnameClean."Row;\n");
    printf (OUT "  theClass->ObitCopy        = NULL;\n");
    printf (OUT "  theClass->ObitClone       = NULL;\n");
    printf (OUT "  theClass->ObitClear       = (ObitClearFP)ObitTable".$ExtnameClean."RowClear;\n");
    printf (OUT "  theClass->ObitInit        = (ObitInitFP)ObitTable".$ExtnameClean."RowInit;\n");
    printf (OUT "\n");
    printf (OUT "} /* end ObitTable".$ExtnameClean."RowClassDefFn */\n");
    printf (OUT "\n");
} # end WriteTableRowPrivate

# Write Table Private routines
sub WriteTablePrivate {
    local $k, $temp, $TC, $want, $count, $sl, $colname;
    printf (OUT "/*------------------  Table".$ExtnameClean."  ------------------------*/\n");
    printf (OUT "\n");
    printf (OUT "/**\n");
    printf (OUT " * Creates empty member objects, initialize reference count.\n");
    printf (OUT " * Parent classes portions are (recursively) initialized first\n");
    printf (OUT " * \x5cparam inn Pointer to the object to initialize.\n");
    printf (OUT " */\n");
    printf (OUT "void ObitTable".$ExtnameClean."Init  (gpointer inn)\n");
    printf (OUT "{\n");
    printf (OUT "  ObitClassInfo *ParentClass;\n");
    printf (OUT "  ObitTable".$ExtnameClean." *in = inn;\n");
    printf (OUT "\n");
    printf (OUT "  /* error checks */\n");
    printf (OUT "  g_assert (in != NULL);\n");
    printf (OUT "\n");
    printf (OUT "  /* recursively initialize parent class members */\n");
    printf (OUT "  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);\n");
    printf (OUT "  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) \n");
    printf (OUT "    ParentClass->ObitInit (inn);\n");
    printf (OUT "\n");
    printf (OUT "  /* set members in this class */\n");
    printf (OUT "\n");
    printf (OUT "} /* end ObitTable".$ExtnameClean."Init */\n");
    printf (OUT "\n");
    printf (OUT "/**\n");
    printf (OUT " * Deallocates member objects.\n");
    printf (OUT " * Does (recursive) deallocation of parent class members.\n");
    printf (OUT " * For some reason this wasn't build into the GType class.\n");
    printf (OUT " * \x5cparam  inn Pointer to the object to deallocate.\n");
    printf (OUT " *           Actually it should be an ObitTable".$ExtnameClean."* cast to an Obit*.\n");
    printf (OUT " */\n");
    printf (OUT "void ObitTable".$ExtnameClean."Clear (gpointer inn)\n");
    printf (OUT "{\n");
    printf (OUT "  ObitClassInfo *ParentClass;\n");
    printf (OUT "  ObitTable".$ExtnameClean." *in = inn;\n");
    printf (OUT "\n");
    printf (OUT "  /* error checks */\n");
    printf (OUT "  g_assert (ObitIsA(in, &myClassInfo));\n");
    printf (OUT "\n");
    printf (OUT "  /* delete this class members */\n");
    printf (OUT "  \n");
    printf (OUT "  /* unlink parent class members */\n");
    printf (OUT "  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);\n");
    printf (OUT "  /* delete parent class members */\n");
    printf (OUT "  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) \n");
    printf (OUT "    ParentClass->ObitClear (inn);\n");
    printf (OUT "  \n");
    printf (OUT "} /* end ObitTable".$ExtnameClean."Clear */\n");
    printf (OUT "\n");

    printf (OUT "/**\n");
    printf (OUT " * Initialize global ClassInfo Structure.\n");
    printf (OUT " */\n");
    printf (OUT "void ObitTable".$ExtnameClean."ClassInit (void)\n");
    printf (OUT "{\n");
    printf (OUT "  if (myClassInfo.initialized) return;  /* only once */\n");
    printf (OUT "  \n");
    printf (OUT "  /* Set name and parent for this class */\n");
    printf (OUT "  myClassInfo.ClassName   = g_strdup(myClassName);\n");
    printf (OUT "  myClassInfo.ParentClass = ObitParentGetClass();\n");
    printf (OUT "\n");
    printf (OUT "  /* Set function pointers */\n");
    printf (OUT "  ObitTable".$ExtnameClean."ClassInfoDefFn ((gpointer)&myClassInfo);\n");
    printf (OUT " \n");
    printf (OUT "  myClassInfo.initialized = TRUE; /* Now initialized */\n");
    printf (OUT " \n");
    printf (OUT "} /* end ObitTable".$ExtnameClean."ClassInit */\n");
    printf (OUT "\n");
    printf (OUT "/**\n");
    printf (OUT " * Initialize global ClassInfo Function pointers.\n");
    printf (OUT " */\n");
    printf (OUT "static void ObitTable".$ExtnameClean."ClassInfoDefFn (gpointer inClass)\n");
    printf (OUT "{\n");
    printf (OUT "  ObitTable".$ExtnameClean."ClassInfo *theClass = (ObitTable".$ExtnameClean."ClassInfo*)inClass;\n");
    printf (OUT "  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;\n");
    printf (OUT "\n");
    printf (OUT "  if (theClass->initialized) return;  /* only once */\n");
    printf (OUT "\n");
    printf (OUT "  /* Check type of inClass */\n");
    printf (OUT "  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));\n");
    printf (OUT "\n");
    printf (OUT "  /* Initialize (recursively) parent class first */\n");
    printf (OUT "  if ((ParentClass!=NULL) && \n");
    printf (OUT "      (ParentClass->ObitClassInfoDefFn!=NULL))\n");
    printf (OUT "    ParentClass->ObitClassInfoDefFn(theClass);\n");
    printf (OUT "\n");
    printf (OUT "  /* function pointers defined or overloaded this class */\n");
    printf (OUT "  theClass->ObitClassInit = (ObitClassInitFP)ObitTable".$ExtnameClean."ClassInit;\n");
    printf (OUT "  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitTable".$ExtnameClean."ClassInfoDefFn;\n");
    printf (OUT "  theClass->ObitGetClass  = (ObitGetClassFP)ObitTable".$ExtnameClean."GetClass;\n");
    printf (OUT "  theClass->newObit       = (newObitFP)newObitTable".$ExtnameClean.";\n");
    printf (OUT "  theClass->ObitCopy      = (ObitCopyFP)ObitTable".$ExtnameClean."Copy;\n");
    printf (OUT "  theClass->ObitClone     = (ObitCloneFP)ObitTableClone;\n");
    printf (OUT "  theClass->ObitClear     = (ObitClearFP)ObitTable".$ExtnameClean."Clear;\n");
    printf (OUT "  theClass->ObitInit      = (ObitInitFP)ObitTable".$ExtnameClean."Init;\n");
    printf (OUT "  theClass->ObitTableConvert = (ObitTableConvertFP)ObitTable".$ExtnameClean."Convert;\n");
    printf (OUT "  theClass->ObitTableOpen    = (ObitTableOpenFP)ObitTable".$ExtnameClean."Open;\n");
    printf (OUT "  theClass->ObitTableClose   = (ObitTableCloseFP)ObitTable".$ExtnameClean."Close;\n");
    printf (OUT "  theClass->ObitTableRead    = (ObitTableReadFP)ObitTableRead;\n");
    printf (OUT "  theClass->ObitTableReadSelect = \n");
    printf (OUT "    (ObitTableReadSelectFP)ObitTableReadSelect;\n");
    printf (OUT "  theClass->ObitTableWrite = (ObitTableWriteFP)ObitTableWrite;\n");
    printf (OUT "  theClass->ObitTableReadRow = \n");
    printf (OUT "    (ObitTableReadRowFP)ObitTable".$ExtnameClean."ReadRow;\n");
    printf (OUT "  theClass->ObitTableWriteRow = \n");
    printf (OUT "    (ObitTableWriteRowFP)ObitTable".$ExtnameClean."WriteRow;\n");
    printf (OUT "\n");
    printf (OUT "} /* end ObitTable".$ExtnameClean."ClassDefFn */\n");
    printf (OUT "\n");
    printf (OUT "/**\n");
    printf (OUT " * Get table specific information from the infolist or descriptor\n");
    printf (OUT " * \x5cparam info Table to update\n");
    printf (OUT " * \x5cparam err  ObitErr for reporting errors.\n");
    printf (OUT " */\n");
    printf (OUT "static void ObitTable".$ExtnameClean."Update (ObitTable".$ExtnameClean." *in, ObitErr *err)\n");
    printf (OUT "{\n");
    printf (OUT "  olong i;\n");
   if ($NumKey>0) { # any keywords?
	printf (OUT "  ObitInfoType type;\n");
	printf (OUT "  gint32 dim[MAXINFOELEMDIM];\n");
    }
    printf (OUT "  ObitTableDesc *desc;\n");
    # Is a union needed?
    $temp = " ";
    for ($k=0; $k<$NumKey; $k++) {
	# For real values , be able to handle either float or double
	if (($KeyType[$k] eq "E") || ($KeyType[$k] eq "D")) {
	    $temp = "  union ObitInfoListEquiv InfoReal;";
	}
    }
    printf (OUT " $temp \n");
    printf (OUT "\n");
    printf (OUT " /* error checks */\n");
    printf (OUT "   g_assert(ObitErrIsA(err));\n");
    printf (OUT "  if (err->error) return;\n");
    printf (OUT "  g_assert (ObitIsA(in, &myClassInfo));\n");
    printf (OUT "\n");
    printf (OUT "  /* Get Keywords */\n");

    # Get Keywords
    for ($k=0; $k<$NumKey; $k++) {
	# For real values , be able to handle either float or double
	if (($KeyType[$k] ne "E") && ($KeyType[$k] ne "D")) {
	    # non real
	    printf (OUT "   /* $KeyName[$k] */\n");
	    # use default is keyword missing and default given
	    #if (($KeyDefV[$k] ne "") && ($KeyType[$k] eq "J")) {
	    if ($KeyDefV[$k] ne "")  {
		if ($KeyType[$k] eq "A") {
		    printf (OUT "  strncpy (in->$KeyName[$k], $KeyDefV[$k], MAXKEYCHARTABLE$ExtnameClean); \n");
		} else {
		    printf (OUT "  in->$KeyName[$k] = $KeyDefV[$k]; \n");
		}
		printf (OUT "  ObitInfoListGetTest(in->myDesc->info, $KeyWord[$k], &type, dim, \n");
		printf (OUT "		       (gpointer)&in->$KeyName[$k]);\n");
	    } else {
		printf (OUT "  if (!ObitInfoListGet(in->myDesc->info, $KeyWord[$k], &type, dim, \n");
		printf (OUT "		       (gpointer)&in->$KeyName[$k], err)) return;\n");
	    }
	} else {
	    # a real value use union and allow value in info list to be float or double
	    if ($KeyDefV[$k] ne "")  {
		printf (OUT "  in->$KeyName[$k] = $KeyDefV[$k]; \n");
		printf (OUT "  InfoReal.dbl = $KeyDefV[$k]; \n");
		printf (OUT "  type = OBIT_double;\n");
		printf (OUT "  ObitInfoListGetTest(in->myDesc->info, $KeyWord[$k], &type, dim, \n");
		printf (OUT "		       (gpointer)&InfoReal);\n");
	    } else {
		printf (OUT "   /* $KeyName[$k] */\n");
		printf (OUT "  if (!ObitInfoListGet(in->myDesc->info, $KeyWord[$k], &type, dim, \n");
		printf (OUT "		       (gpointer)&InfoReal, err)) return;\n");
	    }
	    printf (OUT "  if (type==OBIT_double) in->$KeyName[$k] = InfoReal.dbl;\n");
	    printf (OUT "  if (type==OBIT_float)  in->$KeyName[$k] = InfoReal.flt;\n");
	}
    }
    printf (OUT "\n");
    printf (OUT "  /* initialize column numbers/offsets */\n");

    # Initialize columns
    # Unindexed columns
    for ($k=0; $k<$NumCol; $k++) {
	if ($ColSuff[$k] eq "") {
	    printf (OUT "  in->$ColVName[$k]Off = -1;\n");
	    printf (OUT "  in->$ColVName[$k]Col = -1;\n");
	}
    }

    # indexed - loop over keywords
    for ($k=0; $k<$NumKey; $k++) {
	# it must have a range to be useful
	($lo, $hi) = split(',', substr($KeyRang[$k], 1));
	if (($lo>=1) && ($hi>=$lo)) {
	    # Loop over index values
	    for ($j=$lo; $j<=$hi; $j++) {

		# Loop over columns
		for ($i=0; $i<$NumCol; $i++) {
		    $colname =  $ColVName[$i]."$j";
		    $count = GetColumnCount($ColDim[$i]);
		    #print (" col $i $colname $ColSuff[$i] $KeyWord[$k] \n");

		    # do we want this one? indexed by keyword $k
		    $want = $ColSuff[$i] eq NoQuote($KeyWord[$k]);
		    
		    if ($want) {
			printf (OUT "  in->$colname\Off = -1;\n");
			printf (OUT "  in->$colname\Col = -1;\n");
		    }
		} # end loop over columns
		
	    } # end loop over indices
	}
    } # end loop over keywords
    

    printf (OUT "  /* Find columns and set offsets */\n");
    printf (OUT "  desc = in->myDesc;\n");
    printf (OUT "  if (desc->FieldName) {\n");
    printf (OUT "    for (i=0; i<desc->nfield; i++) {\n");

    # Unindexed columns
    for ($k=0; $k<$NumCol; $k++) {
	if ($ColSuff[$k] eq "") {
	    $sl = length($ColName[$k])-2; # Length of string minus quotes
	    printf (OUT "      if (!strncmp (desc->FieldName[i], $ColName[$k], $sl)) {\n");
	    printf (OUT "	 in->$ColVName[$k]Off = desc->offset[i];\n");
	    printf (OUT " 	 in->$ColVName[$k]Col = i;\n");
	    printf (OUT "      }\n");
	}
    }

    # indexed - loop over keywords
    for ($k=0; $k<$NumKey; $k++) {
	# it must have a range to be useful
	($lo, $hi) = split(',', substr($KeyRang[$k], 1));
	if (($lo>=1) && ($hi>=$lo)) {
	    # Loop over index values
	    for ($j=$lo; $j<=$hi; $j++) {

		# Loop over columns
		for ($i=0; $i<$NumCol; $i++) {
		    $colname =  $ColName[$i]."$j\"";
		    $count = GetColumnCount($ColDim[$i]);
		    #print (" col $i $colname $ColSuff[$i] $KeyWord[$k] \n");

		    # do we want this one? indexed by keyword $k
		    $want = $ColSuff[$i] eq NoQuote($KeyWord[$k]);
		    
		    if ($want) {
			$sl = length($colname)-2; # Length of string minus quotes
			printf (OUT "      if (!strncmp (desc->FieldName[i], $colname, $sl)) {\n");
			printf (OUT "	     in->$ColVName[$i]$j\Off = desc->offset[i];\n");
			printf (OUT " 	     in->$ColVName[$i]$j\Col = i;\n");
			printf (OUT "      }\n");
		    }
		} # end loop over columns
		
	    } # end loop over indices
	}
    } # end loop over keywords

    printf (OUT "     }\n");
    printf (OUT "  }\n");
    printf (OUT "\n");
    printf (OUT "  /* Check required columns */\n");
    # Make sure everything found - unindexed columns - last is status
    for ($k=0; $k<($NumCol-1); $k++) {
	# That which is not optional is required
	if (($ColSuff[$k] eq "") && (index($ColComm[$k],"[OPTIONAL]")<0)) {
	    printf (OUT "  Obit_return_if_fail((in->$ColVName[$k]Off > -1), err,\n");
	    printf (OUT "       \"ObitTable".$ExtnameClean."Update: Could not find column ".$ColVName[$k]."\");\n");
	}
    } # end loop checking if found
    # Make sure everything found - indexed columns
    for ($k=0; $k<$NumKey; $k++) {
 	# it must have a range to be useful
	($lo, $hi) = split(',', substr($KeyRang[$k], 1));
	if (($lo>=1) && ($hi>=$lo)) {
	    # Loop over index values
	    for ($j=$lo; $j<=$hi; $j++) {

		# Make sure range appropriate
		printf (OUT "  if (in->".$KeyName[$k].">=$j) {\n");

		# Loop over columns
		for ($i=0; $i<($NumCol-1); $i++) {
		    $colname =  $ColVName[$i]."$j";
		    $count = GetColumnCount($ColDim[$i]);

		    # do we want this one? indexed by keyword $k
		    $want = $ColSuff[$i] eq NoQuote($KeyWord[$k]);
		    # That which is not optional is required
		    $want = $want && (index($ColComm[$i],"[OPTIONAL]")<0);
		    
		    if ($want) {
			printf (OUT "    Obit_return_if_fail((in->$colname\Off > -1), err,\n");
			printf (OUT "       \"ObitTable".$ExtnameClean."Update: Could not find column ".$colname."\");\n");
		    }
		} # end loop over columns
		
		printf (OUT "  } /* end of if index in range */\n");
	    } # end loop over indices
	}
    } # end loop over keywords checking if found
    printf (OUT "} /* end ObitTable".$ExtnameClean."Update */\n");
    printf (OUT "\n");
    printf (OUT "/**\n");
    printf (OUT " * Copy table specific (keyword) information  to infolist.\n");
    printf (OUT " * \x5cparam info Table to update\n");
    printf (OUT " * \x5cparam err  ObitErr for reporting errors.\n");
    printf (OUT " */\n");
    printf (OUT "static void ObitTable".$ExtnameClean."DumpKey (ObitTable".$ExtnameClean." *in, ObitErr *err)\n");
    printf (OUT "{\n");
    printf (OUT "  ObitInfoList *info=NULL;\n");
    if ($NumKey>0) { # any keywords?
	printf (OUT "  ObitInfoType type;\n");
	printf (OUT "  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};\n");
    }
    printf (OUT "\n");
    printf (OUT " /* error checks */\n");
    printf (OUT "   g_assert(ObitErrIsA(err));\n");
    printf (OUT "  if (err->error) return;\n");
    printf (OUT "  g_assert (ObitIsA(in, &myClassInfo));\n");
    printf (OUT "\n");
    printf (OUT "  /* Set Keywords */\n");
    printf (OUT "  if (in->myIO!=NULL) info = ((ObitTableDesc*)(in->myIO->myDesc))->info;\n");
    printf (OUT "  else info = in->myDesc->info;\n");
    for ($k=0; $k<$NumKey; $k++) {
	printf (OUT "  /* $KeyName[$k] */\n");
	$TC = GetTypeCode($KeyType[$k]);
	printf (OUT "  type  = $TC;\n");
	if ($KeyType[$k] eq "A") { # special case for strings
	    printf (OUT "  dim[0] = MAXKEYCHARTABLE".$ExtnameClean.";\n");
	} else {
	    printf (OUT "  dim[0] = 1;\n");
	}
	printf (OUT "  ObitInfoListAlwaysPut(info, $KeyWord[$k], type, dim, \n");
	printf (OUT "		  (gpointer)&in->$KeyName[$k]);\n");
    }
    printf (OUT "   \n");
    printf (OUT "} /* end ObitTable".$ExtnameClean."DumpKey */\n");
} # end WriteTablePrivate

# Write Class function definition file
sub WriteSource {
    Copyright;
    WritePreamble;
    WriteTableRow;
    WriteConstructor;
    WriteCopy;
    WriteOpen;
    WriteReadRow;
    WriteSetRow;
    WriteWriteRow;
    WriteTableRowPrivate;
    WriteTablePrivate;
} # end WriteSource

# Write Class Python interface definition file
sub WritePy {
    Copyright;
    printf (OUT "\%{\n");
    printf (OUT "\#include \"Obit.h\"\n");
    printf (OUT "\#include \"ObitData.h\"\n");
    printf (OUT "\#include \"ObitTable".$ExtnameClean.".h\"\n");
    printf (OUT "\%}\n");
    printf (OUT " \n");
    printf (OUT "\%");  # Keep python from trashing this line
    printf (OUT "inline \%{\n");
    printf (OUT " \n");
    printf (OUT "extern ObitTable* Table".$ExtnameClean." (ObitData *inData, long *tabVer,\n");
    printf (OUT " 	                   int access,\n");
    printf (OUT " 	                   char *tabName,\n");

    # Add structural keywords to call sequence
    $temp = " ";
    for ($k=0; $k<$NumKey; $k++) {
	if ($KeyRang[$k] ne "") {
# redefine types to basal C types
	    $TC = DclTypeCode($KeyType[$k]);
	    if    ($KeyType[$k] eq "I"){$TC = "int";}
	    elsif ($KeyType[$k] eq "J"){$TC = "int";}
	    elsif ($KeyType[$k] eq "K"){$TC = "int";}
	    $temp = $temp." $TC $KeyName[$k],";
	}
    }
    printf (OUT "                         $temp\n");
    printf (OUT "                           ObitErr *err)\n");
    printf (OUT " {\n");
    printf (OUT "   ObitIOAccess laccess;\n");
    printf (OUT "   /* Cast structural keywords to correct type */\n");
    for ($k=0; $k<$NumKey; $k++) {
	if ($KeyRang[$k] ne "") {
	    $TC = DclTypeCode($KeyType[$k]);
	    printf (OUT "   $TC l$KeyName[$k] = ($TC)$KeyName[$k];\n");
	}
    }
    printf (OUT "   olong ltabVer = (olong)*tabVer;\n");
    printf (OUT "   ObitTable *outTable=NULL;\n");
    printf (OUT "   laccess = OBIT_IO_ReadOnly;\n");
    printf (OUT "   if (access==2) laccess = OBIT_IO_WriteOnly;\n");
    printf (OUT "   else if (access==3) laccess = OBIT_IO_ReadWrite;\n");
    printf (OUT "   outTable = (ObitTable*)newObitTable".$ExtnameClean."Value ((gchar*)tabName, inData, (olong*)&ltabVer,\n");
    printf (OUT "   			   laccess, \n");
    # Add structural keywords to call sequence
    $temp = " ";
    for ($k=0; $k<$NumKey; $k++) {
	if ($KeyRang[$k] ne "") {
	    $temp = $temp." l$KeyName[$k],";
	}
    }
    printf (OUT "                         $temp\n");
    printf (OUT "                           err);\n");
    printf (OUT "   *tabVer = (long)ltabVer;\n");
    printf (OUT "   return outTable;\n");
    printf (OUT "   }\n");
    printf (OUT " \n");

# TableXXGetHeadKeys
    printf (OUT "extern PyObject* Table".$ExtnameClean."GetHeadKeys (ObitTable *inTab) {\n");
    printf (OUT "  PyObject *outDict=PyDict_New();\n");
    printf (OUT "  ObitTable".$ExtnameClean." *lTab = (ObitTable".$ExtnameClean."*)inTab;\n");
    # Do structural keywords - assume all integers
    for ($k=0; $k<$NumKey; $k++) {
	if ($KeyRang[$k] ne "") {
	    printf (OUT "  PyDict_SetItemString(outDict, \"$KeyName[$k]\",  PyInt_FromLong((long)lTab->$KeyName[$k]));\n");
	}
    }

    # Do nonstructural keywords
    for ($k=0; $k<$NumKey; $k++) {
	if (($KeyRang[$k] eq "") && ($KeyDefV[$k] ne "")){
	    # By data type
	    if ($KeyType[$k] eq "A") {
		# String
		printf (OUT "  PyDict_SetItemString(outDict, \"$KeyName[$k]\", PyString_InternFromString(lTab->$KeyName[$k]));\n");
	    } elsif ($KeyType[$k] eq "K") {
		# Integer
		printf (OUT "  PyDict_SetItemString(outDict, \"$KeyName[$k]\",  PyInt_FromLong((long)lTab->$KeyName[$k]));\n");
	    } elsif ($KeyType[$k] eq "J") {
		# Integer
		printf (OUT "  PyDict_SetItemString(outDict, \"$KeyName[$k]\",  PyInt_FromLong((long)lTab->$KeyName[$k]));\n");
	    } elsif ($KeyType[$k] eq "D") {
		# Double
		printf (OUT "  PyDict_SetItemString(outDict, \"$KeyName[$k]\",  PyFloat_FromDouble((double)lTab->$KeyName[$k]));\n");
	    } elsif ($KeyType[$k] eq "E") {
		# Float
		printf (OUT "  PyDict_SetItemString(outDict, \"$KeyName[$k]\",  PyFloat_FromDouble((double)lTab->$KeyName[$k]));\n");
	    } elsif ($KeyType[$k] eq "L") {
		# Boolean
		printf (OUT " PyDict_SetItemString(outDict, \"$KeyName[$k]\",  PyInt_FromLong((long)lTab->$KeyName[$k]));\n");
	    }
	}
    }
    printf (OUT "\n");
    printf (OUT "  return outDict;\n");
    printf (OUT "} \n");
    printf (OUT "\n");


# TableXXSetHeadKeys
    printf (OUT "extern void Table".$ExtnameClean."SetHeadKeys (ObitTable *inTab, PyObject *inDict) {\n");
    printf (OUT "  ObitTable".$ExtnameClean." *lTab = (ObitTable".$ExtnameClean."*)inTab;\n");
    printf (OUT "  char *tstr;\n");
    printf (OUT "  int lstr=MAXKEYCHARTABLE".$ExtnameClean.";\n");
    printf (OUT "\n");
    # Do nonstructural keywords only
    for ($k=0; $k<$NumKey; $k++) {
	if (($KeyRang[$k] eq "") && ($KeyDefV[$k] ne "")){
	    # By data type
	    if ($KeyType[$k] eq "A") {
		# String
 		printf (OUT "  tstr = PyString_AsString(PyDict_GetItemString(inDict, \"$KeyName[$k]\"));\n");
		printf (OUT "  strncpy (lTab->$KeyName[$k], tstr, lstr); lTab->$KeyName[$k]\[lstr-1]=0;\n");
	    } elsif ($KeyType[$k] eq "K") {
		# Integer
 		printf (OUT "  lTab->$KeyName[$k] = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, \"$KeyName[$k]\"));\n");
	    } elsif ($KeyType[$k] eq "J") {
		# Integer
		printf (OUT "  lTab->$KeyName[$k] = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, \"$KeyName[$k]\"));\n");
	    } elsif ($KeyType[$k] eq "D") {
		# Double
		printf (OUT "  lTab->$KeyName[$k] = (odouble)PyFloat_AsDouble(PyDict_GetItemString(inDict, \"$KeyName[$k]\"));\n");
	    } elsif ($KeyType[$k] eq "E") {
		# Float
		printf (OUT "  lTab->$KeyName[$k] = (ofloat)PyFloat_AsDouble(PyDict_GetItemString(inDict, \"$KeyName[$k]\"));\n");
	    } elsif ($KeyType[$k] eq "L") {
		# Boolean
		printf (OUT "  lTab->$KeyName[$k] = (oint)PyInt_AsLong(PyDict_GetItemString(inDict, \"$KeyName[$k]\"));\n");
	    }
	}
    }
    printf (OUT "\n");
    printf (OUT "  if ((lTab->myDesc->access==OBIT_IO_ReadWrite) || (lTab->myDesc->access==OBIT_IO_WriteOnly)) \n");
    printf (OUT "    lTab->myStatus = OBIT_Modified;\n");
    printf (OUT "} \n");
    printf (OUT "\n");
    printf (OUT "\%}\n");
} # end WritePy
