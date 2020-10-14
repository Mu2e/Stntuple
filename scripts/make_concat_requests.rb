#!/usr/bin/env ruby
#------------------------------------------------------------------------------
#  concatenation input for STNTUPLE files
#  example: 
#  --------
#-----------------------------------------------------------------------
#  last check: 2004-08-01 - works out of the box .P.Murat
#
#  make_concat_requests.rb -i [host:]input_dir -b book -d dataset -p name_pattern \
#                          -o output_dir -t output_tcl -f mode [-m maxsize] [-l list_of_files]
#
# example:
# --------
#  cdfopr/scripts/make_concat_requests.rb \
#        -i fcdfdata034.fnal.gov:/cdf/scratch/cdfopr/datasets/cdfopr/sopr00/tmp \
#        -o ftp://ewk@fcdfdata030.fnal.gov:/cdf/scratch/cdfopr/datasets/cdfopr/sopr00/tmp \
#        -p express -b cdfopr -d sopr00
#
#  last check: 2004-08-19
#  defaults:
#  ---------
#  output_dir: fcdflnx3:/cdf/home/www/usr/cdfopr/datasets/$BOOK/$DATASET
#              (output tcl files go there)
#  
#  need default for user in output_dir URL
#  comment: it is important to remember the user name in the output_dir URL
#  modes: "stntuple", "dst"
#
#  parameters are parsed independently, they are globals, not data members of 
#  ConcatenationRequest objet
#------------------------------------------------------------------------------
puts "starting---"

require 'find'
require 'fileutils'
require 'getoptlong'

puts " emoe"
#-----------------------------------------------------------------------
def usage
  puts "usage: make_concat_requests -f format -h host -i input_dir "
  puts "                            -p name_pattern -o output_dir  "
  puts "                            -b book -d dataset -t output_tcl"
  exit(-1)
end
#------------------------------------------------------------------------------
# specify defaults for the global variables and parse command line options
#------------------------------------------------------------------------------
$input_dir      = ""
$output_dir     = ""
$output_tcl     = ""
$book           = "cdfpewk"
$dataset        = ""
$host           = "fcdfdata122.fnal.gov"
$pattern        = ""
$format         = "dst"
$iuser          = `whoami`.strip
$list_of_files  = nil
$max_file_size  = 1100000000
$ouser          = `whoami`.strip
$output_dir     = "fcdfsgi2"
$output_tcl_dir = ""
$ihost          = `hostname -f`.strip
$ohost          = `hostname -f`.strip
$local_host     = `hostname -f`.strip
$verbose        = nil

opts = GetoptLong.new(
  [ "--dataset"       , "-d",        GetoptLong::REQUIRED_ARGUMENT ],
  [ "--book"          , "-b",        GetoptLong::REQUIRED_ARGUMENT ],
  [ "--format"        , "-f",        GetoptLong::REQUIRED_ARGUMENT ],
  [ "--input_dir"     , "-i",        GetoptLong::REQUIRED_ARGUMENT ],
  [ "--list_of_files" , "-l",        GetoptLong::REQUIRED_ARGUMENT ], 
  [ "--max_file_size" , "-m",        GetoptLong::REQUIRED_ARGUMENT ],   # in megabytes
  [ "--output_tcl"    , "-t",        GetoptLong::REQUIRED_ARGUMENT ],
  [ "--output_dir"    , "-o",        GetoptLong::REQUIRED_ARGUMENT ],
  [ "--pattern"       , "-p",        GetoptLong::REQUIRED_ARGUMENT ],
  [ "--verbose"       , "-v",        GetoptLong::NO_ARGUMENT       ]
)

opts.each do |opt, arg|
  if    (opt == "--dataset"       ) ; $dataset        = arg
  elsif (opt == "--book"          ) ; $book           = arg
  elsif (opt == "--format"        ) ; $format         = arg
  elsif (opt == "--pattern"       ) ; $pattern        = arg
  elsif (opt == "--verbose"       ) ; $verbose        = 1
  elsif (opt == "--input_dir"     ) 
#-----------------------------------------------------------------------
#  input directory
#-----------------------------------------------------------------------
    if ( arg.split("@")[1] != nil )
      $iuser     = arg.split("@")[0];
      x          = arg.split("@")[1];
    else
      x          = arg;
    end
    if ( x.split(":")[1] != nil )   ; #  host and directory are specified
      $ihost     = x.split(":")[0].strip;
      $input_dir = x.split(":")[1].strip;
    else                            ; #  host is not specified
      $input_dir = x.strip;
    end
  elsif (opt == "--list_of_files" ) then
    $list_of_files = arg
  elsif (opt == "--max_file_size" ) ; then 
    $max_file_size = arg.to_i*1000000
  elsif (opt == "--output_dir"    ) 
#-----------------------------------------------------------------------
#  output directory
#-----------------------------------------------------------------------
    if ( arg.split("@")[1] != nil )
      $ouser     = arg.split("@")[0].strip;
      x          = arg.split("@")[1].strip;
    else
      x          = arg.strip;
    end
    if ( x.split(":")[1] != nil )   ; #  host and directory are specified
      $ohost      = x.split(":")[0].strip;
      $output_dir = x.split(":")[1].strip;
    else                            ; #  host is not specified
      $output_dir = x.strip;
    end
  elsif (opt == "--output_tcl") 
    $output_tcl=arg
  end

  if ($verbose != 0) ; puts "Option: #{opt}, arg #{arg.inspect}" ; end
end

usage if ($list_of_files == nil) && ($pattern == "")
#-----------------------------------------------------------------------
# derived global variables
#-----------------------------------------------------------------------
if ( $output_tcl.split("/")[1] == nil )
  $output_tcl = "cdfopr@fcdflnx4.fnal.gov:/cdf/home/www/usr/cdfopr/datasets/#{$book}/#{$dataset}/#{$output_tcl}"
end

if ($list_of_files == nil) then
  if ( $input_dir == "" ) 
    $input_dir = "/cdf/scratch/cdfopr/datasets/#{$book}/#{$dataset}"
  end

  if    ( $output_dir == "fcdfdata131" ) 
    $output_dir="murat@fcdfdata131.fnal.gov/export/data4/ewk/murat/datasets/#{$dataset}"
  end
end


#------------------------------------------------------------------------------
class ConcatenationRequest
  def initialize(format,max_size)

    if ($verbose) then ; puts "format : #{format} max_size=#{max_size}" ; end

    @format             = format;
    @pid              = `echo $$`.strip();
    @tmp_fn           = "/tmp/#{$pid}.tmp"
    @tmp_file         = File.open(@tmp_fn,"w");
    @output_file_name = `date +%Y_%m_%d.%H_%M_%S`.strip
    @max_size         = max_size
    @list             = Array.new
    @total_size       = 0
    @index            = 0
#-----------------------------------------------------------------------
#  here the execution starts: determine command to get directory catalog
#  and get it
#-----------------------------------------------------------------------
    puts "hostname = #{`hostname`.strip} ihost = #{$ihost}"
    if     ( $ihost == `hostname`.strip )                         # local host
      cmd = "ls -l #{$input_dir} | grep #{$pattern}"
    elsif  ( $ihost[0..7] == "fcdfdata" )                   # CAF fileserver
#       cmd="cafhostdir #{$ihost} #{$input_dir} | grep #{$pattern} | grep -v .log" ;

      cmd="ssh -l #{$iuser} #{$ihost} ls -l #{$input_dir} | grep #{$pattern} | grep -v .log" ;
    else  
      cmd="rsh -l #{$iuser} #{$ihost} ls -l #{$input_dir} | grep #{$pattern} | grep -v .log" ;
    end

    puts " [ConcatenationRequest::initialize] cmd = #{cmd}"

    if ($list_of_files) then
      if ($verbose); then puts 'list of files comes from '+$list_of_files; end
      f = File.new($list_of_files)
      @list_of_files = f.readlines();
      f.close()
    else
      @list_of_files = `#{cmd}`.split("\n");  
    end

    if ($verbose) then; puts @list_of_files ; end

    puts "@output_file_name = #{@output_file_name}"
    puts "@format           = #{@format}"
    puts "$ihost            = #{$ihost}"
    puts "$iuser            = #{$iuser}"
    puts "$ouser            = #{$ouser}"
    puts "$output_dir       = #{$output_dir}"
  end

#------------------------------------------------------------------------------
  def write_stntuple_header() 
    @tmp_file.puts "//----------------------------------------------"
    @tmp_file.puts "// NJOBS           __NJOBS__"
    @tmp_file.puts "// DATA_SERVER     root://#{$ihost}"
    @tmp_file.puts "// OUTPUT_DIR      ftp://#{$ouser}@#{$ohost}#{$output_dir}"
    @tmp_file.puts "// DATASET         #{$dataset}"
    @tmp_file.puts "// BOOK            #{$book}"
    @tmp_file.puts "//----------------------------------------------"
    @tmp_file.puts "void n_concatenation_jobs(int* N) {"
    @tmp_file.puts "  *N = __NJOBS__;"
    @tmp_file.puts "}"
    @tmp_file.puts "//----------------------------------------------"
    @tmp_file.puts "int init_chain(TChain* Chain, int JobNumber, TString* OutputDir, TString* Book, TString* Dataset) {"
  end

#------------------------------------------------------------------------------
  def write_stntuple_trailer() 
    @tmp_file.puts "  return 0;"
    @tmp_file.puts "}"
  end

#------------------------------------------------------------------------------
  def write_production_trailer() 
    @tmp_file.puts "  show include"
    @tmp_file.puts "exit"
    @tmp_file.puts ""
    @tmp_file.puts "path create NULL_PATH"
    @tmp_file.puts "talk FileOutput"
    @tmp_file.puts "  dhCache set      KAHUNA"
    @tmp_file.puts "  output  stream   AA #{$dataset}"
    @tmp_file.puts "  AA"
    @tmp_file.puts "    fileSize         set 1500000"
    @tmp_file.puts "    dataSetId        set #{$dataset}"
    @tmp_file.puts "    dataSetBook      set enemoe"
    @tmp_file.puts "    abortOnDBFailure set false"
    @tmp_file.puts "    compression      set t"
    @tmp_file.puts "    show"
    @tmp_file.puts "  exit"
    @tmp_file.puts "  output paths       AA NULL_PATH"
    @tmp_file.puts "exit"
    @tmp_file.puts ""
  end

#------------------------------------------------------------------------------
  def write_production_header() 
    @tmp_file.puts "#----------------------------------------------"
    @tmp_file.puts "# NJOBS           __NJOBS__"
    @tmp_file.puts "# DATA_SERVER     root://#{$ihost}"
    @tmp_file.puts "# OUTPUT_DIR      #{$output_dir}"
    @tmp_file.puts "# DATASET         #{$dataset}"
    @tmp_file.puts "# BOOK            #{$book}"
    @tmp_file.puts "#----------------------------------------------"
    @tmp_file.puts "module input DHInput"
    @tmp_file.puts "talk DHInput"
  end

#-----------------------------------------------------------------------
  def write_stntuple_request()
#    puts " --------------- write_request -----------------"
    @index = @index+1

    @tmp_file.puts "  if ( JobNumber == #{@index} ) {"
    @tmp_file.puts "//----------------------------------------------"
    @tmp_file.puts "// total size:     #{@total_size}"
    @tmp_file.puts "//----------------------------------------------"
    @tmp_file.puts "    *Dataset   = \"#{$dataset}\";"
    @tmp_file.puts "    *Book      = \"#{$book}\";"
    @tmp_file.puts "    *OutputDir = \"#{$output_dir}\";"
 
    nf = @list.length;
    for i in 0...nf
      @tmp_file.puts "    Chain->AddFile(\"#{@list[i]}\",TChain::kBigNumber);"
    end
    @tmp_file.puts "  }"

    @total_size = 0;
    @list.clear();
#    puts " --- exit from write_stntuple_request ---"
  end

#-----------------------------------------------------------------------
  def write_production_request()
#    puts " --------------- write_request -----------------"
    @index = @index+1

    @tmp_file.puts "  if { $env(JOB_NUMBER) == #{@index} } {"
    @tmp_file.puts "#----------------------------------------------"
    @tmp_file.puts "# total size:     #{@total_size}"
    @tmp_file.puts "#----------------------------------------------"
 
    nf = @list.length;
    for i in 0...nf
      if ($ihost.strip == $local_host)
        @tmp_file.puts "    include file #{@list[i]}"
      else 
        @tmp_file.puts "    include file #{$ihost}/#{@list[i]}"
      end
    end
    @tmp_file.puts "  }"

    @total_size = 0;
    @list.clear();
#    puts " --- exit from write_stntuple_request ---"
  end


#-----------------------------------------------------------------------
  def write_request()
    if (@format == "stntuple") 
      write_stntuple_request()
    elsif (@format == "dst") 
      write_production_request()
    end
  end


#-----------------------------------------------------------------------
  def write_header()
    if (@format == "stntuple") 
      write_stntuple_header()
    elsif (@format == "dst") 
      write_production_header()
    end
  end


#-----------------------------------------------------------------------
  def write_trailer()
    if (@format == "stntuple") 
      write_stntuple_trailer()
    elsif (@format == "dst") 
      write_production_trailer()
    end
  end


#------------------------------------------------------------------------------
  def make_request_file()
    if ($verbose) then ; puts "[make_request_file] : $verbose = #{$verbose} START" ; end

    write_header();

    if ($verbose) then ; puts "[make_request_file] : after write_header" ; end

                                           # loop over the files and do the job
    snew=0
    nfiles = @list_of_files.length
    for i in 0...nfiles
      if ($verbose) then ; puts @list_of_files[i]; end

      word = @list_of_files[i].strip.split(" ");
#  puts " word = #{word} #{word.length} #{word[4]}"
      name = word[8]
      size = word[4].to_i

      if ($input_dir != "") then 
        filename = "#{$input_dir}/#{name}"
      else
        filename = name
      end

#      if ($verbose) then puts " --- #{filename}" ; end

      if (size > @max_size) 
#-----------------------------------------------------------------------
#  large file
#-----------------------------------------------------------------------
        snew = @total_size + size;
        if (snew > @max_size)
#-----------------------------------------------------------------------
#  total size of the output chunk is above 1.5 GBytes, write 
#  out previous files
#-----------------------------------------------------------------------
          write_request();
          @list.push(filename);
          @total_size=size;
        else
#-----------------------------------------------------------------------
#  total size below $max_file_size GBytes, write out everything
#-----------------------------------------------------------------------
          @list.push(filename);
          @total_size=snew
          write_request();
        end
      else
#-----------------------------------------------------------------------
#  small file
#-----------------------------------------------------------------------
        snew = @total_size+size
        if (snew > @max_size) 
#-----------------------------------------------------------------------
#  small file, write request and put new file first into the list
#-----------------------------------------------------------------------
          write_request()
          @list.push(filename);
          @total_size=size
        else
          @list.push(filename);
          @total_size=snew;
        end
      end
    end

    if (@total_size != 0) ; write_request() ; end
    write_trailer();

    @tmp_file.close()
#-----------------------------------------------------------------------
#  make sure output directory exists and write njobs
#-----------------------------------------------------------------------
    @tmp_file   = File.open(@tmp_fn,"r");
    output_file = File.open(@tmp_fn+".1","w");

    puts "NJOBS:#{@index}"
    @tmp_file.each_line { |line|
      new_line = line.gsub("__NJOBS__","#{@index}");
      output_file.puts new_line
    }

    output_file.close();
    @tmp_file.close();
#-----------------------------------------------------------------------
#  copy output file to its final destination
#-----------------------------------------------------------------------
   puts "$output_tcl = #{$output_tcl}"
   if ($output_tcl.split(":")[1] != nil) 
#-----------------------------------------------------------------------
#  copy to a remote node
#-----------------------------------------------------------------------
     puts "cmd =  rcp  #{@tmp_fn}.1 #{$output_tcl}"
     rc        = `rcp  #{@tmp_fn}.1 #{$output_tcl}`
   else
#-----------------------------------------------------------------------
#  local file
#-----------------------------------------------------------------------
     dir = `dirname #{@output_file_name}`.strip
 
#     puts "..dir = #{dir}"

     if (! FileTest.exists?(dir) ) 
       puts " creating #{dir}"
       rc = `mkdir -p #{dir}`
       puts " rc = #{rc}"
     end

     rc = `mv #{@tmp_fn}.1 #{$output_tcl}`

   end

  end
#------------------------------------------------------------------------------
#  end of class definition
#------------------------------------------------------------------------------
end


if __FILE__ == $PROGRAM_NAME
  puts "emoe! #{$PROGRAM_NAME}"
  if ($verbose) then puts(">>>>>>>>>>>>>> creating concatenation request"); end

  req = ConcatenationRequest.new($format,1500000000);

  req.make_request_file();
end

exit(0)
