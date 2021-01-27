use strict;
use warnings;
use LWP::Simple;

my $infile = $ARGV[0];

my @aux = split(/\//, $infile);

my $last_piece = pop @aux;

my $dir_name = join("_", $last_piece, "out_prot");

my $out_dir_path = join("/", @aux, $dir_name);

#if (-d $out_dir_path) {
#  die($!);
#} else {
#  system("mkdir $out_dir_path")
#}

chomp $infile;

open(IN, $infile);

my @ids;

my $i = 0;

my %ids;

while (my $line = <IN>) {
  chomp $line;
  my @aux = split (/\t/, $line);
  my $id = $aux[$#aux];
  next if (defined $ids{$id});
  $ids{$id} = 1;
#  $ids[$i] = $line;
#  $i++;
}

close IN;

foreach my $key(keys %ids) {
  $ids[$i] = $key;
  $i++;
}


#my $first = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=";

my $first = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=";

#https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=525912259&rettype=gb&retmode=text

my $base = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?";

my $last = "&rettype=fasta&retmode=text";

#my $ids = join(",", @ids);

my $tmp_file_name_index = 0;

my $root = "_tmp_";

my $tmp_file_name = $last_piece.$root.$tmp_file_name_index;

for (my $i = 0; $i <= $#ids; $i = $i + 199) {
  my @curr_ids = @ids[$i..($i+199)];
#  print "@curr_ids\n";
#  my $a = <STDIN>;
  my $ids = join(",", @curr_ids);
  if (-e "$last_piece.fa") {
    print "file $last_piece.fa already downloaded\n";
    last;
  } else {
    my $url = join ("", $first, $ids, $last);
    my $command = "wget -O \"$tmp_file_name\" \"$url\"";
    system $command;
    print $command."\n";
    $tmp_file_name_index++;
    $tmp_file_name = $last_piece.$root.$tmp_file_name_index;
    sleep 1;
    $i++;
  #  system ("mv $last.fa");
#  sleep 1;
  }
}

print("Joinning temporary files...\n");

my $command = ("cat $last_piece$root* > $last_piece\.aa.fa");
system $command;

my $command = ("rm $last_piece$root*");
system $command;

#https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=protein&term=291498256+OR+171850157&usehistory=y

#my @acc_array = split(/,/, $acc_list);

#assemble the esearch URL
#my $base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';

#assemble the efetch URL
#my $url = $base . "efetch.fcgi?db=protein&query_key=$key&WebEnv=$web";
#$url .= "&rettype=fasta&retmode=text";

#post the efetch URL
#my $fasta = get($url);
#print "$fasta";
