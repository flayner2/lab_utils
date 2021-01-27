use strict;
use warnings;
use Bio::SeqIO;

my $seqio_obj = Bio::SeqIO->new(  -format => "fasta",
                                   -file => "$ARGV[0]"
                                );

my %gene_data; #gene-centered data
#my %locus_id_data;

while (my $seq_obj = $seqio_obj->next_seq) {
  my $seq = $seq_obj->seq();
  my $acc = $seq_obj->id();
  my @aux = split(/\|\|/, $acc);
  my $uniq_id;
  my $locus_id;
  my $gene_id;
  my $protein_id;
  my $flag = 0;
  foreach my $element (@aux) {
    if (($element =~ /^GCA_/)||($element =~ /^GCF_/)) {
#      print "Here 1\n";
      $uniq_id = $element;
      next;
    }
    if ($element =~ /gene:/) {
#      print "Here 2\n";
      $gene_id = $element;
      next;
    }
    elsif ($element =~ /locus_tag:/) {
#      print "Here 3\n";
      $locus_id = $element;
      next;
    }
    elsif ($element =~ /protein_id:/) {
#      print "Here 4\n";
      $protein_id = $element;
      next;
    }
    else {
#      print "Here 5\n";
      $flag = 1;
    }
  }
  if ($flag == 1) {
    next;
    print $acc."\n";
  }
  next if (! defined $protein_id);
  next if ((! defined $gene_id)&&(! defined $locus_id));
  if (defined $gene_id) {
    if (defined $gene_data{$gene_id}{gene_id}{sequence}) {
      my $actual_length = length($gene_data{$gene_id}{gene_id}{sequence});
      my $new_length = length($seq);
      if ($new_length > $actual_length) {
        $gene_data{$gene_id}{gene_id}{sequence} = $seq;
        $gene_data{$gene_id}{gene_id}{protein_id} = $protein_id;
      }
    }
    else {
      $gene_data{$gene_id}{gene_id}{sequence} = $seq;
      $gene_data{$gene_id}{gene_id}{protein_id} = $protein_id;
    }
  }
  if (defined $locus_id) {
    if (defined $gene_data{$locus_id}{locus_id}{sequence}) {
      my $new_length = length($seq);
      my $actual_length = length($gene_data{$locus_id}{locus_id}{sequence});
      if ($new_length > $actual_length) {
        $gene_data{$locus_id}{locus_id}{sequence} = $seq;
        $gene_data{$locus_id}{locus_id}{protein_id} = $protein_id;
      }
    }
    else {
      $gene_data{$locus_id}{locus_id}{sequence} = $seq;
      $gene_data{$locus_id}{locus_id}{protein_id} = $protein_id;
    }
  }
  
#  print "$acc\n";
}

my %print_flag;

foreach my $key (keys %gene_data) {
  my $tmp_prot_id;
  if (defined $gene_data{$key}{gene_id}{protein_id}) {
    $tmp_prot_id = $gene_data{$key}{gene_id}{protein_id};
  }
  if (defined $gene_data{$key}{locus_id}{protein_id}) {
    $tmp_prot_id = $gene_data{$key}{locus_id}{protein_id};
  }
  if (defined $print_flag{$tmp_prot_id}) {
    next;
  }
  else {
    if ((defined $gene_data{$key}{gene_id}{sequence})&&(!defined $print_flag{$tmp_prot_id})) {
      print ">$key|$gene_data{$key}{gene_id}{protein_id}\n$gene_data{$key}{gene_id}{sequence}\n";
      $print_flag{$tmp_prot_id} = 1;
      next;
    }
    if ((defined $gene_data{$key}{locus_id}{sequence})&&(!defined $print_flag{$tmp_prot_id})) {
      print ">$key|$gene_data{$key}{locus_id}{protein_id}\n$gene_data{$key}{locus_id}{sequence}\n";
      $print_flag{$tmp_prot_id} = 1;
      next;
    }
#    $print_flag{$key} = 1;
  }
}
