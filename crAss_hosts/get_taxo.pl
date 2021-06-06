use strict;
use warnings;
use LWP::Protocol::https;
use Mozilla::CA;
use Bio::DB::GenBank;
use Bio::Annotation::Collection;
use Bio::SeqIO;
use Data::Dumper;


my @IDS = @ARGV;


my $db = Bio::DB::GenBank -> new();


print "genome_id\tNCBI_taxid\ttaxonomy\n";

while (@IDS) {

    my @ids = splice @IDS, 0, 50;

	my $seqio = $db -> get_Stream_by_version( [@ids] );

	while( my $entry = $seqio -> next_seq ) {

		my $genome_id = $entry -> id . "." . $entry -> version;


		my $taxonomy = join ";", reverse($entry -> species -> classification);


        my $NCBI_taxid = "";

        for my $feat_object ($entry->get_SeqFeatures) {

            if (($feat_object -> primary_tag eq "source") && ($feat_object -> has_tag("db_xref"))) {

                my @arr = $feat_object -> get_tag_values("db_xref");

                $NCBI_taxid = join(";", @arr);

            }

        }


        print $genome_id . "\t" . $NCBI_taxid . "\t" . $taxonomy . "\n";

	}

}
