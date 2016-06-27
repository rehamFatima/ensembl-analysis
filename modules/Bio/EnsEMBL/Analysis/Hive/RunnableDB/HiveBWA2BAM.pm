# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016] EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
=pod

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::BWA




=head1 SYNOPSIS

my $runnableDB =  Bio::EnsEMBL::Analysis::RunnableDB::BWA->new( );

$runnableDB->fetch_input();
$runnableDB->run();

=head1 DESCRIPTION

This module uses BWA to process the alignment from the .sai files to create
a BAM file recording the pairing information if present

=head1 CONTACT

Post general queries to B<http://lists.ensembl.org/mailman/listinfo/dev>

=head1 APPENDIX
=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBWA2BAM;

use warnings ;
use strict;
use Bio::EnsEMBL::Analysis::Runnable::BWA2BAM;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


=head2 fetch_input

 Arg [1]    : None
 Description: Retrieve the .sai files generated by BWA, choose the algorithm based on the 'is_paired' key and
              create a Bio::EnsEMBL::Analysis::Runnable::BWA2BAM object to run bwa sampe/samse
 Returntype : None
 Exceptions : None

=cut

sub fetch_input {
  my ($self) = @_;
  my $program = $self->param('wide_short_read_aligner');
  $self->throw("BWA program not defined in analysis\n") unless (defined $program);
  my $fastqfile;
  my $fastqpair;

  my $method = $self->param('is_paired') ? ' sampe '.$self->param('sampe_options') : ' samse '.$self->param('samse_options');
  foreach my $fastq (@{$self->param('fastq')}) {
      my $abs_filename = $self->param('wide_input_dir').'/'.$fastq->{filename};
      $self->throw("Fastq file $abs_filename not found\n") unless (-e $abs_filename);
      if ($fastq->{is_mate_1} == 1) {
          $fastqfile = $abs_filename;
      }
      else {
          $fastqpair = $abs_filename;
      }
  }
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::BWA2BAM->new
    (
     -analysis   => $self->create_analysis,
     -program    => $program,
     -fastq      => $fastqfile,
     -fastqpair  => $fastqpair,
     -options    => $method,
     -outdir     => $self->param('wide_output_dir'),
     -genome     => $self->param('wide_genome_file'),
	 -samtools => $self->param('wide_samtools'),
     -header => $self->param('header_file'),
     -min_mapped => $self->param('min_mapped'),
     -min_paired => $self->param('min_paired'),
    );
    if ($self->param_is_defined('bam_prefix')) {
        $runnable->bam_prefix($self->param($self->param('bam_prefix')));
    }
  $self->runnable($runnable);
}


=head2 write_output

 Arg [1]    : None
 Description: Dataflow the absolute name of the BAM file, accessible via $self->param('filename') on branch 1
 Returntype : None
 Exceptions : None

=cut

sub write_output {
    my $self = shift;

    $self->dataflow_output_id({filename => $self->output->[0]}, 1);
}

1;