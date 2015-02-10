# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

=head1 NAME

Bio::EnsEMBL::Analysis::Tools::LayerFilter::AbstractLayerFilter

=head1 DESCRIPTION

This module is an "abstract" module to filter low quality gene models against
higher quality genes models. NEVER instanciate this module but either use one
of the existing module like GenomeOverlapFilter, AllExonOverlapFilter or
CodingExonOverlapFilter. At the moment CodingExonOverlapFilter is the prefered
filter.

To create a new filter module, you need to
 - inherit from this module: Bio::EnsEMBL::Analysis::Tools::LayerFilter::AbstractLayerFilter
 - implement has_overlap which will define the filtering.
    the function return 1 if the model has to be reject and 0 if the model has to be kept

Previously the filter was rejecting all models overlapping higher layer models.
Now we accept the models from lower layer overlapping higher layer models when the exons
boundaries are the same.

=head1 CONTACT

Post questions to the Ensembl development list: http://lists.ensembl.org/mailman/listinfo/dev

=cut

package Bio::EnsEMBL::Analysis::Tools::LayerFilter::AbstractLayerFilter;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );



sub new{
  my ($class, @args) = @_;
  my $self = bless {},$class;

  if (scalar(@args)) {
    throw("$class should have no args in new");
  }

  return $self;
}

=head2 has_overlap

 Arg [1]    : $upper, object to filtered against
 Arg [1]    : $lower, object to be filtered
 Example    : $self->has_overlap($upper, $lower);
 Description: All inheriting module must implement has_overlap. It determines if the object from the lower layer is
              overlapping objects from the upper layer. It can call the is_not_equal function which will allow the
              lower layer object to be carried on when it is similar to the upper layer object.
 Returntype : Integer, 1 when the model overlaps the other. 0 if the model does not overlap. Any other value
              for a model that overlaps but has to be kept
 Exceptions : None


=cut

sub has_overlap {
    throw("You need to implement a method has_overlap!");
}

=head2 is_not_equal

 Arg [1]    : $upper, array ref of object to filtered against
 Arg [1]    : $lower, array ref object to be filtered
 Example    : $self->is_not_equal($upper, $lower);
 Description: It compares the start and end of 2 array refs of features to determine if the objects have the same boundaries.
              The objects need ta have the seq_region_start and seq_region_end functions. Mostly it will be Bio::EnsEMBL:Exon objects.
 Returntype : Integer, 1 if they differ, 2 if they are equal
 Exceptions : None


=cut

sub is_not_equal {
    my ($self, $objects, $discarded_objects) = @_;

    # If the number of objects differs, they are not equal
    return 1 if (scalar(@$objects) != scalar(@$discarded_objects));
    my @sorted_objects = sort {$a->seq_region_start <=> $b->seq_region_start} @$objects;
    my @sorted_disc_objects = sort {$a->seq_region_start <=> $b->seq_region_start} @$discarded_objects;
    for (my $i = 0; $i < @sorted_objects; $i++) {
    # If the boundaries for one object differs, they are not equal
        return 1 if ($sorted_objects[$i]->seq_region_start != $sorted_disc_objects[$i]->seq_region_start or
                     $sorted_objects[$i]->seq_region_end != $sorted_disc_objects[$i]->seq_region_end);
    }
    return 2;
}

=head2 filter

 Arg [1]    : $lower_layer_genes, Array ref of Bio::EnsEMBL::Gene
 Arg [2]    : $upper_layer_genes, Array ref of Bio::EnsEMBL::Gene
 Example    : $self->filter($lower_layer_genes, $upper_layer_genes);
 Description: Filter the genes from the lower layer against the genes from the upper layer. This function will be called
              in the RunnableDB and returns the genes that can be kept. When creating a new filter module, the has_overlap
              function need to be implemented. The use of is_not_equal is prefered as it deepens the number of supporting
              evidences.
 Returntype : Array ref of Bio::EnsEMBL::Gene
 Exceptions : None


=cut

sub filter {
  my ($self, $lower_layer_genes, $upper_layer_genes) = @_;

  # assumption is that @upper_layer_genes is sorted by gene start

  my @filtered;

  my $cur_idx = 0;

  foreach my $lower_layer_gene (@$lower_layer_genes) {
    my (@genomic_overlap, $left_bound);


    for(my $i=$cur_idx; $i < @$upper_layer_genes; $i++) {
      my $upper_layer_gene = $upper_layer_genes->[$i];

      if ($upper_layer_gene->end >= $lower_layer_gene->start and not defined $left_bound) {
        $left_bound = $i;
      }

      if ($upper_layer_gene->end < $lower_layer_gene->start) {
        next;
      } elsif ($upper_layer_gene->start > $lower_layer_gene->end) {
        last;
      } else {
        push @genomic_overlap, $upper_layer_gene;
      }
    }

    $cur_idx = $left_bound if defined $left_bound;

    my $exon_overlap = 0;
    if (@genomic_overlap) {
        foreach my $upper_layer_gene (@genomic_overlap) {
            $exon_overlap = $self->has_overlap($upper_layer_gene, $lower_layer_gene);
            last if ($exon_overlap);
        }
    }

    # If there is no overlap, $exon_overlap should be zero. If there is an overlap and we want to discard the model,
    # exon_overlap should be 1. If there is an overlap but we want to keep the model, exon_overlap should be > 1
    if ($exon_overlap != 1) {
      push @filtered, $lower_layer_gene;
    }
  }

  return \@filtered;
}

1;
