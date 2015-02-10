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

Bio::EnsEMBL::Analysis::Tools::LayerFilter::AllExonOverlapFilter

=head1 DESCRIPTION

This module low quality gene models against higher quality genes models based on all their
exons.
At the moment CodingExonOverlapFilter is the prefered filter.

To create a new filter module, you need to
 - inherit from this module: Bio::EnsEMBL::Analysis::Tools::LayerFilter::AbstractLayerFilter
 - implement has_overlap which will define the filtering.

Previously the filter was rejecting all models overlapping higher layer models.
Now we accept the models from lower layer overlapping higher layer models when the exons
boundaries are the same.

=head1 CONTACT

Post questions to the Ensembl development list: http://lists.ensembl.org/mailman/listinfo/dev

=cut

package Bio::EnsEMBL::Analysis::Tools::LayerFilter::AllExonOverlapFilter;

use strict;
use warnings;
use vars qw(@ISA);

use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Analysis::Tools::LayerFilter::AbstractLayerFilter;

@ISA = ('Bio::EnsEMBL::Analysis::Tools::LayerFilter::AbstractLayerFilter');

=head2 has_overlap

 Arg [1]    : $upper, object to filtered against
 Arg [1]    : $lower, object to be filtered
 Example    : $self->has_overlap($upper, $lower);
 Description: Filtering on exon overlap. All models from a lower level which has at least one exon overlapping
              an exon from a higher layer model will be discarded unless it is similar to the higher model
              (same coding exon boundaries).
 Returntype : Integer, 1 when the model overlaps the other. 0 if the model does not overlap. Any other value
              for a model that overlaps but has to be kept
 Exceptions : None


=cut

sub has_overlap {
    my ($self, $upper_layer_gene, $lower_layer_gene, $add_supporting_evidence) = @_;

    my $exons = $upper_layer_gene->get_all_Exons;
    foreach my $oe (@{$lower_layer_gene->get_all_Exons}) {
        foreach my $e (@$exons) {
            if ( $oe->end >= $e->start and $oe->start <= $e->end) {
                return $self->is_not_equal($exons, $lower_layer_gene->get_all_Exons);
            }
        }
    }
    return 0;
}

1;
