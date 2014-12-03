# Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

package Bio::EnsEMBL::Analysis::Tools::LayerFilter::AllExonOverlapFilter;

use strict;
use warnings;
use vars qw(@ISA);

use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Analysis::Tools::LayerFilter::AbstractLayerFilter;

@ISA = ('Bio::EnsEMBL::Analysis::Tools::LayerFilter::AbstractLayerFilter');

sub has_overlap {
    my ($self, $upper_layer_gene, $lower_layer_gene, $add_supporting_evidence) = @_;

    my @exons = @{$upper_layer_gene->get_all_Exons};
    foreach my $oe (@{$lower_layer_gene->get_all_Exons}) {
        foreach my $e (@exons) {
            if ( $oe->end >= $e->start and $oe->start <= $e->end) {
                return 1;
            }
        }
    }
    return 0;
}

1;
